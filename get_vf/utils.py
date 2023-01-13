import argparse
import sys
import gzip
import os
import shutil
import logging
import pandas as pd
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull
from get_vf import __version__
from itertools import chain
from pathlib import Path
from get_vf.defaults import vfdbs, derep_bins

# from get_vf.defaults import markers, aln_formats, gtdb_releases
from Bio import SeqIO
import datetime


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


filters = ["breadth", "depth", "depth_evenness", "breadth_expected_ratio"]

# From https://stackoverflow.com/a/59617044/15704171
def convert_list_to_str(lst):
    n = len(lst)
    if not n:
        return ""
    if n == 1:
        return lst[0]
    return ", ".join(lst[:-1]) + f" or {lst[-1]}"


def check_filter_values(val, parser, var):
    value = str(val)
    if value in filters:
        return value
    else:
        parser.error(
            f"argument {var}: Invalid value {value}. Filter has to be one of {convert_list_to_str(filters)}"
        )


def check_values(val, minval, maxval, parser, var):
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            f"argument {var}: Invalid value value. Range has to be between {minval} and {maxval}!"
        )
    return value


def get_date_dir():
    current_date = datetime.datetime.now()
    folder_name = current_date.strftime("%Y%m%d_%H-%M-%S")
    return folder_name


# From: https://note.nkmk.me/en/python-check-int-float/
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


# function to check if the input value has K, M or G suffix in it
def check_suffix(val, parser, var):
    units = ["K", "M", "G"]
    unit = val[-1]
    value = int(val[:-1])

    if is_integer(value) & (unit in units) & (value > 0):
        return val
    else:
        parser.error(
            "argument %s: Invalid value %s. Memory has to be an integer larger than 0 with the following suffix K, M or G"
            % (var, val)
        )


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        parser.error("argument %s: The file %s does not exist!" % (var, arg))
    else:
        return arg


def is_valid_dir(parser, arg, var):
    if not os.path.exists(arg):
        parser.error("argument %s: The folder %s does not exist!" % (var, arg))
    else:
        return arg


def is_executable(parser, arg, var):
    """Check whether `name` is on PATH and marked as executable."""
    if shutil.which(arg):
        return arg
    else:
        parser.error("argument %s: %s not found in path." % (var, arg))


def is_in_vfdb(parser, arg, var):
    if arg in vfdbs:
        return arg
    else:
        parser.error(
            "argument %s: Invalid value %s.\nDB has to be one of %s"
            % (var, arg, convert_list_to_str(vfdbs))
        )


def is_in_derep_bins(parser, arg, var):
    derep_exec = Path(arg).name
    if derep_exec in derep_bins:
        is_executable(parser, arg, var)
        return arg
    else:
        parser.error(
            "argument %s: Invalid value %s.\De-replication binary has to be one of %s"
            % (var, arg, convert_list_to_str(derep_bins))
        )


def create_folder_and_delete(path):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        logging.info(f"Folder {path} already exists. Deleting.")
        try:
            shutil.rmtree(path)
            os.makedirs(path)
        except OSError as e:
            print("Error: %s : %s" % (path, e.strerror))


def delete_folder(path):
    if os.path.exists(path):
        logging.info(f"Deleting folder {path}.")
        try:
            shutil.rmtree(path)
        except OSError as e:
            print("Error: %s : %s" % (path, e.strerror))
    else:
        logging.info(f"Path {path} does not exist.")


def create_folder(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as e:
            print("Error: %s : %s" % (path, e.strerror))
    else:
        logging.info(f"Folder {path} already exists. Skipping creation.")


defaults = {
    "db": "core",
    "threads": 1,
    "outdir": "get-ancient-vf-output-fetch",
    "mmseqs_bin": "mmseqs",
    "extend_bin": "tadpole.sh",
    "extend_length": 100,
    "extend_memory": None,
    "extend_k": 17,
    "tmp": None,
    "derep_bin": "seqkit",
    "derep_min_length": 30,
    "createdb_release": "latest",
    "output": "get-ancient-vf-output-search",
    "db_dir": "get-ancient-vf-output-fetch",
    "mmseqs2_bin": "mmseqs",
    "mmseqs2_min_length": 15,
    "mmseqs2_min_seqid": 0.6,
    "mmseqs2_cov": 0.6,
    "mmseqs2_cov_mode": 2,
    "mmseqs2_evalue": 1e-5,
    "bitscore": 60,
    "evalue": 1e-5,
    "breadth": 0.5,
    "breadth_expected_ratio": 0.5,
    "depth": 0.1,
    "depth_evenness": 1.0,
    "prefix": None,
    "sort_memory": "1G",
    "mapping_file": None,
    "iters": 25,
    "scale": 0.9,
    "filter": "depth_evenness",
    "xfilter_bin": "xFilter",
    "extract_bin": "filterbyname.sh",
    "vfdb": "core",
}

help_msg = {
    "input": "A FASTA file containing the query sequences",
    "threads": "Number of threads to use",
    "outdir": "Output folder",
    "local_db": "Path to a local getVF DB",
    "vfdb": "Which VFDB to use",
    "hmmalign_bin": "Path to the hmmalign executable",
    "mmseqs_bin": "Path to the mmseqs2 executable",
    "no_extend": "Disable read extension",
    "extend_bin": "Path to the the executable for the extension step",
    "extend_length": "How much to extend in both ends",
    "extend_memory": "How much memory to use for the extension",
    "extend_k": "K-mer length to use for the extension",
    "createdb_release": "Release where to get the data",
    "createdb_output": "Output folder where to store the DB",
    "recreate": "Remove folders if they exist",
    "output": "Output folder to write the search results",
    "prefix": "Prefix used for the output files",
    "tmp": "Path to the temporary directory",
    "db_dir": "Folder with the VFDB data generated by the createdb subcommand",
    "no_derep": "Disable de-replication of identical reads",
    "derep_bin": "Path to the the executable for the de-replication step",
    "derep_min_length": "Minimum read length for the de-replication step",
    "mmseqs2_bin": "mmseqs",
    "ancient": "Use mmseqs2 aDNA optimized parameters",
    "mmseqs2_min_length": "Minimum length of the predicted ORF",
    "mmseqs2_min_seqid": "List matches above this sequence identity",
    "mmseqs2_cov": "List matches above this fraction of aligned (covered) residues",
    "mmseqs2_cov_mode": "Type of mmseq2 coverage calculation",
    "mmseqs2_evalue": " List matches below this E-value",
    "no_filter": "Disable x-filter filtering",
    "bitscore": "Bitscore where to filter the results",
    "evalue": "Evalue where to filter the results",
    "filter": "Which filter to use. Possible values are: breadth, depth, depth_evenness, breadth_expected_ratio",
    "xfilter_bin": "Path to the the executable for the x-filter step",
    "breadth": "Breadth of the coverage",
    "breadth_expected_ratio": "Expected breath to observed breadth ratio (scaled)",
    "depth": "Depth to filter out",
    "depth_evenness": "Reference with higher evenness will be removed",
    "mapping_file": "File with mappings to genes for aggregation",
    "iters": "Number of iterations for the FAMLI-like filtering",
    "scale": "Scale to select the best weithing alignments",
    "agg": "Aggregate the results by virulence factor",
    "help": "Help message",
    "debug": f"Print debug messages",
    "version": f"Print program version",
    "trim": f"Deactivate the trimming for the coverage calculations",
    "extract": f"Extract the reads that match the VFDB",
    "extract_bin": "Path to the the executable for the read extraction step",
    "keep": "Keep temporary files",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to filter BLASTx m8 files using the FAMLI algorithm",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )
    parser.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )

    # Same subparsers as usual
    sub_parsers = parser.add_subparsers(
        help="positional arguments",
        dest="action",
    )

    # Create parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = argparse.ArgumentParser(add_help=False)
    optional = parent_parser._action_groups.pop()

    optional.add_argument(
        "--threads",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
        ),
        dest="threads",
        default=1,
        help=help_msg["threads"],
    )
    # Create the parser sub-command for db creation
    parser_createdb = sub_parsers.add_parser(
        "createdb",
        help="Gather data from GTDB and create a DB",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # create the parser sub-command for search
    parser_search = sub_parsers.add_parser(
        "search",
        help="Search short reads against VFDB",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # createdb_required_args = parser_createdb.add_argument_group("required arguments")
    createdb_required_args = parser_createdb.add_argument_group(
        "Createdb required arguments"
    )
    createdb_optional_args = parser_createdb.add_argument_group(
        "Createdb optional arguments"
    )

    search_required_args = parser_search.add_argument_group("Search required arguments")
    search_optional_args = parser_search.add_argument_group("Search optional arguments")
    extend_args = parser_search.add_argument_group("Read extension arguments")
    derep_args = parser_search.add_argument_group("Read dereplication arguments")
    search_args = parser_search.add_argument_group("MMseqs2 search arguments")
    filter_args = parser_search.add_argument_group("X-filter processing arguments")
    extract_args = parser_search.add_argument_group(
        "Extract VFDB mapping reads arguments"
    )

    # createdb_required_args.add_argument(
    #     "-r",
    #     "--release",
    #     metavar="STR",
    #     type=lambda x: is_valid_gtdb_release(x, parser, "--release"),
    #     required=True,
    #     help=help_msg["createdb_release"],
    #     dest="gtdb_release",
    #     default=defaults["createdb_release"],
    # )
    createdb_output = get_date_dir()
    createdb_optional_args.add_argument(
        "--output",
        type=str,
        required=True,
        # default=createdb_output,
        metavar="STR",
        dest="createdb_output",
        help=help_msg["createdb_output"],
    )
    createdb_optional_args.add_argument(
        "--tmp",
        type=str,
        required=False,
        metavar="STR",
        default=None,
        help="Temporary directory location",
        dest="createdb_tmp",
    )
    createdb_optional_args.add_argument(
        "--recreate",
        dest="recreate",
        action="store_true",
        help=help_msg["recreate"],
        default=argparse.SUPPRESS,
    )
    createdb_optional_args.add_argument(
        "--mmseqs-bin",
        default=defaults["mmseqs_bin"],
        metavar="STR",
        type=lambda x: is_executable(parser, x, "--mmseqs-bin"),
        dest="mmseqs_bin",
        help=help_msg["mmseqs_bin"],
    )
    search_required_args.add_argument(
        "--input",
        required=True,
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser, x, "--input"),
        help=help_msg["input"],
    )
    search_required_args.add_argument(
        "--vfdb",
        required=True,
        metavar="STR",
        default=defaults["vfdb"],
        type=lambda x: is_in_vfdb(parser, x, "--vfdb"),
        dest="db",
        help=help_msg["vfdb"],
    )
    search_required_args.add_argument(
        "--db-dir",
        required=True,
        metavar="STR",
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser, x, "--db-dir"),
        dest="db_dir",
        help=help_msg["db_dir"],
    )
    search_optional_args.add_argument(
        "--output",
        type=str,
        default=defaults["output"],
        metavar="STR",
        dest="output",
        help=help_msg["output"],
    )
    search_optional_args.add_argument(
        "--tmp",
        type=lambda x: is_valid_file(parser, x, "--tmp"),
        default=defaults["tmp"],
        metavar="STR",
        dest="tmp",
        help=help_msg["tmp"],
    )
    search_optional_args.add_argument(
        "--prefix",
        type=str,
        metavar="STR",
        default=defaults["prefix"],
        dest="prefix",
        help=help_msg["prefix"],
    )
    search_optional_args.add_argument(
        "--no-extend",
        dest="no_extend",
        action="store_true",
        help=help_msg["no_extend"],
        default=argparse.SUPPRESS,
    )
    extend_args.add_argument(
        "--extend-bin",
        default=defaults["extend_bin"],
        metavar="STR",
        type=lambda x: is_executable(parser, x, "--extend-bin"),
        dest="extend_bin",
        help=help_msg["extend_bin"],
    )
    extend_args.add_argument(
        "--extend-length",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--extend-length")
        ),
        metavar="INT",
        default=defaults["extend_length"],
        dest="extend_length",
        help=help_msg["extend_length"],
    )
    extend_args.add_argument(
        "--extend-k",
        type=lambda x: int(
            check_values(x, minval=9, maxval=121, parser=parser, var="--extend-k")
        ),
        metavar="INT",
        default=defaults["extend_k"],
        dest="extend_k",
        help=help_msg["extend_k"],
    )
    extend_args.add_argument(
        "--extend-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--extend-memory"),
        default=defaults["extend_memory"],
        dest="extend_memory",
        metavar="INT",
        help=help_msg["extend_memory"],
    )
    search_optional_args.add_argument(
        "--no-derep",
        dest="no_derep",
        action="store_true",
        default=argparse.SUPPRESS,
        help=help_msg["no_derep"],
    )
    derep_args.add_argument(
        "--derep-bin",
        default=defaults["derep_bin"],
        metavar="STR",
        type=lambda x: is_in_derep_bins(parser, x, "--derep-bin"),
        dest="derep_bin",
        help=help_msg["derep_bin"],
    )
    # derep_args.add_argument(
    #     "--derep-min-length",
    #     type=lambda x: int(
    #         check_values(
    #             x, minval=1, maxval=100000, parser=parser, var="--derep-min-length"
    #         )
    #     ),
    #     default=defaults["derep_min_length"],
    #     dest="derep_min_length",
    #     help=help_msg["derep_min_length"],
    # )
    search_args.add_argument(
        "--mmseqs2-bin",
        default=defaults["mmseqs2_bin"],
        metavar="STR",
        type=lambda x: is_executable(parser, x, "--mmseqs2-bin"),
        dest="mmseqs_bin",
        help=help_msg["mmseqs_bin"],
    )
    search_args.add_argument(
        "--mmseq2-min-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=100000, parser=parser, var="--mmseqs2-min-length"
            )
        ),
        metavar="INT",
        default=defaults["mmseqs2_min_length"],
        dest="mmseqs_min_length",
        help=help_msg["mmseqs2_min_length"],
    )
    search_args.add_argument(
        "--mmseqs2-evalue",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1e6, parser=parser, var="--mmseqs2-evalue")
        ),
        metavar="FLOAT",
        default=defaults["mmseqs2_evalue"],
        dest="mmseqs_evalue",
        help=help_msg["evalue"],
    )
    search_args.add_argument(
        "--mmseq2-min-seqid",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1, parser=parser, var="--mmseqs2-min-seqid"
            )
        ),
        metavar="FLOAT",
        default=defaults["mmseqs2_min_seqid"],
        dest="mmseqs_min_seqid",
        help=help_msg["mmseqs2_min_seqid"],
    )
    search_args.add_argument(
        "--mmseq2-cov",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1, parser=parser, var="--mmseqs2-cov")
        ),
        metavar="FLOAT",
        default=defaults["mmseqs2_cov"],
        dest="mmseqs_cov",
        help=help_msg["mmseqs2_cov"],
    )
    search_args.add_argument(
        "--mmseq2-cov-mode",
        type=lambda x: int(
            check_values(x, minval=0, maxval=5, parser=parser, var="--mmseqs2-cov-mode")
        ),
        metavar="INT",
        default=defaults["mmseqs2_cov_mode"],
        dest="mmseqs_cov_mode",
        help=help_msg["mmseqs2_cov_mode"],
    )
    search_optional_args.add_argument(
        "--ancient",
        dest="ancient",
        action="store_true",
        help=help_msg["ancient"],
        default=argparse.SUPPRESS,
    )
    search_optional_args.add_argument(
        "--keep",
        dest="no_keep",
        action="store_true",
        help=help_msg["keep"],
        default=argparse.SUPPRESS,
    )
    search_optional_args.add_argument(
        "--no-filter",
        dest="no_filter",
        action="store_true",
        help=help_msg["no_filter"],
        default=argparse.SUPPRESS,
    )
    filter_args.add_argument(
        "--x-filter-bin",
        default=defaults["xfilter_bin"],
        metavar="STR",
        type=lambda x: is_executable(parser, x, "--x-filter-bin"),
        dest="xfilter_bin",
        help=help_msg["xfilter_bin"],
    )
    filter_args.add_argument(
        "--n-iters",
        type=lambda x: int(
            check_values(x, minval=1, maxval=100000, parser=parser, var="--n-iters")
        ),
        default=defaults["iters"],
        metavar="INT",
        dest="iters",
        help=help_msg["iters"],
    )
    filter_args.add_argument(
        "--evalue",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1e6, parser=parser, var="--evalue")
        ),
        default=defaults["evalue"],
        metavar="FLOAT",
        dest="evalue",
        help=help_msg["evalue"],
    )
    filter_args.add_argument(
        "--scale",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1, parser=parser, var="--scale")
        ),
        default=defaults["scale"],
        metavar="FLOAT",
        dest="scale",
        help=help_msg["scale"],
    )
    filter_args.add_argument(
        "--bitscore",
        type=lambda x: int(
            check_values(x, minval=0, maxval=1e6, parser=parser, var="--bitscore")
        ),
        default=defaults["bitscore"],
        metavar="INT",
        dest="bitscore",
        help=help_msg["bitscore"],
    )
    filter_args.add_argument(
        "--filter",
        type=lambda x: str(check_filter_values(x, parser=parser, var="--filter")),
        default=defaults["filter"],
        metavar="STR",
        dest="filter",
        help=help_msg["filter"],
    )
    filter_args.add_argument(
        "--breadth",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1, parser=parser, var="--breadth")
        ),
        default=defaults["breadth"],
        metavar="FLOAT",
        dest="breadth",
        help=help_msg["breadth"],
    )
    filter_args.add_argument(
        "--breadth-expected-ratio",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1, parser=parser, var="--breadth-expected-ratio"
            )
        ),
        default=defaults["breadth_expected_ratio"],
        metavar="FLOAT",
        dest="breadth_expected_ratio",
        help=help_msg["breadth_expected_ratio"],
    )
    filter_args.add_argument(
        "--depth",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1e6, parser=parser, var="--depth")
        ),
        default=defaults["depth"],
        metavar="FLOAT",
        dest="depth",
        help=help_msg["depth"],
    )
    filter_args.add_argument(
        "--depth-evenness",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1e6, parser=parser, var="--depth-evenness")
        ),
        default=defaults["depth_evenness"],
        metavar="FLOAT",
        dest="depth_evenness",
        help=help_msg["depth_evenness"],
    )
    filter_args.add_argument(
        "--no-trim",
        dest="trim",
        action="store_false",
        help=help_msg["trim"],
        default=argparse.SUPPRESS,
    )
    filter_args.add_argument(
        "--no-aggregate",
        dest="agg",
        action="store_true",
        help=help_msg["agg"],
        default=argparse.SUPPRESS,
    )
    extract_args.add_argument(
        "--extract",
        dest="extract",
        action="store_true",
        help=help_msg["extract"],
        default=argparse.SUPPRESS,
    )
    extract_args.add_argument(
        "--extract-bin",
        default=defaults["extract_bin"],
        metavar="STR",
        type=lambda x: is_executable(parser, x, "--extract-bin"),
        dest="extract_bin",
        help=help_msg["extract_bin"],
    )
    args = parser.parse_args(None if sys.argv[1:] else ["-h"])

    if not hasattr(args, "no_derep"):
        args.no_derep = False
    if not hasattr(args, "no_extend"):
        args.no_extend = False
    if not hasattr(args, "no_filter"):
        args.no_filter = False
    if not hasattr(args, "trim"):
        args.trim = True
    if not hasattr(args, "ancient"):
        args.ancient = False
    if not hasattr(args, "no_keep"):
        args.no_keep = True
    if not hasattr(args, "recreate"):
        args.recreate = False
    if not hasattr(args, "local_db"):
        args.local_db = False
    if not hasattr(args, "agg"):
        args.agg = True
    if not hasattr(args, "extract"):
        args.extract = False

    if args.action is not None and len(sys.argv) == 2:
        args = parser.parse_args([args.action, "-h"])

    return args


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def concat_df(frames):
    COLUMN_NAMES = frames[0].columns
    df_dict = dict.fromkeys(COLUMN_NAMES, [])
    for col in COLUMN_NAMES:
        extracted = (frame[col] for frame in frames)
        # Flatten and save to df_dict
        df_dict[col] = fast_flatten(extracted)
    df = pd.DataFrame.from_dict(df_dict)[COLUMN_NAMES]
    return df


def is_fastq(filename):
    with get_open_func(filename)(filename, "rt") as handle:
        fastq = SeqIO.parse(handle, "fastq")
        return any(fastq)  # False when `fasta` is empty, i.e. wasn't a FASTA file


def create_output_files(prefix, input, db, fastq):

    if prefix is None:
        prefix = Path(input).resolve().stem.split(".")[0]
    if fastq:
        extension = "fq"
    else:
        extension = "fa"

    # create output files
    out_files = {
        "extended_reads": f"{prefix}.extended.fq.gz",
        "extend_reads_log": f"{prefix}.extend.log",
        "derep_reads": f"{prefix}.derep.fq.gz",
        "derep_reads_log": f"{prefix}.derep.log",
        "results": f"{prefix}.results.{db}.tsv",
        "results_gz": f"{prefix}.results.{db}.tsv.gz",
        "results_log": f"{prefix}.results.{db}.log",
        "results_filtered_prefix": f"{prefix}.results-filtered.{db}",
        "results_filtered_log": f"{prefix}.results-filtered.{db}.log",
        "results_filtered_mm": f"{prefix}.results-filtered.{db}_no-multimap.tsv.gz",
        "results_filtered_cov": f"{prefix}.results-filtered.{db}_cov-stats.tsv.gz",
        "results_filtered_group_agg": f"{prefix}.results-filtered.{db}_group-abundances-agg.tsv.gz",
        "results_filtered_group": f"{prefix}.results-filtered.{db}_group-abundances.tsv.gz",
        "reads_db_log": f"{prefix}.reads.{db}.log",
        "reads_db": f"{prefix}.reads.{db}.{extension}.gz",
    }
    return out_files
