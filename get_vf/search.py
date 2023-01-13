import sys
import logging
import pathlib
import pandas as pd
import os
import requests
import http.client
from Bio import SeqIO, AlignIO, Seq
import numpy as np
from scipy import stats
import subprocess
import re
import io
from multiprocessing import Pool
from get_vf.utils import (
    create_folder,
    create_folder_and_delete,
    is_fastq,
    create_output_files,
    get_date_dir,
    delete_folder,
)
from get_vf.defaults import VFDB, vfdbs
import gzip
import shutil
from contextlib import ExitStack

log = logging.getLogger("my_logger")


def extend_reads(
    reads,
    output,
    extend_bin,
    extend_k,
    extend_length,
    extend_memory,
    tmp_dir,
    threads,
    db,
    log_file,
):
    if extend_memory is None:
        extend_memory = " "
    else:
        extend_memory = f"-Xmx{extend_memory}"

    output = pathlib.Path(tmp_dir, output)
    if output.exists():
        logging.info("Extending reads already found. Skipping.")
        return output
    logging.info(f"Extending reads [k:{extend_k}; length:{extend_length}]")
    proc = subprocess.Popen(
        [
            extend_bin,
            extend_memory,
            f"in={reads}",
            f"out={output}",
            f"k={extend_k}",
            f"mode=extend",
            f"el={extend_length}",
            f"er={extend_length}",
            f"overwrite=true",
            f"threads={threads}",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    stdout, stderr = proc.communicate()

    with open(log_file, "w") as f:
        for row in stderr.split("\n"):
            f.write(row + "\n")

    if proc.returncode != 0:
        # rgx = re.compile("ERROR")
        # for line in stdout.splitlines():
        #     if rgx.search(line):
        #         error = line.replace("ERROR: ", "")
        logging.error(f"Error extending reads for DB: {db}")
        logging.error(stderr)
        exit(1)
    return output


def dereplicate_reads(
    derep_bin,
    reads,
    output,
    tmp_dir,
    log_file,
    db,
    # derep_min_length,
    threads,
):
    # vsearch --derep_fulllength ${i} --output - --minseqlength 30 --strand both
    output = pathlib.Path(tmp_dir, output)
    if output.exists():
        logging.info("Dereplicated reads already found. Skipping.")
        return output
    if derep_bin == "vsearch":
        logging.info(f"Dereplicating reads with {derep_bin}")
        proc = subprocess.Popen(
            [
                derep_bin,
                "--fastx_uniques",
                reads,
                "--fastqout",
                output,
                "--strand",
                "both",
                # "--minseqlength",
                # str(derep_min_length),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf-8",
        )
        stdout, stderr = proc.communicate()

    elif derep_bin == "seqkit":
        logging.info(f"Dereplicating reads with {derep_bin}")
        proc = subprocess.Popen(
            [
                derep_bin,
                "rmdup",
                "-s",
                reads,
                "-j",
                str(threads),
                "-o",
                output,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf-8",
        )
        stdout, stderr = proc.communicate()

    with open(log_file, "w") as f:
        for row in stderr.split("\n"):
            f.write(row + "\n")
    if proc.returncode != 0:
        # rgx = re.compile("ERROR")
        # for line in stdout.splitlines():
        #     if rgx.search(line):
        #         error = line.replace("ERROR: ", "")
        logging.error(f"Error dereplicating reads for DB: {db}")
        logging.error(stderr)
        exit(1)
    return output


def search_mmseqs(
    reads,
    db,
    vfdb_db,
    output,
    tmp_dir,
    log_file,
    mmseqs_bin,
    mmseqs_min_length,
    mmseqs_evalue,
    mmseqs_min_seqid,
    mmseqs_cov,
    mmseqs_cov_mode,
    ancient,
    threads,
):
    tmp = pathlib.Path(tmp_dir, "mmseqs", db)
    create_folder(tmp)
    if output.exists():
        output.unlink()
    logging.info(f"Searching reads against DB {db}")
    if ancient:
        logging.info("::: Using MMseqs2 parameters for aDNA")
        cmd = [
            mmseqs_bin,
            "easy-search",
            reads,
            vfdb_db,
            output,
            tmp,
            "--min-length",
            str(mmseqs_min_length),
            "-e",
            str(mmseqs_evalue),
            "--min-seq-id",
            str(mmseqs_min_seqid),
            "-c",
            str(mmseqs_cov),
            "--cov-mode",
            str(mmseqs_cov_mode),
            "--format-mode",
            str(2),
            "--comp-bias-corr",
            str(0),
            "--mask",
            str(0),
            "--exact-kmer-matching",
            str(1),
            "--sub-mat",
            "PAM30.out",
            "-s",
            str(5),
            "-k",
            str(6),
            "--spaced-kmer-pattern",
            "11011101",
            "--threads",
            str(threads),
        ]
    else:
        cmd = [
            mmseqs_bin,
            "easy-search",
            reads,
            vfdb_db,
            output,
            tmp,
            "--min-length",
            str(mmseqs_min_length),
            "-e",
            str(mmseqs_evalue),
            "--min-seq-id",
            str(mmseqs_min_seqid),
            "-c",
            str(mmseqs_cov),
            "--cov-mode",
            str(mmseqs_cov_mode),
            "--format-mode",
            str(2),
            "--threads",
            str(threads),
        ]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    stdout, stderr = proc.communicate()

    with open(log_file, "w") as f:
        for row in stderr.split("\n"):
            f.write(row + "\n")

    if proc.returncode != 0:
        # rgx = re.compile("ERROR")
        # for line in stdout.splitlines():
        #     if rgx.search(line):
        #         error = line.replace("ERROR: ", "")
        logging.error(f"Error searching reads for DB: {db}")
        logging.error(stderr)
        exit(1)
    delete_folder(tmp)
    return output


def filter_results(
    xfilter_bin,
    results,
    prefix,
    db,
    log_file,
    n_iters,
    evalue,
    scale,
    bitscore,
    fltr,
    breadth,
    breadth_expected_ratio,
    depth,
    depth_evenness,
    trim,
    threads,
    output,
    metadata,
    aggregate,
):
    logging.info(f"Filtering VFDB {db} BLASTx results")
    if trim:
        if aggregate:
            cmd = [
                xfilter_bin,
                "--input",
                results,
                "--prefix",
                prefix,
                "--n-iters",
                str(n_iters),
                "--evalue",
                str(evalue),
                "--scale",
                str(scale),
                "--bitscore",
                str(bitscore),
                "--filter",
                str(fltr),
                "--breadth",
                str(breadth),
                "--breadth-expected-ratio",
                str(breadth_expected_ratio),
                "--depth",
                str(depth),
                "--depth-evenness",
                str(depth_evenness),
                "--threads",
                str(threads),
                "-m",
                metadata,
            ]
        else:
            cmd = [
                xfilter_bin,
                "--input",
                results,
                "--prefix",
                prefix,
                "--n-iters",
                str(n_iters),
                "--evalue",
                str(evalue),
                "--scale",
                str(scale),
                "--bitscore",
                str(bitscore),
                "--filter",
                str(fltr),
                "--breadth",
                str(breadth),
                "--breadth-expected-ratio",
                str(breadth_expected_ratio),
                "--depth",
                str(depth),
                "--depth-evenness",
                str(depth_evenness),
                "--threads",
                str(threads),
            ]
    else:
        if aggregate:
            cmd = [
                xfilter_bin,
                "--input",
                results,
                "--prefix",
                prefix,
                "--n-iters",
                str(n_iters),
                "--evalue",
                str(evalue),
                "--scale",
                str(scale),
                "--bitscore",
                str(bitscore),
                "--filter",
                str(fltr),
                "--breadth",
                str(breadth),
                "--breadth-expected-ratio",
                str(breadth_expected_ratio),
                "--depth",
                str(depth),
                "--depth-evenness",
                str(depth_evenness),
                "--threads",
                str(threads),
                "--no-trim",
                "-m",
                metadata,
            ]
        else:
            cmd = [
                xfilter_bin,
                "--input",
                results,
                "--prefix",
                prefix,
                "--n-iters",
                str(n_iters),
                "--evalue",
                str(evalue),
                "--scale",
                str(scale),
                "--bitscore",
                str(bitscore),
                "--filter",
                str(fltr),
                "--breadth",
                str(breadth),
                "--breadth-expected-ratio",
                str(breadth_expected_ratio),
                "--depth",
                str(depth),
                "--depth-evenness",
                str(depth_evenness),
                "--threads",
                str(threads),
                "--no-trim",
            ]
    proc = subprocess.Popen(
        cmd,
        cwd=output,
        stdout=subprocess.PIPE,
        # stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    stdout, stderr = proc.communicate()

    # with open(log_file, "w") as f:
    #     for row in stderr.split("\n"):
    #         f.write(row + "\n")

    if proc.returncode != 0:
        # rgx = re.compile("ERROR")
        # for line in stdout.splitlines():
        #     if rgx.search(line):
        #         error = line.replace("ERROR: ", "")
        logging.error(f"Error filtering BLASTx results for DB: {db}")
        # logging.error(stderr)
        exit(1)
    # gzip output file if exists


def extract_reads(results, db, output, threads, reads, extract_bin, tmp_dir, log_file):
    logging.info(f"Extracting reads for DB {db}")
    results = pd.read_csv(results, sep="\t")
    name_lst = pathlib.Path(tmp_dir, f"{db}-names.tsv")
    results[["queryId"]].to_csv(name_lst, sep="\t", index=False, header=False)
    cmd = [
        extract_bin,
        f"in={reads}",
        f"out={output}",
        f"names={name_lst}",
        f"overwrite=true",
        f"threads={threads}",
        f"include=t",
    ]

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )

    stdout, stderr = proc.communicate()

    with open(log_file, "w") as f:
        for row in stderr.split("\n"):
            f.write(row + "\n")

    if proc.returncode != 0:
        # rgx = re.compile("ERROR")
        # for line in stdout.splitlines():
        #     if rgx.search(line):
        #         error = line.replace("ERROR: ", "")
        logging.error(f"Error extracting reads for DB: {db}")
        logging.error(stderr)
        exit(1)
    return output


def search_db(args):

    if args.no_derep:
        derep = False
    else:
        derep = True

    if args.no_extend:
        extend = False
    else:
        extend = True

    if args.no_filter:
        fltr = False
    else:
        fltr = True

    vfdb_db = pathlib.Path(
        args.db_dir,
        args.db,
        "mmseqs",
        f"{args.db}-db",
    )

    metadata = pathlib.Path(args.db_dir, args.db, VFDB[args.db]["metadata"]).resolve()
    if not os.path.exists(vfdb_db):
        logging.error(
            f"The DB for the {args.db} VFDB does not exist. Please run the createdb subcommand to create it."
        )
        exit(1)
    else:
        logging.info(f"DB for the {args.db} VFDB found")

    logging.info(f"Creating output and tmp folders")
    date_dir = get_date_dir()

    output = pathlib.Path(args.output, args.db)
    create_folder(output)
    if args.tmp is None:
        tmp_dir = pathlib.Path(output, "tmp")
        create_folder(tmp_dir)
    else:
        tmp_dir = pathlib.Path(args.tmp, args.db, date_dir, "tmp")
        create_folder(tmp_dir)

    fastq = is_fastq(args.input)
    output_files = create_output_files(
        prefix=args.prefix, input=args.input, fastq=fastq, db=args.db
    )

    log_dir = pathlib.Path(output, "logs")
    create_folder(log_dir)
    extend_log_file = pathlib.Path(log_dir, output_files["extend_reads_log"])

    extended_reads = pathlib.Path(tmp_dir, output_files["extended_reads"])
    if extend:
        if not os.path.exists(extended_reads):
            # extend
            extend_reads(
                reads=args.input,
                output=output_files["extended_reads"],
                extend_bin=args.extend_bin,
                extend_k=args.extend_k,
                extend_length=args.extend_length,
                extend_memory=args.extend_memory,
                threads=args.threads,
                tmp_dir=tmp_dir,
                db=args.db,
                log_file=extend_log_file,
            )
        else:
            logging.info("Extended reads already exists. Skipping.")
    # dereplicate
    derep_log_file = pathlib.Path(log_dir, output_files["derep_reads_log"])
    derep_reads = pathlib.Path(tmp_dir, output_files["derep_reads"])

    if extend and os.path.exists(extended_reads) and derep:
        dereplicate_reads(
            reads=extended_reads,
            db=args.db,
            output=output_files["derep_reads"],
            derep_bin=args.derep_bin,
            tmp_dir=tmp_dir,
            log_file=derep_log_file,
            # derep_min_length=args.derep_min_length,
            threads=args.threads,
        )
    elif derep:
        dereplicate_reads(
            reads=args.input,
            db=args.db,
            output=output_files["derep_reads"],
            derep_bin=args.derep_bin,
            tmp_dir=tmp_dir,
            log_file=derep_log_file,
            # derep_min_length=args.derep_min_length,
            threads=args.threads,
        )
    # search
    if os.path.exists(derep_reads):
        reads = derep_reads
    elif os.path.exists(extended_reads) and not os.path.exists(derep_reads):
        reads = extended_reads
    else:
        reads = args.input

    results = pathlib.Path(output, output_files["results"])
    results_gz = pathlib.Path(output, output_files["results_gz"])
    search_log_file = pathlib.Path(log_dir, output_files["results_log"])
    if not os.path.exists(results_gz):
        search_mmseqs(
            reads=reads,
            db=args.db,
            vfdb_db=vfdb_db,
            output=results,
            tmp_dir=tmp_dir,
            log_file=search_log_file,
            mmseqs_bin=args.mmseqs_bin,
            mmseqs_min_length=args.mmseqs_min_length,
            mmseqs_evalue=args.mmseqs_evalue,
            mmseqs_min_seqid=args.mmseqs_min_seqid,
            mmseqs_cov=args.mmseqs_cov,
            mmseqs_cov_mode=args.mmseqs_cov_mode,
            ancient=args.ancient,
            threads=args.threads,
        )
    else:
        logging.info("Results already exists. Skipping.")

    if os.path.exists(results):
        with ExitStack() as stack:
            f_in = stack.enter_context(open(results, "rb"))
            f_out = stack.enter_context(gzip.open(results_gz, "wb"))
            shutil.copyfileobj(f_in, f_out)
        results.unlink()
    # filter
    extract_log_file = pathlib.Path(log_dir, output_files["reads_db_log"])
    reads_db = pathlib.Path(output, output_files["reads_db"])
    if fltr and os.path.exists(results_gz):
        filter_log_file = pathlib.Path(log_dir, output_files["results_filtered_log"])
        filter_results(
            results=output_files["results_gz"],
            log_file=filter_log_file,
            n_iters=args.iters,
            evalue=args.evalue,
            scale=args.scale,
            bitscore=args.bitscore,
            fltr=args.filter,
            breadth=args.breadth,
            breadth_expected_ratio=args.breadth_expected_ratio,
            depth=args.depth,
            depth_evenness=args.depth_evenness,
            trim=args.trim,
            db=args.db,
            threads=args.threads,
            xfilter_bin=args.xfilter_bin,
            prefix=output_files["results_filtered_prefix"],
            output=output,
            metadata=metadata,
            aggregate=args.agg,
        )
        cov_file = pathlib.Path(output, output_files["results_filtered_cov"]).resolve()
        mm_file = pathlib.Path(output, output_files["results_filtered_mm"]).resolve()
        # group_file = pathlib.Path(output, args.db, output_files["results_filtered_group"])
        # group_file_agg = pathlib.Path(output, args.db, output_files["results_filtered_group_agg"])

        if not os.path.exists(cov_file) or not os.path.exists(mm_file):
            logging.error("Cannot find filtered results. Exiting.")
            exit(1)
        if args.extract:
            if os.path.exists(mm_file):
                extract_reads(
                    results=mm_file,
                    db=args.db,
                    output=reads_db,
                    threads=args.threads,
                    reads=args.input,
                    extract_bin=args.extract_bin,
                    tmp_dir=tmp_dir,
                    log_file=extract_log_file,
                )
            else:
                if os.path.exists(results_gz):
                    extract_reads(
                        results=results_gz,
                        db=args.db,
                        output=reads_db,
                        threads=args.threads,
                        reads=args.input,
                        extract_bin=args.extract_bin,
                        tmp_dir=tmp_dir,
                        log_file=extract_log_file,
                    )

    # Clean up temporary
    if args.no_keep and os.path.exists(tmp_dir):
        logging.info("Removing temporary files.")
        shutil.rmtree(tmp_dir)
