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
from ete3 import Tree
import subprocess
import re
import io
from multiprocessing import Pool
from get_vf.defaults import gtdb_releases
import shutil

http.client.HTTPConnection.debuglevel = 0

import tqdm
from get_vf.utils import (
    get_arguments,
    create_folder_and_delete,
    create_folder,
    fast_flatten,
)


log = logging.getLogger("my_logger")

# Function to connect to the KEGG API and return the genomes in TSV format
def download(url: str, dest_folder: str):
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)  # create folder if it does not exist

    filename = os.path.basename(url)
    file_path = os.path.join(dest_folder, filename)
    chunk_size = 1024
    filesize = int(requests.head(url).headers["Content-Length"])
    with requests.get(url, stream=True) as r, open(file_path, "wb") as f, tqdm.tqdm(
        unit="B",  # unit string to be displayed.
        unit_scale=True,  # let tqdm to determine the scale in kilo, mega..etc.
        unit_divisor=1024,  # is used when unit_scale is true
        total=filesize,  # the total iteration.
        file=sys.stdout,  # default goes to stderr, this is the display on console.
        desc=filename,  # prefix to be displayed on progress bar.
        leave=False,
    ) as progress:
        if r.status_code == 200:
            for chunk in r.iter_content(chunk_size=chunk_size):
                # download the file chunk by chunk
                datasize = f.write(chunk)
                # on each chunk update the progress bar.
                progress.update(datasize)
        else:
            log.error(f"Error {r.status_code} when downloading {url}")
            os.remove(file_path)
            exit(1)


def retrieve_data(outdir, url, filetype, filename):
    file_dir = pathlib.Path(outdir, filetype)
    file_name = pathlib.Path(file_dir, filename)
    # check if metadata file exists
    if not os.path.exists(file_name):
        logging.info(f"Retrieving {filename}")
        create_folder(file_dir)
        download(f"{url}/{filetype}/{filename}", dest_folder=file_dir)
    else:
        logging.info(f"File {filename} already exists. Skipping download.")
    return file_name


def filter_faa(marker, marker_faa_file, marker_faa_file_filt, zscore):
    # Get seq statistics from the marker faa file
    logging.info(f"Calculating statistics for {marker} [faa]")
    seqs = {}
    for seq in SeqIO.parse(marker_faa_file, "fasta"):
        seqs[seq.id] = len(seq.seq)

    faa_lens = pd.DataFrame(seqs.items(), columns=["genome", "length"])
    faa_len_mean = faa_lens["length"].mean()
    faa_len_std = faa_lens["length"].std()

    logging.info(
        f"Fitering sequences with aa length z-score <= {zscore} [Mean: {faa_len_mean:.2f}, Std: {faa_len_std:.2f}]"
    )
    faa_ids = faa_lens[(np.abs(stats.zscore(faa_lens["length"])) <= zscore)][
        "genome"
    ].values
    faa_filt = (
        record
        for record in SeqIO.parse(marker_faa_file, "fasta")
        if record.id in faa_ids
    )

    logging.info(f"Writing aa filtered sequences to {marker_faa_file_filt}")
    SeqIO.write(faa_filt, marker_faa_file_filt, "fasta")
    return faa_ids


def prune_gtdb_tree(tree_file, marker, filt_tree_file, faa_ids):
    logging.info(f"Reading GTDB bac120 tree")
    t = Tree(str(tree_file), quoted_node_names=True, format=1)
    logging.info(f"Pruning tree for {marker} [leaves: {len(t.get_leaves()):,}]")
    t.prune(faa_ids)
    logging.info(f"Writing pruned tree for {marker} [leaves: {len(t.get_leaves()):,}]")
    t.write(outfile=str(filt_tree_file), quoted_node_names=True, format=1)


def align_job(faa, hmm, hmmalign_bin, marker):
    proc = subprocess.Popen(
        [
            hmmalign_bin,
            "--outformat",
            "PFAM",
            hmm,
            faa,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        rgx = re.compile("ERROR")
        for line in stdout.splitlines():
            if rgx.search(line):
                error = line.replace("ERROR: ", "")
        logging.error(f"Error aligning marker {marker} with hmmalign")
        logging.error(error)
        exit(1)

    alignment = AlignIO.read(io.StringIO(stdout), "stockholm")
    aln = []
    for record in alignment:
        record.seq = record.seq.upper()
        record.desc = ""
        aln.append(record)
    return aln


# From https://stackoverflow.com/a/20887930/15704171
def gapsFromPeptide(peptide_seq, nucleotide_seq):
    """Transfers gaps from aligned peptide seq into codon partitioned nucleotide seq (codon alignment)
    - peptide_seq is an aligned peptide sequence with gaps that need to be transferred to nucleotide seq
    - nucleotide_seq is an un-aligned dna sequence whose codons translate to peptide seq"""

    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i : i + n]

    codons = [
        codon for codon in chunks(nucleotide_seq, 3)
    ]  # splits nucleotides into codons (triplets)

    gappedCodons = []
    codonCount = 0
    for aa in peptide_seq:  # adds '---' gaps to nucleotide seq corresponding to peptide
        if aa != "-":
            gappedCodons.append(codons[codonCount])
            codonCount += 1
        else:
            gappedCodons.append("---")
    return "".join(gappedCodons)


def check_order(pro_align, nucl_seqs):
    nucl_id = set(nucl_seqs.keys())
    pro_id = {i.id for i in pro_align}
    # check if there is pro_id that does not have a nucleotide match
    if pro_id - nucl_id:
        diff = pro_id - nucl_id
        raise ValueError(
            f"Protein Record {', '.join(diff)} cannot find a "
            "nucleotide sequence match, please check the id"
        )
    else:
        pro_nucl_pair = []
        for pro_rec in pro_align:
            pro_nucl_pair.append((pro_rec, nucl_seqs[pro_rec.id]))
    return pro_nucl_pair


def aa2codon(aa, nt):
    # try:
    nt.seq = Seq.Seq(gapsFromPeptide(str(aa.seq), str(nt.seq)))
    # except IndexError as e:
    #     print(f"{e} -- {aa.id} -- {nt.id}")
    #     exit(1)
    nt.seq = nt.seq.upper()
    nt.description = ""
    return nt


def aa2codon_star(args):
    return aa2codon(*args)


def align_marker(
    faa, hmm, marker, aln_file_aa, aln_file_nt, aln_format, hmmalign_bin, fna, threads
):
    if not os.path.exists(aln_file_aa):
        logging.info(f"Aligning {marker} [faa] with hmmalign. Be patient...")
        aln = align_job(faa=faa, hmm=hmm, hmmalign_bin=hmmalign_bin, marker=marker)

        logging.info(f"Writing aa alignment for {marker} to {aln_file_aa}")
        output_handle = open(aln_file_aa, "w")
        AlignIO.write(AlignIO.MultipleSeqAlignment(aln), output_handle, aln_format)
        output_handle.close()
    else:
        logging.info(f"Aa alignment for {marker} already exists. Skipping alignment.")

    logging.info(f"Getting codon alignment for {marker} [faa]")
    alignment = AlignIO.read(aln_file_aa, aln_format)

    pro_nucl_pair = zip(list(alignment), list(fna))
    aln_nt = []
    p = Pool(threads)
    # for pair in pro_nucl_pair:
    #     print(f"{len(str(pair[0].seq).replace('-', ''))} -- {len(str(pair[1].seq))/3}")
    #     pair[1].seq = Seq.Seq(gapsFromPeptide(str(pair[0].seq), str(pair[1].seq)))
    #     pair[1].seq = pair[1].seq.upper()
    #     pair[1].description = ""
    #     aln_nt.append(pair[1])
    aln_nt = list(
        tqdm.tqdm(
            p.imap(aa2codon_star, pro_nucl_pair),
            total=len(list(alignment)),
            leave=False,
            desc=f"Alignments processed",
        )
    )

    logging.info(f"Writing nt alignment for {marker} to {aln_file_nt}")
    output_handle = open(aln_file_nt, "w")
    AlignIO.write(AlignIO.MultipleSeqAlignment(aln_nt), output_handle, aln_format)
    output_handle.close()


def create_mmseqs_db(faa, mmseqs_bin, mmseqs_db, marker):
    logging.info(f"Creating aa MMseqs2 DB for {marker}")
    proc = subprocess.Popen(
        [
            mmseqs_bin,
            "createdb",
            faa,
            mmseqs_db,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        # rgx = re.compile("ERROR")
        # for line in stdout.splitlines():
        #     if rgx.search(line):
        #         error = line.replace("ERROR: ", "")
        logging.error(f"Error creating mmseqs DB for {marker}")
        logging.error(stderr)
        exit(1)


def fetch_data(args):
    if args.local_db:
        logging.info("Fetching data from local DB")
        # Check that we have the files we need in the local db
        version_file = pathlib.Path(
            f"{args.local_db}/{gtdb_releases['latest']['version']}"
        )
        tree_file = pathlib.Path(
            f"{args.local_db}/tree/{gtdb_releases['latest']['tree']}"
        )
        metadata_file = pathlib.Path(
            f"{args.local_db}/metadata/{gtdb_releases['latest']['info_out']}"
        )
        faa_dir = pathlib.Path(f"{args.local_db}/faa")
        fna_dir = pathlib.Path(f"{args.local_db}/fna")
        hmm_dir = pathlib.Path(f"{args.local_db}/hmm")

        if (
            not os.path.exists(version_file)
            or not os.path.exists(tree_file)
            or not os.path.exists(metadata_file)
            or not os.path.exists(faa_dir)
            or not os.path.exists(fna_dir)
            or not os.path.exists(hmm_dir)
        ):
            logging.error(
                f"Files are missing. Please make sure you have the latest version of the findMarker database"
            )
            exit(1)

    logging.info(f"Creating output data {args.outdir}")
    create_folder(args.outdir)
    # tmp_dir = f"{args.outdir}/tmp"
    # create_folder_and_delete(tmp_dir)

    if args.local_db:
        # copy metadata file to outdir
        create_folder(f"{args.outdir}/metadata")
        shutil.copy(
            metadata_file,
            f"{args.outdir}/metadata/{gtdb_releases['latest']['info_out']}",
        )
        create_folder(f"{args.outdir}/tree")
        shutil.copy(tree_file, f"{args.outdir}/tree/{gtdb_releases['latest']['tree']}")

    else:
        BASE_URL = gtdb_releases["latest"]["fm_url"]
        # get list of marker genes
        metadata_file = retrieve_data(
            outdir=args.outdir,
            url=BASE_URL,
            filetype="metadata",
            filename=gtdb_releases["latest"]["info_out"],
        )
        # get GTDB tree
        tree_file = retrieve_data(
            outdir=args.outdir,
            url=BASE_URL,
            filetype="tree",
            filename=gtdb_releases["latest"]["tree"],
        )
    # get list of MSA genes
    markers_df = pd.read_csv(metadata_file, sep="\t")
    marker_data = markers_df[markers_df["marker"] == args.marker].to_dict(
        orient="records"
    )[0]
    if not marker_data:
        logging.error(f"Marker {args.marker} not found")
        sys.exit(1)
    marker_dir = pathlib.Path(args.outdir, "markers", args.marker, str(args.zscore))
    create_folder(marker_dir)

    if args.local_db:
        create_folder(f"{args.outdir}/faa")
        create_folder(f"{args.outdir}/fna")
        create_folder(f"{args.outdir}/hmm")

        marker_faa_file_src = pathlib.Path(
            f"{args.local_db}/faa/{marker_data['marker']}.faa"
        )
        marker_fna_file_src = pathlib.Path(
            f"{args.local_db}/fna/{marker_data['marker']}.fna"
        )
        marker_hmm_file_src = pathlib.Path(
            f"{args.local_db}/hmm/{marker_data['marker']}.hmm"
        )

        marker_faa_file = pathlib.Path(f"{args.outdir}/faa/{marker_data['marker']}.faa")
        marker_fna_file = pathlib.Path(f"{args.outdir}/fna/{marker_data['marker']}.fna")
        marker_hmm_file = pathlib.Path(f"{args.outdir}/hmm/{marker_data['marker']}.hmm")

        shutil.copy(marker_faa_file_src, marker_faa_file)
        shutil.copy(marker_fna_file_src, marker_fna_file)
        shutil.copy(marker_hmm_file_src, marker_hmm_file)

    else:
        marker_faa_file = retrieve_data(
            outdir=pathlib.Path(marker_dir),
            url=BASE_URL,
            filetype="faa",
            filename=f"{marker_data['marker']}.faa",
        )
        marker_fna_file = retrieve_data(
            outdir=pathlib.Path(marker_dir),
            url=BASE_URL,
            filetype="fna",
            filename=f"{marker_data['marker']}.fna",
        )
        marker_hmm_file = retrieve_data(
            outdir=pathlib.Path(marker_dir),
            url=BASE_URL,
            filetype="hmm",
            filename=f"{marker_data['marker']}.hmm",
        )
    marker_faa_file_filt = pathlib.Path(
        str(marker_faa_file).replace(".faa", ".filt.faa")
    )
    if not os.path.exists(marker_faa_file) or not os.path.exists(marker_faa_file_filt):
        faa_ids = filter_faa(
            marker=marker_data["marker"],
            marker_faa_file=marker_faa_file,
            marker_faa_file_filt=marker_faa_file_filt,
            zscore=args.zscore,
        )
    else:
        logging.info("Aa filtered file already exists. Skipping filter.")
        faa_ids = []
        for seq in SeqIO.parse(marker_faa_file_filt, "fasta"):
            faa_ids.append(seq.id)
    # prune tree
    tree_marker_dir = pathlib.Path(marker_dir, "tree")
    filt_tree_file = pathlib.Path(tree_marker_dir, f"{marker_data['marker']}.tree")
    if not os.path.exists(filt_tree_file):
        create_folder(tree_marker_dir)
        prune_gtdb_tree(
            tree_file=tree_file,
            marker=marker_data["marker"],
            filt_tree_file=filt_tree_file,
            faa_ids=faa_ids,
        )
    else:
        logging.info("Tree already pruned. Skipping pruning.")

    fna_filt = (
        record
        for record in SeqIO.parse(marker_fna_file, "fasta")
        if record.id in faa_ids
    )
    # align amino acids with hmmalign
    aln_dir = pathlib.Path(marker_dir, "aln")
    aln_file_faa = pathlib.Path(aln_dir, f"{marker_data['marker']}.filt.aa.aln")
    aln_file_fna = pathlib.Path(aln_dir, f"{marker_data['marker']}.filt.nt.aln")
    if not os.path.exists(aln_file_faa) or not os.path.exists(aln_file_fna):
        create_folder(aln_dir)
        align_marker(
            faa=marker_faa_file_filt,
            hmm=marker_hmm_file,
            marker=marker_data["marker"],
            aln_file_aa=aln_file_faa,
            aln_file_nt=aln_file_fna,
            aln_format=args.aln_format,
            hmmalign_bin=args.hmmalign_bin,
            fna=fna_filt,
            threads=args.threads,
        )
    else:
        logging.info("Sequences already aligned. Skipping alignment step.")

    # Create mmseqs2 db
    mmseqs_dir = pathlib.Path(marker_dir, "mmseqs")
    mmseqs_db = pathlib.Path(mmseqs_dir, f"{marker_data['marker']}-db")
    if not os.path.exists(mmseqs_db):
        create_folder(mmseqs_dir)
        create_mmseqs_db(
            faa=marker_faa_file_filt,
            marker=marker_data["marker"],
            mmseqs_bin=args.mmseqs_bin,
            mmseqs_db=mmseqs_db,
        )
    else:
        logging.info("MMseqs2 DB already built. Skipping creating DB.")
    logging.info(f"All files downloaded")
