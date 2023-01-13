from operator import gt
from unittest import getTestCaseNames
import pandas as pd
import numpy as np
import requests
import io
import os
import sys
import tqdm
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from multiprocessing import Pool
from functools import partial
import re
import time
import tarfile
import shutil
import pathlib
from get_vf.defaults import VFDB, vfdbs
import logging
from datetime import date
import gzip
from Bio import SeqIO
import subprocess


def create_mmseqs_db(db_dir, db, faa, mmseqs_bin):
    # Create mmseqs2 db
    mmseqs_dir = pathlib.Path(db_dir, "mmseqs")
    mmseqs_db = pathlib.Path(mmseqs_dir, f"{db}-db")
    if not os.path.exists(mmseqs_db):
        create_folder(mmseqs_dir)
        create_mmseqs_db_worker(
            faa=faa,
            mmseqs_bin=mmseqs_bin,
            mmseqs_db=mmseqs_db,
            db=db,
        )
    else:
        logging.info("MMseqs2 DB already built. Skipping creating DB.")


def create_mmseqs_db_worker(faa, mmseqs_bin, mmseqs_db, db):
    logging.info(f"Creating aa MMseqs2 DB for {db} VFDB")
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
        logging.error(f"Error creating mmseqs DB for {db}")
        logging.error(stderr)
        exit(1)


def parse_header(header):
    # Use a regular expression to extract the gene name, description, virulence factor, and organism
    match = re.match(r"(\S+)\s(\S.+)\[(\S.+)\]\s\[(\S.+)\]", header)

    # If the regular expression matched, return the extracted information
    if match:
        return (match.group(1), match.group(2), match.group(3), match.group(4))
    # If the regular expression did not match, return None
    else:
        return None


# def create_output_files(outdir):
#     # create output files
#     out_files = {
#         "core": {
#             "out_tsv": f"{outdir}/core/VFDB.metadata.tsv.gz",
#             "out_fasta": f"{outdir}/core/VFDB.fasta.gz",
#         },
#         "full": {
#             "out_tsv": f"{outdir}/full/VFDB.metadata.tsv.gz",
#             "out_fasta": f"{outdir}/full/VFDB.full.fasta.gz",
#         },
#     }
#     return out_files


def create_db_files(input, faa, metadata):
    with open(input, "rb") as handle:
        # Read the first two bytes of the file
        magic = handle.read(2)

        # If the first two bytes are the GZIP magic number, the file is compressed
        if magic == b"\x1f\x8b":
            mode = "rt"
        else:
            mode = "r"

        # Seek back to the beginning of the file
        handle.seek(0)

        # Open the file with gzip.open, using the determined mode
        with gzip.open(handle, mode) as handle:
            # Open the output file
            with gzip.open(metadata, "wt") as output_tsv, gzip.open(
                faa, "wt"
            ) as fasta_output:
                # Loop through the sequences in the FASTA file
                for record in SeqIO.parse(handle, "fasta"):
                    # Parse the header and extract the information
                    info = parse_header(record.description)
                    # If the header was successfully parsed, write the extracted information to the output file
                    if info:
                        string = f"{info[0]}\t{info[2]}\n"
                        record.id = info[0]
                        record.description = ""
                        output_tsv.write(string)
                        SeqIO.write(record, handle=fasta_output, format="fasta")


def get_data(url):
    # Get the data from the API
    retry_strategy = Retry(
        total=10,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504, 403],
        allowed_methods=["HEAD", "GET", "OPTIONS"],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    http.mount("http://", adapter)
    r = http.get(url, timeout=10)
    # Parse the data
    if r.status_code != 200:
        r = io.StringIO(f"Error: {r.status_code}")
    else:
        r = io.StringIO(r.content.decode("utf-8"))
    return r


def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        log.info(f"Folder {path} already exists. Skipping...")


def create_folder_and_delete(path):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        log.info(f"Folder {path} already exists. Deleting...")
        try:
            shutil.rmtree(path)
            os.makedirs(path)
        except OSError as e:
            log.error(f"{path} : {e.strerror}")


def is_folder_empty(path):
    with os.scandir(path) as it:
        if any(it):
            return False
        else:
            return True


def delete_folder(path):
    if os.path.exists(path):
        log.info(f"Deleting folder {path}...")
        try:
            shutil.rmtree(path)
        except OSError as e:
            log.error(f"{path} : {e.strerror}")


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
    ) as progress:
        for chunk in r.iter_content(chunk_size=chunk_size):
            # download the file chunk by chunk
            datasize = f.write(chunk)
            # on each chunk update the progress bar.
            progress.update(datasize)


log = logging.getLogger("my_logger")


def createdb(args):
    # check if output directory exists if not create it else delete it
    outdir = f"{args.createdb_output}/DB"
    core_dir = f"{outdir}/core"
    full_dir = f"{outdir}/full"
    # out_files = create_output_files(outdir)

    # check if output directory exists if not create it
    create_folder(outdir)
    # check if the subfolders are in place, if they exist remove
    # create_folder_and_delete(faa_dir)
    # create_folder_and_delete(fna_dir)
    # create_folder_and_delete(hmm_dir)
    if args.recreate:
        delete_folder(core_dir)
        delete_folder(full_dir)
        create_folder(core_dir)
        create_folder(full_dir)
    else:
        create_folder(core_dir)
        create_folder(full_dir)

    # create a temporary directory to store the downloaded files
    if args.createdb_tmp is None:
        tmp_dir = f"{args.createdb_output}/tmp"
    else:
        tmp_dir = args.createdb_tmp
    create_folder(tmp_dir)

    # Get the date of the latest release
    today_date = date.today()
    version_file = pathlib.Path(f"{outdir}/version")
    with open(version_file, "w") as f:
        f.write(f"{today_date}")
    # Download the files to the tmp directory
    for db in vfdbs:
        log.info(f"Downloading VFDB {db} database")
        download(VFDB[db]["url"], tmp_dir)

    # Create metadata file
    log.info("Processing files")
    for db in vfdbs:
        create_db_files(
            input=f"{tmp_dir}/{VFDB[db]['filename']}",
            faa=pathlib.Path(outdir, db, VFDB[db]["faa"]),
            metadata=pathlib.Path(outdir, db, VFDB[db]["metadata"]),
        )
        create_mmseqs_db(
            db_dir=core_dir,
            db=db,
            faa=pathlib.Path(outdir, db, VFDB[db]["faa"]),
            mmseqs_bin=args.mmseqs_bin,
        )

    # Delete the temporary directory
    log.info("Removing temporary directory")
    delete_folder(tmp_dir)

    log.info("Done")
