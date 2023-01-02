from distutils.log import info
from operator import gt
from unittest import getTestCaseNames
import pandas as pd
import numpy as np
import requests
import io
import os, sys
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
from get_vf.defaults import gtdb_releases, info_file_columns
import logging


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
    faa_dir = f"{outdir}/faa"
    fna_dir = f"{outdir}/fna"
    hmm_dir = f"{outdir}/hmm"
    metadata_dir = f"{outdir}/metadata"
    tree_dir = f"{outdir}/tree"
    # check if output directory exists if not create it
    create_folder(outdir)
    # check if the subfolders are in place, if they exist remove
    # create_folder_and_delete(faa_dir)
    # create_folder_and_delete(fna_dir)
    # create_folder_and_delete(hmm_dir)
    if args.recreate:
        create_folder_and_delete(metadata_dir)
        create_folder_and_delete(tree_dir)
        delete_folder(faa_dir)
        delete_folder(fna_dir)
        delete_folder(hmm_dir)
    else:
        create_folder_and_delete(metadata_dir)
        create_folder_and_delete(tree_dir)
        # create_folder(faa_dir)
        # create_folder(fna_dir)

    # create a temporary directory to store the downloaded files
    if args.createdb_tmp is None:
        tmp_dir = f"{args.createdb_output}/tmp"
    else:
        tmp_dir = args.createdb_tmp
    create_folder(tmp_dir)

    # list files in tmp dir
    files = os.listdir(tmp_dir)

    # Define URLs
    base_url = gtdb_releases["latest"]["url"]
    version_url = f"{base_url}/VERSION"
    tree_url = f"{base_url}/{gtdb_releases['latest']['tree']}"
    marker_url = f"{base_url}/genomic_files_reps/{gtdb_releases['latest']['markers']}"
    gtdbtk_url = f"{base_url}/auxillary_files/{gtdb_releases['latest']['gtdbtk']}"
    info_url = f"{base_url}/auxillary_files/{gtdb_releases['latest']['info']}"
    # Get the version number of latest release
    log.info("Retrieving GTDB version information")
    version_file = pathlib.Path(f"{outdir}/version")
    r = get_data(version_url)
    rgx = re.compile("^v")
    with open(version_file, "w") as f:
        for row in r:
            if rgx.search(row):
                gtdb_version = row.replace("\n", "")
                gtdb_version = gtdb_version.replace("v", "")
            f.write(row)
    logging.info(f"Using GTDB version {gtdb_version}")

    # Download GTDB tree
    log.info(f"Retrieving GTDB v{gtdb_version} tree...")
    tree_file = pathlib.Path(f"{tree_dir}/{gtdb_releases['latest']['tree']}")

    if os.path.exists(tree_file):
        log.info(f"{tree_file} already exists. Skipping...")
    else:
        download(tree_url, tree_dir)

    # Download the marker genes from the GTDB website
    log.info(f"Downloading marker genes...")
    marker_file = pathlib.Path(f"{tmp_dir}/{gtdb_releases['latest']['markers']}")
    marker_dir = pathlib.Path(f"{tmp_dir}/bac120_marker_genes_reps_r{gtdb_version}")
    # Check if file exists if not download it
    if not os.path.exists(marker_file):
        download(marker_url, dest_folder=tmp_dir)
    else:
        log.info(
            f"File {gtdb_releases['latest']['markers']} already exists. Skipping..."
        )

    # if os.path.exists(fna_dir) and args.recreate:
    #     log.info(f"{fna_dir} already exists. Deleting...")
    #     shutil.rmtree(fna_dir)

    # if os.path.exists(faa_dir) and args.recreate:
    #     log.info(f"{faa_dir} already exists. Deleting...")
    #     shutil.rmtree(faa_dir)

    # check if {tmp_dir}/bac120_marker_genes_reps_r{gtdb_version} folder exists if not delete it

    if os.path.exists(marker_dir):
        log.info(f"{marker_file} already exists. Deleting...")
        shutil.rmtree(marker_dir)
    # Extract the marker genes and MSA files
    if os.path.exists(faa_dir) and os.path.exists(fna_dir):
        log.info(f"{faa_dir} already exists. Skipping...")
    else:
        log.info(f"Extracting {tmp_dir}/{gtdb_releases['latest']['markers']}")
        # delete folder
        if os.path.exists(faa_dir):
            shutil.rmtree(faa_dir)
        if os.path.exists(fna_dir):
            shutil.rmtree(fna_dir)
        with tarfile.open(marker_file) as tar:
            for member in tar.getmembers():
                if member.name.endswith(".faa"):
                    tar.extract(member, path=tmp_dir)
                elif member.name.endswith(".fna"):
                    tar.extract(member, path=tmp_dir)
                else:
                    continue
        shutil.move(
            f"{tmp_dir}/bac120_marker_genes_reps_r{gtdb_version}/fna", f"{outdir}"
        )
        shutil.move(
            f"{tmp_dir}/bac120_marker_genes_reps_r{gtdb_version}/faa", f"{outdir}"
        )

    # Download the HMM files from the GTDB website
    gtdbtk_file = pathlib.Path(f"{tmp_dir}/{gtdb_releases['latest']['gtdbtk']}")
    gtdbtk_dir = pathlib.Path(f"{tmp_dir}/release{gtdb_version}_v2")
    # Check if file exists if not download it
    if not os.path.exists(gtdbtk_file):
        log.info(f"Downloading {gtdb_releases['latest']['gtdbtk']}...")
        download(gtdbtk_url, dest_folder=tmp_dir)
    else:
        log.info(
            f"File {gtdb_releases['latest']['gtdbtk']} already exists. Skipping..."
        )

    if os.path.exists(hmm_dir) and is_folder_empty(hmm_dir):
        log.info(f"{hmm_dir} already exists but is empty. Deleting...")
        delete_folder(hmm_dir)

    if not os.path.exists(hmm_dir):
        create_folder(hmm_dir)
        if os.path.exists(gtdbtk_dir):
            log.info(f"{gtdbtk_dir} already exists. Deleting...")
            shutil.rmtree(gtdbtk_dir)
        # Extract the HMM files
        log.info(f"Extracting HMMs from {tmp_dir}/{gtdb_releases['latest']['gtdbtk']}")
        with tarfile.open(gtdbtk_file) as tar:
            for member in tar.getmembers():
                if member.name.endswith(".hmm") and "individual_hmms" in member.name:
                    tar.extract(member, path=tmp_dir)
                    shutil.move(
                        f"{tmp_dir}/{member.name}",
                        f"{hmm_dir}",
                    )
                elif member.name.endswith(".HMM") and "individual_hmms" in member.name:
                    tar.extract(member, path=tmp_dir)
                    fname = pathlib.Path(member.name).with_suffix(".hmm").name
                    fname = os.path.join(hmm_dir, fname)
                    shutil.move(
                        f"{tmp_dir}/{member.name}",
                        fname,
                    )
                else:
                    continue
    else:
        log.info(f"{hmm_dir} already exists. Skipping...")

    # Generate info file
    log.info(f"Downloading info file...")
    # Download the marker genes from the GTDB website
    info_file = pathlib.Path(f"{tmp_dir}/{gtdb_releases['latest']['info']}")
    info_file_out = pathlib.Path(
        f"{metadata_dir}/{gtdb_releases['latest']['info_out']}"
    )
    # Check if file exists if not download it
    if not os.path.exists(info_file_out):
        download(info_url, dest_folder=tmp_dir)
    else:
        log.info(
            f"File {gtdb_releases['latest']['markers']} already exists. Skipping..."
        )
    # read the info file
    info_df = pd.read_csv(info_file, sep="\t")
    info_df.columns = info_file_columns
    # split column marker into db and marker
    dbs = info_df["marker"].str.split("_", expand=True)[0]
    info_df["marker"] = info_df["marker"].str.split("_", expand=True)[1]
    info_df.insert(0, "db", dbs)
    log.info(f"Writing info file to {info_file_out}")
    info_df.to_csv(info_file_out, sep="\t", index=False)
