
# getVF: recruit reads mapping to Virulence Factor databases


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/get-ancient-vf?include_prereleases&label=version)](https://github.com/genomewalker/get-ancient-vf/releases) [![get-ancient-vf](https://github.com/genomewalker/get-ancient-vf/workflows/getVF_ci/badge.svg)](https://github.com/genomewalker/get-ancient-vf/actions) [![PyPI](https://img.shields.io/pypi/v/get-ancient-vf)](https://pypi.org/project/get-ancient-vf/) [![Conda](https://img.shields.io/conda/v/genomewalker/get-ancient-vf)](https://anaconda.org/genomewalker/get-ancient-vf)

A simple tool to find reads that might map to virulence factors. 

# Installation

We recommend having [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First, we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/get-ancient-vf/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install get-ancient-vf
```

### Using conda

```bash
conda install -c conda-forge -c bioconda -c genomewalker get-ancient-vf
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/genomewalker/get-ancient-vf.git
```

By cloning in a dedicated conda environment

```bash
git clone https://github.com/genomewalker/get-ancient-vf.git
cd get-ancient-vf
conda env create -f environment.yml
conda activate get-ancient-vf
pip install -e .
```


# Usage

getVF has two subcommands `createdb` and `search`. 


## Create the database structure needed for getVF

The subcommand `createdb` creates the database structure needed for getVF. It will download all the necessary files from [VFDB](http://www.mgc.ac.cn/VFs/download.htm) and create the basic data for `search` to work. For a complete list of options:

```bash
$ getVF createdb --help
usage: getVF createdb [-h] [--threads THREADS] --output STR [--tmp STR] [--recreate] [--mmseqs-bin STR]

optional arguments:
  -h, --help         show this help message and exit
  --threads THREADS  Number of threads to use (default: 1)

Createdb optional arguments:
  --output STR       Output folder where to store the DB (default: None)
  --tmp STR          Temporary directory location (default: None)
  --recreate         Remove folders if they exist
  --mmseqs-bin STR   Path to the mmseqs2 executable (default: mmseqs)
```


One would run the `createdb` subcommand as:

```bash
getVF createdb --output VFDB --recreate
```

This command will retrieve the data from VFDB and will process it. The output will be stored in the `VFDB` folder. It will get the amino acid sequences and the metadata for the `core` and `full` DBs. The generated folder can be used with the `search` subcommand.

It will generate the following files:

```
VFDB
└── DB
    ├── core
    │   ├── VFDB.core.fasta.gz
    │   ├── VFDB.core.metadata.tsv.gz
    │   └── mmseqs
    │       ├── core-db
    │       ├── core-db.dbtype
    │       ├── core-db.index
    │       ├── core-db.lookup
    │       ├── core-db.source
    │       ├── core-db_h
    │       ├── core-db_h.dbtype
    │       └── core-db_h.index
    ├── full
    │   ├── VFDB.full.fasta.gz
    │   ├── VFDB.full.metadata.tsv.gz
    │   └── mmseqs
    │       ├── full-db
    │       ├── full-db.dbtype
    │       ├── full-db.index
    │       ├── full-db.lookup
    │       ├── full-db.source
    │       ├── full-db_h
    │       ├── full-db_h.dbtype
    │       └── full-db_h.index
    └── version
```

## Search reads against the marker DB

With the subcommand search one can search the reads against the marker DB. Here we can decide between different approaches. By default, the short reads will be extended by a gentle assembly on both ends, then the extended reads will be de-replicated at 100% identity and length using [seqkit]() or [vsearch](); and finally, the reads will be mapped against the VFDB using MMseqs2. The results from the BLASTx search will be filtered using [x-filter](https://github.com/genomewalker/x-filter) to identify the genes with the highest likelihood of being in the sample. Once all steps are done, a fastQ or a fastA file will be generated with the reads that map against the VFDB. 


For a complete list of options:

```
usage: getVF search [-h] [--threads THREADS] --input INPUT --vfdb STR --db-dir STR [--output STR] [--tmp STR]
                    [--prefix STR] [--no-extend] [--extend-bin STR] [--extend-length INT] [--extend-k INT]
                    [--extend-memory INT] [--no-derep] [--derep-bin STR] [--mmseqs2-bin STR]
                    [--mmseq2-min-length INT] [--mmseqs2-evalue FLOAT] [--mmseq2-min-seqid FLOAT]
                    [--mmseq2-cov FLOAT] [--mmseq2-cov-mode INT] [--ancient] [--keep] [--no-filter]
                    [--x-filter-bin STR] [--n-iters INT] [--evalue FLOAT] [--scale FLOAT] [--bitscore INT]
                    [--filter STR] [--breadth FLOAT] [--breadth-expected-ratio FLOAT] [--depth FLOAT]
                    [--depth-evenness FLOAT] [--no-trim] [--no-aggregate] [--extract] [--extract-bin STR]

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads to use (default: 1)

Search required arguments:
  --input INPUT         A FASTA file containing the query sequences
  --vfdb STR            Which VFDB to use (default: core)
  --db-dir STR          Folder with the VFDB data generated by the createdb subcommand

Search optional arguments:
  --output STR          Output folder to write the search results (default: get-ancient-vf-output-search)
  --tmp STR             Path to the temporary directory (default: None)
  --prefix STR          Prefix used for the output files (default: None)
  --no-extend           Disable read extension
  --no-derep            Disable de-replication of identical reads
  --ancient             Use mmseqs2 aDNA optimized parameters
  --keep                Keep temporary files
  --no-filter           Disable x-filter filtering

Read extension arguments:
  --extend-bin STR      Path to the the executable for the extension step (default: tadpole.sh)
  --extend-length INT   How much to extend in both ends (default: 100)
  --extend-k INT        K-mer length to use for the extension (default: 17)
  --extend-memory INT   How much memory to use for the extension (default: None)

Read dereplication arguments:
  --derep-bin STR       Path to the the executable for the de-replication step (default: seqkit)

MMseqs2 search arguments:
  --mmseqs2-bin STR     Path to the mmseqs2 executable (default: mmseqs)
  --mmseq2-min-length INT
                        Minimum length of the predicted ORF (default: 15)
  --mmseqs2-evalue FLOAT
                        Evalue where to filter the results (default: 1e-05)
  --mmseq2-min-seqid FLOAT
                        List matches above this sequence identity (default: 0.6)
  --mmseq2-cov FLOAT    List matches above this fraction of aligned (covered) residues (default: 0.6)
  --mmseq2-cov-mode INT
                        Type of mmseq2 coverage calculation (default: 2)

X-filter processing arguments:
  --x-filter-bin STR    Path to the the executable for the x-filter step (default: xFilter)
  --n-iters INT         Number of iterations for the FAMLI-like filtering (default: 25)
  --evalue FLOAT        Evalue where to filter the results (default: 1e-05)
  --scale FLOAT         Scale to select the best weithing alignments (default: 0.9)
  --bitscore INT        Bitscore where to filter the results (default: 60)
  --filter STR          Which filter to use. Possible values are: breadth, depth, depth_evenness,
                        breadth_expected_ratio (default: depth_evenness)
  --breadth FLOAT       Breadth of the coverage (default: 0.5)
  --breadth-expected-ratio FLOAT
                        Expected breath to observed breadth ratio (scaled) (default: 0.5)
  --depth FLOAT         Depth to filter out (default: 0.1)
  --depth-evenness FLOAT
                        Reference with higher evenness will be removed (default: 1.0)
  --no-trim             Deactivate the trimming for the coverage calculations
  --no-aggregate        Aggregate the results by virulence factor

Extract VFDB mapping reads arguments:
  --extract             Extract the reads that match the VFDB
  --extract-bin STR     Path to the the executable for the read extraction step (default: filterbyname.sh)
```

One would run the following command to search the reads against the VFDB:

```bash
getVF search --input test.fq.gz --vfdb core --db-dir VFDB/DB/ --threads 8 --ancient --filter depth_evenness --depth-evenness 1.0
```

The output folder will contain the following files:

```
get-ancient-vf-output-search/
└── core
    ├── logs
    │   ├── test.derep.log
    │   ├── test.extend.log
    │   ├── test.reads.core.log
    │   └── test.results.core.log
    ├── test.reads.core.fq.gz
    ├── test.results.core.tsv.gz
    ├── test.results-filtered.core_cov-stats.tsv.gz
    ├── test.results-filtered.core_group-abundances-agg.tsv.gz
    ├── test.results-filtered.core_group-abundances.tsv.gz
    └── test.results-filtered.core_no-multimap.tsv.gz
```