VFDB = {
    "core": {
        "url": "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz",
        "filename": "VFDB_setA_pro.fas.gz",
        "metadata": "VFDB.core.metadata.tsv.gz",
        "faa": "VFDB.core.fasta.gz",
    },
    "full": {
        "url": "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz",
        "filename": "VFDB_setB_pro.fas.gz",
        "metadata": "VFDB.full.metadata.tsv.gz",
        "faa": "VFDB.full.fasta.gz",
    },
}


vfdbs = ["core", "full"]

derep_bins = ["vsearch", "seqkit"]
