"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import logging
from get_vf.createdb import createdb
from get_vf.utils import (
    get_arguments,
)
from get_vf.search import search_db

log = logging.getLogger("my_logger")


def main():

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )
    if args.action == "createdb":
        createdb(args)
    elif args.action == "search":
        search_db(args)


if __name__ == "__main__":
    main()


# import logging
# import pathlib
# import pandas as pd
# import os
# from get_vf.utils import (
#     get_arguments,
#     download,
#     create_folder_and_delete,
#     create_folder,
# )


# log = logging.getLogger("my_logger")


# def main():

#     logging.basicConfig(
#         level=logging.DEBUG, format="%(levelname)s ::: %(asctime)s ::: %(message)s"
#     )

#     args = get_arguments()

#     args = get_arguments()
#     logging.getLogger("my_logger").setLevel(
#         logging.DEBUG if args.debug else logging.INFO
#     )
#     if args.action == "split":
#         split_contigs(args)
#     elif args.action == "merge":
#         merge_contigs(args)

#     logging.getLogger("my_logger").setLevel(
#         logging.DEBUG if args.debug else logging.INFO
#     )
#     create_folder(args.outdir)
#     tmp_dir = f"{args.outdir}/tmp"
#     create_folder_and_delete(tmp_dir)

#     BASE_URL = f"http://files.metagenomics.eu//molecular_clock/DB/latest/"
#     # get list of marker genes
#     # check if metadata folder exists if not create it else delete it
#     metadata_dir = f"{args.outdir}/metadata"
#     metadata_file = pathlib.Path(metadata_dir, "bac120_marker_info.tsv")
#     # check if metadata file exists

#     if not os.path.exists(metadata_file):
#         create_folder(metadata_dir)
#         download(
#             f"{BASE_URL}/metadata/bac120_marker_info.tsv", dest_folder=metadata_dir
#         )
#     # get list of MSA genes
#     markers_df = pd.read_csv(metadata_file, sep="\t")
#     print(markers_df)


# if __name__ == "__main__":
#     main()
