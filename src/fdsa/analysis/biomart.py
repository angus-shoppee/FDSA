# Feature-Directed Splice Analysis
# Copyright (C) 2024 Angus Shoppee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# TRANSCRIPT ID (ENSEMBL, REFSEQ) & GENE NAME CONVERSION

from typing import List, Dict
import logging
from pybiomart import Dataset as BiomartDataset
import os
import json


logger = logging.getLogger(__name__)


# TODO: Replace biomart with NCBI's "gene2ensembl" for name lookup
#       https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
# TODO: Add option in build config to specify path to gene2ensembl file and avoid downloading
#       (if auto-downloaded, should be stored in root of build folder for use in future build steps)
# TODO: Remove redundant conversions from name lookup - only need ensembl -> refseq & ensembl -> geneID
# TODO: (redundant once implementation upgraded?) fix values stored against "None" as key
# Possible TODO: For human, allow use of the MANE_Select tag to filter transcripts?


class NameLookup:

    def __init__(
        self,
        name_lookup_dict: Dict[str, Dict[str, Dict[str, str]]]
    ):

        self._dict = name_lookup_dict

    def get_all_gene_names(self) -> List[str]:

        # TODO: Refactor to get_all_gene_ids

        return [gene_name for gene_name in set(self._dict["ensembl"]["gene"].values()) if gene_name]

    def convert(
        self,
        value: str,
        from_format: str,
        to_format: str,
    ) -> str:

        if from_format == "refseq":
            if to_format == "ensembl" or to_format == "gene":
                return self._dict[from_format][to_format][value]
            elif to_format == "refseq":
                return value
            else:
                raise ValueError(f"to_format must be either refseq, ensembl or gene. Got: '{to_format}'")

        elif from_format == "ensembl":
            if to_format == "refseq" or to_format == "gene":
                return self._dict[from_format][to_format][value]
            elif to_format == "ensembl":
                return value
            else:
                raise ValueError(f"to_format must be either refseq, ensembl or gene. Got: '{to_format}'")

        else:
            raise ValueError(f"from_format must be either refseq or ensembl. Got: {from_format}")


def create_and_save_name_lookup(
    database_name: str,
    host_url: str,
    output_path: str,
    use_first_id_match = True
) -> NameLookup:

    biomart_dataset = BiomartDataset(
        name=database_name,
        host=host_url
    )

    query_attributes = [
        "ensembl_transcript_id_version",
        "refseq_mrna",
        "external_gene_name"
    ]

    query_filters = {}

    logger.info("Obtaining info from Biomart...")

    id_df = biomart_dataset.query(
        attributes=query_attributes,
        filters=query_filters
    )

    logger.info("...done\n")
    logger.info("Creating lookup object from Biomart data")

    name_lookup_dict = {
        "refseq": {
            "ensembl": {},
            "gene": {},
        },
        "ensembl": {
            "refseq": {},
            "gene": {}
        }
    }

    _i = 0
    _t = len(id_df)

    for rec in id_df.iterrows():

        if _i % 1000 == 0:
            logger.info(f"Progress: ({_i}/{_t})")
        _i += 1

        if use_first_id_match:
            name_lookup_dict["refseq"]["ensembl"].setdefault(
                rec[1]["RefSeq mRNA ID"],
                str(rec[1]["Transcript stable ID version"]).replace("nan", "")
            )
            name_lookup_dict["refseq"]["gene"].setdefault(
                rec[1]["RefSeq mRNA ID"],
                str(rec[1]["Gene name"]).replace("nan", "")
            )
            name_lookup_dict["ensembl"]["refseq"].setdefault(
                rec[1]["Transcript stable ID version"],
                str(rec[1]["RefSeq mRNA ID"]).replace("nan", "")
            )
            name_lookup_dict["ensembl"]["gene"].setdefault(
                rec[1]["Transcript stable ID version"],
                str(rec[1]["Gene name"]).replace("nan", "")
            )
        else:
            # name_lookup_dict["refseq"]["ensembl"][rec[1]["RefSeq mRNA ID"]] = str(
            #     rec[1]["Transcript stable ID version"]
            # ).replace("nan", "")
            # name_lookup_dict["refseq"]["gene"][rec[1]["RefSeq mRNA ID"]] = str(
            #     rec[1]["Gene name"]
            # ).replace("nan", "")
            # name_lookup_dict["ensembl"]["refseq"][rec[1]["Transcript stable ID version"]] = str(
            #     rec[1]["RefSeq mRNA ID"]
            # ).replace("nan", "")
            # name_lookup_dict["ensembl"]["gene"][rec[1]["Transcript stable ID version"]] = str(
            #     rec[1]["Gene name"]
            # ).replace("nan", "")
            raise NotImplementedError("Only first-match is currently supported for multiple ID mappings")

    logger.info("...done\n")
    logger.info("Saving to file...")

    with open(os.path.join(output_path), "w") as f:
        json.dump(name_lookup_dict, f)

    logger.info("...done\n")

    return NameLookup(name_lookup_dict)


def load_name_lookup_from_file(path: str) -> NameLookup:

    with open(path, "r") as f:

        name_lookup_dict = json.load(f)

    return NameLookup(name_lookup_dict)
