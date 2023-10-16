
# TRANSCRIPT ID (ENSEMBL, REFSEQ) & GENE NAME CONVERSION

from typing import List, Dict, Any
from pybiomart import Dataset as BiomartDataset
import os
import json


class NameLookup:

    def __init__(self, name_lookup_dict: Dict[str, Dict[str, Dict[str, str]]]):

        self._dict = name_lookup_dict

    def get_all_gene_names(self) -> List[str]:

        return [gene_name for gene_name in self._dict["ensembl"]["gene"].values() if gene_name]

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
            if to_format == "refseq" or "gene":
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
    verbose: bool = False
) -> NameLookup:

    def print_if_verbose(s: Any):
        if verbose:
            print(s)

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

    print_if_verbose("Obtaining info from Biomart...")

    id_df = biomart_dataset.query(
        attributes=query_attributes,
        filters=query_filters
    )

    print_if_verbose("...done\n")
    print_if_verbose("Creating lookup object from Biomart data")

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

        if verbose:
            if _i % 1000 == 0:
                print(f"Progress: ({_i}/{_t})")
            _i += 1

        name_lookup_dict["refseq"]["ensembl"][rec[1]["RefSeq mRNA ID"]] = str(
            rec[1]["Transcript stable ID version"]
        ).replace("nan", "")
        name_lookup_dict["refseq"]["gene"][rec[1]["RefSeq mRNA ID"]] = str(
            rec[1]["Gene name"]
        ).replace("nan", "")
        name_lookup_dict["ensembl"]["refseq"][rec[1]["Transcript stable ID version"]] = str(
            rec[1]["RefSeq mRNA ID"]
        ).replace("nan", "")
        name_lookup_dict["ensembl"]["gene"][rec[1]["Transcript stable ID version"]] = str(
            rec[1]["Gene name"]
        ).replace("nan", "")

    print_if_verbose("...done\n")
    print_if_verbose("Saving to file...")

    with open(os.path.join(output_path), "w") as f:
        json.dump(name_lookup_dict, f)

    print_if_verbose("...done\n")

    return NameLookup(name_lookup_dict)


def load_name_lookup_from_file(path: str) -> NameLookup:

    with open(path, "r") as f:

        name_lookup_dict = json.load(f)

    return NameLookup(name_lookup_dict)
