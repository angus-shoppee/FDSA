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


from typing import List, Tuple, Dict, Any, Union, Optional
import csv
import pandas as pd

from config import output_column_names as cols
from analysis.experiment import Sample
from reporting.parse_stringtie_combined import calculate_fraction_lacking_feature


FDSA_RESULT_FREQUENCY_COLUMN_PREFIX = "Percent "


class QuantifiedSpliceJunctionLocus:
    start: int
    end: int
    n: int

    def __init__(
        self,
        locus: str,
        n: int
    ) -> None:

        locus_components = locus.split("-")

        self.start = int(locus_components[0])
        self.end = int(locus_components[1])
        self.n = n

    def __eq__(self, other: Any) -> bool:
        """start and end, but not necessarily n, must be equal"""
        if isinstance(other, QuantifiedSpliceJunctionLocus):
            if self.start == other.start and self.end == other.end:
                return True
        return False


class FdsaResult:
    transcript_id: str
    gene_id: str
    gene_name: str
    chromosome: str
    strand: str
    exon_positions: List[Tuple[int, int]]
    feature_qualifiers: str
    feature_region: Tuple[int, int]
    feature_number: int
    total_features_in_transcript: int
    all: Dict[str, List[QuantifiedSpliceJunctionLocus]]
    overlapping: Dict[str, List[QuantifiedSpliceJunctionLocus]]
    frequencies: Dict[str, float]
    # frequencies: List[Dict[str, float]]
    stringtie_frequencies: Optional[Dict[str, float]]  # TODO: Give frequencies and stringtie_frequencies better names

    def __init__(
        self,
        transcript_id: str,
        gene_id: str,
        gene_name: str,
        chromosome: str,
        strand: str,
        exon_positions_string: str,
        feature_qualifiers: str,
        feature_region_string: str,
        feature_number: int,
        total_features_in_transcript: int,
        junction_definition: str,
        all_junction_dump_string: str,
        overlapping_junction_dump_string: str,
        # TODO: Confirm intended type is Dict[str, float] and not List[Dict[str, float]]
        frequencies: Dict[str, float],
        # frequencies: List[Dict[str, float]]
        stringtie_frequencies: Optional[Dict[str, float]]

    ) -> None:
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chromosome = chromosome
        self.strand = strand

        _positions = [position.split("-") for position in exon_positions_string.split(" ")]
        self.exon_positions = [(int(position[0]), int(position[1])) for position in _positions if len(position) == 2]

        self.feature_qualifiers = feature_qualifiers

        _feature_region = feature_region_string.split("-")
        self.feature_region = (
            int(_feature_region[0]), int(_feature_region[1])
        ) if len(_feature_region) == 2 else (
            int(_feature_region[0]), int(_feature_region[0])
        )

        self.feature_number = feature_number
        self.total_features_in_transcript = total_features_in_transcript

        junction_locus_by_id = parse_junction_locus_by_id_from_string(junction_definition)

        self.all = parse_junctions_from_string(all_junction_dump_string, junction_locus_by_id)
        self.overlapping = parse_junctions_from_string(overlapping_junction_dump_string, junction_locus_by_id)

        self.frequencies = frequencies

        self.stringtie_frequencies = stringtie_frequencies

    # Unused
    def total_number_occurrences_across_all_samples(self) -> int:

        return sum([sum([locus.n for locus in loci] for loci in self.overlapping.values())])

    def get_exon_positions_string(self) -> str:

        return " ".join([f"{e[0]}-{e[1]}" for e in self.exon_positions])

    def get_feature_region_string(self) -> str:

        return f"{self.feature_region[0]}-{self.feature_region[1]}"

    @classmethod
    def get_serialized_header_for_info_cols(cls) -> List[str]:

        return [
            cols.GENE_ID, cols.GENE_NAME, cols.CHROMOSOME, cols.STRAND,
            cols.EXON_POSITIONS, cols.FEATURE_QUALIFIERS, cols.FEATURE_REGION,
            cols.N_FEATURES_IN_TRANSCRIPT, cols.FEATURE_NUMBER
        ]

    def serialize_info_cols(self) -> List[str]:

        return [
            self.gene_id, self.gene_name, self.chromosome, self.strand,
            self.get_exon_positions_string(), self.feature_qualifiers, self.get_feature_region_string(),
            str(self.total_features_in_transcript), str(self.feature_number)
        ]


def load_fdsa_results_as_df(
    fdsa_results_path: str,
    samples: Dict[str, Sample],
    rank_by: str = "frequency",   # TODO: Parse default rank_by from config (add stringtie_frequency option?)
    force_gene_names: Union[None, List[str]] = None,
    min_total_number_occurrences_across_all_samples: int = 1,
    min_per_sample_occurrences_number_occurrences: int = 0,
    min_per_sample_occurrences_in_at_least_n_samples: int = 0,
    fdsa_result_frequency_column_prefix: str = FDSA_RESULT_FREQUENCY_COLUMN_PREFIX
) -> pd.DataFrame:

    if force_gene_names is None:
        force_gene_names = list()

    with open(fdsa_results_path, "r") as f:

        results = csv.reader(f)

        header = next(results)

        raw_number_column_names = [sample_name for sample_name in samples.keys()]
        frequency_column_names = [fdsa_result_frequency_column_prefix + sample_name for sample_name in samples.keys()]

        gene_name_column_index: Union[None, int] = None
        info_column_indexes, raw_number_column_indexes, frequency_column_indexes = [], [], []
        for i, value in enumerate(header):

            if value == cols.GENE_NAME:
                gene_name_column_index = i

            if value in raw_number_column_names:
                raw_number_column_indexes.append(i)
            elif value in frequency_column_names:
                frequency_column_indexes.append(i)
            else:
                info_column_indexes.append(i)

        if gene_name_column_index is None:
            raise ValueError(f"{cols.GENE_NAME} column was not found in header row")

        if not all([raw_number_column_indexes, frequency_column_indexes]):
            raise ValueError(
                "Unable to find matching sample data in FDSA results. Check that BAM files match FDSA output.\n" +
                f"Sample names from BAM files: {samples.keys()}\n" +
                f"Column names in FDSA output: {header}"
            )

        if force_gene_names:

            rows_to_import = [row for row in results if row[gene_name_column_index].lower() in force_gene_names]

        else:

            rows_to_import = [row for row in results if all([
                sum(
                    [int(n) for n in row[raw_number_column_indexes[0]:frequency_column_indexes[0]]]
                ) >= min_total_number_occurrences_across_all_samples,
                len(
                    [True for n in row[raw_number_column_indexes[0]:frequency_column_indexes[0]] if
                     int(n) >= min_per_sample_occurrences_number_occurrences]
                ) >= min_per_sample_occurrences_in_at_least_n_samples
            ])]

        fdsa_results_df = pd.DataFrame(
            rows_to_import,
            columns=header
        )

        # fdsa_results_df = pd.DataFrame(
        #     [row for row in results if row[gene_name_column_index] in force_gene_names or all([
        #         sum(
        #             [int(n) for n in row[raw_number_column_indexes[0]:frequency_column_indexes[0]]]
        #         ) >= min_total_number_occurrences_across_all_samples,
        #         len(
        #             [True for n in row[raw_number_column_indexes[0]:frequency_column_indexes[0]] if
        #              int(n) >= min_per_sample_occurrences_number_occurrences]
        #         ) >= min_per_sample_occurrences_in_at_least_n_samples
        #     ])],
        #     columns=header
        # )

        if len(fdsa_results_df) == 0:
            raise ValueError("No results were loaded - check that specified criteria are not too stringent")

        # Convert integer columns
        fdsa_results_df[cols.TRANSCRIPT_START] = fdsa_results_df[cols.TRANSCRIPT_START].astype(int)
        fdsa_results_df[cols.N_FEATURES_IN_TRANSCRIPT] = fdsa_results_df[cols.N_FEATURES_IN_TRANSCRIPT].astype(int)
        fdsa_results_df[cols.FEATURE_NUMBER] = fdsa_results_df[cols.FEATURE_NUMBER].astype(int)
        for column_name in raw_number_column_names:
            fdsa_results_df[column_name] = fdsa_results_df[column_name].astype(int)

        # Convert float columns
        fdsa_results_df[cols.AVG_NUMBER] = fdsa_results_df[cols.AVG_NUMBER].astype(float)
        fdsa_results_df[cols.AVG_FREQUENCY] = fdsa_results_df[cols.AVG_FREQUENCY].astype(float)
        for column_name in frequency_column_names:
            fdsa_results_df[column_name] = fdsa_results_df[column_name].astype(float)

        # Order results by either number or avg frequency, depending on report config
        if rank_by == "frequency":
            fdsa_results_df = fdsa_results_df.sort_values(
                by=[cols.AVG_FREQUENCY, cols.AVG_NUMBER],
                ascending=[False, False]
            )
        elif rank_by == "number":
            fdsa_results_df = fdsa_results_df.sort_values(
                by=[cols.AVG_NUMBER, cols.AVG_FREQUENCY],
                ascending=[False, False]
            )
        else:
            raise ValueError(f"Invalid value for rank_by specified (\"{rank_by}\"). Options are frequency/number")

        return fdsa_results_df


def convert_fdsa_results_df_to_objects(
    fdsa_results_df: pd.DataFrame,
    stringtie_results_df: Optional[pd.DataFrame]
) -> List[FdsaResult]:

    zipped_info_columns = list(zip(
        fdsa_results_df[cols.TRANSCRIPT_ID],
        fdsa_results_df[cols.GENE_ID],
        fdsa_results_df[cols.GENE_NAME],
        fdsa_results_df[cols.CHROMOSOME],
        fdsa_results_df[cols.STRAND],
        fdsa_results_df[cols.EXON_POSITIONS],
        fdsa_results_df[cols.FEATURE_QUALIFIERS],
        fdsa_results_df[cols.FEATURE_REGION],
        fdsa_results_df[cols.FEATURE_NUMBER],
        fdsa_results_df[cols.N_FEATURES_IN_TRANSCRIPT],
        fdsa_results_df[cols.JUNCTION_DEFINITION],
        fdsa_results_df[cols.ALL_JUNCTIONS],
        fdsa_results_df[cols.OVERLAPPING_JUNCTIONS]
    ))

    # Keep indexes updated if modifying column zip!
    _gene_id_index_in_zipped = 1
    _feature_number_index_in_zipped = 8

    frequency_column_names = fdsa_results_df.columns[
        fdsa_results_df.columns.str.contains(FDSA_RESULT_FREQUENCY_COLUMN_PREFIX)]

    zipped_frequency_columns = list(zip(
        *[fdsa_results_df[column_name] for column_name in frequency_column_names]
    ))

    # TODO: Investigate type warning (expected str got dict)
    return [
        FdsaResult(
            *[row[col_i] for col_i in range(len(zipped_info_columns[0]))],
            {
                frequency_column_names[freq_i].replace(FDSA_RESULT_FREQUENCY_COLUMN_PREFIX, ""): frequency
                for freq_i, frequency in enumerate(zipped_frequency_columns[row_i])
            },
            None if stringtie_results_df is None else calculate_fraction_lacking_feature(
                stringtie_results_df,
                row[_gene_id_index_in_zipped],
                int(row[_feature_number_index_in_zipped]),
                scale_to_percentage=True
            )
        )
        for row_i, row in enumerate(zipped_info_columns)
    ]


def load_fdsa_results(
    fdsa_results_path: str,
    samples: Dict[str, Sample],
    rank_by: str = "frequency",
    force_gene_names: Union[None, List[str]] = None,
    min_total_number_occurrences_across_all_samples: int = 1,
    min_per_sample_occurrences_number_occurrences: int = 0,
    min_per_sample_occurrences_in_at_least_n_samples: int = 0,
    fdsa_result_frequency_column_prefix: str = FDSA_RESULT_FREQUENCY_COLUMN_PREFIX,
    stringtie_results_path: Optional[str] = None,
    stringtie_remove_sample_name_ending: Optional[str] = None
) -> List[FdsaResult]:

    fdsa_results_df = load_fdsa_results_as_df(
        fdsa_results_path=fdsa_results_path,
        samples=samples,
        rank_by=rank_by,
        force_gene_names=force_gene_names,
        min_total_number_occurrences_across_all_samples=min_total_number_occurrences_across_all_samples,
        min_per_sample_occurrences_number_occurrences=min_per_sample_occurrences_number_occurrences,
        min_per_sample_occurrences_in_at_least_n_samples=min_per_sample_occurrences_in_at_least_n_samples,
        fdsa_result_frequency_column_prefix=fdsa_result_frequency_column_prefix
    )

    stringtie_results_df = None if stringtie_results_path is None else pd.read_csv(stringtie_results_path)

    if stringtie_results_df is not None:

        suffix = stringtie_remove_sample_name_ending \
            if "." not in stringtie_remove_sample_name_ending \
            else stringtie_remove_sample_name_ending[:stringtie_remove_sample_name_ending.index(".")]

        stringtie_results_df.rename(
            columns={
                col: col.replace(suffix, '')
                for col in [_col for _col in stringtie_results_df.columns if suffix in _col]
            },
            inplace=True
        )

    return convert_fdsa_results_df_to_objects(
        fdsa_results_df,
        stringtie_results_df
    )


def parse_junction_locus_by_id_from_string(
    junction_definition: str
) -> Dict[str, str]:
    junction_locus_by_id = {}
    for definition in junction_definition.split(" "):
        definition_components = definition.split("=")
        junction_locus_by_id[definition_components[0]] = definition_components[1][1:-1]  # Strip square brackets

    return junction_locus_by_id


def _get_sample_name_and_junction_string(s: str) -> Tuple[str, str]:

    components = [component.strip() for component in s.split(":")]
    components_length = len(components)

    if components_length > 2 or components_length < 1:
        raise ValueError(f"Failed to process junction string: {s}")

    junction_string = "" if components_length == 1 else components[1]

    return components[0], junction_string


def parse_junctions_from_string(
    junction_dump_string: str,
    junction_locus_by_id: Dict[str, str]
) -> Dict[str, List[QuantifiedSpliceJunctionLocus]]:

    junction_string_by_sample_name = {sample_name: junction_string for sample_name, junction_string in [
        _get_sample_name_and_junction_string(s) for s in junction_dump_string.split(" | ")
    ]}

    junctions_by_sample_name = {sample_name: [] for sample_name in junction_string_by_sample_name.keys()}
    for sample_name, junction_string in junction_string_by_sample_name.items():
        for quantified_junction in junction_string.split(" "):
            if not quantified_junction:
                continue
            quantified_junction_components = quantified_junction.split("*")
            junctions_by_sample_name[sample_name].append(
                QuantifiedSpliceJunctionLocus(
                    junction_locus_by_id[quantified_junction_components[1]],
                    int(quantified_junction_components[0])
                )
            )

    return junctions_by_sample_name
