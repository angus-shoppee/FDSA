# Feature-directed Analysis of Splice Events (FASE)
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


from typing import List, Tuple, Dict, Any, Union
import csv
import pandas as pd

from config import output_column_names as cols
from analysis.experiment import Sample


FASE_RESULT_FREQUENCY_COLUMN_PREFIX = "Percent "


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


class FaseResult:
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
        frequencies: Dict[str, float]
        # frequencies: List[Dict[str, float]]
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

    def total_number_occurrences_across_all_samples(self) -> int:

        return sum([sum([locus.n for locus in loci] for loci in self.overlapping.values())])


def load_fase_results_as_df(
    fase_results_path: str,
    samples: Dict[str, Sample],
    rank_by: str = "frequency",
    force_gene_names: Union[None, List[str]] = None,
    min_total_number_occurrences_across_all_samples: int = 1,
    min_per_sample_occurrences_number_occurrences: int = 0,
    min_per_sample_occurrences_in_at_least_n_samples: int = 0,
    fase_result_frequency_column_prefix: str = FASE_RESULT_FREQUENCY_COLUMN_PREFIX
) -> pd.DataFrame:

    if force_gene_names is None:
        force_gene_names = list()

    with open(fase_results_path, "r") as f:

        results = csv.reader(f)

        header = next(results)

        raw_number_column_names = [sample_name for sample_name in samples.keys()]
        frequency_column_names = [fase_result_frequency_column_prefix + sample_name for sample_name in samples.keys()]

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
                "Unable to find matching sample data in FASE results. Check that BAM files match FASE output.\n" +
                f"Sample names from BAM files: {samples.keys()}\n" +
                f"Column names in FASE output: {header}"
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

        fase_results_df = pd.DataFrame(
            rows_to_import,
            columns=header
        )

        # fase_results_df = pd.DataFrame(
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

        if len(fase_results_df) == 0:
            raise ValueError("No results were loaded - check that specified criteria are not too stringent")

        # Convert integer columns
        fase_results_df[cols.TRANSCRIPT_START] = fase_results_df[cols.TRANSCRIPT_START].astype(int)
        fase_results_df[cols.N_FEATURES_IN_TRANSCRIPT] = fase_results_df[cols.N_FEATURES_IN_TRANSCRIPT].astype(int)
        fase_results_df[cols.FEATURE_NUMBER] = fase_results_df[cols.FEATURE_NUMBER].astype(int)
        for column_name in raw_number_column_names:
            fase_results_df[column_name] = fase_results_df[column_name].astype(int)

        # Convert float columns
        fase_results_df[cols.AVG_NUMBER] = fase_results_df[cols.AVG_NUMBER].astype(float)
        fase_results_df[cols.AVG_FREQUENCY] = fase_results_df[cols.AVG_FREQUENCY].astype(float)
        for column_name in frequency_column_names:
            fase_results_df[column_name] = fase_results_df[column_name].astype(float)

        # Order results by either number or avg frequency, depending on report config
        if rank_by == "frequency":
            fase_results_df = fase_results_df.sort_values(
                by=[cols.AVG_FREQUENCY, cols.AVG_NUMBER],
                ascending=[False, False]
            )
        elif rank_by == "number":
            fase_results_df = fase_results_df.sort_values(
                by=[cols.AVG_NUMBER, cols.AVG_FREQUENCY],
                ascending=[False, False]
            )
        else:
            raise ValueError(f"Invalid value for rank_by specified (\"{rank_by}\"). Options are frequency/number")

        return fase_results_df


def convert_fase_results_df_to_objects(fase_results_df: pd.DataFrame) -> List[FaseResult]:

    zipped_info_columns = list(zip(
        fase_results_df[cols.TRANSCRIPT_ID],
        fase_results_df[cols.GENE_ID],
        fase_results_df[cols.GENE_NAME],
        fase_results_df[cols.CHROMOSOME],
        fase_results_df[cols.STRAND],
        fase_results_df[cols.EXON_POSITIONS],
        fase_results_df[cols.FEATURE_QUALIFIERS],
        fase_results_df[cols.FEATURE_REGION],
        fase_results_df[cols.FEATURE_NUMBER],
        fase_results_df[cols.N_FEATURES_IN_TRANSCRIPT],
        fase_results_df[cols.JUNCTION_DEFINITION],
        fase_results_df[cols.ALL_JUNCTIONS],
        fase_results_df[cols.OVERLAPPING_JUNCTIONS]
    ))

    frequency_column_names = fase_results_df.columns[
        fase_results_df.columns.str.contains(FASE_RESULT_FREQUENCY_COLUMN_PREFIX)]

    zipped_frequency_columns = list(zip(
        *[fase_results_df[column_name] for column_name in frequency_column_names]
    ))

    # TODO: Investigate type warning (expected str got dict)
    return [
        FaseResult(
            *[row[col_i] for col_i in range(len(zipped_info_columns[0]))],
            {
                frequency_column_names[freq_i].replace(FASE_RESULT_FREQUENCY_COLUMN_PREFIX, ""): frequency
                for freq_i, frequency in enumerate(zipped_frequency_columns[row_i])
            }
        ) for row_i, row in enumerate(zipped_info_columns)
    ]


def load_fase_results(
    fase_results_path: str,
    samples: Dict[str, Sample],
    rank_by: str = "frequency",
    force_gene_names: Union[None, List[str]] = None,
    min_total_number_occurrences_across_all_samples: int = 1,
    min_per_sample_occurrences_number_occurrences: int = 0,
    min_per_sample_occurrences_in_at_least_n_samples: int = 0,
    fase_result_frequency_column_prefix: str = FASE_RESULT_FREQUENCY_COLUMN_PREFIX
) -> List[FaseResult]:

    fase_results_df = load_fase_results_as_df(
        fase_results_path=fase_results_path,
        samples=samples,
        rank_by=rank_by,
        force_gene_names=force_gene_names,
        min_total_number_occurrences_across_all_samples=min_total_number_occurrences_across_all_samples,
        min_per_sample_occurrences_number_occurrences=min_per_sample_occurrences_number_occurrences,
        min_per_sample_occurrences_in_at_least_n_samples=min_per_sample_occurrences_in_at_least_n_samples,
        fase_result_frequency_column_prefix=fase_result_frequency_column_prefix
    )

    return convert_fase_results_df_to_objects(fase_results_df)


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
