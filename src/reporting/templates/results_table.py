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


from typing import List, Union
import os

from src.reporting.process_results import FaseResult
from src.config.parse_config import FaseRunConfig


RESULTS_TABLE_TEMPLATE = "results_table.html"
RESULTS_TABLE_CSS = "results_table.css"


def _generate_results_criteria_text(run_config: FaseRunConfig) -> str:

    # _s1 and _s2 control plurality of term "splice event/events" in rendered text
    _s1 = "" if run_config.report_min_total_n_occurrences_across_all_samples == 1 else "s"
    _s2 = "" if run_config.report_min_n_occurrences_in_sample == 1 else "s"

    return (
        f"{run_config.report_min_total_n_occurrences_across_all_samples} total splice event{_s1} across all samples, " +
        f"{run_config.report_min_n_occurrences_in_sample} splice event{_s2} per sample in at least " +
        f"{run_config.report_occurrences_in_at_least_n_samples} samples."
    )


def _get_results_row(fase_result: FaseResult) -> List[Union[str, int, float]]:

    numbers = [
        sum([junction.n for junction in junctions])
        for junctions in fase_result.overlapping.values()
    ]
    frequencies = list(fase_result.frequencies.values())

    average_number = round(sum(numbers) / len(numbers), 2)
    average_frequency = round(sum(frequencies) / len(frequencies), 2)

    return [
        fase_result.gene_name, fase_result.feature_number, fase_result.total_features_in_transcript,
        average_number, average_frequency
    ] + numbers + frequencies


def _generate_table_contents(
    table_header: List[str],
    table_data: List[List[str]]
) -> str:

    return (
        "<tr>" +
        f"<th class=\"gene-name-column\">{table_header[0]}</th>" +
        "".join([f"<th>{column_label}</th>" for column_label in table_header[1:]]) +
        "</tr>\n" +
        "\n".join([
            "<tr>" +
            f"<td class=\"gene-name-column\">{row[0]}</td>" +
            "".join([f"<td>{value}</td>" for value in row[1:]]) +
            "</tr>"
            for row in table_data
        ])
    )


def results_table(
    run_config: FaseRunConfig,
    fase_results: List[FaseResult]
) -> str:

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), RESULTS_TABLE_TEMPLATE), "r") as f:
        template = f.read()

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), RESULTS_TABLE_CSS), "r") as f:
        style = "<style>\n" + f.read() + "\n</style>"

    results_criteria_text = _generate_results_criteria_text(run_config)

    table_header = [
       "Gene name", "Feature #", "Total features", "Avg. N", "Avg. %"
    ] + list(fase_results[0].overlapping.keys()) + [
       f"% {sample_name}"
       for sample_name in fase_results[0].frequencies.keys()
    ]
    table_data = [
        _get_results_row(fase_result)
        for fase_result in fase_results
    ]

    table_contents = _generate_table_contents(table_header, table_data)

    return template.format(
        style=style,
        results_criteria_text=results_criteria_text,
        table_contents=table_contents
    )
