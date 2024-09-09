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


from typing import Dict, List, Tuple, Optional
import os
import re
from dataclasses import dataclass

from config.parse_config import ProgramRunConfig
from reporting.process_results import FdsaResult


TRANSCRIPT_SECTIONS_TEMPLATE = "transcript_sections.html"
TRANSCRIPT_SECTION_UNIT_TEMPLATE = "transcript_section_unit.html"
TRANSCRIPT_SECTIONS_CSS = "transcript_sections.css"

MAX_PLOTTED_TEXT = """ The number of transcripts to visualise has been limited to {max_plotted}."""
MAX_PER_GROUP_TEXT = """ The number of splice plots per experimental group has been limited to {max_per_group}."""
JUNCTION_LIMIT_TEXT = """ Splice plots have been limited to show non-feature-overlapping junctions with a minimum 
                          of {junction_limit} occurrences."""


@dataclass
class TranscriptSectionRenderInfo:
    section_id: str
    section_title: str
    gene_name: str
    gene_id: str
    n_exons: int
    feature_name: str
    feature_qualifiers: str
    feature_no: int
    n_features: int
    expression_plot_uri: Optional[str]
    splice_plot_uri: str


def transcript_unit_section(
    html_template: str,
    render_info: TranscriptSectionRenderInfo
) -> str:

    formatted_html = html_template.format(**vars(render_info))

    # Delete expression plot container div if image is null
    null_src_regex_pattern = r"""<div class="expression-plot-container">.*?<img src=""[^>]*>.*?</div>"""
    return re.sub(null_src_regex_pattern, '', formatted_html, flags=re.DOTALL)


def transcript_sections(
    run_config: ProgramRunConfig,
    fdsa_results: List[FdsaResult],
    plots: Dict[str, Dict[str, Optional[str]]],
    svg_scheme: str
) -> Tuple[str, List[TranscriptSectionRenderInfo]]:

    all_section_render_info = [
        TranscriptSectionRenderInfo(
            section_id=f"transcript-{fdsa_result.transcript_id}-{fdsa_result.feature_number}",
            section_title=f"{fdsa_result.gene_name} " +
                          f"({fdsa_result.feature_number} of {fdsa_result.total_features_in_transcript})",
            gene_name=fdsa_result.gene_name,
            gene_id=fdsa_result.gene_id,
            n_exons=len(fdsa_result.exon_positions),
            feature_name=run_config.feature_name,
            feature_qualifiers=fdsa_result.feature_qualifiers,
            feature_no=fdsa_result.feature_number,
            n_features=fdsa_result.total_features_in_transcript,
            expression_plot_uri="" if plots[
                    f"{fdsa_result.transcript_id}-{fdsa_result.feature_number}"
                ]["expression"] is None else (
                svg_scheme + plots[
                    f"{fdsa_result.transcript_id}-{fdsa_result.feature_number}"
                ]["expression"]
            ),
            splice_plot_uri=svg_scheme + plots[
                f"{fdsa_result.transcript_id}-{fdsa_result.feature_number}"
            ]["splice"]
        )
        for fdsa_result in fdsa_results if plots.get(
            f"{fdsa_result.transcript_id}-{fdsa_result.feature_number}", None
        ) is not None
    ]

    max_plotted_text = "" if run_config.report_max_n_plotted is None else \
        MAX_PLOTTED_TEXT.format(max_plotted=run_config.report_max_n_plotted)

    max_per_group_text = "" if run_config.report_transcript_plot_max_samples_per_group is None else \
        MAX_PER_GROUP_TEXT.format(max_per_group=run_config.report_transcript_plot_max_samples_per_group)

    junction_limit_text = "" if run_config.report_draw_junctions_min_count in (0, 1) else \
        JUNCTION_LIMIT_TEXT.format(junction_limit=run_config.report_draw_junctions_min_count)

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), TRANSCRIPT_SECTIONS_CSS), "r") as f:
        style = "<style>\n" + f.read() + "\n</style>"

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), TRANSCRIPT_SECTIONS_TEMPLATE), "r") as f:
        template = f.read()

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), TRANSCRIPT_SECTION_UNIT_TEMPLATE), "r") as f:
        section_unit_template = f.read()

    sections_html = (
        template.format(
            style=style,
            rank_by=run_config.rank_results_by,
            max_plotted_text=max_plotted_text,
            max_per_group_text=max_per_group_text,
            junction_limit_text=junction_limit_text,
            transcript_unit_sections="\n".join([
                transcript_unit_section(section_unit_template, render_info)
                for render_info in all_section_render_info
            ])
        )
    )

    return sections_html, all_section_render_info
