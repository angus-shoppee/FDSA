
from typing import List, Dict, Union
import os

from src.reporting.process_results import FaseResult
from src.reporting.templates.transcript_section import transcript_section, TranscriptSectionRenderInfo
from src.reporting.templates.table_of_contents import table_of_contents


REPORT_TEMPLATE = "templates/main_report.html"
REPORT_CSS = "templates/main_report.css"

SVG_SCHEME = "data:image/svg+xml;base64,"

JUMP_TO_RESULTS_TABLE_ID = "results_table_container"


def _get_results_row(fase_result: FaseResult) -> List[Union[str, int, float]]:

    numbers = [
        sum([junction.n for junction in junctions])
        for junctions in fase_result.overlapping.values()
    ]
    frequencies = list(fase_result.frequencies.values())

    return [
        fase_result.gene_name, fase_result.feature_number, fase_result.total_features_in_transcript
    ] + numbers + frequencies


def generate_html_report(
    run_name: str,
    feature_name: str,
    fase_results: List[FaseResult],
    plots: Dict[str, Dict[str, str]],
    svg_scheme: str = SVG_SCHEME
) -> str:

    results_table_header = [
        "Gene name", "Feature #", "Total features"
    ] + list(fase_results[0].overlapping.keys()) + [
        f"% {sample_name}"
        for sample_name in fase_results[0].frequencies.keys()
    ]
    results_table_data = [
        _get_results_row(fase_result)
        for fase_result in fase_results
    ]

    print(results_table_header)
    print(results_table_data)

    all_section_render_info = [
        TranscriptSectionRenderInfo(
            section_id=f"transcript-{fase_result.transcript_id}-{fase_result.feature_number}",
            section_title=f"{fase_result.gene_name} " +
                          f"({fase_result.feature_number} of {fase_result.total_features_in_transcript})",
            gene_name=fase_result.gene_name,
            transcript_id=fase_result.transcript_id,
            n_exons=len(fase_result.exon_positions),
            feature_name=feature_name,
            feature_no=fase_result.feature_number,
            n_features=fase_result.total_features_in_transcript,
            expression_plot_uri=svg_scheme+plots[
                f"{fase_result.transcript_id}-{fase_result.feature_number}"
            ]["expression"],
            splice_plot_uri=svg_scheme+plots[
                f"{fase_result.transcript_id}-{fase_result.feature_number}"
            ]["splice"]
        )
        for fase_result in fase_results if plots.get(
            f"{fase_result.transcript_id}-{fase_result.feature_number}", None
        ) is not None
    ]

    toc_html = table_of_contents(
        JUMP_TO_RESULTS_TABLE_ID,
        [render_info.section_id for render_info in all_section_render_info],
        [render_info.section_title for render_info in all_section_render_info]
    )

    sections_html = "\n".join(transcript_section(render_info) for render_info in all_section_render_info)

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), REPORT_TEMPLATE), "r") as f:
        template = f.read()

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), REPORT_CSS), "r") as f:
        style = "<style>\n" + f.read() + "\n</style>"

    return template.format(
        run_name=run_name,
        style=style,
        toc_html=toc_html,
        sections_html=sections_html
    )
