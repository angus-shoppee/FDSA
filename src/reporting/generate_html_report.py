
from typing import List, Dict
import os

from src.reporting.process_results import FaseResult
from src.reporting.templates.transcript_section import transcript_section, TranscriptSectionRenderInfo
from src.reporting.templates.table_of_contents import table_of_contents


REPORT_TEMPLATE = "templates/main_report.html"
REPORT_CSS = "templates/main_report.css"

SVG_SCHEME = "data:image/svg+xml;base64,"

JUMP_TO_RESULTS_TABLE_ID = "results_table_container"


def generate_html_report(
    run_name: str,
    feature_name: str,
    fase_results: List[FaseResult],
    plots: Dict[str, Dict[str, str]],
    svg_scheme: str = SVG_SCHEME
) -> str:

    all_render_info = [
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
        [render_info.section_id for render_info in all_render_info],
        [render_info.section_title for render_info in all_render_info]
    )

    sections_html = "\n".join(transcript_section(render_info) for render_info in all_render_info)

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
