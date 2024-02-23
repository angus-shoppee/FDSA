
from typing import List, Dict
import os

from src.config.parse_config import FaseRunConfig
from src.reporting.process_results import FaseResult
from src.reporting.templates.results_table import results_table
from src.reporting.templates.transcript_sections import transcript_sections
from src.reporting.templates.table_of_contents import table_of_contents


REPORT_TEMPLATE = "templates/main_report.html"
REPORT_CSS = "templates/main_report.css"

SVG_SCHEME = "data:image/svg+xml;base64,"

JUMP_TO_RESULTS_TABLE_ID = "results-table-section"
JUMP_TO_TRANSCRIPTS_ID = "transcript-sections-intro"


def generate_html_report(
    # run_name: str,
    # feature_name: str,
    run_config: FaseRunConfig,
    fase_results: List[FaseResult],
    plots: Dict[str, Dict[str, str]],
    svg_scheme: str = SVG_SCHEME
) -> str:

    results_table_html = results_table(run_config, fase_results)

    sections_html, all_section_render_info = transcript_sections(
        run_config,
        fase_results,
        plots,
        svg_scheme
    )

    toc_html = table_of_contents(
        JUMP_TO_RESULTS_TABLE_ID,
        JUMP_TO_TRANSCRIPTS_ID,
        [render_info.section_id for render_info in all_section_render_info],
        [render_info.section_title for render_info in all_section_render_info]
    )

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), REPORT_TEMPLATE), "r") as f:
        template = f.read()

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), REPORT_CSS), "r") as f:
        style = "<style>\n" + f.read() + "\n</style>"

    return template.format(
        report_name=run_config.report_name,
        style=style,
        toc_html=toc_html,
        results_table_html=results_table_html,
        sections_html=sections_html
    )
