
from typing import Dict, List, Tuple
import os
from dataclasses import dataclass

from src.config.parse_config import FaseRunConfig
from src.reporting.process_results import FaseResult


TRANSCRIPT_SECTION_TEMPLATE = "transcript_section.html"
TRANSCRIPT_SECTIONS_CSS = "transcript_sections.css"


@dataclass
class TranscriptSectionRenderInfo:
    section_id: str
    section_title: str
    gene_name: str
    transcript_id: str
    n_exons: int
    feature_name: str
    feature_no: int
    n_features: int
    expression_plot_uri: str
    splice_plot_uri: str


def transcript_section(render_info: TranscriptSectionRenderInfo) -> str:

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), TRANSCRIPT_SECTION_TEMPLATE), "r") as f:
        template = f.read()

    return template.format(**vars(render_info))


def transcript_sections(
    run_config: FaseRunConfig,
    fase_results: List[FaseResult],
    plots: Dict[str, Dict[str, str]],
    svg_scheme: str
) -> Tuple[str, List[TranscriptSectionRenderInfo]]:

    all_section_render_info = [
        TranscriptSectionRenderInfo(
            section_id=f"transcript-{fase_result.transcript_id}-{fase_result.feature_number}",
            section_title=f"{fase_result.gene_name} " +
                          f"({fase_result.feature_number} of {fase_result.total_features_in_transcript})",
            gene_name=fase_result.gene_name,
            transcript_id=fase_result.transcript_id,
            n_exons=len(fase_result.exon_positions),
            feature_name=run_config.feature_name,
            feature_no=fase_result.feature_number,
            n_features=fase_result.total_features_in_transcript,
            expression_plot_uri=svg_scheme + plots[
                f"{fase_result.transcript_id}-{fase_result.feature_number}"
            ]["expression"],
            splice_plot_uri=svg_scheme + plots[
                f"{fase_result.transcript_id}-{fase_result.feature_number}"
            ]["splice"]
        )
        for fase_result in fase_results if plots.get(
            f"{fase_result.transcript_id}-{fase_result.feature_number}", None
        ) is not None
    ]

    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), TRANSCRIPT_SECTIONS_CSS), "r") as f:
        style = "<style>\n" + f.read() + "\n</style>"

    sections_html = style + "\n".join([transcript_section(render_info) for render_info in all_section_render_info])

    return sections_html, all_section_render_info
