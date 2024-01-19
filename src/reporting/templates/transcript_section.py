
import os
from dataclasses import dataclass


TRANSCRIPT_SECTION_TEMPLATE = "transcript_section.html"


@dataclass
class TranscriptSectionRenderInfo:
    section_id: str
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
