
from templates.transcript_section import transcript_section, TranscriptSectionRenderInfo


SVG_SCHEME = "data:image/svg+xml;base64,"


if __name__ == "__main__":

    html_snippet = transcript_section(TranscriptSectionRenderInfo(
        section_id="transcript",
        gene_name="Abc1",
        transcript_id="ENSMUS00000000001",
        n_exons=3,
        feature_name="transmembrane region",
        feature_no=1,
        n_features=1,
        expression_plot_uri=SVG_SCHEME,
        splice_plot_uri=SVG_SCHEME
    ))

    print(html_snippet)

