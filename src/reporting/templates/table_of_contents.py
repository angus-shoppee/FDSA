
import os
from typing import List


TABLE_OF_CONTENTS_TEMPLATE = "table_of_contents.html"


def table_of_contents(
    transcript_section_ids: List[str],
    transcript_section_titles: List[str]
) -> str:
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), TABLE_OF_CONTENTS_TEMPLATE), "r") as f:
        template = f.read()

    transcript_entries = (
        "\n".join([
            f"""<li><a href="#{section_id}">{section_title}</a></li>"""
            for section_id, section_title in zip(transcript_section_ids, transcript_section_titles)
        ])
    )

    return template.format(
        transcript_entries=transcript_entries
    )
