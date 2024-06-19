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


import os
from typing import List


TABLE_OF_CONTENTS_TEMPLATE = "table_of_contents.html"


def table_of_contents(
    results_table_id: str,
    transcripts_sections_intro_id: str,
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
        results_table_id=results_table_id,
        transcripts_sections_intro_id=transcripts_sections_intro_id,
        transcript_entries=transcript_entries
    )
