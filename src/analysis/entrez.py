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


# ENTREZ API

from typing import Dict, List, Any
import os
import json
import entrezpy.base.result
import entrezpy.base.analyzer
import entrezpy.conduit
from socket import timeout as socket_timeout
from xml.etree import cElementTree
from xml.etree.ElementTree import Element

from src.analysis.annotation import GBSeq, GBFeature, GBQual


class FetchResult(entrezpy.base.result.EutilsResult):
    xml: str
    fetched: List[Any]  # Inherited attribute

    def __init__(self, response, request):
        super().__init__(request.eutil, request.query_id, request.db)
        self.xml = ""

    def size(self):
        return 1 if self.xml else 0

    def isEmpty(self):
        return True if self.xml else False

    def get_link_parameter(self, reqnum=0):
        print("{} has no elink capability".format(self))
        return {}

    def dump(self):
        return {
            self: {
                'dump:': {
                    'xml': self.xml,
                    'query_id': self.query_id,
                    'db': self.db,
                    'eutil': self.function
                }
            }
        }

    def add_record(self, rec):
        self.fetched.append(rec)


class RecordAnalyzer(entrezpy.base.analyzer.EutilsAnalyzer):

    def __init__(self):
        super().__init__()

    def init_result(self, response, request):
        if self.result is None:
            # Create FetchResult instance to handle/store results
            self.result = FetchResult(response, request)

    def analyze_error(self, response, request):
        print(json.dumps({__name__: {'Response': {'dump': request.dump(), 'error': response}}}))

    def analyze_result(self, response, request):
        self.init_result(response, request)
        self.result.xml = response.getvalue()


# TODO: Refactor use of XML file format to genbank format (>10x more storage efficient, faster parsing with gtfparse)
def download_and_save_feature_annotation_xml(
    refseq_ids_to_search: List[str],
    output_path: str,
    user_email: str,
    chunk_size: int = 100,
    max_retry_attempts_per_chunk: int = 5,
    force_restart_download: bool = False
) -> None:

    progress_file_path = os.path.join(os.path.dirname(output_path), ".download_progress.json")

    if os.path.exists(progress_file_path) and not force_restart_download:
        # Resume previously started download
        with open(progress_file_path, "r") as progress_file:
            print("Resuming previously download")
            progress = json.load(progress_file)
            chunk_size = progress["chunk_size"]  # Override chunk_size to remain consistent
    else:
        # Download all data
        progress = {
            "current_chunk_no": 1,
            "chunk_size": chunk_size
        }
        # Wipe previous output file if it exists
        open(output_path, "w").close()
        with open(progress_file_path, "w") as progress_file:
            json.dump(
                progress,
                progress_file
            )

    # Split input into chunks of specified size
    batch_size = len(refseq_ids_to_search)
    chunks = []
    _chunk_offset = 0
    while _chunk_offset < batch_size:
        chunks.append(refseq_ids_to_search[_chunk_offset:min(batch_size, _chunk_offset + chunk_size)])
        _chunk_offset += chunk_size

    number_of_chunks = len(chunks)
    chunks_iterable = iter(chunks)

    with open(output_path, "a") as f:

        chunk = next(chunks_iterable, [])
        current_chunk_no = 1
        retry_attempt = 0

        while chunk:

            if progress["current_chunk_no"] > current_chunk_no:
                chunk = next(chunks_iterable)
                current_chunk_no += 1
                continue

            if retry_attempt > max_retry_attempts_per_chunk:
                raise RuntimeError("Could not download chunk from entrez - max attempts reached")

            print(f"Downloading chunk {current_chunk_no} of {number_of_chunks}")
            print(f"{chunk[:5]} ...")

            c = entrezpy.conduit.Conduit(user_email)

            fetch_features = c.new_pipeline()

            query = ", ".join([refseq for refseq in chunk])

            sid = fetch_features.add_search({
                "db": "nuccore",
                "term": query
            })

            fetch_features.add_fetch(
                {
                    "rettype": "genbank",
                    # "retmode": "gbwithparts"
                    "retmode": "xml"
                },
                dependency=sid,
                analyzer=RecordAnalyzer()
            )

            # If query successfully returns a result, write to file; otherwise retry
            try:
                a = c.run(fetch_features)

                # Append data for this chunk to XML file, trimming to leave only one opening and closing GBSet tag
                xml_string_buffer = a.result.xml
                if current_chunk_no > 1:
                    xml_string_buffer = xml_string_buffer[xml_string_buffer.index("<GBSet>") + len("<GBSet>"):]
                if current_chunk_no < number_of_chunks:
                    xml_string_buffer = xml_string_buffer[:xml_string_buffer.index("</GBSet>")]

                f.write(xml_string_buffer)

            except (AttributeError, TypeError, socket_timeout) as e:
                retry_attempt += 1
                print(
                    f"(Caught exception: {e})\n" +
                    f"Query did not return XML data." +
                    f"Retrying (attempt {retry_attempt} of {max_retry_attempts_per_chunk})"
                )
                continue

            # Proceed to next chunk, or break if final chunk reached
            try:
                chunk = next(chunks_iterable)
                current_chunk_no += 1
                retry_attempt = 0
                with open(progress_file_path, "w") as progress_file:
                    json.dump(
                        {
                            "current_chunk_no": current_chunk_no,
                            "chunk_size": chunk_size
                        },
                        progress_file
                    )

            except StopIteration:
                break

        # Remove download progress file
        os.remove(progress_file_path)
        print("DEBUG: Removed progress file")


def _get_unique_child(
    elem: Element,
    child_name: str
) -> Element:

    x = [c for c in elem.iter(child_name)]
    if not x:
        raise ValueError("No matching child elements found")
    if len(x) > 1:
        raise ValueError("Multiple matching child elements found")
    return x[0]


def get_gbseq_from_xml(
    xml_file_path: str
) -> Dict[str, GBSeq]:

    """
    Creates GBSeq objects from feature annotation XML.
    Returns a dictionary of GBSeq objects against refseq IDs.
    """

    print("Reading feature annotation data from XML file...")

    tree = cElementTree.parse(xml_file_path)
    root = tree.getroot()

    print("...done\n")
    print("Creating feature annotation records...")

    gbseq_by_refseq = {}
    skipped_feature_invalid_region, skipped_feature_missing_gbqual = 0, 0
    _i = 0
    _t = len(root)

    for gbseq in root:

        if _i % 1000 == 0:
            print(f"Progress: ({_i}/{_t})")
        _i += 1

        refseq_id = _get_unique_child(gbseq, "GBSeq_locus").text

        features = []

        for gbfeature in gbseq.iter("GBFeature"):

            quals = []

            for gbqual in gbfeature.iter("GBQualifier"):

                try:
                    name = _get_unique_child(gbqual, "GBQualifier_name").text
                    value = _get_unique_child(gbqual, "GBQualifier_value").text
                except ValueError:
                    skipped_feature_missing_gbqual += 1
                    continue

                if name == "gene" or name == "gene_synonym":
                    continue

                quals.append(GBQual(name=name, value=value))

            loc = _get_unique_child(gbfeature, "GBFeature_location").text
            loc = loc.replace("<", "").replace(">", "")

            try:
                # Check that n is numeric - int conversion will fail with ValueError if not
                if ".." in loc:
                    for n in loc.split(".."):
                        int(n)
                else:
                    int(loc)

                features.append(GBFeature(
                    _get_unique_child(gbfeature, "GBFeature_key").text,
                    loc,
                    quals
                ))

            except ValueError:
                skipped_feature_invalid_region += 1

        gbseq_by_refseq[refseq_id] = GBSeq(
            refseq_id,
            features
        )

    print(
        f"...finished creating {len(gbseq_by_refseq)} records " +
        f"({skipped_feature_invalid_region} features skipped due to invalid region, "
        f"{skipped_feature_missing_gbqual} features skipped due to missing GBQualifier)\n"
    )

    return gbseq_by_refseq
