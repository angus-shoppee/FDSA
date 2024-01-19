
# FEATURE ANNOTATION DATA (VIA GENBANK)

from dataclasses import dataclass
from typing import List, Dict


def get_unique_child(elem, child_name):
    """Find a child element of a given element in an XML tree, intended for one-off children"""

    x = [c for c in elem.iter(child_name)]
    if not x:
        raise ValueError("No matching children found")
    if len(x) > 1:
        raise ValueError("Multiple matching children found")
    return x[0]


@dataclass
class RefseqExonFeature:
    start: int
    end: int


class RefseqExon:
    start: int
    end: int
    length: int

    def __init__(self, start: int, end: int):

        self.start = start
        self.end = end
        self.length = 1 + end - start


@dataclass
class GBQual:
    name: str
    value: str


class GBFeature:
    key: str
    start: int
    end: int
    quals: List[GBQual]

    def __init__(
            self,
            feature_key: str,
            feature_location: str,
            feature_quals: List[GBQual]
    ):
        self.key = feature_key
        self.start = int(feature_location) if ".." not in feature_location else int(feature_location.split("..")[0])
        self.end = int(feature_location) if ".." not in feature_location else int(feature_location.split("..")[1])
        self.quals = feature_quals

    def has_qual_value_containing(
            self,
            value_substring: str
    ) -> bool:
        for qual in self.quals:
            if not qual.value:
                return False
            if value_substring in qual.value:
                return True
        return False


class GBSeq:
    refseq_id: str
    features: List[GBFeature]
    exons: List[RefseqExon]
    _analysis_features: dict

    def __init__(
            self,
            refseq_id: str,
            features: List[GBFeature]
    ):
        self.refseq_id = refseq_id
        self.features = features
        self._analysis_features = {}

        self.exons = [RefseqExon(
            start=feature.start,
            end=feature.end
        ) for feature in self.features if feature.key == "exon"]

    def get_features_by_key(
            self,
            key: str
    ) -> List[GBFeature]:
        return [feature for feature in self.features if feature.key == key]

    def get_analysis_features(self) -> Dict:
        return self._analysis_features

    def set_analysis_feature(
            self,
            feature_id: str,
            gbfeature: GBFeature
    ) -> None:

        # Create feature list if it doesn't exist yet
        if feature_id not in self._analysis_features.keys():
            self._analysis_features[feature_id] = []

        # Append to feature list
        self._analysis_features[feature_id].append(
            RefseqExonFeature(
                start=gbfeature.start,
                end=gbfeature.end
            )
        )
