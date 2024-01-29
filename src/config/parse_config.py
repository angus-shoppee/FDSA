
from typing import Union, List, Dict, Any
from functools import reduce
import os
import configparser


RANK_BY_ALLOWED_VALUES = ["frequency", "number"]


def flatten_nested_lists(nested_lists: List[List[Any]]) -> List[Any]:
    return reduce(lambda a, b: a + b, nested_lists)


def parse_bool(value: Union[bool, str]) -> bool:

    if isinstance(value, bool):
        return value

    if value.lower() == "true":
        return True
    elif value.lower() == "false":
        return False

    raise ValueError(f"Could not parse bool from input value: {value}. Allowed values are True/False.")


class FaseInternalConfig:
    allowed_species: List[str]
    biomart_name_for_species: Dict[str, str]
    default_mapq_for_unique_mapping: int
    default_feature_junction_overlap_threshold: float
    default_max_n_features_in_transcript: int
    default_only_use_longest_annotated_transcript: bool
    default_skip_transcripts_with_redundant_feature_annotation: bool
    # ^ skip_transcripts_with_redundant_feature_annotation is overridden by only_use_longest_annotated_transcript=True
    default_include_all_junctions_in_output: bool
    default_generate_report: bool
    default_point_color_in_plot: str
    default_feature_counts_primary_alignment_only: bool

    def __init__(self, internal_config_path: str) -> None:

        # WARNING: Internal config is exposed for advanced user control. Deletions will trigger a ValueError here,
        # but improperly formatted changes to values may not be caught and will break program behaviour

        _e = f"Error while parsing internal config at {internal_config_path}: "

        config = configparser.ConfigParser()
        config.read(internal_config_path)

        # [BUILD]

        allowed_species_string = config.get("BUILD", "allowedSpecies", fallback=None)
        if allowed_species_string is None:
            raise ValueError(_e + f"Missing mandatory parameter \"allowedSpecies\" in section BUILD")
        self.allowed_species = allowed_species_string.split(" ")

        biomart_name_for_species_string = config.get("BUILD", "biomartNameForSpecies", fallback=None)
        if biomart_name_for_species_string is None:
            raise ValueError(_e + f"Missing mandatory parameter \"biomartNameForSpecies\" in section BUILD")
        self.biomart_name_for_species = {
            key: value
            for key, value in [
                s.split(":") for s in biomart_name_for_species_string.split(" ")
            ]
        }
        for species_name in self.allowed_species:
            if species_name not in self.biomart_name_for_species.keys():
                raise ValueError(_e + f"The species \"{species_name}\" specified in allowedSpecies is missing a "
                                      f"corresponding key:value entry in biomartNameForSpecies")

        # [ANALYSIS]

        default_mapq_for_unique_mapping = config.get(
            "ANALYSIS", "defaultMapqForUniqueMapping", fallback=None
        )
        if default_mapq_for_unique_mapping is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultMapqForUniqueMapping\" in section ANALYSIS")
        self.default_mapq_for_unique_mapping = int(default_mapq_for_unique_mapping)

        default_feature_junction_overlap_threshold = config.get(
            "ANALYSIS", "defaultFeatureJunctionOverlapThreshold", fallback=None
        )
        if default_feature_junction_overlap_threshold is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultFeatureJunctionOverlapThreshold\" in section "
                                  f"ANALYSIS")
        self.default_feature_junction_overlap_threshold = float(default_feature_junction_overlap_threshold)

        default_max_n_features_in_transcript = config.get(
            "ANALYSIS", "defaultMaxNFeaturesInTranscript", fallback=None
        )
        if default_max_n_features_in_transcript is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultMaxNFeaturesInTranscript\" in section ANALYSIS")
        self.default_max_n_features_per_transcript = int(default_max_n_features_in_transcript)

        default_only_use_longest_annotated_transcript = config.get(
            "ANALYSIS", "defaultOnlyUseLongestAnnotatedTranscript", fallback=None
        )
        if default_only_use_longest_annotated_transcript is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultOnlyUseLongestAnnotatedTranscript\" in section "
                                  f"ANALYSIS")
        self.default_only_use_longest_annotated_transcript = parse_bool(
            default_only_use_longest_annotated_transcript
        )

        default_skip_transcripts_with_redundant_feature_annotation = config.get(
            "ANALYSIS", "defaultSkipTranscriptsWithRedundantFeatureAnnotation", fallback=None
        )
        if default_skip_transcripts_with_redundant_feature_annotation is None:
            raise ValueError(_e + f"Missing mandatory parameter "
                                  f"\"defaultSkipTranscriptsWithRedundantFeatureAnnotation\" in section ANALYSIS")
        self.default_skip_transcripts_with_redundant_feature_annotation = parse_bool(
            default_skip_transcripts_with_redundant_feature_annotation
        )

        default_include_all_junctions_in_output = config.get(
            "ANALYSIS", "defaultIncludeAllJunctionsInOutput", fallback=None
        )
        if default_include_all_junctions_in_output is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultIncludeAllJunctionsInOutput\" in section "
                                  f"ANALYSIS")
        self.default_include_all_junctions_in_output = parse_bool(default_include_all_junctions_in_output)

        # [REPORT]

        default_generate_report = config.get(
            "REPORT", "defaultGenerateReport", fallback=None
        )
        if default_generate_report is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultGenerateReport\" in section REPORT")
        self.default_generate_report = parse_bool(default_generate_report)

        default_point_color_in_plot = config.get(
            "REPORT", "defaultPointColorInPlot", fallback=True
        )
        if default_point_color_in_plot is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultPointColorInPlot\" in section REPORT")
        self.default_point_color_in_plot = default_point_color_in_plot

        default_feature_counts_primary_alignment_only = config.get(
            "REPORT", "defaultFeatureCountsPrimaryAlignmentOnly", fallback=None
        )
        if default_feature_counts_primary_alignment_only is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultFeatureCountsPrimaryAlignmentOnly\" in "
                                  f"section REPORT")
        self.default_feature_counts_primary_alignment_only = parse_bool(
            default_feature_counts_primary_alignment_only
        )


class FaseUserConfig:

    pass


class FaseRunConfig:
    run_name: str
    feature_name: str
    output_path: str
    input_path: str
    species: str
    genome: str
    n_threads: int
    max_n_features_in_transcript: int
    only_use_longest_annotated_transcript: bool
    skip_transcripts_with_redundant_feature_annotation: bool  # overridden by only_use_longest_annotated_transcript=True
    generate_report: bool
    report_max_n_plotted: Union[None, int]
    report_min_total_number_occurrences_across_all_samples: int
    report_min_n_occurrences_in_sample: int
    report_occurrences_in_at_least_n_samples: int
    primary_alignment_only: bool
    mapq_for_unique_mapping: int
    sample_groups: Union[None, Dict[str, List[str]]]
    group_name_by_sample_name: Union[None, Dict[str, str]]  # Interpolated, not defined by the user
    color_by_group_name: Union[None, Dict[str, str]]

    def __init__(
        self,
        run_config_path: str,
        internal_config: FaseInternalConfig
    ) -> None:

        _e = f"Error while parsing run config at {run_config_path}: "

        run_config = configparser.ConfigParser()
        run_config.read(run_config_path)

        if "RUN" not in run_config:
            raise ValueError(_e + f"Missing mandatory section RUN (run settings)")

        # [RUN] - Mandatory parameters

        run_name = run_config.get("RUN", "name", fallback=None)
        if run_name is None:
            raise ValueError(_e + f"Missing mandatory parameter \"name\" (name for this run) in section RUN")
        self.run_name = run_name

        feature_name = run_config.get("RUN", "feature", fallback=None)
        if feature_name is None:
            raise ValueError(_e + f"Missing mandatory parameter \"feature\" (feature annotation name/substring to "
                                  f"match) in section RUN")
        self.feature_name = feature_name

        output_path = run_config.get("RUN", "output", fallback=None)
        if output_path is None:
            raise ValueError(_e + f"Missing mandatory parameter \"output\" (output path for results) in section RUN")
        if os.path.exists(output_path) and not os.path.isdir(output_path):
            raise ValueError(_e + f"The specified output path exists but it is not a directory: {output_path}")
        self.output_path = output_path

        input_path = run_config.get("RUN", "input", fallback=None)
        if input_path is None:
            raise ValueError(_e + f"Missing mandatory parameter \"input\" (path to folder containing BAM files) in "
                                  f"section RUN")
        if not os.path.isdir(output_path):
            raise ValueError(_e + f"The specified input path is invalid: {input_path}")
        self.input_path = input_path

        bam_ending = run_config.get("RUN", "bamEnding", fallback=run_config.get("RUN", "ending", fallback=None))
        if bam_ending is None:
            raise ValueError(_e + f"Missing mandatory parameter \"bamEnding\" in section RUN")
        self.bam_ending = bam_ending

        # Get BAM file names in specified input folder for validation and auto-generation of [SAMPLES] section
        _ending_length = len(bam_ending)
        bam_file_absolute_paths = [os.path.join(
            input_path, p
        ) for p in os.listdir(input_path) if p[-_ending_length:] == bam_ending]
        sample_names = [
            (os.path.basename(p)).replace(bam_ending, "") for p in bam_file_absolute_paths
        ]
        if len(sample_names) == 0:
            raise ValueError(_e + f"No BAM files with the specified ending ({bam_ending}) were found in the specified "
                                  f"location: {input_path}")

        species = run_config.get("RUN", "species", fallback=None)
        if species is None:
            raise ValueError(_e + f"Missing mandatory parameter \"species\" in section RUN")
        if species not in internal_config.allowed_species:
            raise ValueError(_e + f"An invalid species ({species}) was specified. "
                                  f"Allowed species: {internal_config.allowed_species}")
        self.species = species

        genome = run_config.get(
            "RUN", "genome", fallback=run_config.get(
                "RUN", "referenceGenome", fallback=None
            )
        )
        if genome is None:
            raise ValueError(_e + f"Missing mandatory parameter \"genome\" (path to reference genome GTF file) in "
                                  f"section RUN")
        self.genome = genome

        # [RUN] - Optional parameters

        n_threads = run_config.get("RUN", "threads", fallback=1)  # Default to one thread
        self.n_threads = int(n_threads)

        max_n_features_in_transcript = run_config.get(
            "RUN", "maxNumberFeaturesInTranscript", fallback=run_config.get(
                "RUN", "maxNFeaturesInTranscript", fallback=run_config.get(
                    "RUN", "maxFeaturesInTranscript", fallback=0  # Default to 0 (interpreted as no limit)
                )
            )
        )
        self.max_n_features_in_transcript = int(max_n_features_in_transcript)

        generate_report = run_config.get("RUN", "report", fallback=internal_config.default_generate_report)
        self.generate_report = parse_bool(generate_report)

        mapq_for_unique_mapping = run_config.get(
            "RUN", "mapqForUniqueMapping", fallback=internal_config.default_mapq_for_unique_mapping
        )
        self.mapq_for_unique_mapping = int(mapq_for_unique_mapping)

        include_all_junctions_in_output = run_config.get(
            "RUN", "include_all_junctions_in_output", fallback=internal_config.default_include_all_junctions_in_output
        )
        self.include_all_junctions_in_output = parse_bool(include_all_junctions_in_output)

        only_use_longest_annotated_transcript = run_config.get(
            "RUN", "onlyUseLongestAnnotatedTranscript",
            fallback=internal_config.default_only_use_longest_annotated_transcript
        )
        self.only_use_longest_annotated_transcript = parse_bool(only_use_longest_annotated_transcript)

        skip_transcripts_with_redundant_feature_annotation = run_config.get(
            "RUN", "skipTranscriptsWithRedundantFeatureAnnotation",
            fallback=internal_config.default_skip_transcripts_with_redundant_feature_annotation
        )
        self.skip_transcripts_with_redundant_feature_annotation = parse_bool(
            skip_transcripts_with_redundant_feature_annotation
        )

        # Optional section: [REPORT]

        report_max_plotted = run_config.get(
            "REPORT", "maxPlotted", fallback=run_config.get(
                "REPORT", "maxNPlotted", fallback=None  # Default to None (no limit)
            )
        )
        self.report_max_plotted = None if report_max_plotted is None else int(report_max_plotted)

        report_min_total_number_occurrences_across_all_samples = run_config.get(
            "REPORT", "minTotalNumberOccurrencesAcrossAllSamples", fallback=run_config.get(
                "REPORT", "minTotalNOccurrencesAcrossAllSamples", fallback=run_config.get(
                    "REPORT", "minTotalOccurrencesAcrossAllSamples", fallback=1  # Default to 1 (excludes no activity)
                )
            )
        )
        self.report_min_total_number_occurrences_across_all_samples = int(
            report_min_total_number_occurrences_across_all_samples
        )

        report_min_n_occurrences_in_sample = run_config.get(
            "REPORT", "minNumberOccurrencesInSample", fallback=run_config.get(
                "REPORT", "minNOccurrencesInSample", fallback=run_config.get(
                    "REPORT", "minOccurrencesInSample", fallback=None  # Check for paired use, then default to 0
                )
            )
        )

        report_occurrences_in_at_least_n_samples = run_config.get(
            "REPORT", "occurrencesInAtLeastNSamples", fallback=run_config.get(
                "REPORT", "inAtLeastNSamples", fallback=None  # Check for paired use, then default to 0
            )
        )

        # Disallow using either report_min_n_occurrences_in_sample or report_occurrences_in_at_least_n_samples
        # alone, without both being set
        if any([
            report_min_n_occurrences_in_sample is None and report_occurrences_in_at_least_n_samples is not None,
            report_min_n_occurrences_in_sample is not None and report_occurrences_in_at_least_n_samples is None
        ]):
            raise ValueError(_e + f"Parameters minNumberOccurrencesInSample and occurrencesInAtLeastNSamples must be "
                                  f"used together. Either define both, or exclude both")

        self.report_min_n_occurrences_in_sample = 0 if report_min_n_occurrences_in_sample is None else \
            int(report_min_n_occurrences_in_sample)
        if self.report_min_n_occurrences_in_sample < 0:
            raise ValueError(_e + f"The parameter minNumberOccurrencesInSample cannot be less than zero")

        if report_occurrences_in_at_least_n_samples in ("all", "All", "ALL", "*"):
            # Special case: Set value to total number of samples
            self.report_occurrences_in_at_least_n_samples = len(sample_names)
        else:
            self.report_occurrences_in_at_least_n_samples = 0 if report_occurrences_in_at_least_n_samples is None else \
                int(report_occurrences_in_at_least_n_samples)
        if self.report_occurrences_in_at_least_n_samples < 0:
            raise ValueError(_e + f"The parameter occurrencesInAtLeastNSamples cannot be less than zero")

        if self.generate_report:

            # Auto-generated if required: [SAMPLES]

            if "SAMPLES" in run_config:
                sample_groups = {
                    key: value.split(" ")
                    for key, value in run_config.items("SAMPLES") if value
                }
                _flattened_sample_groups_values = flatten_nested_lists(list(sample_groups.values()))
                if not all([
                    sample_name in sample_names
                    for sample_name in _flattened_sample_groups_values
                ]):
                    raise ValueError(_e + "Sample names in run config do not match sample names parsed from BAM files")
                if not all([
                    sample_name in _flattened_sample_groups_values
                    for sample_name in sample_names
                ]):
                    raise ValueError(_e + f"Sample names parsed from BAM files are not all present in the SAMPLES "
                                          f"section of run config. Each sample name must be assigned to a group.")
            else:
                sample_groups = {"all": sample_names}
            self.sample_groups = sample_groups

            group_name_by_sample_name = {}
            for group_name, sample_names in sample_groups.items():
                for sample_name in sample_names:
                    group_name_by_sample_name[sample_name] = group_name
            self.group_name_by_sample_name = group_name_by_sample_name

            # Auto-generated if required: [COLORS]

            # Allow "COLOURS" as alternative spelling
            if "COLOURS" in run_config:
                colors_alt_spelling_section_contents = run_config.items("COLOURS")
                run_config.add_section("COLORS")
                for key, value in colors_alt_spelling_section_contents:
                    run_config.set("COLORS", key, value)
                run_config.remove_section("COLOURS")

            # Use defined values if both COLORS and SAMPLES sections are present, otherwise auto-generate
            if "COLORS" in run_config and "SAMPLES" in run_config:
                color_by_group_name = dict(run_config.items("COLORS"))
            else:
                color_by_group_name = {
                    group_name: internal_config.default_point_color_in_plot
                    for group_name in sample_groups.keys()
                }
            self.color_by_group_name = color_by_group_name

        else:

            self.sample_groups = None
            self.group_name_by_sample_name = None
            self.color_by_group_name = None


# FOR TESTING - run as standalone
if __name__ == "__main__":

    internal_config_relative_path = "internal.config"

    frc = FaseRunConfig(
        "/Users/aasho2/Projects/FASE_V1/fase_run.config",
        FaseInternalConfig(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), internal_config_relative_path)
        )
    )

    print("\n".join([f"{a}: {type(getattr(frc, a))} {getattr(frc, a)}" for a in dir(frc) if a[:2] != "__"]))
