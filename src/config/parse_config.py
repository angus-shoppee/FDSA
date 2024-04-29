
from typing import Union, List, Dict, Any
from functools import reduce
import csv
import os
import configparser

from src.config.marker_aliases import convert_marker_alias


RANK_RESULTS_BY_ALLOWED_VALUES = ["frequency", "number"]


def _remove_lowercase_duplicates(values: List[str]) -> List[str]:

    # Ugly function to help bypass configparser duplicating options and their corresponding values when case
    # Sensitivity is disabled (desired behaviour here is to retain only the mixed case option if applicable)

    all_lowercase = [value for value in values if value == value.lower()]
    lowercase_with_pair = [lowercase_value for lowercase_value in all_lowercase if len([
        value for value in values if value.lower() == lowercase_value
    ]) > 1]
    modified_values = values.copy()
    for value in lowercase_with_pair:
        modified_values.remove(value)
    return modified_values


def remove_quotes(value: str) -> str:
    return value.replace("\"", "").replace("\'", "")


def flatten_nested_lists(nested_lists: List[List[Any]]) -> List[Any]:
    return reduce(lambda a, b: a + b, nested_lists)


def parse_bool(value: Union[bool, str]) -> bool:

    if isinstance(value, bool):
        return value

    # Allow the user to use Yes/No for boolean flags if they prefer
    if value.lower() == "true" or value.lower() == "yes":
        return True
    elif value.lower() == "false" or value.lower() == "no":
        return False

    raise ValueError(f"Could not parse bool from input value: {value}. Allowed values are True/False or Yes/No "
                     f"(not case-sensitive).")


# TODO: Refactor names for default report parameters to be consistently named between internal and user config classes


class FaseInternalConfig:
    allowed_species: List[str]
    biomart_name_for_species: Dict[str, str]
    default_biomart_mirror: str
    default_mapq_for_unique_mapping: int
    default_rank_results_by: str
    default_feature_junction_overlap_threshold: float
    default_max_n_features_in_transcript: int
    default_only_use_longest_annotated_transcript: bool
    default_skip_transcripts_with_redundant_feature_annotation: bool
    # ^ skip_transcripts_with_redundant_feature_annotation is overridden by only_use_longest_annotated_transcript=True
    default_include_all_junctions_in_output: bool
    default_paired_end_reads: bool
    default_generate_report: bool
    default_save_normalized_gene_counts: bool
    default_draw_junctions_min_count: int
    default_point_color_in_plot: str
    default_point_shape_in_plot: str
    default_feature_counts_primary_alignment_only: bool
    # Possible TODO: Generalise default_feature_counts_primary_alignment_only to default_primary_alignment_only

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
        self.allowed_species = [s.lower() for s in allowed_species_string.split(" ")]

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

        default_biomart_mirror = config.get("BUILD", "defaultBiomartMirror", fallback=None)
        if default_biomart_mirror is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultBiomartMirror\" in section BUILD")
        self.default_biomart_mirror = default_biomart_mirror

        # [RUN]

        default_mapq_for_unique_mapping = config.get(
            "RUN", "defaultMapqForUniqueMapping", fallback=None
        )
        if default_mapq_for_unique_mapping is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultMapqForUniqueMapping\" in section RUN")
        self.default_mapq_for_unique_mapping = int(default_mapq_for_unique_mapping)

        default_rank_results_by = config.get("REPORT", "defaultRankResultsBy", fallback=None)
        if default_rank_results_by is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultRankResultsBy\" in section REPORT")
        if default_rank_results_by not in RANK_RESULTS_BY_ALLOWED_VALUES:
            raise ValueError(_e + f"Invalid value \"{default_rank_results_by}\" specified for parameter "
                                  f"\"defaultRankResultsBy\" in section RUN. Allowed options are: "
                                  f"{RANK_RESULTS_BY_ALLOWED_VALUES}")
        self.default_rank_results_by = default_rank_results_by

        default_feature_junction_overlap_threshold = config.get(
            "RUN", "defaultFeatureJunctionOverlapThreshold", fallback=None
        )
        if default_feature_junction_overlap_threshold is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultFeatureJunctionOverlapThreshold\" in section "
                                  f"RUN")
        self.default_feature_junction_overlap_threshold = float(default_feature_junction_overlap_threshold)

        default_max_n_features_in_transcript = config.get(
            "RUN", "defaultMaxNFeaturesInTranscript", fallback=None
        )
        if default_max_n_features_in_transcript is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultMaxNFeaturesInTranscript\" in section RUN")
        self.default_max_n_features_per_transcript = int(default_max_n_features_in_transcript)

        default_only_use_longest_annotated_transcript = config.get(
            "RUN", "defaultOnlyUseLongestAnnotatedTranscript", fallback=None
        )
        if default_only_use_longest_annotated_transcript is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultOnlyUseLongestAnnotatedTranscript\" in "
                                  f"section RUN")
        self.default_only_use_longest_annotated_transcript = parse_bool(
            default_only_use_longest_annotated_transcript
        )

        default_skip_transcripts_with_redundant_feature_annotation = config.get(
            "RUN", "defaultSkipTranscriptsWithRedundantFeatureAnnotation", fallback=None
        )
        if default_skip_transcripts_with_redundant_feature_annotation is None:
            raise ValueError(_e + f"Missing mandatory parameter "
                                  f"\"defaultSkipTranscriptsWithRedundantFeatureAnnotation\" in section RUN")
        self.default_skip_transcripts_with_redundant_feature_annotation = parse_bool(
            default_skip_transcripts_with_redundant_feature_annotation
        )

        default_include_all_junctions_in_output = config.get(
            "RUN", "defaultIncludeAllJunctionsInOutput", fallback=None
        )
        if default_include_all_junctions_in_output is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultIncludeAllJunctionsInOutput\" in section "
                                  f"RUN")
        self.default_include_all_junctions_in_output = parse_bool(default_include_all_junctions_in_output)

        default_paired_end_reads = config.get(
            "RUN", "defaultPairedEndReads", fallback=None
        )
        if default_paired_end_reads is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultPairedEndReads\" in section RUN")
        self.default_paired_end_reads = parse_bool(default_paired_end_reads)

        # [REPORT]

        default_generate_report = config.get(
            "REPORT", "defaultGenerateReport", fallback=None
        )
        if default_generate_report is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultGenerateReport\" in section REPORT")
        self.default_generate_report = parse_bool(default_generate_report)

        default_save_normalized_gene_counts = config.get(
            "REPORT", "defaultSaveNormalizedGeneCounts", fallback=None
        )
        if default_save_normalized_gene_counts is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultSaveNormalizedGeneCounts\" in section REPORT")
        self.default_save_normalized_gene_counts = parse_bool(default_save_normalized_gene_counts)

        default_point_color_in_plot = config.get(
            "REPORT", "defaultPointColorInPlot", fallback=None
        )
        if default_point_color_in_plot is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultPointColorInPlot\" in section REPORT")
        self.default_point_color_in_plot = default_point_color_in_plot

        default_point_shape_in_plot = config.get(
            "REPORT", "defaultPointShapeInPlot", fallback=None
        )
        if default_point_shape_in_plot is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultPointShapeInPlot\" in section REPORT")
        self.default_point_shape_in_plot = convert_marker_alias(default_point_shape_in_plot)

        default_feature_counts_primary_alignment_only = config.get(
            "REPORT", "defaultFeatureCountsPrimaryAlignmentOnly", fallback=None
        )
        if default_feature_counts_primary_alignment_only is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultFeatureCountsPrimaryAlignmentOnly\" in "
                                  f"section REPORT")
        self.default_feature_counts_primary_alignment_only = parse_bool(
            default_feature_counts_primary_alignment_only
        )

        default_draw_junctions_min_count = config.get(
            "REPORT", "defaultDrawJunctionsMinCount", fallback=None
        )
        if default_draw_junctions_min_count is None:
            raise ValueError(_e + f"Missing mandatory parameter \"defaultDrawJunctionsMinCount\" in section REPORT")
        self.default_draw_junctions_min_count = int(default_draw_junctions_min_count)

    def __setattr__(self, name, value):
        # If attribute is a string, strip all quotation marks before applying
        if isinstance(value, str):
            value = remove_quotes(value)
        super().__setattr__(name, value)


class FaseUserConfig:
    email: str
    user_default_biomart_mirror: Union[None, str]
    user_default_generate_report: Union[None, bool]
    user_default_n_threads: Union[None, int]
    user_default_max_n_features_in_transcript: Union[None, int]
    user_default_only_use_longest_annotated_transcript: Union[None, int]
    user_default_skip_transcripts_with_redundant_feature_annotation: Union[None, int]
    user_default_feature_junction_overlap_threshold: Union[None, float]
    user_default_primary_alignment_only: Union[None, bool]
    user_default_mapq_for_unique_mapping: Union[None, int]
    user_default_include_all_junctions_in_output: Union[None, bool]
    user_default_report_feature_counts_executable: Union[None, str]
    user_default_report_save_normalized_gene_counts: Union[None, bool]
    user_default_report_n_threads = Union[None, int]
    user_default_report_rank_results_by: Union[None, str]
    user_default_report_max_n_plotted: Union[None, int]
    user_default_report_min_total_n_occurrences_across_all_samples: Union[None, int]  # Chef's kiss variable name
    user_default_report_min_n_occurrences_in_sample: Union[None, int]
    user_default_report_occurrences_in_at_least_n_samples: Union[None, int, str]
    user_default_report_draw_junctions_min_count: Union[None, int]
    # ^ If * is set as the default for this parameter, the final number will be set at runtime when n samples is known

    def __init__(
        self,
        user_config_path: str,
    ) -> None:

        _e = f"Error while parsing user config at {user_config_path}: "

        user_config = configparser.ConfigParser()
        user_config.read(user_config_path)

        # [USER]

        email = user_config.get(
            "BUILD", "emailAddress", fallback=user_config.get(
                "BUILD", "email", fallback=None
            )
        )
        if email is None:
            raise ValueError(_e + f"Missing mandatory parameter \"email\" in section BUILD. Please note that this "
                                  f"email address is used only for external academic API calls during the build "
                                  f"process.")
        # Possible TODO: Validity check on provided email address (or later in main.build)?
        self.email = email

        user_default_biomart_mirror = user_config.get(
            "BUILD", "biomartMirror", fallback=user_config.get(
                "BUILD", "biomart", fallback=None
            )
        )
        self.user_default_biomart_mirror = user_default_biomart_mirror

        # [DEFAULT RUN]

        user_default_generate_report = user_config.get("DEFAULT RUN", "report", fallback=None)
        self.user_default_generate_report = None if user_default_generate_report is None else \
            parse_bool(user_default_generate_report)

        user_default_n_threads = user_config.get("DEFAULT RUN", "threads", fallback=None)
        self.user_default_n_threads = None if user_default_n_threads is None else int(user_default_n_threads)

        user_default_max_n_features_in_transcript = user_config.get(
            "RUN", "maxNumberFeaturesInTranscript", fallback=user_config.get(
                "RUN", "maxNFeaturesInTranscript", fallback=user_config.get(
                    "RUN", "maxFeaturesInTranscript", fallback=None
                )
            )
        )
        self.user_default_max_n_features_in_transcript = None \
            if user_default_max_n_features_in_transcript is None \
            else int(user_default_max_n_features_in_transcript)

        user_default_only_use_longest_annotated_transcript = user_config.get(
            "DEFAULT RUN", "onlyUseLongestAnnotatedTranscript", fallback=None
        )
        self.user_default_only_use_longest_annotated_transcript = None \
            if user_default_only_use_longest_annotated_transcript is None \
            else parse_bool(user_default_only_use_longest_annotated_transcript)

        user_default_skip_transcripts_with_redundant_feature_annotation = user_config.get(
            "DEFAULT RUN", "skipTranscriptsWithRedundantFeatureAnnotation", fallback=None
        )
        self.user_default_skip_transcripts_with_redundant_feature_annotation = None \
            if user_default_skip_transcripts_with_redundant_feature_annotation is None \
            else parse_bool(user_default_skip_transcripts_with_redundant_feature_annotation)

        user_default_feature_junction_overlap_threshold = user_config.get(
            "RUN", "featureJunctionOverlapThreshold", fallback=user_config.get(
                "RUN", "overlapThreshold", fallback=None
            )
        )
        self.user_default_feature_junction_overlap_threshold = None \
            if user_default_feature_junction_overlap_threshold is None \
            else float(user_default_feature_junction_overlap_threshold)

        user_default_primary_alignment_only = user_config.get(
            "DEFAULT RUN", "primaryAlignmentOnly", fallback=user_config.get(
                "DEFAULT RUN", "primaryAlignment", fallback=None
            )
        )
        self.user_default_primary_alignment_only = None if user_default_primary_alignment_only is None \
            else parse_bool(user_default_primary_alignment_only)

        user_default_mapq_for_unique_mapping = user_config.get("DEFAULT RUN", "mapqForUniqueMapping", fallback=None)
        self.user_default_mapq_for_unique_mapping = None if user_default_mapq_for_unique_mapping is None \
            else int(user_default_mapq_for_unique_mapping)

        user_default_include_all_junctions_in_output = user_config.get(
            "DEFAULT RUN", "includeAllJunctionsInOutput", fallback=None
        )
        self.user_default_include_all_junctions_in_output = None \
            if user_default_include_all_junctions_in_output is None \
            else parse_bool(user_default_include_all_junctions_in_output)

        # [DEFAULT REPORT]

        user_default_report_feature_counts_executable = user_config.get(
            "DEFAULT REPORT", "featureCountsExecutableLocation", fallback=user_config.get(
                "DEFAULT REPORT", "featureCountsExecutablePath", fallback=user_config.get(
                    "DEFAULT REPORT", "featureCountsExecutable", fallback=user_config.get(
                        "DEFAULT REPORT", "featureCounts", fallback=None
                    )
                )
            )
        )
        self.user_default_report_feature_counts_executable = user_default_report_feature_counts_executable

        user_default_report_save_normalized_gene_counts = user_config.get(
            "DEFAULT REPORT", "saveNormalizedGeneCounts", fallback=user_config.get(
                "DEFAULT REPORT", "saveNormGeneCounts", fallback=None
            )
        )
        self.user_default_report_save_normalized_gene_counts = None \
            if user_default_report_save_normalized_gene_counts is None \
            else parse_bool(user_default_report_save_normalized_gene_counts)

        user_default_report_n_threads = user_config.get(
            "DEFAULT REPORT", "threads", fallback=self.user_default_n_threads
        )
        self.user_default_report_n_threads = user_default_report_n_threads

        user_default_report_rank_results_by = user_config.get(
            "DEFAULT REPORT", "rankResultsBy", fallback=user_config.get(
                "DEFAULT REPORT", "rankBy", fallback=None
            )
        )
        self.user_default_report_rank_results_by = None if user_default_report_rank_results_by is None \
            else user_default_report_rank_results_by
        if (
            self.user_default_report_rank_results_by is not None
            and self.user_default_report_rank_results_by not in RANK_RESULTS_BY_ALLOWED_VALUES
        ):
            raise ValueError(_e + f"Invalid value \"{user_default_report_rank_results_by}\" specified for parameter "
                                  f"\"rankResultsBy\" in section REPORT DEFAULT. Allowed options are: "
                                  f"{RANK_RESULTS_BY_ALLOWED_VALUES}")

        user_default_report_max_n_plotted = user_config.get(
            "DEFAULT REPORT", "maxNumberPlotted", fallback=user_config.get(
                "DEFAULT REPORT", "maxNPlotted", fallback=user_config.get(
                    "DEFAULT REPORT", "maxPlotted", fallback=None
                )
            )
        )
        self.user_default_report_max_n_plotted = None if user_default_report_max_n_plotted is None \
            else int(user_default_report_max_n_plotted)

        user_default_report_min_total_n_occurrences_across_all_samples = user_config.get(
            "DEFAULT REPORT", "minTotalNumberOccurrencesAcrossAllSamples", fallback=user_config.get(
                "DEFAULT REPORT", "minTotalNOccurrencesAcrossAllSamples", fallback=user_config.get(
                    "DEFAULT REPORT", "minTotalOccurrencesAcrossAllSamples", fallback=None
                )
            )
        )
        self.user_default_report_min_total_n_occurrences_across_all_samples = None \
            if user_default_report_min_total_n_occurrences_across_all_samples is None \
            else int(user_default_report_min_total_n_occurrences_across_all_samples)

        user_default_report_min_n_occurrences_in_sample = user_config.get(
            "DEFAULT REPORT", "minNumberOccurrencesInSample", fallback=user_config.get(
                "DEFAULT REPORT", "minNOccurrencesInSample", fallback=user_config.get(
                    "DEFAULT REPORT", "minOccurrencesInSample", fallback=None
                )
            )
        )

        user_default_report_occurrences_in_at_least_n_samples = user_config.get(
            "DEFAULT REPORT", "occurrencesInAtLeastNSamples", fallback=user_config.get(
                "DEFAULT REPORT", "inAtLeastNSamples", fallback=None
            )
        )

        if any([
            user_default_report_min_n_occurrences_in_sample is None
            and user_default_report_occurrences_in_at_least_n_samples is not None,
            user_default_report_min_n_occurrences_in_sample is not None
            and user_default_report_occurrences_in_at_least_n_samples is None
        ]):
            raise ValueError(_e + f"Parameters minNumberOccurrencesInSample and occurrencesInAtLeastNSamples must be "
                                  f"used together. Either define both, or exclude both")

        self.user_default_report_min_n_occurrences_in_sample = 0 \
            if user_default_report_min_n_occurrences_in_sample is None \
            else int(user_default_report_min_n_occurrences_in_sample)
        if self.user_default_report_min_n_occurrences_in_sample < 0:
            raise ValueError(_e + f"The parameter minNumberOccurrencesInSample cannot be less than zero")

        if user_default_report_occurrences_in_at_least_n_samples in ("all", "All", "ALL", "*"):
            # Special case: Set value to total number of samples
            self.user_default_report_occurrences_in_at_least_n_samples = "*"
        else:
            self.user_default_report_occurrences_in_at_least_n_samples = 0 \
                if user_default_report_occurrences_in_at_least_n_samples is None \
                else int(user_default_report_occurrences_in_at_least_n_samples)
        if self.user_default_report_occurrences_in_at_least_n_samples < 0:
            raise ValueError(_e + f"The parameter occurrencesInAtLeastNSamples cannot be less than zero")

        user_default_report_draw_junctions_min_count = user_config.get(
            "DEFAULT REPORT", "drawJunctionsMinCount", fallback=user_config.get(
                "DEFAULT REPORT", "drawJunctionsMinN", fallback=user_config.get(
                    "DEFAULT REPORT", "drawJunctionsMin", fallback=None
                )
            )
        )
        self.user_default_report_draw_junctions_min_count = None \
            if user_default_report_draw_junctions_min_count is None \
            else int(user_default_report_draw_junctions_min_count)

    def __setattr__(self, name, value):
        # If attribute is a string, strip all quotation marks before applying
        if isinstance(value, str):
            value = remove_quotes(value)
        super().__setattr__(name, value)


class FaseRunConfig:
    run_name: str
    feature_name: str
    output_path: str
    input_path: str
    bam_ending: str
    paired_end_reads: bool
    species: str
    genome: str
    genes: Union[None, Dict[str, bool]]  # Check .genes.get(gene_name, False) for speed improvement over list membership
    n_threads: int
    rank_results_by: str
    max_n_features_in_transcript: int
    only_use_longest_annotated_transcript: bool
    skip_transcripts_with_redundant_feature_annotation: bool  # overridden by only_use_longest_annotated_transcript=True
    feature_junction_overlap_threshold: float
    generate_report: bool
    report_name: str
    report_genes: Union[None, Dict[str, bool]]  # As with .genes, check .report_genes.get(gene_name, False)
    report_n_threads: int
    report_feature_counts_executable: Union[None, str]
    report_save_normalized_gene_counts: bool
    report_gene_count_matrix: Union[None, str]
    report_max_n_plotted: Union[None, int]
    report_min_total_n_occurrences_across_all_samples: int
    report_min_n_occurrences_in_sample: int
    report_occurrences_in_at_least_n_samples: int
    report_draw_junctions_min_count: int
    report_transcript_plot_max_samples_per_group: Union[None, int]
    primary_alignment_only: bool
    mapq_for_unique_mapping: int
    sample_groups: Dict[str, List[str]]
    group_name_by_sample_name: Dict[str, str]  # Interpolated, not defined by the user
    color_by_group_name: Dict[str, str]
    shape_by_group_name: Dict[str, str]

    def __init__(
        self,
        run_config_path: str,
        internal_config: FaseInternalConfig,
        user_config: Union[None, FaseUserConfig]
    ) -> None:

        if not os.path.exists(run_config_path):
            raise ValueError(f"The specified file does not exist: {run_config_path}")

        _e = f"Error while parsing run config at {run_config_path}: "

        run_config = configparser.ConfigParser()
        run_config.read(run_config_path)

        if user_config is None:
            user_config = FaseUserConfig(os.path.join("src", "config", ".null_user.config"))

        if "RUN" not in run_config:
            raise ValueError(_e + f"Missing mandatory section RUN (run settings)")

        # [RUN] - Mandatory parameters

        run_name = run_config.get(
            "RUN", "runName", fallback=run_config.get(
                "RUN", "name", fallback=None
            )
        )
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
        self.input_path = input_path
        if not os.path.isdir(self.input_path):
            raise ValueError(_e + f"The specified input path is invalid: {self.input_path}")

        bam_ending = run_config.get("RUN", "bamEnding", fallback=run_config.get("RUN", "ending", fallback=None))
        if bam_ending is None:
            raise ValueError(_e + f"Missing mandatory parameter \"bamEnding\" in section RUN")
        self.bam_ending = bam_ending

        # Get BAM file names in specified input folder for validation and auto-generation of [SAMPLES] section
        _ending_length = len(self.bam_ending)
        bam_file_absolute_paths = [os.path.join(
            self.input_path, p
        ) for p in os.listdir(self.input_path) if p[-_ending_length:] == self.bam_ending]
        sample_names = [
            (os.path.basename(p)).replace(self.bam_ending, "") for p in bam_file_absolute_paths
        ]
        if len(sample_names) == 0:
            raise ValueError(_e + f"No BAM files with the specified ending ({self.bam_ending}) were found in the "
                                  f"specified location: {self.input_path}")

        # Allow this parameter to be set in REPORT section rather than RUN (since it is technically for report only)
        paired_end_reads = run_config.get(
            "RUN", "pairedEndReads", fallback=run_config.get(
                "RUN", "pairedEnd", fallback=run_config.get(
                    "REPORT", "pairedEndReads", fallback=run_config.get(
                        "REPORT", "pairedEnd", fallback=internal_config.default_paired_end_reads
                    )
                )
            )
        )
        self.paired_end_reads = parse_bool(paired_end_reads)

        species = run_config.get("RUN", "species", fallback=None)
        if species is None:
            raise ValueError(_e + f"Missing mandatory parameter \"species\" in section RUN")
        species = species.lower()
        self.species = species
        if self.species not in internal_config.allowed_species:
            raise ValueError(_e + f"An invalid species ({self.species}) was specified. "
                                  f"Allowed species: {internal_config.allowed_species}")

        genome = run_config.get(
            "RUN", "genome", fallback=run_config.get(
                "RUN", "referenceGenome", fallback=None
            )
        )
        if genome is None:
            raise ValueError(_e + f"Missing mandatory parameter \"genome\" (path to reference genome GTF file) in "
                                  f"section RUN")
        self.genome = genome

        genes = run_config.get(
            "RUN", "geneList", fallback=run_config.get(
                "RUN", "geneSet", fallback=run_config.get(
                    "RUN", "genes", fallback=None
                )
            )
        )
        # If gene names are supplied (either inline or in a single-column, no-header CSV file), map gene name in
        # lowercase to True in self.genes for fast checking.
        # No case sensitivity to avoid issues with lack of consistency in capitalisation within names of human
        # ORF genes, which are often spelled in partial lowercase
        if genes is None:
            self.genes = None
        elif ".csv" in genes.lower():
            # Handle input given as path to a CSV file
            if not os.path.exists(genes):
                raise ValueError(
                    "A CSV file was specified for the parameter genes in section RUN, but its path is not valid: " +
                    genes
                )
            with open(genes, "r") as gene_set_file:
                self.genes = {str(row[0]).lower(): True for row in csv.reader(gene_set_file)}
        else:
            # Allow inline gene set to be delimited by either spaces or commas
            self.genes = {remove_quotes(gene_name).lower(): True for gene_name in genes.replace(",", " ").split()}

        # [RUN] - Optional parameters

        n_threads = run_config.get("RUN", "threads", fallback=(
            user_config.user_default_n_threads
            if user_config.user_default_n_threads is not None
            else 1
        ))  # Default to one thread
        self.n_threads = int(n_threads)

        max_n_features_in_transcript = run_config.get(
            "RUN", "maxNumberFeaturesInTranscript", fallback=run_config.get(
                "RUN", "maxNFeaturesInTranscript", fallback=run_config.get(
                    "RUN", "maxFeaturesInTranscript", fallback=(
                        user_config.user_default_max_n_features_in_transcript
                        if user_config.user_default_max_n_features_in_transcript is not None
                        else 0  # Default to 0 (interpreted as no limit)
                    )
                )
            )
        )
        if max_n_features_in_transcript in ("*", "ANY", "Any", "any", "ALL", "All", "all"):
            self.max_n_features_in_transcript = 0
        else:
            self.max_n_features_in_transcript = int(max_n_features_in_transcript)

        generate_report = run_config.get("RUN", "report", fallback=(
            user_config.user_default_generate_report
            if user_config.user_default_generate_report is not None
            else internal_config.default_generate_report
        ))
        self.generate_report = parse_bool(generate_report)

        primary_alignment_only = run_config.get(
            "RUN", "primaryAlignmentOnly", fallback=run_config.get(
                "RUN", "primaryAlignment", fallback=(
                    user_config.user_default_primary_alignment_only
                    if user_config.user_default_primary_alignment_only is not None
                    else internal_config.default_feature_counts_primary_alignment_only
                )
            )
        )
        self.primary_alignment_only = primary_alignment_only

        mapq_for_unique_mapping = run_config.get(
            "RUN", "mapqForUniqueMapping", fallback=(
                user_config.user_default_mapq_for_unique_mapping
                if user_config.user_default_mapq_for_unique_mapping is not None
                else internal_config.default_mapq_for_unique_mapping
            )
        )
        self.mapq_for_unique_mapping = int(mapq_for_unique_mapping)

        include_all_junctions_in_output = run_config.get(
            "RUN", "includeAllJunctionsInOutput", fallback=(
                user_config.user_default_include_all_junctions_in_output
                if user_config.user_default_include_all_junctions_in_output is not None
                else internal_config.default_include_all_junctions_in_output
            )
        )
        self.include_all_junctions_in_output = parse_bool(include_all_junctions_in_output)

        only_use_longest_annotated_transcript = run_config.get(
            "RUN", "onlyUseLongestAnnotatedTranscript",
            fallback=(
                user_config.user_default_only_use_longest_annotated_transcript
                if user_config.user_default_only_use_longest_annotated_transcript is not None
                else internal_config.default_only_use_longest_annotated_transcript
            )
        )
        self.only_use_longest_annotated_transcript = parse_bool(only_use_longest_annotated_transcript)

        skip_transcripts_with_redundant_feature_annotation = run_config.get(
            "RUN", "skipTranscriptsWithRedundantFeatureAnnotation",
            fallback=(
                user_config.user_default_skip_transcripts_with_redundant_feature_annotation
                if user_config.user_default_skip_transcripts_with_redundant_feature_annotation is not None
                else internal_config.default_skip_transcripts_with_redundant_feature_annotation
            )
        )
        self.skip_transcripts_with_redundant_feature_annotation = parse_bool(
            skip_transcripts_with_redundant_feature_annotation
        )

        feature_junction_overlap_threshold = run_config.get(
            "RUN", "featureJunctionOverlapThreshold", fallback=run_config.get(
                "RUN", "overlapThreshold", fallback=(
                    user_config.user_default_feature_junction_overlap_threshold
                    if user_config.user_default_feature_junction_overlap_threshold is not None
                    else internal_config.default_feature_junction_overlap_threshold
                )
            )
        )
        self.feature_junction_overlap_threshold = float(feature_junction_overlap_threshold)

        # Optional section: [REPORT]

        report_name = run_config.get(
            "REPORT", "reportName", fallback=run_config.get(
                "REPORT", "name", fallback=self.run_name  # Use run name if no report name specified
            )
        )
        self.report_name = report_name

        report_genes = run_config.get(
            "REPORT", "geneList", fallback=run_config.get(
                "REPORT", "geneSet", fallback=run_config.get(
                    "REPORT", "genes", fallback=None
                )
            )
        )
        # Same process as for .genes
        if report_genes is None:
            self.report_genes = None
        elif ".csv" in report_genes.lower():
            # Handle input given as path to a CSV file
            if not os.path.exists(report_genes):
                raise ValueError(
                    "A CSV file was specified for the parameter genes in section RUN, but its path is not valid: " +
                    report_genes
                )
            with open(report_genes, "r") as gene_set_file:
                self.report_genes = {str(row[0]).lower(): True for row in csv.reader(gene_set_file)}
        else:
            # Allow inline gene set to be delimited by either spaces or commas
            self.report_genes = {
                remove_quotes(gene_name).lower(): True
                for gene_name in report_genes.replace(",", " ").split()
            }

        report_n_threads = run_config.get("REPORT", "threads", fallback=user_config.user_default_report_n_threads)
        self.report_n_threads = report_n_threads

        report_feature_counts_executable = run_config.get(
            "REPORT", "featureCountsExecutableLocation", fallback=run_config.get(
                "REPORT", "featureCountsExecutablePath", fallback=run_config.get(
                    "REPORT", "featureCountsExecutable", fallback=run_config.get(
                        "REPORT", "featureCounts", fallback=user_config.user_default_report_feature_counts_executable
                    )
                )
            )
        )
        self.report_feature_counts_executable = report_feature_counts_executable

        report_save_normalized_gene_counts = run_config.get(
            "REPORT", "saveNormalizedGeneCounts", fallback=run_config.get(
                "REPORT", "saveNormGeneCounts", fallback=(
                    user_config.user_default_report_save_normalized_gene_counts
                    if user_config.user_default_report_save_normalized_gene_counts is not None
                    else internal_config.default_save_normalized_gene_counts
                )
            )
        )
        self.report_save_normalized_gene_counts = parse_bool(report_save_normalized_gene_counts)

        report_gene_count_matrix = run_config.get(
            "REPORT", "featureCountsOutput", fallback=run_config.get(
                "REPORT", "geneCountsMatrix", fallback=run_config.get(
                    "REPORT", "geneCounts", fallback=None
                )
            )
        )
        self.report_gene_count_matrix = report_gene_count_matrix

        rank_results_by = run_config.get(
            "REPORT", "rankResultsBy", fallback=run_config.get(
                "REPORT", "rankBy", fallback=(
                    user_config.user_default_report_rank_results_by
                    if user_config.user_default_report_rank_results_by is not None
                    else internal_config.default_rank_results_by
                )
            )
        )
        self.rank_results_by = rank_results_by
        if self.rank_results_by not in RANK_RESULTS_BY_ALLOWED_VALUES:
            raise ValueError(_e + f"Invalid value \"{self.rank_results_by}\" specified for parameter "
                                  f"\"rankResultsBy\" in section REPORT. Allowed options are: "
                                  f"{RANK_RESULTS_BY_ALLOWED_VALUES}")

        report_max_n_plotted = run_config.get(
            "REPORT", "maxNumberPlotted", fallback=run_config.get(
                "REPORT", "maxNPlotted", fallback=run_config.get(
                    "REPORT", "maxPlotted", fallback=(
                        user_config.user_default_report_max_n_plotted  # Default to None (no limit)
                    )
                )
            )
        )
        # Allow */all to represent no limit on number of plots (same as omitting the parameter)
        self.report_max_n_plotted = None if (
            report_max_n_plotted is None or report_max_n_plotted in ("all", "All", "ALL", "*")
        ) else int(report_max_n_plotted)

        report_min_total_n_occurrences_across_all_samples = run_config.get(
            "REPORT", "minTotalNumberOccurrencesAcrossAllSamples", fallback=run_config.get(
                "REPORT", "minTotalNOccurrencesAcrossAllSamples", fallback=run_config.get(
                    "REPORT", "minTotalOccurrencesAcrossAllSamples", fallback=(
                        user_config.user_default_report_min_total_n_occurrences_across_all_samples
                        if user_config.user_default_report_min_total_n_occurrences_across_all_samples is not None
                        else 1  # Default to 1 (excludes no activity)
                    )
                )
            )
        )
        self.report_min_total_n_occurrences_across_all_samples = int(
            report_min_total_n_occurrences_across_all_samples
        )

        report_min_n_occurrences_in_sample = run_config.get(
            "REPORT", "minNumberOccurrencesInSample", fallback=run_config.get(
                "REPORT", "minNOccurrencesInSample", fallback=run_config.get(
                    "REPORT", "minOccurrencesInSample", fallback=(
                        user_config.user_default_report_min_n_occurrences_in_sample
                        # Default to None - Check for paired use below, then set stored value to 0
                    )
                )
            )
        )

        report_occurrences_in_at_least_n_samples = run_config.get(
            "REPORT", "occurrencesInAtLeastNSamples", fallback=run_config.get(
                "REPORT", "inAtLeastNSamples", fallback=(
                    user_config.user_default_report_occurrences_in_at_least_n_samples
                    # Default to None - Check for paired use below, then set stored value to 0
                )
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

        report_draw_junctions_min_count = run_config.get(
            "REPORT", "drawJunctionsMinCount", fallback=run_config.get(
                "REPORT", "drawJunctionsMinN", fallback=run_config.get(
                    "REPORT", "drawJunctionsMin", fallback=(
                        user_config.user_default_report_draw_junctions_min_count
                        if user_config.user_default_report_draw_junctions_min_count is not None
                        else internal_config.default_draw_junctions_min_count
                    )
                )
            )
        )
        self.report_draw_junctions_min_count = int(report_draw_junctions_min_count)

        report_transcript_plot_max_samples_per_group = run_config.get(
            "REPORT", "transcriptPlotMaxNSamplesPerGroup", fallback=run_config.get(
                "REPORT", "transcriptPlotMaxSamplesPerGroup", fallback=None
            )
        )
        if report_transcript_plot_max_samples_per_group in (None, "*", "ALL", "All", "all"):
            self.report_transcript_plot_max_samples_per_group = None
        else:
            self.report_transcript_plot_max_samples_per_group = int(report_transcript_plot_max_samples_per_group)

        # Ensure valid parameters are set to enable featureCounts to run if report is enabled
        self.check_feature_counts()

        # Enable case-sensitivity for option names in subsequent (SAMPLES, COLORS) sections
        run_config.optionxform = str
        run_config.read(run_config_path)

        # Auto-generated if required: [SAMPLES]

        if "SAMPLES" in run_config:
            _non_redundant_keys = _remove_lowercase_duplicates(run_config.options("SAMPLES"))
            sample_groups = {
                key: [remove_quotes(s.strip()) for s in value.split(",")]
                for key, value in run_config.items("SAMPLES") if all([
                    value,
                    key in _non_redundant_keys
                ])
            }
            _flattened_sample_groups_values = flatten_nested_lists(list(sample_groups.values()))
            if not all([
                sample_name in sample_names
                for sample_name in _flattened_sample_groups_values
            ]):
                raise ValueError(_e + f"Sample names in run config do not match sample names parsed from BAM files.\n\n"
                                      f"Sample names parsed from config file: {_flattened_sample_groups_values}\n"
                                      f"Sample names parsed from BAM files: {sample_names}")
            if not all([
                sample_name in _flattened_sample_groups_values
                for sample_name in sample_names
            ]):
                raise ValueError(_e + f"Sample names parsed from BAM files are not all present in the SAMPLES "
                                      f"section of run config. Each sample name must be assigned to a group.\n\n"
                                      f"Sample names parsed from config file: {_flattened_sample_groups_values}\n"
                                      f"Sample names parsed from BAM files: {sample_names}"
                                 )
        else:
            sample_groups = {"All": sample_names}
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
            # color_by_group_name = dict(run_config.items("COLORS"))
            color_by_group_name = {k: remove_quotes(v) for k, v in run_config.items("COLORS")}
        else:
            color_by_group_name = {
                group_name: internal_config.default_point_color_in_plot
                for group_name in sample_groups.keys()
            }
        self.color_by_group_name = color_by_group_name

        # Auto-generated if required: [SHAPES]
        if "SHAPES" in run_config and "SAMPLES" in run_config:
            shape_by_group_name = {k: convert_marker_alias(remove_quotes(v)) for k, v in run_config.items("SHAPES")}
        else:
            shape_by_group_name = {
                group_name: internal_config.default_point_shape_in_plot
                for group_name in sample_groups.keys()
            }
        self.shape_by_group_name = shape_by_group_name

    def __setattr__(self, name, value):
        # If attribute is a string, strip all quotation marks before applying
        if isinstance(value, str):
            value = remove_quotes(value)
        super().__setattr__(name, value)

    def check_feature_counts(self) -> None:

        if all([
            self.generate_report,
            self.report_feature_counts_executable is None,
            self.report_gene_count_matrix is None
        ]):
            raise ValueError(f"Report generation requires that gene counts are either calculated or supplied. "
                             f"If the geneCounts parameter is not set in the REPORT section of run config "
                             f"file, then the featureCountsExecutable parameter must be set either in the "
                             f"REPORT section of run config, or in the DEFAULT REPORT section of user config.")
