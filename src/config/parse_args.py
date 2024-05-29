
from argparse import ArgumentParser


def get_mode_parser(program_description: str) -> ArgumentParser:

    mode_parser = ArgumentParser(
        description=program_description,
        add_help=False
    )

    mode_parser.add_argument(
        "mode",
        type=str,
        choices=["user", "build", "run", "report", "filter", "quant"],
        # TODO: Format mode help text
        help="\n".join([
            "user - (Required only once) Sets the user config file containing build information and optional default "
            "parameters.",
            "build - (Required once per species, before first run) Parses transcripts from the provided reference "
            "genome GTF, downloads feature annotation data, and links the two to enable directed splice analysis.",
            "run - Performs splice analysis according to parameters set in the specified run config file.",
            "report - Generates a graphical report from the output of \"fase run\".",
            "filter - Generates filtered BAM files containing only reads aligned to genes that have splice events in "
            "the output of \"fase run\" (thresholds are configurable via the [FILTER] section in run config). "
            "quant - Wraps stringtie: Stringtie is run individually for each supplied BAM file, then results are "
            "merged, and stringtie is re-run using the merged GTF as a reference; finally, estimated transcript counts"
            "are extracted using stringtie's prepDE.py3 script."
        ])
    )

    return mode_parser


def get_build_parser() -> ArgumentParser:

    build_parser = ArgumentParser()

    build_parser.add_argument(
        "run_config_path",
        type=str,
        nargs="?",
        default=None,
        help=(
            "Path to a run config file containing a [RUN] section with the \"species\" and \"genome\" "
            "parameters set. If not specified, the --species and --genome options required."
        )
    )
    build_parser.add_argument(
        "-s",
        "--species",
        type=str,
        required=False,
        help="Species for which transcript annotation data will be downloaded. Must match the species of the "
             "provided reference genome GTF."
    )
    build_parser.add_argument(
        "-g",
        "--genome",
        type=str,
        required=False,
        help="Reference genome GTF from which transcripts will be parsed."
    )
    build_parser.add_argument(
        "-e",
        "--email",
        type=str,
        required=False,
        help="An identifying email address to be attached to GenBank API queries. Note: no emails will be "
             "sent - used for identification/anti-spam purposes only."
    )
    build_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        help="Number of threads to use during build process for speed improvements"
    )
    build_parser.add_argument(
        "--regenerate-lookup",
        action="store_true",
        help="Re-downloads transcript information from BioMart and regenerates the locally stored "
             "lookup used to convert between transcript IDs, gene IDs and gene names."
    )
    build_parser.add_argument(
        "--regenerate-transcripts",
        action="store_true",
        help="Re-parses transcripts from the provided reference genome GTF to create a stored transcript "
             "library."
    )
    build_parser.add_argument(
        "--regenerate-features",
        action="store_true",
        help="Re-downloads transcript feature annotation data from GenBank and stores the information locally."
    )
    build_parser.add_argument(
        "--redo-annotation",
        action="store_true",
        help="Re-applies feature annotation to the stored transcript library."
    )
    build_parser.add_argument(
        "--regenerate-all",
        action="store_true",
        help="Shorthand for --regenerate-lookup --regenerate-transcripts "
             "--regenerate-features --redo-annotation "
    )
    build_parser.add_argument(
        "--regenerate-partial-features-download",
        action="store_true",
        help="When combined with the --regenerate-features flag, forces any data partially downloaded data from "
             "GenBank to be discarded, restarting the download from scratch."
    )

    return build_parser


def get_run_parser() -> ArgumentParser:

    run_parser = ArgumentParser()

    run_parser.add_argument(
        "run_config_path",
        nargs="?",
        default=None,  # Allow no value so missing value error can be manually handled
        type=str,
        help=(
            "Path to a run config file containing a [RUN] section with the mandatory \"name\" (run name), \"feature\", "
            "\"species\", \"input\", \"bamEnding\", and \"output\" parameters set."
        )
    )
    run_parser.add_argument(
        "--report",
        action="store_true",
        help="Generates a visual report for the output of FASE's analysis upon completion."
    )
    run_parser.add_argument(
        "--no-report",
        action="store_true",
        help="Disables automatic report generation, overriding default behaviour and behaviour set in run config."
    )
    run_parser.add_argument(
        "--filter",
        action="store_true",
        help=(
            "Generates filtered BAM files containing only reads aligning to genes that have splice events in FASE's "
            "output upon completion."
        )
    )
    run_parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Disables filtered BAM file generation, overriding default behaviour and behaviour set in run config."
    )

    return run_parser


def get_report_parser() -> ArgumentParser:

    report_parser = ArgumentParser()

    report_parser.add_argument(
        "run_config_path",
        nargs="?",
        default=None,  # Allow no value so missing value error can be manually handled
        type=str,
        help=(
            "Path to a run config file containing a valid [RUN] section (see \"fase run --help\") and optional "
            "[REPORT], [SAMPLES], and [COLORS] sections."
            ""
        )
    )

    return report_parser


def get_filter_parser() -> ArgumentParser:

    filter_parser = ArgumentParser()

    filter_parser.add_argument(
        "run_config_path",
        nargs="?",
        default=None,  # Allow no value so missing value error can be manually handled
        type=str,
        help=(
            "Path to a run config file containing a valid [RUN] section (see \"fase run --help\") and optional "
            "[FILTER] section."
            ""
        )
    )

    return filter_parser


def get_quant_parser() -> ArgumentParser:

    quant_parser = ArgumentParser()

    quant_parser.add_argument(
        "run_config_path",
        nargs="?",
        default=None,  # Allow no value so missing value error can be manually handled
        type=str,
        help="Path to a run config file containing a valid [RUN] section (see \"fase run --help\"). Can be omitted if "
             "all of the -g, -i and -o parameters are set."
    )

    quant_parser.add_argument(
        "-s",
        "--stringtie",
        type=str,
        required=False,
        help="(If not defined in user config or run config) Path to the stringtie executable"
    )

    quant_parser.add_argument(
        "-p",
        "--prep-de",
        type=str,
        required=False,
        help="(If not defined in user config or run config) Path to the prepDE.py3 script provided with stringtie"
    )

    quant_parser.add_argument(
        "-g",
        "--genome",
        type=str,
        required=False,
        help="(If not defined in run config) Path to the reference genome annotation (GTF) file"
    )

    quant_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        help="(If not defined in run config) Input directory containing filtered BAM files (must end in \".bam\")"
    )

    quant_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="(If not defined in run config) Output directory where stringtie results will be saved"
    )

    quant_parser.add_argument(
        "-l",
        "--read-length",
        type=int,
        required=False,
        help="(If not defined in run config, optional) Input read length (stringtie default = 75)"
    )

    quant_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        help="(If not defined in run config, optional) Number of threads for stringtie to use"
    )

    return quant_parser
