
from typing import Any
import os
from gtfparse import read_gtf

from transcript import (
    create_and_save_transcript_library, annotate_and_save_transcript_library, load_transcript_library_from_file
)
from entrez import download_and_save_feature_annotation_xml, get_gbseq_from_xml
from biomart import create_and_save_name_lookup, load_name_lookup_from_file
from experiment import Sample
from counts import run_feature_counts, get_gene_counts_from_tsv
from core import FaseConfig, perform_splice_analysis

# ==================================================================================================================== #

# GLOBAL SETTINGS


# Provide user identity to make entrez queries
EMAIL = "angus.shoppee@monash.edu"

# NEW IN V0.4: Species from which sequencing data originates
SPECIES = "Mouse"

# NEW IN V0.4: Path to single-column csv file / Leave empty to default to whole-genome
# NOTE: Currently unused
GENES_TO_ANALYZE = ""

# REFACTORED IN V0.4: Which annotated features do we want to analyze splice events for?
# * values: the strings that will be searched for in the genbank feature's qualifiers
FEATURES_TO_ANALYZE = [
    "transmembrane region"
]

# NEW (VERSION 0.3) - How much of the feature (as fraction, 0.0 to 1.0)
#                     should overlap with the splice junction?
FEATURE_JUNCTION_OVERLAP_THRESHOLD = 0.5

EXPERIMENTAL_ALLOW_JUNCTION_START_OUTSIDE_TRANSCRIPT = False
EXPERIMENTAL_ALLOW_JUNCTION_END_OUTSIDE_TRANSCRIPT = False

# Where to save output files
OUTPUT_DIR = "/Users/aasho2/Projects/FASE_V1/OUTPUT/V0_4"

# Sorted, indexed BAM files to be analyzed
BAM_FILES_DIR = "/Users/aasho2/PHD/Bioinformatics/STAR/runs/hons_PD1KO/sorted"

# Everything that comes after sample identifier in the bam paths, including file extension
BAM_SUFFIX = "_Aligned_Sorted.out.bam"

# Reference gtf file - make sure these are the same assembly your BAMs were aligned to
REFERENCE_GENOME_GTF_PATH = "/Users/aasho2/PHD/Bioinformatics/STAR/genomes/GRCm39/gencode_M27_primary/gencode.vM27.primary_assembly.annotation.gtf"

# Provide the full path to your featureCounts executable
FEATURECOUNTS_EXECUTABLE = "/Users/aasho2/opt/anaconda3/envs/bbmap/bin/featureCounts"

# NEW IN V0.4: Indicate to featureCounts whether reads are paired-end
PAIRED_END_READS = True

# NEW IN V0.4: Number of threads to use for featureCounts and splice analysis
NUMBER_OF_THREADS = 12

# REFACTORED IN V0.4: Path to existing featureCounts output - leave empty to run featureCounts
COUNTS_FILE = ""

# Optional: Choose a biomart mirror
BIOMART_MIRROR = 'http://asia.ensembl.org'

# NEW IN V0.4
FORCE_REGENERATE_WHOLE_GENOME_LOOKUP = False

# NEW IN V0.4
FORCE_REGENERATE_TRANSCRIPT_LIBRARY = False

# NEW IN V0.4
FORCE_REGENERATE_GENBANK_FEATURE_XML = False


# ==================================================================================================================== #

def main(verbose: bool = False) -> None:

    def print_if_verbose(s: Any):
        if verbose:
            print(s)

    print_if_verbose("Initializing...")

    base_dir = os.path.realpath(os.path.dirname(__file__))

    # Load config
    fase_config = FaseConfig(
        os.path.join(base_dir, "config", "internal.config")
    )

    # Validate user settings (not exhaustive)
    # TODO: Complete
    if SPECIES not in fase_config.allowed_species:
        raise ValueError(
            f"Provided species '{SPECIES}' is not supported. Currently supported species: {fase_config.allowed_species}"
        )

    # Create data directory for this species if required
    species_specific_data_dir = os.path.join(base_dir, "data", SPECIES)
    if not os.path.isdir(species_specific_data_dir):
        os.makedirs(species_specific_data_dir)

    name_lookup_path = os.path.join(species_specific_data_dir, "name_lookup.json")

    # Create name lookup if required
    if FORCE_REGENERATE_WHOLE_GENOME_LOOKUP or not os.path.exists(name_lookup_path):
        name_lookup = create_and_save_name_lookup(
            fase_config.biomart_name_for_species[SPECIES],
            BIOMART_MIRROR,
            name_lookup_path,
            verbose=verbose
        )

    # Load name lookup
    else:
        name_lookup = load_name_lookup_from_file(name_lookup_path)

    print_if_verbose("...done\n")

    annotated_transcript_library_path = os.path.join(species_specific_data_dir, "annotated_transcript_library.object")

    # Create annotated transcript library if required
    if any([
        FORCE_REGENERATE_TRANSCRIPT_LIBRARY,
        FORCE_REGENERATE_GENBANK_FEATURE_XML,
        not os.path.exists(annotated_transcript_library_path)
    ]):

        print_if_verbose("Annotated transcript library has not been built and will be initialized this run.\n")

        transcript_library_path = os.path.join(species_specific_data_dir, "transcript_library.object")

        # Create transcript library if required
        if FORCE_REGENERATE_TRANSCRIPT_LIBRARY or not os.path.exists(transcript_library_path):
            print_if_verbose("Importing reference GTF...")
            ref_gtf = read_gtf(REFERENCE_GENOME_GTF_PATH)
            print_if_verbose("...done\n")

            transcript_library = create_and_save_transcript_library(
                SPECIES,
                transcript_library_path,
                name_lookup.get_all_gene_names(),
                ref_gtf,
                verbose=verbose
            )
            del ref_gtf

        # Load transcript library
        else:
            transcript_library = load_transcript_library_from_file(transcript_library_path)

        genbank_feature_xml_path = os.path.join(species_specific_data_dir, "genbank_feature_annotation.xml")

        # Create genbank feature XML if required
        if FORCE_REGENERATE_GENBANK_FEATURE_XML or not os.path.exists(genbank_feature_xml_path):
            print_if_verbose("Downloading genbank feature annotation data...")

            transcripts_with_refseq = []
            for t in transcript_library.get_all_transcripts():
                try:
                    refseq = name_lookup.convert(t.transcript_id, "ensembl", "refseq")
                    if refseq:
                        transcripts_with_refseq.append(t.transcript_id)
                except KeyError:
                    pass

            download_and_save_feature_annotation_xml(
                list(transcripts_with_refseq),
                genbank_feature_xml_path,
                EMAIL,
                verbose=verbose
            )

            del transcripts_with_refseq

            print_if_verbose("...done\n")

        # Load genbank feature XML
        gbseq_by_refseq = get_gbseq_from_xml(
            genbank_feature_xml_path,
            verbose=verbose
        )

        # Apply feature annotation to transcript library
        annotated_transcript_library = annotate_and_save_transcript_library(
            transcript_library,
            name_lookup,
            gbseq_by_refseq,
            annotated_transcript_library_path,
            verbose=verbose
        )

        del transcript_library
        del gbseq_by_refseq

        print_if_verbose("Annotated transcript library has been created and saved.\n")

    # Load annotated transcript library
    else:
        print_if_verbose("Loading transcript library...")
        annotated_transcript_library = load_transcript_library_from_file(annotated_transcript_library_path)
        print_if_verbose("...done\n")

    # # START OF DEBUG BLOCK
    # print(annotated_transcript_library)
    # print("Number of transcripts:", len(annotated_transcript_library.get_all_transcripts()))
    #
    # for t in annotated_transcript_library.get_all_transcripts()[10:20]:
    # 	print(t.gene_name, "\t", t.transcript_id, "\t", t.seqname, t.strand, t.start, t.end)
    # 	print([(e.start, e.end) for e in t.exons])
    # 	print([(e.start, e.end) for e in t.gbseq.exons])
    # 	print(t.exon_map)
    # 	try:
    # 		print("local pos 85 --> genomic pos", t.genomic_position_from_local(85))
    # 	except ValueError:
    # 		print("local pos 85 couldn't be converted to genomic pos")
    # 	print()
    # # END OF DEBUG BLOCK

    print_if_verbose(f"Enabling screening for features containing term(s): \"{', '.join(FEATURES_TO_ANALYZE)}\"...")

    for feature_substring in FEATURES_TO_ANALYZE:
        for transcript in annotated_transcript_library.get_all_transcripts():
            for gbfeature in transcript.gbseq.features:
                if gbfeature.has_qual_value_containing(feature_substring):
                    transcript.gbseq.set_analysis_feature(feature_substring, gbfeature)

    print_if_verbose("... done\n")

    if not os.path.isdir(OUTPUT_DIR):
        if os.path.exists(OUTPUT_DIR):
            raise ValueError(f"The specified output directory ({OUTPUT_DIR}) already exists as a file.")
        os.makedirs(OUTPUT_DIR)

    feature_counts_output_path = COUNTS_FILE if COUNTS_FILE else os.path.join(OUTPUT_DIR, "gene_counts.tsv")

    # Run featureCounts if counts file has not been supplied and featureCounts has not been run in a previous analysis
    if not os.path.exists(feature_counts_output_path):

        print_if_verbose("Running featureCounts...")
        run_feature_counts(
            FEATURECOUNTS_EXECUTABLE,
            BAM_FILES_DIR,
            BAM_SUFFIX,
            REFERENCE_GENOME_GTF_PATH,
            feature_counts_output_path,
            paired_end_reads=PAIRED_END_READS,
            threads=NUMBER_OF_THREADS
        )
        print_if_verbose("...done\n")

    print_if_verbose(f"Importing gene count data...")
    gene_counts = get_gene_counts_from_tsv(feature_counts_output_path)
    print_if_verbose("...done\n")

    # # START OF DEBUG BLOCK
    # print(index_by_sample_name)
    # print(counts_by_gene_id["ENSMUSG00000016496.8"])
    # # END OF DEBUG BLOCK

    _suffix_length = len(BAM_SUFFIX)
    bam_file_absolute_paths = [os.path.join(
        BAM_FILES_DIR, p
    ) for p in os.listdir(BAM_FILES_DIR) if p[-_suffix_length:] == BAM_SUFFIX]

    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, BAM_SUFFIX)
        samples[sample.name] = sample

    # Run analysis
    perform_splice_analysis(
        FEATURES_TO_ANALYZE,
        FEATURE_JUNCTION_OVERLAP_THRESHOLD,
        samples,
        annotated_transcript_library,
        gene_counts,
        OUTPUT_DIR,
        verbose=verbose
    )


if __name__ == "__main__":

    main(verbose=True)