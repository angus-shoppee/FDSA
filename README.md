
# **Feature-Directed Splice Analysis (FDSA)**

## **Overview**
The Feature-Directed Splice Analysis (FDSA) pipeline is designed for analyzing splice events affecting specific transcript features in sequencing data. The pipeline relies on configuration files and command-line arguments for flexibility and reproducibility.

---

## **Runtime Modes**
FDSA comprises six runtime modes, divided into **setup** and **analysis** categories:

### **Setup Modes**
1. **`user`**: Specifies the location of the user config file, used for setting default parameters and paths to optional dependencies.  
   **Usage**: Run once after installation.  
   **Command**:  
   ```bash
   fdsa user $user_config_path
   ```
    <br>

2. **`build`**: Prepares the program for a specific species by processing a supplied reference genome and downloading 
                feature annotations.  
   **Usage**: Run once per species.  
   **Command**:  
   ```bash
   fdsa build --species $species_name --genome $reference_genome_gtf_path --threads $n_threads
   ```
   **Or, supplying a run config file with a valid [RUN] section specifying "species" and "genome"**:  
   ```bash
   fdsa build --threads $n_threads $run_config_path
   ```
   <br>

---

### **Analysis Modes**
3. **`run`**: Core analysis mode for detecting splice events.  
   **Usage**: Requires a run config file.  
   **Command**:  
   ```bash
   fdsa run $run_config_path
   ```
   **Mandatory Run Config Section**: `[RUN]`.

   <br>


4. **`filter`**: Produces filtered BAM files containing only reads aligned to genes with relevant splice activity.  
   **Usage**: Requires a run config file.  
   **Command**:  
   ```bash
   fdsa filter RUN_CONFIG_PATH
   ```
   **Mandatory Run Config Section**: `[RUN]`.

   **Optional Run Config Section**: `[FILTER]`.

   <br>


5. **`quant`**: Performs transcript quantification using StringTie, annotates predicted transcript structures, and 
   calculates an estimated percentage of transcripts lacking the designated feature(s).

   **Usage**: Can be supplied a run config file, or run independently if required arguments are set.  
   **Command**:  
   ```bash
   fdsa quant $run_config_path
   ```
   **Mandatory Run Config Section**: `[RUN]`.

   **Optional Run Config Section**: `[QUANT]`.
   <br><br>

   **Or, if running independently**:
   ```bash
   fdsa quant --stringtie $stringtie_path --prep-de $prep_de_path --input $bam_files_dir --output $output_dir --read-length $read_length --threads $n_threads
   ```
   *Note: If run independently, transcripts will not be annotated unless an FDSA results file is supplied (see below)*

   <br>


6. **`report`**: Generates a visual report from the output of the `run` and `quant` modes.  
   **Usage**: Requires a run config file with optional sections for customization.

   **Command**:  
   ```bash
   fdsa report RUN_CONFIG_PATH
   ```
   **Mandatory Run Config Section**: `[RUN]`.

   **Optional Sections**: `[REPORT]`, `[SAMPLES]`, `[COLORS]`, `[SHAPES]`.

   <br>


### **Example Setup & Usage**
   ```bash
   fdsa user $path_to_user_config_file
   fdsa build --species human --genome $path_to_genome_gtf --threads $n_threads
   fdsa run --filter --quant --report $path_to_run_config_file
   ```

---

## **Configuration Files**

Parameters may either be strings, integers, decimal numbers, or booleans. For example:
- `species = "human"` (STRING)
- `threads = 12` (INTEGER)
- `featureJunctionOverlapThreshold = 0.5` (DECIMAL)
- `enableReport = False` (BOOLEAN)

Some parameters (indicated below) can also be delimited. For example:
- `geneList = "FAS", "BTLA", "CTLA4"`

String parameters should ideally be surrounded with quotation marks (`"`), though they are not strictly required in many cases.

Certain integer parameters (indicated below) relating to threshold/maximum amounts can be set to the wildcard character (`*`) to allow any amount.

Boolean parameters can be specified as either `True`/`False` or `Yes`/`No` and are case-insensitive.

### **User Config**
- **Purpose**: Sets default parameters and paths to optional dependencies.
- **Mandatory Section**:
  - `[BUILD]`
    - `email` / `emailAddress`: Email address to attach to GenBank API queries.
    - `biomart` / `biomartMirror`: (Optional) Ensembl base URL for API queries.
    <br>Default = "http://asia.ensembl.org"

- **Optional Sections**:
  - `[DEFAULT RUN]`
    <br>*For all parameters, see Run Config below*
  - `[DEFAULT FILTER]`
    - `samtools` / `samtoolsExecutablePath`: Path to the samtools executable file.
    <br>*For all parameters, see Run Config below*
  - `[DEFAULT QUANT]`
    - `stringtie / stringtieExecutablePath`: Path to the stringtie executable file.
    - `prepDE / prepDEScriptPath`: Path to the prepDE.py3 script distributed alongside stringtie.
    <br>*For all parameters, see Run Config below*
  - `[DEFAULT REPORT]`
    - `featureCounts` / `featureCountsExecutablePath`: Path to the featureCounts executable file.
    <br>*For all parameters, see Run Config below*

### **Run Config**
- **Purpose**: Specifies parameters for each analysis mode.
- **Mandatory Sections**:
  - `[RUN]`
    - (Mandatory) `name` / `runName`: Name to use when generating output directories and files.
    - (Mandatory) `feature`: Exact text to be matched with feature annotation records to identify transcript features. For example: `feature = "transmembrane region"`
    - (Mandatory) `input` / `inputPath`: Directory containing indexed BAM files.
    - (Mandatory) `bamEnding`: Everything following the sample identifying part of the input BAM file names. For example: `bamEnding = "_Aligned.sortedByCoord.out.bam"`
    - (Mandatory) `pairedEnd` / `pairedEndReads`: Indicates whether input reads are paired-end type (True/False).
    - (Mandatory) `species`: Species of input data.
    - (Mandatory) `genome` / `referenceGenome`: Path to the reference genome GTF file used to align input BAM files. **IMPORTANT: Input read files must be aligned to the same reference genome assembly as the reference GTF used to run `fdsa build` for the corresponding species**.
    - (Mandatory) `output` / `outputPath`: Location to generate output directories and files.
    - `threads`: Number of threads to use for each enabled runtime mode.
    <br>Default = 1
    - `maxFeaturesInTranscript` / `maxNumberFeaturesInTranscript`: Maximum number of the specified feature annotation allowed in a given transcript for it to be included in the analysis. For example, combining `feature = "transmembrane region"` and
      `maxNumberFeaturesInTranscript = 1` will exclude multi-pass receptors from the analysis.
    <br>Default = *
    - `overlapThreshold` / `featureJunctionOverlapThreshold`: Minimum exonic overlap (fraction, 0-1) of a splice junction with a feature locus required to count a junction as overlapping the feature.
    <br>Default = 0.5
    - `mapqForUniqueMapping`: Map quality score used by the aligner that generated input BAM files to indicate a unique alignment. For example, STAR aligner uses a mapq of 255 to indicate uniquely aligned reads.
    <br>Default = 255
    - `includeAllJunctionsInOutput`: Include a column in the run results file listing all splice junctions for each gene (Required for downstream runtime modes).
    <br>Default = True
    - `useOnlyLongestAnnotatedTranscript`: If there are multiple reference transcripts for a gene, only use the longest transcript with feature annotation information. Disabling this option will produce redundant output, but may catch edge cases where annotation of reference transcript variants is incomplete.
    <br>Default = True
    - `skipTranscriptsWithRedundantFeatureAnnotation` If there are multiple reference transcripts for a gene, skip generating results for multiple transcripts with duplicate feature annotations. Disabling this option will produce redundant output, but may catch edge cases where annotation of reference transcript variants is incomplete.
    <br>Default = True
    - `genes` / `geneList`: If desired, perform the analysis for a supplied list of genes instead of screening.
    <br>*Space- or comma-delimited*
    <br>*Alternatively, a path to a headerless single-column .csv file may be supplied*
    <br>Default = None
    - `filter` / `enableFilter`: Enable `filter` mode after run completion. Analogous to using the "--filter" command line flag.
    <br>Default = False
    - `quant` / `enableQuant`: Enable `quant` mode after run completion. Analogous to using the "--quant" command line flag.
    <br>Default = False
    - `report` / `enableReport`: Enable `report` mode after run completion. Analogous to using the "--report" command line flag.
    <br>Default = False
- **Optional Sections**:
  - `[FILTER]`
    - `unique` / `uniqueMappingOnly`: Specifies whether only uniquely mapped reads should be included in the filtered BAM output.
    <br>Default = True
    - `minTotalOccurrencesAcrossAllSamples` / `minTotalNumberOccurrencesAcrossAllSamples`: Minimum overall number of splice events across all samples required for a gene to be included in the filtered BAM output.
    <br>Default = 1
    - *The following two options are to be used as a pair:*
      - A^ `minOccurrencesInSample` / `minNumberOccurrencesInSample`: Only genes with a minimum of A^ splice events in at least B^ samples will be included in the filtered BAM output.
      <br>Default = 1
      - B^ `inAtLeastNSamples` / `occurrencesInAtLeastNSamples`: Only genes with a minimum of A^ splice events in at least B^ samples will be included in the filtered BAM output.
      <br>*The wildcard character `*` may be used to specify that there must be at least A^ reads in ALL samples*
      <br>Default = 1
  - `[QUANT]`
    - `readLength` / `inputReadLength`: Length of reads in input BAM files.
    <br>Default = 75
    - `assignReferenceGene`: Enables inference of reference gene IDs and reference gene names for stringtie
                        results without these values. When enabled, if a given stringtie gene
                        has at least one transcript with an associated reference gene ID and
                        reference gene name, `fdsa quant` will assign these values to all other
                        transcripts for that gene, as long as there is no ambiguity
                        (transcripts from a single stringtie gene associated with multiple
                        reference gene IDs).
    <br>Default = True
  - `[REPORT]`
    - `name` / `reportName`: Name to use when creating the HTML report file. Overrides run name. If not specified, the run name will be used. Setting this parameter allows multiple reports with different configurations to be created without overwriting files.
    <br>Default = None
    - `genes` / `geneList`: If desired, the report will be generated for a supplied list of genes, rather than genes which meet splice activity criteria.
    <br>*Space- or comma-delimited*
    <br>*Alternatively, a path to a headerless single-column .csv file may be supplied*
    <br>Default = None
    - `rankBy` / `rankResultsBy`: ("number" OR "frequency") Specifies whether results included in the report should be ordered by the frequency or overall number of splice events.
    <br>Default = "frequency"
    - `maxPlotted` / `maxNumberPlotted`: Maximum number of genes for which plots should be drawn.
    - `minTotalOccurrencesAcrossAllSamples` / `minTotalNumberOccurrencesAcrossAllSamples`: Minimum overall number of splice events across all samples required for a gene to be included in the report.
    <br>Default = 1
    - *The following two options are to be used as a pair:*
      - A^ `minOccurrencesInSample` / `minNumberOccurrencesInSample`: Only genes with a minimum of A^ splice events in at least B^ samples will be included in the report.
      <br>Default = 1
      - B^ `inAtLeastNSamples` / `occurrencesInAtLeastNSamples`: Only genes with a minimum of A^ splice events in at least B^ samples will be included in the report.
      <br>*The wildcard character `*` may be used to specify that there must be at least A^ reads in ALL samples*
      <br>Default = 1
    - `drawJunctionsMin` / `drawJunctionsMinCount`: Minimum number of junction counts (only applies to junctions NOT overlapping features) required for a junction to be included in a splice graph.
    <br>Default = 2
    - `transcriptPlotMaxSamplesPerGroup`: Maximum number of samples per experimental group to include in the splice graph section for a gene. Useful to prevent bloat in analyses where there are many samples per experimental group.
    <br>Default = *
    - `geneCounts` / `featureCountsOutput`: If desired, a path to a gene count matrix can be supplied to bypass the need to run featureCounts. If featureCounts has already been run by `fdsa report`, its output will be detected automatically.
    <br>Default = None
    - `saveNormGeneCounts` / `saveNormalizedGeneCounts`: Write TMM normalized gene counts to a file instead of performing normalization each time `fdsa report` is run.
    <br>Default = True
- **Report Customization Sections**:
  - `[SAMPLES]`
    - *Designate experimental groups in the following format:*
      - `group_x = "SRR1000001", "SRR1000002", "SRR1000003"`
      - `group_y = "SRR1000004", "SRR1000005", "SRR1000006"`
  - `[COLORS]`
    - *Designate colours to associate with experimental groups in the following format:*
      - `group_x = "red"`
      - `group_x = "#0000FF"`
      <br>*Either hex values or color names recognized by matplotlib may be used*
  - `[SHAPES]`
    - *Designate point shapes to associate with experimental groups in the following format:*
      - `group_x = "circle"`
      - `group_x = "square"`
      <br>*Either matplotlib shape codes or aliases such as "circle", "triangle", "star" may be used*
      <br>*All available aliases can be found in src/config/marker_aliases.py*

---

## **Command-Line Arguments**
For each mode, required arguments and options are parsed via the command line. Below is a summary of notable arguments for each mode:

### **User**
- `user_config_path`: Path to a config file with a valid `[BUILD]` section, and optional `[DEFAULT RUN]`,
                      `[DEFAULT FILTER]`, `[DEFAULT QUANT]`, and `[DEFAULT REPORT]` sections.

### **Build**

Positional argument:
  - `run_config_path`: (Optional) Path to a run config file containing a valid `[RUN]` section with "species" and "genome"
  parameters set.

Options (email only required if not defined in user config file; species & genome only required if no run config file is supplied):
 - `-s SPECIES, --species SPECIES`: Specifies the species. Options: mouse/human (*See note 1).
 - `-g GENOME, --genome GENOME`: Path to the reference genome GTF file.
 - `-e EMAIL, --email EMAIL`: (If not already specified in user config) Identifying email for GenBank API queries.
 - `-t THREADS, --threads THREADS`: Number of threads to use for parallel processing.
 - `--regenerate-lookup`: Re-downloads transcript information from BioMart and regenerates the
                    locally stored lookup used to convert between transcript IDs, gene IDs
                    and gene names.
 - `--regenerate-transcripts`: Re-parses transcripts from the provided reference genome GTF to create
                    a stored transcript library.
 - `--regenerate-features`: Re-downloads transcript feature annotation data from GenBank and
                    stores the information locally.
 - `--redo-annotation`: Re-applies feature annotation to the stored transcript library.
 - `--regenerate-all`: Shorthand for --regenerate-lookup --regenerate-transcripts
                    --regenerate-features --redo-annotation
 - `--regenerate-partial-features-download`: When combined with the --regenerate-features flag, forces any data
                    partially downloaded data from GenBank to be discarded, restarting the
                    download from scratch.

### **Run**

Positional argument:
  - `run_config_path`: Path to a config file with a valid `[RUN]` section, and optional `[FILTER]`,
                       `[QUANT]`, `[REPORT]`, `[SAMPLES]`, `[COLORS]`, and `[SHAPES]` sections.

Options:
  - `--filter`: Chains results into filter mode.
  - `--quant`: Chains results into quant mode
  - `--report`: Generates a HTML report upon completion.

### **Filter**

Positional argument:
  - `run_config_path`

### **Quant**

Positional argument:
  - `run_config_path`: (Optional)

Options (only required if no run config file is supplied):
  - `-s STRINGTIE_PATH, --stringtie STRINGTIE_PATH`: Path to the StringTie executable.
  - `-p PREP_DE_PATH, --prep-de PREP_DE_PATH`: Path to the prepDE.py3 script distributed together with StringTie.
  - `-i INPUT, --input INPUT`: Input directory containing BAM files.
  - `-o OUTPUT, --output OUTPUT`: Output base directory.
  - `-l READ_LENGTH, --read-length READ_LENGTH`: Read length in input files (optional; StringTie default = 75)
  - `-f FDSA_RESULTS, --fdsa-results FDSA_RESULTS`: Path to fdsa run output file (optional; quantified transcripts will not be annotated with feature information if not set)

### **Report**

Positional argument:
  - `run_config_path`

---

## **Dependencies**
- **Mandatory**:
  - Python 3.8+
  - Core libraries (e.g., argparse)
  - Automatically installed with pip installation: conorm, entrezpy, gtfparse, matplotlib, numpy, pandas, pybiomart, pysam, setuptools
- **Optional**:
  - `Samtools`: Used in `filter` mode.
  - `StringTie`: Used in `quant` mode.
  - `FeatureCounts`: Used in `report` mode.


---

## **Installation**
1. Install program and required dependencies:
   ```bash
   pip install .
   ```
2. Configure user settings:
   ```bash
   fdsa user $user_config_path
   ```

---

## **Notes**
1. Currently, only mouse and human are officially supported. However, this is simply because only these species have
   been tested and confirmed to work as intended. This program is species-agnostic, and any species should be usable as
   long as it has a reference genome and Genbank transcript annotations. To add support for a species:
   - Edit the [BUILD] section in src/config/internal.config
   - Add a user-facing species name to `allowedSpecies`, e.g. "human"
   - Add a mapping from user-facing species name to BioMart gene name to `biomartNameForSpecies`, e.g.
     "human:hsapiens_gene_ensembl"

---
