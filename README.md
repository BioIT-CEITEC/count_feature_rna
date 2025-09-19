# Count Feature RNA Workflow

This repository contains a Snakemake workflow for RNA-seq quantification using multiple tools including featureCount, HTSeq-count, RSEM, Kallisto, and Salmon.

## Features
- **Modular Design:** Each step is implemented as a separate rule and wrapper, allowing easy customization.
- **Conda Environments:** Each rule uses its own environment for reproducibility.
- **BioRoot Utilities:** Integrates with external modules for sample and reference management.
- **Configurable:** All parameters (paths, organism, alignment options, etc.) are set via `workflow.config.json`.
  
## Workflow Overview

- **featureCount** (Subread package) - Gene/feature counting from BAM files
- **HTSeq-count** - Python-based read counting for features
- **RSEM** - RNA-Seq by Expectation Maximization for transcript quantification
- **Kallisto** - Fast transcript quantification using pseudoalignment
- **Salmon** - Fast and accurate transcript quantification
  - Salmon alignment mode (BAM-based)
  - Salmon mapping mode (FASTQ-based)

## Directory Structure
- `Snakefile`: Main workflow file.
- `workflow.config.json`: Configuration file.
- `rules/`: Snakemake rule files for each workflow step.
- `wrappers/`: Scripts and conda environments for each step.
- `qc_reports/`: Output directory for QC results and reports.
- `mapped/`: Input directory for alignment BAM files.
- `logs/`: Log files for each step.

### Key Capabilities
- Multiple feature types support (exons, genes, transcripts, 3'UTR, 5'UTR)
- Strand-specific and unstranded RNA-seq data support
- Single-end and paired-end data compatibility
- Comprehensive quality control reporting via MultiQC
- Flexible configuration through JSON config files
- Integration with BioRoot utilities for reference management

## Usage
1. **Configure the workflow:**
   - Edit `workflow.config.json` to specify sample information, reference genome, and parameters.
2. **Run the workflow:**
   ```bash
   snakemake --cores <N>
   ```
   Replace `<N>` with the number of CPU cores to use.
3. **Outputs:**
   - Count files for each tool per samples and summary HTML files.

## Requirements
- [Snakemake >=5.18.0](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/)
- Python 3

## Customization
- Modify rules or wrapper scripts to adapt to specific project needs.
- Add or remove steps by editing the `rules/` and `wrappers/` directories.

## Contact
For questions or contributions, please contact the BioIT-CEITEC team.
