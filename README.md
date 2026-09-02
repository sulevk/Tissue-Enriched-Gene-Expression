# Tissue-enriched RNA-seq analysis workflow

This folder contains three R-based workflow scripts for tissue-specific expression analyses:

- `analyze_gtex_profiles_median_of_ratios.R`
- `DESeq2tissueRNAseq05July2026WithGreatPlot.Rmd`
- `TissueBoxplotsTop9genes.Rmd`

## Overview

1. `analyze_gtex_profiles_median_of_ratios.R`
   - Loads GTEx-style expression data and one-hot tissue metadata.
   - Filters low-count genes.
   - Normalizes counts using DESeq2 size factors.
   - Computes tissue means, variability, and tissue-vs-rest significance tests.
   - Writes summary tables and TSV/CSV outputs to a project-local output directory.

2. `DESeq2tissueRNAseq05July2026WithGreatPlot.Rmd`
   - Loads RNA-seq quantification files from a tissue directory.
   - Builds a DESeq2 dataset for tissue-specific differential expression.
   - Extracts per-tissue results using a single helper function.
   - Saves intermediate checkpoints and normalized count outputs.

3. `TissueBoxplotsTop9genes.Rmd`
   - Uses normalized DESeq2 counts to generate tissue-specific boxplots.
   - Uses a reusable plotting helper and a named gene-list map.
   - Saves plot PDFs to a `plots` folder.

## How to run

### Local machine

Open the script or notebook in RStudio or VS Code and run it from the project folder.
The scripts now automatically resolve the active project directory by checking, in order:

- `WORKSPACE_DIR`
- `PROJECT_ROOT`
- `PWD`
- working directory (`getwd()`)
- common cluster roots if needed

This means you do not need to hard-code a machine-specific absolute path.

Example:

```bash
cd /path/to/TissueSpecific
Rscript analyze_gtex_profiles_median_of_ratios.R
```

For the Rmd notebooks, open them in RStudio and run the chunks in order.

### Cluster or remote environment

Set one of the environment variables before running the scripts:

```bash
export PROJECT_ROOT=/path/to/project
# or
export WORKSPACE_DIR=/path/to/project
```

Then run the same commands as above.

## Expected input files

The scripts assume that the project directory contains the relevant data files, such as:

- `GTExWithTissues.RData`
- `matrix1_onehot_tissue.txt`
- `HumanTissueSamples205.xlsx`
- quantification folders under `tissue205quant/.../quant.sf`

If a required file is missing, the script will stop with a clear error.

## Output folders

The workflow writes outputs to folders created automatically under the active project directory:

- `gtex_profile_analysis_median_of_ratios/`
- `results/`
- `plots/`

These outputs include TSV tables, summary files, checkpoint `.rds` files, and plot PDFs.

## Best practices

- Keep the working directory at the project root or set `PROJECT_ROOT`/`WORKSPACE_DIR`.
- Do not edit the data paths manually unless the file layout changes.
- Re-run from a clean project directory if you want a fresh output set.
- If a large run fails mid-way, the checkpoint files can help recover intermediate objects.

## Notes

These scripts are designed to be reproducible and easier to maintain than the earlier fully copy-pasted versions. They aim to keep the same biological logic while reducing duplication, improving readability, and making execution more robust in both local and cluster environments.
