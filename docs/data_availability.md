# Data Availability

This repository includes compact result tables, selected figure outputs, and a final PDF report. It does not include raw sequencing files, Salmon index files, large model binaries, or archived development workspaces.

## Included

- Curated MATLAB and Python source code
- Original GEMsplice reference code and license
- Summary TSV comparison tables
- Selected PNG figures and figure text summaries
- Compact CSV pathway and summary tables
- One final PDF report

## Not Included

- FASTQ/BAM/CRAM files
- Salmon index directories
- Full `quant.sf` output trees
- Large `.mat` model/result binaries
- Human-GEM model binaries
- Old solver, notebook, and local environment trial workspaces
- Intermediate scratch outputs

## Why Large Data Are Excluded

The excluded files are either too large for a clean GitHub portfolio repository, are generated intermediates, or belong to local computational environments. The purpose of this repository is to present the final analysis logic and selected interpretable outputs.

## Reproducing The Full Analysis

To rerun the full workflow, the user must provide compatible RNA-seq quantification files, Human-GEM/GEMsplice model resources, transcript-to-gene mapping files, and MATLAB-compatible model structures. The included code documents the workflow structure but does not package all raw inputs.
