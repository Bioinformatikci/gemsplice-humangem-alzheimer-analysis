# Workflow

This project summarizes a curated Alzheimer/control analysis workflow built around transcript-level quantification, Human-GEM reaction mapping, and GEMsplice/iMAT-style model analysis.

## High-Level Steps

1. Quantify transcript abundance from RNA-seq samples with Salmon.
2. Map transcript-level expression to gene and reaction-level model features.
3. Prepare Alzheimer and control patient structures for GEMsplice-compatible analysis.
4. Build iMAT input structures using log2-transformed expression values.
5. Load the curated HumanGEM/GEMsplice model state.
6. Generate Alzheimer/control comparison plots and summary tables.
7. Compare tx2gene, HumanGEM-rule, and Alzheimer-objective result sets.
8. Compile final figures, tables, and interpretation into a single report.

## Code Entry Points

- `src/python/build_isoform_mapping.py`: builds transcript/gene/reaction mapping tables.
- `src/matlab/build_patients_from_salmon_quants.m`: converts Salmon quantification outputs into patient-level structures.
- `src/matlab/build_imat_dati_micro.m`: prepares iMAT expression inputs.
- `src/matlab/load_imat_gemsplice_model.m`: loads the project model state.
- `src/matlab/plot_GEMsplice_iMAT_alz_vs_control.m`: generates Alzheimer/control comparison figures.
- `src/matlab/write_imat_gemsplice_summary.m`: writes compact summary outputs.
- `src/python/build_final_master_report.py`: assembles the final report.

## Result Sets

- `results/figures/humangem_tx2gene/`: HumanGEM tx2gene log2 result plots and tables.
- `results/figures/imat_humangem_rules/`: iMAT results using HumanGEM reaction rules.
- `results/figures/imat_ad_objectives/`: iMAT results using Alzheimer-focused objective comparisons.

## Excluded Work

Exploratory solver trials, local environment setup attempts, old notebooks, and archived development folders were intentionally excluded from this curated repository.
