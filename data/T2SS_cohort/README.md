# Reproducible data for the T2SS cohort figures

This directory contains the input data used to reproduce:

1. `T2SS_abundance_boxplot_Gp.pdf`
2. `T2SS_abundance_boxplot_Gv.pdf`
3. `gardnerella_abundance_filtered.pdf`

The corresponding R scripts are located in `../../code/`, and the reproduced outputs are written to `../output/`.


## Data notes

- The Gp and Gv CSV files contain samples passing the original T2SS completeness filter (`Completeness > 0.6`).
- The two TSV files contain per-sample relative abundance estimates.
- The combined abundance figure retains only species/sample pairs present in the corresponding filtered T2SS sample list.
- Sample identifiers and all quantitative columns are preserved from the analysis inputs used for the original figures.


