# Gardnerella T2SS: reproducible code and data

This package contains the minimum direct inputs and derived plotting tables for the TerL cladogram (panel E) and the Gv142-2_P1 genome map (panel F). It also preserves the available inputs for rerunning the TerL analysis. 

The `code/` and `data/T2SS_cohort/` directories contain the R scripts and input data required to reproduce the *G. piotii* and *G. vaginalis* T2SS abundance boxplots and the filtered *Gardnerella* relative-abundance figure. The corresponding outputs are provided in `data/output/`.

## Reproduce the T2SS cohort figures

The three cohort figures are generated with R. Install the required packages:

```r
install.packages(c(
  "ggplot2",
  "ggpubr",
  "dplyr",
  "readr",
  "tidyr",
  "patchwork",
  "ragg",
  "png"
))
```
## Reproduce panels E and F

Tested with Python 3.9.6 on macOS ARM64. Install into an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/reproduce.py
```

The command validates the supplied data, exports tab-separated source-data tables, and redraws both panels. Outputs are written to `generated/figure_E/` and `generated/figure_F/`. TSV tables in `source_data/` can be imported into Excel (tab delimiter; identifiers and sequence columns as text). They are the canonical tabular source data; native Newick, FASTA and GenBank files are retained to avoid losing biological structure. No Excel software is required.

Use `python scripts/export_source_data.py` for data extraction and validation only. 

From the repository root, run:

```bash
Rscript code/run.R
```
The reproduced figures and associated summary files are written to
data/output/.

## Data and provenance

| Panel | Direct plotting input | Preserved upstream inputs |
|---|---|---|
| E | `data/figure_E/TerL_final.treefile` and label/category definitions in `scripts/plot_figure_E.py` | 19 original protein sequences, original headers, MAFFT alignment, trimmed alignment, IQ-TREE report and log |
| F | `data/figure_F/Gv142-2_P1.gbk` | Complete nucleotide sequence and the 71 CDS annotations in the same GenBank file; a FASTA export is also provided |


## Panel E interpretation and reanalysis

There are 19 taxa (5 study prophages, 2 published/reconstructed references and 12 public homologues). The supplied trimmed alignment contains 552 columns. Node labels are SH-aLRT/UFBoot. The plot displays labels only when both values are at least 80. The original tree representation is retained without midpoint or outgroup rooting. Horizontal line lengths in the cladogram represent topology steps, not inferred evolutionary distances. The apparent basal position of a tip in this drawing must not be interpreted as evidence of an evolutionary root.

Original recorded analysis: IQ-TREE 3.1.3; ModelFinder (`-m MFP`) selected Q.PFAM+G4 by BIC; 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates; seed 886838. The original command used `-T AUTO`. 

For an optional fresh alignment and tree analysis, install MAFFT, trimAl and IQ-TREE 3, then run:

```bash
bash scripts/reinfer_tree.sh
```

## Panel F interpretation

The genome is 37,522 bp with GC 43.6597196311% and 71 CDS features. 

GC tracks use 500-bp windows, 100-bp steps, truncated terminal windows and omission of windows shorter than 50 bp: 375 windows in total. The gray track, labelled “GC content” in the original figure, actually plots window GC fraction minus genome GC fraction. GC skew is `(G-C)/(G+C)`, with zero when the denominator is zero. These windows do not wrap across the sequence end.

The GenBank record declares a **linear** sequence. The circular drawing is a visualization layout and does not demonstrate a biologically circular molecule.

## Source tables

* `Figure_E_taxa_sequences.tsv`: tip order, IDs, display names, groups, original/aligned/trimmed protein sequences.
* `Figure_E_nodes.tsv`: node and parent IDs, child ordering, branch lengths, both support statistics and descendant tip IDs. Together these preserve tree structure, not just a list of support values.
* `Figure_E_newick.tsv`: complete original Newick string.
* `Figure_F_CDS.tsv`: all 71 CDS coordinates, strands, products, categories, colors and plot labels. Start/end columns use 1-based inclusive coordinates; the separate BioPython location string uses 0-based half-open coordinates.
* `Figure_F_GC_tracks.tsv`: window bounds (0-based, end-exclusive), midpoint, G/C counts, GC fraction, genome GC, GC deviation and skew. Fractions are on a 0–1 scale.
* `Figure_F_genome.tsv` and `Figure_F_sequence_chunks.tsv`: metadata and ordered nucleotide sequence chunks with 1-based inclusive coordinates.


