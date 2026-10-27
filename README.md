# Gv manuscript: reproducibility of panels E and F

This package contains the minimum direct inputs and derived plotting tables for the TerL cladogram (panel E) and the Gv142-2_P1 genome map (panel F). It also preserves the available inputs for rerunning the TerL analysis. This is a figure-level package, not the complete study or a raw-read-to-assembly pipeline.

## Reproduce the figures

Tested with Python 3.9.6 on macOS ARM64. Install into an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/reproduce.py
```

The command validates the supplied data, exports tab-separated source-data tables, and redraws both panels. Outputs are written to `generated/figure_E/` and `generated/figure_F/`. TSV tables in `source_data/` can be imported into Excel (tab delimiter; identifiers and sequence columns as text). They are the canonical tabular source data; native Newick, FASTA and GenBank files are retained to avoid losing biological structure. No Excel software is required.

Use `python scripts/export_source_data.py` for data extraction and validation only. The scripts locate the repository relative to their own file paths, so they may be invoked from another working directory.

## Data and provenance

| Panel | Direct plotting input | Preserved upstream inputs |
|---|---|---|
| E | `data/figure_E/TerL_final.treefile` and label/category definitions in `scripts/plot_figure_E.py` | 19 original protein sequences, original headers, MAFFT alignment, trimmed alignment, IQ-TREE report and log |
| F | `data/figure_F/Gv142-2_P1.gbk` | Complete nucleotide sequence and the 71 CDS annotations in the same GenBank file; a FASTA export is also provided |

`provenance/source_manifest.json` maps every copied original to its path relative to the original Gv folder and records its SHA-256 hash. Original scripts are preserved unchanged in `provenance/original_scripts/`. They are archival context, not the entry point; some contain commands for other analyses and require the original folder layout. Run only the adapted entry point above for these two panels.

The original collected GenBank file is byte-identical to `results/10_pharokka/Gv142-2_P1/Gv142-2_P1.gbk` in the source project. The reference figures are the original standalone PNG/PDF exports, before manuscript layout editing. Files under `generated/` are local verification results and are ignored by Git.

## Panel E interpretation and reanalysis

There are 19 taxa (5 study prophages, 2 published/reconstructed references and 12 public homologues). The supplied trimmed alignment contains 552 columns. Node labels are SH-aLRT/UFBoot. The plot displays labels only when both values are at least 80. The original tree representation is retained without midpoint or outgroup rooting. Horizontal line lengths in the cladogram represent topology steps, not inferred evolutionary distances. The apparent basal position of a tip in this drawing must not be interpreted as evidence of an evolutionary root.

Original recorded analysis: IQ-TREE 3.1.3; ModelFinder (`-m MFP`) selected Q.PFAM+G4 by BIC; 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates; seed 886838. The original command used `-T AUTO`. See `provenance/TerL_final.iqtree` and `.log` for details and software references.

For an optional fresh alignment and tree analysis, install MAFFT, trimAl and IQ-TREE 3, then run:

```bash
bash scripts/reinfer_tree.sh
```

This workflow writes separate results in `generated/reinferred_tree/`. It uses the original alignment/trimming options and recorded IQ-TREE seed, with one thread for portability. Original MAFFT/trimAl versions were not established from the supplied files. This optional workflow has not been run during package preparation; an identical inferred tree or support values are not guaranteed across versions and thread settings. The archived final tree remains the source for the manuscript plot. To rerun IQ-TREE alone with the exact archived trimmed alignment, use that FASTA as `-s`, the recorded model-search/support options, and a separate output prefix.

## Panel F interpretation

The genome is 37,522 bp with GC 43.6597196311% and 71 CDS features. CDS colors are derived from product-name keyword rules in the original script, not directly from a PHROG category field. The exported CDS table records both the product and its final plotting category. Only the first integrase encountered receives the integrase text label, as in the original script.

GC tracks use 500-bp windows, 100-bp steps, truncated terminal windows and omission of windows shorter than 50 bp: 375 windows in total. The gray track, labelled “GC content” in the original figure, actually plots window GC fraction minus genome GC fraction. GC skew is `(G-C)/(G+C)`, with zero when the denominator is zero. These windows do not wrap across the sequence end.

The GenBank record declares a **linear** sequence. The circular drawing is a visualization layout and does not demonstrate a biologically circular molecule. The screenshot uses “Gp142-2_P1”; the source data and original exports use “Gv142-2_P1”. Resolve this display-name inconsistency in the final manuscript without silently renaming sequence identifiers.

## Source tables

* `Figure_E_taxa_sequences.tsv`: tip order, IDs, display names, groups, original/aligned/trimmed protein sequences.
* `Figure_E_nodes.tsv`: node and parent IDs, child ordering, branch lengths, both support statistics and descendant tip IDs. Together these preserve tree structure, not just a list of support values.
* `Figure_E_newick.tsv`: complete original Newick string.
* `Figure_F_CDS.tsv`: all 71 CDS coordinates, strands, products, categories, colors and plot labels. Start/end columns use 1-based inclusive coordinates; the separate BioPython location string uses 0-based half-open coordinates.
* `Figure_F_GC_tracks.tsv`: window bounds (0-based, end-exclusive), midpoint, G/C counts, GC fraction, genome GC, GC deviation and skew. Fractions are on a 0–1 scale.
* `Figure_F_genome.tsv` and `Figure_F_sequence_chunks.tsv`: metadata and ordered nucleotide sequence chunks with 1-based inclusive coordinates.

## Scope of reproduction

The package reproduces the biological content of the standalone plots. Manuscript-specific text placement, label leaders, figure letters, line styling, legend placement and the Gp/Gv display-name edit are not reconstructed. The adapted E script replaces Arial with bundled DejaVu Sans to avoid a proprietary-font dependency. Other E changes are input/output paths. F changes are input/output paths and restriction to the one genome shown. Original scripts remain available for comparison.

## Publication and citation

Before public release, supply the real authors, manuscript title, repository URL, release identifier and author-approved licenses. No authorship, DOI, repository URL or reuse license is invented in this package. Third-party reference sequences remain subject to their source terms. See `docs/publication_notes_zh.md` for the remaining publication metadata and English availability-statement templates. A GitHub release archived in Zenodo provides a fixed, citable DOI; cite the specific archived version used in the paper.

Public accession IDs are retained in the sequence labels. Accession IDs for the newly identified prophage nucleotide sequences, the reconstruction coordinates/source for pGamma, and a complete extraction provenance for every TerL sequence are not established by this figure-only package. Add these from the study records before describing the entire study as independently reproducible. Raw reads, full database installations, unrelated BLAST outputs and other figure analyses are outside this package.
