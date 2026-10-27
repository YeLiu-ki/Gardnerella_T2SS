#!/usr/bin/env python3

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from pycirclize import Circos


# ============================================================
# Paths
# ============================================================

GBKDIR = ROOT / "data/figure_F"
OUTDIR = ROOT / "generated/figure_F"

OUTDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Genome order
# ============================================================

GENOMES = ["Gv142-2_P1"]


# ============================================================
# Display names
# ============================================================

DISPLAY_NAMES = {
    "vB_Gva_AB1_MW387018.1": "vB_Gva_AB1",
    "Gv88_P1": "Gv88_P1",
    "Gv144-3_P1": "Gv144-3_P1",
    "Gv44089_P1": "Gv44089_P1",
    "pGamma_DNF00622A": "pGamma",
    "Gv142-2_P1": "Gv142-2_P1",
    "Gv88_P2": "Gv88_P2",
}


# ============================================================
# Functional categories
# ============================================================

CATEGORY_COLORS = {
    "DNA packaging": "#E69F00",
    "Head / capsid": "#56B4E9",
    "Tail": "#009E73",
    "Lysis": "#D55E00",
    "Integration": "#CC79A7",
    "Replication / DNA": "#0072B2",
    "Other known": "#999999",
    "Hypothetical": "#DDDDDD",
}


# ============================================================
# Functional classification
# ============================================================

def classify_product(product):

    p = product.lower()

    # --------------------------------------------------------
    # DNA packaging
    # --------------------------------------------------------

    if any(x in p for x in [
        "terminase",
        "dna packaging",
        "packaging protein",
    ]):
        return "DNA packaging"

    # --------------------------------------------------------
    # Head / capsid
    # --------------------------------------------------------

    if any(x in p for x in [
        "capsid",
        "major head",
        "head protein",
        "portal",
        "scaffold",
        "head maturation",
        "head-tail connector",
        "head tail connector",
    ]):
        return "Head / capsid"

    # --------------------------------------------------------
    # Tail
    # --------------------------------------------------------

    if any(x in p for x in [
        "tail protein",
        "major tail",
        "minor tail",
        "tape measure",
        "tail fiber",
        "tail fibre",
        "tail assembly",
        "distal tail",
        "baseplate",
    ]):
        return "Tail"

    # --------------------------------------------------------
    # Lysis
    # --------------------------------------------------------

    if any(x in p for x in [
        "holin",
        "endolysin",
        "lysin",
        "lysozyme",
        "lysis protein",
    ]):
        return "Lysis"

    # --------------------------------------------------------
    # Integration / excision
    # --------------------------------------------------------

    if any(x in p for x in [
        "integrase",
        "recombinase",
        "excisionase",
        "integration protein",
        "excision protein",
    ]):
        return "Integration"

    # --------------------------------------------------------
    # Replication / DNA metabolism
    # --------------------------------------------------------

    if any(x in p for x in [
        "dna polymerase",
        "rna polymerase",
        "helicase",
        "primase",
        "nuclease",
        "exonuclease",
        "endonuclease",
        "dna-binding",
        "dna binding",
        "replication",
        "recombination",
        "single-stranded dna",
        "ssb",
    ]):
        return "Replication / DNA"

    # --------------------------------------------------------
    # Hypothetical
    # --------------------------------------------------------

    if any(x in p for x in [
        "hypothetical",
        "unknown",
        "uncharacterized",
        "unannotated",
    ]):
        return "Hypothetical"

    return "Other known"


# ============================================================
# Key gene labels
#
# Only biologically informative marker genes are labelled.
# TerS is intentionally NOT labelled to reduce crowding.
# ============================================================

def get_short_label(product):

    p = product.lower()

    if "terminase large" in p:
        return "TerL"

    if "integrase" in p:
        return "Integrase"

    if "portal" in p:
        return "Portal"

    if "major capsid" in p:
        return "Major capsid"

    if "major head" in p:
        return "Major head"

    if "tape measure" in p:
        return "Tape measure"

    if "major tail" in p:
        return "Major tail"

    if "endolysin" in p:
        return "Endolysin"

    if "holin" in p:
        return "Holin"

    # Avoid generic "lysin" if it is already endolysin
    if "lysin" in p:
        return "Lysin"

    return None


# ============================================================
# GC content / GC skew
# ============================================================

def sliding_gc(seq, window=500, step=100):

    seq = str(seq).upper()

    genome_gc = gc_fraction(seq)

    positions = []
    gc_deviation = []
    gc_skew = []

    for start in range(0, len(seq), step):

        end = min(
            start + window,
            len(seq)
        )

        fragment = seq[start:end]

        if len(fragment) < 50:
            continue

        g = fragment.count("G")
        c = fragment.count("C")

        gc = (
            (g + c) / len(fragment)
            if len(fragment) > 0
            else 0
        )

        skew = (
            (g - c) / (g + c)
            if (g + c) > 0
            else 0
        )

        positions.append(
            (start + end) / 2
        )

        gc_deviation.append(
            gc - genome_gc
        )

        gc_skew.append(
            skew
        )

    return (
        np.array(positions),
        np.array(gc_deviation),
        np.array(gc_skew),
        genome_gc,
    )


# ============================================================
# Circular genome plotting
# ============================================================

def plot_genome(gbk_file):

    record = SeqIO.read(
        gbk_file,
        "genbank"
    )

    genome_name = gbk_file.stem

    display_name = DISPLAY_NAMES.get(
        genome_name,
        genome_name
    )

    genome_size = len(record.seq)

    print()
    print("======================================")
    print(display_name)
    print(f"Genome length: {genome_size:,} bp")
    print("======================================")

    # --------------------------------------------------------
    # Initialize Circos
    # --------------------------------------------------------

    circos = Circos(
        sectors={
            genome_name: genome_size
        },
        space=1,
    )

    sector = circos.sectors[0]


    # ========================================================
    # Track 1
    # Coordinate axis
    # ========================================================

    axis_track = sector.add_track(
        (94, 100)
    )

    axis_track.axis(
        fc="white",
        ec="black",
        lw=0.7,
    )

    # 5 kb ticks are suitable for ~40–50 kb prophages
    axis_track.xticks_by_interval(
        5000,
        label_formatter=lambda v:
            f"{int(v / 1000)} kb",
        label_size=7,
        label_orientation="vertical",
        tick_length=1.5,
    )


    # ========================================================
    # Track 2
    # CDS arrows
    # ========================================================

    cds_track = sector.add_track(
        (72, 92)
    )

    cds_track.axis(
        fc="white",
        ec="none"
    )

    label_positions = []

    n_cds = 0

    category_counts = {
        x: 0
        for x in CATEGORY_COLORS
    }

    # Only one representative integrase gets a text label
    integrase_labelled = False


    for feature in record.features:

        if feature.type != "CDS":
            continue

        n_cds += 1

        start = int(
            feature.location.start
        )

        end = int(
            feature.location.end
        )

        strand = feature.location.strand

        product = feature.qualifiers.get(
            "product",
            ["hypothetical protein"]
        )[0]

        category = classify_product(
            product
        )

        category_counts[category] += 1

        color = CATEGORY_COLORS[
            category
        ]


        # ----------------------------------------------------
        # Forward-strand CDS
        # ----------------------------------------------------

        if strand == 1:

            cds_track.arrow(
                start,
                end,
                r_lim=(82, 91),
                fc=color,
                ec="white",
                lw=0.35,
                head_length=250,
                shaft_ratio=0.65,
            )


        # ----------------------------------------------------
        # Reverse-strand CDS
        # ----------------------------------------------------

        elif strand == -1:

            cds_track.arrow(
                end,
                start,
                r_lim=(73, 82),
                fc=color,
                ec="white",
                lw=0.35,
                head_length=250,
                shaft_ratio=0.65,
            )


        # ----------------------------------------------------
        # Key gene labels
        # ----------------------------------------------------

        short_label = get_short_label(
            product
        )

        if short_label is None:
            continue


        # Label only one representative integrase
        if short_label == "Integrase":

            if integrase_labelled:
                continue

            integrase_labelled = True


        midpoint = (
            start + end
        ) / 2

        label_positions.append(
            (
                midpoint,
                short_label,
            )
        )


    print(f"CDS: {n_cds}")

    for category, count in category_counts.items():

        print(
            f"  {category:<20} {count}"
        )


    # ========================================================
    # Track 3
    # GC content deviation
    # ========================================================

    (
        x,
        gc_dev,
        gc_skew,
        genome_gc,
    ) = sliding_gc(
        record.seq
    )

    print(
        f"GC content: "
        f"{genome_gc * 100:.2f}%"
    )


    gc_track = sector.add_track(
        (55, 65)
    )

    gc_track.axis(
        fc="#F8F8F8",
        ec="none"
    )


    max_gc = max(
        abs(gc_dev.min()),
        abs(gc_dev.max()),
        0.01,
    )


    gc_track.fill_between(
        x,
        gc_dev,
        0,
        vmin=-max_gc,
        vmax=max_gc,
        color="#666666",
        alpha=0.75,
    )


    # ========================================================
    # Track 4
    # GC skew
    # ========================================================

    skew_track = sector.add_track(
        (42, 52)
    )

    skew_track.axis(
        fc="#F8F8F8",
        ec="none"
    )


    max_skew = max(
        abs(gc_skew.min()),
        abs(gc_skew.max()),
        0.01,
    )


    positive_skew = np.maximum(
        gc_skew,
        0
    )

    negative_skew = np.minimum(
        gc_skew,
        0
    )


    skew_track.fill_between(
        x,
        positive_skew,
        0,
        vmin=-max_skew,
        vmax=max_skew,
        color="#4C78A8",
        alpha=0.80,
    )


    skew_track.fill_between(
        x,
        negative_skew,
        0,
        vmin=-max_skew,
        vmax=max_skew,
        color="#E45756",
        alpha=0.80,
    )


    # ========================================================
    # Center genome name
    # ========================================================

    sector.text(
        (
            f"{display_name}\n"
            f"{genome_size / 1000:.1f} kb\n"
            f"GC {genome_gc * 100:.1f}%"
        ),
        r=15,
        size=12,
        weight="bold",
    )


    # ========================================================
    # Draw Circos
    # ========================================================

    fig = circos.plotfig()

    ax = fig.axes[0]


    # ========================================================
    # Key-gene text labels
    # ========================================================

    for position, label in label_positions:

        theta = (
            2
            * np.pi
            * position
            / genome_size
        )

        angle = np.degrees(
            theta
        )

        # Make labels less likely to appear upside down
        rotation = angle - 90

        ha = "left"

        if 90 < angle < 270:

            rotation += 180

            ha = "right"


        ax.text(
            theta,
            102,
            label,
            fontsize=7.5,
            rotation=rotation,
            rotation_mode="anchor",
            ha=ha,
            va="center",
        )


    # ========================================================
    # Legend
    # ========================================================

    legend_elements = [

        Patch(
            facecolor=color,
            edgecolor="none",
            label=category,
        )

        for category, color
        in CATEGORY_COLORS.items()

    ]


    # Add GC tracks to legend
    legend_elements.extend([

        Patch(
            facecolor="#666666",
            edgecolor="none",
            label="GC content",
        ),

        Patch(
            facecolor="#4C78A8",
            edgecolor="none",
            label="GC skew (+)",
        ),

        Patch(
            facecolor="#E45756",
            edgecolor="none",
            label="GC skew (−)",
        ),

    ])


    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            -0.02
        ),
        ncol=4,
        frameon=False,
        fontsize=8,
        columnspacing=1.3,
        handlelength=1.6,
    )


    # Leave room for legend
    fig.subplots_adjust(
        bottom=0.14
    )


    # ========================================================
    # Save
    # ========================================================

    png = (
        OUTDIR
        / f"{genome_name}_circular.png"
    )

    pdf = (
        OUTDIR
        / f"{genome_name}_circular.pdf"
    )


    fig.savefig(
        png,
        dpi=600,
        bbox_inches="tight",
    )


    fig.savefig(
        pdf,
        bbox_inches="tight",
    )


    plt.close(fig)


    print(
        f"Saved: {png}"
    )

    print(
        f"Saved: {pdf}"
    )


# ============================================================
# Run all genomes
# ============================================================

print(
    "======================================"
)

print(
    "Circular functional genome maps"
)

print(
    "======================================"
)


n_success = 0


for genome in GENOMES:

    gbk = (
        GBKDIR
        / f"{genome}.gbk"
    )


    if not gbk.exists():

        print()
        print(
            "WARNING: GenBank file "
            f"not found: {gbk}"
        )

        continue


    plot_genome(
        gbk
    )

    n_success += 1


# ============================================================
# Finish
# ============================================================

print()
print(
    "======================================"
)

print(
    "Finished"
)

print(
    "======================================"
)

print(
    f"Successfully plotted: "
    f"{n_success}/{len(GENOMES)}"
)

print(
    f"Output directory: "
    f"{OUTDIR}"
)