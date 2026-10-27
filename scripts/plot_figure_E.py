#!/usr/bin/env python3

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re


# ============================================================
# 1. Paths
# ============================================================

TREE_FILE = Path(
    str(ROOT / "data/figure_E/TerL_final.treefile")
)

OUTDIR = Path(
    str(ROOT / "generated/figure_E")
)

OUTDIR.mkdir(parents=True, exist_ok=True)

PNG_OUT = OUTDIR / "TerL_cladogram_rectangular.png"
PDF_OUT = OUTDIR / "TerL_cladogram_rectangular.pdf"


# ============================================================
# 2. Sequence categories
# ============================================================

# Five prophages identified in this study
CANDIDATES = {
    "Gv88_P1",
    "Gv88_P2",
    "Gv142-2_P1",
    "Gv144-3_P1",
    "Gv44089_P1",
}

# Published / reconstructed references
REFERENCES = {
    "vB_Gva_AB1_MW387018.1",
    "pGamma_DNF00622A",
}


# ============================================================
# 3. Display labels
# ============================================================

LABELS = {

    # --------------------------------------------------------
    # Candidate prophages
    # --------------------------------------------------------

    "Gv88_P1":
        "Gv88_P1",

    "Gv88_P2":
        "Gv88_P2",

    "Gv142-2_P1":
        "Gv142-2_P1",

    "Gv144-3_P1":
        "Gv144-3_P1",

    "Gv44089_P1":
        "Gv44089_P1",

    # --------------------------------------------------------
    # Reference sequences
    # --------------------------------------------------------

    "vB_Gva_AB1_MW387018.1":
        "vB_Gva_AB1 (MW387018.1)",

    "pGamma_DNF00622A":
        "pGamma (DNF00622A)",

    # --------------------------------------------------------
    # Public TerL homologues
    # --------------------------------------------------------

    "WP_410311434.1":
        "G. vaginalis | WP_410311434.1",

    "WP_481632921.1":
        "G. vaginalis | WP_481632921.1",

    "WP_337978145.1":
        "S. anginosus | WP_337978145.1",

    "WP_075038564.1":
        "G. vaginalis | WP_075038564.1",

    "WP_064340359.1":
        "Gardnerella spp. | WP_064340359.1",

    "WP_048730474.1":
        "G. vaginalis | WP_048730474.1",

    "WP_086114633.1":
        "A. macacae | WP_086114633.1",

    "WP_071414771.1":
        "Alloscardovia sp. | WP_071414771.1",

    "WP_335340446.1":
        "A. omnicolens | WP_335340446.1",

    "WP_241513427.1":
        "B. amazonense | WP_241513427.1",

    "MDU6662587.1":
        "Actinomyces sp. | MDU6662587.1",

    "DAU07354.1":
        "Caudoviricetes sp. | DAU07354.1",
}


# ============================================================
# 4. Read IQ-TREE tree
# ============================================================

tree = Phylo.read(
    TREE_FILE,
    "newick"
)

# IMPORTANT:
#
# The topology inferred by IQ-TREE is retained.
#
# No midpoint rooting.
# No outgroup rooting.
#
# Original evolutionary branch lengths are NOT used
# in this topology-only cladogram.


# ============================================================
# 5. Parse SH-aLRT / UFBoot support
# ============================================================

def get_support(clade):

    if clade.is_terminal():
        return None

    if clade.name is None:
        return None

    text = str(clade.name).strip()

    match = re.fullmatch(
        r"([0-9.]+)\s*/\s*([0-9.]+)",
        text
    )

    if match is None:
        return None

    sh_alrt = float(match.group(1))
    ufboot = float(match.group(2))

    return sh_alrt, ufboot


def fmt_support(value):

    if abs(value - round(value)) < 1e-8:
        return str(int(round(value)))

    return f"{value:.1f}"


# ============================================================
# 6. Calculate topology depth
#
# Root/representation point = depth 0
# Every internal edge = one topological step.
#
# This is NOT evolutionary distance.
# ============================================================

topological_depth = {}


def assign_depth(clade, depth=0):

    topological_depth[clade] = depth

    for child in clade.clades:

        assign_depth(
            child,
            depth + 1
        )


assign_depth(tree.root)


# Maximum internal depth
max_depth = max(
    topological_depth.values()
)


# ============================================================
# 7. X coordinates
#
# Internal nodes:
#     topology depth
#
# Terminal nodes:
#     all aligned to the same right-hand position
# ============================================================

TIP_X = max_depth + 1.5

xpos = {}


for clade in tree.find_clades():

    if clade.is_terminal():

        xpos[clade] = TIP_X

    else:

        xpos[clade] = topological_depth[clade]


# ============================================================
# 8. Y coordinates
#
# Tips are evenly spaced.
# Internal nodes are centered between first and last child.
# ============================================================

terminals = tree.get_terminals()

N = len(terminals)

print(f"Number of terminal sequences: {N}")


# Top-to-bottom order
ypos = {
    terminal: i
    for i, terminal in enumerate(terminals)
}


def assign_internal_y(clade):

    if clade.is_terminal():

        return ypos[clade]

    child_y = [
        assign_internal_y(child)
        for child in clade.clades
    ]

    y = (
        child_y[0]
        + child_y[-1]
    ) / 2.0

    ypos[clade] = y

    return y


assign_internal_y(tree.root)


# ============================================================
# 9. Matplotlib settings
# ============================================================

plt.rcParams.update({

    "font.family": "DejaVu Sans",

    "font.size": 9,

    # Editable text in Illustrator
    "pdf.fonttype": 42,

    "ps.fonttype": 42,

})


fig, ax = plt.subplots(
    figsize=(10.5, 7.5)
)


# ============================================================
# 10. Draw topology manually
#
# All branches black.
# No branch-length scaling.
# ============================================================

def draw_clade(clade):

    if clade.is_terminal():
        return

    parent_x = xpos[clade]

    children = clade.clades

    child_ys = [
        ypos[child]
        for child in children
    ]


    # --------------------------------------------------------
    # Vertical connector
    # --------------------------------------------------------

    ax.plot(
        [parent_x, parent_x],
        [
            min(child_ys),
            max(child_ys)
        ],
        color="black",
        linewidth=0.9,
        solid_capstyle="round",
        zorder=1,
    )


    # --------------------------------------------------------
    # Horizontal child branches
    # --------------------------------------------------------

    for child in children:

        child_x = xpos[child]

        child_y = ypos[child]

        ax.plot(
            [
                parent_x,
                child_x
            ],
            [
                child_y,
                child_y
            ],
            color="black",
            linewidth=0.9,
            solid_capstyle="round",
            zorder=1,
        )

        draw_clade(child)


draw_clade(tree.root)


# ============================================================
# 11. Draw terminal labels
# ============================================================

LABEL_OFFSET = 0.12


for terminal in terminals:

    raw_name = terminal.name

    label = LABELS.get(
        raw_name,
        raw_name
    )

    x = xpos[terminal]
    y = ypos[terminal]


    # --------------------------------------------------------
    # Candidate prophages
    # --------------------------------------------------------

    if raw_name in CANDIDATES:

        display = label

        fontsize = 9.5

        fontweight = "bold"


    # --------------------------------------------------------
    # Reference sequences
    # --------------------------------------------------------

    elif raw_name in REFERENCES:

        display =  label

        fontsize = 9.5

        fontweight = "bold"


    # --------------------------------------------------------
    # Public homologues
    # --------------------------------------------------------

    else:

        display = label

        fontsize = 8.5

        fontweight = "normal"


    ax.text(
        x + LABEL_OFFSET,
        y,
        display,
        ha="left",
        va="center",
        fontsize=fontsize,
        fontweight=fontweight,
        color="black",
        zorder=3,
    )


# ============================================================
# 12. Draw support labels
#
# Only show nodes where:
#
# SH-aLRT >= 80
# AND
# UFBoot >= 80
# ============================================================

for clade in tree.get_nonterminals():

    support = get_support(clade)

    if support is None:
        continue

    sh_alrt, ufboot = support


    # Strict support threshold
    if not (
        sh_alrt >= 80
        and
        ufboot >= 80
    ):
        continue


    x = xpos[clade]

    y = ypos[clade]


    support_text = (
        f"{fmt_support(sh_alrt)}/"
        f"{fmt_support(ufboot)}"
    )


    # Place support slightly above/right of node
    ax.text(
        x + 0.08,
        y - 0.18,
        support_text,
        ha="left",
        va="bottom",
        fontsize=7,
        color="black",
        zorder=4,

        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.90,
            "pad": 0.3,
        },
    )


# ============================================================
# 13. Optional symbols at terminal branch ends
#
# Small symbols improve identification without using colors.
# ============================================================

for terminal in terminals:

    raw_name = terminal.name

    x = xpos[terminal]
    y = ypos[terminal]


    if raw_name in CANDIDATES:

        ax.scatter(
            x,
            y,
            s=16,
            marker="o",
            facecolor="black",
            edgecolor="black",
            linewidth=0.5,
            zorder=5,
        )


    elif raw_name in REFERENCES:

        ax.scatter(
            x,
            y,
            s=18,
            marker="D",
            facecolor="black",
            edgecolor="black",
            linewidth=0.5,
            zorder=5,
        )


# ============================================================
# 14. Clean axes
# ============================================================

ax.set_xlabel("")
ax.set_ylabel("")

ax.set_xticks([])
ax.set_yticks([])


for spine in [
    "top",
    "right",
    "bottom",
    "left",
]:

    ax.spines[spine].set_visible(False)


# ============================================================
# 15. Orientation
#
# First terminal appears at top.
# ============================================================

ax.invert_yaxis()


# ============================================================
# 16. Limits
#
# Leave enough room on right for long labels.
# ============================================================

ax.set_xlim(
    -0.3,
    TIP_X + 5.2
)

ax.set_ylim(
    N - 0.3,
    -0.8
)


# ============================================================
# 17. Legend
# ============================================================

legend_handles = [

    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="black",
        markeredgecolor="black",
        markersize=6,
        label="Prophage identified in this study",
    ),

    Line2D(
        [0],
        [0],
        marker="D",
        linestyle="None",
        markerfacecolor="black",
        markeredgecolor="black",
        markersize=5,
        label="Published/reconstructed reference",
    ),
]


ax.legend(
    handles=legend_handles,

    loc="lower left",

    bbox_to_anchor=(
        0.0,
        -0.065
    ),

    frameon=False,

    fontsize=8.5,

    ncol=2,

    handletextpad=0.6,

    columnspacing=1.5,
)


# ============================================================
# 18. No title inside plot
#
# Better for manuscript multi-panel figures.
# ============================================================

ax.set_title("")


# ============================================================
# 19. Layout
# ============================================================

plt.subplots_adjust(
    left=0.035,
    right=0.985,
    top=0.985,
    bottom=0.11,
)


# ============================================================
# 20. Save
# ============================================================

fig.savefig(
    PNG_OUT,
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    PDF_OUT,
    bbox_inches="tight",
)

plt.close(fig)


# ============================================================
# 21. Report
# ============================================================

print()
print("DONE")
print()
print(f"PNG: {PNG_OUT}")
print(f"PDF: {PDF_OUT}")
print()
print(
    "Topology-only rectangular cladogram generated."
)
print(
    "All terminal taxa are right-aligned."
)
print(
    "Original IQ-TREE branch lengths were not used."
)
print(
    "Only nodes with SH-aLRT >=80 AND UFBoot >=80 "
    "are labelled."
)