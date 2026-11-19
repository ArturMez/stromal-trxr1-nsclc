"""
The house look of every panel.

WHAT: the lineage and histology colors, the axis labels the panels share, the matplotlib defaults,
      and one saver that writes a panel as PNG and PDF.
FROM: nothing; these are constants.
TO:   work_root/output/chapter<N>/<name>.png and .pdf.
WHY:  colors carry meaning across several figures, so they are defined once, under lineage names.
      Used by chapter1.py, chapter3.py, chapter4.py, tumor_vs_stroma.py and chapter2/06, 07, 09
      and 11. Panel files are individual, one per panel; the assembled multi-panel figures are
      not part of the package.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (the backend must be chosen before pyplot loads)

from . import paths  # noqa: E402

# Lineage colors, keyed by the lineage names of cells.LINEAGES.
LC = {"T": "#1baf7a", "Mac": "#c98f3e", "Fib": "#4a3aa7", "Other": "#9e9d97", "CK": "#d9d8d2"}

# Histology colors: blue squamous, orange adenocarcinoma.
HC = {"LUSC": "#2a78d6", "LUAD": "#eb6834"}

# Histology markers, so that histology does not rest on color alone: squares squamous, circles
# adenocarcinoma. A pooled "combined" estimate keeps its diamond.
HM = {"LUSC": "s", "LUAD": "o"}

# How a cohort is named in a shown label. Code, file names and table columns keep BOMI1 and BOMI2.
COHORT_LABEL = {"BOMI1": "Uppsala I", "BOMI2": "Uppsala II"}

# Ink, secondary ink, muted rules and grid lines.
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"

# The gray of a shuffled or expected value, against which an observed value is drawn.
SHUF = "#b9b8b2"

# How a lineage is named in a label or a legend. US English, spelled out, lower case.
LABEL = {"T": "T cells", "Mac": "macrophages", "Fib": "fibroblasts",
         "Other": "no lineage marker", "CK": "cancer cells"}


def style():
    """The matplotlib defaults of the panels (7 pt text). Call once before drawing a panel."""
    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "text.color": INK,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def clean(ax, left=False):
    """Drop the frame down to the axes a panel needs: bottom, and left only when it carries a scale."""
    for side in ("top", "right") + (() if left else ("left",)):
        ax.spines[side].set_visible(False)


def save(fig, chapter, name):
    """Write a panel as PNG at 300 dpi and as PDF, and return the two paths.

    The axis labels are passed as extra artists of the tight crop: matplotlib measures an axis for
    the crop without the width of a label that is wider than the axes, so such a label would be cut
    at the edge of the file (for example Supplementary Figure 3b and Figure 5h). A panel whose
    labels fit inside the crop is unaffected.
    """
    labels = [lab for ax in fig.axes for lab in (ax.xaxis.label, ax.yaxis.label)
              if lab.get_visible() and lab.get_text()]
    extra = fig.get_default_bbox_extra_artists() + labels   # a given list replaces the default
    written = []
    for ext in ("png", "pdf"):
        path = paths.output(chapter, f"{name}.{ext}")
        fig.savefig(path, dpi=300, facecolor="white", bbox_inches="tight",
                    bbox_extra_artists=extra)
        written.append(path)
    plt.close(fig)
    return written
