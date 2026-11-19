#!/usr/bin/env python3
"""
Chapter 2, step 05: the panels of Figure 3 and Supplementary Figures 3 and 4, drawn from the tables.

WHAT: fig2_programs            Figure 3: per cell set (fibroblasts, macrophages, T and NK cells), the
                               programs of fig2_rows.csv with the bar (effect), the line (+-1.96
                               permutation SD), the effect, p and "batches in the same direction", and the
                               bracket of the nested inflammatory CAF row
      suppfig3_programs      Supplementary Figure 3: the heatmap of every primary, contrast and
                               technical-control program of each cell set's panel, three columns, the
                               myCAF row labeled "Matrix CAF signature"
      suppfig4_fibroblast_genes  Supplementary Figure 4: the fibroblast genes that pass the scan rule
                               (04_fibroblast_genes.py) and every tested complement member, one dot per batch
ITEMS: Figure 3, Supplementary Figure 3, Supplementary Figure 4 (the gene tables of Supplementary Figure 4 are written by
      04_fibroblast_genes.py).
FROM: work_root/output/chapter2/fig2_rows.csv, program_results.csv, gene_programs_definitions.csv
      (Figure 3 and Supplementary Figure 3); suppfig4_rows.csv, suppfig4_dots.csv, suppfig4_numbers.csv
      (Supplementary Figure 4). Every panel is drawn from these tables re-read from disk.
TO:   work_root/output/chapter2/<panel>.png (300 dpi for Figure 3, 200 dpi for the supplementary figures) and .pdf
WHY:  the figures are drawn from the saved tables only, so that a figure can be redrawn without rerunning a test.
RUN:  python3 chapter2/05_plots.py [--only fig2 suppfig3 suppfig4]
"""
import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402

# SURF is the surface color of Supplementary Figures 3 and 4 (figure, axes, empty heatmap cells and dot outlines); white, so
# the figures have no background.
BLUE, RED, INK, INK2, MUTED, GRID, SURF, AXIS = ("#2a78d6", "#e34948", "#0b0b0b", "#52514e", "#898781", "#e1e0d9",
                                                 "#ffffff", "#c3c2b7")
MINUS = "−"
# Figure 3 bars: pastel colors as in Figure 1 (the tile colors of Figure 1c, where blue marks longer and red shorter
# survival); here blue marks a higher and red a lower program score.
BAR_UP, BAR_DOWN = "#9dbde3", "#e99391"
NOTE = ("Gene lists and references: Supplementary Data 1. Effect = mean paired difference in the program "
        "score (mean within-batch z of its genes), TXNRD1-positive minus depth-matched cells without detected TXNRD1.")


def read(name, **kw):
    """Read an output table of chapter 2 from disk (floats parsed round-trip)."""
    return pd.read_csv(L.out(name), float_precision="round_trip", **kw)


def save(fig, name, dpi=300, **kw):
    """Save a figure as <name>.png and <name>.pdf in the chapter 2 output folder and close it."""
    for ext in ("png", "pdf"):
        fig.savefig(L.out(f"{name}.{ext}"), dpi=dpi, **kw)
    plt.close(fig)
    print(f"wrote {name}.png and .pdf")


# ---------------------------------------------------------------- Figure 3
def fig2():
    """Figure 3 (file fig2_programs): bars, permutation bands, p and batch votes of the programs of fig2_rows.csv."""
    F = read("fig2_rows.csv")
    plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "text.color": INK})
    blocks = [(g.cell_block.iloc[0], int(g.denominator.iloc[0].split(" ")[0].replace(",", "")), g.sort_values("row"))
              for _, g in F.groupby("cell_set", sort=False)]
    # every text is at least 7 pt at the printed width (7.2 in, within the 190 mm of a full-width Elsevier figure;
    # Elsevier asks 7 pt for normal text); positions below are in inches on the figure
    W, ROW, HEAD, GAP, TOP, BOT, NOTEH = 7.2, 0.19, 0.30, 0.02, 0.46, 0.03, 0.09
    H = TOP + BOT + sum(HEAD + ROW / 2 + len(d) * ROW + GAP + (NOTEH if (d.nested_under.fillna("") != "").any() else 0)
                        for _, _, d in blocks)
    XN, XB0, XB1, XV, XP, XS, SC = 2.85, 3.05, 4.95, 5.10, 5.85, 6.55, 0.10
    F_HEAD, F_CELL, F_LAB, F_VAL, F_NOTE = 7.5, 9.0, 8.0, 7.5, 7.0

    def ef(e):
        """x position of an effect e on the bar axis (clipped at +-SC)."""
        return (XB0 + XB1) / 2 + np.clip(e, -SC, SC) / SC * (XB1 - XB0) / 2

    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    rend = fig.canvas.get_renderer()
    ax.text((XB0 + XB1) / 2, TOP - 0.06, "effect under the TXNRD1 split\n(bar, line = ±1.96 permutation SD)",
            ha="center", va="bottom", fontsize=F_HEAD, color=MUTED)
    ax.text(XP, TOP - 0.06, "p", ha="center", va="bottom", fontsize=F_HEAD, color=MUTED)
    ax.text(XS, TOP - 0.06, "batches in the\nsame direction", ha="center", va="bottom", fontsize=F_HEAD, color=MUTED)
    nested = L.panels()["nested"]
    y = TOP
    for cell, pairs, d in blocks:
        y += HEAD / 2
        ax.text(0.05, y, cell, fontsize=F_CELL, fontweight="bold", va="center")
        ax.text(XN, y, f"{pairs:,} pairs", fontsize=F_HEAD, color=MUTED, va="center", ha="right")
        ax.plot([0.05, W - 0.05], [y + HEAD / 2 - 0.03] * 2, color=GRID, lw=0.6)
        y += HEAD / 2 + ROW / 2
        y0 = y
        left = {}
        for r in d.itertuples():
            t = ax.text(XN, y, r.label_as_printed, ha="right", va="center", fontsize=F_LAB)
            left[r.set_id] = (y, XN - t.get_window_extent(rend).width / fig.dpi)
            x0, x1 = ef(0), ef(r.effect)
            ax.add_patch(Rectangle((min(x0, x1), y - 0.05), abs(x1 - x0), 0.10,
                                   color=BAR_UP if r.effect > 0 else BAR_DOWN, lw=0, zorder=2))
            ax.plot([ef(r.band95_lo), ef(r.band95_hi)], [y, y], color=INK, lw=0.7, zorder=3)
            ax.text(XV, y, f"{r.effect:+.3f}", va="center", fontsize=F_VAL, color=INK2)
            ax.text(XP, y, "< 0.001" if r.p_gauss < 0.001 else f"{r.p_gauss:.3f}", ha="center", va="center",
                    fontsize=F_VAL, color=INK2)
            ax.text(XS, y, r.study_vote, ha="center", va="center", fontsize=F_VAL, color=INK2)
            step = ROW
            if r.set_id in nested:
                parent, note = nested[r.set_id]["parent"], nested[r.set_id]["note"]
                ax.text(XN, y + 0.11, note, ha="right", va="center", fontsize=F_NOTE, color=MUTED, style="italic")
                py, px = left[parent]
                xb = min(px, left[r.set_id][1]) - 0.10
                ax.plot([xb, xb], [py, y], color=AXIS, lw=0.8, zorder=1)
                for yy in (py, y):
                    ax.plot([xb, xb + 0.04], [yy, yy], color=AXIS, lw=0.8, zorder=1)
                step += NOTEH
            y += step
        ax.plot([ef(0)] * 2, [y0 - ROW / 2, y - ROW / 2], color=AXIS, lw=0.8, zorder=1)
        y += GAP
    save(fig, "fig2_programs", transparent=True)


# ---------------------------------------------------------------- Supplementary Figure 3
def suppfig3():
    """Supplementary Figure 3 (file suppfig3_programs): heatmap of the program effects per cell set; programs that the
    library-size split reproduces are marked "RNA amount" instead of drawn."""
    plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 9.5, "text.color": INK})
    cfg = L.panels()
    RS = read("program_results.csv")
    D = pd.read_csv(L.out("gene_programs_definitions.csv"), dtype=str).fillna("")
    RS["stem"] = RS.set_id.map(L.stem)
    RS["kind"] = np.select([RS.verdict.str.startswith("specific"), RS.verdict.str.startswith(("withdrawn", "not specific"))],
                           ["specific", "rna"], "none")
    lab2cell = {s["label"]: c for c, s in L.CELL_SETS.items()}
    RS["cell"] = RS.cell_type.map(lab2cell)
    shown = RS[RS.role.isin(cfg["suppfig3_roles"])]
    shown = shown[[s in cfg["panel"][c] for s, c in zip(shown.stem, shown.cell)]]
    cells = list(L.CELL_SETS)
    pairs = {c: int(RS[RS.cell == c].denominator.iloc[0].split(" ")[0].replace(",", "")) for c in cells}
    prim = D[D.role.isin(cfg["suppfig3_roles"])].assign(stem=lambda x: x.set_id.map(L.stem))
    name_of = {s: cfg["printed_labels"].get(i, L.us(p)) for s, i, p in zip(prim.stem, prim.set_id, prim.program)}
    rows = []
    for g in cfg["groups"]:
        stems = list(dict.fromkeys(shown[shown.group == g].sort_values("stem").stem))
        keep = [s for s in stems if (shown[shown.stem == s].kind != "rna").any()]
        if keep:
            rows.append(("h", g))
            rows += [("r", s) for s in keep]
    nB = len(rows)
    fig, ax = plt.subplots(figsize=(10.4, 3.4 + 0.29 * nB))
    fig.patch.set_facecolor(SURF)
    ax.set_facecolor(SURF)
    fig.subplots_adjust(left=0.01, right=0.99, top=1 - 0.95 / fig.get_figheight(), bottom=1.5 / fig.get_figheight())
    ax.set_xlim(0, 10.4)
    ax.set_ylim(nB, -2.6)
    ax.axis("off")
    XN, CW = 4.6, 1.35
    x0 = XN + 0.25
    NORM = Normalize(-0.08, 0.08)
    cmap = LinearSegmentedColormap.from_list("d", [RED, "#f4f3f0", BLUE])
    for j, c in enumerate(cells):
        xc = x0 + j * CW + CW / 2
        ax.text(xc, -1.9, textwrap.fill(L.CELL_SETS[c]["label"], 13), ha="center", va="center", fontsize=9.5,
                fontweight="bold")
        ax.text(xc, -0.85, f"{pairs[c]:,} pairs", ha="center", va="center", fontsize=8, color=MUTED)
    for i, (kind, val) in enumerate(rows):
        if kind == "h":
            ax.text(0.05, i + 0.1, val, fontsize=9.5, fontweight="bold", va="center", color=INK2)
            ax.plot([0.05, x0 + len(cells) * CW], [i - 0.45, i - 0.45], color=GRID, lw=0.8)
            continue
        ax.text(XN, i, name_of.get(val, val), ha="right", va="center", fontsize=9.2)
        for j, c in enumerate(cells):
            xl = x0 + j * CW
            r = shown[(shown.cell == c) & (shown.stem == val)]
            if r.empty:
                ax.add_patch(Rectangle((xl + 0.03, i - 0.42), CW - 0.06, 0.84, facecolor=SURF, edgecolor=GRID, lw=0.6))
                continue
            r = r.iloc[0]
            if r.kind == "rna":
                ax.add_patch(Rectangle((xl + 0.03, i - 0.42), CW - 0.06, 0.84, facecolor="#e8e7e2", edgecolor=GRID, lw=0.6))
                ax.text(xl + CW / 2, i, "RNA amount", ha="center", va="center", fontsize=7, color=AXIS, style="italic")
                continue
            ax.add_patch(Rectangle((xl + 0.03, i - 0.42), CW - 0.06, 0.84, facecolor=cmap(NORM(r.effect)),
                                   edgecolor=GRID, lw=0.6))
            ax.text(xl + CW / 2, i, f"{r.effect:+.3f}".replace("-", MINUS), ha="center", va="center", fontsize=8.5,
                    fontweight="bold" if abs(r.z) >= L.Z95 else "normal",
                    color=INK if abs(NORM(r.effect) - 0.5) < 0.32 else "white")
    fig.text(0.02, 0.985, "Programs across cell types", fontsize=14, fontweight="bold", va="top")
    fig.text(0.02, 0.985 - 0.30 / fig.get_figheight(),
             "Effect of the TXNRD1 split on each program score, in within-batch z. Bold where |z| is 1.96 or more.",
             fontsize=9.5, color=MUTED, va="top")
    cax = fig.add_axes([0.04, 0.86 / fig.get_figheight(), 0.20, 0.11 / fig.get_figheight()])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([-0.08, 0, 0.08])
    cb.set_ticklabels([f"{MINUS}0.08 lower", "0", "+0.08 higher"])
    cb.outline.set_visible(False)
    cax.tick_params(labelsize=8, length=0)
    foot = ["“RNA amount” marks a program the library-size split moves the same way, which is reported rather "
            "than drawn as an effect. An empty cell means the program was scored in that cell type but is not part of its "
            "panel. Color is clipped at ±0.08.", NOTE]
    fig.text(0.02, 0.01, "\n".join(textwrap.fill(t, 165) for t in foot), fontsize=7.6, color=MUTED, va="bottom")
    save(fig, "suppfig3_programs", dpi=200, facecolor=SURF)


# ---------------------------------------------------------------- Supplementary Figure 4
def suppfig4():
    """Supplementary Figure 4 (file suppfig4_fibroblast_genes): the fibroblast genes, one dot per batch, pooled effect and band."""
    plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 9.5, "text.color": INK})
    Rw = read("suppfig4_rows.csv", keep_default_na=False, na_values=[""]).sort_values("row")
    dots = read("suppfig4_dots.csv")
    Nm = read("suppfig4_numbers.csv", keep_default_na=False).set_index("item")
    drawn = dots[dots.drawn].drop_duplicates("block").sort_values("column")
    cols, lab, npairs = list(drawn.block), dict(zip(drawn.block, drawn.label)), dict(zip(drawn.block, drawn.n_pairs))
    ncol, nr1 = len(cols), len(Rw)
    maxn = max(npairs.values())

    def psize(k):
        """Dot area for a batch with k pairs."""
        return 20 + k / maxn * 300

    fig, ax = plt.subplots(figsize=(11.5, 3.9 + 0.36 * nr1))
    fig.patch.set_facecolor(SURF)
    ax.set_facecolor(SURF)
    fig.subplots_adjust(left=0.01, right=0.99, top=1 - 0.95 / fig.get_figheight(), bottom=1.6 / fig.get_figheight())
    X0, dxs = 1.9, 0.62
    XP0 = X0 + ncol * dxs + 0.35
    XP1 = XP0 + 1.8
    XA = XP1 + 0.75
    SCB = 0.2
    ax.set_xlim(0, XA + 0.6)
    ax.set_ylim(nr1 - 0.4, -3.6)
    ax.axis("off")
    cmap = LinearSegmentedColormap.from_list("d", [RED, "#f0efec", BLUE])
    norm_g = Normalize(-0.4, 0.4)
    for j, b in enumerate(cols):
        ax.text(X0 + j * dxs, -0.75, f"{lab[b]} ({npairs[b]})", rotation=45, ha="left", va="bottom", fontsize=8,
                color=INK2, rotation_mode="anchor")

    def efb(e):
        """x position of an effect e on the pooled-effect axis (clipped at the scale ends)."""
        return (XP0 + XP1) / 2 + np.clip(e, -SCB, SCB) / SCB * (XP1 - XP0) / 2

    ax.text((XP0 + XP1) / 2, -0.75, "pooled effect\n(±1.96 permutation SD)", ha="center", va="bottom", fontsize=8.5,
            color=MUTED)
    ax.text(XA, -0.75, "batches\nsame way", ha="center", va="bottom", fontsize=8.5, color=MUTED)
    for i, r in enumerate(Rw.itertuples()):
        if r.kind == "header":
            ax.text(0.05, i, r.section, fontsize=10, fontweight="bold", va="center", color=INK2)
            ax.plot([0.05, XA + 0.5], [i + 0.45, i + 0.45], color=AXIS, lw=0.8)
            continue
        ax.text(X0 - 0.35, i, r.gene, ha="right", va="center", fontsize=10)
        ax.plot([X0 - 0.2, XA + 0.5], [i + 0.5, i + 0.5], color=GRID, lw=0.6, zorder=0)
        vals = dots[dots.gene == r.gene].set_index("block").effect
        for j, b in enumerate(cols):
            v = vals.get(b, np.nan)
            if not np.isfinite(v):
                ax.scatter([X0 + j * dxs], [i], s=4, color=AXIS)
                continue
            ax.scatter([X0 + j * dxs], [i], s=psize(npairs[b]), color=cmap(norm_g(v)), edgecolor=SURF, lw=1, zorder=3)
        x0, x1 = efb(0), efb(r.eff)
        ax.add_patch(Rectangle((min(x0, x1), i - 0.15), abs(x1 - x0), 0.30, color=BLUE if r.eff > 0 else RED, lw=0, zorder=2))
        ax.plot([efb(r.lo), efb(r.hi)], [i, i], color=INK, lw=0.9, zorder=3)
        ax.text(XA, i, f"{int(r.batches_same_way)} of {int(r.batches_drawn)}", ha="center", va="center", fontsize=9.5,
                color=INK2)
    ax.plot([efb(0)] * 2, [0.5, nr1 - 0.5], color=AXIS, lw=1, zorder=1)
    fig.text(0.02, 0.985, "Fibroblast genes, batch by batch", fontsize=14, fontweight="bold", va="top")
    fig.text(0.02, 0.985 - 0.30 / fig.get_figheight(),
             "Genes chosen by rule, not by rank. Each dot is one sequencing batch: color = effect in that batch, size = matched pairs.",
             fontsize=9.5, color=MUTED, va="top")
    cax = fig.add_axes([0.05, 0.95 / fig.get_figheight(), 0.22, 0.13 / fig.get_figheight()])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm_g, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([-0.4, 0, 0.4])
    cb.set_ticklabels([f"{MINUS}0.4 lower", "0", "+0.4 higher"])
    cb.outline.set_visible(False)
    cax.tick_params(labelsize=8, length=0)
    fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=np.sqrt(psize(k)), color=MUTED, label=f"{k} pairs")
                        for k in (25, 100, 250)],
               loc="lower left", bbox_to_anchor=(0.33, 0.8 / fig.get_figheight()), ncol=3, frameon=False, fontsize=8.3)
    removed = Nm.loc["genes_removed_rna", "genes"].split()
    foot = [f"Design: {int(Nm.loc['pairs_all', 'value']):,} depth-matched pairs from "
            f"{int(Nm.loc['cells', 'value']):,} fibroblasts (at most 2,000 per batch and histology). Gene effects and q from "
            f"the pre-registered gene scan (Benjamini-Hochberg over {int(Nm.loc['genes_scanned', 'value']):,} genes). A gene is "
            "listed when q is 0.20 or below and the scan's library-size split does not reproduce it (same direction, |z| "
            f"1.96 or more). Removed by that rule: {', '.join(removed) or 'none'}.",
            f"Batches with fewer than 20 pairs not drawn: {Nm.loc['batches_not_drawn', 'genes'] or 'none'}. Color clipped at "
            "the ends of each scale. Unit: within-batch z of the count fraction (genes) or of the program score."]
    fig.text(0.02, 0.01, "\n".join(textwrap.fill(t, 178) for t in foot), fontsize=7.6, color=MUTED, va="bottom")
    save(fig, "suppfig4_fibroblast_genes", dpi=200, facecolor=SURF)


def main():
    """Draw the panels named with --only."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["fig2", "suppfig3", "suppfig4"],
                    choices=["fig2", "suppfig3", "suppfig4"])
    a = ap.parse_args()
    for name in a.only:
        {"fig2": fig2, "suppfig3": suppfig3, "suppfig4": suppfig4}[name]()


if __name__ == "__main__":
    main()
