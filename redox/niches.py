"""
TrxR1 niches: what they are, how large the largest one is, and which cells lie in their area.

WHAT: builds the connected groups of TrxR1-positive stromal cells, measures confluence (the niche
      connectivity of the paper: the fraction of positive cells in the largest group), and marks
      the cells that fall inside the niche area.
FROM: positions in micrometers and a TrxR1-positive flag, one core at a time.
TO:   an integer label per cell, a fraction, and a second integer label per cell.
      Used by chapter4.py (Figure 5, Supplementary Table 3 and Supplementary Figure 6).
WHY:  Figure 5 rests on these three operations. Both radii are inclusive: two cells are linked when
      their distance is at most link_um, and a cell lies in the niche area when it is at most
      zone_um from a niche cell. scipy's query_pairs and query both use "at most", so those two
      calls are used rather than a hand-rolled distance test.
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def link(xy_um, positive, link_um=30.0):
    """Connected groups of TrxR1-positive cells.

    Returns one label per cell of the input: 0, 1, 2 ... for the positive cells, one value per
    group, and -1 for every cell that is not positive. Two positive cells are in the same group
    when a chain of positive cells joins them in steps of at most link_um micrometers.
    """
    positive = np.asarray(positive, dtype=bool)
    labels = np.full(len(positive), -1, dtype=int)
    idx = np.flatnonzero(positive)
    k = len(idx)
    if k == 0:
        return labels
    xy = np.asarray(xy_um, dtype=float)
    pairs = cKDTree(xy[idx]).query_pairs(float(link_um), output_type="ndarray")
    if len(pairs):
        graph = csr_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(k, k))
        labels[idx] = connected_components(graph, directed=False)[1]
    else:
        labels[idx] = np.arange(k)
    return labels


def sizes(labels):
    """Cells per group, indexed by label. An empty array when there is no group."""
    lab = np.asarray(labels)
    lab = lab[lab >= 0]
    return np.bincount(lab) if len(lab) else np.array([], dtype=int)


def confluence(labels, positive):
    """Cells in the largest group divided by all TrxR1-positive cells of the core.

    One when every positive cell belongs to a single group, near zero when the positive cells are
    scattered. Missing when the core has no positive cell.
    """
    positive = np.asarray(positive, dtype=bool)
    k = int(positive.sum())
    if k == 0:
        return float("nan")
    counts = sizes(np.asarray(labels)[positive])
    return float(counts.max() / k) if len(counts) else float("nan")


def area_membership(xy_um, positive, labels, zone_um=15.0):
    """Which niche area each cell of the core lies in.

    Returns one label per cell: the label of the nearest niche cell when that cell is at most
    zone_um micrometers away, and -1 for a cell outside every area. Niche cells are the positive
    cells that carry a group label. TrxR1-negative cells are labeled like any other cell, which is
    what makes the area a piece of tissue rather than a set of positive cells.
    """
    positive = np.asarray(positive, dtype=bool)
    labels = np.asarray(labels)
    out = np.full(len(positive), -1, dtype=int)
    niche = positive & (labels >= 0)
    if not niche.any():
        return out
    xy = np.asarray(xy_um, dtype=float)
    src = np.flatnonzero(niche)
    dist, nearest = cKDTree(xy[src]).query(xy)
    inside = dist <= float(zone_um)
    out[inside] = labels[src[nearest[inside]]]
    return out
