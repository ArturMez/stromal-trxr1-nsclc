"""
The shared module behind every chapter script.

WHAT: one definition each of the paths, the coordinate scale, the clinical schema, the per-cell
      table, the survival models, the niche construction, the small statistics, the TrxR1 staining
      (intensity, slide and background), the rank tests of the driver analyses and the house style.
FROM: the tables under work_root/intermediate, themselves built by mif/ (staining also reads the
      TrxR1 channel of the raw cell exports under raw_root).
TO:   nothing on disk; the chapter scripts write the tables and panels.
WHY:  every chapter script needs the same definitions, so each lives here once. Import from here;
      do not copy.

The submodules are imported lazily by name (from redox import survival), so that a script that only
needs paths does not pull in matplotlib.
"""
__all__ = ["paths", "coords", "clinical", "cells", "survival", "niches", "stats", "staining", "ranktests",
           "plotting"]
