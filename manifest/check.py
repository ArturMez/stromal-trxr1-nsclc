#!/usr/bin/env python3
"""
The drift alarm: does every manuscript value still come out of the package tables?

WHAT: reads manifest.csv, finds each value in the table that holds it, formats it by the rounding
      rule the manuscript used, and reports every row whose formatted value differs from the one
      printed in the manuscript.
FROM: manifest/manifest.csv and the tables under work_root/output named by its "file" column.
TO:   a report on stdout, and exit code 1 when a row mismatches or cannot be read. A row whose
      expected_diff column names a reason is reported as a known difference instead, and does not
      set the exit code.
WHY:  an analysis change that moves a number by a rounding step is invisible in a diff and shows up
      only when somebody rereads the paper. Run it after any change to a chapter script.

manifest.csv columns
  id        a stable identifier, chapter tag and a number
  chapter   the chapter the value belongs to
  location  where the reader finds it: a paragraph, a figure legend letter or a table
  item      what the value is, in words; where the text states a value in words ("about as
            dense"), the item quotes the words and says which number stands for them
  printed   the value as the manuscript prints it, with no unit and no surrounding words
  file      the table that holds it, relative to work_root/output
  selector  which row and column: "key=value,key=value|column"; the keys are columns of the table
  rule      how the value is formatted, one of
              "<n> dp"                    fixed decimals, for example "2 dp"
              "signed <n> dp"             fixed decimals with a leading plus or minus
              "diff <n> dp"               as signed, but a value that rounds to zero is printed
                                          "0.0" without a sign (the supplementary tables' rule)
              "abs <n> dp"                fixed decimals of the absolute value, for a difference
                                          the text states in words ("13.0 points lower")
              "floor <n> dp"              rounded down, for a lower bound the text states ("larger than
                                          for 98 % of the comparison genes" from a percentile of 98.88)
              "ceiling <n> dp"            rounded up, for an upper bound the text states ("all p < 0.04"
                                          from a largest p of 0.0338)
              "integer"                   a count
              "percent 0 dp"              a value already stored in percent, rounded to whole
              "fraction as percent 0 dp"  a value stored as a fraction, printed in whole percent
              "p: 3 dp, < 0.001 below"    the house rule for a p value
  expected_diff
            empty by default. A short reason where the package knowingly differs from the printed
            value (for example a value that depends on a random stream the package seeds
            differently); check.py then reports the row as a known difference, not a mismatch.

RUN: python3 manifest/check.py [--chapter 1] [--quiet]
"""
import argparse
import re
import sys
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from redox import paths  # noqa: E402  (the repository has to be on the path first)

MANIFEST = Path(__file__).resolve().parent / "manifest.csv"

HEADER = ["id", "chapter", "location", "item", "printed", "file", "selector", "rule", "expected_diff"]


def select(table, selector):
    """The one value a selector names. Raises when it matches no row or more than one."""
    where, _, column = selector.partition("|")
    if not column:
        raise ValueError(f"selector {selector!r} names no column")
    keep = table
    for term in [t for t in where.split(",") if t.strip()]:
        key, _, want = term.partition("=")
        key, want = key.strip(), want.strip()
        if key not in keep.columns:
            raise ValueError(f"{key} is not a column of the table")
        keep = keep[keep[key].astype(str) == want]
    if column.strip() not in keep.columns:
        raise ValueError(f"{column.strip()} is not a column of the table")
    if len(keep) != 1:
        raise ValueError(f"selector {selector!r} matches {len(keep)} rows, expected one")
    return keep.iloc[0][column.strip()]


def apply_rule(value, rule):
    """Format one value the way the manuscript prints it."""
    rule = rule.strip()
    v = float(value)
    if rule.startswith("p:"):
        return "< 0.001" if v < 0.001 else f"{v:.3f}"
    if rule == "integer":
        return str(int(round(v)))
    if rule == "percent 0 dp":
        return f"{v:.0f}"
    if rule == "fraction as percent 0 dp":
        return f"{100 * v:.0f}"
    m = re.fullmatch(r"(signed |diff |abs |floor |ceiling )?(\d+) dp", rule)
    if m:
        n = int(m.group(2))
        if m.group(1) in ("floor ", "ceiling "):      # on the decimal string, so that 0.07 stays 0.07
            mode = ROUND_FLOOR if m.group(1) == "floor " else ROUND_CEILING
            return f"{Decimal(repr(v)).quantize(Decimal(1).scaleb(-n), rounding=mode):.{n}f}"
        if m.group(1) == "abs ":
            return f"{abs(v):.{n}f}"
        if m.group(1) == "diff ":
            return f"{0.0:.{n}f}" if round(v, n) == 0 else f"{v:+.{n}f}"
        return f"{v:+.{n}f}" if m.group(1) else f"{v:.{n}f}"
    raise ValueError(f"unknown rule {rule!r}")


def check(chapter=None, quiet=False):
    """Every row of the manifest, or every row of one chapter. Returns the failing rows.

    A row with an expected_diff reason that does not match is listed as a known difference and is
    not returned; one that matches after all is listed too, so that its reason can be removed.
    """
    rows = pd.read_csv(MANIFEST, dtype=str).fillna("")
    missing = [c for c in HEADER if c not in rows.columns]
    if missing:
        raise SystemExit(f"{MANIFEST} is missing the columns {missing}")
    if chapter is not None:
        rows = rows[rows.chapter.astype(str) == str(chapter)]
    tables, failures, known, stale, checked = {}, [], [], [], 0
    for r in rows.itertuples():
        try:
            if r.file not in tables:
                path = paths.WORK_ROOT / "output" / r.file
                tables[r.file] = pd.read_csv(path)
            got = apply_rule(select(tables[r.file], r.selector), r.rule)
        except Exception as exc:                       # a missing file or a selector that moved
            failures.append((r.id, r.item, r.printed, f"{type(exc).__name__}: {exc}"))
            continue
        checked += 1
        reason = r.expected_diff.strip()
        if got != r.printed.strip():
            (known if reason else failures).append((r.id, r.item, r.printed, got, reason))
        elif reason:
            stale.append((r.id, r.item, r.printed, got, reason))
    if not quiet:
        print(f"{len(rows)} rows in the manifest, {checked} read from the tables")
        for i, item, printed, got, *_ in failures:
            print(f"  {i}  {item}: manuscript {printed!r}, tables {got!r}")
        for i, item, printed, got, reason in known:
            print(f"  known difference  {i}  {item}: manuscript {printed!r}, tables {got!r} ({reason})")
        for i, item, printed, got, reason in stale:
            print(f"  note  {i}: matches although expected_diff says {reason!r}; the reason can go")
    print(f"{'PASS' if not failures else 'FAIL'}: {len(failures)} mismatches, "
          f"{len(known)} known differences")
    return failures


def main():
    """Command line: check all rows or one chapter; the exit code is 1 when a row fails."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chapter", default=None, help="check one chapter only")
    ap.add_argument("--quiet", action="store_true", help="print the verdict only")
    a = ap.parse_args()
    return 1 if check(a.chapter, a.quiet) else 0


if __name__ == "__main__":
    sys.exit(main())
