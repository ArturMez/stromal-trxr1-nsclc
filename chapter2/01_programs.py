#!/usr/bin/env python3
"""
Chapter 2, step 01: the gene programs (gene lists) and the gene-program table.

WHAT: s1 builds every program of the catalog (CAT) from its source: MSigDB sets (Reactome, GO,
      Hallmark, KEGG legacy, NABA, Travaglini lung atlas, curated) and paper-derived lists. The catalog
      holds the programs that a script of this folder reads. s1 applies the rule-based edits
      (immunoglobulin genes, collagen chains, B2M), reconciles symbols to the atlas gene axis through
      HGNC (previous and alias symbols, replaced when exactly one match exists), and writes the gmt
      with a definitions table. It also writes the list-in-list check behind the nested programs of
      Figure 3.
      s2 writes the gene-program table behind the Supplementary Figure 3 footnote: program id, label
      as printed, source (identifier, reference, PMID, DOI, database version), every gene, the edits
      made, and the genes tested in each cell set (measured in at least 60 % of the cells and detected
      in at least 2 % of them, TXNRD1 removed). 08_supptable2.py s2 builds the sheet Gene programs of
      Supplementary Data 1 from it.
ITEMS: Figure 3 and Supplementary Figures 3 and 4 (the programs and their labels); Supplementary Data 1
      (sheet Gene programs, through 08_supptable2.py and 10_supplementary_data1.py).
FROM: chapter2/data/sources/ (shipped with the package, see README):
          msigdb_v2026.1.Hs_resource_sets.json.gz   the MSigDB API responses for the database sets
                                                    (release v2026.1.Hs) and the Reactome release 97
                                                    cross-check
          literature_lists.json.gz                   the paper-derived lists: two independent
                                                    extractions of each list from its source paper,
                                                    and an adjudication where they differed
          pubmed_hgnc_responses.json.gz              PubMed and HGNC responses, so no network is needed
      the gene axis of the atlas object (configuration key "atlas_object"); s2 also reads the working
      sets written by 02_working_sets.py.
TO:   work_root/output/chapter2/
          gene_programs.gmt               set id, "label | role | first reference", genes
          gene_programs_definitions.csv   one row per set: role, group, source, reference, PMID, DOI,
                                          version, genes in the source, genes removed with the reason,
                                          genes used, symbols renamed (+ .md)
          gene_programs_nested.csv        the list-in-list check behind the nested programs of Figure 3 (+ .md)
          supptable_gene_programs.csv     s2, every list with its source, genes and tested genes per cell set
                                          (+ .md); 08_supptable2.py s2 builds Supplementary Data 1 (Gene programs) from it
WHY:  the lists are shipped with the package as data, so that the programs are rebuilt without a
      network and without the original downloads. Labels are in US spelling (the ids are the keys that
      every table joins on).
RUN:  python3 chapter2/01_programs.py [--only s1 s2] [--allow-network]
      s1 needs no other step (it reads the gene axis of the atlas object or, failing that, the working sets).
      s2 needs the working sets of 02_working_sets.py. --allow-network fetches PubMed or HGNC entries that
      are missing from the shipped responses.
"""
import argparse
import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ch2lib as L  # noqa: E402
from redox import paths  # noqa: E402

SRC = L.DATA / "sources"
LOG = L.Log()

# the resources' own reference papers (PubMed ids)
RESOURCE_PAPER = {"Reactome": "37941124", "GO": "36866529", "Hallmark": "26771021", "KEGG": "10592173",
                  "NABA": "22159717", "Travaglini": "33208946", "MSigDB": "16199517"}

# ---------------------------------------------------------------- editing rules (applied to the genes of a list)
IG = re.compile(r"^(IG[HKL][VDJ]\S*|IGH[GAMDE]\d*|IGKC|IGLC\d*|IGLL\d*)$")
NONHUMAN = re.compile(r"^([0-9].*|[A-Z]|.*[a-z].*)$")
RULES = {
    "immunoglobulin": (lambda g: bool(IG.match(g)),
                       "immunoglobulin genes: outside B and plasma cells their reads are ambient RNA"),
    "duplicate_locus_symbol": (lambda g: g == "C4B_2", "C4B_2 duplicates the C4B locus"),
    "collagen_chains": (lambda g: bool(re.match(r"^COL\d+A\d+$", g)),
                        "collagen chain genes: the set keeps the genes that modify, process and crosslink collagen; "
                        "the chains are in the core matrisome set"),
    "class_I_light_chain": (lambda g: g == "B2M", "B2M is the MHC class I light chain; it is in the class I set"),
}

# ---------------------------------------------------------------- catalog
# group, set id, label, role, source, editing rules. Labels in US spelling.
# The catalog holds the programs that a script of this folder reads: the primary, contrast and technical-control
# programs, and three lists with another role. e2f_targets (role sensitivity) and matrix.tgfb (role matrix panel)
# are scored by 06_thioredoxin.py for Supplementary Data 1 (sheet Thioredoxin programs); activation.chu is promoted
# to primary by s1 because the first-choice activation list could not be verified. "msigdb:" sources are read from
# the shipped MSigDB responses, "lit:" sources from literature_lists.json.gz.
CAT = [
 ("Fibroblast state", "complement_cascade", "Complement cascade genes", "primary", "msigdb:REACTOME_COMPLEMENT_CASCADE", ["immunoglobulin", "duplicate_locus_symbol"]),
 ("Fibroblast state", "adventitial_fibroblast", "Adventitial fibroblast signature", "primary", "msigdb:TRAVAGLINI_LUNG_ADVENTITIAL_FIBROBLAST_CELL", []),
 ("Fibroblast state", "alveolar_fibroblast", "Alveolar fibroblast signature", "primary", "msigdb:TRAVAGLINI_LUNG_ALVEOLAR_FIBROBLAST_CELL", []),
 ("Fibroblast state", "mycaf.cords", "Myofibroblastic CAF signature", "primary", "lit:cords2023_mCAF", []),
 ("Fibroblast state", "icaf.cords", "Inflammatory CAF signature", "primary", "lit:cords2023_iCAF", []),
 ("Fibroblast state", "core_matrisome", "Core matrisome genes", "primary", "msigdb:NABA_CORE_MATRISOME", []),
 ("Fibroblast state", "collagen_modification", "Collagen modification and crosslinking genes", "primary",
  "msigdb:REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES+msigdb:REACTOME_CROSSLINKING_OF_COLLAGEN_FIBRILS", ["collagen_chains"]),
 ("Redox", "nrf2_targets", "NRF2 target genes", "primary", "msigdb:SINGH_NFE2L2_TARGETS", []),
 ("Antigen presentation", "mhc_class_ii", "MHC class II genes", "primary", "msigdb:GOCC_MHC_CLASS_II_PROTEIN_COMPLEX", ["class_I_light_chain"]),
 ("Antigen presentation", "mhc_class_i", "MHC class I genes", "contrast", "msigdb:GOCC_MHC_CLASS_I_PROTEIN_COMPLEX", []),
 ("Growth and biosynthesis", "cell_cycle", "Cell-cycle genes", "primary", "lit:tirosh2016_S_phase+lit:tirosh2016_G2M_phase", []),
 ("Growth and biosynthesis", "e2f_targets", "E2F target genes", "sensitivity", "msigdb:HALLMARK_E2F_TARGETS", []),  # read by 06
 ("Growth and biosynthesis", "myc_targets", "MYC target genes", "primary", "msigdb:HALLMARK_MYC_TARGETS_V1", []),
 ("Growth and biosynthesis", "mtorc1", "mTORC1 signaling genes", "primary", "msigdb:HALLMARK_MTORC1_SIGNALING", []),
 ("Interferon and inflammation", "ifn_gamma", "Interferon-gamma response genes", "primary", "msigdb:HALLMARK_INTERFERON_GAMMA_RESPONSE", []),
 ("Interferon and inflammation", "inflammatory_response", "Inflammatory response genes", "primary", "msigdb:HALLMARK_INFLAMMATORY_RESPONSE", []),
 ("Macrophage state", "lam", "Lipid-associated macrophage signature", "primary", "lit:jaitin2019_LAM", []),
 ("Macrophage state", "resident_like_macrophage", "LYVE1+ resident-like macrophage markers", "primary", "lit:cheng2021_Macro_LYVE1", []),
 ("Macrophage state", "classical_complement", "Classical complement pathway genes", "primary", "msigdb:GOBP_COMPLEMENT_ACTIVATION_CLASSICAL_PATHWAY", ["immunoglobulin"]),
 ("Macrophage state", "lysosome", "Lysosome genes", "primary", "msigdb:KEGG_LYSOSOME", []),
 ("Macrophage state", "apoptotic_cell_clearance", "Apoptotic cell clearance genes", "primary", "msigdb:GOBP_APOPTOTIC_CELL_CLEARANCE", []),
 ("Macrophage state", "inflammasome", "Inflammasome genes", "primary", "msigdb:REACTOME_INFLAMMASOMES", []),
 ("Macrophage state", "iron_homeostasis", "Intracellular iron homeostasis genes", "primary", "msigdb:GOBP_INTRACELLULAR_IRON_ION_HOMEOSTASIS", []),
 ("Macrophage state", "hypoxia", "Hypoxia response genes", "primary", "msigdb:HALLMARK_HYPOXIA", []),
 ("Macrophage state", "angiogenesis", "Angiogenesis genes", "primary", "msigdb:HALLMARK_ANGIOGENESIS", []),
 ("Macrophage state", "protein_secretion", "Protein secretion genes", "primary", "msigdb:HALLMARK_PROTEIN_SECRETION", []),
 ("T and NK cell state", "exhaustion", "T cell exhaustion markers", "primary", "lit:chu2023_exhaustion", []),
 ("T and NK cell state", "cytotoxic", "Cytotoxic effector genes", "primary", "lit:chu2023_cytotoxicity", []),
 ("T and NK cell state", "treg", "Tumor-infiltrating Treg signature", "primary", "lit:magnuson2018_tumour_treg", []),
 ("T and NK cell state", "naive_memory", "Memory-like CD8 T cell markers", "primary", "lit:sadefeldman2018_memory_like_CD8", []),
 ("T and NK cell state", "activation", "T cell activation signature", "primary", "lit:szabo2019_activation", []),
 ("T and NK cell state", "activation.chu", "T cell activation signature", "sensitivity", "lit:chu2023_activation", []),  # promoted to primary
 ("Technical control", "dissociation_stress", "Dissociation stress score", "technical control", "lit:denisenko2020_stress17", []),
 ("Matrix panel", "matrix.tgfb", "TGF-beta signaling (Hallmark)", "matrix panel", "msigdb:HALLMARK_TGF_BETA_SIGNALING", []),  # read by 06
]


def read_gz_json(name):
    """A gzipped JSON file of chapter2/data/sources."""
    with gzip.open(SRC / name, "rt", encoding="utf-8") as fh:
        return json.load(fh)


class Metadata:
    """PubMed and HGNC responses from the shipped cache; the network only with --allow-network."""

    def __init__(self, allow_network):
        self.cache = read_gz_json("pubmed_hgnc_responses.json.gz")
        self.allow, self.fetched = allow_network, 0

    def get(self, url, accept=None):
        """The response for a URL: from the cache, else from the network (only with --allow-network; 4 tries)."""
        if url in self.cache:
            return self.cache[url]
        if not self.allow:
            raise SystemExit(f"not in the shipped PubMed and HGNC responses: {url}\n"
                             "rerun with --allow-network to fetch it")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", **({"Accept": accept} if accept else {})})
        txt = None
        for k in range(4):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    txt = r.read().decode("utf-8", "replace")
                    break
            except Exception:
                time.sleep(2 + 3 * k)
        self.cache[url] = txt
        self.fetched += 1
        time.sleep(0.12)
        return txt

    def pubmed(self, pmid):
        """PubMed summary of one PMID: doi, year, journal, title and a short first-author citation."""
        t = self.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json")
        r = json.loads(t)["result"][str(pmid)]
        doi = next((a["value"] for a in r.get("articleids", []) if a["idtype"] == "doi"), "")
        a0 = r["authors"][0] if r.get("authors") else {"name": "", "authtype": ""}
        au = a0["name"] if a0.get("authtype") == "CollectiveName" else f"{a0['name'].rsplit(' ', 1)[0]} et al."
        return dict(pmid=str(pmid), doi=doi, year=r["pubdate"][:4], journal=r["source"], title=r["title"].rstrip("."),
                    short=f"{au} {r['pubdate'][:4]}, {r['source']}")

    def pmid_from_doi(self, doi):
        """The PMID that PubMed holds for a DOI, or "" when there is not exactly one."""
        t = self.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term="
                     + urllib.parse.quote(f"{doi}[doi]"))
        ids = json.loads(t)["esearchresult"]["idlist"] if t else []
        return ids[0] if len(ids) == 1 else ""


def gene_axis():
    """Gene symbols of the atlas object (var index), or of the fibroblast working set."""
    named = paths.CFG.get("atlas_object")
    if named and Path(named).expanduser().exists():
        with h5py.File(Path(named).expanduser(), "r") as f:
            gk = f["var"].attrs["_index"]
            gk = gk.decode() if isinstance(gk, bytes) else gk
            v = f[f"var/{gk}"][:]
        return {x.decode() if isinstance(x, bytes) else str(x) for x in v}, "atlas object var index"
    g = L.work_dir("fib") / "genes.npy"
    if g.exists():
        return {str(x) for x in np.load(g, allow_pickle=True)}, "genes.npy of the fibroblast working set"
    raise SystemExit('no gene axis: set "atlas_object" in the configuration or run 02_working_sets.py first')


def s1_programs(a):
    """Build the gmt, the definitions table and the nested-program check from the catalog."""
    RES = read_gz_json("msigdb_v2026.1.Hs_resource_sets.json.gz")
    LIT = read_gz_json("literature_lists.json.gz")
    LITD = {e["id"]: e for grp in LIT["groups"] for e in grp["lists"]}
    OBJ, axis = gene_axis()
    meta = Metadata(a.allow_network)
    MSIG_VER = [v for v in RES["versions"]["msigdb"] if v.endswith(".Hs")][0]
    REACT_VER = f"Reactome release {RES['versions']['reactome']}"
    LOG(f"01_programs.py s1: {len(RES['sets'])} MSigDB sets ({MSIG_VER}), {len(LITD)} paper lists, "
        f"gene axis from the {axis} ({len(OBJ):,} symbols)")

    def resource_part(name):
        """Genes and source description of one MSigDB set, with its resource paper and release."""
        s = RES["sets"][name]
        src = s.get("exactSource") or ""
        if name.startswith("REACTOME_"):
            kind, ident = "curated pathway (Reactome)", f"Reactome {src}; MSigDB {s['systematicName']}"
            paper, ver = RESOURCE_PAPER["Reactome"], f"MSigDB {MSIG_VER}; {REACT_VER}"
            xc = RES.get("reactome_crosscheck", {}).get(name)
            if xc:
                def clean(gs):
                    """Drop immunoglobulin, non-human and duplicate-locus symbols before the comparison."""
                    return {g for g in gs if not IG.match(g) and not NONHUMAN.match(g) and g != "C4B_2"}
                x, y = clean(s["geneSymbols"]), clean(xc["reactome_genes"])
                check = ("MSigDB list equals the human Reactome participants" if x == y else
                         f"MSigDB list differs from Reactome participants: only MSigDB {sorted(x - y)}, only Reactome {sorted(y - x)}")
                check += " (immunoglobulin and non-human participants set aside for the comparison)"
            else:
                check = "fetched from the MSigDB API"
        elif name.startswith("GOBP_") or name.startswith("GOCC_"):
            kind, ident = "ontology term (Gene Ontology)", f"{src}; MSigDB {s['systematicName']}"
            paper, ver, check = RESOURCE_PAPER["GO"], f"MSigDB {MSIG_VER}", "fetched from the MSigDB API"
        elif name.startswith("HALLMARK_"):
            kind, ident = "curated gene set (MSigDB Hallmark)", f"MSigDB {s['systematicName']}"
            paper, ver, check = RESOURCE_PAPER["Hallmark"], f"MSigDB {MSIG_VER}", "fetched from the MSigDB API"
        elif name.startswith("KEGG_"):
            kind, ident = "curated pathway (KEGG, MSigDB legacy)", f"KEGG {src}; MSigDB {s['systematicName']}"
            paper, ver, check = RESOURCE_PAPER["KEGG"], f"MSigDB {MSIG_VER}", "fetched from the MSigDB API"
        elif name.startswith("NABA_"):
            kind, ident = "published compendium (matrisome)", f"MSigDB {s['systematicName']}"
            paper, ver, check = RESOURCE_PAPER["NABA"], f"MSigDB {MSIG_VER}", "fetched from the MSigDB API"
        elif name.startswith("TRAVAGLINI_"):
            kind, ident = "published signature (human lung cell atlas)", f"MSigDB C8 {s['systematicName']}"
            paper, ver, check = RESOURCE_PAPER["Travaglini"], f"MSigDB {MSIG_VER}", "fetched from the MSigDB API"
        else:
            kind, ident = "published signature (MSigDB curated)", f"MSigDB {s['systematicName']}"
            paper, ver, check = str(s.get("pmid") or ""), f"MSigDB {MSIG_VER}", "fetched from the MSigDB API"
        where = f"MSigDB gene set {name}" + (f"; source given by MSigDB: {src}"
                                             if src and not src.startswith(("R-HSA", "GO:", "hsa")) else "")
        return dict(genes=list(s["geneSymbols"]), kind=kind, identifier=ident, paper=paper, version=ver,
                    where=where, species="human", verification=check, url=s.get("msigdbURL", ""))

    def lit_part(lid):
        """Genes and source description of one paper-derived list, or (None, reason) when it could not be verified."""
        e = LITD.get(lid)
        if e is None or e["status"] == "not_found_by_either" or not e.get("final") or not e["final"].get("found"):
            why = "not requested" if e is None else (e["status"] if e["status"] == "not_found_by_either"
                                                     else "adjudicator could not confirm the list")
            return None, why
        f = e["final"]
        genes = [g.strip() for g in (f.get("human_symbols") or []) if g and g.strip()]
        pm = (f.get("pmid") or "").strip() or (meta.pmid_from_doi(f["doi"]) if f.get("doi") else "")
        ver = {"agreed": "two independent readers extracted the identical list from the source",
               "adjudicated": "two readers disagreed; a third reader settled it from the source"}[e["status"]]
        if e["status"] == "adjudicated":
            ver += (f" (only reader A: {', '.join(e.get('only_in_A', [])) or 'none'}; "
                    f"only reader B: {', '.join(e.get('only_in_B', [])) or 'none'})")
        return dict(genes=genes, kind="published signature", identifier=f.get("pmcid") or "", paper=pm,
                    version="as printed in the source", where=f.get("where_in_source", ""), species=f.get("species", ""),
                    verification=ver, url=f.get("access_route", ""), selection=f.get("selection_rule", ""),
                    notes=f.get("notes", "")), ""

    def hgnc_to_object(sym):
        """The atlas symbol for a gene symbol: itself when present, else the one HGNC previous or alias match."""
        if sym in OBJ:
            return sym
        hits = set()
        for field in ("symbol", "prev_symbol", "alias_symbol"):
            t = meta.get(f"https://rest.genenames.org/fetch/{field}/{urllib.parse.quote(sym)}", "application/json")
            if not t:
                continue
            for d in json.loads(t)["response"]["docs"]:
                for s in [d.get("symbol")] + d.get("prev_symbol", []) + d.get("alias_symbol", []):
                    if s in OBJ:
                        hits.add(s)
        return next(iter(hits)) if len(hits) == 1 else sym

    rows, gmt, missing = [], [], []
    for group, sid, name, role, source, rules in CAT:
        parts = []
        for piece in source.split("+"):
            kind, key = piece.split(":", 1)
            if kind == "msigdb":
                parts.append(resource_part(key))
            else:
                p, why = lit_part(key)
                if p is None:
                    missing.append((sid, name, role, key, why))
                    parts = None
                    break
                parts.append(p)
        if not parts:
            continue
        listed = list(dict.fromkeys(g for p in parts for g in p["genes"]))
        removed, used = {}, []
        for g in listed:
            hit = next((r for r in rules if RULES[r][0](g)), None)
            if hit:
                removed.setdefault(hit, []).append(g)
            else:
                used.append(g)
        renamed, final = [], []
        for g in used:
            s2 = hgnc_to_object(g)
            if s2 != g:
                renamed.append(f"{g} -> {s2}")
            final.append(s2)
        final = list(dict.fromkeys(final))
        papers = [meta.pubmed(p["paper"]) for p in parts if p["paper"]]
        ref = "; ".join(dict.fromkeys(f"{m['short']}. {m['title']}" for m in papers))
        rows.append(dict(
            set_id=sid, program=name, role=role, group=group,
            source_type=" + ".join(dict.fromkeys(p["kind"] for p in parts)),
            identifier=" + ".join(dict.fromkeys(p["identifier"] for p in parts if p["identifier"])),
            reference=ref, doi="; ".join(dict.fromkeys(m["doi"] for m in papers if m["doi"])),
            pmid="; ".join(dict.fromkeys(m["pmid"] for m in papers)),
            where_in_source=" + ".join(p["where"] for p in parts),
            selection_in_source=" + ".join(p.get("selection", "") for p in parts if p.get("selection")),
            resource_version=" + ".join(dict.fromkeys(p["version"] for p in parts)),
            species_of_source=" + ".join(dict.fromkeys(p["species"] for p in parts)),
            n_genes_in_source=len(listed),
            genes_removed="; ".join(f"{RULES[r][1]}: {' '.join(gs)}" for r, gs in removed.items()),
            n_genes_used=len(final), genes_used=" ".join(final),
            symbols_renamed_to_object="; ".join(renamed),
            verification=" + ".join(p["verification"] for p in parts),
            reader_notes=" | ".join(p.get("notes", "") for p in parts if p.get("notes")),
            source_url=" ; ".join(p["url"] for p in parts if p["url"])))
        gmt.append(f"{sid}\t{name} | {role} | {ref.split('. ')[0] if ref else ''}\t" + "\t".join(final))

    D = pd.DataFrame(rows)
    # a program whose primary list could not be verified takes its first verified sensitivity list (one case:
    # activation.chu for activation)
    for name in {m[1] for m in missing if m[2] == "primary"}:
        sub = D[D.program == name]
        if not (sub.role == "primary").any() and (sub.role == "sensitivity").any():
            i = sub.index[0]
            D.loc[i, "role"] = "primary"
            D.loc[i, "verification"] += " + promoted to primary because the planned primary list could not be verified"

    def nodash(v):
        """Replace long dashes in a text value by plain ones."""
        return v.replace("\u2014", " - ").replace("\u2013", "-") if isinstance(v, str) else v
    D = D.apply(lambda c: c.map(nodash))
    gmt = [nodash(line) for line in gmt]
    o = L.out("gene_programs_definitions.csv")
    D.to_csv(o, index=False)
    L.sidecar(o, "One row per gene program (set) as built by 01_programs.py s1; the columns of "
              "published definitions table with the label column named program.", [
                  ("set_id", "program id, used by every Chapter 2 table"),
                  ("program", "label (US spelling); the printed label may differ, see supptable_gene_programs.csv"),
                  ("role", "primary (used in text and figures), contrast, technical control; sensitivity (e2f_targets) "
                           "and matrix panel (matrix.tgfb): kept only because 06_thioredoxin.py scores them, not drawn"),
                  ("group", "group of Supplementary Figure 3"), ("source_type", "kind of source"),
                  ("identifier", "database identifier or PMC id"), ("reference", "first author, year, journal and title (PubMed)"),
                  ("doi", "DOI of the reference"), ("pmid", "PubMed id of the reference"),
                  ("where_in_source", "where the list is in the source"), ("selection_in_source", "the source's own selection rule"),
                  ("resource_version", "database release, or 'as printed in the source'"),
                  ("species_of_source", "species of the source list"), ("n_genes_in_source", "genes in the source list"),
                  ("genes_removed", "rule-based edits with the reason and the genes removed"),
                  ("n_genes_used", "genes in the list used"), ("genes_used", "the genes used, space separated"),
                  ("symbols_renamed_to_object", "symbols replaced by the atlas symbol (HGNC, one match)"),
                  ("verification", "how the list was checked"), ("reader_notes", "the readers' notes, verbatim"),
                  ("source_url", "access route")])
    (L.out("gene_programs.gmt")).write_text("\n".join(gmt) + "\n")
    L.sidecar(L.out("gene_programs.gmt"), "The gene programs, one per line, tab separated (gmt format); TXNRD1 "
              "is removed from every list at scoring, not here.",
              [("field 1", "set id"), ("field 2", "label | role | first reference"), ("fields 3 on", "the genes")])
    # the list-in-list check behind the nested programs of Figure 3
    GU = dict(zip(D.set_id, D.genes_used.str.split()))
    nest = []
    for child, spec in L.panels()["nested"].items():
        c, p = GU[child], set(GU[spec["parent"]])
        nest.append(dict(child=child, parent=spec["parent"], n_child=len(c), n_child_in_parent=sum(g in p for g in c),
                         child_genes=" ".join(c)))
    o = L.out("gene_programs_nested.csv")
    pd.DataFrame(nest).to_csv(o, index=False)
    L.sidecar(o, "Programs drawn nested under another in Figure 3: how many of the child's genes are in the parent.",
              [("child", "nested program id"), ("parent", "program it is drawn under"), ("n_child", "genes of the child"),
               ("n_child_in_parent", "of them, in the parent's list"), ("child_genes", "the child's genes")])
    LOG(f"  {len(D)} programs written, {meta.fetched} network requests")
    for r in nest:
        LOG(f"  nested: {r['child']} has {r['n_child']} genes, {r['n_child_in_parent']} of them in {r['parent']}")


def tested_universe(cell):
    """Genes that pass arm presence 0.60 and detection 0.02 in the capped working set of a cell set, and its cell count."""
    S = L.load_set(cell)
    return S["UNI"], S["N"]


def s2_table(a):
    """Write supptable_gene_programs.csv: every program with its source, genes and tested genes per cell set."""
    D = pd.read_csv(L.out("gene_programs_definitions.csv"), dtype=str).fillna("")
    cfg = L.panels()
    T = pd.DataFrame({
        "program_id": D.set_id,
        "label_as_printed": [cfg["printed_labels"].get(s, L.us(p)) for s, p in zip(D.set_id, D.program)],
        "role": D.role, "group": D.group, "source_type": D.source_type, "source_identifier": D.identifier,
        "reference": D.reference, "pmid": D.pmid, "doi": D.doi, "database_version": D.resource_version,
        "where_in_source": D.where_in_source, "n_genes_in_source": D.n_genes_in_source.astype(int),
        "edits": D.genes_removed, "symbols_renamed": D.symbols_renamed_to_object,
        "n_genes_used": D.n_genes_used.astype(int), "genes": D.genes_used})
    for cell, spec in L.CELL_SETS.items():
        uni, n = tested_universe(cell)
        tested = [[g for g in gs.split() if g in uni and g != L.SPLITTER] for gs in D.genes_used]
        T[f"n_tested_{cell}"] = [len(t) for t in tested]
        T[f"scored_{cell}"] = [len(t) >= L.MINSET for t in tested]
        T[f"not_testable_{cell}"] = [" ".join(g for g in gs.split() if g not in uni and g != L.SPLITTER)
                                     for gs in D.genes_used]
        LOG(f"  {spec['label']}: {n:,} cells after the cap, {len(uni):,} genes in the tested universe, "
            f"{int(T[f'scored_{cell}'].sum())} programs with at least {L.MINSET} tested genes")
    o = L.out("supptable_gene_programs.csv")
    T.to_csv(o, index=False)
    cols = [("program_id", "program id, as in every Chapter 2 table"),
            ("label_as_printed", "label as printed in Figure 3 and Supplementary Figure 3"),
            ("role", "primary: used in the text and figures; contrast; technical control; sensitivity (e2f_targets, "
                     "another published definition of proliferation) and matrix panel (matrix.tgfb, of an earlier "
                     "matrix panel): kept only because 06_thioredoxin.py scores them for Supplementary "
                     "Data 1 (Thioredoxin programs), not drawn"),
            ("group", "group of Supplementary Figure 3"), ("source_type", "kind of source"),
            ("source_identifier", "database identifier (Reactome, GO, KEGG, MSigDB systematic name) or PMC id"),
            ("reference", "first author, year, journal and title"), ("pmid", "PubMed id"), ("doi", "DOI"),
            ("database_version", "MSigDB release v2026.1.Hs, Reactome release 97, or 'as printed in the source' "
                                 "for paper lists"),
            ("where_in_source", "where the list is in the source"), ("n_genes_in_source", "genes in the source list"),
            ("edits", "rule-based edits with the reason and the genes removed (none when empty)"),
            ("symbols_renamed", "symbols replaced by the atlas symbol through HGNC"),
            ("n_genes_used", "genes in the list used (TXNRD1 included when listed; it is removed at scoring)"),
            ("genes", "every gene of the list used, space separated")]
    for cell, spec in L.CELL_SETS.items():
        cols += [(f"n_tested_{cell}", f"genes tested in {spec['short']}: detected in at least 2 % of the cells whose "
                                      "arm measures the gene and measured in at least 60 % of the cells, TXNRD1 removed"),
                 (f"scored_{cell}", f"True when at least {L.MINSET} genes are tested, so the program is scored in {spec['short']}"),
                 (f"not_testable_{cell}", f"genes of the list not tested in {spec['short']}")]
    L.sidecar(o, "Every gene program of Chapter 2 with its source, genes and tested genes per cell set; "
              "08_supptable2.py s2 builds Supplementary Data 1 (Gene programs; the table the Supplementary Figure 3 footnote "
              "cites) from its panel programs. Redistribution: lists from MSigDB KEGG legacy sets (lysosome) carry KEGG terms; "
              "see README.", cols)
    LOG(f"wrote {o.name}: {len(T)} programs")


def main():
    """Run s1 and/or s2 and write the log."""
    ap = argparse.ArgumentParser(description=__doc__.split("RUN:")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", default=["s1", "s2"], choices=["s1", "s2"])
    ap.add_argument("--allow-network", action="store_true", help="fetch PubMed or HGNC entries missing from the cache")
    a = ap.parse_args()
    if "s1" in a.only:
        s1_programs(a)
    if "s2" in a.only:
        s2_table(a)
    LOG.write(paths.WORK_ROOT / "logs" / "chapter2" / "01_programs.log")


if __name__ == "__main__":
    (paths.WORK_ROOT / "logs" / "chapter2").mkdir(parents=True, exist_ok=True)
    main()
