# stromal-trxr1-nsclc

## 1. What this is

This repository holds the analysis code of the article "Stromal thioredoxin reductase 1, T-cell-rich stroma and
survival in lung squamous cell carcinoma". The code computes every number, table and figure panel of the article
from three sources: the image-analysis exports of multiplex immunofluorescence on tissue microarrays of two Swedish
NSCLC cohorts, the clinical tables of those cohorts, and public single-cell RNA data that the code assembles into an
atlas of fibroblasts, myeloid cells, T and NK cells and cancer cells.

The code calls the two cohorts BOMI1 and BOMI2. The article calls them **Uppsala I (BOMI1)** and **Uppsala II
(BOMI2)**. The code calls niche connectivity "confluence" (the column `biggest_30`).

Layout:

| folder or file | what it holds |
|---|---|
| `mif/` | three scripts that build the per-cohort imaging tables from the raw exports and clinical files |
| `redox/` | the shared Python module (paths, cell tables, clinical tables, survival models, statistics, plotting) |
| `chapter1.py`, `tumor_vs_stroma.py`, `chapter3.py`, `chapter4.py` | the imaging analyses |
| `atlas/` | the build of the single-cell atlas from public deposits (stages 00 to 06; `atlas/steps/` holds the step scripts) |
| `chapter2/` | the single-cell analyses (twelve numbered scripts, the library `ch2lib.py` and the gene-list sources in `chapter2/data/`) |
| `manifest/` | `manifest.csv`, every value the article prints with the table that holds it, and `check.py`, which compares them |
| `environment/` | the Python and R package versions |

## 2. What runs without patient data

The imaging data and clinical tables of Uppsala I and II are not public (section 9). Without them, the following
parts run from public data alone:

- `atlas/`: downloads the public deposits listed in `atlas/sources.tsv`, checks each by sha256, and builds the atlas
  objects. `atlas/05_obs_summary.py` compares the rebuilt object with `atlas/reference/fullgene_obs_summary.tsv`.
- `chapter2/`: every single-cell result of the article (Figures 2 and 3, Supplementary Figures 2 to 4 and
  Supplementary Data 1) from the atlas objects.
- The manifest rows of the single-cell part (file names under `chapter2/` and `supplementary/`).

The imaging scripts (`mif/`, `chapter1.py`, `tumor_vs_stroma.py`, `chapter3.py`, `chapter4.py`) need the restricted
data. They write the patient-level tables under
`work_root/intermediate/`, which stay with the data and are never shared.

## 3. Requirements and installation

- Python 3.13 for the analyses (the outputs of the article were made with Python 3.13 on an Apple M2). The atlas build
  ran under Python 3.12 with anndata and scanpy; the versions are at the end of `environment/requirements.txt`.
- R 4.5 for `mif/02_cells.R` and three atlas steps.

```bash
pip install -r environment/requirements.txt
Rscript environment/r_packages.R
```

For the atlas build, also:

```bash
pip install "anndata==0.12.*" "scanpy==1.12.*" pyarrow
```

`environment/r_packages.R` installs data.table, readxl, autothresholdr, jsonlite and Matrix and warns when an
installed version differs from the one used for the article.

## 4. Configuration

Copy `config.example.json` to `config.json` and set the paths. Every script reads `config.json` in the repository
root, or the file named by the environment variable `REDOX_CONFIG`.

| key | needed by | what it names |
|---|---|---|
| `raw_root` | all | the folder of the raw data: `Raw_original_BOMI2/` and `Raw_original_BOMI1/` (the inForm cell and tissue exports), `Bomi2/` and `Bomi1/` (clinical files, core keys, mutation summary) and `external_scRNA/` (the public single-cell deposits, as `atlas/00_download.sh` writes them) |
| `work_root` | all | the folder every script writes to: `intermediate/` (patient-level tables), `output/` (the tables and panels of the article) and `atlas/` |
| `input_files` | `mif/` | the file names of three raw inputs: the Uppsala II clinical database sheet in `raw_root/Bomi2/` and the Uppsala I merged cell and tissue exports in `raw_root/Raw_original_BOMI1/` |
| `examples` | `chapter3.py`, `chapter4.py` | `examples.json`, the example cores and patients drawn in Figures 4c, 5a, 5c, 5d, 5f and Supplementary Figure 6; `examples.example.json` shows the format, and its keys use the working names of the panels (`fig3c` is Figure 4c; `fig4a_window` is the window of Figures 5a and 5c; `fig4d` and `fig4f_squares` are Figures 5d and 5f). Without it, the example panels are skipped with a message |
| `pellinen` | `chapter1.py` | optional: the folder with the published fibroblast data of Pellinen and colleagues, used to reconstruct alpha-SMA-positive fibroblasts in Uppsala I, whose panel has no alpha-SMA (the reconstructed rows of Figure 1e and Supplementary Table 1c) |
| `gse81089` | `tumor_vs_stroma.py` | optional: the bulk RNA table of the Uppsala II tumors (GEO GSE81089; one row per patient, one column per Ensembl gene). Without it, the bulk RNA rows and Supplementary Figure 1f are skipped |
| `atlas_object` | `chapter2/` | the atlas object `fullgene_object.h5ad` that `atlas/04_build_objects.py` writes |
| `epi_object` | `chapter2/02_working_sets.py --cells can` | the epithelial-cell object `epi_fullgene_object.h5ad` that `atlas/06_cancer_cell_object.py` writes |
| `scrna_raw`, `atlas_root` | `atlas/` | optional: other places for the public deposits (default `raw_root/external_scRNA`) and for the atlas build (default `work_root/atlas`) |
| `rscript` | `atlas/` | optional: the R interpreter (default: `Rscript` on the PATH) |

Optional keys may be left empty. `config.json` and `examples.json` hold local paths and identifiers and are not part
of the repository.

## 5. Run order

The imaging part and the single-cell part are independent. All commands run from the repository root. Runtimes are
approximate: the lower values were measured on an Apple M2 (8 cores, 24 GB of memory), the upper ones on an
Intel machine.

**Imaging part** (restricted data):

| command | needs | writes | runtime |
|---|---|---|---|
| `python3 mif/01_clinical.py` | raw clinical files | `intermediate/{BOMI2,BOMI1}/clinical.csv` (de-identified) and `clinical_full.csv` | seconds |
| `Rscript mif/02_cells.R --cohort BOMI2` | raw cell and tissue exports | `intermediate/BOMI2/thresholds.csv`, `cells.tsv`, `areas.csv` | 0.5 to 5 min |
| `Rscript mif/02_cells.R --cohort BOMI1` | the same for Uppsala I | `intermediate/BOMI1/...` | 0.5 to 2.5 min |
| `python3 mif/03_percase.py` | the tables of 01 and 02 | `intermediate/{BOMI2,BOMI1}/percase.csv` | 20 s |
| `python3 chapter1.py` | the mif tables, optional `pellinen` | `output/chapter1/` | 1 to 2.5 min |
| `python3 chapter3.py` | the mif tables, `examples` | `output/chapter3/` | 1.5 to 5 min |
| `python3 chapter4.py` | the mif tables, `examples` | `output/chapter4/` | 2 to 8 min |
| `python3 tumor_vs_stroma.py` | the mif tables, raw cell exports, mutation and smoking files, optional `gse81089` | `output/tumor_vs_stroma/` | 1 to 3.5 min |

`mif/02_cells.R --phase-a` (Uppsala II only) also repeats the review of the hand-annotated cores
(`Raw_original_BOMI2/outputPheno/`) from which the three fixed Uppsala II lineage cuts were derived, and writes its
suggestions to `thresholds_phase_a.csv`; it changes no other output. Every analysis script takes `--only` with the
names of its sections (s1, s2, ...); the sections and their outputs are listed in the header of each script.

**Single-cell part** (public data):

| command | needs | writes | runtime |
|---|---|---|---|
| `bash atlas/00_download.sh` | network | the deposits of `atlas/sources.tsv` under `raw_root/external_scRNA`, each checked by sha256 (`--check-only`, `--urls-only`) | about 27 GB of downloads |
| `python atlas/01_arms_luca_lambrechts.py` | stage 00 | gene panel; LuCA, E-MTAB-6149 and E-MTAB-6653 arms; three-arm panel object; block QC; E-MTAB-6653 label transfer | about 80 min |
| `python atlas/02_arm_emtab13526.py` | stages 00, 01; network (reads the object of the original authors by HTTP range requests) | the E-MTAB-13526 arm and its merge | about 1 h |
| `python atlas/03_arms_wave1.py` | stages 00 to 02 | GSE325414, GSE246215 and STTT arms | about 85 min |
| `python atlas/04_build_objects.py` | stages 00 to 03 | full-transcriptome arms, `fullgene_object.h5ad` (set `atlas_object`), `panel_object_v2.h5ad`, block QC on 24 blocks, QC transfer, obs repair, state check | several hours |
| `python atlas/05_obs_summary.py` | stage 04 | the obs summary of the object, compared with `atlas/reference/fullgene_obs_summary.tsv` | 10 s |
| `python atlas/06_cancer_cell_object.py` | stages 00 to 04 | `epi_fullgene_object.h5ad`, the epithelial cells of primary LUSC and LUAD tumors (set `epi_object`) | several hours |
| `python3 chapter2/02_working_sets.py` | `atlas_object` | the fibroblast, macrophage and T/NK working sets | 1 min |
| `python3 chapter2/01_programs.py` | `chapter2/data/` | the gene programs (`gene_programs.gmt`) and their definitions | seconds |
| `python3 chapter2/03_program_tests.py` | 01, 02 | the program tests of Figure 3 | about 10 min with 6 workers |
| `python3 chapter2/04_fibroblast_genes.py` | 02 | the per-gene fibroblast scan of Supplementary Figure 4 | 1.5 min |
| `python3 chapter2/05_plots.py` | 03, 04 | Figure 3, Supplementary Figures 3 and 4 | seconds |
| `python3 chapter2/06_thioredoxin.py` | 01, 02 | Figure 2c (stromal cell sets) and the sheet Thioredoxin programs | 3.5 min |
| `python3 chapter2/07_smoking_atlas.py` | 02, `atlas_object` | Supplementary Figure 2b | 15 s |
| `python3 chapter2/08_supptable2.py` | 01, 03, 06 | the tables of Supplementary Data 1 | under a second |
| `python3 chapter2/09_txnrd1_profile.py` | `atlas_object`, 02 | Figure 2a | 45 s |
| `python3 chapter2/10_supplementary_data1.py` | 08 | Supplementary Data 1 (Excel) | under a second |

Atlas stages take `--list`, `--dry-run`, `--only step ...` and `--from step`, run one process per step and write
one log per step under `work_root/atlas/logs/`. A clean rebuild needs about 80 GB of free disk; the steps stream the
large matrices by row block or library and ran with 24 GB of memory.

Then the leave-out runs. They repeat the single-cell tests without the cells of the pulsed electric field ablation
study (Jimenez_Bruno_2026, GSE325414), as the Methods state, and for T and NK cells without the cells of Guo and
colleagues (Guo_Zhang_2018, PMID 29942094), as the Results state (about 4.5 minutes in all):

```bash
python3 chapter2/03_program_tests.py --only s1 --leave-out-study Jimenez_Bruno_2026
python3 chapter2/04_fibroblast_genes.py --leave-out-study Jimenez_Bruno_2026
python3 chapter2/06_thioredoxin.py --leave-out-study Jimenez_Bruno_2026
python3 chapter2/09_txnrd1_profile.py --only s2 s3 --leave-out-study Jimenez_Bruno_2026
python3 chapter2/03_program_tests.py --only s1 --cells tnk --leave-out-study Guo_Zhang_2018
```

A leave-out run writes its tables with the suffix `_without_<study>` beside the ordinary ones, which it leaves
untouched. The study is given as the prefix of its batch names and is removed before the cap of 2,000 cells per batch
and histology. Because the cap and the pairing draw from one generator in batch order, a leave-out run differs from
the ordinary one by more than the removed cells when the removed study was capped.

Then the cancer cells and the additional tests of the thioredoxin section (about 10 minutes). The cancer cells are a
reference: the cells that the original authors labeled as cancer cells, from the object of stage 06. The atlas object
and the tables of the three stromal cell sets are not changed by these runs.

```bash
python3 chapter2/02_working_sets.py --cells can
python3 chapter2/06_thioredoxin.py --only s1 --cells can
python3 chapter2/06_thioredoxin.py --only s1 --cells fib mac tnk can --split-cycle
python3 chapter2/06_thioredoxin.py --only s1 --cells can --histology LUSC
python3 chapter2/06_thioredoxin.py --only s1 --cells can --histology LUAD
python3 chapter2/06_thioredoxin.py --only s2 s3 s4
python3 chapter2/11_txn_txnip.py
python3 chapter2/08_supptable2.py --only s3
```

`--split-cycle` adds TXN, PRDX1, PRDX2 and SRXN1 as single-gene sets (Supplementary Figure 2a); `--histology` keeps
the cells of one histology within the sample of the main run (Supplementary Figure 2c). `11_txn_txnip.py` splits each
cell set by TXN and ranks the TXNIP difference among the differences given by genes expressed as strongly as TXN
(Supplementary Figure 2d).

Last, the T and NK subtypes and the workbook (about 2 minutes):

```bash
python3 chapter2/12_tnk_subtypes.py
python3 chapter2/10_supplementary_data1.py
```

`12_tnk_subtypes.py` tests the four programs drawn for T and NK cells in Figure 3 within conventional CD4 T cells,
regulatory T cells, CD8 T cells and NK cells, each with all cells and without the cells of Guo and colleagues. Its
first block must equal the T and NK rows of Figure 3, or the run stops.

### The single-cell tests

Each cell that detects TXNRD1 is paired 1:1 with a non-detecting cell of the same batch and histology with the
nearest number of genes detected (ties at random). The effect is the mean paired difference in a program score (the
mean within-batch z of the count fractions of its tested genes, TXNRD1 removed from every list). z and p come from
1,000 permutations of TXNRD1 detection within patient, with the matcher rerun in each. Verdict: no effect when
|z| < 1.96; otherwise withdrawn when the library-size split (control A) moves the score the same way with
|z| >= 1.96 and at least as far, not specific when it does so less far, and specific to TXNRD1 otherwise. These are
the verdict names of the result tables. There, "specific to TXNRD1" means only that the library-size split did not
reproduce the difference. The article uses "specific to TXNRD1" for the comparison with other genes, which is
described next. Figure 3 draws the programs of each cell set's panel (`chapter2/data/figure_panels.json`) whose
verdict is "specific to TXNRD1", that is, the programs with p < 0.05 that the library-size split did not reproduce.
In the thioredoxin section (Figure 2c, Supplementary Figure 2), the TXNRD1 effect is placed among the effects of splits by
the genes detected within 20 % of TXNRD1's rate; an entry is drawn filled at or above the 97.5th percentile with
|z| >= 1.96 and when it does not rest on one gene.

Random numbers: every script that draws them names its seeds in its header. The single-cell tests use `SEED = 6`
(the cap of 2,000 cells per batch and histology, the observed pairing and the permutation stream).
`chapter2/03_program_tests.py` draws all permuted pairings in order from one generator and computes the statistics
in `--workers` processes in the original order, so the results do not depend on the number of workers. The fast
matcher `match_fast` returns the pairs of the reference matcher `match_random_ties` (`chapter2/ch2lib.py`); every run
checks this at its start, and `--slow-matcher` uses the reference matcher throughout.

## 6. Outputs and the items of the article

Every output lies under `work_root/output/`, each table with a sidecar `.md` that describes its columns. The output
file names follow the working order of the analysis, not the figure numbers of the article (for example,
`chapter3/fig3a_clustering` is Figure 4a, and `chapter4/supptable4c.csv` is Supplementary Table 3b). So do column
names, the keys of `examples.json` and the panel labels `Fig3a` and `Fig3b` of `chapter3/fig3_stats.csv` (Figures 4a
and 4b), which other code reads. The sidecar texts and messages use the numbers of the article. The table below is
the reference. Panels are written as PNG and PDF.

| item of the article | output files | script |
|---|---|---|
| Figure 1a (schematic) | the counts in `chapter1/ch1_numbers.csv` and `tumor_vs_stroma/suppfig1_numbers.csv` | `chapter1.py`, `tumor_vs_stroma.py` |
| Figure 1b | `chapter1/fig1b_models.csv`, `fig1b_forest` | `chapter1.py` |
| Figure 1c | `chapter1/supptable1b.csv`, `fig1c_tiles` | `chapter1.py` |
| Figure 1d | `tumor_vs_stroma/lineage_profile_summary.csv`, `fig1d_lineage_profile` | `tumor_vs_stroma.py` |
| Figure 1e | `chapter1/supptable1c.csv`, `fig1e_pies` | `chapter1.py` |
| Figure 2a | `chapter2/txnrd1_profile_summary.csv`, `suppfig2a_txnrd1_profile` | `chapter2/09_txnrd1_profile.py` |
| Figure 2b | schematic, not code | |
| Figure 2c | `chapter2/trx_results.csv`, `trx_results_with_extra.csv`, `trx_numbers.csv`, `suppfig2c_programs_calibrated` | `chapter2/06_thioredoxin.py` |
| Figure 3 | `chapter2/program_results.csv`, `fig2_rows.csv`, `fig2_programs` | `chapter2/03_program_tests.py`, `05_plots.py` |
| Figure 4a | `chapter3/fig3a_moran_percore.csv`, `fig3_stats.csv`, `fig3a_clustering` | `chapter3.py` |
| Figure 4b | `chapter3/fig3b_neighborhood_percase.csv`, `fig3b_fibroblast_center` | `chapter3.py` |
| Figure 4c | `chapter3/fig3c_core_cells.csv`, `fig3c_map`, `fig3c_zoom` | `chapter3.py` |
| Figure 5a to 5h | `chapter4/fig4_stats.csv`, `fig4ef_patients.csv`, `fig4d_examples.csv`, `ch4_niche_survival.csv`; panels `fig4a_niches`, `fig4b_bars`, `fig4c_confluence`, `fig4d_cores`, `fig4e_histology`, `fig4f_fraction`, `fig4g_forest`, `fig4h_sweep` | `chapter4.py` |
| Supplementary Figure 1a to 1f | `tumor_vs_stroma/suppfig1a_density` to `suppfig1f_bulk`, `suppfig1_numbers.csv` | `tumor_vs_stroma.py` |
| Supplementary Figure 1g | micrograph, not code | |
| Supplementary Figure 2a | `chapter2/trx_cycle_genes.csv`, `panel_trx_cycle_genes` | `chapter2/06_thioredoxin.py` |
| Supplementary Figure 2b | `chapter2/smoking_atlas_tests.csv`, `suppfig2d_smoking_atlas` | `chapter2/07_smoking_atlas.py` |
| Supplementary Figure 2c | `chapter2/trx_cancer_histology.csv`, `panel_trx_cancer_histology` | `chapter2/06_thioredoxin.py` |
| Supplementary Figure 2d | `chapter2/txn_txnip.csv`, `panel_txn_txnip` | `chapter2/11_txn_txnip.py` |
| Supplementary Figure 3 | `chapter2/suppfig3_programs` (from the tables of 03) | `chapter2/05_plots.py` |
| Supplementary Figure 4 | `chapter2/suppfig4_rows.csv`, `suppfig4_dots.csv`, `suppfig4_fibroblast_genes` | `chapter2/04_fibroblast_genes.py`, `05_plots.py` |
| Supplementary Figure 5a, 5b, 5c | `chapter3/suppfig5a_same_label`, `suppfig5b_same_type_distance`, `suppfig5c_bomi1_unstained_center` and their tables | `chapter3.py` |
| Supplementary Figure 6 | `chapter4/suppfig6_stats.csv`, `suppfig6_matrix` | `chapter4.py` |
| Supplementary Table 1a | `tumor_vs_stroma/supptable1a.csv`, `supptable1a_celltypes.csv`, `supptable1a_covariates.csv` | `tumor_vs_stroma.py` |
| Supplementary Table 1b | `chapter1/supptable1b.csv` | `chapter1.py` |
| Supplementary Table 1c | `chapter1/supptable1c.csv`, `supptable1c_staining.csv`, `supptable1c_tcell_density.csv` | `chapter1.py` |
| Supplementary Table 2 | `chapter3/supptable3.csv` | `chapter3.py` |
| Supplementary Table 3a | `chapter4/supptable4a.csv`, `supptable4a_nichearea.csv` | `chapter4.py` |
| Supplementary Table 3b | `chapter4/supptable4c.csv`, `supptable4c_density.csv` | `chapter4.py` |
| Supplementary Table 3c | `chapter4/supptable4e.csv` | `chapter4.py` |
| Supplementary Table 3d | `chapter4/supptable3d.csv` (the printed cells), `ch4_niche_survival.csv` (every model) | `chapter4.py` |
| Supplementary Data 1 | `supplementary/Supplementary_Data_1.xlsx`, from `chapter2/supptable2_programs.csv` (sheet Thioredoxin programs), `supptable2_cancer.csv`, `supptable2_cycle_genes.csv` (Thioredoxin cycle genes), `supptable2_txn_txnip.csv` (TXN and TXNIP), `supptable2b_programs.csv` (Gene programs), `tnk_subtypes.csv` (T and NK subtypes) and `supptable2b_genes.csv` (Gene lists) | `chapter2/08_supptable2.py`, `10_supplementary_data1.py`, `12_tnk_subtypes.py` |

Numbers printed only in the text come from tables in the same folders; `manifest/manifest.csv` names the table and
the row of each (section 8).

### The imaging tables of mif/

`mif/` writes, per cohort, the tables that every imaging analysis reads, under `work_root/intermediate/{BOMI2,BOMI1}/`.
After the exclusions, Uppsala II holds 488 cores of 294 patients and Uppsala I 338 cores of 216 patients; the survival
analyses hold 269 Uppsala II patients (182 deaths) and 194 Uppsala I patients (125 deaths).

- `cells.tsv`, one row per cell of the Tumor or Stroma tissue category of a primary core: `core` (the core identifier
  of the export), `Case` (the patient identifier), `tissue` (Tumor or Stroma), `x_raw`, `y_raw` (the position as
  exported: pixels of 0.4964 um in Uppsala II, micrometers in Uppsala I; `redox.cells.load_cells` converts them once),
  `CD3`, `CD68`, `ASMA`, `CK` (0/1 positivity at the cohort threshold; Uppsala I has CD3 and CK only), `lineage` (T,
  Mac, CK, Fib or Other; a cell positive for more than one lineage marker is assigned to the first lineage in the
  priority order CD3, CD68, CK, alpha-SMA) and `TrxR1` (0/1 cytoplasmic TrxR1 positivity; empty when the intensity is
  missing). Channels not analyzed in the study are not read.
- `clinical.csv` (de-identified): `Case`, `histology` (LUSC, LUAD, other), `OS_weeks`, `event` (1 dead, 0 alive),
  `stage` (1 to 4 for I to IV) and `os_source`. `clinical_full.csv` keeps every merged column and is never shared.
- `thresholds.csv` (marker, cut, method, number of cells) and `areas.csv` (core, tissue, area in mm2).
- `percase.csv`, one row per patient: the number of stromal cells and of stromal cells with a TrxR1 call, the
  TrxR1-positive fraction of stromal cells with a call, the proportion of each lineage among all stromal cells, the
  stromal cells of each lineage per mm2 of stroma, and the clinical columns.

## 7. Figures

The figure pages of the article are composed from the panels and tables of this repository by separate layout code,
which sets the page size, letters and fonts and adds no values. The schematics (Figures 1a and 2b) and the micrograph
(Supplementary Figure 1g) are not made by code.

## 8. The value check

`manifest/manifest.csv` lists every value the article prints, one row per value: where it is printed
(`location`), what it is (`item`), the printed value (`printed`), the output table (`file`, relative to
`work_root/output/`), the row and column (`selector`), the rounding rule (`rule`) and, where a difference is expected,
its reason (`expected_diff`). `python3 manifest/check.py` formats every value from the tables and compares it with the
printed value. A row whose `expected_diff` gives a reason is reported as a known difference and does not fail the
check. Most rows read tables of the imaging part, so the full check needs the restricted data; the rows of the
single-cell part run from public data.

On the machine that made the article's outputs the check gives `PASS: 0 mismatches, 1 known differences`. On an
Intel x86_64 machine the label shuffle of `chapter4.py` ranks cells by distance, and equal distances round to
different last bits, so two p values of Supplementary Table 3b can move by 0.001 (manifest rows ch4-216 and ch4-241
carry this as a known difference); tied nearest neighbors can move Moran's I of a few cores of Figure 4a in the fourth
decimal without changing any p; and the example pair of Figure 4c can differ, because numpy builds sort tied values
differently. Every other value differs at most in its last digits.

## 9. Data availability and ethics

Imaging, clinical and mutation data of Uppsala I and Uppsala II are available from the authors on request. Public
data came from GEO (GSE81089, GSE325414), ArrayExpress (E-MTAB-6149, E-MTAB-6653, E-MTAB-13526), the LuCA atlas in CZ
CELLxGENE Discover and the data portal of Zhang and colleagues (PMID 35027529). `atlas/sources.tsv` lists every public
file with its URL, size, sha256, version and license. The repository holds code and small hand-made tables only, no
expression data and no patient data. The data portal of Zhang and colleagues states no license; nothing derived from
it is included here.

Gene-list sources in `chapter2/data/sources/` (compressed JSON):

- `msigdb_v2026.1.Hs_resource_sets.json.gz`: the MSigDB responses (release v2026.1.Hs) for the database sets of the
  catalog, with Reactome participant lists used for a cross-check. MSigDB gene sets are shared under the MSigDB terms
  of use; the KEGG lysosome set (KEGG_LYSOSOME, from the MSigDB "KEGG legacy" collection) is covered by KEGG's own
  license. Every set can be fetched again from gsea-msigdb.org by the systematic names in
  `gene_programs_definitions.csv`.
- `literature_lists.json.gz`: the gene lists taken from publications, each with two independent extractions and an
  adjudication, the gene symbols as printed and where in the publication the list is.
- `pubmed_hgnc_responses.json.gz`: the PubMed and HGNC responses for the citations and the symbol reconciliation, so
  a rebuild needs no network (`chapter2/01_programs.py --allow-network` fetches anything missing).

A rebuilt atlas differs from the object analyzed in the article in one place: the GSE325414 arm reads the first entry
of the deposit matrix (row 19, AAGAB, first cell, count 1), which the analyzed object lacked. That cell is epithelial,
so it is not in `fullgene_object.h5ad`; in `panel_object_v2.h5ad` it changes the AAGAB detection of one cell, which
cannot change the block QC verdict.

The study was conducted in accordance with the Declaration of Helsinki and approved by the Regional Ethical Review
Board in Uppsala (2006/325 for Uppsala I and 2012/532 for Uppsala II), which waived the requirement for written
informed consent.

## 10. License, citation and contact

The code is released under the MIT license (`LICENSE`).

When you use this code, please cite the article: "Stromal thioredoxin reductase 1, T-cell-rich stroma and survival in
lung squamous cell carcinoma" (citation to be added on publication). The code is at https://github.com/ArturMez/stromal-trxr1-nsclc.
`CITATION.cff` gives the citation of the code in machine-readable form, and `CHANGELOG.md` lists the versions.

Contact: the corresponding author, Alfonso Martín-Bernabé, Department of Oncology-Pathology, Karolinska Institutet,
Stockholm, Sweden (alfonso.martin.bernabe@ki.se).
