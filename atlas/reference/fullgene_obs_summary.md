# fullgene_obs_summary.tsv

The obs summary of `output/RNA/meta_dataset/fullbuild/fullgene_object.h5ad`, read with both
corrections of `redox/atlas_fixes.py`, as every analysis reads it: the study label of GSE325414 is
Jimenez_Bruno_2026, and the 50,651 cells of the redistributed copies of eight Lambrechts 2018 patients
(Zhang_Chen_2022_DL1 to DL8) are left out of every section and counted once in `left_out_copies`. Counted:
1,838,259 cells, 368 patients, 24 blocks, 17 studies; 1,722,161 cells and 322 patients pass block QC (the
object itself holds 1,888,910 cells and 376 patient identifiers). The file was written with
`python atlas/05_obs_summary.py --object <object> --write-reference` (h5py, obs only). A rebuilt object
must reproduce every row: `python atlas/05_obs_summary.py`. sha256 of the file: in the atlas README.

Tab separated, one row per group:

| column | meaning |
|---|---|
| section | total, block_id, study, source_key, histology, tissue_class, level1_compartment, level2_celltype (key level1 and level2 joined by a bar), block_qc_pass, chapter2_selection (key histology and level1, or histology and "all"; cells of primary tumor that pass block QC, LUSC or LUAD), chapter2_selection_level2 (key histology, level1 and level2), left_out_copies (key "all"; the cells the other sections leave out) |
| key | the group; the object's value "primary_tumour" is written "primary_tumor" |
| n_cells | cells in the group |
| n_patients | distinct obs patient_id in the group |
| n_qc_pass | cells of the group with block_qc_pass True |

Values that tie to chapter2/, equal to the cells and patients of the working sets in
`work_root/intermediate/chapter2/work_{fib,mac,tnk}/meta_{LUSC,LUAD}.npz`: fibroblasts
6,808 LUSC + 5,168 LUAD = 11,976 (55 + 107 = 162 patients); macrophages 64,239 + 97,125 = 161,364 (67 + 148 =
215 patients); T/NK 155,689 + 215,534 = 371,223 (70 + 155 = 225 patients); the pooled non-cancer set before its
cap 289,299 + 387,179 = 676,478 cells (all compartments of the selection). Without the exclusion of the
Lambrechts copies the table would count them (fibroblasts 12,821, macrophages 164,486, T/NK 386,022,
376 patients).
