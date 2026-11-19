# atlas/data: tables the build reads

Both tables were prepared by hand and have no producing script. They ship here because the build
cannot run without them and they hold no patient-level data. They keep the spelling of the public
deposits (for example the obs value `primary_tumour`), so they use British spellings in places.

## build-plan_obs_schema.csv

The obs schema every arm must fill (the 30 obs columns of the atlas). Read by
`steps/stage1a_arm_emtab6653.py` and `steps/stage2_emtab13526_15_verify_arm.py`, which use only the
columns `column` and `required`.

| column | meaning |
|---|---|
| column | obs column name |
| dtype | intended type |
| allowed_values_or_form | the vocabulary or the form of the value |
| required | yes or no |
| purpose | why the column exists |
| trap_it_defends_against | the error in a public deposit that the column guards against |

## arm_registry.tsv

The arm registry, read by `steps/fullbuild_state.py` (columns arm_id, status, reason, cells_full). The
LuCA row states that the object holds 16 LuCA blocks from 12 studies, and its reason says "12 LuCA
studies reach the object". No column that the build reads for an included arm depends on those
words (fullbuild_state.py prints the reason of excluded arms only).

| column | meaning |
|---|---|
| arm_id | arm name as in `varm/gene_present` of the full-gene object |
| source_key | obs source_key of the arm |
| block_id | blocks the arm contributes (free text) |
| status | included or EXCLUDED |
| reason | why (free text) |
| panel_genes, full_genes | genes on the panel axis and in the full arm |
| cells_full | cells the arm contributes to fullgene_object.h5ad |
| counts | raw, recovered or estimated counts (free text) |
| block_qc | block QC verdict (free text) |
| notes | free text |
