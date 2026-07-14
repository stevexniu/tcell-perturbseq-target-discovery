# Public data pointers

All inputs are public. Nothing in this repo requires credentialed access.

## Primary dataset — Marson 2025 genome-scale CRISPRi Perturb-seq (CD4+ T cells)

- **Dataset page:** https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- **Public S3 bucket (anonymous / requester-pays-free):** `s3://genome-scale-tcell-perturb-seq/marson2025_data/`
  - HTTPS mirror: `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/<file>`

| File | Size | Used for |
|---|---|---|
| `GWCD4i.DE_stats.h5ad` | 16.8 GB | per-perturbation DE statistics (6 layers: log_fc, p_value, adj_p_value, baseMean, lfcSE, zscore). Scored to build both axes; the per-condition extraction reads ~42 rows from it directly. |
| `GWCD4i.pseudobulk_merged.h5ad` | 44.5 GB | non-targeting-control pseudobulk; builds the activation axis (Stim8hr vs Rest CPM log-ratio). Only the 11,018 non-targeting rows are used. |
| `data_sharing_readme.md` | — | authoritative schema for the h5ad files. |

Structure of the DE table: 33,983 obs (perturbation × condition) × 10,282 measured genes.
Key obs columns: `target_contrast` (Ensembl ID), `target_contrast_gene_name` (symbol),
`culture_condition` (Rest / Stim8hr / Stim48hr), `ontarget_significant`,
`guide_correlation_signif`, `n_downstream`, plus off-target/neighboring flags.

## Authors' published supplementary tables (validation ground truth)

Base: `https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/metadata/suppl_tables/`

| File | Used for |
|---|---|
| `polarization_prediction_condition_comparison_regulator_coefficients.csv` | authors' own regulator coefficients — the validation target for both axes (`signature` ∈ {activation, ota}; `known_regulators` flag). |
| `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` | the published Th1/Th2 signature. **Use `contrast == "Th2_vs_Th1 (Ota 2021)"` only** — the file also contains a Hollbacker-2021 contrast, and a `str.contains("Th2_vs_Th1")` filter silently mixes the two (they correlate only r≈0.69). |
| `K562_comparison.suppl_table.csv` | cell-line-shared / generic effect flag (`logfc_pearson_r`). |

## Evidence resources (per-candidate workup)

Public APIs/portals queried per gene (see `scripts/fetch_evidence.py`):
GTEx v8, Human Protein Atlas, ChEMBL, RCSB PDB, AlphaFold DB, UniProt/InterPro,
GWAS Catalog, eQTL Catalogue, Open Targets (Genetics), gnomAD v4,
ClinicalTrials.gov v2, OpenAlex/PubMed, Drugs@FDA.
