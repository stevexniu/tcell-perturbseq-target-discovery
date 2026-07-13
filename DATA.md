# Data sources (all public)

## 1. Perturb-seq differential-expression table (primary input)
- **S3 bucket:** `s3://genome-scale-tcell-perturb-seq/marson2025_data/`
- **DE table:** `GWCD4i.DE_stats.h5ad` (16.8 GB) — anonymous HTTPS mirror:
  `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad`
- **Pseudobulk (for the activation signature):** `GWCD4i.pseudobulk_merged.h5ad` (44.6 GB), same mirror.
- **Dataset page:** https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- **Paper:** Genome-scale perturb-seq in primary human CD4+ T cells, bioRxiv 2025, DOI 10.64898/2025.12.23.696273

The scripts read the DE/pseudobulk files lazily over HTTP range reads (h5py + fsspec) — no full download needed.
NOTE: fsspec's http filesystem must be created with `client_kwargs={"trust_env": True}` behind a proxy.

## 2. Authors' own analysis tables (for Figure 1 validation)
- **Repo:** https://github.com/emdann/GWT_perturbseq_analysis_2025 (branch `master`)
- **Regulator coefficients:** `metadata/suppl_tables/polarization_prediction_condition_comparison_regulator_coefficients.csv`
  (we use rows with `signature == "activation"`, `celltype == "Stim8hr"`).

## 3. Annotation connectors (Figure 3)
- **GTEx v8** tissue median TPM — via the Expression MCP connector (`gtex_expression_summary`).
- **Human Protein Atlas** (release 25.1) — via the Protein Annotation MCP connector (`get_protein_atlas_gene`).
- **GWAS Catalog** — MAP3K1 associations verified live (`gwas_associations_for_gene`).

Pre-computed copies of every table these produce are in `data/` so the figures rebuild offline.
