#!/usr/bin/env bash
# One-command rerun of the reproducible parts of the pipeline.
#
#   ./run.sh          # regenerate figures from the shipped tables (fast, no downloads)
#   ./run.sh full     # also re-extract per-condition scores over S3 + re-fetch GTEx
#   ./run.sh discovery# also re-run the genome-scale dual-axis discovery (needs 61 GB / ~30 min)
#
# Environment: create it once with
#   conda env create -f env/environment.yml && conda activate tcell-target-discovery
# (or: python -m venv .venv && . .venv/bin/activate && pip install -r env/requirements.txt)

set -euo pipefail
cd "$(dirname "$0")"
MODE="${1:-figures}"

echo "== [1/3] regenerating figures from results/tables/ =="
python scripts/make_figures.py

if [[ "$MODE" == "full" || "$MODE" == "discovery" ]]; then
  echo "== [2/3] re-extracting per-condition polarization scores (reads ~42 rows over S3) =="
  python scripts/extract_per_condition.py
  echo "== re-fetching GTEx focus tissues (public API) =="
  python scripts/fetch_evidence.py
  echo "== regenerating figures with refreshed tables =="
  python scripts/make_figures.py
fi

if [[ "$MODE" == "discovery" ]]; then
  echo "== [3/3] genome-scale dual-axis discovery =="
  echo "   downloads GWCD4i.DE_stats.h5ad (16.8 GB) + GWCD4i.pseudobulk_merged.h5ad (44.5 GB)"
  echo "   ~30 min on 8 CPU / 32 GB. See scripts/driver.py header for the Modal recipe."
  python scripts/driver.py
fi

echo "== done. Figures in results/figures/, tables in results/tables/, reports in results/reports/ =="
