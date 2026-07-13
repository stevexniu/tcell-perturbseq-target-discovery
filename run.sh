#!/usr/bin/env bash
set -euo pipefail
echo "[1/2] Rebuilding figures from data/ (offline)..."
python scripts/04_make_figures.py
echo "[2/2] Re-fetching authors' coefficients to reconfirm Fig1 validation..."
python scripts/03_fig1_authors_coef.py || echo "  (network fetch skipped/failed; shipped data/fig1_score_vs_authorcoef.csv still valid)"
echo "Done. See figures/ and summary_report.md"
