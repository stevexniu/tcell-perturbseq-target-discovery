# Novel druggable regulators of CD4+ T-cell activation

Reproducible target-discovery pipeline over the Marson 2025 genome-scale CRISPRi
Perturb-seq dataset in primary human CD4+ T cells. Two discovery axes are built
**from the data itself**, each **validated by reproducing the authors' own
regulator coefficients** before any gene is nominated, and candidates are worked
up on four independent evidence dimensions (novelty, druggability, human
genetics, safety) against public databases.

## One-command rerun

```bash
conda env create -f env/environment.yml
conda activate tcell-target-discovery
./run.sh            # regenerate all figures from the shipped tables (fast, no downloads)
```

Deeper reruns:

```bash
./run.sh full       # + re-extract per-condition scores over S3 and re-fetch GTEx
./run.sh discovery  # + re-run the genome-scale dual-axis discovery (61 GB inputs, ~30 min)
```

(No conda? `python -m venv .venv && . .venv/bin/activate && pip install -r env/requirements.txt`.)

## The result in one paragraph

The **polarization axis** (Th1/Th2) passes validation strongly — Spearman ρ 0.72
against the authors' coefficients, AUROC 0.92 on their known regulators, 100%
sign agreement — so all nominations are drawn from it. The **activation axis
fails** (ρ 0.13, AUROC ≈ chance) and everything from it is exploratory only. The
central finding is a **novelty-versus-actionability tradeoff**: the most
reproducible genuinely-new hits are mostly intracellular scaffolds / RNA-binding
proteins / housekeeping enzymes (hard to drug, broad expression), while the
tractable, safe, genetically-supported candidates are the ones already known and
in the clinic — which the axis re-finds in the correct direction, its proof of
correctness. Best genuinely-new lead: **TATDN2** (Th2-skewing, TatD-family
nuclease active site, most T-cell-specific in its direction; caveat: no genetic
support yet). Second: **JADE2** (Th1-skewing, chromatin reader).

## Repository layout

```
tcell_target_discovery/
├── README.md                 # this file
├── DATA.md                   # public-data pointers (S3 bucket, authors' tables, evidence APIs)
├── run.sh                    # one-command rerun
├── env/                      # environment.yml, requirements.txt, modal_image.txt
├── scripts/
│   ├── discover.py           # validated building blocks: axis construction, projection, validation
│   ├── driver.py             # genome-scale dual-axis discovery run (Modal / any 8-CPU host)
│   ├── extract_per_condition.py  # per-condition polarization scores (reads ~42 rows over S3)
│   ├── fetch_evidence.py     # evidence-source documentation + runnable GTEx/HPA fetch
│   └── make_figures.py       # regenerate all figures from results/tables/
└── results/
    ├── reports/
    │   ├── summary_report.md            # main report: method, nominations, tradeoff, falsifiers
    │   └── candidate_workup_report.md   # detailed per-candidate workup with all identifiers
    ├── tables/               # all nomination + evidence tables (see below)
    └── figures/              # Fig 1–3 + coverage heatmap
```

## Key outputs

| File | What it is |
|---|---|
| `results/reports/summary_report.md` | main deliverable — validated method, evidence-ranked shortlist per direction, novelty/actionability tradeoff, one falsifiable experiment per lead |
| `results/tables/coverage_table.csv` | 15 candidates × 34 evidence columns (every gap explicit) |
| `results/tables/ranked_polarization_Stim8hr_full.csv` | full ranked polarization axis (the trusted one) |
| `results/tables/ranked_activation_Stim8hr_full.csv` | full ranked activation axis (exploratory) |
| `results/tables/{novelty,druggability,genetics,safety}_read.csv` | the four specialist reads |
| `results/figures/validation_scatter.png` | Fig 1 — score vs authors' coefficients (both axes) |
| `results/figures/fig2_ranking.png` | Fig 2 — ranked candidates, nominations highlighted, hubs greyed |
| `results/figures/fig3_leads.png` | Fig 3 — leads' GTEx expression, HPA specificity, per-condition KD effect |
| `results/figures/coverage_heatmap.png` | evidence coverage across candidates × dimensions |

## Method notes

- **Two axes, both from the data.** Activation = non-targeting-control pseudobulk,
  Stim8hr vs Rest CPM log-ratio. Polarization = the authors' published Th1/Th2
  signature (Ota 2021 contrast **only** — see the footgun below). Each perturbation
  is projected onto the standardized signature.
- **Validate before trusting.** Each axis's per-perturbation score is correlated
  against the authors' own regulator coefficients on the overlapping genes, and
  AUROC / sign-agreement are computed on their known regulators. Only the
  polarization axis clears the bar.
- **Nomination cutoff (documented).** QC-pass (on-target KD significant, two-guide
  agreement ≥ 0.70, ≥ 20 downstream genes, no off-target/neighboring flag),
  non-hub, not cell-line-shared (authors' K562 comparison), and absent from the
  authors' regulator table (`novel`).
- **Known vs new is derived, not typed.** Every candidate is marked
  recovered-known / partial / genuinely-new against the authors' own table and a
  literature/ChEMBL/trials search — never from a hand-picked gene list.

### Footgun (documented)

The published polarization signature file contains **two** contrasts —
`Th2_vs_Th1 (Ota 2021)` and `Th2_vs_Th1 (Hollbacker 2021)`. A
`str.contains("Th2_vs_Th1")` filter matches both and, via dict-overwrite,
silently substitutes Hollbacker values for the 11,616 overlapping genes (the two
correlate only r≈0.69). This pipeline filters on the exact Ota string. See
`scripts/discover.py` / `scripts/driver.py`.

## Data

All inputs are public; see `DATA.md`. The heavy discovery step downloads 61 GB
and was run on Modal (8 CPU / 32 GB, ~30 min) — the recipe is in
`env/modal_image.txt` and the `driver.py` header, but it runs on any host with
enough disk.
