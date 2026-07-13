# Novel druggable regulators of T-cell activation — CD4+ Perturb-seq target discovery

Reproducible target-discovery pipeline over the genome-scale CRISPRi Perturb-seq screen in primary
human CD4+ T cells (Marson lab, bioRxiv 2025, DOI 10.64898/2025.12.23.696273). It builds an
activation signature from the dataset's own controls, validates the scoring method against known
T-cell regulators **and** against the authors' own published regulator coefficients, ranks all
perturbations in both therapeutic directions, and nominates one lead per direction with novelty,
druggability, and human-genetics evidence.

## Nominations (summary)

| direction | lead | status | why | main catch |
|---|---|---|---|---|
| **Anti-tumor** (KD raises activation) | **RSBN1L** | genuinely new | most reproducible brake (guide-repro 1.000), right-direction autoimmune genetics, LoF-tolerant | no ligand/structure; mechanism unproven; authors' own coefficient disagrees with the brake call |
| **Autoimmune** (KD lowers activation) | **MAP3K1** (MEKK1) | recovered-known-broad | all three axes align; druggable kinase (Phase-2 inhibitor); genome-wide autoimmune loci | not novel; needs a MEKK1-selective molecule; LoF-intolerant (partial inhibition only) |

Central tradeoff: **novelty and actionability are anti-correlated** — the novel hits lack chemical
matter; the tractable hits are known targets. Full evidence and caveats: [`summary_report.md`](summary_report.md)
and [`target_nominations.md`](target_nominations.md).

## One-command rerun

```bash
conda env create -f environment.yml && conda activate tcell-perturbseq-targets
bash run.sh
```

`run.sh` rebuilds all three figures from the shipped `data/` tables (fully offline, ~10 s) and
re-fetches the authors' regulator coefficients to reconfirm the Figure-1 validation. The heavy
stages (signature build, DE scoring, benchmark) read the 17-45 GB h5ad files lazily over HTTP; their
outputs are checkpointed in `data/` so the figures and report reproduce without the large reads. See
[`DATA.md`](DATA.md) for every data pointer and the methods notes in each script.

## Layout

```
summary_report.md              # main deliverable: method, 2 nominations, tradeoff, experiments
target_nominations.md          # full per-axis findings with source ids
DATA.md                        # public data pointers (S3 + authors' GitHub + connectors)
environment.yml
run.sh                         # one-command figure rebuild + Fig1 re-validation
scripts/
  01_build_signature.py        # activation signature from NTC controls (methods)
  02_benchmark.py              # known-gene AUROC method selection + held-out validation (methods)
  03_fig1_authors_coef.py      # RUNNABLE: fetch authors' coefficients, rebuild Fig1 validation table
  04_make_figures.py           # RUNNABLE: rebuild all 3 figures from data/
figures/  fig1_validation.png fig2_ranking.png fig3_leads.png
data/     signature.csv benchmark_auroc.csv benchmark_validation.json
          shortlist_{autoimmune,antitumor}.csv candidate_evidence_{matrix,detailed}.csv
          validation_criteria.json fig1_score_vs_authorcoef.csv fig2_ranking.csv
          fig3_{activation_shift,gtex,hpa}.csv
```

## Method in brief

1. **Signature** from the dataset's own non-targeting controls (log2 CPM Stim8hr vs Rest).
2. **Method selection** by a known-gene AUROC benchmark (14 drivers, 11 brakes vs a low-expression
   null); z-score/mean_signed chosen. 200 random 50/50 splits re-selecting the method each fold give
   held-out AUROC **drivers 0.935±0.045, brakes 0.738±0.078** — the anti-tumor direction is weaker.
3. **Independent validation**: our per-perturbation score reproduces the authors' own activation
   regulator coefficients (known-gene Pearson **r=0.93**; global r=0.23 as most genes are noise), and
   recovers CD28/CBLB where expected.
4. **Ranking**: QC filter keeps 487/11,415 Stim8hr perturbations; hub/fitness genes flagged and set
   aside from the leads.

## Caveats (carried openly)

- The **anti-tumor/brake direction is weaker** than the autoimmune/driver direction at every level
  (benchmark AUROC, author-coefficient validation, and RSBN1L's own gaps).
- **RSBN1L conflicts with the authors' coefficient** (their −0.022 places it at the driver end); the
  brake call rests on our signature projection + the direct activation-shift assay and must be
  confirmed in an independent arrayed assay.
- No formal **eQTL colocalization** was run for either lead — the highest-value next genetics step.
- "Recovered-known" genes come from the benchmark set and are **not** independent evidence.

Data are public and de-identified; nothing here is clinical advice.
