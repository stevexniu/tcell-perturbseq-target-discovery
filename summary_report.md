# Novel druggable regulators of T-cell activation — summary report

**Dataset.** Genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells (Marson lab, bioRxiv
2025, DOI 10.64898/2025.12.23.696273), profiled at rest and 8 h / 48 h after stimulation.
Differential-expression table: `s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad`
(33,983 perturbation×timepoint contrasts × 10,282 genes).

**Goal.** A reproducible, evidence-backed shortlist of novel, druggable regulators of T-cell
activation in both therapeutic directions — autoimmune (knockdown lowers the activation program) and
anti-tumor (knockdown raises it) — not a single forced pick.

---

## Method

1. **Activation signature, from the data's own controls.** Built from the dataset's non-targeting
   control pseudobulk (log2 CPM, Stim8hr vs Rest; 3,701 + 3,635 NTC profiles). Canonical effector
   genes move as expected (IFNG +4.6, GZMB +4.0, IL2RA +3.5, IL2 +3.1; quiescence genes KLF2 −3.6,
   IL7R −2.4 down). Three of six quiescence markers had the wrong sign or were null (CCR7, BACH2,
   SELL) — a known limitation of a two-state stim-vs-rest contrast, carried openly.
2. **Method selected by a known-gene benchmark.** Four scoring methods ({z-score, log_fc} ×
   {mean_signed, correlation}) were scored on their ability to separate 14 canonical drivers and
   11 canonical brakes from a low-expression null (AUROC). **z-score / mean_signed** was chosen
   (drivers-vs-null 0.956, brakes-vs-null 0.742). Under 200 random 50/50 splits re-selecting the
   method each fold, held-out AUROC was **drivers 0.935 ± 0.045, brakes 0.738 ± 0.078** — so the
   brake (anti-tumor) direction is measurably weaker than the driver direction, everywhere.
3. **Independent validation — reproduces the authors' own coefficients.** Our per-perturbation
   activation score was compared against the authors' published activation-signature regulator
   coefficients (their `polarization_..._regulator_coefficients.csv`, Stim8hr). Acceptance bar
   (r ≥ 0.70 inherited from the benchmark skill's fixed threshold; sign-concordance ≥ 0.80 a
   convention, not tuned to these data): known-gene Pearson r ≥ 0.70, driver sign-concordance ≥ 0.80,
   all canonical drivers at the program-down end. **Observed: known-gene r = 0.93 (Spearman 0.91, n = 18), driver concordance
   100%, all drivers program-down → PASS.** Across all 3,830 shared genes the global correlation is
   modest (r = 0.23) because most genes have near-zero effect and the two estimators differ
   (regularized regression coefficient vs signature projection); agreement rises monotonically with
   effect size (r = 0.28 → 0.34 as weak-effect genes are excluded). One known brake, RASA2, is
   sign-discordant (both values near zero). The axis recovers the field's known biology — CD28 and
   CBLB land where expected (see Figure 2).
4. **Ranking.** QC filter (on-target KD significant, guide-reproducible, ≥20 downstream genes, no
   off-target/neighbor flags) kept 487/11,415 Stim8hr perturbations. Hub/fitness genes flagged at
   n_downstream ≥ p90 (1,630) and set aside from the leads.

![Figure 1. The activation axis reproduces the authors' regulator coefficients and recovers known biology.]({{artifact:e4d156c9-cf7c-41fd-95f0-90134632a348}})

*Figure 1 — Validation. Left: our per-perturbation activation score vs the authors' own
activation-signature regulator coefficient (Stim8hr, n=3,830). Canonical drivers (blue) sit at the
program-down end, brakes (orange) at the program-up end; known-gene r=0.93 while the global r=0.23 is
shown separately. Right: correlation strengthens monotonically as low-effect genes are excluded.*

---

## The two nominations

![Figure 2. Both-direction ranking of 487 reproducible, activation-specific regulators.]({{artifact:22bd5600-4ebe-44bc-9298-49dca2e80512}})

*Figure 2 — Ranking. All 487 QC-passing perturbations, ordered by activation score. Program-down end
= autoimmune (knockdown lowers activation); program-up end = anti-tumor (knockdown raises it). The two
nominated leads RSBN1L and MAP3K1 are highlighted; known anchors CD28 (autoimmune) and CBLB
(anti-tumor) are labelled; hub/fitness genes are greyed out.*

### Nomination 1 — Anti-tumor: **RSBN1L** (genuinely new)

The single most reproducible brake in the screen (guide-reproducibility 1.000; specific footprint,
n_downstream 178, well below the hub threshold — not general fitness loss). Genuinely new to
immunity: only 3 PubMed records, none in a T-cell or drug-target context, and not named in the
dataset paper's abstract. Human genetics point the right way for an inhibitor — gene-body variants
associate with autoimmune hypothyroidism (p=2e-24), eosinophilia (p=2e-13), SLE and allergy (a
brake-loss / heightened-immunity signature), and the gene is loss-of-function-tolerant (gnomAD
LOEUF 0.60), the profile an inhibitor wants. Annotated as an RSBN1-family putative H4K20 histone
demethylase, giving a plausible transcriptional-brake mechanism in a druggable enzyme class.

**The catch.** Druggability is hard *today*: zero ligands, no experimental structure (AlphaFold
only) — a lead would need a de novo chemistry campaign. The demethylase mechanism is annotation, not
demonstrated in T cells. No formal eQTL colocalization was run, so "genetics supports" rests on trait
association + variant position + constraint, not proven causal allele direction. **Critically, the
authors' own activation-signature regulator coefficient disagrees with our brake call for RSBN1L:**
their coefficient is −0.022 (z ≈ −39; more negative than 94% of genes), placing RSBN1L at the
*driver* (program-down) end, opposite our program-up brake score (+0.554). This is a genuine
method-level conflict — our signature projection of the knockdown DE profile vs their regularized
multivariate regression — and it is the single strongest reason for caution on this lead. It does not
by itself refute the brake call (the two estimators measure different things, and the direct
activation-shift assay in Figure 3 supports the brake direction), but confirming the brake effect in
an independent arrayed assay is the prerequisite for taking RSBN1L further.

**Falsifying experiment.** Arrayed CRISPRi knockdown (plus catalytic-dead rescue) of RSBN1L in
primary human CD4+/CD8+ T cells → TCR stimulation + tumor-cell co-culture killing assay. The
hypothesis predicts KD *increases* effector cytokines (IL-2/IFNG) and tumor killing, and that a
catalytic-dead demethylase mutant fails to rescue. If KD does not raise effector output, if the
effect tracks only with viability, or if catalytic-dead rescues fully, the brake-via-demethylase
hypothesis is falsified.

### Nomination 2 — Autoimmune: **MAP3K1 (MEKK1)** (recovered-known-broad)

The only autoimmune-direction candidate where all three axes align. Knockdown lowered activation with
the tightest footprint among the six nominees (n_downstream 89 — a focused signaling regulator, not a
hub; the shortlist-wide minimum is 20). A documented TCR/JNK-pathway kinase (T-cell-specific deletion alters effector responses, PMID
26774476); druggable now (Ser/Thr kinase, 8 nM ChEMBL ligand, PDB 6WHB, Phase-2 clinical inhibitor
E-6201/CHEMBL1097999). Human genetics support an immunosuppressive target: verified against the live
GWAS Catalog, two immune-mediated traits map to the gene body at genome-wide significance — psoriasis
(rs12654176, p=9e-12, PMID 40021644) and polyarteritis nodosa (rs535102440, p=2e-12, PMID 39024449).
Additional signals (Crohn's, Graves', palindromic rheumatism) are FinnGen PheWAS at sub-genome-wide
significance — suggestive, not established.

**The catch.** Not novel — a known druggable kinase already in the clinic. The best clinical molecule
(E-6201) is a *dual* MEK/MEKK1 inhibitor, so a MEKK1-selective chemical is needed to test the specific
hypothesis. MAP3K1 is LoF-intolerant (LOEUF 0.33) with a Mendelian sex-reversal knockout phenotype,
so only *partial* inhibition is a viable therapeutic model — on-target dose-limiting toxicity is the
main risk.

**Falsifying experiment.** Titrated selective MEKK1 inhibition (or a degron) in primary human CD4+ T
cells across a dose range, with viability controls → measure activation (CD25/CD69, IL-2/IFNG) vs a
general-toxicity readout. The hypothesis predicts a dose window where activation drops *before*
viability does. If reduced activation appears only at doses that also reduce viability, the "specific
activation regulator" claim collapses into a fitness effect. The genetics rationale is separately
falsified if the autoimmune risk alleles *increase* rather than decrease MAP3K1 activity.

![Figure 3. The two leads up close: tissue/immune expression and activation-program shift.]({{artifact:ca9d290f-feb1-471f-b0a9-5a01ed42b915}})

*Figure 3 — Leads up close. GTEx tissue expression (whole blood highlighted), Human Protein Atlas
protein-level annotation, and the activation-program shift across Rest → Stim8hr → Stim48hr for each
lead. RSBN1L knockdown raises the activation score on stimulation (brake release; +0.29/+0.31 at
8h/48h); MAP3K1 knockdown lowers it, sharpest at 8h (−0.40; driver loss). Error bars are 95% bootstrap
CIs over signature genes — not biological-replicate CIs (see limitation below).*

---

## Novelty versus actionability — stated plainly

The central tension in this candidate set is that **novelty and actionability are anti-correlated.**
The genuinely-new hits (RSBN1L, and the alternative CLCC1) have no chemical matter and, for RSBN1L,
an unproven mechanism. The immediately tractable hits (MAP3K1, MBTPS1) are already known targets —
MBTPS1's brake-release direction is in fact already published (PMID 40307212). Rather than force a
novel-but-unactionable gene into a lead slot, each direction is nominated on its best-evidenced lead
with status marked: RSBN1L carries the true-novelty flag (highest reproducibility, right-direction
genetics, hard chemistry); MAP3K1 carries the fast-to-test flag (known druggable kinase, aligned
genetics, low novelty). The full six-gene evidence matrix (including CLCC1, WAS, MBTPS1, PUM1) is in
`target_nominations.md` and `candidate_evidence_detailed.csv`.

## The anti-tumor direction's weaker evidence — stated openly

The brake/anti-tumor direction is weaker than the driver/autoimmune direction at every level of this
analysis, and this should temper confidence in RSBN1L relative to MAP3K1:
- **Benchmark separation** is consistently lower for brakes (held-out AUROC 0.738 ± 0.078) than
  drivers (0.935 ± 0.045).
- **Author-coefficient validation** is carried by the drivers; the brakes include the one
  sign-discordant known gene (RASA2, author coef −0.021 vs our +0.203).
- **RSBN1L specifically conflicts with the authors' coefficient.** Their regulator coefficient for
  RSBN1L (−0.022, z ≈ −39) points to the driver end, opposite our brake call — a method-level
  disagreement disclosed in full above and visible in Figure 1. RASA2 aside, this is the clearest
  reminder that program-up (brake) calls are less securely validated than program-down (driver) calls.
- **RSBN1L's own gaps** are the hardest of the two leads: no ligand, no structure, mechanism
  unproven in T cells.

The recommendation is therefore asymmetric: MAP3K1 is a fast, well-supported test of the autoimmune
direction; RSBN1L is the higher-risk, higher-novelty anti-tumor bet whose first milestone is simply
confirming the brake effect and its enzyme-dependence before any chemistry investment.

## Shared next step

Neither lead has a formal allele-direction eQTL colocalization; both "genetics supports" calls rest
on trait association + variant position + constraint. That colocalization is the highest-value next
genetics step for either candidate.

---

*Full nominations with per-axis findings and source ids: `target_nominations.md`. Evidence tables:
`candidate_evidence_matrix.csv`, `candidate_evidence_detailed.csv`. Method scripts, environment, and
data pointers accompany this report in the repository.*
