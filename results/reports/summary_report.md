# Novel druggable regulators of CD4+ T-cell activation — evidence-backed shortlist

**Dataset:** genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells (Marson 2025), profiled at rest and 8 h / 48 h post-stimulation. DE table `GWCD4i.DE_stats.h5ad` (33,983 perturbation×condition observations × 10,282 measured genes). Analysis condition: **Stim8hr**.

**Goal:** identify novel, druggable regulators of T-cell activation worth pursuing as therapeutic targets — in both an immunosuppressive direction (autoimmune) and a brake-release direction (anti-tumor).

---

## 1. Method — axes built from the data, validated against the authors' own numbers

Two discovery axes were built entirely from the dataset itself, never from a hand-typed gene list:

- **Activation axis** — the dataset's own non-targeting-control cells, stimulated (Stim8hr) vs resting, as a CPM log-ratio signature.
- **Polarization axis** — the authors' published Th1/Th2 signature (Ota 2021 contrast only; the discovery skill's default filter silently mixes in a second contrast, which was corrected).

Each per-perturbation score was then **validated by reproducing the authors' own regulator coefficients** (fetched live from their public supplementary tables) on the 3,830 perturbed genes that overlap the authors' table:

| Axis | Spearman ρ vs authors' coef | AUROC (known regulators) | Sign agreement (21 knowns) | Verdict |
|---|---|---|---|---|
| **Polarization** | **0.72** | **0.92** | **100%** | **PASS — trusted** |
| Activation | 0.13 | 0.49 (≈ chance) | 76% | exploratory only |

The polarization axis reproduces the authors' coefficients strongly and orders all 21 known regulators perfectly by direction (Th1 drivers IRF1/JAK2/IFNGR1/STAT1/2 positive; Th2 drivers GATA3/STAT6/IL4R/RARA negative). The activation axis fails formal validation — its coarse extremes still recover the canonical TCR module (CD247, PLCG1, ITK, MALT1), but 4 of those 5 are hub/broad-effect knockdowns, so it captures direction, not activation-specific signal. **All nominations are drawn from the validated polarization axis; activation-axis hits are exploratory support only.**

![Validation — per-perturbation score vs authors' coefficients]({{artifact:964c1e9c-4612-4de0-89cd-62bae8f4bed7}})

*Figure 1. The polarization axis (left) reproduces the authors' regulator coefficients (ρ 0.72) and recovers their known regulators (red) in the correct rank order; the activation axis (right) does not (ρ 0.13).*

**Recovering known biology before trusting anything new:** the top of both polarization directions is populated by the authors' own known regulators (Th1 end: RARA, GATA3, STAT6, CHD4, MTA2, PTPN2; Th2 end: JAK2, IFNGR1/2, STAT2, TRAF3), and the specialist workup independently re-finds clinically-drugged immune targets in the correct direction (below). This is the axis's proof-of-correctness.

---

## 2. Nomination pipeline

From the validated axis, candidates were filtered with an explicit, documented cutoff — QC-pass (on-target KD significant, two-guide agreement ≥ 0.70, ≥ 20 downstream genes, no off-target flag), non-hub, not cell-line-shared, and absent from the authors' regulator table (`novel`) — then worked up on four independent evidence dimensions via live database connectors: **novelty** (OpenAlex/PubMed, ChEMBL, Drugs@FDA, ClinicalTrials.gov), **druggability** (UniProt/InterPro, ChEMBL, PDB/AlphaFold), **human genetics** (GWAS Catalog, eQTL Catalogue, Open Targets), and **safety** (GTEx/HPA expression breadth, gnomAD constraint).

![Ranked candidates with nominations highlighted]({{artifact:79aa7033-f9ad-4789-b452-5fd9c42935b4}})

*Figure 2. All QC-pass perturbations ranked by polarization score. Nominations are highlighted by direction (blue = Th1-skewing, red = Th2-skewing), recovered-known anchors are labelled "(known)", and hub/fitness genes are greyed.*

![Leads up close — expression breadth vs polarization effect]({{artifact:223e8449-0bce-45e2-b26d-7913d442d30d}})

*Figure 3. For each lead: (A) GTEx tissue expression, (B) HPA RNA specificity, (C) knockdown effect on the polarization axis across Rest / Stim8hr / Stim48hr. Recovered-known targets (* ) are immune-restricted (TNFRSF9, CD6); the novel intracellular hits are broadly expressed. Knockdown direction is consistent across conditions for the polarization nominations; among the exploratory activation genes, CD6 flips sign by 48 h and INTS6/TATDN2 decay toward zero.*

---

## 3. The nominations

### The novelty-versus-actionability tradeoff, stated plainly

This is the central finding and it is not hidden. The strongest, most-reproducible **genuinely-new** hits are almost all **intracellular scaffolds, chromatin factors, RNA-binding proteins and housekeeping enzymes** — real regulators, but poor systemic-drug targets (11/15 challenging-to-intractable; 13/15 carry a safety watch or liability from broad expression and/or high loss-of-function intolerance — computed over all 15 worked-up candidates, enumerated in §3 below and in `coverage_table.csv`). The **tractable, safe, genetically-supported** candidates are exactly the ones that are **recovered-known** — cell-surface receptors and an E3 ligase already in the clinic. A genuinely-new *and* cleanly-druggable target is rare in this screen; the honest yield is small, and it is reported as such.

### Th1-skewing direction — inhibitor points toward allergic/atopic disease

| Gene | New vs known | Case | Druggability | Genetics | Safety |
|---|---|---|---|---|---|
| **JADE2** | **genuinely-new** | Best genuinely-new Th1 lead; perfect two-guide reproducibility (guide-r 1.00); PHD-finger chromatin reader in a histone-acetyltransferase complex | challenging — small-molecule (PHD/HAT) or degrader | weak (MS rs2084007, p=2e-13) | watch — broad incl. brain, pLI 1.0 |
| TNFRSF9 (4-1BB) | recovered-known | Costimulatory receptor; axis re-finds it in the KD-lowers-Th1 direction | tractable — cell-surface, agonist/antagonist Ab (urelumab, utomilumab; NCT02658981) | strong (RA p=3e-9, SLE) | clean/immune-restricted (resting-state caveat) |
| RAC2 | recovered-known | Hematopoietic Rho GTPase, established immune biology | challenging — small GTPase | weak, protective-on-LoF (RA) | watch — pLI 1.0; LoF → immunodeficiency |
| INTS6, URM1, BAHD1 | genuinely-new | Real regulators but essential-complex (INTS6=Integrator) or undruggable class (URM1 sulfur carrier; BAHD1 heterochromatin scaffold) | intractable / challenging | none–weak | liability / watch |

### Th2-skewing direction — inhibitor points toward autoimmune (MS/IBD/RA)

| Gene | New vs known | Case | Druggability | Genetics | Safety |
|---|---|---|---|---|---|
| **TATDN2** | **genuinely-new** | **Best genuinely-new lead in the whole screen**; no prior target literature, no ChEMBL target, no trials; highest T-cell specificity of its direction (K562 concordance 0.26); a TatD-family metal-dependent nuclease with a defined active site | challenging — small molecule (metal-dependent active site) | **none** (no immune GWAS) | watch — broad expression, LoF-tolerant |
| IFNAR1 | recovered-known | Type-I IFN receptor; axis re-finds it in the KD-raises-Th2 direction, and its approved antagonist (anifrolumab, NCT03435601) matches that prediction | tractable — cell-surface antagonist biologic | strong, protective-on-LoF (MS OT 0.61, SLE) | watch — ubiquitous; on-target = impaired antiviral immunity |
| TRIT1 | genuinely-new | tRNA isopentenyltransferase with a substrate pocket | challenging | weak (asthma, WBC count) | watch — recessive mito disease liability |
| TWF1, SNRPB2, PAXIP1 | new / partial | Undruggable class (TWF1 actin regulator) or essential complex (SNRPB2 core spliceosome; PAXIP1 MLL3/4 + DDR scaffold) | intractable | weak | watch / liability |

### Activation axis (exploratory — support only)

- **CBLB** — recovered-known canonical T-cell brake, active immuno-oncology target (Cbl-b inhibitors NCT05107674; 566 ChEMBL ligands), tractable E3 ligase. Genetics **risk-on-LoF** (hypothyroidism p=3e-41): loss of function *raises* autoimmune risk — consistent with its brake-release (anti-tumor) direction being desirable but an autoimmune liability. Confirmation of known biology on a failed axis, not a nomination.
- **CD6** — recovered-known (itolizumab Ph3), immune-restricted, strong protective-on-LoF genetics (MS rs4939490 p=2e-29, IBD). A clean known target the exploratory axis re-finds.
- **CLCC1** — genuinely-new chloride-channel-like protein; passes the measured-specificity filter but is the weakest of the measured set (K562 concordance 0.48, highest/least-specific of the 8 measured candidates), with no genetic support, challenging druggability (ion-channel-like membrane protein), and a safety watch on broad expression. Exploratory-axis hit only; not carried forward.

The 15 worked-up candidates are: Th1-skewing (6) TNFRSF9, RAC2, JADE2, BAHD1, INTS6, URM1; Th2-skewing (6) IFNAR1, TATDN2, TRIT1, TWF1, SNRPB2, PAXIP1; activation-exploratory (3) CBLB, CD6, CLCC1. The aggregate figures below are computed over all 15 (tractability: 4 tractable, 7 challenging, 4 intractable → 11 challenging-to-intractable; safety: 2 clean, 10 watch, 3 liability → 13 watch-or-liability).

---

## 4. Recommendation

1. **Primary genuinely-new lead: TATDN2** (Th2-skewing → inhibit for a Th1 shift; autoimmune-relevant). The only candidate that is simultaneously target-novel, most T-cell-specific in its direction, and has a real enzymatic active site. Caveats stated openly: **zero human-genetic support so far** and broad tissue expression.
2. **Second, harder: JADE2** (Th1-skewing; chromatin reader, best reproducibility, degrader modality) as a tool-compound / lead-generation project, contingent on ruling out broad essentiality.
3. **Recovered-known set (TNFRSF9, IFNAR1, CD6, CBLB, RAC2) is not a discovery** — it is the axis's proof-of-correctness and is labelled as such throughout.
4. **Weaker-evidence direction, stated openly:** the entire **activation axis failed validation** (ρ 0.13) and everything drawn from it — including the biologically-compelling CBLB — is exploratory support, not a nomination. Rebuilding that axis from the authors' activation model (rather than a bulk CPM contrast) is the prerequisite for trusting its hits.

Nothing here is called druggable or novel beyond what the connector evidence shows. Every gap is recorded in the coverage table.

---

## 5. Falsifiable experiment per lead

**TATDN2 (primary new lead).** *Claim:* CRISPRi knockdown of TATDN2 in polarizing human CD4+ T cells shifts them toward Th2, so a TATDN2 inhibitor would push toward Th1 and be useful in Th2-driven atopic disease. *Competing explanation not yet ruled out:* the shift is an indirect proliferation/fitness artifact of a broadly-expressed gene. *Falsifier:* knock down TATDN2 (2 guides) in Th1- and Th2-skewing cultures; measure IFN-γ⁺ vs IL-4/IL-13⁺ fractions by flow at 96 h with a proliferation dye in parallel. The claim predicts a **≥1.5-fold rise in the IL-4⁺/IFN-γ⁺ ratio** under knockdown; if the ratio is unchanged (or moves the other way) while the cells simply divided less, the polarization claim is wrong.

**JADE2 (second lead).** *Claim:* JADE2 knockdown skews CD4+ T cells toward Th1, so its inhibition/degradation would be useful where a Th1 shift is undesired (atopy). *Falsifier:* same polarization-flow readout; the claim predicts a **reproducible drop in the IL-4⁺/IFN-γ⁺ ratio** with both guides. Because JADE2 is a broadly-essential chromatin factor (pLI 1.0), the discriminating control is a parallel viability/cell-cycle stain — if the polarization shift tracks a general loss of viable dividing cells rather than a specific cytokine re-balance, the target-specific claim fails.

**CBLB (exploratory, brake-release direction).** *Claim:* CBLB knockdown releases the T-cell activation brake (anti-tumor direction). *Caveat:* the activation axis it comes from failed validation, and CBLB genetics are risk-on-LoF (autoimmune liability). *Falsifier:* CBLB knockdown should **raise CD25/CD69 and IL-2 secretion** upon suboptimal TCR stimulation; if activation markers are unchanged, the brake-release claim fails — but note this tests known biology, not a novel nomination.

---

## Artifacts
- `summary_report.md` — this report.
- `candidate_workup_report.md` — detailed per-candidate workup with all identifiers.
- `coverage_table.csv` — 15 candidates × 34 evidence columns (novelty, druggability, genetics, safety) with every gap.
- `validation_scatter.png` (Fig 1), `fig2_ranking.png` (Fig 2), `fig3_leads.png` (Fig 3).
- `novelty_read.csv`, `druggability_read.csv`, `genetics_read.csv`, `safety_read.csv` — the four specialist reads.
- `ranked_polarization_Stim8hr_full.csv`, `ranked_activation_Stim8hr_full.csv` — full ranked tables.
- Supporting data: `polar_score_by_condition.csv`, `gtex_focus_tpm.csv`, `hpa_summary.csv`.
