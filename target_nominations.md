# T-cell Target Nominations — CD4+ Perturb-seq (Stim8hr)

Two nominations, one per direction. Candidates were drawn from the reproducible, non-hub tail of
each shortlist (hub/fitness genes set aside), then investigated on three axes in parallel:
novelty (prior target/drugging/author-nomination), druggability (class/modality/ligands/structure),
and human genetics (autoimmune GWAS/eQTL, constraint, direction-of-effect for an inhibitor).
Every claim is grounded in a retrieved source id. Where an axis returned no support, it says so.

## Evidence matrix (all six investigated)

| gene | dir | score | guide-repro | downstream | novelty | druggability | genetics |
|---|---|---|---|---|---|---|---|
| CLCC1 | auto | -0.968 | 0.798 | 679 | genuinely-new | hard | neutral |
| MAP3K1 | auto | -0.867 | 0.672 | 89 | recovered-known-broad | tractable-SM | supports (LoF-intolerant) |
| WAS | auto | -0.965 | 0.632 | 494 | recovered-known-broad | hard | contradicts |
| RSBN1L | anti | +0.554 | 1.000 | 178 | genuinely-new | hard | supports |
| MBTPS1 | anti | +0.451 | 0.923 | 279 | recovered-known-broad | tractable-SM | neutral |
| PUM1 | anti | +0.524 | 0.999 | 226 | partially-characterized | hard | neutral (haploinsufficiency flag) |

Key finding: **novelty and actionability are anti-correlated here.** The genuinely-new hits
(RSBN1L, CLCC1) have no chemical matter yet; the tractable hits (MAP3K1, MBTPS1) are already
known targets. So I nominate the best-*evidenced* lead per direction and mark its status plainly,
rather than forcing a novel-but-unactionable gene into the lead slot.

---

# NOMINATION 1 — Anti-tumor: RSBN1L  (genuinely new)

**The case.** RSBN1L was the single most reproducible brake in the screen: knockdown raised the
effector program with guide-reproducibility 1.000 and a specific footprint (n_downstream 178, well
below the hub threshold), so the effect is unlikely to be general fitness loss. It is genuinely new
to immunity — only 3 PubMed records exist for the gene, none in a T-cell or drug-target context,
and it is not named in the dataset paper's abstract. Human genetics point the right way for an
inhibitor: gene-body variants associate with autoimmune hypothyroidism (p=2e-24), eosinophilia
(p=2e-13), SLE and allergy — a heightened-immunity / brake-loss signature — and the gene is
loss-of-function-tolerant (gnomAD LOEUF 0.60), the profile an inhibitor wants. Mechanistically it
is annotated as an RSBN1-family putative H4K20 histone demethylase (a 2-OG-dioxygenase-like
enzyme class), giving a plausible transcriptional-brake mechanism and an attractive druggable
enzyme family.

**Recovered-known vs new.** Genuinely new — no prior target nomination, no drugs, not
author-highlighted (abstract only; full text was 403-blocked). The recovered-known anchor in this
direction is MBTPS1, whose inhibition-enhances-antitumor-immunity direction is already published
(PMID 40307212) — reassuring that the axis finds real brakes, but not novel.

**The catch, plainly.** Two gaps. (1) Druggability is *hard today*: the enzyme class is attractive
but there are zero ligands and no experimental structure — a lead would need a de novo chemistry
campaign (AlphaFold model only). (2) The demethylase mechanism is annotation, not demonstrated in
T cells; the eQTL colocalization to RSBN1L expression was not formally run, so 'genetics supports'
rests on trait association + variant position + constraint, not proven causal direction. (3) The
authors' own activation regulator coefficient for RSBN1L is −0.022 (z ≈ −39, more negative than 94%
of genes), which places it at the *driver* end — discordant with our brake call (+0.554). The two
estimators measure different things (signature projection of the KD profile vs regularized
multivariate regression), and the direct activation-shift assay (Rest→Stim) supports the brake
direction, but this method-level conflict is the primary reason to confirm the brake effect in an
independent arrayed assay before investing further.

**Safety note (therapeutic claim made).** Favorable on constraint: LoF-tolerant (LOEUF 0.60), no
Mendelian phenotype, so partial inhibition is unlikely to hit an essential function — but its
non-immune expression breadth is uncharacterized and would need a tissue-expression/constraint
check before committing.

**One experiment that would falsify it.** Arrayed CRISPRi knockdown (and catalytic-dead rescue) of
RSBN1L in primary human CD4+ and CD8+ T cells, then TCR stimulation + a tumor-cell co-culture
killing assay. The hypothesis predicts KD *increases* effector cytokines (IL-2/IFNG) and tumor
killing, and that a catalytic-dead demethylase mutant fails to rescue (enzyme-dependent brake). If
KD does not raise effector output, or if the effect tracks only with a viability/fitness change, or
if catalytic-dead rescues fully, the brake-via-demethylase hypothesis is falsified.

---

# NOMINATION 2 — Autoimmune: MAP3K1 (MEKK1)  (recovered-known-broad)

**The case.** Among the autoimmune-direction candidates, MAP3K1 is the only one where all three
axes align. Knockdown lowered activation with a highly specific footprint (n_downstream 89 — the
tightest among the six nominees; the shortlist-wide minimum is 20), consistent with a focused
signaling regulator rather than a hub.
It is a documented TCR/JNK-pathway kinase (T-cell-specific deletion alters effector responses,
PMID 26774476), it is druggable now (Ser/Thr kinase; 8 nM ChEMBL ligand, PDB 6WHB,
a Phase-2 clinical inhibitor E-6201 / CHEMBL1097999, Open Targets 'SM Advanced Clinical'), and
human genetics support an immunosuppressive target: MAP3K1 is a genuine inflammatory-disease locus.
Genome-wide (verified against the live GWAS Catalog): psoriasis, rs12654176, p=9e-12, PMID 40021644,
and polyarteritis nodosa (PheCode 446), rs535102440, p=2e-12, PMID 39024449 — both mapped to
MAP3K1/C5orf67, the only two immune-mediated traits mapped to the gene body at genome-wide
significance. Additional autoimmune signals (Crohn's
p=6.6e-6, Graves' ophthalmopathy p=5.0e-6, palindromic rheumatism p=1.8e-6) come from FinnGen R12
PheWAS at *sub*-genome-wide significance, so they are suggestive rather than established (an earlier
draft presented these as equivalent to the psoriasis locus, which was overstated).

**Recovered-known vs new.** Recovered-known-broad — a known druggable kinase already pursued as an
inhibitor. It is *not* a novel target. The genuinely-new autoimmune hit, CLCC1, is an ER chloride
channel with zero T-cell literature — but it fails the actionability bar (no small-molecule pocket,
genetics neutral, broadly-essential ER housekeeping channel), so it is carried as a
new-but-not-yet-actionable alternative, not the lead.

**The catch, plainly.** (1) Selectivity: the best clinical molecule (E-6201) is a dual MEK/MEKK1
inhibitor, so a MEKK1-selective chemical would be needed to test the specific hypothesis. (2)
Constraint: MAP3K1 is LoF-intolerant (LOEUF 0.33) with a Mendelian sex-reversal knockout phenotype,
so only *partial* inhibition is viable — full blockade is not a safe therapeutic model. (3) Novelty
is low by design here; the value is a well-supported, fast-to-test lead, not a first-in-class idea.

**Safety note (therapeutic claim made).** LoF-intolerant (LOEUF 0.33) + Mendelian developmental
knockout phenotype means the therapeutic window depends on partial/tunable inhibition; dose-limiting
on-target toxicity is the main risk. This is the reason to keep CLCC1 and the novel tail in view
despite weaker current evidence.

**One experiment that would falsify it.** Titrated selective MEKK1 inhibition (or a degron) in
primary human CD4+ T cells across a dose range, with viability controls, then measure activation
(CD25/CD69, IL-2/IFNG) vs a general-toxicity readout. The hypothesis predicts a dose window where
activation drops *before* viability does (specific immunosuppression). If reduced activation only
appears at doses that also reduce viability, the 'specific activation regulator' claim collapses
into a fitness effect. Separately, the genetics rationale is falsified if the autoimmune
risk alleles are shown to *increase* rather than decrease MAP3K1 activity.

---

## Coverage & gaps (both nominations)

| axis | RSBN1L (anti-tumor) | MAP3K1 (autoimmune) |
|---|---|---|
| screen reproducibility | covered — repro 1.000, specific (nd 178) | covered — repro 0.672, most-specific (nd 89) |
| novelty | covered — genuinely new | covered — recovered-known-broad |
| prior drugging | covered — none | covered — E-6201 Phase 2 (dual MEK/MEKK1) |
| protein class / structure | partial — enzyme class known, AlphaFold only, no PDB | covered — kinase, PDB 6WHB |
| chemical matter | GAP — zero ligands | covered — 8 nM ligand |
| autoimmune GWAS | covered — hypothyroid/eosinophil/SLE | covered — psoriasis (p=9e-12) + polyarteritis nodosa (p=2e-12), both GWAS; Crohn/Graves FinnGen-suggestive only |
| eQTL colocalization | GAP — not formally run | GAP — not formally run |
| direction-of-effect | supports (brake-loss, LoF-tolerant) | supports (inflammatory locus) |
| constraint/safety | favorable (LOEUF 0.60) | flag (LOEUF 0.33 + Mendelian) |
| mechanism in T cells | GAP — demethylase role unproven | covered — TCR/JNK kinase |

Shared gap: neither has a formal allele-direction eQTL colocalization; both 'genetics supports'
calls rest on trait association + variant position + constraint. That coloc is the highest-value
next genetics step for either.
