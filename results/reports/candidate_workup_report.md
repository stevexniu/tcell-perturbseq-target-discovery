# Evidence-ranked target shortlist — CD4+ T-cell polarization regulators

**Source:** 15 candidates selected from the dual-axis Perturb-seq discovery with a documented cutoff — validated **polarization axis** (Spearman 0.72 vs authors' coefficients) weighted over the **activation axis** (failed validation, ρ 0.13, exploratory support only). Selection rule for polarization: QC-pass, non-hub, not cell-line-shared, guide-reproducibility ≥ 0.70, and absent from the authors' own regulator table (`novel`); top 6 by |score| per direction. Activation-axis genes (CBLB, CD6, CLCC1) are an exploratory judgment subset, not from a single threshold.

Each candidate was worked up on four independent evidence dimensions via live database connectors — **novelty** (OpenAlex/PubMed, ChEMBL, Drugs@FDA, ClinicalTrials.gov), **druggability** (UniProt/InterPro, ChEMBL, PDB/AlphaFold), **human genetics** (GWAS Catalog, eQTL Catalogue, Open Targets), and **safety** (GTEx/HPA expression breadth, gnomAD constraint). All identifiers below come from those queries; gaps are marked, not filled.

![Evidence coverage across candidates × dimensions]({{artifact:f239610f-d120-44be-b14a-50abbb0dde04}})

---

## Headline finding — the novelty/druggability tension is the result

The strongest-scoring, most-reproducible **genuinely-new** hits are almost all **intracellular scaffolds, chromatin factors, RNA-binding proteins and housekeeping enzymes** — real regulators, but poor systemic-drug targets: 11 of 15 are challenging/intractable, and **13 of 15 carry a safety watch/liability** (broad expression and/or high LoF-intolerance). The **tractable, safe, genetically-supported** candidates are exactly the ones that are **recovered-known** (cell-surface receptors TNFRSF9, CD6, IFNAR1; E3 ligase CBLB) — all with drugs already in the clinic. That the validated axis re-finds these known druggable immune targets *in the correct direction* is the strongest confirmation it works; it also means the honest novel-and-druggable yield is small, and I flag it rather than inflate it.

**Gap census across the 15:** no genetic support 4/15 · specificity untested (no K562 data) 7/15 · safety watch-or-liability 13/15 · challenging/intractable druggability 11/15.

---

## Th1-skewing direction (knockdown → more Th1 / less Th2)

**Therapeutic area this points to:** enhancing Th1 / dampening Th2 → **allergic and atopic disease** (asthma, atopic dermatitis) by an inhibitor of the target; the Th1-driving direction is also broadly relevant to anti-tumor immunity.

| Rank | Gene | Novelty | Druggability | Genetics | Safety | Verdict |
|---|---|---|---|---|---|---|
| 1 | **TNFRSF9** (4-1BB/CD137) | **recovered-known** — agonist Abs urelumab (Ph2), utomilumab (Ph3); ChEMBL CHEMBL3712857; NCT00351325, NCT02658981 | tractable — cell-surface TNFR, biologic (agonist/antagonist Ab), PDB available | strong — RA (p=3e-9), autoimmune thyroid, SLE (Open Targets) | clean (resting-state); immune-restricted, LoF-tolerant (pLI 0.003) — caveat: activation-inducible, GTEx resting-only | Known IO target; recovered, not new |
| 2 | **RAC2** | **recovered-known** — hematopoietic Rho GTPase, established immune-cell biology; ChEMBL CHEMBL5581 | challenging — small GTPase, no easy pocket | weak, but **protective-on-LoF** for RA; celiac (rs9610686) | **watch** — pLI 1.0, LOEUF 0.20; germline LoF → immunodeficiency/neutrophil dysfunction | Known biology; essentiality liability |
| 3 | **INTS6** | genuinely-new (as a *target*) | challenging → **no obvious modality** (Integrator scaffold) | none (immune) | **liability** — core Integrator complex, broadly essential, pLI 1.0 | New but essential-complex; drop |
| 4 | **JADE2** | **genuinely-new** | challenging — PHD/HAT E3-like; degrader speculative | weak — MS (rs2084007, p=2e-13), WBC counts | watch — broad incl. brain, pLI 1.0 LOEUF 0.24 | New; best-reproduced (guide-r 1.00) but chromatin-essential |
| 5 | **BAHD1** | **genuinely-new** | challenging — heterochromatin reader/scaffold, degrader-only speculative | weak — UC; TSH levels (rs1659563) | watch — broad, pLI 0.998 | New; hard modality, broad expression |
| 6 | **URM1** | **genuinely-new** | **intractable** — ubiquitin-like sulfur carrier, no pocket | none | watch — housekeeping tRNA-thiolation | New but essentially undruggable |

**Th1 read:** the direction is validated and its top of list is a bona-fide druggable IO target (TNFRSF9), recovered not discovered. Among genuinely-new hits, **JADE2** is the most interesting lead-generation prospect (perfect guide reproducibility, chromatin-reader with degrader potential) but carries a broad-essentiality caveat; INTS6/URM1 should be dropped as essential/undruggable.

---

## Th2-skewing direction (knockdown → more Th2 / less Th1)

**Therapeutic area this points to:** enhancing Th1 by *inhibiting* one of these Th2-promoting genes → **autoimmune / Th1-Th17 disease** (MS, IBD, RA); conversely the Th2-promoting direction itself is relevant where a Th2 shift is desired.

| Rank | Gene | Novelty | Druggability | Genetics | Safety | Verdict |
|---|---|---|---|---|---|---|
| 1 | **IFNAR1** | **recovered-known** — antagonist Ab anifrolumab (Ph4, approved SLE); ChEMBL CHEMBL1887; NCT03435601 | tractable — cell-surface type-I IFN receptor, antagonist biologic | **strong** — MS (Open Targets 0.61), SLE 0.58, **protective-on-LoF** | watch — ubiquitous surface receptor; on-target = impaired antiviral immunity (anifrolumab precedent) | Known target; direction matches drug (antagonist) |
| 2 | **TATDN2** | **genuinely-new** — no target pubs, no ChEMBL target, no trials | challenging — TatD metal-dependent nuclease, a *defined active site* (better than a scaffold) | none | watch — broad, LoF-tolerant (pLI 0), no strong essentiality | **Most genuinely-new + most T-cell-specific (K562 0.26); best new lead** |
| 3 | **TRIT1** | **genuinely-new** | challenging — tRNA isopentenyltransferase, substrate pocket | weak — asthma, WBC count (quantitative) | watch — cyto+mito enzyme; **recessive mito disease** liability | New; enzyme pocket but housekeeping/mito risk |
| 4 | **TWF1** | genuinely-new | **intractable** — twinfilin actin regulator, no modality | weak — IBD (suggestive p=6e-6) | watch — ubiquitous cytoskeletal, pLI 0.998 | New; undruggable class |
| 5 | **SNRPB2** | genuinely-new | **intractable** — core U2 snRNP | weak — food-allergy (sub-genome-wide) | **liability** — core spliceosome, broadly essential | New but essential; drop |
| 6 | **PAXIP1** | partial | **intractable** — BRCT scaffold, MLL3/4 + DDR | weak — colitis (suggestive) | **liability** — pLI 1.0 LOEUF 0.20, broad incl. brain | Essential chromatin/genome-maintenance; drop |

**Th2 read:** IFNAR1 anchors the direction as a recovered-known target whose approved antagonist (anifrolumab) points the way our knockdown predicts — strong axis validation. The standout **genuinely-new** candidate is **TATDN2**: no prior target literature, highest T-cell specificity of the direction (K562 concordance 0.26), a metal-dependent nuclease active site (a real, if challenging, small-molecule opportunity), and no strong essentiality signal. It is the single best new lead in the whole set, with the honest caveats that it has *no* human-genetic support yet and broad expression.

---

## Activation axis (exploratory — axis failed validation, support only)

- **CBLB** — **recovered-known**: canonical T-cell brake, active IO target (Cbl-b inhibitors NCT05107674, NCT05662397; ChEMBL CHEMBL4879459, **566 ligands**), tractable E3 ligase. Genetics **strong but risk-on-LoF** (hypothyroidism p=3e-41, rs1921445) — loss of function *raises* autoimmune risk, exactly matching its brake-release (immune-enhancing/anti-tumor) direction and warning against it for autoimmune indications. Recovered, not new; the axis's failure to validate means this is confirmation-of-known, not a nomination.
- **CD6** — **recovered-known**: itolizumab (Ph3), Oncolysin-CD6 (Ph2); cell-surface, immune-restricted, strong genetics (MS rs4939490 p=2e-29, IBD, **protective-on-LoF**). A clean known target the exploratory axis re-finds.
- **CLCC1** — genuinely-new but ER chloride channel, no immune genetics, recessive-disease (retinitis pigmentosa) liability; not pursued.

---

## Recommendation

1. **If a genuinely-new lead is the goal:** **TATDN2** (Th2-skewing → inhibit for a Th1 shift, autoimmune-relevant) is the only genuinely-new candidate that is simultaneously target-novel, most T-cell-specific in its direction, and has a real (if challenging) enzymatic pocket. Prioritize it — with eyes open about zero current genetic support and broad expression.
2. **Second tier, new but harder: JADE2** (Th1-skewing, chromatin reader, best reproducibility, degrader modality) as a lead-generation / tool-compound project, contingent on ruling out broad essentiality.
3. **Everything else genuinely-new is either essential-complex (INTS6, SNRPB2, PAXIP1 — drop) or class-undruggable (URM1, TWF1).**
4. **The recovered-known set (TNFRSF9, IFNAR1, CD6, CBLB, RAC2) is not a discovery** but is the axis's proof-of-correctness: it re-finds clinically-drugged immune targets in the right direction. Do not present them as novel.

Nothing here is called druggable or novel beyond what the connector evidence shows; the coverage table records every gap.

---

## One falsifying experiment per lead

**TATDN2 (primary new lead).** *Claim:* CRISPRi knockdown of TATDN2 in polarizing human CD4+ T cells shifts them toward Th2 (raises the Th2-vs-Th1 axis), so a TATDN2 inhibitor would push toward Th1 and be useful in Th2-driven atopic disease. *Bet:* TATDN2's nuclease activity supports a transcriptional program permissive for Th1; competing explanation not yet ruled out — the effect is an indirect fitness/proliferation artifact of a broadly-expressed gene rather than a specific polarization role. *Falsifier:* knock down TATDN2 (2 guides) in Th1- and Th2-skewing cultures and measure IFN-γ⁺ vs IL-4/IL-13⁺ fractions by flow at 96 h. The claim predicts a **≥1.5-fold rise in IL-4⁺/IFN-γ⁺ ratio** under knockdown; if the ratio is unchanged (or moves the other way) while a proliferation dye shows the cells simply divided less, the polarization claim is wrong.

**JADE2 (second-tier new lead).** *Claim:* JADE2 knockdown skews CD4+ T cells toward Th1. *Falsifier:* same polarization-flow readout; the claim predicts a **reproducible drop in the IL-4⁺/IFN-γ⁺ ratio** with both guides. Because JADE2 is a broadly-essential chromatin factor (pLI 1.0), the discriminating control is a viability/cell-cycle stain run in parallel — if the polarization shift tracks a general loss of viable dividing cells rather than a specific cytokine re-balance, the target-specific claim fails.

## Artifacts
- `coverage_table.csv` — all 15 candidates × 34 evidence columns (novelty, druggability, genetics, safety) with every identifier and gap.
- `novelty_read.csv`, `druggability_read.csv`, `genetics_read.csv`, `safety_read.csv` — the four specialist reads.
- `coverage_heatmap.png` — the coverage figure above.
- `candidates_for_workup.csv` — the 15-candidate set with documented selection basis.
