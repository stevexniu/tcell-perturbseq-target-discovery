#!/usr/bin/env python
"""Dual-axis discovery driver for Modal. Reuses the validated building blocks in discover.py,
but (1) builds the polarization axis from Ota 2021 ONLY (the skill's contains('Th2_vs_Th1')
filter contaminates it with Hollbacker 2021 via dict-overwrite), and (2) dumps the FULL per-gene
ranked table and the FULL validation merge for BOTH axes at Stim8hr, not just the top-40 JSON.
Outputs land under ./out/ for harvest."""
import os, json
import numpy as np, pandas as pd, h5py
from scipy.stats import spearmanr

import discover as D  # validated building blocks: read_col, fetch, build_activation_signature, project, load_authors_reg, load_k562

OUT = "out"; os.makedirs(OUT, exist_ok=True)
DATA_DIR = "/tmp/tcell_data"
CONDITION = "Stim8hr"

POLAR_SIG_URL = ("https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/"
                 "metadata/suppl_tables/Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv")

def build_polarization_signature_ota(genes):
    """Authors' published Th1/Th2 polarization signature, Ota 2021 contrast ONLY (Th2-up = +).
    Matches the 'ota' validation coefficient column; avoids Hollbacker contamination."""
    print(f"polarization signature (Ota 2021 only) <- {POLAR_SIG_URL}", flush=True)
    t = pd.read_csv(POLAR_SIG_URL)
    t = t[t["contrast"].astype(str) == "Th2_vs_Th1 (Ota 2021)"]
    print(f"  Ota rows: {len(t)}", flush=True)
    wmap = dict(zip(t["variable"].astype(str), t["zscore"].astype(float)))
    return np.array([wmap.get(g, 0.0) for g in genes], np.float32)

def run_axis(axis, genes, o_arrays, Z, k562_df, areg):
    """Score, validate, rank both directions, flag hub + cell-line-shared. Returns (full_df, validation, thresholds)."""
    P, gsig, gall, ndown, ncells, nguides, single, ontgt, offd, neigh = o_arrays

    if axis == "activation":
        pb = D.fetch(DATA_DIR, D.PB_FILE)
        w = D.build_activation_signature(pb, genes)
    else:
        w = build_polarization_signature_ota(genes)
    s = D.project(Z, w)  # per-observation score (this condition only)

    # per-gene mean score across guides/obs of this (gene, condition)
    per_gene = pd.DataFrame({"gene": P, "score": s}).groupby("gene", as_index=False)["score"].mean()

    # ---- VALIDATION vs authors' own coefficients ----
    val_sig = D.VALIDATION_SIG[axis]
    act = areg[(areg["signature"] == val_sig) & (areg["celltype"] == CONDITION)][
        ["regulator", "coef_mean", "known_regulators"]].dropna(subset=["coef_mean"])
    merged = per_gene.merge(act, left_on="gene", right_on="regulator", how="inner")
    rho = spearmanr(merged["score"], merged["coef_mean"])[0] if len(merged) > 20 else None
    known = set(act.loc[act["known_regulators"].astype(str).str.lower().isin(["true","1","1.0"]), "regulator"])
    act = act.assign(abscoef=act["coef_mean"].abs())
    authors_top = set(act.sort_values("abscoef", ascending=False).head(300)["regulator"])
    # dump the FULL validation merge (for AUROC + scatter locally)
    merged["is_known"] = merged["regulator"].isin(known)
    merged.to_csv(f"{OUT}/validation_merge_{axis}.csv", index=False)
    validation = {"axis": axis, "condition": CONDITION, "validation_signature": val_sig,
                  "authors_coef_spearman": None if rho is None else round(float(rho),4),
                  "n_matched": int(len(merged)), "n_authors_known_regulators": int(len(known)),
                  "n_authors_top300": int(len(authors_top))}
    print("VALIDATION:", validation, flush=True)

    # ---- per-gene table over ALL perturbations (QC flags kept as columns, nothing dropped) ----
    df = pd.DataFrame({"gene": P, "score": s, "guide_signif": gsig, "guide_all": gall,
                       "n_downstream": ndown, "n_cells_target": ncells, "n_guides": nguides,
                       "single_guide": single, "ontarget_signif": ontgt,
                       "distal_offtarget": offd, "neighboring_KD": neigh})
    # aggregate to one row per gene (median for QC metrics, any() for on-target)
    agg = df.groupby("gene", as_index=False).agg(
        score=("score","mean"), guide_signif=("guide_signif","median"),
        guide_all=("guide_all","median"), n_downstream=("n_downstream","median"),
        n_cells_target=("n_cells_target","median"), n_guides=("n_guides","median"),
        single_guide=("single_guide","max"), ontarget_signif=("ontarget_signif","max"),
        distal_offtarget=("distal_offtarget","max"), neighboring_KD=("neighboring_KD","max"))

    # QC-pass gate (same as skill): on-target worked, guides agree, real effect, no cis/distal off-target
    agg["qc_pass"] = (agg["ontarget_signif"].astype(bool) & (agg["guide_signif"] >= 0.5)
                      & (agg["n_downstream"] >= 20) & ~agg["distal_offtarget"].astype(bool)
                      & ~agg["neighboring_KD"].astype(bool)
                      & ~agg["gene"].isin(["NTC","None","non-targeting"]))

    # status vs authors' own table
    def status(g):
        if g in known: return "authors-known"
        if g in authors_top: return "authors-nominated"
        return "novel"
    agg["status"] = agg["gene"].map(status)

    # hub / broad-effect: top-decile downstream count among QC-pass genes
    hub_thr = float(agg.loc[agg["qc_pass"], "n_downstream"].quantile(0.9))
    agg["hub_flag"] = np.where(agg["n_downstream"] >= hub_thr, "hub/broad-effect", "")

    # cell-line-shared / generic: high concordance with same KD in K562 (authors' comparison)
    k = k562_df[k562_df["condition"] == CONDITION].dropna(subset=["logfc_pearson_r"])
    kr = k.groupby("target_contrast_gene_name")["logfc_pearson_r"].mean()
    agg["k562_concordance"] = agg["gene"].map(kr)
    qc_k = agg.loc[agg["qc_pass"], "k562_concordance"].dropna()
    k_thr = float(qc_k.quantile(0.67)) if len(qc_k) else None
    agg["specificity_flag"] = np.where(
        agg["k562_concordance"].notna() & (agg["k562_concordance"] >= (k_thr if k_thr is not None else np.inf)),
        "cell-line-shared/generic", "")
    agg["k562_concordance"] = agg["k562_concordance"].round(4)

    agg["axis"] = axis; agg["condition"] = CONDITION
    agg = agg.sort_values("score").reset_index(drop=True)
    agg.to_csv(f"{OUT}/per_gene_{axis}.csv", index=False)
    thresholds = {"hub_threshold_n_downstream": hub_thr, "k562_concordance_threshold": k_thr}
    validation.update(thresholds)
    return agg, validation

def main():
    de = D.fetch(DATA_DIR, D.DE_FILE)
    with h5py.File(de, "r") as f:
        genes = np.array([D._dec(x) for x in f["var"]["gene_name"][:]], dtype=object)
        o = f["obs"]
        pert  = D.read_col(o, "target_contrast_gene_name").astype(str)
        cond  = D.read_col(o, "culture_condition").astype(str)
        gsig  = D.read_col(o, "guide_correlation_signif").astype(float)
        gall  = D.read_col(o, "guide_correlation_all").astype(float)
        ndown = D.read_col(o, "n_downstream").astype(float)
        ncells= D.read_col(o, "n_cells_target").astype(float)
        nguides=D.read_col(o, "n_guides").astype(float)
        tobool= lambda col: np.array([str(x).lower() in ("true","1","1.0") for x in D.read_col(o, col)])
        single= tobool("single_guide_estimate")
        ontgt = tobool("ontarget_significant")
        offd  = tobool("distal_offtarget_flag")
        neigh = tobool("neighboring_gene_KD")
        sel = cond == CONDITION
        Z = f["layers"]["zscore"][:][sel].astype(np.float32)
    o_arrays = tuple(a[sel] for a in (pert, gsig, gall, ndown, ncells, nguides, single, ontgt, offd, neigh))

    areg = D.load_authors_reg()
    k562 = D.load_k562()

    all_val = []
    for axis in ["activation", "polarization"]:
        _, val = run_axis(axis, genes, o_arrays, Z, k562, areg)
        all_val.append(val)
    json.dump({"condition": CONDITION, "axes": all_val}, open(f"{OUT}/validation_summary.json","w"), indent=2, default=str)
    print("DONE. wrote per_gene_{activation,polarization}.csv, validation_merge_*.csv, validation_summary.json", flush=True)

if __name__ == "__main__":
    main()
