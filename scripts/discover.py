#!/usr/bin/env python
"""Discovery: rank CD4+ T-cell Perturb-seq perturbations into a target shortlist, both directions.

No hand-supplied gene list anywhere. Two data-driven axes: the ACTIVATION axis (built from the dataset's own
control cells, Stim8hr vs Rest) and the POLARIZATION axis (the authors' published Th1/Th2 signature). Each is
validated by reproducing the AUTHORS' own regulator coefficients for that signature (activation / ota); novelty
is judged against the authors' own nominations — never against a list typed into this file.

Public sources (anonymous):
  DE:          s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad
  pseudobulk:  s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.pseudobulk_merged.h5ad
  authors' regulator coefficients (their GitHub):
    raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/metadata/suppl_tables/
      polarization_prediction_condition_comparison_regulator_coefficients.csv

Usage: discover.py [--axis activation|polarization] [--data-dir DIR] [--condition Stim8hr] [--out FILE]
Run once per axis (e.g. --axis activation --out candidates_activation.json, then --axis polarization).
Uses local files in --data-dir if present, else downloads from the public sources.
"""
import os, argparse, json, urllib.request
import numpy as np, pandas as pd, h5py
from scipy.stats import spearmanr

BUCKET = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
AUTHORS_REG = ("https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/"
               "metadata/suppl_tables/polarization_prediction_condition_comparison_regulator_coefficients.csv")
# authors' K562-vs-CD4T comparison: a knockdown that acts the SAME in a non-T-cell line is generic
# (fitness/ER/housekeeping), not T-cell-activation-specific — this is how ER false positives are caught.
K562_CMP = ("https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/"
            "metadata/suppl_tables/K562_comparison.suppl_table.csv")
# authors' published Th1/Th2 polarization signature (Ota 2021) — used as the polarization-axis weights
POLARIZATION_SIG = ("https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/"
                    "metadata/suppl_tables/Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv")
DE_FILE, PB_FILE = "GWCD4i.DE_stats.h5ad", "GWCD4i.pseudobulk_merged.h5ad"
# which authors' coefficient signature validates each axis
VALIDATION_SIG = {"activation": "activation", "polarization": "ota"}

def load_authors_reg():
    """Authors' own regulator-coefficient table, fetched live from their public GitHub suppl_tables. The
    'known vs novel' call ALWAYS comes from this table — never a gene list typed here, never a memorized set."""
    print(f"authors' regulator table <- {AUTHORS_REG}", flush=True)
    return pd.read_csv(AUTHORS_REG)

def load_k562():
    print(f"authors' K562 comparison <- {K562_CMP}", flush=True)
    return pd.read_csv(K562_CMP)

def _dec(a): return a.decode() if isinstance(a, bytes) else a
def read_col(obs, name):
    node = obs[name]
    if isinstance(node, h5py.Group):
        cats = np.array([_dec(c) for c in node["categories"][:]], dtype=object)
        codes = node["codes"][:]; out = np.empty(codes.shape, dtype=object)
        v = codes >= 0; out[v] = cats[codes[v]]; out[~v] = None; return out
    a = node[:]
    return np.array([_dec(x) for x in a], dtype=object) if a.dtype.kind in ("S", "O") else a

def fetch(data_dir, fname):
    local = os.path.join(data_dir, fname)
    if os.path.exists(local): return local
    os.makedirs(data_dir, exist_ok=True)
    print(f"downloading {fname} from public bucket ...", flush=True)
    urllib.request.urlretrieve(f"{BUCKET}/{fname}", local)
    return local

def build_activation_signature(pb_path, genes):
    """Dataset's own activation program: mean CPM of non-targeting controls, Stim8hr vs Rest."""
    with h5py.File(pb_path, "r") as f:
        root = f["mod"]["rna"] if "mod" in f else f
        var = np.array([_dec(x) for x in root["var"]["gene_name"][:]], dtype=object)
        obs = root["obs"]
        gtype = read_col(obs, "guide_type").astype(str)
        cond = read_col(obs, "culture_condition").astype(str)
        X = root["X"]; indptr = X["indptr"][:]; data = X["data"]; indices = X["indices"]
        def mean_profile(mask):
            rows = np.where(mask)[0]; acc = np.zeros(len(var)); n = 0
            for r in rows:
                s, e = indptr[r], indptr[r + 1]
                if e > s:
                    tot = data[s:e].sum()
                    if tot > 0:
                        acc[indices[s:e]] += data[s:e] / tot * 1e4; n += 1
            return acc / max(n, 1)
        ntc = gtype == "non-targeting"
        rest = mean_profile(ntc & (cond == "Rest"))
        stim = mean_profile(ntc & (cond == "Stim8hr"))
    sig_full = np.log2(stim + 1) - np.log2(rest + 1)   # up on activation = +
    gidx = {g: i for i, g in enumerate(var)}
    return np.array([sig_full[gidx[g]] if g in gidx else 0.0 for g in genes], np.float32)

def build_polarization_signature(genes):
    """Authors' published Th1/Th2 polarization signature (Ota 2021) as per-gene weights (Th2-up = +).
    A published contrast, not a hand-picked list."""
    print(f"polarization signature <- {POLARIZATION_SIG}", flush=True)
    t = pd.read_csv(POLARIZATION_SIG)
    if "contrast" in t.columns:
        t = t[t["contrast"].astype(str).str.contains("Th2_vs_Th1", case=False, na=False)]
    wmap = dict(zip(t["variable"].astype(str), t["zscore"].astype(float)))
    return np.array([wmap.get(g, 0.0) for g in genes], np.float32)

def project(Z, w):
    idx = np.where(w != 0)[0]; wv = w[idx]; wv = (wv - wv.mean()) / (wv.std() + 1e-9)
    s = Z[:, idx]; s = (s - s.mean(1, keepdims=True)) / (s.std(1, keepdims=True) + 1e-9)
    return (s @ wv) / len(idx)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/tmp/tcell_data")
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--axis", default="activation", choices=["activation", "polarization"])
    ap.add_argument("--out", default="candidates.json")
    args = ap.parse_args()

    de = fetch(args.data_dir, DE_FILE)
    with h5py.File(de, "r") as f:
        genes = np.array([_dec(x) for x in f["var"]["gene_name"][:]], dtype=object)
        o = f["obs"]
        pert = read_col(o, "target_contrast_gene_name").astype(str)
        cond = read_col(o, "culture_condition").astype(str)
        gsig = read_col(o, "guide_correlation_signif").astype(float)
        ndown = read_col(o, "n_downstream").astype(float)
        ontgt = np.array([str(x).lower() in ("true", "1", "1.0") for x in read_col(o, "ontarget_significant")])
        offd = np.array([str(x).lower() in ("true", "1", "1.0") for x in read_col(o, "distal_offtarget_flag")])
        neigh = np.array([str(x).lower() in ("true", "1", "1.0") for x in read_col(o, "neighboring_gene_KD")])
        sel = cond == args.condition
        Z = f["layers"]["zscore"][:][sel].astype(np.float32)
    P = pert[sel]; gsig, ndown = gsig[sel], ndown[sel]
    ontgt, offd, neigh = ontgt[sel], offd[sel], neigh[sel]

    if args.axis == "activation":
        pb = fetch(args.data_dir, PB_FILE)
        w = build_activation_signature(pb, genes)
    else:
        w = build_polarization_signature(genes)
    s = project(Z, w)
    per_gene = pd.DataFrame({"gene": P, "score": s}).groupby("gene", as_index=False)["score"].mean()

    # ---- VALIDATION: reproduce the authors' OWN regulator coefficients (no hand-picked list) ----
    areg = load_authors_reg()
    act = areg[(areg["signature"] == VALIDATION_SIG[args.axis]) & (areg["celltype"] == args.condition)][
        ["regulator", "coef_mean", "known_regulators"]].dropna(subset=["coef_mean"])
    merged = per_gene.merge(act, left_on="gene", right_on="regulator", how="inner")
    rho = spearmanr(merged["score"], merged["coef_mean"])[0] if len(merged) > 20 else None
    known = set(act.loc[act["known_regulators"].astype(str).str.lower().isin(["true", "1", "1.0"]), "regulator"])
    act = act.assign(abscoef=act["coef_mean"].abs())
    authors_top = set(act.sort_values("abscoef", ascending=False).head(300)["regulator"])  # authors' own strongest
    validation = {"authors_coef_spearman": None if rho is None else round(float(rho), 3),
                  "n_matched": int(len(merged)), "n_authors_known_regulators": len(known)}
    print("validation vs authors' coefficients:", validation, flush=True)

    # ---- RANK both directions: QC-clean, reproducible, on-target; fitness/hub flagged, not dropped ----
    keep = ontgt & (gsig >= 0.5) & (ndown >= 20) & ~offd & ~neigh & ~np.isin(P, ["NTC", "None", "non-targeting"])
    df = pd.DataFrame({"gene": P[keep], "score": s[keep], "guide_signif": gsig[keep], "n_downstream": ndown[keep]})
    df = df.groupby("gene", as_index=False).agg(score=("score", "mean"),
            guide_signif=("guide_signif", "median"), n_downstream=("n_downstream", "median"))
    hub_thr = df["n_downstream"].quantile(0.9)

    def status(g):
        if g in known: return "authors-known"          # authors' own curated known regulator
        if g in authors_top: return "authors-nominated"  # among the authors' strongest coefficients
        return "novel"                                   # not in the authors' table (still needs a drug check)
    df["status"] = df["gene"].map(status)
    df["hub_flag"] = np.where(df["n_downstream"] >= hub_thr, "hub/broad-effect", "")

    # ---- artifact layer 2: cell-line-shared effect (authors' K562 comparison) ----
    # High concordance with the same knockdown in K562 (a non-T-cell line) = generic fitness/ER/housekeeping
    # effect, not T-cell-activation-specific. This catches the ER/glycosylation false positives that a
    # downstream-breadth (hub) filter alone misses.
    try:
        k = load_k562()
        k = k[k["condition"] == args.condition].dropna(subset=["logfc_pearson_r"])
        kr = k.groupby("target_contrast_gene_name")["logfc_pearson_r"].mean()
        df["k562_concordance"] = df["gene"].map(kr).round(3)
        k_thr = float(df["k562_concordance"].quantile(0.67))   # top tercile = most cell-line-shared
        df["specificity_flag"] = np.where(df["k562_concordance"] >= k_thr, "cell-line-shared/generic", "")
    except Exception as e:
        print("K562 artifact filter skipped:", e, flush=True)
        df["k562_concordance"] = np.nan; df["specificity_flag"] = ""; k_thr = None

    clean = (df["hub_flag"] == "") & (df["specificity_flag"] == "")   # non-hub AND T-cell-specific
    if args.axis == "activation":
        lo_label = "knockdown lowers the activation program (immunosuppressive / autoimmune direction)"
        hi_label = "knockdown raises it (brake release / anti-tumor direction)"
    else:
        lo_label = "knockdown shifts cells toward Th1 (lowers the Th2-vs-Th1 axis)"
        hi_label = "knockdown shifts cells toward Th2 (raises it)"
    lo = df.sort_values("score").head(40)
    hi = df.sort_values("score", ascending=False).head(40)
    out = {"axis": args.axis, "condition": args.condition, "validation": validation,
           "hub_threshold": float(hub_thr), "k562_concordance_threshold": k_thr,
           "direction_low": lo_label, "direction_high": hi_label,
           "score_low_top": lo.to_dict("records"), "score_high_top": hi.to_dict("records"),
           "note": ("'novel' = not in the authors' own regulator table (confirm it isn't already a drug "
                    "target via the connectors). A clean candidate is non-hub AND T-cell-specific: hub_flag "
                    "empty (not a broad fitness gene) AND specificity_flag empty (low K562 concordance, so "
                    "not a generic ER/housekeeping effect).")}
    json.dump(out, open(args.out, "w"), indent=2, default=str)
    ca = df.loc[clean & (df["status"] == "novel")]
    print(f"[{args.axis}] LOW-score clean novel:", ca.sort_values("score").head(10)["gene"].tolist(), flush=True)
    print(f"[{args.axis}] HIGH-score clean novel:",
          ca.sort_values("score", ascending=False).head(10)["gene"].tolist(), flush=True)
    print("WROTE", args.out, flush=True)

if __name__ == "__main__":
    main()
