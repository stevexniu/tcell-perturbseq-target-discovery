#!/usr/bin/env python
"""Extract per-condition polarization scores for the nominated genes.

Reads only the ~42 relevant zscore rows from the 16.8 GB DE h5ad directly over
S3 (anonymous) — no full download needed — and projects each (gene, condition)
onto the Ota-2021 Th1/Th2 polarization signature using the exact same
projection math as discover.py / run_discovery.py.

Output: results/tables/polar_score_by_condition.csv

Run:  python scripts/extract_per_condition.py
"""
import os
import numpy as np
import pandas as pd
import h5py
import s3fs

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(HERE, "results", "tables")
os.makedirs(TAB, exist_ok=True)

BUCKET = "genome-scale-tcell-perturb-seq"
KEY = "marson2025_data/GWCD4i.DE_stats.h5ad"
POLAR_URL = ("https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/"
             "master/metadata/suppl_tables/Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv")

GENES = ["TATDN2", "JADE2", "URM1", "BAHD1", "INTS6", "RAC2",
         "IFNAR1", "TNFRSF9", "CD6", "CBLB", "TRIT1", "SNRPB2", "PAXIP1", "TWF1"]
CONDS = ["Rest", "Stim8hr", "Stim48hr"]


def _dec(x):
    return x.decode() if isinstance(x, bytes) else x


def read_col(obs, name):
    node = obs[name]
    if isinstance(node, h5py.Group):
        cats = np.array([_dec(c) for c in node["categories"][:]], dtype=object)
        codes = node["codes"][:]; out = np.empty(codes.shape, dtype=object)
        v = codes >= 0; out[v] = cats[codes[v]]; out[~v] = None; return out
    a = node[:]
    return np.array([_dec(x) for x in a], dtype=object) if a.dtype.kind in ("S", "O") else a


def project_row(z, w):
    idx = np.where(w != 0)[0]; wv = w[idx]; wv = (wv - wv.mean()) / (wv.std() + 1e-9)
    s = z[idx]; s = (s - s.mean()) / (s.std() + 1e-9)
    return float((s @ wv) / len(idx))


def main():
    fs = s3fs.S3FileSystem(anon=True)
    with fs.open(f"{BUCKET}/{KEY}", "rb") as fo:
        f = h5py.File(fo, "r"); obs = f["obs"]
        tgt = read_col(obs, "target_contrast_gene_name").astype(str)
        cond = read_col(obs, "culture_condition").astype(str)
        var_names = np.array([_dec(x) for x in f["var"]["gene_name"][:]], dtype=object)
        idx_map = {}
        for g in GENES:
            for c in CONDS:
                hit = np.where((tgt == g) & (cond == c))[0]
                if len(hit):
                    idx_map[(g, c)] = int(hit[0])
        Zlayer = f["layers"]["zscore"]
        rows = {i: Zlayer[i, :].astype(np.float32) for i in sorted(set(idx_map.values()))}
        f.close()

    t = pd.read_csv(POLAR_URL)
    t = t[t["contrast"].astype(str) == "Th2_vs_Th1 (Ota 2021)"]
    wmap = dict(zip(t["variable"].astype(str), t["zscore"].astype(float)))
    w = np.array([wmap.get(g, 0.0) for g in var_names], np.float32)

    recs = [{"gene": g, "condition": c, "polar_score": project_row(rows[i], w)}
            for (g, c), i in idx_map.items()]
    pc = pd.DataFrame(recs).pivot(index="gene", columns="condition", values="polar_score")[CONDS]
    pc.to_csv(os.path.join(TAB, "polar_score_by_condition.csv"))
    print(f"wrote polar_score_by_condition.csv  ({pc.shape[0]} genes x {pc.shape[1]} conditions)")


if __name__ == "__main__":
    main()
