#!/usr/bin/env python
"""Method selection by the known-gene benchmark + held-out validation.
Reproduces data/benchmark_auroc.csv and data/benchmark_validation.json.
Requires the Stim8hr DE subset (zscore + log_fc); see README for how it is read over HTTP.
"""
import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

DRIVERS=["CD3D","CD3E","CD3G","CD247","ZAP70","LCK","LAT","PLCG1","VAV1","ITK","LCP2","FYN","RASGRP1","THEMIS"]
BRAKES =["CBLB","RASA2","SOCS1","CISH","PTPN2","TNFAIP3","SOCS3","DGKA","DGKZ","UBASH3A","TNFAIP8L2"]

def score_perturbations(Z, LFC, sig_w):
    """Four candidate scores: {z,log_fc} x {mean_signed, correlation} of the perturbation
    DE profile against the activation-signature weight vector sig_w (aligned to columns)."""
    out={}
    for name, M in [("z",Z),("log_fc",LFC)]:
        ms=(M*sig_w).mean(1)                     # mean_signed
        cor=np.array([np.corrcoef(M[i],sig_w)[0,1] for i in range(M.shape[0])])  # correlation
        out[(name,"mean_signed")]=ms; out[(name,"correlation")]=cor
    return out

def auroc(scores, genes, pos, null):
    idx={g:i for i,g in enumerate(genes)}
    p=[i for g in pos if (i:=idx.get(g)) is not None]
    n=[i for g in null if (i:=idx.get(g)) is not None]
    y=np.r_[np.ones(len(p)),np.zeros(len(n))]; s=np.r_[scores[p],scores[n]]
    return roc_auc_score(y, s)

if __name__=="__main__":
    print("See README: this script documents the benchmark. Precomputed results in "
          "data/benchmark_auroc.csv and data/benchmark_validation.json.")
