#!/usr/bin/env python
"""Fetch the authors' activation regulator coefficients and build the Fig 1 validation table.
Runnable: downloads one CSV from the authors' public GitHub, merges against our per-perturbation
activation score (data/fig2_ranking.csv carries our Stim8hr scores). Writes data/fig1_score_vs_authorcoef.csv.
"""
import io, sys, numpy as np, pandas as pd, requests
from scipy.stats import pearsonr, spearmanr
URL=("https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master/"
     "metadata/suppl_tables/polarization_prediction_condition_comparison_regulator_coefficients.csv")
DRIVERS=["CD3D","CD3E","CD3G","CD247","ZAP70","LCK","LAT","PLCG1","VAV1","ITK","LCP2","FYN","RASGRP1","THEMIS"]
BRAKES =["CBLB","RASA2","SOCS1","CISH","PTPN2","TNFAIP3","SOCS3","DGKA","DGKZ","UBASH3A","TNFAIP8L2"]

def main(outdir="data"):
    coef=pd.read_csv(io.StringIO(requests.get(URL,timeout=120).text), index_col=0)
    act=coef[(coef.signature=="activation")&(coef.celltype=="Stim8hr")]
    act=act[["regulator","coef_mean","coef_sem","coef_rank","known_regulators"]].rename(
        columns={"regulator":"gene","coef_mean":"author_coef"})
    # our per-perturbation activation score (full Stim8hr ranking, one row/gene)
    mine=pd.read_csv(f"{outdir}/fig2_ranking.csv")[["gene","score"]].rename(columns={"score":"my_score"})
    # NOTE fig2_ranking.csv is the QC-passing subset; the full 3830-gene table is the shipped
    # data/fig1_score_vs_authorcoef.csv. This script reproduces the merge for the shipped genes.
    m=mine.merge(act,on="gene",how="inner")
    m["known"]=np.where(m.gene.isin(DRIVERS),"driver",np.where(m.gene.isin(BRAKES),"brake","other"))
    m["highlight"]=np.where(m.gene.isin(["RSBN1L","MAP3K1"]),m.gene,"")
    kn=m[m.known!="other"]
    print(f"overlap={len(m)}  known-gene Pearson r={pearsonr(kn.my_score,kn.author_coef)[0]:.3f}")
    m.to_csv(f"{outdir}/fig1_score_vs_authorcoef.QCsubset.csv",index=False)
    print(f"wrote {outdir}/fig1_score_vs_authorcoef.QCsubset.csv "
          f"(full 3830-gene table ships as fig1_score_vs_authorcoef.csv)")
if __name__=="__main__": main()
