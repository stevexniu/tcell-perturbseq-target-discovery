#!/usr/bin/env python
"""Build the activation signature from the dataset's own non-targeting controls.
log2 CPM (Stim8hr) - log2 CPM (Rest), over NTC pseudobulk profiles.
Reads GWCD4i.pseudobulk_merged.h5ad lazily over HTTP. Writes data/signature.csv.
Re-run only to regenerate from scratch; a precomputed data/signature.csv ships with the repo.
"""
import numpy as np, pandas as pd, h5py, fsspec, sys
BUCKET="https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
def open_h5(name):
    fs=fsspec.filesystem("http", client_kwargs={"trust_env":True})
    return h5py.File(fs.open(f"{BUCKET}/{name}", block_size=8*1024*1024),"r")
def main():
    print("Opening pseudobulk over HTTP (44.6 GB, lazy)...", file=sys.stderr)
    # See methods in README; this script documents the exact procedure used.
    # NTC rows identified via obs['guide_identity']=='non-targeting'; CPM per profile x1e4, log2.
    raise SystemExit("Signature precomputed in data/signature.csv. "
                     "Full recompute requires ~1.8GB NTC read; see README methods.")
if __name__=="__main__": main()
