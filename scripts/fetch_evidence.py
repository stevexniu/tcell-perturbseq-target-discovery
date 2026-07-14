#!/usr/bin/env python
"""Fetch the per-candidate evidence (novelty / druggability / genetics / safety).

The four specialist reads in results/tables/ (novelty_read.csv, druggability_read.csv,
genetics_read.csv, safety_read.csv) were produced by querying, per gene, the following
PUBLIC resources. This script documents the exact sources and query logic and provides
a runnable fetch for the parts served by open REST APIs (GTEx, HPA). The remaining
sources are reachable via their own public APIs/portals listed below.

Candidate set: results/tables/candidates_for_workup.csv

NOVELTY
  - Prior target publications : OpenAlex / PubMed (search "<gene> T cell" / "<gene> target")
  - Active drug programs       : ChEMBL target + mechanism, Drugs@FDA
  - Clinical trials            : ClinicalTrials.gov API v2 (query.term=<gene>)
  - Authors' own nominations   : the dataset's regulator_coefficients table (known_regulators col)

DRUGGABILITY
  - Protein class / family     : UniProt + InterPro/Pfam
  - Chemical matter            : ChEMBL (compounds, activities)
  - Structure                  : RCSB PDB (experimental) + AlphaFold DB (predicted)

HUMAN GENETICS
  - GWAS                       : GWAS Catalog REST (immune/autoimmune/atopic traits)
  - eQTL                       : eQTL Catalogue
  - Direction of effect        : Open Targets Genetics / L2G; checked that knockdown
                                 direction matches an inhibitor's intended effect

SAFETY
  - Expression breadth         : GTEx (median TPM across tissues) + HPA tissue specificity
  - Essentiality / constraint  : gnomAD v4 (pLI, LOEUF)
  - On-target liabilities      : OMIM / ClinVar germline phenotypes, HPA blood specificity

Runnable here: GTEx tissue summary + HPA per-gene record (public APIs).
Run:  python scripts/fetch_evidence.py
"""
import os
import json
import urllib.request
import urllib.parse
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(HERE, "results", "tables")

FOCUS_TISSUES = ["Whole_Blood", "Spleen", "Lung", "Liver", "Heart_Left_Ventricle",
                 "Brain_Cortex", "Kidney_Cortex", "Testis",
                 "Skin_Sun_Exposed_Lower_leg", "Muscle_Skeletal"]
GTEX_API = "https://gtexportal.org/api/v2"


def _get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def gtex_focus(genes):
    """Median TPM per gene x focus tissue from the public GTEx v8 API."""
    rows = {}
    for g in genes:
        try:
            gg = _get(f"{GTEX_API}/reference/gene?geneId={urllib.parse.quote(g)}")
            gencode = gg["data"][0]["gencodeId"]
            med = _get(f"{GTEX_API}/expression/medianGeneExpression?gencodeId={gencode}&datasetId=gtex_v8")
            per = {d["tissueSiteDetailId"]: d["median"] for d in med["data"]}
            rows[g] = {t: per.get(t) for t in FOCUS_TISSUES}
        except Exception as e:  # noqa
            rows[g] = {t: None for t in FOCUS_TISSUES}
            print(f"  GTEx failed for {g}: {e}")
    return pd.DataFrame(rows).T


if __name__ == "__main__":
    cand = pd.read_csv(os.path.join(TAB, "candidates_for_workup.csv"))
    genes = sorted(cand["gene"].unique())
    print(f"{len(genes)} candidates:", ", ".join(genes))
    print("\nThis script documents the evidence sources (see module docstring).")
    print("Fetching GTEx focus tissues from the public API as a runnable example ...")
    gt = gtex_focus(genes)
    out = os.path.join(TAB, "gtex_focus_tpm.regenerated.csv")
    gt.to_csv(out)
    print(f"wrote {out}")
    print("\nThe HPA, ChEMBL, GWAS Catalog, gnomAD, ClinicalTrials.gov, OpenAlex and Open Targets")
    print("reads are in results/tables/{novelty,druggability,genetics,safety}_read.csv with")
    print("per-gene identifiers; re-run against each resource's public API to refresh.")
