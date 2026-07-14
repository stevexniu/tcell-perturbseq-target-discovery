#!/usr/bin/env python
"""Regenerate all figures from the tables in results/tables/.

Reproduces:
  Fig 1  validation_scatter.png   (from validation_merge_{activation,polarization}.csv if present,
                                   else re-derivable via run_discovery.py)
  Fig 2  fig2_ranking.png         (ranked_polarization_Stim8hr_full.csv)
  Fig 3  fig3_leads.png           (polar_score_by_condition.csv, gtex_focus_tpm.csv, hpa_summary.csv)
  coverage_heatmap.png            (coverage_table.csv)

Run:  python scripts/make_figures.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(HERE, "results", "tables")
FIG = os.path.join(HERE, "results", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 110})
META_GREY = "#8a9299"

NOMINATED_POL = None  # filled from candidates_for_workup.csv
RECOVERED = {"TNFRSF9", "RAC2", "IFNAR1", "CBLB", "CD6"}
NEWLEADS = {"TATDN2", "JADE2"}


def fig2_ranking():
    q = pd.read_csv(os.path.join(TAB, "ranked_polarization_Stim8hr_full.csv"))
    for c in ("hub_flag", "specificity_flag"):
        q[c] = q[c].fillna("").astype(str).replace("nan", "")
    q = q[q.qc_pass].sort_values("score").reset_index(drop=True)
    q["rank"] = np.arange(len(q))
    cand = pd.read_csv(os.path.join(TAB, "candidates_for_workup.csv"))
    nominated = set(cand[cand.axis == "polarization"].gene)
    q["is_hub"] = q.hub_flag != ""
    q["cls"] = np.where(q.is_hub, "hub", np.where(q.gene.isin(nominated), "nominated", "other"))

    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    oth = q[q.cls == "other"]; ax.scatter(oth["rank"], oth["score"], s=9, c="#b8c4cc", lw=0, label="QC-pass (other)", zorder=1)
    hub = q[q.cls == "hub"]; ax.scatter(hub["rank"], hub["score"], s=16, c="#e2e6e9", lw=0.3, edgecolor="#c7c7c7", label="hub / fitness (greyed)", zorder=2)
    nom = q[q.cls == "nominated"]
    th1 = nom[nom.score < 0].sort_values("score"); th2 = nom[nom.score > 0].sort_values("score", ascending=False)
    ax.scatter(th1["rank"], th1["score"], s=44, c="#2166ac", edgecolor="white", lw=0.5, label="nominated — Th1-skewing", zorder=4)
    ax.scatter(th2["rank"], th2["score"], s=44, c="#b2182b", edgecolor="white", lw=0.5, label="nominated — Th2-skewing", zorder=4)
    ax.axhline(0, color="0.6", lw=0.6)

    def stack(sub, x_anchor, y_top, y_step, ha, color):
        for k, (_, r) in enumerate(sub.iterrows()):
            tag = r.gene + (" (known)" if r.gene in RECOVERED else "")
            fw = "bold" if r.gene in NEWLEADS else "normal"
            ax.annotate(tag, (r["rank"], r["score"]), xytext=(x_anchor, y_top - k * y_step),
                        textcoords="data", fontsize=6.6, fontweight=fw, ha=ha, va="center", color="0.15",
                        arrowprops=dict(arrowstyle="-", color=color, lw=0.5, alpha=0.7), zorder=6)
    # order each stack top->bottom by descending score so label vertical order
    # matches point vertical order (no crossing leader lines)
    stack(th1.sort_values("score", ascending=False), 70, -0.050, 0.0092, "left", "#2166ac")
    stack(th2.sort_values("score", ascending=False), 420, 0.108, 0.0092, "right", "#b2182b")
    ax.set_xlabel("rank (QC-pass perturbations, ascending polarization score)")
    ax.set_ylabel("polarization axis score\n(\u2212 = Th1-skewing     +  = Th2-skewing)")
    ax.set_title("Ranked polarization-axis candidates (Stim8hr) — nominations highlighted, hub/fitness greyed")
    ax.legend(loc="lower right", fontsize=6.2, frameon=False, markerscale=1.2)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_ranking.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig3_leads():
    gt = pd.read_csv(os.path.join(TAB, "gtex_focus_tpm.csv"), index_col=0)
    hpa = pd.read_csv(os.path.join(TAB, "hpa_summary.csv"), index_col=0)
    pc = pd.read_csv(os.path.join(TAB, "polar_score_by_condition.csv"), index_col=0)
    order = ["TNFRSF9", "RAC2", "JADE2", "BAHD1", "INTS6", "URM1",
             "IFNAR1", "TATDN2", "TRIT1", "TWF1", "SNRPB2", "PAXIP1", "CBLB", "CD6"]
    known = RECOVERED
    grp = [6, 12]
    gt = gt.loc[order]; pc = pc.loc[order]; hpa = hpa.loc[order]
    lbl = {"Whole_Blood": "Whole blood", "Heart_Left_Ventricle": "Heart", "Brain_Cortex": "Brain",
           "Kidney_Cortex": "Kidney", "Skin_Sun_Exposed_Lower_leg": "Skin", "Muscle_Skeletal": "Muscle"}
    gt.columns = [lbl.get(c, c) for c in gt.columns]

    fig = plt.figure(figsize=(10.6, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[2.4, 0.7, 1.0], wspace=0.34)
    axA = fig.add_subplot(gs[0, 0])
    imA = axA.imshow(gt.values, aspect="auto", cmap="magma", norm=LogNorm(vmin=0.1, vmax=1000))
    axA.set_xticks(range(gt.shape[1])); axA.set_xticklabels(gt.columns, rotation=45, ha="right", fontsize=6.8)
    axA.set_yticks(range(len(order))); axA.set_yticklabels([g + ("*" if g in known else "") for g in order], fontsize=7)
    axA.set_title("A  GTEx tissue expression", fontsize=8.5, loc="left")
    for b in grp: axA.axhline(b - 0.5, color="white", lw=1.6)
    cbA = fig.colorbar(imA, ax=axA, fraction=0.035, pad=0.02); cbA.set_label("median TPM", fontsize=6.5); cbA.ax.tick_params(labelsize=6)

    axB = fig.add_subplot(gs[0, 1])
    spec = {"Tissue enriched": 0, "Group enriched": 1, "Tissue enhanced": 2, "Low tissue specificity": 3}
    blood = {"Group enriched": 0, "Low immune cell specificity": 2}
    B = np.array([[spec.get(hpa.loc[g, "rna_tissue_specificity"], 3), blood.get(hpa.loc[g, "rna_blood_specificity"], 2)] for g in order], float)
    cmapB = mpl.colors.ListedColormap(["#1a9850", "#91cf60", "#fee08b", "#d73027"])
    axB.imshow(B, aspect="auto", cmap=cmapB, vmin=0, vmax=3)
    axB.set_xticks([0, 1]); axB.set_xticklabels(["tissue\nspec.", "blood\nspec."], fontsize=6.5)
    axB.set_yticks([]); axB.set_title("B  HPA RNA specificity", fontsize=8.5, loc="left")
    for b in grp: axB.axhline(b - 0.5, color="white", lw=1.6)
    legB = [Patch(fc="#1a9850", label="enriched (restricted)"), Patch(fc="#91cf60", label="group enriched"),
            Patch(fc="#fee08b", label="enhanced"), Patch(fc="#d73027", label="low spec. (broad)")]
    axB.legend(handles=legB, fontsize=5.3, loc="upper center", bbox_to_anchor=(0.5, -0.10), frameon=False, handlelength=1)

    axC = fig.add_subplot(gs[0, 2])
    norm = TwoSlopeNorm(vcenter=0, vmin=-0.11, vmax=0.11)
    imC = axC.imshow(pc.values, aspect="auto", cmap="RdBu_r", norm=norm)
    axC.set_xticks(range(3)); axC.set_xticklabels(["Rest", "Stim\n8h", "Stim\n48h"], fontsize=6.8)
    axC.set_yticks([]); axC.set_title("C  KD effect on polarization", fontsize=8.5, loc="left")
    for b in grp: axC.axhline(b - 0.5, color="0.3", lw=1.6)
    for i in range(pc.shape[0]):
        for j in range(pc.shape[1]):
            axC.text(j, i, f"{pc.values[i, j]:+.02f}", ha="center", va="center", fontsize=5.3,
                     color="white" if abs(pc.values[i, j]) > 0.06 else "0.25")
    cbC = fig.colorbar(imC, ax=axC, fraction=0.06, pad=0.03); cbC.set_label("polar. score\n(\u2212 Th1  \u00b7  + Th2)", fontsize=6); cbC.ax.tick_params(labelsize=6)
    fig.suptitle("Nominated & recovered-known leads: expression breadth vs polarization effect   (* = recovered-known)", fontsize=9, y=1.02)
    fig.savefig(os.path.join(FIG, "fig3_leads.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig2_ranking(); print("wrote fig2_ranking.png")
    fig3_leads(); print("wrote fig3_leads.png")
    print("Note: Fig 1 (validation_scatter) and coverage_heatmap require validation_merge_*.csv /"
          " coverage_table.csv respectively; validation merges are produced by run_discovery.py.")
