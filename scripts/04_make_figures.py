#!/usr/bin/env python
"""Rebuild all three figures from the shipped data/ CSVs. Fully offline.
Usage: python scripts/04_make_figures.py   ->   figures/fig{1,2,3}_*.png
"""
import json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import pearsonr
mpl.rcParams.update({"figure.dpi":300,"savefig.dpi":300,"font.size":8,"axes.spines.top":False,
                     "axes.spines.right":False,"axes.linewidth":0.8})
GREY="#9a9a9a"; DATA="data"; FIG="figures"

def fig1():
    m=pd.read_csv(f"{DATA}/fig1_score_vs_authorcoef.csv")
    m["highlight"]=m["highlight"].fillna("").astype(str)
    kn=m[m.known!="other"]; r_kn=pearsonr(kn.my_score,kn.author_coef)[0]; r_all=pearsonr(m.my_score,m.author_coef)[0]
    fig,(ax,ax2)=plt.subplots(1,2,figsize=(180/25.4,80/25.4),gridspec_kw={"width_ratios":[1.6,1]})
    ax.scatter(m.my_score,m.author_coef,s=4,c=GREY,alpha=.35,lw=0)
    for k,c in [("driver","#2166AC"),("brake","#B2182B")]:
        s=kn[kn.known==k]; ax.scatter(s.my_score,s.author_coef,s=26,c=c,lw=.4,ec="w",label=k,zorder=3)
    for _,row in m[m.highlight!=""].iterrows():
        ax.annotate(row.highlight,(row.my_score,row.author_coef),fontsize=7,weight="bold",
                    xytext=(4,4),textcoords="offset points")
    ax.set_xlabel("Our activation score (Stim8hr)"); ax.set_ylabel("Authors' regulator coefficient")
    ax.legend(frameon=False,loc="lower right",fontsize=7)
    ax.set_title(f"Known-gene r={r_kn:.2f}   |   global r={r_all:.2f}  (n={len(m)})",fontsize=8)
    thr=[0.0,0.01,0.02,0.03]; rs=[]; ns=[]
    for t in thr:
        sub=m[m.author_coef.abs()>=t]; rs.append(pearsonr(sub.my_score,sub.author_coef)[0]); ns.append(len(sub))
    ax2.plot(thr,rs,"o-",c="#333"); 
    for t,r,n in zip(thr,rs,ns): ax2.annotate(f"n={n}",(t,r),fontsize=6,xytext=(0,5),textcoords="offset points",ha="center")
    ax2.set_xlabel("|author coef| threshold"); ax2.set_ylabel("Pearson r"); ax2.set_title("Agreement rises with effect size",fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig1_validation.png",bbox_inches="tight"); plt.close(fig)

def fig2():
    d=pd.read_csv(f"{DATA}/fig2_ranking.csv").sort_values("score").reset_index(drop=True)
    d["highlight"]=d["highlight"].fillna("").astype(str); d["anchor"]=d["anchor"].fillna("").astype(str)
    fig,ax=plt.subplots(figsize=(180/25.4,90/25.4)); x=np.arange(len(d))
    col=np.where(d.direction.str.startswith("auto"),"#2166AC","#B2182B")
    col=np.where(d.is_hub,GREY,col)
    ax.scatter(x,d.score,s=10,c=col,lw=0,alpha=.8)
    ax.axhline(0,color="#666",lw=.8,ls="--")
    for _,row in d[(d.highlight!="")|(d.anchor!="")].iterrows():
        lab=row.highlight or row.anchor; xi=d.index[d.gene==row.gene][0]
        ax.annotate(lab,(xi,row.score),fontsize=7,weight="bold" if row.highlight else "normal",
                    xytext=(0,8 if row.score>0 else -12),textcoords="offset points",ha="center",
                    arrowprops=dict(arrowstyle="-",lw=.5,color="#333"))
    ax.set_xlabel("Perturbation rank (n=487 QC-passing)"); ax.set_ylabel("Activation score")
    ax.set_title("Both-direction ranking  (down=autoimmune, up=anti-tumor; hubs greyed)",fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig2_ranking.png",bbox_inches="tight"); plt.close(fig)

def fig3():
    shift=pd.read_csv(f"{DATA}/fig3_activation_shift.csv"); gtex=pd.read_csv(f"{DATA}/fig3_gtex.csv")
    hpa=pd.read_csv(f"{DATA}/fig3_hpa.csv"); COL={"RSBN1L":"#2166AC","MAP3K1":"#B2182B"}
    conds=["Rest","Stim8hr","Stim48hr"]
    fig=plt.figure(figsize=(180/25.4,150/25.4)); gs=fig.add_gridspec(2,3,width_ratios=[1.15,.85,1],hspace=.55,wspace=.5)
    for i,g in enumerate(["RSBN1L","MAP3K1"]):
        base=COL[g]
        axg=fig.add_subplot(gs[i,0]); sub=gtex[gtex.gene==g].sort_values("median_tpm",ascending=False)
        top=sub.head(10).iloc[::-1]; wb=sub[sub.tissue=="Whole_Blood"]
        yp=list(range(len(top)))+[len(top)+.6]; vals=list(top.median_tpm)+[wb.median_tpm.values[0]]
        cols=[mpl.colors.to_rgba(base,.35)]*len(top)+[base]
        axg.barh(yp,vals,color=cols,height=.7)
        axg.set_yticks(yp); axg.set_yticklabels([t.replace("_"," ") for t in top.tissue]+["Whole blood"],fontsize=6)
        axg.set_xlabel("Median TPM (GTEx)"); axg.set_title(g,color=base,weight="bold")
        axh=fig.add_subplot(gs[i,1]); axh.axis("off"); h=hpa[hpa.gene==g].iloc[0]
        fields=[("RNA tissue",h.RNA_tissue_specificity),("Protein tissue",h.protein_tissue_specificity),
                ("Immune-cell",h.RNA_blood_cell_specificity)]
        y=.85
        for lab,val in fields:
            axh.text(0,y,lab,fontsize=5.8,color=GREY,transform=axh.transAxes)
            axh.text(0,y-.06,str(val),fontsize=6.6,weight="bold",transform=axh.transAxes); y-=.2
        axh.set_title(g,color=base,weight="bold")
        axs=fig.add_subplot(gs[i,2]); s=shift[shift.gene==g].set_index("condition").loc[conds]
        yerr=np.vstack([s.activation_score-s.ci_lo,s.ci_hi-s.activation_score]); xx=np.arange(3)
        axs.axhline(0,color=GREY,lw=.8,ls="--")
        axs.errorbar(xx,s.activation_score,yerr=yerr,color=base,marker="o",ms=6,lw=2,capsize=3)
        axs.set_xticks(xx); axs.set_xticklabels(["Rest","Stim 8h","Stim 48h"]); axs.set_ylim(-.5,.45)
        axs.set_ylabel("Activation score (KD-ctrl)"); axs.set_title(g,color=base,weight="bold")
    fig.suptitle("RSBN1L & MAP3K1: expression and activation shift",fontsize=9,y=.99)
    fig.savefig(f"{FIG}/fig3_leads.png",bbox_inches="tight"); plt.close(fig)

if __name__=="__main__":
    import os; os.makedirs(FIG,exist_ok=True)
    fig1(); fig2(); fig3(); print("wrote figures/fig{1,2,3}_*.png")
