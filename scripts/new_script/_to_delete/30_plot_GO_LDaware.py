#!/usr/bin/env python3
"""Figures for class-stratified GO-BP enrichment with LD-aware validation.

One figure per trait layer, three panels sharing the GO-term axis:
  A  fold enrichment            (dot size = candidate genes carrying the term)
  B  independent LD blocks      (how much of the support is distributed)
  C  significance before/after  (gene-level q -> LD-aware permutation q)

Single-locus terms (1 LD block) are drawn as hollow markers throughout, so
"this rests on one tandem array" is encoded by shape as well as by panel B.
"""
import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import pathlib
REPO = str(pathlib.Path(__file__).resolve().parents[2])
GO   = f"{REPO}/data/LD_mapped/go_enrichment"
FIG  = f"{REPO}/fig/go_ldaware"; os.makedirs(FIG, exist_ok=True)

SURF, INK, INK2, MUT = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8985"
COL  = {"CTL": "#2a78d6", "LIN": "#eb6834"}          # validated categorical slots 1,2

def theme(ax, xlab=None, title=None):
    ax.set_facecolor(SURF)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUT); ax.spines[s].set_linewidth(0.8)
    ax.grid(axis="x", color="#e6e5e1", lw=0.6, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3, width=0.8)
    if xlab:  ax.set_xlabel(xlab, fontsize=9.5, color=INK, fontweight="bold")
    if title: ax.set_title(title, fontsize=10.5, color=INK, fontweight="bold", loc="left", pad=8)

def build(layer, title, outfile, onto="BP"):
    rows = list(csv.DictReader(open(f"{GO}/GO_{onto}_enrichment_by_lipid_class_{layer}_LDaware.tsv"), delimiter="\t"))
    for r in rows:
        r["fold"]   = float(r["fold_enrichment"]); r["genes"] = int(r["term_test_count"])
        r["blocks"] = int(r["blocks_carrying_term"])
        r["q_gene"] = float(r["p_adj_bh"]);        r["q_perm"] = float(r["q_perm_bh"])
    # order: condition, then lipid class, then fold
    rows.sort(key=lambda r: (r["condition"], r["lipid_class"], -r["fold"]))
    rows = rows[::-1]                                    # bottom-up plotting
    y = np.arange(len(rows))
    sizes = np.interp([r["genes"] for r in rows], [3, 31], [34, 300])

    h = max(4.2, 0.30 * len(rows) + 2.0)
    fig, axes = plt.subplots(1, 3, figsize=(14.6, h), sharey=True,
                             gridspec_kw={"width_ratios": [2.5, 1.0, 1.35], "wspace": 0.06})
    fig.patch.set_facecolor(SURF)

    # ---- A  fold enrichment
    ax = axes[0]
    for i, r in enumerate(rows):
        c = COL[r["condition"]]; solid = r["blocks"] >= 3
        ax.plot([0, r["fold"]], [y[i], y[i]], color=c, lw=1.4, alpha=.45, zorder=2, solid_capstyle="round")
        ax.scatter(r["fold"], y[i], s=sizes[i], zorder=3, linewidths=1.3,
                   facecolor=c if solid else SURF, edgecolor=c)
    ax.axvline(1, color=MUT, ls=(0, (4, 3)), lw=1.0, zorder=1)
    ax.set_xscale("log"); ax.set_xlim(0.8, max(r["fold"] for r in rows) * 1.9)
    ax.set_xticks([1, 2, 5, 10, 20, 40]); ax.set_xticklabels(["1", "2", "5", "10", "20", "40"])
    ax.set_yticks(y); ax.set_yticklabels([f'{r["lipid_class"]:>5}  {r["go_term"][:55]}' for r in rows],
                       fontsize=8.2, color=INK, fontfamily="DejaVu Sans Mono")
    theme(ax, "Fold enrichment (log scale)",
          "A   Enrichment of GO %s terms" % ("biological-process" if onto=="BP" else "molecular-function"))
    ax.set_ylim(-0.8, len(rows) - 0.2)

    # ---- B  independent LD blocks
    ax = axes[1]
    for i, r in enumerate(rows):
        c = COL[r["condition"]]; solid = r["blocks"] >= 3
        ax.plot([0, r["blocks"]], [y[i], y[i]], color=c, lw=1.4, alpha=.45, zorder=2, solid_capstyle="round")
        ax.scatter(r["blocks"], y[i], s=58, zorder=3, linewidths=1.3,
                   facecolor=c if solid else SURF, edgecolor=c)
    ax.axvline(3, color=MUT, ls=(0, (4, 3)), lw=1.0, zorder=1)
    ax.text(3, len(rows) - 0.1, " 3 intervals", fontsize=7.6, color=INK2, va="top")
    ax.set_xlim(0, max(r["blocks"] for r in rows) + 1.6)
    theme(ax, "Independent 250 kb intervals", "B   Distributed support")

    # ---- C  gene-level q vs LD-aware q
    ax = axes[2]
    for i, r in enumerate(rows):
        a, b = -np.log10(r["q_gene"]), -np.log10(r["q_perm"])
        ax.plot([a, b], [y[i], y[i]], color=MUT, lw=1.1, alpha=.55, zorder=2, solid_capstyle="round")
        ax.scatter(a, y[i], s=30, marker="o", facecolor=SURF, edgecolor=MUT, linewidths=1.1, zorder=3)
        ax.scatter(b, y[i], s=52, marker="D", facecolor=COL[r["condition"]],
                   edgecolor=COL[r["condition"]], linewidths=1.0, zorder=4)
    ax.axvline(-np.log10(0.05), color="#c0392b", ls=(0, (4, 3)), lw=1.1, zorder=1)
    ax.text(-np.log10(0.05), len(rows) - 0.1, "  q = 0.05", fontsize=7.6, color="#c0392b", va="top")
    theme(ax, r"$-\log_{10}(q)$", "C   Before and after the LD-aware null")

    # separator between the CTL and LIN groups
    conds = [r["condition"] for r in rows]
    for i in range(1, len(conds)):
        if conds[i] != conds[i-1]:
            for a in axes:
                a.axhline(y[i] - 0.5, color=MUT, lw=0.8, alpha=.55, zorder=1)

    handles = [Line2D([], [], marker="o", ls="", mfc=COL["CTL"], mec=COL["CTL"], ms=8, label="CTL"),
               Line2D([], [], marker="o", ls="", mfc=COL["LIN"], mec=COL["LIN"], ms=8, label="LIN"),
               Line2D([], [], marker="o", ls="", mfc=SURF, mec=INK2, ms=8, label="carried by fewer than 3 independent 250 kb intervals"),
               Line2D([], [], marker="o", ls="", mfc=SURF, mec=MUT, ms=6, label="gene-level q"),
               Line2D([], [], marker="D", ls="", mfc=MUT, mec=MUT, ms=6, label="LD-aware q")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8.6,
               bbox_to_anchor=(0.5, -0.004))
    fig.suptitle(title, fontsize=12.5, fontweight="bold", color=INK, x=0.005, ha="left",
                 y=1 - 0.30 / h)
    fig.text(0.005, 1 - 0.62 / h, f"Dot area in panel A is the number of candidate genes carrying the term "
             f"({min(int(r['term_test_count']) for r in rows)}–{max(int(r['term_test_count']) for r in rows)}). "
             "Hollow markers rest on fewer than three independent 250 kb intervals.",
             fontsize=8.4, color=INK2, ha="left")
    fig.tight_layout(rect=[0, 0.045, 1, 1 - 0.80 / h])
    fig.savefig(outfile, dpi=300, facecolor=SURF, bbox_inches="tight")
    plt.close(fig)
    print(f"{outfile}  ({len(rows)} terms)")


ONTO = {"BP": ("biological-process", "GO-BP"), "MF": ("molecular-function", "GO-MF")}
for onto, (longname, short) in ONTO.items():
    src = f"{GO}/GO_{onto}_enrichment_by_lipid_class_individual_LDaware.tsv"
    if not os.path.exists(src):
        print(f"skip {onto}: no results on disk"); continue
    build("individual",
          f"{short} enrichment of individual-lipid GWAS candidates, by lipid class",
          f"{FIG}/SuppFig_GO{onto}_LDaware_Individual.png", onto=onto)
    build("sumratio",
          f"{short} enrichment of class-sum and class-ratio GWAS candidates, by lipid class",
          f"{FIG}/SuppFig_GO{onto}_LDaware_SumRatio.png", onto=onto)
