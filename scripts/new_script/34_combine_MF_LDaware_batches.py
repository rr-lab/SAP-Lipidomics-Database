#!/usr/bin/env python3
"""Combine independent permutation batches of 33_GO_MF_LDaware_permutation.py.

Each batch contributes N draws from the same LD-aware null with a different seed,
so exceedance counts add and the pooled p is (sum_exceed + 1) / (sum_perm + 1).
Used to reach 8,000 permutations in 2,000-draw chunks. BH is applied per layer.
"""
import csv, glob, json, pathlib, sys

REPO = pathlib.Path(__file__).resolve().parents[2]
GO = REPO / "data/LD_mapped/go_enrichment"
BATCHDIR = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/mfb")

def bh(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i]); q = [0.0]*len(ps); run = 1.0
    for r, i in enumerate(reversed(order), 1):
        run = min(run, ps[i]*len(ps)/(len(ps)-r+1)); q[i] = min(run, 1.0)
    return q

for L, pat in (("individual", "ind_*.json"), ("sumratio", "sr_*.json")):
    files = sorted(BATCHDIR.glob(pat))
    if not files: continue
    tot = {}
    for f in files:
        for k, (exceed, nperm, nullsum) in json.load(open(f)).items():
            e, n, s = tot.get(k, (0, 0, 0.0))
            tot[k] = (e + exceed, n + nperm, s + nullsum)
    path = GO / f"GO_MF_enrichment_by_lipid_class_{L}_LDaware.tsv"
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows:
        k = f"{L}|{r['condition']}|{r['lipid_class']}|{r['go_term']}"
        e, n, s = tot[k]
        r["p_perm"] = round((e + 1) / (n + 1), 6)
        r["null_mean"] = round(s / n, 2)
        r["n_permutations"] = n
    for r, q in zip(rows, bh([r["p_perm"] for r in rows])):
        r["q_perm_bh"] = round(q, 6)
    rows.sort(key=lambda r: r["p_perm"])
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    keep = [r for r in rows if r["q_perm_bh"] < 0.05]
    print(f"{L}: {len(files)} batches, {rows[0]['n_permutations']} permutations, "
          f"{len(rows)} terms -> {len(keep)} survive")
    print(f"  -> {path}")
