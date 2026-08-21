#!/usr/bin/env python3

"""GO-BP enrichment stratified by lipid class for individual-lipid GWAS.

Candidate genes are assigned to a class when at least one of their associated
individual lipid phenotypes belongs to that class. CTL and LIN are analyzed
separately. Class sums and class ratios are intentionally excluded because a
ratio belongs to two classes rather than one.

All eligible class-by-term tests are adjusted together with BH correction.
"""

import csv
import math
import os
import re
from collections import defaultdict

REPO = "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database"
ANNOT = f"{REPO}/data/annotation/gene_annotation.txt"
MASTER = f"{REPO}/data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv"
OUT = f"{REPO}/data/LD_mapped/go_enrichment"
os.makedirs(OUT, exist_ok=True)


def parse_go_bp(value):
    if not value:
        return []
    terms = []
    for token in value.split(";"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r"^(.*)\((GO:[0-9]+)\)$", token)
        terms.append((match.group(1).strip(), match.group(2)) if match else (token, ""))
    return terms


def phenotype_class(phenotype):
    match = re.match(r"^([A-Za-z0-9]+)\(", phenotype.strip())
    return match.group(1) if match else None


def bh_adjust(pvalues):
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted = [0.0] * len(pvalues)
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = len(pvalues) - reverse_rank + 1
        running = min(running, pvalues[index] * len(pvalues) / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def fisher_greater(a, n_tested, term_bg, n_background):
    """One-sided Fisher exact p-value using the hypergeometric tail."""
    upper = min(n_tested, term_bg)
    denominator = math.comb(n_background, n_tested)
    return sum(
        math.comb(term_bg, overlap) *
        math.comb(n_background - term_bg, n_tested - overlap) /
        denominator
        for overlap in range(a, upper + 1)
    )


# GO background and term membership.
gene_terms = defaultdict(set)
with open(ANNOT, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        gene = row["GeneID"].strip()
        gene_terms[gene].update(parse_go_bp(row.get("GO_BP", "")))

background = {gene for gene, terms in gene_terms.items() if terms}
n_background = len(background)
term_genes = defaultdict(set)
term_ids = {}
for gene in background:
    for term, go_id in gene_terms[gene]:
        term_genes[term].add(gene)
        if go_id:
            term_ids[term] = go_id

# Build candidate genes per condition and lipid class from individual traits.
class_genes = defaultdict(lambda: defaultdict(set))
class_phenotypes = defaultdict(lambda: defaultdict(set))

with open(MASTER, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["layer"].strip().lower() != "individual":
            continue
        condition = row["condition"].strip()
        gene = row["GeneID"].strip()
        phenotypes = [p.strip() for p in row.get("Phenotypes", "").split(";") if p.strip()]
        for phenotype in phenotypes:
            lipid_class = phenotype_class(phenotype)
            if not lipid_class:
                continue
            class_genes[condition][lipid_class].add(gene)
            class_phenotypes[condition][lipid_class].add(phenotype)

all_rows = []
summary = []

for condition in sorted(class_genes):
    for lipid_class in sorted(class_genes[condition]):
        candidates = class_genes[condition][lipid_class]
        testable = candidates & background
        summary.append({
            "condition": condition,
            "lipid_class": lipid_class,
            "individual_phenotypes": len(class_phenotypes[condition][lipid_class]),
            "candidate_genes": len(candidates),
            "GO_BP_testable_genes": len(testable),
        })

        n_testable = len(testable)
        for term, term_set in term_genes.items():
            term_bg = len(term_set)
            overlap = sorted(testable & term_set)
            n_overlap = len(overlap)
            if term_bg < 5 or n_overlap < 3:
                continue

            a = n_overlap
            pvalue = fisher_greater(a, n_testable, term_bg, n_background)
            expected = n_testable * term_bg / n_background
            all_rows.append({
                "condition": condition,
                "lipid_class": lipid_class,
                "go_term": term,
                "go_id": term_ids.get(term, ""),
                "term_bg_count": term_bg,
                "term_test_count": n_overlap,
                "expected": expected,
                "fold_enrichment": n_overlap / expected,
                "p_value": pvalue,
                "overlap_genes": ";".join(overlap),
            })

# One global correction across every class, condition, and GO term test.
adjusted = bh_adjust([row["p_value"] for row in all_rows])
for row, p_adj in zip(all_rows, adjusted):
    row["p_adj_bh"] = p_adj

all_rows.sort(key=lambda row: (row["p_adj_bh"], row["p_value"], -row["fold_enrichment"]))
result_file = f"{OUT}/GO_BP_enrichment_by_lipid_class_individual_LD.tsv"
with open(result_file, "w", newline="") as handle:
    fields = [
        "condition", "lipid_class", "go_term", "go_id", "term_bg_count",
        "term_test_count", "expected", "fold_enrichment", "p_value",
        "p_adj_bh", "overlap_genes",
    ]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for row in all_rows:
        writer.writerow(row)

summary_file = f"{OUT}/GO_BP_enrichment_by_lipid_class_individual_LD_summary.tsv"
with open(summary_file, "w", newline="") as handle:
    fields = [
        "condition", "lipid_class", "individual_phenotypes", "candidate_genes",
        "GO_BP_testable_genes", "GO_terms_tested", "BH_significant_terms",
    ]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for row in summary:
        matching = [
            result for result in all_rows
            if result["condition"] == row["condition"]
            and result["lipid_class"] == row["lipid_class"]
        ]
        row["GO_terms_tested"] = len(matching)
        row["BH_significant_terms"] = sum(result["p_adj_bh"] < 0.05 for result in matching)
        writer.writerow(row)

print(f"GO background genes: {n_background}")
print(f"Total class-by-term tests: {len(all_rows)}")
print(f"Global BH-significant tests: {sum(row['p_adj_bh'] < 0.05 for row in all_rows)}")
for row in summary:
    print(
        f"{row['condition']} {row['lipid_class']}: "
        f"{row['candidate_genes']} candidate genes, "
        f"{row['GO_BP_testable_genes']} GO-testable, "
        f"{row['BH_significant_terms']} significant terms"
    )
print(f"Results: {result_file}")
print(f"Summary: {summary_file}")
