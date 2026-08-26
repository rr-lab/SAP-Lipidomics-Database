#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from copy import deepcopy
from pathlib import Path


ROOT = Path("/Users/nirwantandukar/Documents/Github/SoLD_paper")
THESIS_ROOT = Path(
    "/Users/nirwantandukar/Documents/Github/Thesis_NirwanTandukar_GeneticsGenomics/ncsuthesis-0.6"
)


def fmt_tag(p: float) -> str:
    return f"{p:.3e}".replace("+", "")


def fmt_pval_list(p: float) -> str:
    return f"{p:.1e}"


CONFIGS = [
    {
        "label": "CTL individual",
        "threshold": 0.05 / 6106358,
        "table_path": THESIS_ROOT
        / "Supplementary_tables/Chapter2/SuppTable_S7_CTL_GWAS_candidate_genes_for_individual_lipid_traits.tsv",
        "gene_by_path": ROOT / "results/gene_count/individual/control/gene_by_phenotype_bestP_p<=1e-07.tsv",
        "gene_summary_out": ROOT
        / f"results/gene_count/individual/control/control_individual_gene_summary_across_phenotypes_p<={fmt_tag(0.05 / 6106358)}.tsv",
        "gene_by_out": ROOT
        / f"results/gene_count/individual/control/gene_by_phenotype_bestP_p<={fmt_tag(0.05 / 6106358)}.tsv",
        "phenotype_summary_out": ROOT
        / f"results/gene_count/individual/control/phenotype_counts_sigGenes_p<={fmt_tag(0.05 / 6106358)}.tsv",
    },
    {
        "label": "CTL sum_ratio",
        "threshold": 0.05 / 6106358,
        "table_path": THESIS_ROOT
        / "Supplementary_tables/Chapter2/SuppTable_S8_CTL_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv",
        "gene_by_path": ROOT / "results/gene_count/sum_ratio/control/gene_by_phenotype_bestP_p<=1e-07.tsv",
        "gene_summary_out": ROOT
        / f"results/gene_count/sum_ratio/control/control_sum_ratio_gene_summary_across_phenotypes_p<={fmt_tag(0.05 / 6106358)}.tsv",
        "gene_by_out": ROOT
        / f"results/gene_count/sum_ratio/control/gene_by_phenotype_bestP_p<={fmt_tag(0.05 / 6106358)}.tsv",
        "phenotype_summary_out": ROOT
        / f"results/gene_count/sum_ratio/control/phenotype_counts_sigGenes_p<={fmt_tag(0.05 / 6106358)}.tsv",
    },
    {
        "label": "LIN individual",
        "threshold": 0.05 / 6085245,
        "table_path": THESIS_ROOT
        / "Supplementary_tables/Chapter2/SuppTable_S9_LIN_GWAS_candidate_genes_for_individual_lipid_traits.tsv",
        "gene_by_path": ROOT / "results/gene_count/individual/lowinput/gene_by_phenotype_bestP_p<=1e-07.tsv",
        "gene_summary_out": ROOT
        / f"results/gene_count/individual/lowinput/lowinput_individual_gene_summary_across_phenotypes_p<={fmt_tag(0.05 / 6085245)}.tsv",
        "gene_by_out": ROOT
        / f"results/gene_count/individual/lowinput/gene_by_phenotype_bestP_p<={fmt_tag(0.05 / 6085245)}.tsv",
        "phenotype_summary_out": ROOT
        / f"results/gene_count/individual/lowinput/phenotype_counts_sigGenes_p<={fmt_tag(0.05 / 6085245)}.tsv",
    },
    {
        "label": "LIN sum_ratio",
        "threshold": 0.05 / 6085245,
        "table_path": THESIS_ROOT
        / "Supplementary_tables/Chapter2/SuppTable_S10_LIN_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv",
        "gene_by_path": ROOT / "results/gene_count/sum_ratio/lowinput/gene_by_phenotype_bestP_p<=1e-07.tsv",
        "gene_summary_out": ROOT
        / f"results/gene_count/sum_ratio/lowinput/sum_ratio_gene_summary_across_phenotypes_p<={fmt_tag(0.05 / 6085245)}.tsv",
        "gene_by_out": ROOT
        / f"results/gene_count/sum_ratio/lowinput/gene_by_phenotype_bestP_p<={fmt_tag(0.05 / 6085245)}.tsv",
        "phenotype_summary_out": ROOT
        / f"results/gene_count/sum_ratio/lowinput/phenotype_counts_sigGenes_p<={fmt_tag(0.05 / 6085245)}.tsv",
    },
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ordered_unique(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def rebuild_one(cfg: dict[str, object]) -> dict[str, object]:
    threshold = float(cfg["threshold"])
    base_rows = read_tsv(Path(cfg["table_path"]))
    base_by_gene = {row["GeneID"]: row for row in base_rows}
    base_fields = list(base_rows[0].keys())

    gene_by_rows = read_tsv(Path(cfg["gene_by_path"]))
    kept_rows = [row for row in gene_by_rows if float(row["Best_P_Value"]) <= threshold]

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    pheno_grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in kept_rows:
        grouped[row["GeneID"]].append(row)
        pheno_grouped[row["Phenotype"]].append(row)

    rebuilt_rows: list[dict[str, str]] = []
    gene_summary_rows: list[dict[str, str]] = []

    for gene_id, hits in grouped.items():
        hits_by_p = sorted(hits, key=lambda r: (float(r["Best_P_Value"]), r["Phenotype"], r["Best_SNP"]))
        best_hit = hits_by_p[0]
        phenotypes = sorted({row["Phenotype"] for row in hits})
        min_p_by_pheno = {
            phenotype: min(float(row["Best_P_Value"]) for row in hits if row["Phenotype"] == phenotype)
            for phenotype in phenotypes
        }
        best_snp_rows = ordered_unique([row["Best_SNP"] for row in hits_by_p])

        base = deepcopy(base_by_gene[gene_id])
        base["Best_SNP"] = best_hit["Best_SNP"]
        base["Best_P_Value"] = best_hit["Best_P_Value"]
        base["N_Phenotypes"] = str(len(phenotypes))
        base["Phenotypes"] = "; ".join(phenotypes)
        base["All_SNPs"] = "; ".join(best_snp_rows)
        base["N_SNPs"] = str(len(best_snp_rows))
        base["P_Values"] = "; ".join(fmt_pval_list(min_p_by_pheno[p]) for p in phenotypes)
        rebuilt_rows.append(base)

        gene_summary_rows.append(
            {
                "GeneID": gene_id,
                "Chromosome": base["Chromosome"],
                "Gene_Start": base["Gene_Start"],
                "Gene_End": base["Gene_End"],
                "Best_SNP": best_hit["Best_SNP"],
                "Best_SNP_Position": base["Best_SNP_Position"],
                "Best_P_Value": best_hit["Best_P_Value"],
                "N_Phenotypes": str(len(phenotypes)),
                "Phenotypes": "; ".join(phenotypes),
                "All_SNPs": "; ".join(best_snp_rows),
                "N_SNPs": str(len(best_snp_rows)),
                "P_Values": "; ".join(fmt_pval_list(min_p_by_pheno[p]) for p in phenotypes),
            }
        )

    rebuilt_rows.sort(key=lambda r: (float(r["Best_P_Value"]), -int(r["N_Phenotypes"]), r["GeneID"]))
    gene_summary_rows.sort(key=lambda r: (float(r["Best_P_Value"]), -int(r["N_Phenotypes"]), r["GeneID"]))

    phenotype_summary_rows = []
    for phenotype, rows in sorted(pheno_grouped.items()):
        best_row = min(rows, key=lambda r: float(r["Best_P_Value"]))
        phenotype_summary_rows.append(
            {
                "Phenotype": phenotype,
                "N_Genes": str(len({row["GeneID"] for row in rows})),
                "N_SNPs": str(len({row["Best_SNP"] for row in rows})),
                "Min_P": best_row["Best_P_Value"],
                "Top_Gene": best_row["GeneID"],
            }
        )
    phenotype_summary_rows.sort(key=lambda r: (-int(r["N_Genes"]), float(r["Min_P"]), r["Phenotype"]))

    write_tsv(Path(cfg["table_path"]), rebuilt_rows, base_fields)
    write_tsv(
        Path(cfg["gene_summary_out"]),
        gene_summary_rows,
        [
            "GeneID",
            "Chromosome",
            "Gene_Start",
            "Gene_End",
            "Best_SNP",
            "Best_SNP_Position",
            "Best_P_Value",
            "N_Phenotypes",
            "Phenotypes",
            "All_SNPs",
            "N_SNPs",
            "P_Values",
        ],
    )
    write_tsv(
        Path(cfg["gene_by_out"]),
        sorted(kept_rows, key=lambda r: (r["Phenotype"], float(r["Best_P_Value"]), r["GeneID"])),
        list(gene_by_rows[0].keys()),
    )
    write_tsv(
        Path(cfg["phenotype_summary_out"]),
        phenotype_summary_rows,
        ["Phenotype", "N_Genes", "N_SNPs", "Min_P", "Top_Gene"],
    )

    return {
        "label": cfg["label"],
        "threshold": threshold,
        "rows": len(rebuilt_rows),
        "gene_pheno_rows": len(kept_rows),
        "phenotypes": len(phenotype_summary_rows),
    }


def main() -> None:
    summaries = [rebuild_one(cfg) for cfg in CONFIGS]
    print("Rebuilt GWAS candidate tables using condition-specific Bonferroni thresholds:")
    for summary in summaries:
        print(
            f"  - {summary['label']}: {summary['rows']} genes, "
            f"{summary['gene_pheno_rows']} gene-phenotype rows, "
            f"{summary['phenotypes']} phenotypes at p <= {summary['threshold']:.6e}"
        )


if __name__ == "__main__":
    main()
