#!/usr/bin/env bash
# =============================================================================
# Build genes.range (chr, start, end, gene_id) from the Ensembl GFF3.
#
# One-off, runs in seconds on the login node. genes.range is the interval file
# every later step joins SNP positions against.
#
# Only nuclear chromosomes 1-10 are kept, because the GWAS chr column is 1-10
# and scaffolds cannot be matched to it. Gene ids are taken from gene_id= when
# present, otherwise from ID= with the "gene:" prefix stripped, which is how
# Ensembl writes SORBI_ ids.
#
# Usage   bash 10_build_genes_range.sh
# =============================================================================
set -euo pipefail

GFF="${GFF:-/rsstu/users/r/rrellan/DOE_CAREER/SAP/ref/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3}"
OUT_DIR="${OUT_DIR:-/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/gene_annotation_final}"
OUT="$OUT_DIR/genes.range"

[[ -f "$GFF" ]] || { echo "ERROR: GFF3 not found: $GFF" >&2; exit 1; }
mkdir -p "$OUT_DIR"

awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ { next }
  $3 == "gene" {
    chr = $1
    sub(/^[Cc]hromosome_?/, "", chr); sub(/^[Cc]hr_?/, "", chr)
    sub(/^0+/, "", chr)                       # 01 -> 1, the HapMap padding style
    if (chr !~ /^([1-9]|10)$/) { skipped++; next }
    id = ""
    n = split($9, a, ";")
    for (i = 1; i <= n; i++) {
      if (a[i] ~ /^gene_id=/) { id = substr(a[i], 9); break }
    }
    if (id == "") {
      for (i = 1; i <= n; i++) {
        if (a[i] ~ /^ID=/) { id = substr(a[i], 4); sub(/^gene:/, "", id); break }
      }
    }
    if (id == "") { noid++; next }
    print chr, $4, $5, id
    kept++
  }
  END {
    printf "genes kept (chr 1-10) : %d\n", kept+0   > "/dev/stderr"
    printf "genes on scaffolds    : %d\n", skipped+0 > "/dev/stderr"
    printf "gene lines with no id : %d\n", noid+0    > "/dev/stderr"
  }' "$GFF" | sort -k1,1n -k2,2n > "$OUT"

echo
echo "Saved: $OUT  ($(wc -l < "$OUT") genes)"
echo
echo "chromosome spread"
cut -f1 "$OUT" | sort -n | uniq -c | awk '{printf "  chr%-3s %s genes\n", $2, $1}'
echo
echo "first rows"
head -3 "$OUT"

# sanity: ids should look like SORBI_3001G103800
bad=$(awk -F'\t' '$4 !~ /^SORBI_/' "$OUT" | wc -l)
[[ "$bad" -eq 0 ]] && echo "all gene ids are SORBI_" \
  || echo "WARNING: $bad gene ids are not SORBI_, check the GFF3 attribute field" >&2
