#!/usr/bin/env bash
# =============================================================================
# Rebuild the high-heterozygosity SNP removal list and its counts.
#
# Reproduces filter_high_het_snps.R exactly
#   het codes           Y R S W K M
#   genotype columns    12 .. NF
#   het fraction        (het calls) / (number of genotype columns)
#                       -- the denominator is ALL genotype columns, so missing
#                          calls sit in the denominator and are not counted as
#                          heterozygous, which is what the R version does since
#                          NA %in% het_codes is FALSE
#   rule                remove when fraction > 0.80, keep when <= 0.80
#
# What it adds. The R version printed removed SNP names to the console and
# never wrote them out, so the removal list was lost and the exact count could
# not be recovered. This writes the list, the per-SNP het fraction, and the
# per-chromosome counts, so the Methods can state a number.
#
# awk rather than R because each HapMap is 550-790 MB and there are ten of
# them. This is one streaming pass per chromosome at constant memory. Reading a
# whole HapMap with fread and then apply() over rows needs many GB and hours.
#
# Usage   bash 03_het_filter_audit.sh
#         HMP_DIR=/path/to/hapmap OUT_DIR=/path/to/out bash 03_het_filter_audit.sh
# =============================================================================
set -euo pipefail

HMP_DIR="${HMP_DIR:-/share/maize/ntanduk/SoLD/genotype/hapmap}"
OUT_DIR="${OUT_DIR:-/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/sum_ratio_lipids_final/genotype_qc}"
THRESH="${THRESH:-0.80}"

mkdir -p "$OUT_DIR"
REMOVED="$OUT_DIR/removed_high_het_SNPs_80perc.tsv"
NAMES="$OUT_DIR/removed_high_het_SNPs_80perc_names_only.txt"
SUMMARY="$OUT_DIR/het_filter_summary_by_chr.tsv"

printf 'rs\tchr\tpos\thet_fraction\tn_het\tn_genotype_cols\n' > "$REMOVED"
printf 'chr\tn_snps_total\tn_removed\tn_kept\tpct_removed\n'  > "$SUMMARY"

tot_all=0; rem_all=0
for chr in $(seq 1 10); do
  f="$HMP_DIR/SAP_only_samples_bialleles_MAF_0.05_chr${chr}.hmp.txt"
  if [[ ! -f "$f" ]]; then
    echo "WARNING: missing, skipping -> $f" >&2
    continue
  fi
  echo "chr${chr} ..." >&2

  read -r tot rem < <(
    awk -F'\t' -v thr="$THRESH" -v OFS='\t' '
      NR == 1 {
        if (NF < 12) { print "ERROR: fewer than 12 columns, not a HapMap" > "/dev/stderr"; exit 1 }
        next
      }
      {
        n = 0; het = 0
        for (i = 12; i <= NF; i++) {
          n++
          v = $i
          if (v=="Y"||v=="R"||v=="S"||v=="W"||v=="K"||v=="M") het++
        }
        total++
        if (n > 0) {
          frac = het / n
          if (frac > thr) {
            removed++
            printf "%s\t%s\t%s\t%.6f\t%d\t%d\n", $1, $3, $4, frac, het, n >> REMOVED_FILE
          }
        }
      }
      END { print total+0, removed+0 }
    ' REMOVED_FILE="$REMOVED" "$f"
  )

  kept=$(( tot - rem ))
  pct=$(awk -v r="$rem" -v t="$tot" 'BEGIN{ printf "%.4f", (t>0 ? 100*r/t : 0) }')
  printf '%s\t%s\t%s\t%s\t%s\n' "$chr" "$tot" "$rem" "$kept" "$pct" >> "$SUMMARY"
  echo "  chr${chr}: $tot SNPs, removed $rem, kept $kept (${pct}%)" >&2

  tot_all=$(( tot_all + tot )); rem_all=$(( rem_all + rem ))
done

cut -f1 "$REMOVED" | tail -n +2 > "$NAMES"

pct_all=$(awk -v r="$rem_all" -v t="$tot_all" 'BEGIN{ printf "%.4f", (t>0 ? 100*r/t : 0) }')
printf 'ALL\t%s\t%s\t%s\t%s\n' "$tot_all" "$rem_all" "$(( tot_all - rem_all ))" "$pct_all" >> "$SUMMARY"

echo
echo "total SNPs before filter : $tot_all"
echo "removed (>${THRESH} het)  : $rem_all  (${pct_all}%)"
echo "kept                     : $(( tot_all - rem_all ))"
echo
echo "removed list  : $REMOVED"
echo "names only    : $NAMES"
echo "per-chr counts: $SUMMARY"
echo
echo "The kept count is the number to use as the Bonferroni denominator, and the"
echo "removed count is the number the Methods should state."
