#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32GB]"
#BSUB -J grmfilt
#BSUB -o stdout.%J
#BSUB -e stderr.%J
# =============================================================================
# Rebuild the genomic relationship matrix with the high-heterozygosity markers
# excluded, as a drop-in replacement for data/kinship/sap_grm.rel.
#
# The shipped GRM predates the heterozygosity filter, so it was computed on all
# markers including the 83,640 removed at >80% heterozygous calls. The same GRM
# feeds the heritability estimates, so both inherit that.
#
# Expect a small change. 83,640 of roughly 9 million markers is under 1%, and a
# relationship matrix is an average over markers, so the two GRMs should agree
# closely. The point is to know that rather than assume it, and to be able to
# say so in the Methods. 21_heritability_compare_grm.R measures the difference.
#
# Output is written in exactly the format the heritability scripts read:
# a square .rel with a matching .rel.id, so nothing downstream changes except
# the path.
#
# Usage   bsub < 20_build_filtered_grm.sh
# =============================================================================
set -euo pipefail

PLINK="${PLINK:-/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/plink}"
BED_DIR="${BED_DIR:-/rsstu/users/r/rrellan/DOE_CAREER/SorghumGEA/data/SAP/bed_files}"
BFILE="${BFILE:-$BED_DIR/SAP_ldpanel}"
CHR_PREFIX="${CHR_PREFIX:-$BED_DIR/SAP_ldpanel_Chr}"
HETLIST="${HETLIST:-/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/sum_ratio_lipids_final/genotype_qc/removed_high_het_SNPs_80perc.tsv}"
OUT_DIR="${OUT_DIR:-/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/kinship_hetfiltered}"
WORK="$OUT_DIR/work"
mkdir -p "$WORK"

[[ -x "$PLINK"    ]] || { echo "ERROR: plink not executable at $PLINK" >&2; exit 1; }
[[ -f "$HETLIST"  ]] || { echo "ERROR: het removal list not found at $HETLIST" >&2; exit 1; }

# ---- 1. assemble one genome-wide bfile --------------------------------------
# The LD panel ships per chromosome. A GRM must be computed across all markers
# at once, not per chromosome, so merge first if the whole-genome file is absent.
if [[ -f "$BFILE.bed" ]]; then
  GENO="$BFILE"
  echo "using genome-wide bfile: $GENO"
else
  echo "merging per-chromosome bfiles ..."
  : > "$WORK/merge_list.txt"
  for c in 02 03 04 05 06 07 08 09 10; do
    [[ -f "${CHR_PREFIX}${c}.bed" ]] && echo "${CHR_PREFIX}${c}" >> "$WORK/merge_list.txt"
  done
  "$PLINK" --bfile "${CHR_PREFIX}01" --merge-list "$WORK/merge_list.txt" \
    --allow-extra-chr --make-bed --out "$WORK/merged" > "$WORK/merge.stdout" 2>&1
  GENO="$WORK/merged"
fi
echo "markers available: $(wc -l < "$GENO.bim")"

# ---- 2. resolve the high-het markers to this bfile's ids --------------------
# Keyed on chromosome and position, because the HapMap the list came from spells
# chromosomes 01..10 while a .bim may use 1..10 or Chr01. nc() normalises both.
NC='function nc(x){ sub(/^[Cc]hr_?/,"",x); return x+0 }'
awk "$NC"'NR>1{ print nc($2)":"($3+0) }' "$HETLIST" | sort -u > "$WORK/het_chrbp.txt"
awk -v FS='[ \t]+' "$NC"'
  FNR==NR { bad[$0]; next }
  (nc($1)":"($4+0)) in bad { print $2 }
' "$WORK/het_chrbp.txt" "$GENO.bim" | sort -u > "$WORK/het_exclude_ids.txt"

n_list=$(wc -l < "$WORK/het_chrbp.txt"); n_hit=$(wc -l < "$WORK/het_exclude_ids.txt")
echo "high-het markers on the list : $n_list"
echo "             found in the bfile: $n_hit"
if [[ "$n_hit" -eq 0 ]]; then
  echo "ERROR: none matched on chr:pos. Refusing to build an unchanged GRM and call it filtered." >&2
  echo "  list key : $(head -1 "$WORK/het_chrbp.txt")" >&2
  echo "  bfile key: $(awk -v FS='[ \t]+' "$NC"'NR==1{print nc($1)":"($4+0)}' "$GENO.bim")" >&2
  exit 1
fi

# ---- 3. build both GRMs, so the comparison is like for like ------------------
# Same samples, same options, only the marker set differs.
"$PLINK" --bfile "$GENO" --allow-extra-chr --make-rel square \
  --out "$OUT_DIR/sap_grm_allmarkers" > "$WORK/grm_all.stdout" 2>&1
"$PLINK" --bfile "$GENO" --allow-extra-chr --exclude "$WORK/het_exclude_ids.txt" \
  --make-rel square --out "$OUT_DIR/sap_grm_hetfiltered" > "$WORK/grm_filt.stdout" 2>&1

for tag in allmarkers hetfiltered; do
  f="$OUT_DIR/sap_grm_${tag}.rel"
  [[ -s "$f" ]] || { echo "ERROR: plink produced no $f, see $WORK/grm_${tag}.stdout" >&2; exit 1; }
  echo "  sap_grm_${tag}.rel : $(wc -l < "$f") individuals"
done

echo
echo "markers used, all       : $(wc -l < "$GENO.bim")"
echo "markers used, filtered  : $(( $(wc -l < "$GENO.bim") - n_hit ))"
echo "output: $OUT_DIR"
echo
echo "Next: Rscript 21_heritability_compare_grm.R"
