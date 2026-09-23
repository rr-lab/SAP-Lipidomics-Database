#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J ldannot
#BSUB -o stdout.%J
#BSUB -e stderr.%J
# =============================================================================
# LD-aware gene annotation of the lipid GWAS, both layers, both trials.
#
#   raw_gwas/*.assoc.txt
#     -> Bonferroni-significant SNPs                 (trait_snp.tsv)
#     -> PLINK ids via the .bim, joined on chr:bp    (sig_snp_map.tsv)
#     -> strong-LD partners, r2 >= 0.4 within 250 kb (sig_strongLD.ld)
#     -> genes whose body contains a partner         (snp_gene_LD.tsv)
#
# Then run 12_make_annotation_tables.R for the four output tables.
#
# Chromosome naming. Every source spells the chromosome differently. The
# HapMap-derived het list uses 01..10, the assoc files use 1..10, and a .bim
# may use 1, 01 or Chr01. nc() strips any Chr prefix and coerces to a number,
# and every key in this script goes through it. Without that the joins fail
# silently and the analysis quietly runs unfiltered.
#
# LD is computed per chromosome when the per-chromosome bfiles are present,
# because r2 is only ever needed inside a 250 kb window, so a whole-genome
# --r2 run wastes memory for no extra pairs.
#
# Usage   bsub < 11_LD_annotate_genes.sh
#         BED_DIR=... LAYERS=sum_ratio_lipids_final bash 11_LD_annotate_genes.sh
# =============================================================================
set -euo pipefail

BLUP="${BLUP:-/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP}"
OUT_DIR="${OUT_DIR:-$BLUP/gene_annotation_final}"
WORK="${WORK:-$OUT_DIR/work}"
PLINK="${PLINK:-/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/plink}"
BED_DIR="${BED_DIR:-/rsstu/users/r/rrellan/DOE_CAREER/SorghumGEA/data/SAP/bed_files}"
BFILE="${BFILE:-$BED_DIR/SAP_ldpanel}"                 # genome-wide set
CHR_PREFIX="${CHR_PREFIX:-$BED_DIR/SAP_ldpanel_Chr}"   # per-chromosome sets
GENES="${GENES:-$OUT_DIR/genes.range}"
HETLIST="${HETLIST:-$BLUP/sum_ratio_lipids_final/genotype_qc/removed_high_het_SNPs_80perc.tsv}"
LAYERS="${LAYERS:-individual_lipids_final sum_ratio_lipids_final}"
R2="${R2:-0.4}"
KB="${KB:-250}"
BS=100000

mkdir -p "$WORK" "$OUT_DIR"
[[ -f "$GENES" ]] || { echo "ERROR: run 10_build_genes_range.sh first, missing $GENES" >&2; exit 1; }
[[ -x "$PLINK" ]] || { echo "ERROR: plink not executable at $PLINK" >&2; exit 1; }

NC='function nc(x){ sub(/^[Cc]hr_?/,"",x); return x+0 }'

# ---- 0. which bfiles exist ---------------------------------------------------
CHRS=()
for c in 01 02 03 04 05 06 07 08 09 10; do
  [[ -f "${CHR_PREFIX}${c}.bim" ]] && CHRS+=("$c")
done
BIM_ALL="$WORK/all.bim"
if [[ ${#CHRS[@]} -gt 0 ]]; then
  echo "per-chromosome bfiles found: ${#CHRS[@]}  (${CHR_PREFIX}NN)"
  cat "${CHR_PREFIX}"{01,02,03,04,05,06,07,08,09,10}.bim 2>/dev/null > "$BIM_ALL" || true
  [[ -s "$BIM_ALL" ]] || { for c in "${CHRS[@]}"; do cat "${CHR_PREFIX}${c}.bim"; done > "$BIM_ALL"; }
elif [[ -f "$BFILE.bim" ]]; then
  echo "using the genome-wide bfile: $BFILE"
  cp "$BFILE.bim" "$BIM_ALL"
else
  echo "ERROR: no bfile found. Looked for ${CHR_PREFIX}NN.bim and $BFILE.bim" >&2
  echo "       set BED_DIR, or BFILE / CHR_PREFIX directly" >&2
  exit 1
fi
echo "markers in the LD panel: $(wc -l < "$BIM_ALL")"
echo "first .bim row         : $(head -1 "$BIM_ALL")"
echo "derived chromosome key : $(awk -v FS='[ \t]+' "$NC"'NR==1{print nc($1)}' "$BIM_ALL")"

# ---- 1. high-heterozygosity markers to exclude --------------------------------
HET="$WORK/het_exclude_chrbp.txt"
HETID="$WORK/het_exclude_plink_ids.txt"
if [[ -f "$HETLIST" ]]; then
  awk -v FS='[ \t]+' "$NC"'NR>1{ print nc($2)":"($3+0) }' "$HETLIST" | sort -u > "$HET"
  awk -v FS='[ \t]+' "$NC"'
    FNR==NR { bad[$0]; next }
    (nc($1)":"($4+0)) in bad { print $2 }
  ' "$HET" "$BIM_ALL" | sort -u > "$HETID"
  n_het=$(wc -l < "$HET"); n_hetid=$(wc -l < "$HETID")
  echo "high-het markers on the removal list : $n_het"
  echo "                 present in the panel : $n_hetid"
  if [[ "$n_het" -gt 0 && "$n_hetid" -eq 0 ]]; then
    echo "ERROR: none of the $n_het high-het markers matched the LD panel on chr:pos." >&2
    echo "       removal list key : $(head -1 "$HET")" >&2
    echo "       panel key        : $(awk -v FS='[ \t]+' "$NC"'NR==1{print nc($1)":"($4+0)}' "$BIM_ALL")" >&2
    echo "       Refusing to run an unfiltered analysis." >&2
    exit 1
  fi
else
  : > "$HET"; : > "$HETID"
  echo "NOTE: no high-het removal list at $HETLIST, proceeding without that filter" >&2
fi

# ---- 2. Bonferroni-significant SNPs ------------------------------------------
TS="$WORK/trait_snp.tsv"
printf 'Layer\tCondition\tTrait\tChr\tBP\tP\n' > "$TS"
for layer in $LAYERS; do
  for trial in control lowinput; do
    cond=$([[ $trial == control ]] && echo CTL || echo LIN)
    dir="$BLUP/$layer/$trial/raw_gwas"
    [[ -d "$dir" ]] || { echo "skip, no such dir: $dir" >&2; continue; }
    mapfile -t files < <(find "$dir" -maxdepth 1 -type f -name '*.assoc.txt' | sort)
    [[ ${#files[@]} -gt 0 ]] || { echo "skip, no assoc files in $dir" >&2; continue; }

    ntest=$(awk "$NC"'FNR==NR{bad[$0]; next} FNR>1 && !((nc($1)":"($3+0)) in bad)' "$HET" "${files[0]}" | wc -l)
    [[ "$ntest" -gt 0 ]] || { echo "ERROR: no markers left after the het filter in $dir" >&2; exit 1; }
    thr=$(awk -v n="$ntest" 'BEGIN{ printf "%.6e", 0.05/n }')
    echo "$layer / $cond : ${#files[@]} traits, $ntest markers tested, Bonferroni $thr"

    for f in "${files[@]}"; do
      trait=$(basename "$f"); trait="${trait%%_mod_sub_*}"; trait="${trait%.assoc.txt}"
      awk -v L="$layer" -v C="$cond" -v T="$trait" -v thr="$thr" -v OFS='\t' "$NC"'
        NR==1 { for(i=1;i<=NF;i++) h[$i]=i
                if(!("p_wald" in h)){ print "ERROR: no p_wald column" > "/dev/stderr"; exit 1 }
                pc=h["p_wald"]; cc=h["chr"]; bc=h["ps"]; next }
        $pc+0 <= thr { print L, C, T, nc($cc), $bc+0, $pc }
      ' "$f" >> "$TS"
    done
  done
done
if [[ -s "$HET" ]]; then
  head -1 "$TS" > "$TS.tmp"
  awk -F'\t' 'NR==FNR{bad[$0]; next} FNR>1 && !((($4+0)":"($5+0)) in bad)' "$HET" "$TS" >> "$TS.tmp"
  mv "$TS.tmp" "$TS"
fi
echo "significant trait-SNP rows: $(( $(wc -l < "$TS") - 1 ))"

# ---- 3. resolve PLINK ids on chr:bp ------------------------------------------
MAP="$WORK/sig_snp_map.tsv"
awk -F'\t' 'NR>1{ print ($4+0)"\t"($5+0) }' "$TS" | sort -u > "$WORK/sig_chrbp.tsv"
awk -v FS='[ \t]+' -v OFS='\t' "$NC"'
  FNR==NR { want[($1+0)":"($2+0)]; next }
  (nc($1)":"($4+0)) in want { print $2, nc($1), $4+0 }
' "$WORK/sig_chrbp.tsv" "$BIM_ALL" | sort -u > "$MAP"

n_want=$(wc -l < "$WORK/sig_chrbp.tsv"); n_got=$(wc -l < "$MAP")
echo "distinct significant positions: $n_want   present in the LD panel: $n_got"
[[ "$n_got" -gt 0 ]] || { echo "ERROR: no significant position is in the LD panel. Wrong bfile." >&2; exit 1; }
[[ "$n_got" -lt "$n_want" ]] && echo "NOTE: $(( n_want - n_got )) significant positions are not in the LD panel" >&2
cut -f1 "$MAP" | sort -u > "$WORK/sig_snps.txt"

# ---- 4. strong-LD partners ----------------------------------------------------
LD="$WORK/sig_strongLD.ld"
: > "$LD"
run_plink () {  # $1 bfile prefix   $2 snp list   $3 out prefix
  local excl=()
  [[ -s "$HETID" ]] && excl=(--exclude "$HETID")
  "$PLINK" --bfile "$1" --allow-extra-chr --r2 \
    "${excl[@]}" --ld-snp-list "$2" \
    --ld-window-kb "$KB" --ld-window 99999 --ld-window-r2 "$R2" \
    --out "$3" >> "$WORK/plink.stdout" 2>&1
}
: > "$WORK/plink.stdout"
echo "running plink --r2, r2 >= $R2, window ${KB} kb${HETID:+, excluding $(wc -l < "$HETID") high-het markers}"

if [[ ${#CHRS[@]} -gt 0 ]]; then
  for c in "${CHRS[@]}"; do
    cnum=$((10#$c))
    awk -F'\t' -v c="$cnum" '$2+0==c{print $1}' "$MAP" > "$WORK/sig_chr${c}.txt"
    if [[ ! -s "$WORK/sig_chr${c}.txt" ]]; then echo "  chr${c}: no significant SNPs, skipped"; continue; fi
    if run_plink "${CHR_PREFIX}${c}" "$WORK/sig_chr${c}.txt" "$WORK/ld_chr${c}"; then
      tail -n +2 "$WORK/ld_chr${c}.ld" >> "$LD"
      echo "  chr${c}: $(wc -l < "$WORK/sig_chr${c}.txt") lead SNPs, $(( $(wc -l < "$WORK/ld_chr${c}.ld") - 1 )) LD pairs"
    else
      echo "  chr${c}: plink failed, see $WORK/plink.stdout and $WORK/ld_chr${c}.log" >&2; exit 1
    fi
  done
else
  run_plink "$BFILE" "$WORK/sig_snps.txt" "$WORK/ld_all" || {
    echo "ERROR: plink failed, see $WORK/plink.stdout" >&2; exit 1; }
  tail -n +2 "$WORK/ld_all.ld" >> "$LD"
fi
echo "LD pairs total: $(wc -l < "$LD")"

# ---- 5. genes containing a partner, plus each lead SNP's own gene ------------
SG="$WORK/snp_gene_LD.tsv"
awk -v BS="$BS" -v OFS='\t' "$NC"'
  FNR==NR { c=nc($1); s=$2+0; e=$3+0; g=$4
            for (b=int(s/BS); b<=int(e/BS); b++){ k=c":"b; n[k]++; G[k,n[k]]=g; S[k,n[k]]=s; E[k,n[k]]=e }
            next }
  $1=="CHR_A" { next }
  { lead=$3; c=nc($4); bp=$5+0; r2=$7
    k=c":"int(bp/BS); m=n[k]
    for (i=1;i<=m;i++) if (S[k,i]<=bp && E[k,i]>=bp) print lead, G[k,i], r2 }
' "$GENES" "$LD" > "$SG"

awk -v BS="$BS" -v OFS='\t' "$NC"'
  FNR==NR { c=nc($1); s=$2+0; e=$3+0; g=$4
            for (b=int(s/BS); b<=int(e/BS); b++){ k=c":"b; n[k]++; G[k,n[k]]=g; S[k,n[k]]=s; E[k,n[k]]=e }
            next }
  { id=$1; c=$2+0; bp=$3+0
    k=c":"int(bp/BS); m=n[k]
    for (i=1;i<=m;i++) if (S[k,i]<=bp && E[k,i]>=bp) print id, G[k,i], "1" }
' "$GENES" "$MAP" >> "$SG"

sort -u "$SG" -o "$SG"
echo "SNP-gene pairs: $(wc -l < "$SG")"
echo
echo "work files in $WORK"
echo "now run  Rscript 12_make_annotation_tables.R"
