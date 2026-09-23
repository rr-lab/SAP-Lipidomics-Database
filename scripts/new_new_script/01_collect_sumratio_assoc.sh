#!/usr/bin/env bash
# =============================================================================
# Collect the GEMMA .assoc.txt files from the vcf2gwas output tree into the
# project results folder, split by trial.
#
# vcf2gwas writes one directory per trait, and inside it one timestamped
# directory per run. Both trials land under the same trait directory, so the
# trial is read from the FILENAME (Final_control_ / Final_lowinput_), never
# from the timestamp, which is not a reliable key.
#
# Destination layout matches individual_lipids_final
#     <DEST>/control/raw_gwas/   and   <DEST>/lowinput/raw_gwas/
# with the original vcf2gwas filenames kept, since label_from_file() parses the
# trait off the front of the name. The sibling annotation/ manhattan/ and
# table_annotation/ folders are created empty so the two trees match.
#
# Usage   bash 01_collect_sumratio_assoc.sh
#         MODE=move bash 01_collect_sumratio_assoc.sh
# =============================================================================
set -euo pipefail

SRC="${SRC:-/share/maize/ntanduk/SoLD/ratios/Output/Linear_Mixed_Model}"
DEST="${DEST:-/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/sum_ratio_lipids_final}"
MODE="${MODE:-copy}"          # copy (default) or move
EXPECTED_PER_TRIAL="${EXPECTED_PER_TRIAL:-153}"

[[ -d "$SRC" ]] || { echo "ERROR: source not found: $SRC" >&2; exit 1; }
for t in control lowinput; do
  mkdir -p "$DEST/$t/raw_gwas" "$DEST/$t/annotation" "$DEST/$t/manhattan" "$DEST/$t/table_annotation"
done
MANIFEST="$DEST/assoc_manifest.tsv"
printf 'trial\ttrait\tsrc_path\tdest_path\tbytes\tn_snps\n' > "$MANIFEST"

n_ctl=0; n_lin=0; n_skip=0

# -L so a symlinked output tree still works. Only .assoc.txt, nothing else.
while IFS= read -r f; do
  b=$(basename "$f")

  case "$b" in
    *Final_control_*)  trial=CTL; sub=control/raw_gwas  ;;
    *Final_lowinput_*) trial=LIN; sub=lowinput/raw_gwas ;;
    *) echo "SKIP (no trial token): $f" >&2; n_skip=$((n_skip+1)); continue ;;
  esac

  # trait name is everything before _mod_sub_, the same rule label_from_file() uses
  trait="${b%%_mod_sub_*}"
  [[ "$trait" != "$b" ]] || { echo "SKIP (no _mod_sub_): $f" >&2; n_skip=$((n_skip+1)); continue; }

  out="$DEST/$sub/$b"
  if [[ -e "$out" ]]; then
    echo "ERROR: destination already exists, refusing to overwrite: $out" >&2
    exit 1
  fi

  if [[ "$MODE" == "move" ]]; then mv "$f" "$out"; else cp -p "$f" "$out"; fi

  nsnp=$(( $(wc -l < "$out") - 1 ))
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$trial" "$trait" "$f" "$out" "$(stat -c%s "$out")" "$nsnp" >> "$MANIFEST"
  if [[ "$trial" == CTL ]]; then n_ctl=$((n_ctl+1)); else n_lin=$((n_lin+1)); fi
done < <(find -L "$SRC" -mindepth 2 -type f -name '*.assoc.txt' | sort)

echo
echo "MODE      : $MODE"
echo "CTL files : $n_ctl   -> $DEST/control/raw_gwas"
echo "LIN files : $n_lin   -> $DEST/lowinput/raw_gwas"
echo "skipped   : $n_skip"
echo "manifest  : $MANIFEST"

# ---- checks -----------------------------------------------------------------
fail=0
for t in CTL LIN; do
  n=$([[ $t == CTL ]] && echo $n_ctl || echo $n_lin)
  if [[ "$n" -ne "$EXPECTED_PER_TRIAL" ]]; then
    echo "WARNING: $t has $n files, expected $EXPECTED_PER_TRIAL" >&2; fail=1
  fi
done

# every trait should appear exactly once per trial
dup=$(cut -f1,2 "$MANIFEST" | tail -n +2 | sort | uniq -d)
if [[ -n "$dup" ]]; then
  echo "WARNING: trait appears more than once within a trial:" >&2
  echo "$dup" >&2; fail=1
fi

# traits present in one trial only
only=$(awk -F'\t' 'NR>1{c[$2]=c[$2]$1}END{for(k in c) if(length(c[k])!=6) print k" -> "c[k]}' "$MANIFEST")
if [[ -n "$only" ]]; then
  echo "WARNING: traits not present in both trials:" >&2
  echo "$only" >&2; fail=1
fi

# all assoc files should have the same SNP count within a trial
awk -F'\t' 'NR>1{print $1"\t"$6}' "$MANIFEST" | sort -u | awk -F'\t' '{c[$1]++} END{for(k in c) if(c[k]>1) print "WARNING: "k" has "c[k]" different SNP counts across traits"}' >&2

[[ $fail -eq 0 ]] && echo "checks passed" || echo "checks raised warnings, see above"
