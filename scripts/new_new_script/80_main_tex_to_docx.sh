#!/usr/bin/env bash
# ==============================================================================
# main.tex -> Word, for editors who want .docx.
#
#   bash scripts/new_new_script/80_main_tex_to_docx.sh [out.docx]
#
# Pandoc reads the LaTeX directly. Three things have to be handled first,
# because the rilabRxiv class puts them where pandoc does not look.
#
#   1. The abstract and the keywords sit in the PREAMBLE, before
#      \begin{document}. Pandoc's LaTeX reader skips the preamble, so both were
#      dropped silently. They are moved into the body after \maketitle.
#
#   2. $^{\S}$ and $^{\dagger}$ are footnote marks written as math. Pandoc's
#      math reader rejects \S inside math and leaves the raw TeX in the Word
#      file. They are rewritten as the characters themselves.
#
#   3. Citations are natbib \citep. Pandoc resolves them with --citeproc against
#      SoLD.bib and writes a real reference list, which needs a heading of its
#      own since the LaTeX one comes from \bibliography.
#
# The figures are embedded, so the file is large. Tables come through as real
# Word tables; the booktabs rules become plain borders, which is what a copy
# editor wants anyway.
#
# WHAT DOES NOT SURVIVE, and has to be checked by hand
#   - two-column layout, page geometry and the running heads
#   - \ref cross-references to figures and tables render as the label text,
#     not as live Word fields
#   - landscape tables lose their rotation
#
# Requires pandoc 3. Prints the citation keys it could not resolve.
# ==============================================================================
set -euo pipefail

REPO="${SOLD_REPO:-.}"
OUT="${1:-$REPO/final_submission/SoLD_manuscript.docx}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

cd "$REPO"
[ -f main.tex ] || { echo "no main.tex in $REPO" >&2; exit 1; }
command -v pandoc >/dev/null || { echo "pandoc is not installed" >&2; exit 1; }

python3 - "$TMP/body.tex" <<'PY'
import re, sys
out = sys.argv[1]
s = open("main.tex").read()

abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', s, re.S)
keywords = re.search(r'\\keywords\{(.*?)\}', s, re.S)
front = ""
if abstract:
    front += "\n\\section*{Abstract}\n" + abstract.group(1).strip() + "\n"
    s = s.replace(abstract.group(0), "")
if keywords:
    front += "\n\\section*{Keywords}\n" + keywords.group(1).strip() + "\n"

s = s.replace("\\maketitle", "\\maketitle\n" + front, 1)

# ---- resolve \ref -----------------------------------------------------------
# Pandoc leaves \ref{} as the raw label, so the Word file reads
# "Supplementary [fig:suppfig_s1_workflow]". LaTeX numbers floats by the order
# their environments appear in the source, and \beginsupplement resets both
# counters and prefixes them with S, so the same walk reproduces the numbers.
body_at = s.index("\\begin{document}")
supp_at = s.index("\\beginsupplement", body_at)
num, fig, tab = {}, [0, 0], [0, 0]
for m in re.finditer(r'\\begin\{(figure\*?|table\*?)\}(.*?)\\end\{\1\}', s[body_at:], re.S):
    lab = re.search(r'\\label\{([^}]+)\}', m.group(2))
    if not lab: continue
    supp = (body_at + m.start()) > supp_at
    ctr = fig if m.group(1).startswith("figure") else tab
    ctr[supp] += 1
    kind = "Fig." if m.group(1).startswith("figure") else "Table"
    num[lab.group(1)] = ("S%d %s" % (ctr[1], "Fig" if kind == "Fig." else "Table")) if supp \
                        else str(ctr[0])
# "Supplementary~\ref{...}" already says Supplementary, so S1 Fig reads right there
s = re.sub(r'\\ref\{([^}]+)\}', lambda m: num.get(m.group(1), m.group(1)), s)
# main.tex cites one supplementary figure as Fig.~\ref{} where every other one
# uses Supplementary~\ref{}, and \thefigure already ends in "Fig", so that one
# reads "Fig. S9 Fig". It reads the same way in the compiled PDF; here the
# duplicate word is dropped so the Word file is clean either way.
s = re.sub(r'(?:Fig\.|Figure)[~ ]?(S\d+) Fig', r'Supplementary \1 Fig', s)
s = re.sub(r'Table[~ ]?(S\d+) Table', r'Supplementary \1 Table', s)
print("  cross-references resolved:", len(num))
s = s.replace("$^{\\S}$", "\\textsuperscript{§}")
s = s.replace("$^{\\dagger}$", "\\textsuperscript{†}")
s = s.replace("$^{\\ddagger}$", "\\textsuperscript{‡}")
open(out, "w").write(s)
print("  abstract moved into the body:", bool(abstract))
print("  keywords moved into the body:", bool(keywords))
PY

pandoc "$TMP/body.tex" \
  --from=latex --to=docx \
  --citeproc --bibliography=SoLD.bib \
  --metadata reference-section-title="References" \
  --resource-path=".:fig" \
  -o "$OUT" 2> "$TMP/warn.txt" || { cat "$TMP/warn.txt" >&2; exit 1; }

grep -o "citation [A-Za-z0-9_]* not found" "$TMP/warn.txt" | sort -u | sed 's/^/  UNRESOLVED /' || true
grep -v "citation .* not found" "$TMP/warn.txt" | grep -v '^\s*$' | head -5 || true

python3 - "$OUT" <<'PY'
import sys, zipfile, re
z = zipfile.ZipFile(sys.argv[1])
d = z.read("word/document.xml").decode("utf8")
print("\n  %s" % sys.argv[1])
print("  images %d   tables %d   paragraphs %d"
      % (len([n for n in z.namelist() if n.startswith("word/media/")]),
         d.count("<w:tbl>"), d.count("<w:p ") + d.count("<w:p>")))
for probe in ("Sorghum Lipid Database", "Tandukar", "Abstract", "Keywords", "References"):
    print("    %-24s %s" % (probe, "present" if probe in d else "ABSENT"))
PY
