#!/usr/bin/env python3
"""
main.tex -> a single-column Word manuscript, body through the reference list.

  python3 build_docx.py  [out.docx]

WHAT IT DOES THAT PANDOC ALONE DOES NOT

  1. FLOATS ARE MOVED TO THE TEXT. LaTeX floats sit wherever the page breaks
     allow. Each main figure and table is lifted out of the source and
     re-inserted directly after the paragraph that first cites it, which is
     what a manuscript submission wants.

  2. NUMBERS ARE BAKED IN. Pandoc leaves \ref{} as the raw label. Floats are
     numbered in source order and every \ref is replaced by its number, so the
     Word file reads "Fig. 3" rather than "Fig. [fig:gwas_manhattan]". The
     captions get a "Figure N." / "Table N." prefix for the same reason.

  3. THE SUPPLEMENT IS DROPPED BUT ITS POINTERS SURVIVE. Everything from
     \section*{Supplement} on is cut. The supplementary figures are numbered
     first, so "Supplementary~\ref{...}" still resolves to "Supplementary S4
     Fig" instead of a dangling label.

  4. FRONT MATTER IS REBUILT. The rilabRxiv class puts the title, authors,
     affiliations, abstract and keywords in the preamble, which pandoc's reader
     skips. They are rewritten into the body with the affiliation superscripts
     preserved.

  5. MATH-MODE FOOTNOTE MARKS ARE DEMOTED. $^{\S}$ and friends break pandoc's
     math reader and come through as raw TeX; they are rewritten as characters.

Citations are resolved by citeproc against SoLD.bib. Images are downsampled to
a sensible print width so the file is not 200 MB.
"""
import re, os, shutil, subprocess, sys
import importlib.util as _il, os as _os
_sp = _il.spec_from_file_location("mathfix", _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "82a_mathfix.py"))
mathfix = _il.module_from_spec(_sp); _sp.loader.exec_module(mathfix)
from PIL import Image

SRC, BIB, WORK = "main.tex", "SoLD.bib", "work"
OUT = sys.argv[1] if len(sys.argv) > 1 else "SoLD_manuscript.docx"
IMG_W = 1950                       # px; 6.5 in at 300 dpi
os.makedirs(WORK, exist_ok=True)

s = open(SRC, encoding="utf8").read()

# ---------------------------------------------------------------- front matter
def braced(txt, i):
    """Return the content of the {...} group starting at or after index i."""
    i = txt.index("{", i); d, k = 1, i + 1
    while d:
        if txt[k] == "{": d += 1
        elif txt[k] == "}": d -= 1
        k += 1
    return txt[i+1:k-1], k

title = re.search(r'\\title\{', s)
title = braced(s, title.start())[0] if title else "Untitled"
authors = [(braced(s, m.end()-1)[0], m.group(1)) for m in
           re.finditer(r'\\author\[([^\]]*)\]\s*(?=\{)', s)]
affils  = [(m.group(1), braced(s, m.end()-1)[0]) for m in
           re.finditer(r'\\affil\[([^\]]*)\]\s*(?=\{)', s)]
abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', s, re.S)
abstract = abstract.group(1).strip() if abstract else ""
kw = re.search(r'\\keywords\{', s)
kw = braced(s, kw.start())[0].strip() if kw else ""

def marks(spec):
    """[$1$,$2$,*] -> 1,2,*"""
    return ",".join(x.strip().strip("$") for x in spec.split(",") if x.strip())

front = ["\\section*{%s}" % title, ""]
front.append(", ".join("%s\\textsuperscript{%s}" % (n, marks(a)) for n, a in authors))
front.append("")
for a, text in affils:
    front.append("\\textsuperscript{%s}%s\n" % (marks(a), text))
front += ["", "\\subsection*{Abstract}", abstract, ""]
if kw:
    front += ["\\subsection*{Keywords}", kw, ""]
front = "\n".join(front)

# ---------------------------------------------------------------------- body
body_at = s.index("\\begin{document}")
cut = min(i for i in (s.find("\\section*{Supplement}", body_at),
                      s.find("\\beginsupplement", body_at + 400)) if i > 0)
supp_src, body = s[cut:], s[body_at:cut]

# supplementary figures are numbered from the part we are cutting, so that
# "Supplementary~\ref{...}" in the retained text still resolves
supp_num, n = {}, 0
for m in re.finditer(r'\\begin\{(figure\*?)\}(.*?)\\end\{\1\}', supp_src, re.S):
    lab = re.search(r'\\label\{([^}]+)\}', m.group(2))
    if lab: n += 1; supp_num[lab.group(1)] = "S%d Fig" % n

# ------------------------------------------------- pull the floats out of body
def caption_of(blk):
    m = re.search(r'\\caption', blk)
    if not m: return None, None, None
    i = m.end()
    if blk[i] == "[":                       # \caption[short]{long}
        d, k = 1, i + 1
        while d:
            if blk[k] == "[": d += 1
            elif blk[k] == "]": d -= 1
            k += 1
        i = k
    txt, end = braced(blk, i)
    return txt, m.start(), end

floats, fig_n, tab_n = [], 0, 0
for m in re.finditer(r'\\begin\{(figure\*?|table\*?)\}(?:\[[^\]]*\])?(.*?)\\end\{\1\}', body, re.S):
    blk = m.group(2)
    isfig = m.group(1).startswith("figure")
    lab = re.search(r'\\label\{([^}]+)\}', blk)
    cap, _, _ = caption_of(blk)
    if isfig: fig_n += 1; num, kind = fig_n, "Figure"
    else:     tab_n += 1; num, kind = tab_n, "Table"
    floats.append(dict(raw=m.group(0), inner=blk, isfig=isfig, num=num, kind=kind,
                       label=lab.group(1) if lab else None, cap=cap, span=m.span()))

num_of = {f["label"]: str(f["num"]) for f in floats if f["label"]}
for f in floats: body = body.replace(f["raw"], "\n\n@@FLOAT%d%s@@\n\n" % (f["num"], "F" if f["isfig"] else "T"))

# ------------------------------------------------------------- resolve \ref{}
def ref_sub(m):
    k = m.group(1)
    if k in num_of:  return num_of[k]
    if k in supp_num: return supp_num[k]
    return k
body = re.sub(r'\\ref\{([^}]+)\}', ref_sub, body)
# captions were lifted out before this point, so they need the same pass
for f in floats:
    if f["cap"]:
        f["cap"] = re.sub(r'\\ref\{([^}]+)\}', ref_sub, f["cap"])
        f["cap"] = mathfix.apply(f["cap"])[0]
# "Supplementary~S9 Fig" reads fine; "Fig.~S9 Fig" would not, but that was fixed upstream
body = re.sub(r'(?:Fig\.|Figure)[~ ]?(S\d+) Fig', r'Supplementary \1 Fig', body)

# -------------------------------------------------------- rebuild float blocks
def shrink(path):
    dst = os.path.join(WORK, os.path.basename(path))
    im = Image.open(path)
    if im.width > IMG_W:
        im = im.resize((IMG_W, round(im.height * IMG_W / im.width)), Image.LANCZOS)
    im.save(dst, dpi=(300, 300))
    return dst

def render(f):
    cap = "\\textbf{%s %d.} %s" % (f["kind"], f["num"], f["cap"] or "")
    if f["isfig"]:
        g = re.search(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', f["inner"])
        return "\\begin{figure}\n\\includegraphics{%s}\n\\caption{%s}\n\\end{figure}" % (
                shrink(g.group(1)), cap)
    inner = f["inner"]
    # \cmidrule and \addlinespace have no Word equivalent and pandoc emits
    # their arguments as cell text ("3-4(lr)5-6"), so they go before conversion
    inner = re.sub(r'\\cmidrule(?:\([a-z]+\))?\{[^}]*\}', '', inner)
    inner = inner.replace("\\addlinespace", "")
    inner = mathfix.apply(inner)[0]
    i = inner.index("\\begin{tabular}")
    j = inner.index("\\end{tabular}") + len("\\end{tabular}")
    return "\\begin{table}\n\\caption{%s}\n%s\n\\end{table}" % (cap, inner[i:j])

# ---------------------------------------- move each float after its first cite
# Remove every placeholder FIRST. Removing them one at a time while also
# recording insertion indices shifts the list underneath the indices already
# recorded, which silently drops floats.
paras = [p for p in re.split(r'\n\s*\n', body)
         if not re.fullmatch(r'@@FLOAT\d+[FT]@@', p.strip())]

place, unplaced = {}, []
for f in floats:
    # a panel letter may follow the number ("Fig.~2A"), so the guard is
    # "not another digit" rather than a word boundary
    pat = (r'(?:Fig\.|Figure)[~ ]*%d(?![0-9])' if f["isfig"]
           else r'Table[~ ]*%d(?![0-9])') % f["num"]
    hit = next((i for i, p in enumerate(paras) if re.search(pat, p)), None)
    if hit is None: unplaced.append(f); hit = len(paras) - 1
    place.setdefault(hit, []).append(f)

out = []
for idx, p in enumerate(paras):
    out.append(p)
    for f in sorted(place.get(idx, []), key=lambda x: (not x["isfig"], x["num"])):
        out.append(render(f))
body = "\n\n".join(out)
if unplaced:
    print("  NOT CITED, appended at the end:",
          ", ".join("%s %d" % (f["kind"], f["num"]) for f in unplaced))

# --------------------------------------------------------------- tidy up TeX
# inline math -> real Word runs (see mathfix for why)
body, n = mathfix.apply(body)
print("inline math rewritten as text: %d  |  left as equations: %d" % tuple(n))
body = body.replace("$^{\\S}$", "\\textsuperscript{\\S}")
body = body.replace("$^{\\dagger}$", "\\textsuperscript{\\dag}")
body = body.replace("$^{\\ddagger}$", "\\textsuperscript{\\ddag}")
for junk in (r'\\begin\{document\}', r'\\maketitle', r'\\bibliography\{[^}]*\}',
             r'\\onecolumn', r'\\twocolumn', r'\\pagebreak', r'\\clearpage',
             r'\\setboolean\{[^}]*\}\{[^}]*\}', r'\\linenumbers', r'\\thispagestyle\{[^}]*\}'):
    body = re.sub(junk, "", body)
body = re.sub(r'\n{3,}', '\n\n', body)

open(os.path.join(WORK, "body.tex"), "w", encoding="utf8").write(front + "\n\n" + body)

# ------------------------------------------------------------------- pandoc
cmd = ["pandoc", os.path.join(WORK, "body.tex"), "--from=latex", "--to=docx",
       "--citeproc", "--bibliography=" + BIB,
       "--metadata", "reference-section-title=References",
       "--resource-path=.:fig:" + WORK, "-o", OUT]
r = subprocess.run(cmd, capture_output=True, text=True)
for line in r.stderr.splitlines():
    if "not found" in line or "Error" in line: print("  pandoc:", line)
if r.returncode: sys.exit("pandoc failed:\n" + r.stderr)

print("floats placed:")
for f in floats:
    print("   %-8s %-2d  %s" % (f["kind"], f["num"], (f["cap"] or "")[:58].replace("\n"," ")))
print("\nwrote", OUT, "%.1f MB" % (os.path.getsize(OUT)/1e6))
