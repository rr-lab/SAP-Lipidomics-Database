"""
Inline math -> real Word text runs.

Pandoc turns every $...$ into an OMML equation object. Word renders those, but
they are brittle: converters drop them inside table cells, so a table of
q-values comes out with blank columns, and a manuscript full of equation
objects is awkward for a copy editor.

Almost all the inline math here is a variable, a relation and a number. Those
are rewritten as ordinary text with \textsuperscript / \textsubscript, which
pandoc turns into normal runs with superscript formatting. Anything with real
structure (\frac, \binom, \bar, a display equation) is left as math so nothing
is silently mangled.
"""
import re

SYM = {r'\times':'×', r'\leq':'≤', r'\le':'≤', r'\geq':'≥', r'\ge':'≥',
       r'\sim':'\\textasciitilde{}', r'\rightarrow':'→', r'\leftrightarrow':'↔', r'\approx':'≈',
       r'\varepsilon':'ε', r'\epsilon':'ε', r'\beta':'β', r'\alpha':'α',
       r'\mu':'μ', r'\sigma':'σ', r'\ell':'ℓ', r'\cap':'∩', r'\cup':'∪',
       r'\pm':'±', r'\S':'§', r'\dagger':'†', r'\ddagger':'‡',
       r'\ln':'ln', r'\log':'log', r'\max':'max', r'\min':'min',
       r'\,':'\u2009', r'\;':' ', r'\ ':' ', r'\%':'\\%'}
KEEP = (r'\frac', r'\binom', r'\bar', r'\sqrt', r'\sum', r'\int', r'\begin', r'\mathcal')
# single letters that are variables and get italics; everything else stays upright
VARS = set("abcefghijknpqrstuvwxyzABCDEFGHIJKLMNPQRSTUVXYZ")

def _grp(s, i):
    """(content, next index) for the {...} at i, or the single next character"""
    if i < len(s) and s[i] == '{':
        d, k = 1, i + 1
        while d < 100 and k < len(s):
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
            if d == 0: break
        return s[i+1:k-1], k
    return (s[i], i + 1) if i < len(s) else ('', i)

def convert(m):
    """LaTeX inline math -> text-mode LaTeX, or None if it should stay math."""
    if any(k in m for k in KEEP): return None
    m = m.replace('{>}', '>').replace('{<}', '<').replace('~', ' ')
    m = re.sub(r'\\math(?:rm|bf|it|sf)\{([^{}]*)\}', r'\1', m)
    for k in sorted(SYM, key=len, reverse=True):
        m = m.replace(k, SYM[k])
    # what is left may only be the escapes this function itself emitted
    if re.sub(r'\\\\(?:textasciitilde\\{\\}|%)', '', m).count('\\\\'):
        return None                                # an unknown macro survived

    out, i = [], 0
    while i < len(m):
        ch = m[i]
        if ch in '^_':
            g, i = _grp(m, i + 1)
            g = re.sub(r'(?<![0-9A-Za-z])-', '\u2212', g)   # true minus
            out.append(('\\textsuperscript{%s}' if ch == '^' else '\\textsubscript{%s}') % g)
        elif ch.isalpha():
            j = i
            while j < len(m) and m[j].isalpha(): j += 1
            w = m[i:j]
            out.append('\\emph{%s}' % w if len(w) == 1 and w in VARS else w)
            i = j
        else:
            if ch in '=<>≤≥≈±' and out and not out[-1].endswith(' '):
                out.append(' ')
                out.append(ch); out.append(' ')
            else:
                out.append(ch)
            i += 1
    txt = ''.join(out)
    txt = re.sub(r' {2,}', ' ', txt).strip()
    txt = txt.replace('\\textasciitilde{} ', '\\textasciitilde{}')
    txt = re.sub(r'(?<=[0-9]) ?× ?(?=10)', ' × ', txt)
    txt = re.sub(r'^-(?=[A-Za-z\\])', '\u2212', txt)   # leading minus, not a hyphen
    return txt

def apply(tex):
    """Rewrite every inline $...$ in tex. Display math is untouched."""
    done = [0, 0]
    def one(mo):
        c = convert(mo.group(1).strip())
        if c is None:
            done[1] += 1
            return mo.group(0)
        done[0] += 1
        return c
    out = re.sub(r'(?<!\$)\$([^$]{1,120})\$(?!\$)', one, tex)
    return out, done
