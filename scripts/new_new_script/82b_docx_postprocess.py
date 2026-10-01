#!/usr/bin/env python3
"""
Give every table a full border grid and force the whole document to black.

Pandoc's docx writer emits tables in its own "Table" style, which draws a rule
under the header row and nothing else, and it colours hyperlinks and a few
other runs. Neither is reachable from the command line, so both are set here
by editing the package XML directly.
"""
import re, shutil, sys, zipfile

SRC = sys.argv[1]
DST = sys.argv[2] if len(sys.argv) > 2 else SRC
W = 'w:'
EDGES = ("top", "left", "bottom", "right", "insideH", "insideV")
BORDERS = ("<w:tblBorders>" +
           "".join('<w:%s w:val="single" w:sz="8" w:space="0" w:color="000000"/>' % e
                   for e in EDGES) +
           "</w:tblBorders>")

def grid(xml):
    """Insert a full border set into every table, replacing any it already has."""
    xml = re.sub(r'<w:tblBorders>.*?</w:tblBorders>', '', xml, flags=re.S)
    # tblStyle/tblW/... must stay in schema order; tblBorders follows tblW
    def ins(m):
        inner = m.group(1)
        k = inner.find("<w:tblW")
        if k >= 0:
            k = inner.index("/>", k) + 2 if inner[k:k+40].find("/>") >= 0 else len(inner)
            return "<w:tblPr>" + inner[:k] + BORDERS + inner[k:] + "</w:tblPr>"
        return "<w:tblPr>" + inner + BORDERS + "</w:tblPr>"
    return re.sub(r'<w:tblPr>(.*?)</w:tblPr>', ins, xml, flags=re.S)

TEXTW = 9360          # 6.5 in in twips

def columns(xml):
    """Drop all-empty columns and give the rest content-proportional widths.

    Pandoc derives column widths from the LaTeX column spec, which leaves a
    long first column one letter wide, and \\multicolumn headers can leave a
    trailing column that is empty in every row. Widths are set explicitly
    rather than left to autofit so Word and everything else agree.
    """
    from xml.etree import ElementTree as ET
    W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    ET.register_namespace("w", W[1:-1])

    NS = ('xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
          'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
          'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" '
          'xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math" '
          'xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture" '
          'xmlns:v="urn:schemas-microsoft-com:vml" '
          'xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" '
          'xmlns:wps="http://schemas.microsoft.com/office/word/2010/wordprocessingShape" '
          'xmlns:w14="http://schemas.microsoft.com/office/word/2010/wordml"')

    def fix(m):
        # a bare <w:tbl> fragment has no namespace declarations of its own
        try:
            tbl = ET.fromstring("<root %s>%s</root>" % (NS, m.group(0))).find(W + "tbl")
        except ET.ParseError:
            return m.group(0)
        if tbl is None: return m.group(0)
        rows = tbl.findall(W + "tr")
        if not rows: return m.group(0)
        ncol = max(len(r.findall(W + "tc")) for r in rows)
        text = [[""] * ncol for _ in rows]
        for i, r in enumerate(rows):
            for j, c in enumerate(r.findall(W + "tc")):
                if j < ncol:
                    text[i][j] = "".join(x.text or "" for x in c.iter(W + "t")).strip()
        keep = [j for j in range(ncol) if any(text[i][j] for i in range(len(rows)))]
        if not keep: return m.group(0)
        if len(keep) < ncol:
            for r in rows:
                cells = r.findall(W + "tc")
                if len(cells) == ncol:
                    for j, c in enumerate(cells):
                        if j not in keep: r.remove(c)
            text = [[row[j] for j in keep] for row in text]
        n = len(keep)

        # A column must at least fit its longest unbreakable token, or numbers
        # wrap mid-value ("0.6 / 25"). Whatever is left over is shared out by
        # total content length, damped so one prose column cannot take it all.
        CH, PAD = 118, 260                      # twips per character, cell padding
        need, want = [], []
        for j in range(n):
            # cap the demand: a long prose word may wrap, a number may not
            tok = min(10, max((len(w) for row in text for w in row[j].split()), default=3))
            need.append(CH * tok + PAD)
            want.append(max(3, max(len(row[j]) for row in text)) ** 0.6)
        if sum(need) > TEXTW:
            k = TEXTW / sum(need)
            wid = [max(500, int(w * k)) for w in need]
        else:
            slack = TEXTW - sum(need)
            tot = sum(want)
            wid = [need[j] + int(slack * want[j] / tot) for j in range(n)]
        wid[-1] += TEXTW - sum(wid)

        grid = tbl.find(W + "tblGrid")
        if grid is not None:
            for gc in list(grid): grid.remove(gc)
            for w in wid:
                ET.SubElement(grid, W + "gridCol").set(W + "w", str(w))
        for r in rows:
            for j, c in enumerate(r.findall(W + "tc")):
                pr = c.find(W + "tcPr")
                if pr is None: continue
                for e in pr.findall(W + "tcW"): pr.remove(e)
                e = ET.Element(W + "tcW"); e.set(W + "w", str(wid[min(j, n-1)]))
                e.set(W + "type", "dxa"); pr.insert(0, e)
        pr = tbl.find(W + "tblPr")
        if pr is not None:
            for e in pr.findall(W + "tblW") + pr.findall(W + "tblLayout"): pr.remove(e)
            e = ET.Element(W + "tblW"); e.set(W + "w", str(TEXTW)); e.set(W + "type", "dxa")
            pr.insert(0, e)
            e = ET.Element(W + "tblLayout"); e.set(W + "type", "fixed"); pr.insert(1, e)
        out = ET.tostring(tbl, encoding="unicode")
        return re.sub(r'\s+xmlns:[a-z0-9]+="[^"]*"', '', out)

    return re.sub(r'<w:tbl>.*?</w:tbl>', fix, xml, flags=re.S)

def blacken(xml):
    # Theme attributes go FIRST. A tag like
    #   <w:color w:val="4F81BD" w:themeColor="accent1"/>
    # does not match a w:val-then-close pattern, so stripping the theme
    # attributes afterwards would leave the colour behind for a second pass.
    for a in ("themeColor", "themeTint", "themeShade", "themeFill",
              "themeFillTint", "themeFillShade"):
        xml = re.sub(r'\s+w:%s="[^"]*"' % a, '', xml)
    xml = re.sub(r'<w:color\b[^>]*?w:val="(?!000000)[^"]*"[^>]*?/>',
                 '<w:color w:val="000000"/>', xml)
    xml = re.sub(r'<w:color\b[^>]*?w:val="(?!000000)[^"]*"[^>]*?>',
                 '<w:color w:val="000000"/>', xml)
    return xml

zin = zipfile.ZipFile(SRC)
names = zin.namelist()
tmp = DST + ".tmp"
with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
    for n in names:
        data = zin.read(n)
        if n.endswith(".xml") and ("document" in n or "styles" in n or "numbering" in n):
            x = data.decode("utf8")
            if "document" in n: x = columns(grid(x))
            x = blacken(x)
            data = x.encode("utf8")
        zout.writestr(n, data)
zin.close()
shutil.move(tmp, DST)

z = zipfile.ZipFile(DST)
d = z.read("word/document.xml").decode("utf8")
st = z.read("word/styles.xml").decode("utf8")
print("  tables            %d" % d.count("<w:tbl>"))
print("  with full borders %d" % d.count("<w:tblBorders>"))
bad = set(re.findall(r'<w:color w:val="([^"]+)"', d + st)) - {"000000"}
print("  non-black colours %s" % (sorted(bad) or "none"))
print("  images            %d" % d.count("<a:blip"))
print("  fixed-width tables %d" % d.count('<w:tblLayout w:type="fixed"'))
