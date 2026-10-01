import zipfile, re
from xml.etree import ElementTree as ET
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def read_xlsx(path, sheet=0):
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")):
            ss.append("".join(t.text or "" for t in si.iter(NS+"t")))
    names = sorted(n for n in z.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", n))
    root = ET.fromstring(z.read(names[sheet]))
    rows = []
    for r in root.iter(NS+"row"):
        cells = {}
        for c in r.iter(NS+"c"):
            ref = c.get("r"); col = re.match(r"([A-Z]+)", ref).group(1)
            v = c.find(NS+"v"); t = c.get("t")
            if t == "inlineStr":
                isn = c.find(NS+"is"); val = "".join(x.text or "" for x in isn.iter(NS+"t")) if isn is not None else ""
            elif v is None: val = None
            elif t == "s": val = ss[int(v.text)]
            elif t in ("str", "e"): val = v.text
            elif t == "b": val = (v.text == "1")
            else:
                try:
                    f = float(v.text)
                    val = int(f) if f.is_integer() and abs(f) < 2**53 else f
                except (TypeError, ValueError):
                    val = v.text
            cells[col] = val
        rows.append(cells)
    if not rows: return [], []
    hdr_cells = rows[0]
    cols = sorted(hdr_cells, key=lambda c: (len(c), c))
    hdr = [hdr_cells[c] for c in cols]
    out = []
    for r in rows[1:]:
        if not any(v is not None for v in r.values()): continue
        out.append({h: r.get(c) for h, c in zip(hdr, cols)})
    return hdr, out
