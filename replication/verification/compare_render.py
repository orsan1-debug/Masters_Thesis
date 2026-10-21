"""compares the verify_render.R PDF w/ the stored one. display numbers + all numbers + page words
    python3 replication/verification/compare_render.py [stored.pdf] [render.pdf]"""

import csv, os, re, subprocess, sys, unicodedata
from collections import Counter

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(os.path.dirname(here))
stored = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, "output", "rendered", "thesis", "AAAAA-Working-Thesis.pdf")
render = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "thesis_render.pdf")

# NFKC + dashes to hyphens. page numbers dropped as the PDFs can differ in line and page breaks
dashes = str.maketrans({"ꞏ": "·", "–": "-", "—": "-", "−": "-", "‐": "-", "‑": "-"})
caption = re.compile(r"^(Appendix )?(Table|Figure) (\d+)[:.]")

def lines(pdf):
    """page and line for every line but the page numbers. blank lines kept as empty strings"""
    text = subprocess.run(["pdftotext", "-layout", "-enc", "UTF-8", pdf, "-"], capture_output=True,
                          text=True, check=True).stdout
    out = []
    for page, page_text in enumerate(text.split("\f"), start=1):
        for line in unicodedata.normalize("NFKC", page_text).translate(dashes).splitlines():
            line = " ".join(line.split())
            if not re.fullmatch(r"\d+", line):
                out.append((page, line))
    return out

def numbers(line):
    return [w.strip("([{").rstrip(".,;:)]}") for w in line.split() if re.search(r"\d", w)]

def is_row(line):
    """a line mostly of numbers. a sentence like the end of a caption is not one"""
    n = numbers(line)
    return len(n) >= 2 and len(n) >= 0.4 * len(line.split()) and not line.endswith(".")

def block(region):
    """first run of rows in region through gaps of up to 2 other lines
    a row repeated on a continuation page is kept once"""
    rows, gap = [], 0
    for _, line in region:
        if not line:
            continue
        if is_row(line):
            if line not in rows:
                rows.append(line)
            gap = 0
        elif rows:
            gap += 1
            if gap > 2:
                break
    return [x for line in rows for x in numbers(line)]

def displays(doc):
    """page and row numbers in order for each table and figure. a row is a line mostly of numbers
    table rows follow the caption and figure labels come before it"""
    caps = [(k, page, m.group(0)[:-1]) for k, (page, line) in enumerate(doc) for m in [caption.match(line)] if m]
    out = {}
    for j, (k, page, name) in enumerate(caps):
        if "Table" in name:
            end = caps[j + 1][0] if j + 1 < len(caps) else len(doc)
            out[name] = (page, block(doc[k + 1:min(end, k + 120)]))
        else:
            start = caps[j - 1][0] + 1 if j else 0
            out[name] = (page, block([x for x in doc[start:k] if x[0] == page][::-1]))
    return out

a, b = lines(stored), lines(render)
da, db = displays(a), displays(b)
rows = []
for name, (page, na) in da.items():
    pb, nb = db.get(name, ("", []))
    first = next((i for i, (x, y) in enumerate(zip(na, nb)) if x != y), min(len(na), len(nb)))
    same = na == nb
    rows.append({"display": name, "page": page, "render_page": pb, "numbers": len(na), "render_numbers": len(nb),
                 "identical": same, "first_difference": "" if same else
                 f"{' '.join(na[first:first + 6])} | {' '.join(nb[first:first + 6])}"})
with open(os.path.join(here, "render_compare.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

ca = Counter(x for _, line in a for x in numbers(line))
cb = Counter(x for _, line in b for x in numbers(line))
pages = lambda doc, x: " ".join(sorted({str(p) for p, line in doc if x in numbers(line)}, key=int))
diff = [{"number": x, "in": "stored PDF", "times": n, "pages": pages(a, x)} for x, n in (ca - cb).items()] + \
       [{"number": x, "in": "render", "times": n, "pages": pages(b, x)} for x, n in (cb - ca).items()]
with open(os.path.join(here, "render_numbers.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["number", "in", "times", "pages"])
    w.writeheader()
    w.writerows(diff)

def words(pdf):
    """words of each page in reading order. words split at a line end are joined and hyphens
    between letters dropped so a word reads the same wherever the line breaks"""
    text = subprocess.run(["pdftotext", "-enc", "UTF-8", pdf, "-"], capture_output=True,
                          text=True, check=True).stdout
    pages = []
    for page_text in text.split("\f")[:-1]:
        joined = []
        for line in unicodedata.normalize("NFKC", page_text).translate(dashes).splitlines():
            line = " ".join(line.split())
            if not line or re.fullmatch(r"\d+", line):
                continue
            if joined and re.search(r"[A-Za-z]-$", joined[-1]) and re.match(r"[a-z]", line):
                joined[-1] += line
            else:
                joined.append(line)
        pages.append([re.sub(r"(?<=[A-Za-z])-(?=[A-Za-z])", "", w) for w in " ".join(joined).split()])
    return pages

def same_letters(x, y):
    return sorted("".join(x)) == sorted("".join(y))

# reordered = same words in other order (eg figure labels). URL split = URL or DOI broken elsewhere
# page break elsewhere = same letters once the next or previous page is added
wa, wb = words(stored), words(render)
text_rows = []
for p in range(max(len(wa), len(wb))):
    x, y = (wa[p] if p < len(wa) else []), (wb[p] if p < len(wb) else [])
    near = [q for q in (p - 1, p + 1) if 0 <= q < min(len(wa), len(wb))]
    verdict = ("identical" if x == y else "reordered" if Counter(x) == Counter(y) else
               "URL split elsewhere" if same_letters(x, y) else
               "page break elsewhere" if any(same_letters(x + wa[q], y + wb[q]) for q in near) else "differs")
    text_rows.append({"page": p + 1, "words": len(x), "render_words": len(y), "verdict": verdict,
                      "only_stored": " ".join((Counter(x) - Counter(y)).elements()),
                      "only_render": " ".join((Counter(y) - Counter(x)).elements())})
with open(os.path.join(here, "render_text.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(text_rows[0]))
    w.writeheader()
    w.writerows(text_rows)

print(f"displays: {len(rows)}, numbers identical in {sum(r['identical'] for r in rows)}")
for r in rows:
    if not r["identical"]:
        print(f"  {r['display']} (p. {r['page']}): {r['first_difference']}")
print(f"numbers in the document: stored PDF {sum(ca.values())}, render {sum(cb.values())}. "
      f"{sum((ca - cb).values())} only in the stored PDF, {sum((cb - ca).values())} only in the render")
print(f"pages: stored PDF {len(wa)}, render {len(wb)}. " +
      ", ".join(f"{n} {v}" for v, n in Counter(r["verdict"] for r in text_rows).items()))
