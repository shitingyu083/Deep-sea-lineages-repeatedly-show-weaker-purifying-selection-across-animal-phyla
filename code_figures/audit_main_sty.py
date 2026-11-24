"""Audit condensed main.md: word counts per section, legend lengths, reference
first-citation order, figure / Extended Data citation order. Writes a report file."""
import sys, re
src, out = sys.argv[1], sys.argv[2]
T = open(src, encoding="utf-8").read()
lines = T.split("\n")

def words(s):
    s = re.sub(r"\^[^^]*\^", "", s)          # superscript citations
    s = re.sub(r"[*_`#]", "", s)
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9.,'’/−–\-]*", s))

sec, cur = {}, None
for ln in lines:
    if ln.startswith("## "):
        cur = ln[3:].strip()
        sec.setdefault(cur, [])
    elif cur:
        sec[cur].append(ln)
R = []
title = lines[0].lstrip("# ").strip()
R.append("title chars %d: %s" % (len(title), title))
for k in ("Summary", "Abstract", "Main", "Results", "Discussion", "Methods"):
    body = [l for l in sec.get(k, []) if not l.startswith("#")]
    R.append("%s words %d" % (k, words(" ".join(body))))
mt = sum(words(" ".join(l for l in sec.get(k, []) if not l.startswith("#"))) for k in ("Main", "Results", "Discussion"))
R.append("MAIN TEXT (Main+Results+Discussion) words %d" % mt)
for ln in lines:
    if ln.startswith("### ") and lines.index(ln) < lines.index("## Methods"):
        R.append("subheading chars %d: %s" % (len(ln[4:]), ln[4:]))
for ln in lines:
    m = re.match(r"\*\*((Extended Data )?Fig\. \d+) \|", ln)
    if m:
        R.append("legend %s words %d" % (m.group(1), words(ln)))

# reference order in text before References
body = T.split("## References")[0]
seen = []
for m in re.finditer(r"\^([0-9,–\-]+)\^", body):
    for part in m.group(1).split(","):
        if re.match(r"^\d+[–\-]\d+$", part):
            a, b = re.split(r"[–\-]", part)
            nums = range(int(a), int(b) + 1)
        elif part.isdigit():
            nums = [int(part)]
        else:
            continue
        for n in nums:
            if n not in seen:
                seen.append(n)
for m in re.finditer(r"ref\. (\d+)", body):
    pass
R.append("ref first-appearance order: %s" % seen)
R.append("ref order OK: %s (n=%d)" % (seen == sorted(seen) and seen == list(range(1, len(seen) + 1)), len(seen)))
nref = len(re.findall(r"^\d+\. ", T.split("## References")[1].split("## Acknowledgements")[0], re.M))
R.append("reference list entries %d" % nref)

def order(pat):
    o = []
    body2 = re.sub(r"Extended Data Figs?\. \d", "", body) if pat.startswith("Figs") else body
    for m in re.finditer(pat, body2):
        k = int(m.group(1))
        if k not in o:
            o.append(k)
    return o
R.append("Fig order: %s" % order(r"Figs?\. (\d)"))
R.append("ED order: %s" % order(r"Extended Data Figs?\. (\d)"))
R.append("SuppTable order: %s" % order(r"Supplementary Tables? (\d+)"))
R.append("VERIFY count %d; AUTHOR TO CONFIRM count %d" % (T.count("VERIFY"), T.count("AUTHOR TO CONFIRM")))
open(out, "w", encoding="utf-8").write("\n".join(R) + "\n")
