#!/usr/bin/env python3
"""Sanitize sequence lines of final proteomes / CDS that DIAMOND (and other tools) reject.
2026-10-19: OrthoFinder Annelida add_v2 died in `diamond makedb` on Alvpom ("Invalid character in
sequence: '.'"; the first attempt deadlocked instead of erroring). A scan found one internal '.'
(Alvpom_004721; the matching CDS codon is TGA, i.e. a stop) and spaces inside sequence lines of the
lab proteomes Hha / Rve.
Fix: drop whitespace inside sequence lines; protein '.' -> '*' (same convention as 06b_prep_extra.py,
"some releases use '.' as stop"). Headers, IDs and order are unchanged. The original is copied to
backups/<file>.<ts>; a symlinked file is replaced by a sanitized regular file (the target, which may
sit in the read-only ~/deepsea tree, is never touched). The new file keeps the old mtime: the change is
content-equivalent for compleasm/hmmsearch, so 14_protqc_loop does not re-queue QC.
Usage: python 06c_sanitize_fasta.py [--dry]"""
import glob, os, re, shutil, sys, time

R = os.path.expanduser("~/deepsea_conv")
DRY = "--dry" in sys.argv
TS = time.strftime("%Y%m%d_%H%M%S")
BAD_P = re.compile(r"[^A-Za-z*]")
BAD_N = re.compile(r"[^A-Za-z]")

def fix(path, prot):
    bad = BAD_P if prot else BAD_N
    lines = open(path).read().split("\n")
    out, n = [], 0
    for l in lines:
        if l.startswith(">") or not l:
            out.append(l); continue
        if bad.search(l):
            n += 1
            l = "".join(l.split())
            if prot:
                l = l.replace(".", "*")
        out.append(l)
    left = sum(bool(bad.search(l)) for l in out if l and not l.startswith(">"))
    link = os.readlink(path) if os.path.islink(path) else ""
    if n and not DRY:
        st = os.stat(path)                       # follows the link: mtime the QC loop compares against
        b = f"{R}/backups/{os.path.basename(path)}.{TS}"
        if link:
            open(b + ".symlink", "w").write(link + "\n")
        else:
            shutil.copy2(path, b)
        tmp = path + ".tmp_sanitize"
        with open(tmp, "w") as o:
            o.write("\n".join(out))
        os.utime(tmp, (st.st_atime, st.st_mtime))
        os.replace(tmp, path)                    # replaces the link itself, not its target
    return n, left, link

for pat, prot in ((f"{R}/data/proteomes/*.faa", True), (f"{R}/data/outgroups/*/*.faa", True),
                  (f"{R}/data/cds/*.fna", False), (f"{R}/data/outgroups/*/*.fna", False)):
    for f in sorted(glob.glob(pat)):
        n, left, link = fix(f, prot)
        if n:
            print(f"{'DRY ' if DRY else ''}{os.path.relpath(f, R)}\tfixed_lines={n}\tstill_bad={left}"
                  + (f"\tsymlink->{link}" if link else ""))
print("scan done")
