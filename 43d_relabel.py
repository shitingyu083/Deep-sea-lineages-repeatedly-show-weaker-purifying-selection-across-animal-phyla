#!/usr/bin/env python3
"""Draft -> final species tree for a clade whose selection run (24_selection_run.sh) was labelled on the draft tree
(2026-10-20). Per HOG both trees are pruned to the alignment's species as 23_label_trees.py does and compared by a
canonical rooted signature (tip names and nesting, unary nodes included; child order and branch lengths ignored).
The labels (FG/BG branches, deep and shallow lineages, null sets) depend on that rooted topology only, and
HyPhy and the IQ-TREE run inside CSUBST re-estimate branch lengths, so a HOG with an identical induced topology keeps
its label files and results; a HOG whose induced topology changed is relabelled on FINAL and re-run.
Default (read-only): OUT/diff.tsv (per HOG), OUT/changed.txt, summary on stdout (OUT = work/codon/<clade>/relabel43d,
  or --out DIR). Exit 0.
--verify: also replays 23 on DRAFT in full (its order, its one RNG along the sorted HOGs, --null-from-shallow 5) and
  compares the result with trees/ and lineages.tsv byte for byte (OUT/verify_bad.tsv; exit 2 on a mismatch): the
  files on disk are what 23 wrote from DRAFT with today's habitat labels.
--apply (implies --verify; safe while 24 runs on the clade; refused on a verify mismatch, if DRAFT is not the md5 in
  trees/SOURCE_TREE.txt, if FINAL lacks a species of an alignment, or if OUT/APPLYING or OUT/APPLIED exists).
  For the changed HOGs:
  0. their new label files are written from FINAL into OUT/tmp/ with 23's functions (null sets drawn with a per-HOG
     RNG seeded 20261020 + crc32(hog): 23's single RNG cannot be continued for a subset);
  1. all their trees/ entries (rooted.nwk and the hyphy.nwk / fg.txt of every set) move to OUT/stale/trees/: a job
     dispatched from now on finds no labels and returns before writing anything, a running one fails fast;
  2. their in-flight jobs (24 wrappers: own uid, environment W = this clade) are waited for; after --wait-min minutes
     the wrappers' descendants (HyPhy / CSUBST / IQ-TREE) are killed, never a wrapper (xargs would abort), then
     stragglers (own non-shell processes working in their cs/ dirs); if anything is left 10 min later, step 1 is
     undone and nothing else changes (OUT/APPLY_FAILED);
  3. their outputs (sel/<hog>*.json, sel_logs/<hog>*.log, cs/<hog>{,.nullK}) move to OUT/stale/ (completed ones are
     draft-tree results, kept for a species-tree sensitivity comparison);
  4. the new files move into trees/ (renames; per HOG rooted.nwk, then the null sets, the real set last, fg.txt
     before hyphy.nwk), so a job dispatched meanwhile sees no labels or a complete new set;
  5. lineages.tsv rows replaced (backup OUT/lineages.tsv.pre43d); OUT/SOURCE.txt and trees/SOURCE_TREE.txt (FINAL:
     43b never re-runs 23 on the clade) written; APPLYING renamed APPLIED. Every move is in OUT/manifest.tsv.
  The skipped and moved jobs are re-run by a second 24 pass (43e_relabel_chain.sh: TAG=re43d, same tier lists; the
  HyPhy analysis set is hylist & FINAL-eligible, see diff.tsv).
--wdir DIR: run on a copy of work/codon/<clade> (tests).
Usage: python 43d_relabel.py CLADE DRAFT FINAL [--verify] [--apply] [--wait-min 20] [--out DIR] [--wdir DIR]"""
import fcntl
import hashlib
import itertools
import math
import os
import random
import re
import shutil
import signal
import sys
import time
import zlib

from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv")
sys.path.insert(0, f"{R}/scripts")
from importlib import import_module  # noqa: E402

L23 = import_module("23_label_trees")
LAB = L23.LAB
KNULL = 5  # 43b runs 23 with --null-from-shallow 5
UID = os.getuid()
WRAP = 'one "$0" "$1"'  # 24's xargs command: bash -c 'one "$0" "$1"' HOG SET


def log(msg):
    print(f"{time.strftime('%F_%T')} {msg}", flush=True)


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def opts():
    a = sys.argv[1:]
    if len(a) < 3 or a[0].startswith("-"):
        sys.exit(__doc__)
    val = lambda k, d: a[a.index(k) + 1] if k in a else d  # noqa: E731
    return a[0], a[1], a[2], {"verify": "--verify" in a, "apply": "--apply" in a,
                              "wait": float(val("--wait-min", 20)), "out": val("--out", None),
                              "wdir": val("--wdir", None)}


def load(p):
    t = Tree(p, format=1)
    leaves = set(t.get_leaf_names())
    inner = {n.name for n in t.traverse() if not n.is_leaf() and n.name}
    assert len(leaves) == len(t), f"duplicate tip names in {p}"
    # 23 selects the species with t.search_nodes(name=s): same as `s in leaves` unless an internal node has the name
    assert not inner & leaves, f"internal node named like a tip in {p}: {sorted(inner & leaves)[:5]}"
    return t, leaves


def species(W, hog):
    return [l[1:].strip() for l in open(f"{W}/aln/{hog}.fas") if l.startswith(">")]  # as 23


def pruned(full, leaves, sp):
    t = full.copy()
    t.prune([s for s in sp if s in leaves], preserve_branch_length=True)  # as 23
    return t


def label_sets(t, rng, K=KNULL):
    """23's labels: states, deep lineages (FG), shallow lineages, null picks (distinct sets of len(fg) shallow
    lineages: all combinations if there are <= K, else K random draws)."""
    st = L23.states(t, LAB)
    fg = L23.maximal(t, st, "D")
    sh = L23.maximal(t, st, "S")
    picks = []
    if len(sh) >= len(fg) and K:
        if math.comb(len(sh), len(fg)) <= K:
            picks = [list(c) for c in itertools.combinations(sh, len(fg))]
        else:
            seen = set()
            while len(picks) < K:
                p = rng.sample(range(len(sh)), len(fg))
                key = frozenset(p)
                if key not in seen:
                    seen.add(key)
                    picks.append([sh[i] for i in p])
    return st, fg, sh, picks


def write_hog(od, hog, t, st, fg, sh, picks):
    """23's files for one HOG into od; returns its lineages.tsv row and the file names."""
    t.write(outfile=f"{od}/{hog}.rooted.nwk", format=5)
    L23.write_set(f"{od}/{hog}", t, st, fg, "S")
    for k, pick in enumerate(picks, 1):
        L23.write_set(f"{od}/{hog}.null{k}", t, st, pick, "S")
    names = [f"{hog}.rooted.nwk", f"{hog}.hyphy.nwk", f"{hog}.fg.txt"]
    names += [f"{hog}.null{k}.{x}" for k in range(1, len(picks) + 1) for x in ("hyphy.nwk", "fg.txt")]
    return (f"{hog}\t{len(t)}\t{len(fg)}\t{len(sh)}\t" +
            ";".join(",".join(n.get_leaf_names()) for n in fg) + "\n"), names


def hog_rng(hog):
    return random.Random(20261020 + zlib.crc32(hog.encode()))


def sig(t):
    """Canonical rooted topology: tip names and nesting, unary nodes included, child order and lengths ignored."""
    s = {}
    for n in t.traverse("postorder"):
        s[n] = n.name if n.is_leaf() else "(" + ",".join(sorted(s[c] for c in n.children)) + ")"
    return s[t]


# entries per HOG (24's and 23's names; a HOG name never contains .nullK or these suffixes)
PAT = {"trees": r"(\.null\d+)?\.(rooted\.nwk|hyphy\.nwk|fg\.txt)", "sel": r"(\.null\d+)?\.(relax|bph)\.json",
       "sel_logs": r"(\.null\d+)?\.(relax|bph)\.log", "cs": r"(\.null\d+)?"}


def index(W, hogs, dirs=tuple(PAT)):
    """{dir: {hog: sorted entry names}} of the HOGs in hogs."""
    out = {}
    for d, suf in ((d, PAT[d]) for d in dirs):
        rx, ix = re.compile(r"^(.+?)" + suf + "$"), {}
        for fn in sorted(os.listdir(f"{W}/{d}")) if os.path.isdir(f"{W}/{d}") else []:
            m = rx.match(fn)
            if m and m.group(1) in hogs:
                ix.setdefault(m.group(1), []).append(fn)
        out[d] = ix
    return out


def scan(W, hogs, draft, final, ix, lin, verify, tmp, vtmp):
    """One pass in 23's order. verify: 23 replayed on DRAFT (one RNG along the sorted HOGs) into vtmp and compared
    with trees/ + lineages.tsv byte for byte. Always: FINAL labels per HOG (nulls from hog_rng); the files of a HOG
    whose induced topology changed are written to tmp. Returns per-HOG rows, new lineages rows, mismatches."""
    (td, ld), (tf, lf) = draft, final
    rng = random.Random(20261017)
    rows, new, bad = [], {}, []
    for i, hog in enumerate(hogs):
        sp = species(W, hog)
        a, b = pruned(td, ld, sp), pruned(tf, lf, sp)
        if verify:
            st, fg, sh, picks = label_sets(a, rng)
            row, _ = write_hog(vtmp, hog, a, st, fg, sh, picks)
            got = sorted(os.listdir(vtmp))
            if got != ix["trees"].get(hog, []):
                bad.append((hog, f"{len(got)} files vs {len(ix['trees'].get(hog, []))} in trees/"))
            else:
                for fn in got:
                    if open(f"{vtmp}/{fn}", "rb").read() != open(f"{W}/trees/{fn}", "rb").read():
                        bad.append((hog, fn))
                        break
            if lin.get(hog) != row:
                bad.append((hog, "lineages.tsv row"))
            for fn in got:
                os.remove(f"{vtmp}/{fn}")
        ch = len(b) < len(sp) or len(a) < len(sp) or sig(a) != sig(b)
        st, fg, sh, picks = label_sets(b, hog_rng(hog))
        old = lin[hog].rstrip("\n").split("\t")
        fo = {frozenset(x.split(",")) for x in old[4].split(";") if x}
        fn_ = {frozenset(n.get_leaf_names()) for n in fg}
        if not ch:
            assert fo == fn_ and int(old[3]) == len(sh), f"{hog}: same topology, different labels"
        n1 = f"{W}/trees/{hog}.null1.fg.txt"
        el_d = int(old[2]) >= 2 and os.path.exists(n1) and os.path.getsize(n1) > 0  # 43c's rule
        el_f = len(fg) >= 2 and len(picks) > 0
        if ch:
            new[hog] = write_hog(tmp, hog, b, st, fg, sh, picks) if tmp else None
        rows.append([hog, len(sp), len(sp) - len(b), int(ch), int(old[2]), len(fg), int(fo == fn_), int(old[3]),
                     len(sh), int(el_d), int(el_f), len(ix["sel"].get(hog, [])), len(ix["cs"].get(hog, [])),
                     len(sp) - len(a)])
        if (i + 1) % 2000 == 0:
            log(f"scan {i + 1}/{len(hogs)}: changed {sum(r[3] for r in rows)}, mismatches {len(bad)}")
    return rows, new, bad


def wrappers(W, hogs):
    """24's job wrappers (own uid, argv: bash -c 'one "$0" "$1"' HOG SET, environment W) of the HOGs in hogs."""
    Wr, out = os.path.realpath(W), {}
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            if os.stat(f"/proc/{p}").st_uid != UID:
                continue
            a = open(f"/proc/{p}/cmdline", "rb").read().split(b"\0")
            if (len(a) < 5 or not a[0].endswith(b"bash") or a[1] != b"-c" or a[2] != WRAP.encode()
                    or a[3].decode() not in hogs):
                continue
            env = dict(x.split(b"=", 1) for x in open(f"/proc/{p}/environ", "rb").read().split(b"\0") if b"=" in x)
            if os.path.realpath(env.get(b"W", b"").decode() or "/nonexistent") == Wr:
                out[int(p)] = (a[3].decode(), a[4].decode())
        except OSError:
            continue
    return out


def stragglers(W, hogs):
    """Own processes (shells excluded: 24's wrappers cd there) whose working directory is a CSUBST dir of a HOG in
    hogs: an IQ-TREE whose CSUBST was killed is no wrapper descendant any more and keeps running there."""
    base, rx, out = os.path.realpath(f"{W}/cs") + "/", re.compile(r"^(.+?)" + PAT["cs"] + "$"), {}
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            if os.stat(f"/proc/{p}").st_uid != UID:
                continue
            cwd = os.readlink(f"/proc/{p}/cwd")
            if not cwd.startswith(base):
                continue
            m = rx.match(cwd[len(base):].split("/")[0])
            a0 = os.path.basename(open(f"/proc/{p}/cmdline", "rb").read().split(b"\0")[0])
            if m and m.group(1) in hogs and not a0.endswith(b"sh"):
                out[int(p)] = ("straggler", m.group(0))
        except OSError:
            continue
    return out


def kill(q, why, termed):
    """SIGTERM own process q; SIGKILL from 30 s after the first SIGTERM on. Logs the first of each."""
    try:
        if os.stat(f"/proc/{q}").st_uid != UID:
            return
        t = termed.setdefault(q, [time.time(), 0])
        sg = signal.SIGKILL if time.time() - t[0] > 30 else signal.SIGTERM
        cmd = open(f"/proc/{q}/cmdline", "rb").read().replace(b"\0", b" ").decode(errors="replace")
        os.kill(q, sg)
        if t[1] < (2 if sg == signal.SIGKILL else 1):
            t[1] = 2 if sg == signal.SIGKILL else 1
            log(f"   {sg.name} {q} ({why}): {cmd[:150]}")
    except OSError:
        pass


def kill_desc(wr, termed):
    """kill() the descendants (HyPhy / CSUBST / IQ-TREE ...) of the wrappers, never a wrapper itself: GNU xargs
    aborts when a child of its own dies of a signal (a wrapper whose job was killed exits with 143/137: xargs goes on)."""
    ch = {}
    for p in os.listdir("/proc"):
        if p.isdigit():
            try:
                s = open(f"/proc/{p}/stat").read()
                ch.setdefault(int(s[s.rindex(")") + 2:].split()[1]), []).append(int(p))
            except (OSError, ValueError, IndexError):
                pass
    for w, job in wr.items():
        todo = list(ch.get(w, []))  # the whole subtree as it is now, so children of a killed process are hit too
        while todo:
            q = todo.pop()
            todo.extend(ch.get(q, []))
            kill(q, f"job {job[0]} {job[1]}", termed)


def move(src, dst, man, done):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.lexists(dst):
        raise FileExistsError(dst)
    os.replace(src, dst)  # a rename on the same file system, files and directories alike
    man.write(f"{src}\t{dst}\n")
    man.flush()
    done.append((src, dst))


def order(h, fn):
    """rooted.nwk first, then the null sets, the real set last; fg.txt before hyphy.nwk within a set."""
    r = fn[len(h):]
    return (0 if r == ".rooted.nwk" else 1 if r.startswith(".null") else 2, r.endswith(".hyphy.nwk"), r)


def drain(W, cs, wait):
    """Step 2: wait until no wrapper of a HOG in cs runs (descendants killed after wait min), then until no
    straggler remains (killed at once). False if anything is left after wait + 10 min."""
    t0, termed, last = time.time(), {}, None
    while True:
        wr = wrappers(W, cs)
        sg = {} if wr else stragglers(W, cs)
        el = (time.time() - t0) / 60
        if (len(wr), len(sg)) != last:
            log(f"step 2 at {el:.1f} min: {len(wr)} in-flight jobs {sorted(wr.values())[:6]}, {len(sg)} stragglers")
            last = (len(wr), len(sg))
        if not wr and not sg:
            return True
        if el > wait + 10:
            return False
        if wr and el > wait:
            kill_desc(wr, termed)
        for q, why in sg.items():
            kill(q, " ".join(why), termed)
        time.sleep(10)


def apply(W, out, changed, new, wait, src_note):
    """Steps 1-5 (docstring). True when applied; False when step 2 failed and step 1 was undone."""
    tr, st, cs = f"{W}/trees", f"{out}/stale", set(changed)
    with open(f"{out}/APPLYING", "w") as f:
        f.write(f"pid {os.getpid()} started {time.strftime('%F_%T')}\n")
    man = open(f"{out}/manifest.tsv", "a", newline="\n")
    man.write(f"# 43d apply {time.strftime('%F_%T')}: src\tdst\n")
    s1, ix1 = [], index(W, cs, ("trees",))["trees"]
    for h in changed:  # 1. label files away: jobs dispatched from now on return at once
        for fn in ix1.get(h, []):
            move(f"{tr}/{fn}", f"{st}/trees/{fn}", man, s1)
    log(f"step 1: {len(s1)} trees/ entries of {len(changed)} changed HOGs moved to {st}/trees")
    if not drain(W, cs, wait):  # 2. in-flight jobs
        for src, dst in reversed(s1):
            os.replace(dst, src)
            man.write(f"{dst}\t{src}\tundo\n")
        man.close()
        msg = (f"{time.strftime('%F_%T')} step 2 failed: jobs of changed HOGs still running; step 1 undone, nothing "
               f"else changed. Jobs dispatched for changed HOGs meanwhile were skipped: re-run 24 to fill them in.\n")
        open(f"{out}/APPLY_FAILED", "a").write(msg)
        os.remove(f"{out}/APPLYING")
        log(msg.split(" ", 1)[1].strip())
        return False
    ix = index(W, cs)
    assert not any(ix["trees"].get(h) for h in changed), "trees/ entries of changed HOGs reappeared"
    s3 = []
    for d in ("sel", "sel_logs", "cs"):  # 3. outputs away (complete ones = draft-tree results)
        for h in changed:
            for fn in ix[d].get(h, []):
                move(f"{W}/{d}/{fn}", f"{st}/{d}/{fn}", man, s3)
    log(f"step 3: {len(s3)} sel/ sel_logs/ cs/ entries moved to {st}/")
    s4 = []
    for h in changed:  # 4. new labels in: rooted, null sets, real set; fg.txt before hyphy.nwk
        for fn in sorted(new[h][1], key=lambda x: order(h, x)):
            move(f"{out}/tmp/{fn}", f"{tr}/{fn}", man, s4)
    log(f"step 4: {len(s4)} label files from FINAL in trees/")
    shutil.copy2(f"{W}/lineages.tsv", f"{out}/lineages.tsv.pre43d")  # 5. bookkeeping
    lines = open(f"{W}/lineages.tsv").readlines()
    body = [new[h][0] if h in new else ln for ln, h in ((ln, ln.split("\t", 1)[0]) for ln in lines[1:])]
    with open(f"{W}/lineages.tsv.tmp43d", "w", newline="\n") as f:
        f.write(lines[0] + "".join(body))
    os.replace(f"{W}/lineages.tsv.tmp43d", f"{W}/lineages.tsv")
    with open(f"{out}/SOURCE.txt", "w", newline="\n") as f:
        f.write(src_note["text"])
    with open(f"{tr}/SOURCE_TREE.txt.tmp43d", "w", newline="\n") as f:
        f.write(f"{src_note['final']}\t{src_note['fm']}\n")
    os.replace(f"{tr}/SOURCE_TREE.txt.tmp43d", f"{tr}/SOURCE_TREE.txt")
    man.close()
    os.replace(f"{out}/APPLYING", f"{out}/APPLIED")
    msg = (f"applied {time.strftime('%F_%T')}: {len(changed)} HOGs relabelled on FINAL; entries moved in steps "
           f"1/3/4: {len(s1)}/{len(s3)}/{len(s4)} (manifest.tsv)\n")
    open(f"{out}/APPLIED", "a").write(msg)
    log("step 5: lineages.tsv, SOURCE.txt, trees/SOURCE_TREE.txt written; " + msg.strip())
    return True


def summary(rows, hy, nh):
    """Log lines on what the switch to FINAL changes."""
    pc = lambda a, b: f"{a}/{b} ({100 * a / b:.1f}%)" if b else f"{a}/0"  # noqa: E731
    C = [r for r in rows if r[3]]
    L = [f"changed induced topology: {pc(len(C), len(rows))} HOGs; in hylist {pc(sum(r[0] in hy for r in C), len(hy))}"
         f"; in nullhy {pc(sum(r[0] in nh for r in C), len(nh))}",
         f"  of the changed: same deep-lineage tip sets {sum(r[6] for r in C)}, different n_fg {sum(r[4] != r[5] for r in C)}"
         f"; their outputs so far: {sum(r[11] for r in C)} sel json, {sum(r[12] for r in C)} cs dirs",
         f"eligible on DRAFT {sum(r[9] for r in rows)}, on FINAL {sum(r[10] for r in rows)} (lost "
         f"{sum(r[9] > r[10] for r in rows)}, gained {sum(r[10] > r[9] for r in rows)}); hylist & FINAL-eligible "
         f"{sum(r[0] in hy and r[10] for r in rows)}/{len(hy)}; nullhy & FINAL-eligible "
         f"{sum(r[0] in nh and r[10] for r in rows)}/{len(nh)}",
         f"HOGs with alignment species missing from FINAL: {sum(r[2] > 0 for r in rows)} (from DRAFT: "
         f"{sum(r[13] > 0 for r in rows)})"]
    for x in L:
        log(x)
    return L


COLS = ("hog n_sp n_missing_final changed n_fg_draft n_fg_final fg_tips_same n_sh_draft n_sh_final elig_draft "
        "elig_final n_sel_json n_cs_dirs n_missing_draft in_hylist in_nullhy").split()


def main():
    clade, dp, fp, opt = opts()
    W = os.path.realpath(opt["wdir"] or f"{R}/work/codon/{clade}")
    out = os.path.realpath(opt["out"] or f"{W}/relabel43d")
    os.makedirs(out, exist_ok=True)
    lock = open(f"{out}/.lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit(f"another 43d runs on {out}")
    dm, fm = md5(dp), md5(fp)
    src = open(f"{W}/trees/SOURCE_TREE.txt").read().rstrip("\n").split("\t")
    log(f"43d {clade} {'apply' if opt['apply'] else 'verify' if opt['verify'] else 'diff'}: W {W}; OUT {out}")
    log(f"DRAFT {dp} {dm}; FINAL {fp} {fm}; trees/SOURCE_TREE.txt {src}")
    if opt["apply"]:
        why = [f"{out}/{x} exists" for x in ("APPLIED", "APPLYING") if os.path.exists(f"{out}/{x}")]
        if len(src) < 2 or src[1] != dm:
            why.append("DRAFT md5 is not the one in trees/SOURCE_TREE.txt (the tree 23 labelled)")
        if os.stat(out).st_dev != os.stat(W).st_dev:
            why.append("OUT and W on different file systems (the moves must be renames)")
        if why:
            log("apply refused: " + "; ".join(why))
            return 1
    hogs = [fn[:-4] for fn in sorted(os.listdir(f"{W}/aln")) if fn.endswith(".fas")]
    lin = {ln.split("\t", 1)[0]: ln for ln in open(f"{W}/lineages.tsv").readlines()[1:]}
    assert set(lin) == set(hogs), f"lineages.tsv has {len(lin)} HOGs, aln/ {len(hogs)}"
    rd = lambda p: set(open(p).read().split()) if os.path.exists(p) else set()  # noqa: E731
    hy, nh = rd(f"{W}/tiers/hylist.txt"), rd(f"{W}/tiers/nullhy.txt")
    tmp = vtmp = None
    if opt["apply"]:
        tmp = f"{out}/tmp"
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
    if opt["verify"] or opt["apply"]:
        vtmp = f"{out}/verify_tmp"
        shutil.rmtree(vtmp, ignore_errors=True)
        os.makedirs(vtmp)
    t0 = time.time()
    rows, new, bad = scan(W, hogs, load(dp), load(fp), index(W, set(hogs)), lin, vtmp is not None, tmp, vtmp)
    log(f"scan of {len(hogs)} HOGs: {(time.time() - t0) / 60:.1f} min")
    with open(f"{out}/diff.tsv", "w", newline="\n") as f:
        f.write("\t".join(COLS) + "\n")
        f.writelines("\t".join(map(str, r + [int(r[0] in hy), int(r[0] in nh)])) + "\n" for r in rows)
    changed = [r[0] for r in rows if r[3]]
    with open(f"{out}/changed.txt", "w", newline="\n") as f:
        f.write("".join(h + "\n" for h in changed))
    summ = summary(rows, hy, nh)
    if vtmp:
        with open(f"{out}/verify_bad.tsv", "w", newline="\n") as f:
            f.writelines(f"{h}\t{w}\n" for h, w in bad)
        log(f"verify (23 replayed on DRAFT vs trees/ and lineages.tsv): {len(bad)} mismatches"
            + (f", first {bad[:3]}" if bad else ""))
        os.rmdir(vtmp)
    if not opt["apply"]:
        return 2 if bad else 0
    miss = sum(r[2] > 0 for r in rows)
    if bad or miss:
        log(f"apply refused: {len(bad)} verify mismatches; {miss} HOGs with alignment species missing from FINAL "
            f"(those need 43f_sync_aln.py first)")
        return 1
    text = (f"clade\t{clade}\ndraft\t{dp}\t{dm}\nfinal\t{fp}\t{fm}\n"
            f"changed\t{len(changed)} of {len(hogs)} HOGs (changed.txt; per-HOG comparison diff.tsv)\n"
            f"date\t{time.strftime('%F_%T')}\n" + "".join(f"summary\t{x}\n" for x in summ) +
            "note\tA HOG whose alignment-induced rooted topology is the same on DRAFT and FINAL keeps the label files "
            "23_label_trees.py wrote from DRAFT (same FG/BG labels and null sets; DRAFT branch lengths are starting "
            "values only: HyPhy and the IQ-TREE run inside CSUBST re-estimate them). A changed HOG was relabelled on "
            "FINAL with 23's functions (null sets from a per-HOG RNG seeded 20261020 + crc32(hog)); its draft-tree "
            "outputs are in stale/ (stale CSUBST dirs link the rooted tree by its trees/ path, which now holds the "
            "FINAL tree: map their branch ids with stale/trees/<hog>.rooted.nwk).\n")
    try:
        ok = apply(W, out, changed, new, opt["wait"], {"final": fp, "fm": fm, "text": text})
    except Exception:
        import traceback
        msg = traceback.format_exc()
        open(f"{out}/APPLY_ERROR", "a").write(f"{time.strftime('%F_%T')}\n{msg}")
        log(f"apply ERROR, APPLYING kept (recover with {out}/manifest.tsv):\n{msg}")
        return 1
    if ok and not os.listdir(tmp):
        os.rmdir(tmp)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
