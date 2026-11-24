#!/usr/bin/env python3
"""Autogroup nice of MY OWN long-running background sessions, so the final species tree (19h_final_tree.sh,
2026-10-19) gets the free CPU first. Under Linux autogroup scheduling every session is one scheduling entity
(weight 1024 at autogroup nice 0, 110 at nice 10), whatever the number of threads in it; raising the autogroup nice
of my daemons (OrthoFinder, HyPhy, eggNOG, genome QC) hands most of their share to the tree sessions without
stopping or slowing them when the CPU is free. A session is touched only if every process in it is mine (uid), runs
from or for ~/deepsea_conv and matches none of the protected patterns (the user's BCCMC / SDM / environment jobs);
other users' jobs never (the kernel would refuse anyway). Non-root writes are rate-limited by the kernel (one per
100 ms), hence the pause between writes.
  show SID...            members of each session, their autogroup and whether the session may be touched
  set W NICE SID...      set the autogroup nice (one write per autogroup), verify, log to W/autogroup.tsv
  restore W              write the logged old nice back (only to autogroups that still exist, same id)"""
import os, sys, time

ME = os.getuid()
ROOT = os.path.realpath(os.path.expanduser("~/deepsea_conv"))
PROTECT = ("bccmc", "sdm", "biotime", "maxent", "run_pipeline", "production_v5", "41_score", "40_bccmc",
           "00_make_env", "00b_fix_env", "marine_dd")


def procs():
    out = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            st = open(f"/proc/{d}/stat").read()
            f = st[st.rindex(")") + 2:].split()
            out[int(d)] = {"sid": int(f[3]), "uid": os.stat(f"/proc/{d}").st_uid}
        except (OSError, ValueError, IndexError):
            pass
    return out


def info(pid):
    def rd(fn):
        try:
            with open(fn, "rb") as fh:
                return fh.read().replace(b"\0", b" ").decode(errors="replace").strip()
        except OSError:
            return ""
    try:
        cwd = os.readlink(f"/proc/{pid}/cwd")
    except OSError:
        cwd = ""
    ag = rd(f"/proc/{pid}/autogroup")              # "/autogroup-12345 nice 0"
    return rd(f"/proc/{pid}/cmdline"), cwd, ag


def parse_ag(ag):
    try:
        a, _, n = ag.split()
        return a, int(n)
    except ValueError:
        return None, None


def members(sid, ps):
    return sorted(p for p, v in ps.items() if v["sid"] == sid)


def check(sid, ps):
    """(ok, reason, rows) for one session; rows = (pid, cmdline, cwd, autogroup id, nice)."""
    ms = members(sid, ps)
    if not ms:
        return False, "no processes", []
    rows, bad = [], []
    for p in ms:
        cmd, cwd, ag = info(p)
        a, n = parse_ag(ag)
        rows.append((p, cmd, cwd, a, n))
        if ps[p]["uid"] != ME:
            bad.append(f"{p}: uid {ps[p]['uid']}")
        elif any(k in cmd.lower() or k in cwd.lower() for k in PROTECT):
            bad.append(f"{p}: protected pattern ({cmd[:80]})")
        elif not (os.path.realpath(cwd or "/").startswith(ROOT) or "deepsea_conv" in cmd):
            bad.append(f"{p}: not a deepsea_conv process (cwd {cwd})")
        elif a is None:
            bad.append(f"{p}: no autogroup")
    return (not bad), "; ".join(bad) or "ok", rows


def show(sids):
    ps = procs()
    for s in sids:
        ok, why, rows = check(s, ps)
        ags = sorted({(r[3], r[4]) for r in rows if r[3]})
        print(f"session {s}: {len(rows)} processes, autogroups {ags}, {'OK' if ok else 'REFUSED'}: {why}")
        for r in rows[:6]:
            print(f"   {r[0]} {r[3]} nice {r[4]} cwd {r[2]} :: {r[1][:150]}")
        if len(rows) > 6:
            print(f"   ... {len(rows) - 6} more")


def write_nice(pid, a, nice):
    with open(f"/proc/{pid}/autogroup", "w") as fh:
        fh.write(str(nice))
    time.sleep(0.2)
    a2, n2 = parse_ag(info(pid)[2])
    return a2 == a and n2 == nice


def set_(W, nice, sids):
    ps = procs()
    log = open(f"{W}/autogroup.tsv", "a")
    if log.tell() == 0:
        log.write("sid\tpid\tautogroup\told_nice\tnew_nice\ttime\tcmdline\n")
    for s in sids:
        ok, why, rows = check(s, ps)
        if not ok:
            print(f"session {s} refused: {why}")
            continue
        done = set()
        for p, cmd, cwd, a, n in rows:
            if a in done:
                continue
            done.add(a)
            if n == nice:
                print(f"session {s} {a}: already nice {nice}")
                continue
            try:
                good = write_nice(p, a, nice)
            except OSError as e:
                print(f"session {s} {a} pid {p}: write failed ({e})")
                continue
            print(f"session {s} {a} (pid {p}): nice {n} -> {nice} {'verified' if good else 'NOT verified'}")
            if good:
                log.write(f"{s}\t{p}\t{a}\t{n}\t{nice}\t{time.strftime('%F %T')}\t{cmd[:200]}\n")
                log.flush()
    log.close()


def restore(W):
    fn = f"{W}/autogroup.tsv"
    if not os.path.exists(fn):
        print("nothing to restore")
        return
    rows = [l.rstrip("\n").split("\t") for l in list(open(fn))[1:] if l.strip() and not l.startswith("#")]
    old = {}
    for r in rows:                       # first logged old nice of each autogroup = the original value
        old.setdefault(r[2], (int(r[0]), int(r[3])))
    ps = procs()
    live = {}
    for p, v in ps.items():
        if v["uid"] != ME:
            continue
        a, n = parse_ag(info(p)[2])
        if a in old and a not in live:
            live[a] = (p, n)
    for a, (s, n0) in old.items():
        if a not in live:
            print(f"{a} (session {s}): gone")
            continue
        p, n = live[a]
        if n == n0:
            print(f"{a} (session {s}): already nice {n0}")
            continue
        try:
            good = write_nice(p, a, n0)
        except OSError as e:
            print(f"{a} (session {s}) pid {p}: write failed ({e})")
            continue
        print(f"{a} (session {s}, pid {p}): nice {n} -> {n0} {'verified' if good else 'NOT verified'}")
    with open(fn, "a") as log:
        log.write(f"#restored\t-\t-\t-\t-\t{time.strftime('%F %T')}\t-\n")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "show":
        show([int(x) for x in sys.argv[2:]])
    elif cmd == "set":
        set_(sys.argv[2], int(sys.argv[3]), [int(x) for x in sys.argv[4:]])
    elif cmd == "restore":
        restore(sys.argv[2])
    else:
        sys.exit(f"unknown command {cmd}")
