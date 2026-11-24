"""Habitat (ambient) water temperature per species: OBIS records x WOA23 (2026-10-20).
Depth, pressure and temperature co-vary, so the deep-vs-shallow contrasts are re-tested with habitat
temperature as a covariate and against cold-water (polar) shallow species.
  extract (py310env, duckdb): OBIS dump records of the panel species (exact scientificname match to
          meta/species_codes.tsv 'species', CSV read as in 02_obis_depth.py), deduplicated to unique
          (species, 0.01 deg lat/lon, 1 m depth) rows with record count n and median OBIS bathymetry
          -> work/habitat/obis_panel_records.tsv
  temp    (envs/dc, h5py): WOA23 1991-2020 annual objectively analysed in-situ temperature t_an (1 deg,
          102 standard levels 0-5500 m; data/woa23/woa23_decav91C0_t00_01.nc). Each level is first
          extended into masked cells (land, or below the 1-deg cell-mean bottom) by up to 3 rings of 3x3
          means of valid cells, so coastal records and slope records deeper than the cell-mean bottom take
          the nearby open-water value at their depth; then linear interpolation between the bracketing
          levels. Below the deepest level with a value (incl. > 5500 m) that value is used (flagged).
          Record depth = OBIS depth if 0 <= depth < 11500 (as in 02_obis_depth.py); records without a
          depth are used only for T_lab (at the label depth = meta/species_codes.tsv depth_median).
Output meta/habitat_temp.tsv (one row per species code):
  T_obs_med/q10/q90, n_obs  over units = species x 1-deg cell x 10-m depth bin (records with a depth),
                            each unit weighted equally (limits intensively sampled sites); d_obs_med =
                            median unit depth (compare with depth_median)
  T_obs_med_bqc, n_obs_bqc  same, without records deeper than OBIS bathymetry + max(100 m, 20 %)
  frac_below, frac_ext      share of units below the WOA column / using extension-ring values
  T_lab_med, n_loc          over the unique 1-deg cells of all records, at the label depth
  T_glob                    global median of t_an at the label depth (original WOA cells)
  T_hab = T_obs_med if n_obs >= 5, else T_lab_med if n_loc >= 5, else T_glob if label depth >= 1000 m,
          else NA; T_src = obs / lab / glob_deep / none
  lat_med, abslat_med (1-deg cell centres); polar_shallow = label depth <= 200 m and
  (abslat_med >= 60 or T_hab < 4)
plus work/habitat/flags.tsv (row / record counts per case). Ambient water temperature only, not the
microhabitat of vents and seeps.
Usage: ~/anaconda3/envs/py310env/bin/python scripts/42_habitat_temp.py extract
       envs/dc/bin/python                   scripts/42_habitat_temp.py temp
"""
import csv
import os
import sys

SRC = os.path.expanduser("~/Global_Marine_Biodiversity_Records_All.csv")
WOA = "data/woa23/woa23_decav91C0_t00_01.nc"
CODES = "meta/species_codes.tsv"
OUTD = "work/habitat"
REC = f"{OUTD}/obis_panel_records.tsv"
OUT = "meta/habitat_temp.tsv"
FLAGS = f"{OUTD}/flags.tsv"


def extract():
    import duckdb
    os.makedirs(OUTD, exist_ok=True)
    names = sorted({r["species"].strip() for r in csv.DictReader(open(CODES), delimiter="\t")
                    if r["species"].strip()})
    con = duckdb.connect()
    con.execute(f"SET threads TO 24; SET memory_limit='64GB'; SET temp_directory='{OUTD}/duck_tmp';")
    con.execute("CREATE TABLE sp(name VARCHAR)")
    con.executemany("INSERT INTO sp VALUES (?)", [(n,) for n in names])
    con.execute(f"""
    COPY (
      WITH r AS (
        SELECT scientificname AS sp,
               TRY_CAST(latitude AS DOUBLE) AS lat, TRY_CAST(longitude AS DOUBLE) AS lon,
               TRY_CAST(depth AS DOUBLE) AS d, TRY_CAST(bathymetry AS DOUBLE) AS b
        FROM read_csv('{SRC}', header=true, all_varchar=true, parallel=true,
                      ignore_errors=true, strict_mode=false)
        WHERE scientificname IN (SELECT name FROM sp)
      )
      SELECT sp, round(lat, 2) AS lat, round(lon, 2) AS lon,
             CASE WHEN d >= 0 AND d < 11500 THEN round(d) END AS d,
             count(*) AS n, median(b) FILTER (WHERE b > 0) AS bathy
      FROM r WHERE lat BETWEEN -90 AND 90 AND lon BETWEEN -180 AND 360
      GROUP BY 1, 2, 3, 4 ORDER BY 1, 2, 3, 4
    ) TO '{REC}' (FORMAT csv, DELIMITER '\t', HEADER true);
    """)
    ns, nr, nn = con.execute(f"SELECT count(DISTINCT sp), count(*), sum(n) FROM "
                             f"read_csv('{REC}', delim='\t', header=true)").fetchone()
    print(f"extract done: {len(names)} panel names, {ns} with located records; {nr} unique rows, {nn} records")


def shift_lat(a, di):
    """a shifted along latitude (axis 1) without wrapping: out[:, i] = a[:, i - di]; edge rows 0 / False."""
    import numpy as np
    out = np.zeros_like(a)
    if di == 1:
        out[:, 1:] = a[:, :-1]
    else:
        out[:, :-1] = a[:, 1:]
    return out


def woa_fields():
    """z; tf = t_an extended by <= 3 rings per level, then filled down (and up) each column (NaN where a
    column has no value); ring (0 WOA value, 1-3 extension ring, -1 none); kmax (deepest level with a value)."""
    import h5py
    import numpy as np
    with h5py.File(WOA, "r") as f:
        assert "scale_factor" not in f["t_an"].attrs
        t = f["t_an"][0].astype("float64")
        z = f["depth"][:].astype("float64")
        lat, lon = f["lat"][:], f["lon"][:]
    assert t.shape == (len(z), 180, 360) and lat[0] == -89.5 and lon[0] == -179.5, (t.shape, lat[0], lon[0])
    t[~(t < 1e30)] = np.nan
    ring = np.where(np.isnan(t), -1, 0).astype(np.int8)
    for r in (1, 2, 3):
        v = ~np.isnan(t)
        x = np.where(v, t, 0.0)
        s = np.zeros_like(x)
        n = np.zeros_like(x)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di or dj:
                    xs, vs = np.roll(x, dj, axis=2), np.roll(v, dj, axis=2)   # longitude wraps
                    if di:
                        xs, vs = shift_lat(xs, di), shift_lat(vs, di)
                    s += xs
                    n += vs
        new = ~v & (n > 0)
        t[new] = s[new] / n[new]
        ring[new] = r
    v = ~np.isnan(t)
    kmax = np.where(v.any(axis=0), t.shape[0] - 1 - np.argmax(v[::-1], axis=0), -1)
    for k in range(1, t.shape[0]):            # below the deepest level with a value: that value
        m = np.isnan(t[k]) & ~np.isnan(t[k - 1])
        t[k][m] = t[k - 1][m]
    for k in range(t.shape[0] - 2, -1, -1):   # masked top levels (not expected): the shallowest value
        m = np.isnan(t[k]) & ~np.isnan(t[k + 1])
        t[k][m] = t[k + 1][m]
    return z, t, ring, kmax


def grid_ij(lat, lon):
    import numpy as np
    i = np.clip(np.floor(np.asarray(lat, float) + 90).astype(int), 0, 179)
    j = np.clip(np.floor(np.mod(np.asarray(lon, float) + 180, 360)).astype(int), 0, 359)
    return i, j


def interp(z, tf, ring, kmax, i, j, d):
    """t at depth d (m) in cells (i, j): linear between standard levels, depth clamped to 0-5500 m.
    Flags: below = deeper than the deepest level with a value; ext = an extension-ring value was used."""
    import numpy as np
    d = np.asarray(d, float)
    dc = np.clip(d, 0, z[-1])
    k = np.clip(np.searchsorted(z, dc, side="right") - 1, 0, len(z) - 2)
    w = (dc - z[k]) / (z[k + 1] - z[k])
    a, b = tf[k, i, j], tf[k + 1, i, j]
    val = a + w * (b - a)
    km = kmax[i, j]
    below = (km >= 0) & (d > z[np.maximum(km, 0)])
    ext = ~below & ((ring[k, i, j] > 0) | ((w > 0) & (ring[k + 1, i, j] > 0)))
    return val, below, ext


def temp():
    import numpy as np
    import pandas as pd
    z, tf, ring, kmax = woa_fields()
    for la, lo, d in ((0.3, -150.2, 0), (0.3, -150.2, 4000), (-65.2, 0.4, 0), (36.0, -30.0, 1000)):
        print(f"WOA check lat {la} lon {lo} {d} m: %.2f C" % interp(z, tf, ring, kmax, *grid_ij([la], [lo]), [d])[0][0])
    codes = pd.read_csv(CODES, sep="\t", dtype=str)
    codes["depth_median"] = pd.to_numeric(codes["depth_median"], errors="coerce")
    codes["species"] = codes["species"].str.strip()
    rec = pd.read_csv(REC, sep="\t", dtype={"sp": str})
    rec["i"], rec["j"] = grid_ij(rec["lat"], rec["lon"])
    dep = rec[rec["d"].notna()].copy()
    dep["T"], dep["below"], dep["ext"] = interp(z, tf, ring, kmax, dep["i"].to_numpy(), dep["j"].to_numpy(),
                                                dep["d"].to_numpy())
    b, dd = dep["bathy"].to_numpy(float), dep["d"].to_numpy(float)
    dep["bqc"] = ~((b > 0) & (dd > b + np.maximum(100.0, 0.2 * b)))
    nocell = dep["T"].isna().to_numpy()
    flags = {"rows": len(rec), "records": int(rec["n"].sum()), "species_with_rows": rec["sp"].nunique(),
             "rows_depth": len(dep), "records_depth": int(dep["n"].sum()),
             "rows_no_woa_value": int(nocell.sum()), "rows_below_column": int(dep["below"].sum()),
             "rows_gt5500m": int((dd > z[-1]).sum()), "rows_extended": int(dep["ext"].sum()),
             "rows_deeper_than_bathy": int((~dep["bqc"]).sum())}
    dep = dep[~nocell].copy()
    dep["dbin"] = (dep["d"] // 10).astype(int)

    def unit_t(x):
        return x.groupby(["sp", "i", "j", "dbin"], as_index=False).agg(
            T=("T", "mean"), d=("d", "mean"), below=("below", "max"), ext=("ext", "max"))
    units, units_q = unit_t(dep), unit_t(dep[dep["bqc"]])
    cells = rec[["sp", "i", "j"]].drop_duplicates()
    flags.update(units=len(units), units_bqc=len(units_q), cells=len(cells))
    gi, gj = [a.ravel() for a in np.meshgrid(np.arange(180), np.arange(360), indexing="ij")]
    tg = {}

    def tglob(d):   # global median at depth d over cells with WOA values at both bracketing levels
        if d not in tg:
            val = interp(z, tf, ring, kmax, gi, gj, np.full(gi.size, d))[0]
            k = int(np.clip(np.searchsorted(z, min(max(d, 0), z[-1]), side="right") - 1, 0, len(z) - 2))
            ok = (ring[k, gi, gj] == 0) & (ring[k + 1, gi, gj] == 0)
            tg[d] = float(np.median(val[ok])) if ok.any() else np.nan
        return tg[d]

    ug = {k: g for k, g in units.groupby("sp")}
    uq = {k: g for k, g in units_q.groupby("sp")}
    cg = {k: g for k, g in cells.groupby("sp")}
    out = []
    for r in codes.itertuples(index=False):
        sp, dl = r.species, float(r.depth_median)
        row = {"code": r.code, "species": sp, "depth_median": dl, "n_obs": 0, "n_obs_bqc": 0, "n_loc": 0}
        u = ug.get(sp)
        if u is not None:
            q = np.quantile(u["T"], [0.1, 0.5, 0.9])
            row.update(n_obs=len(u), T_obs_q10=q[0], T_obs_med=q[1], T_obs_q90=q[2], d_obs_med=u["d"].median(),
                       frac_below=u["below"].mean(), frac_ext=u["ext"].mean())
        u = uq.get(sp)
        if u is not None:
            row.update(n_obs_bqc=len(u), T_obs_med_bqc=float(np.median(u["T"])))
        c = cg.get(sp)
        if c is not None:
            lat_c = c["i"].to_numpy() - 89.5
            row.update(n_loc=len(c), lat_med=float(np.median(lat_c)), abslat_med=float(np.median(np.abs(lat_c))))
            if np.isfinite(dl):
                tl, bl, _ = interp(z, tf, ring, kmax, c["i"].to_numpy(), c["j"].to_numpy(), np.full(len(c), dl))
                okl = np.isfinite(tl)
                if okl.any():
                    row.update(T_lab_med=float(np.median(tl[okl])), frac_lab_below=float(bl[okl].mean()))
        row["T_glob"] = tglob(dl) if np.isfinite(dl) else np.nan
        if row["n_obs"] >= 5:
            row.update(T_hab=row["T_obs_med"], T_src="obs")
        elif row["n_loc"] >= 5 and np.isfinite(row.get("T_lab_med", np.nan)):
            row.update(T_hab=row["T_lab_med"], T_src="lab")
        elif np.isfinite(dl) and dl >= 1000:
            row.update(T_hab=row["T_glob"], T_src="glob_deep")
        else:
            row.update(T_hab=np.nan, T_src="none")
        out.append(row)
    o = pd.DataFrame(out)
    o["polar_shallow"] = ((o["depth_median"] <= 200) & ((o["abslat_med"] >= 60) | (o["T_hab"] < 4))).astype(int)
    cols = ["code", "species", "depth_median", "T_hab", "T_src", "n_obs", "T_obs_med", "T_obs_q10", "T_obs_q90",
            "d_obs_med", "n_obs_bqc", "T_obs_med_bqc", "frac_below", "frac_ext", "n_loc", "T_lab_med",
            "frac_lab_below", "T_glob", "lat_med", "abslat_med", "polar_shallow"]
    o.reindex(columns=cols).to_csv(OUT, sep="\t", index=False, float_format="%.6g", na_rep="NA")
    flags.update({f"T_src_{k}": int(v) for k, v in o["T_src"].value_counts().items()})
    flags["polar_shallow"] = int(o["polar_shallow"].sum())
    with open(FLAGS, "w") as f:
        f.write("".join(f"{k}\t{v}\n" for k, v in flags.items()))
    print("flags:", flags)
    deep = o["depth_median"] > 200
    print("median T_hab: deep %.2f C (n %d), shallow %.2f C (n %d)" % (
        o.loc[deep, "T_hab"].median(), deep.sum(), o.loc[~deep, "T_hab"].median(), (~deep).sum()))
    ok = o["T_obs_med"].notna() & o["T_lab_med"].notna()
    print("T_obs_med vs T_lab_med: r %.3f (n %d)" % (np.corrcoef(o.loc[ok, "T_obs_med"], o.loc[ok, "T_lab_med"])[0, 1], ok.sum()))
    x = (o["T_obs_med"] - o["T_lab_med"]).abs()
    print("largest |T_obs_med - T_lab_med|:")
    print(o.assign(dx=x).nlargest(8, "dx")[["code", "depth_median", "d_obs_med", "n_obs", "T_obs_med", "T_lab_med"]].to_string(index=False))
    x = (o["T_obs_med"] - o["T_obs_med_bqc"]).abs()
    print("|T_obs_med - T_obs_med_bqc|: median %.3f, max %.2f, > 1 C in %d species" % (x.median(), x.max(), (x > 1).sum()))
    print("wrote", OUT, len(o), "rows")


if __name__ == "__main__":
    {"extract": extract, "temp": temp}[sys.argv[1]]()
