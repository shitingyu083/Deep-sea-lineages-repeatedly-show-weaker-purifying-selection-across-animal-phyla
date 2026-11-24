#!/usr/bin/env Rscript
# A6 timing (PREREG_addon.md): did the extant deep-sea (> 200 m) lineages colonise the deep sea mainly after the last
# global deep-water anoxic / extinction events (recruitment from shallow refugia; Jacobs & Lindberg 1998, PNAS 95:9396:
# after the Cenomanian-Turonian OAE2 ~94 Ma and the latest-Paleocene event ~56 Ma), after the mid-Cenozoic cooling
# (~34 Ma, the young-fauna variant tested by Lins et al. 2012), or throughout the Phanerozoic (ancient origin; Thuy et
# al. 2012, Lins et al. 2012)?  Usage: Rscript 58_origin_timing.R [TREE] [OUTDIR]
#   (defaults work/final_v3/date/lsd.timetree.nwk, work/addon/A6)
# env LABEL_COL (asr), MA (1 = TREE is the LSD2 tree, rescaled to Ma by addon_lib load_tree), WIDTH (10), EPOCHS
#   (94,56,33.9), NBOOT (200), NPOW (50; 0 = no power module), PWIN (event ages of the planted windows, default all 9),
#   PMULT (5,10,30), PEMULT (2,4,8), NMAP (200), NDRAW (10000), SEED (20261025), NCORES (8), ROBUST (multi-tree newick
#   in Ma, default work/mcmctree_bb/interp/draws.nwk, skipped if absent), NROB (0 = all), POINT (comma list of further
#   Ma trees, default interp/point_correlated.nwk, skipped if absent), OUTGROUP (Monbre,Salros)
# Labels: meta/labels_v2.tsv LABEL_COL, D / S tips only (outgroups dropped); Mk2 = addon_lib (ARD, ER fallback), per Myr.
# (1) Decision tests, directional (1 df): the S->D rate inside an interval is free, D->S and the S->D rate outside are
#     shared as in the homogeneous model (ER null: outside S->D = D->S). Families: windows [E, E - WIDTH] Ma at the 9
#     events; epochs [B, 0] Ma (extant deep fauna recruited after B) at EPOCHS. z = sign(log ratio) x sqrt(LR), p1 =
#     one-sided p for raised colonisation, Holm within family. Parametric bootstrap (NBOOT labellings under the
#     homogeneous fit, deep count within 10 %, all models refitted): p_boot per interval and p_fw = single-step max-z
#     within family. Effect = Holm p < 0.05 and p_fw < 0.05 on the main tree; robust = Holm p < 0.05 in >= 80 % of the
#     ROBUST trees; rates on an optimiser bound are flagged.
# (1a) Attribution (descriptive, per family): the effect interval with the largest z is the attributed interval; every
#     other effect interval of the family is refitted conditional on it (own S->D rate on its part outside the
#     attributed interval, compared with the background, or, if nested inside it, on itself, compared with the rest of
#     the attributed interval; 1 df, one-sided asymptotic p, Holm within family) and reported as a separate effect only
#     if its conditional Holm p < 0.05, else as explained by the attributed interval (nested epochs [B, 0] all contain
#     a younger rise; adjacent windows share long branches). frac_attr = share of the ROBUST trees in which the
#     attributed interval has the largest z of its family.
# (1b) Secondary (ED): 4-parameter era model (addon_lib mk2_era_fit, own q_SD and q_DS in the window), asymptotic LRT
#     df = 4 - npar, Holm over the 9 windows; flat likelihoods can give degenerate optima (reported as they are).
# (1c) Power and localisation (prespecified, descriptive): per planted interval x multiplier (x PMULT in a PWIN window,
#     x PEMULT after an epoch boundary) NPOW labellings with the background S->D rate rescaled so that E[n_D] = observed,
#     conditioned on the deep count as in the bootstrap; detection = rule (1) for the planted interval (p_fw against
#     the bootstrap null of the data); localisation = the planted interval has the largest z of its family. An effect
#     is attributed to its interval only if det_loc (detected and localised) >= 0.8 at the multiplier nearest the
#     observed ratio, else it is reported as not separable from the adjacent intervals (tip states cannot place a
#     change within a branch, and intervals share long branches).
# (2) NMAP stochastic maps under the homogeneous fit (node states from the joint posterior by pre-order sampling, edge
#     histories conditional on the end states by rejection sampling with the first change forced when they differ,
#     Nielsen 2002): per window S->D changes, shallow time, expected S->D changes = q_SD x shallow time. Check: over
#     the whole tree changes / (rate x time) = 1 at the ML rates (EM fixed point), up to Monte Carlo error.
# (3) Origins = maximal all-D clades (addon_lib deep_clades), origin age uniform on the stem branch (NDRAW draws). Per
#     window: origins, expected = N x shallow time in the window / total shallow time (constant colonisation per
#     shallow lineage-time; shallow = outside the D clades and the stem above the drawn age), ratio R = origins /
#     expected. Random-window null: per draw 9 non-overlapping WIDTH windows, older bounds uniform in
#     [WIDTH, min(541, root age)] Ma; p = (1 + #R_random >= R) / (1 + n) per window (the 9 random windows of the same
#     draw) and pooled over the 9 windows. Predate table: per epoch boundary B, deep clades whose crown is older than B
#     (supported = P(D) at the crown >= 0.95 in the maps; a minimum count of deep lineages older than B), crown counts
#     also on every ROBUST tree. Descriptive (no decision role).
# (4) Robustness: (1) without the bootstrap on every ROBUST and POINT tree.
# Caveats: from extant tips a raised recent S->D rate cannot be separated from extinction of old deep lineages
# (state-dependent extinction); both read as recent recruitment of the extant deep fauna. > 200 m includes the upper
# bathyal. Outputs (plain tsv): dir_main, dir_boot, era4_main, power, power_summary, smap_windows, origins,
# origin_windows, predate, dir_robust, dir_point, attribution, decision (.tsv), report.txt
fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(Sys.getenv("ADDON_LIB", file.path(dirname(normalizePath(sub("^--file=", "", fa[1]))), "addon_lib.R")))
suppressPackageStartupMessages(library(parallel))
a <- commandArgs(TRUE); ev <- function(k, d) Sys.getenv(k, d)
TREE <- if (length(a) >= 1) a[1] else file.path(W, "date/lsd.timetree.nwk")
OUT <- if (length(a) >= 2) a[2] else file.path(R, "work/addon/A6")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LC <- ev("LABEL_COL", "asr"); WIDTH <- as.numeric(ev("WIDTH", "10")); NBOOT <- as.integer(ev("NBOOT", "200"))
NMAP <- as.integer(ev("NMAP", "200")); NDRAW <- as.integer(ev("NDRAW", "10000")); SEED <- as.integer(ev("SEED", "20261025"))
NC <- as.integer(ev("NCORES", "8")); NROB <- as.integer(ev("NROB", "0")); BB <- file.path(R, "work/mcmctree_bb/interp")
MA <- ev("MA", "1") == "1"; OG <- strsplit(ev("OUTGROUP", "Monbre,Salros"), ",")[[1]]
ROB <- ev("ROBUST", file.path(BB, "draws.nwk")); PTS <- strsplit(ev("POINT", file.path(BB, "point_correlated.nwk")), ",")[[1]]
EVT <- data.table(event = c("end-Permian", "end-Triassic", "Toarcian OAE", "OAE1a", "OAE2", "K-Pg", "PETM",
                            "Eocene-Oligocene", "mid-Miocene"), E = c(252, 201.4, 183, 120, 94, 66, 56, 33.9, 14))
nv <- function(k, d) as.numeric(strsplit(ev(k, d), ",")[[1]])
EPO <- nv("EPOCHS", "94,56,33.9"); NPOW <- as.integer(ev("NPOW", "50")); PWIN <- nv("PWIN", paste(EVT$E, collapse = ","))
PMULT <- nv("PMULT", "5,10,30"); PEMULT <- nv("PEMULT", "2,4,8")
if (!all(vapply(PWIN, function(v) any(abs(EVT$E - v) < 1e-6), NA))) stop("PWIN must be event ages of EVT")
TST <- rbind(EVT[, .(family = "window", event, E, width = WIDTH)],
             data.table(family = "epoch", event = paste0("after ", EPO, " Ma"), E = EPO, width = EPO))
lg <- function(...) { cat(format(Sys.time(), "%F %T"), ..., "\n"); flush(stdout()) }
mcl <- function(X, f) if (NC > 1 && .Platform$OS.type != "windows") mclapply(X, f, mc.cores = NC, mc.preschedule = FALSE) else lapply(X, f)
lb <- fread(file.path(R, "meta/labels_v2.tsv"), sep = "\t", quote = ""); stopifnot(!anyDuplicated(lb$code))
if (!LC %in% names(lb)) stop("meta/labels_v2.tsv has no column ", LC)
lab <- setNames(as.character(lb[[LC]]), lb$code)
prep <- function(tr) {                                          # D / S Metazoa tips only; x = 1 for D
  tr <- keep.tip(tr, setdiff(intersect(tr$tip.label, names(lab)[lab %in% c("D", "S")]), OG))
  if (!is.rooted(tr)) stop("tree is not rooted")
  d <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]               # ages are read as root age - depth
  if (diff(range(d)) > 1e-3 * max(d)) stop("tree is not ultrametric (tip depths ", signif(min(d), 6), " to ", signif(max(d), 6), ")")
  tr$edge.length[tr$edge.length <= 0] <- 1e-6
  list(tr = tr, x = setNames(as.integer(lab[tr$tip.label] == "D"), tr$tip.label))
}
# ---- (1) directional tests ------------------------------------------------------------------------------------------
on_bd <- function(q, bd) any(log(q) < bd[["lo"]] + 0.01 | log(q) > bd[["hi"]] - 0.01, na.rm = TRUE)
# own S->D rate inside [E, E - width] Ma; D->S and the outside S->D rate shared (ER null: outside S->D = D->S); 1 df.
# The m = 1 start is the null itself, so LR >= 0 up to the optimiser tolerance
dir_fit <- function(pp, bd, E, width, f0) {
  seg <- era_segments(pp$tr, E, width); wb <- sum(seg[, 2]); lo <- bd[["lo"]]; hi <- bd[["hi"]]
  if (!(wb > 0)) return(data.table(win_bl = wb, q_SD_bg = f0$q_SD, q_DS = f0$q_DS, q_SD_in = NA_real_, ratio = NA_real_,
                                   LR = 0, z = 0, p1 = 0.5, at_bound = on_bd(c(f0$q_SD, f0$q_DS), bd), conv = 0L))
  er <- f0$model == "ER"; rt <- function(l) if (er) c(l[1], l[1], l[2]) else l
  nll <- function(l) { q <- exp(rt(l)); v <- -mk2_era_ll(pp, seg, q[1], q[2], q[3], q[2]); if (is.finite(v)) v else 1e300 }
  b0 <- if (er) log(f0$q_SD) else log(c(f0$q_SD, f0$q_DS)); best <- NULL
  for (m in c(1, 3, 10, 30, 100, 1 / 3)) {
    o <- tryCatch(optim(pmin(pmax(c(b0, log(f0$q_SD * m)), lo), hi), nll, method = "L-BFGS-B", lower = lo, upper = hi),
                  error = function(e) NULL)
    if (!is.null(o) && o$value < 1e300 && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(data.table(win_bl = wb, q_SD_bg = NA_real_, q_DS = NA_real_, q_SD_in = NA_real_, ratio = NA_real_,
                                       LR = NA_real_, z = NA_real_, p1 = NA_real_, at_bound = NA, conv = NA_integer_))
  q <- exp(rt(best$par)); lr <- max(0, 2 * (-best$value - f0$loglik)); z <- sign(log(q[3] / q[1])) * sqrt(lr)
  data.table(win_bl = wb, q_SD_bg = q[1], q_DS = q[2], q_SD_in = q[3], ratio = q[3] / q[1], LR = lr, z = z,
             p1 = pnorm(z, lower.tail = FALSE), at_bound = on_bd(q, bd), conv = best$convergence)
}
dir_all <- function(tr, x, f0, par = FALSE) {           # every row of TST; Holm within family
  pp <- mk2_prep(tr, x); bd <- mk2_bounds(tr)
  f <- function(k) dir_fit(pp, bd, TST$E[k], TST$width[k], f0)
  r <- cbind(TST, rbindlist(if (par) mcl(seq_len(nrow(TST)), f) else lapply(seq_len(nrow(TST)), f)))
  r[, p_holm := p.adjust(p1, "holm"), by = family][]
}
# (1a) conditional fit of interval o given the attributed interval a: S->D rate classes 1 = background, 2 = inside a,
# 3 = the part of o outside a (o not nested in a; reference class 1) or o itself (o nested in a; reference class 2);
# null = classes 1 and 2 only (the dir_fit model of a, started at its optimum s0), D->S shared (ER: class 1 = D->S)
cond_fit <- function(pp, bd, Ea, wa, Eo, wo, f0, s0) {
  dep <- node.depth.edgelength(pp$tr); H <- max(dep[seq_len(pp$n)]); t0 <- dep[pp$tr$edge[, 1]]; t1 <- dep[pp$tr$edge[, 2]]
  cu <- sort(unique(pmin(pmax(c(H, Ea, Ea - wa, Eo, Eo - wo, 0), 0), H)), decreasing = TRUE); J <- length(cu) - 1L
  S <- vapply(seq_len(J), function(j) pmax(0, pmin(t1, H - cu[j + 1L]) - pmax(t0, H - cu[j])), numeric(length(t0)))
  mid <- (cu[-1] + cu[-length(cu)]) / 2; ia <- mid <= Ea & mid >= Ea - wa; io <- mid <= Eo & mid >= Eo - wo
  nest <- !any(io & !ia); fr <- if (nest) io & ia else io & !ia; cl0 <- ifelse(ia, 2L, 1L); cl1 <- replace(cl0, fr, 3L)
  er <- f0$model == "ER"; o <- if (er) 1L else 2L; lo <- bd[["lo"]]; hi <- bd[["hi"]]; ir <- if (nest) o + 1L else 1L
  out <- data.table(nested = nest, free_bl = sum(S[, fr]), q_ref = NA_real_, q_free = NA_real_, ratio = NA_real_,
                    ll0 = NA_real_, ll1 = NA_real_, LR = NA_real_, z = NA_real_, p1 = NA_real_, at_bound = NA, conv = NA_integer_)
  if (!(out$free_bl > 0)) return(out[, `:=`(LR = 0, z = 0, p1 = 0.5)][])
  nll <- function(l, cl) {
    qs <- exp(c(l[1], l[o + 1L], if (length(l) > o + 1L) l[o + 2L] else NA_real_)); b <- exp(l[o])
    P <- mk2_P(S[, 1], qs[cl[1]], b); for (j in seq_len(J)[-1]) P <- mm2(P, mk2_P(S[, j], qs[cl[j]], b))
    v <- -mk2_ll(pp, P = P); if (is.finite(v)) v else 1e300
  }
  fit <- function(starts, cl) {
    best <- NULL
    for (s in starts) {
      r <- tryCatch(optim(pmin(pmax(s, lo), hi), function(l) nll(l, cl), method = "L-BFGS-B", lower = lo, upper = hi),
                    error = function(e) NULL)
      if (!is.null(r) && r$value < 1e300 && (is.null(best) || r$value < best$value)) best <- r
    }
    best
  }
  b0 <- function(q) log(if (er) q[c(1, 3)] else q)
  n0 <- fit(list(b0(c(s0$q_SD_bg, s0$q_DS, s0$q_SD_in)), b0(c(f0$q_SD, f0$q_DS, f0$q_SD))), cl0)
  if (is.null(n0)) return(out)
  n1 <- fit(lapply(c(1, 3, 10, 30, 100, 1 / 3, 1 / 10), function(m) c(n0$par, n0$par[ir] + log(m))), cl1)
  if (is.null(n1)) return(out)
  lr <- max(0, 2 * (n0$value - n1$value)); q <- exp(n1$par); rt <- q[o + 2L] / q[ir]
  out[, `:=`(q_ref = q[ir], q_free = q[o + 2L], ratio = rt, ll0 = -n0$value, ll1 = -n1$value, LR = lr,
             z = sign(log(rt)) * sqrt(lr), p1 = pnorm(sign(log(rt)) * sqrt(lr), lower.tail = FALSE), at_bound = on_bd(q, bd),
             conv = n1$convergence)][]
}
# parametric bootstrap under the homogeneous fit: every labelling refitted (homogeneous + all directional tests);
# p_boot per interval, p_fw = single-step max-z within family
boot_dir <- function(tr, x, f0, obs) {
  sm <- mk2_sim_cond(tr, f0$q_SD, f0$q_DS, NBOOT, sum(x), seed = SEED + 1000L)
  lg("bootstrap:", ncol(sm$sims), "labellings, deep-count tolerance", sm$tol, ", acceptance", signif(sm$acc_rate, 3))
  B <- rbindlist(Filter(is.data.table, mcl(seq_len(ncol(sm$sims)), function(j) {
    xs <- setNames(as.integer(sm$sims[, j]), rownames(sm$sims))
    f <- tryCatch(mk2_fit(tr, xs, "ARD"), error = function(e) NULL)
    r <- if (is.null(f)) NULL else tryCatch(dir_all(tr, xs, f), error = function(e) NULL)
    if (is.null(r) || nrow(r) != nrow(TST) || anyNA(r$z)) NULL else r[, .(sim = j, family, event, E, z, LR, at_bound, model_0 = f$model)]
  })))
  o <- list(B = B, n = 0L, tol = sm$tol, acc = sm$acc_rate, mx = NULL, p_boot = rep(NA_real_, nrow(obs)), p_fw = rep(NA_real_, nrow(obs)))
  if (!nrow(B)) return(o)
  o$n <- ns <- uniqueN(B$sim); o$mx <- B[, .(mz = max(z)), by = .(sim, family)]
  fam <- obs$family; EE <- obs$E; zz <- obs$z
  o$p_boot <- vapply(seq_along(zz), function(k) (1 + sum(B$z[B$family == fam[k] & abs(B$E - EE[k]) < 1e-9] >= zz[k] - 1e-9)) / (1 + ns), 0)
  o$p_fw <- vapply(seq_along(zz), function(k) (1 + sum(o$mx$mz[o$mx$family == fam[k]] >= zz[k] - 1e-9)) / (1 + ns), 0)
  o
}
era4_all <- function(tr, x, f0) {                       # secondary (ED): 4-parameter era model, asymptotic LRT
  bd <- mk2_bounds(tr)
  r <- rbindlist(Filter(is.data.table, mcl(EVT$E, function(E) tryCatch(mk2_era_fit(tr, x, E, WIDTH, f0), error = function(e) NULL))))
  if (!nrow(r)) return(NULL)
  r[win_bl <= 0, c("q_SD_win", "q_DS_win", "ratio_SD", "ratio_DS") := NA_real_]
  r[, at_bound := vapply(seq_len(.N), function(k) on_bd(c(q_SD_bg[k], q_DS_bg[k], q_SD_win[k], q_DS_win[k]), bd), NA)]
  r[, `:=`(event = EVT$event[match(E, EVT$E)], p_holm = p.adjust(p, "holm"))]; setcolorder(r, "event"); r[]
}
# ---- (1c) power and localisation -------------------------------------------------------------------------------------
# labellings (tips x nsim, 1 = D) with S->D = aw inside [E, E - width] Ma, (a, b) elsewhere; tc in cladewise order;
# root from the stationary distribution of the background
era_sim <- function(tc, a, b, aw, E, width, nsim) {
  n <- Ntip(tc); seg <- era_segments(tc, E, width)
  P <- mm2(mm2(mk2_P(seg[, 1], a, b), mk2_P(seg[, 2], aw, b)), mk2_P(seg[, 3], a, b))
  X <- matrix(0L, n + tc$Nnode, nsim); X[n + 1L, ] <- as.integer(runif(nsim) < a / (a + b))
  for (i in seq_len(nrow(tc$edge))) { p <- tc$edge[i, 1]; cn <- tc$edge[i, 2]
    X[cn, ] <- as.integer(runif(nsim) < ifelse(X[p, ] == 1L, P$DD[i], P$SD[i])) }
  out <- X[seq_len(n), , drop = FALSE]; rownames(out) <- tc$tip.label; out
}
# background scale s (S->D = s x a outside, s x a x mult inside) with E[n_D] = nD (common random numbers, so the
# mean deep count is monotone in s)
calib <- function(tc, a, b, mult, E, width, nD, nsim = 2000L, seed = 1L) {
  f <- function(ls) { set.seed(seed); mean(colSums(era_sim(tc, exp(ls) * a, b, exp(ls) * a * mult, E, width, nsim))) - nD }
  lo <- log(1e-3); hi <- log(10); flo <- f(lo); fhi <- f(hi)
  if (flo >= 0) return(exp(lo))
  if (fhi <= 0) return(exp(hi))
  exp(uniroot(f, c(lo, hi), f.lower = flo, f.upper = fhi, tol = 0.02)$root)
}
pow_cond <- function(tc, a, b, aw, E, width, nsim, nD, tol = 0.1, tol2 = 0.2, batch = 2000L, max_draw = 2e5, seed = 1L) {
  set.seed(seed); n <- Ntip(tc)
  okf <- function(X, tl) { k <- colSums(X); abs(k - nD) <= tl * nD & k >= 5 & n - k >= 5 }
  X <- era_sim(tc, a, b, aw, E, width, batch); nd <- batch; tl <- if (mean(okf(X, tol)) < 0.01) tol2 else tol
  keep <- X[, okf(X, tl), drop = FALSE]
  while (ncol(keep) < nsim && nd < max_draw) { X <- era_sim(tc, a, b, aw, E, width, batch); nd <- nd + batch
    keep <- cbind(keep, X[, okf(X, tl), drop = FALSE]) }
  list(sims = keep[, seq_len(min(nsim, ncol(keep))), drop = FALSE], tol = tl, acc_rate = ncol(keep) / nd)
}
# per planted cell: detection by rule (1) (p_fw against the bootstrap max-z null of the data), localisation = the
# planted interval has the largest z of its family; argmax kept for the confusion matrix (ED)
power_run <- function(tr, x, f0, bt) {
  tc <- reorder(tr, "cladewise"); nD <- sum(x)
  cells <- rbind(CJ(family = "window", E = PWIN, mult = PMULT), CJ(family = "epoch", E = EPO, mult = PEMULT))
  cells[, width := fifelse(family == "window", WIDTH, E)]
  mxw <- bt$mx[family == "window", mz]; mxe <- bt$mx[family == "epoch", mz]
  rbindlist(lapply(seq_len(nrow(cells)), function(k) {
    cc <- cells[k]; s <- calib(tc, f0$q_SD, f0$q_DS, cc$mult, cc$E, cc$width, nD, seed = SEED + 3000L + k)
    ps <- pow_cond(tc, s * f0$q_SD, f0$q_DS, s * f0$q_SD * cc$mult, cc$E, cc$width, NPOW, nD, seed = SEED + 4000L + k)
    lg("power", cc$family, cc$E, "x", cc$mult, ": background scale", signif(s, 3), ",", ncol(ps$sims), "labellings, acceptance",
       signif(ps$acc_rate, 3))
    if (!ncol(ps$sims)) return(NULL)
    rbindlist(Filter(is.data.table, mcl(seq_len(ncol(ps$sims)), function(j) {
      xs <- setNames(as.integer(ps$sims[, j]), rownames(ps$sims))
      f <- tryCatch(mk2_fit(tr, xs, "ARD"), error = function(e) NULL)
      r <- if (is.null(f)) NULL else tryCatch(dir_all(tr, xs, f), error = function(e) NULL)
      if (is.null(r) || anyNA(r$z)) return(NULL)
      fam <- r$family; zz <- r$z; ph <- r$p_holm
      pfw <- vapply(seq_along(zz), function(m) { mx <- if (fam[m] == "window") mxw else mxe
        (1 + sum(mx >= zz[m] - 1e-9)) / (1 + length(mx)) }, 0)
      inf <- which(fam == cc$family); i0 <- inf[abs(r$E[inf] - cc$E) < 1e-6]; am <- inf[which.max(zz[inf])]
      data.table(family = cc$family, E = cc$E, mult = cc$mult, bg_scale = s, tol = ps$tol, acc_rate = ps$acc_rate, sim = j,
                 nD = sum(xs), z = zz[i0], p_holm = ph[i0], p_fw = pfw[i0], detected = ph[i0] < 0.05 & pfw[i0] < 0.05,
                 argmax = r$event[am], local = am == i0, any_sig = any(ph[inf] < 0.05 & pfw[inf] < 0.05))
    })))
  }))
}
pw_sum <- function(P) P[, .(n = .N, nD_mean = mean(nD), bg_scale = bg_scale[1], acc_rate = acc_rate[1], tol = tol[1],
                           power = mean(detected), loc = mean(local), det_loc = mean(detected & local), any_sig = mean(any_sig)),
                        by = .(family, E, mult)]
# ---- (2) stochastic maps ------------------------------------------------------------------------------------------
mk2_partials <- function(pp, a, b) {                     # conditional likelihoods (postorder), rows scaled to max 1
  P <- mk2_P(pp$tr$edge.length, a, b); L <- pp$L0; e <- pp$tr$edge
  for (i in seq_len(nrow(e))) {
    l <- L[e[i, 2], ] / max(L[e[i, 2], ]); p <- e[i, 1]
    L[p, ] <- L[p, ] * c(P$SS[i] * l[1] + P$SD[i] * l[2], P$DS[i] * l[1] + P$DD[i] * l[2])
  }
  L / pmax(L[, 1], L[, 2])
}
# jump times (from the top of the edge) of a path from state i to state j (0 = S, 1 = D) over time t: rejection
# sampling with the first change forced when i != j (Nielsen 2002)
path_jumps <- function(i, j, t, a, b, cap = 1e5) {
  q <- c(a, b)
  for (k in seq_len(cap)) {
    s <- i; tt <- 0; out <- numeric(0)
    if (i != j) { tt <- -log1p(-runif(1) * (-expm1(-q[i + 1] * t))) / q[i + 1]; out <- tt; s <- 1L - i }
    repeat { tt <- tt + rexp(1, q[s + 1]); if (tt >= t) break; out <- c(out, tt); s <- 1L - s }
    if (s == j) return(out)
  }
  stop("path_jumps: no path accepted in ", cap, " attempts")
}
# nmap maps: node states from the joint posterior (root prior 1/2, 1/2 as mk2_ll; pre-order sampling), then edge
# histories; seg = (hi, lo, state) age spans, jmp = (age, from-state) of every change
smap <- function(tr, x, f0, nmap, seed) {
  set.seed(seed); a <- f0$q_SD; b <- f0$q_DS
  pp <- mk2_prep(tr, x); tp <- pp$tr; n <- pp$n; e <- tp$edge; bl <- tp$edge.length; ne <- nrow(e)
  L <- mk2_partials(pp, a, b); P <- mk2_P(bl, a, b); dep <- node.depth.edgelength(tp); age <- max(dep[seq_len(n)]) - dep
  ST <- matrix(0L, n + tp$Nnode, nmap); r <- n + 1L
  ST[r, ] <- as.integer(runif(nmap) < L[r, 2] / (L[r, 1] + L[r, 2]))
  for (i in rev(seq_len(ne))) {                                           # reverse postorder = pre-order
    up <- ST[e[i, 1], ] == 1L; cn <- e[i, 2]
    w1 <- ifelse(up, P$DD[i], P$SD[i]) * L[cn, 2]; w0 <- ifelse(up, P$DS[i], P$SS[i]) * L[cn, 1]
    ST[cn, ] <- as.integer(runif(nmap) < w1 / (w0 + w1))
  }
  A0 <- age[e[, 1]]; A1 <- age[e[, 2]]
  maps <- lapply(seq_len(nmap), function(m) {
    s0 <- ST[e[, 1], m]; J <- Map(path_jumps, s0, ST[e[, 2], m], bl, MoreArgs = list(a = a, b = b)); k <- lengths(J)
    list(seg = cbind(hi = unlist(Map(function(A, j) c(A, A - j), A0, J)), lo = unlist(Map(function(A, j, B) c(A - j, B), A0, J, A1)),
                     s = unlist(Map(function(s, kk) (s + 0:kk) %% 2L, s0, k))),
         jmp = matrix(as.numeric(c(unlist(Map(`-`, A0, J)), unlist(Map(function(s, kk) (s + seq_len(kk) - 1L) %% 2L, s0, k)))),
                      ncol = 2, dimnames = list(NULL, c("age", "from"))))
  })
  list(tr = tp, ST = ST, maps = maps)
}
qn <- function(v, p) { v <- v[is.finite(v)]; if (length(v)) unname(quantile(v, p)) else NA_real_ }
# per map and window (age span (lo, hi]): time in S / D and S->D / D->S changes; "all" = whole tree
WIN <- data.table(window = c(EVT$event, "all"), E = c(EVT$E, NA), lo = c(EVT$E - WIDTH, -Inf), hi = c(EVT$E, Inf))
win_tab <- function(sm) rbindlist(lapply(seq_along(sm$maps), function(m) {
  g <- sm$maps[[m]]$seg; j <- sm$maps[[m]]$jmp; s0 <- g[, "s"] == 0
  rbindlist(lapply(seq_len(nrow(WIN)), function(k) {
    o <- pmax(0, pmin(g[, "hi"], WIN$hi[k]) - pmax(g[, "lo"], WIN$lo[k])); w <- j[, "age"] > WIN$lo[k] & j[, "age"] <= WIN$hi[k]
    data.table(map = m, window = WIN$window[k], T_S = sum(o[s0]), T_D = sum(o[!s0]), N_SD = sum(w & j[, "from"] == 0),
               N_DS = sum(w & j[, "from"] == 1))
  }))
}))
smap_summary <- function(M, f0) {                       # expected changes = rate x time in the source state
  M <- copy(M)[, `:=`(X_SD = f0$q_SD * T_S, X_DS = f0$q_DS * T_D)]
  M <- rbind(M, M[window %in% EVT$event, .(window = "pooled_9", T_S = sum(T_S), T_D = sum(T_D), N_SD = sum(N_SD),
                                            N_DS = sum(N_DS), X_SD = sum(X_SD), X_DS = sum(X_DS)), by = map], use.names = TRUE)
  M[, .(E = WIN$E[match(window[1], WIN$window)], N_SD = mean(N_SD), T_S = mean(T_S), exp_SD = mean(X_SD),
        ratio_SD = sum(N_SD) / sum(X_SD), ratio_SD_lo = qn(N_SD / X_SD, 0.025), ratio_SD_hi = qn(N_SD / X_SD, 0.975),
        P_any_SD = mean(N_SD > 0), N_DS = mean(N_DS), T_D = mean(T_D), exp_DS = mean(X_DS), ratio_DS = sum(N_DS) / sum(X_DS)),
    by = window]
}
# ---- (3) origins ----------------------------------------------------------------------------------------------------
# windows (lo, hi] in Ma; per draw: origins C and expected X = q_hat x shallow lineage time (shallow = edges whose child
# clade is not all-D, plus the stem above the drawn origin age; q_hat = N / total shallow time); rand: also length(lo)
# random non-overlapping WIDTH windows per draw, older bounds uniform in [WIDTH, min(541, root age)]
origins_null <- function(tr, x, ndraw, seed, lo = EVT$E - WIDTH, hi = EVT$E, rand = TRUE) {
  set.seed(seed); n <- Ntip(tr); dep <- node.depth.edgelength(tr); H <- max(dep[seq_len(n)]); age <- H - dep; e <- tr$edge
  od <- deep_clades(tr, x); N <- nrow(od); if (!N) return(NULL)
  po <- reorder(tr, "postorder")$edge; nt <- c(rep(1L, n), integer(tr$Nnode)); nd <- c(as.integer(x[tr$tip.label]), integer(tr$Nnode))
  for (i in seq_len(nrow(po))) { nt[po[i, 1]] <- nt[po[i, 1]] + nt[po[i, 2]]; nd[po[i, 1]] <- nd[po[i, 1]] + nd[po[i, 2]] }
  sh <- nd[e[, 2]] < nt[e[, 2]]; ya <- age[e[sh, 2]]; yb <- age[e[sh, 1]]
  br <- sort(unique(c(0, ya, yb))); Ff <- approxfun(br, vapply(br, function(u) sum(pmax(0, pmin(u, yb) - ya)), 0), rule = 2)
  cr <- od$crown_age; st <- od$stem_age; o <- matrix(runif(N * ndraw, cr, st), N)
  qh <- N / (Ff(H) + colSums(st - o)); cs <- function(v) colSums(matrix(v, N))
  cnt <- function(l, h) {                               # l, h: scalars or one value per draw
    L <- matrix(l, N, ndraw, byrow = TRUE); U <- matrix(h, N, ndraw, byrow = TRUE)
    list(C = colSums(o > L & o <= U), X = qh * (Ff(h) - Ff(l) + cs(pmax(0, pmin(U, st) - pmax(o, L)))))
  }
  k <- length(lo); w9 <- lapply(seq_len(k), function(i) cnt(lo[i], hi[i]))
  r <- list(od = od, N = N, H = H, qh = qh, C = do.call(rbind, lapply(w9, `[[`, "C")), X = do.call(rbind, lapply(w9, `[[`, "X")))
  if (rand) {
    r$Hm <- Hm <- min(541, H); if (Hm - (k - 1) * WIDTH <= WIDTH) stop("tree too young for ", k, " random windows")
    Hr <- matrix(apply(matrix(runif(k * ndraw, WIDTH, Hm - (k - 1) * WIDTH), k), 2, sort), k) + (seq_len(k) - 1) * WIDTH
    wr <- lapply(seq_len(k), function(i) cnt(Hr[i, ] - WIDTH, Hr[i, ]))
    r$C2 <- do.call(rbind, lapply(wr, `[[`, "C")); r$X2 <- do.call(rbind, lapply(wr, `[[`, "X"))
  }
  r
}
orig_summary <- function(r, nm = EVT$event, E = EVT$E) {
  R1 <- r$C / r$X; R1[!is.finite(R1)] <- NA; R2 <- r$C2 / r$X2; R2[!is.finite(R2)] <- NA
  one <- function(lab, E, C, X, R, g) { ok <- !is.na(g)
    data.table(window = lab, E = E, C_mean = mean(C), X_mean = mean(X), R_mean = mean(R, na.rm = TRUE), R_med = qn(R, 0.5),
               R_lo = qn(R, 0.025), R_hi = qn(R, 0.975), p_rand = (1 + sum(g[ok])) / (1 + sum(ok))) }
  P1 <- colSums(r$C) / colSums(r$X); P2 <- colSums(r$C2) / colSums(r$X2)
  rbind(rbindlist(lapply(seq_along(nm), function(i) one(nm[i], E[i], r$C[i, ], r$X[i, ], R1[i, ], sweep(R2, 2, R1[i, ], ">=")))),
        one("pooled_9", NA_real_, colSums(r$C), colSums(r$X), P1, P2 >= P1))
}
# deep clades whose crown is older than each epoch boundary: deep lineages that predate it (against complete post-event
# replacement; a minimum count of colonisations before it); supported = P(D) at the crown >= 0.95 in the stochastic
# maps (node numbers of sm$tr = those of tr)
predate_tab <- function(od, sm, B = EPO) {
  pD <- rowMeans(sm$ST[od$node, , drop = FALSE] == 1L)
  rbindlist(lapply(B, function(b) data.table(boundary = b, n_clades = nrow(od), n_older = sum(od$crown_age > b),
    n_supported = sum(od$crown_age > b & pD >= 0.95), n_tips_older = sum(od$n_D[od$crown_age > b]), oldest_crown = max(od$crown_age))))
}
# ---- (4) robustness -------------------------------------------------------------------------------------------------
rd_trees <- function(f) {
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  t <- read.tree(f); if (inherits(t, "phylo")) list(t) else lapply(seq_along(t), function(i) t[[i]])
}
dir_tree <- function(t) tryCatch({
  m <- prep(t); f <- mk2_fit(m$tr, m$x, "ARD"); r <- dir_all(m$tr, m$x, f); od <- deep_clades(m$tr, m$x)
  r[, `:=`(model_h = f$model, q_SD_h = f$q_SD, q_DS_h = f$q_DS, root_age = max(node.depth.edgelength(m$tr)), n_tip = Ntip(m$tr),
           n_crown_older = vapply(E, function(e) sum(od$crown_age > e), 0L))][]
}, error = function(e) { lg("robustness tree skipped:", conditionMessage(e)); NULL })
# ---- main -----------------------------------------------------------------------------------------------------------
options(width = 250)
lg("58 timing:", TREE, "| MA", MA, "| labels", LC, "| width", WIDTH, "| epochs", paste(EPO, collapse = ","), "| dropped",
   paste(OG, collapse = ","))
m0 <- prep(load_tree(TREE, ma = MA)); tr <- m0$tr; x <- m0$x; H <- max(node.depth.edgelength(tr)[seq_len(Ntip(tr))])
lg(Ntip(tr), "tips (D", sum(x), "/ S", sum(x == 0L), "), root age", round(H, 2), "Ma")
f0 <- mk2_fit(tr, x, "ARD"); f0b <- on_bd(c(f0$q_SD, f0$q_DS), mk2_bounds(tr))
lg("homogeneous", f0$model, ": q_SD", signif(f0$q_SD, 4), "q_DS", signif(f0$q_DS, 4), "per Myr, logLik", round(f0$loglik, 3),
   if (f0b) "(a rate on an optimiser bound)" else "")
dm <- dir_all(tr, x, f0, par = TRUE)
lg("directional fits done; min Holm p", signif(min(dm$p_holm, na.rm = TRUE), 3), "(", dm$event[which.min(dm$p_holm)], ")")
bt <- if (NBOOT > 0) boot_dir(tr, x, f0, dm) else NULL
dm[, `:=`(p_boot = if (is.null(bt)) NA_real_ else bt$p_boot, p_fw = if (is.null(bt)) NA_real_ else bt$p_fw)]
fwrite(dm, file.path(OUT, "dir_main.tsv"), sep = "\t")
if (!is.null(bt) && nrow(bt$B)) fwrite(cbind(bt$B, tol = bt$tol, acc_rate = bt$acc), file.path(OUT, "dir_boot.tsv"), sep = "\t")
e4 <- era4_all(tr, x, f0); if (!is.null(e4)) fwrite(e4, file.path(OUT, "era4_main.tsv"), sep = "\t")
lg("era model (secondary):", if (is.null(e4)) "not fitted" else paste("min Holm p", signif(min(e4$p_holm), 3)))
pw <- if (NPOW > 0 && !is.null(bt) && !is.null(bt$mx)) power_run(tr, x, f0, bt) else NULL
ps <- if (!is.null(pw) && nrow(pw)) pw_sum(pw) else NULL
if (!is.null(ps)) { fwrite(pw, file.path(OUT, "power.tsv"), sep = "\t"); fwrite(ps, file.path(OUT, "power_summary.tsv"), sep = "\t") }
lg("power module:", if (is.null(ps)) "not run" else paste(nrow(ps), "cells,", nrow(pw), "labellings"))
sm <- smap(tr, x, f0, NMAP, SEED); ss <- smap_summary(win_tab(sm), f0)
fwrite(ss, file.path(OUT, "smap_windows.tsv"), sep = "\t")
lg("stochastic maps:", NMAP, "| whole-tree S->D / expected", round(ss[window == "all", ratio_SD], 3), "(1 expected)")
og <- origins_null(tr, x, NDRAW, SEED + 2000L); ow <- NULL; pd <- NULL
if (!is.null(og)) {
  od <- copy(og$od); nk <- nrow(EVT)
  ms <- matrix(vapply(seq_len(nk), function(k) pmax(0, pmin(od$stem_age, EVT$E[k]) - pmax(od$crown_age, EVT$E[k] - WIDTH)) /
                        (od$stem_age - od$crown_age), numeric(nrow(od))), nrow(od), dimnames = list(NULL, paste0("w", EVT$E)))
  od <- cbind(od, as.data.table(ms)); ow <- orig_summary(og); pd <- predate_tab(og$od, sm)
  fwrite(od, file.path(OUT, "origins.tsv"), sep = "\t"); fwrite(ow, file.path(OUT, "origin_windows.tsv"), sep = "\t")
  fwrite(pd, file.path(OUT, "predate.tsv"), sep = "\t")
  lg(og$N, "origins (maximal all-D clades); pooled R", round(ow[window == "pooled_9", R_mean], 3), "p_rand", signif(ow[window == "pooled_9", p_rand], 3),
     "| crowns older than", paste(pd$boundary, "Ma:", pd$n_older, collapse = ", "))
} else lg("no deep clades: origin analysis skipped")
RB <- rd_trees(ROB); if (!is.null(RB) && NROB > 0) RB <- RB[seq_len(min(NROB, length(RB)))]
dr <- if (length(RB)) rbindlist(Filter(is.data.table, mcl(seq_along(RB), function(i) {
  r <- dir_tree(RB[[i]]); if (!is.null(r)) r[, tree := i]; r }))) else data.table()
if (nrow(dr)) fwrite(dr, file.path(OUT, "dir_robust.tsv"), sep = "\t")
lg("robustness:", length(RB), "trees read,", if (nrow(dr)) uniqueN(dr$tree) else 0L, "fitted")
pt <- rbindlist(lapply(PTS[nzchar(PTS) & file.exists(PTS)], function(f) {
  tl <- rd_trees(f); rbindlist(lapply(seq_along(tl), function(i) { r <- dir_tree(tl[[i]]); if (!is.null(r)) r[, `:=`(file = basename(f), tree = i)]; r }))
}))
if (nrow(pt)) fwrite(pt, file.path(OUT, "dir_point.tsv"), sep = "\t")
# ---- decision and report ----------------------------------------------------------------------------------------------
dc <- dm[, .(family, event, E, width, LR, z, p1, p_holm, p_boot, p_fw, ratio, q_SD_bg, q_SD_in, at_bound)]
dc[, `:=`(sig_main = p_holm < 0.05 & p_fw < 0.05, n_robust = 0L, frac_robust = NA_real_, loc_mult = NA_real_, loc_power = NA_real_)]
if (nrow(dr)) {
  fr <- dr[, .(n_robust = .N, frac_robust = mean((p_holm < 0.05) %in% TRUE)), by = .(family, E)]
  dc[fr, on = .(family, E), `:=`(n_robust = i.n_robust, frac_robust = i.frac_robust)]
}
if (!is.null(ps)) for (k in which(is.finite(dc$ratio) & dc$ratio > 0)) {       # power cell nearest the observed ratio
  mu <- if (dc$family[k] == "window") PMULT else PEMULT; mm <- mu[which.min(abs(log(mu) - log(dc$ratio[k])))]
  v <- ps[family == dc$family[k] & abs(E - dc$E[k]) < 1e-6 & mult == mm, det_loc]
  if (length(v) == 1L) set(dc, k, c("loc_mult", "loc_power"), list(mm, v))
}
dc[, robust := frac_robust >= 0.8]
# (1a) attribution: per family the effect interval with the largest z; the other effect intervals conditional on it
dc[, `:=`(attr = FALSE, frac_attr = NA_real_, cond_given = NA_character_, cond_ratio = NA_real_, cond_z = NA_real_,
          cond_p_holm = NA_real_)]
pp0 <- mk2_prep(tr, x); bd0 <- mk2_bounds(tr); cd <- list()
for (fm in unique(dc$family)) {
  ix <- which(dc$family == fm & dc$sig_main %in% TRUE); if (!length(ix)) next
  ka <- ix[which.max(dc$z[ix])]; set(dc, ka, "attr", TRUE)
  if (nrow(dr)) { mz <- dr[family == fm, .SD[which.max(z)], by = tree]; set(dc, ka, "frac_attr", mean(abs(mz$E - dc$E[ka]) < 1e-6)) }
  ko <- setdiff(ix, ka); if (!length(ko)) next
  cf <- rbindlist(lapply(ko, function(k) cond_fit(pp0, bd0, dc$E[ka], dc$width[ka], dc$E[k], dc$width[k], f0, dm[ka])))
  cf[, `:=`(family = fm, attr_event = dc$event[ka], event = dc$event[ko], E = dc$E[ko], p_holm = p.adjust(p1, "holm"))]
  set(dc, ko, c("cond_given", "cond_ratio", "cond_z", "cond_p_holm"), list(dc$event[ka], cf$ratio, cf$z, cf$p_holm))
  cd[[fm]] <- cf
}
cd <- rbindlist(cd)
if (nrow(cd)) { setcolorder(cd, c("family", "attr_event", "event", "E")); fwrite(cd, file.path(OUT, "attribution.tsv"), sep = "\t") }
lg("attribution:", if (nrow(cd)) paste(nrow(cd), "conditional fits") else "no family with two or more effect intervals")
call_of <- function(k) {
  if (is.na(dc$sig_main[k])) return("untested (no bootstrap)")
  if (!dc$sig_main[k]) return("no effect")
  fr <- dc$frac_robust[k]; lp <- dc$loc_power[k]; nb <- paste0(dc$family[k], "s")
  s <- if (is.na(fr)) "effect, dating robustness not assessed" else if (fr >= 0.8) "effect, robust to dating" else "effect, not robust to dating"
  lq <- if (is.na(lp)) ", localisation not assessed" else if (lp < 0.8) paste(", not separable from adjacent", nb) else ""
  if (dc$attr[k]) return(paste0(s, ", attributed (largest z of the ", nb, if (!is.na(dc$frac_attr[k]))
    sprintf("; in %.0f %% of the ROBUST trees", 100 * dc$frac_attr[k]), ")", lq))
  if (is.na(dc$cond_p_holm[k])) return(paste0(s, "; conditional test given ", dc$cond_given[k], " failed"))
  if (dc$cond_p_holm[k] < 0.05) return(sprintf("%s, separate effect given %s (conditional Holm p %.2g)%s", s, dc$cond_given[k],
                                               dc$cond_p_holm[k], lq))
  sprintf("%s; explained by %s (conditional Holm p %.2g)", s, dc$cond_given[k], dc$cond_p_holm[k])
}
dc[, call := vapply(seq_len(.N), call_of, "")]
fwrite(dc, file.path(OUT, "decision.tsv"), sep = "\t")
pr <- function(d) capture.output(print(d, digits = 3, row.names = FALSE))
rp <- c(paste("58_origin_timing", format(Sys.time(), "%F %T")),
        sprintf("tree %s (MA %s), labels %s, dropped %s", TREE, MA, LC, paste(OG, collapse = ",")),
        sprintf("%d tips (D %d, S %d), root age %.2f Ma; homogeneous %s: q_SD %.4g, q_DS %.4g per Myr, logLik %.3f%s",
                Ntip(tr), sum(x), sum(x == 0L), H, f0$model, f0$q_SD, f0$q_DS, f0$loglik, if (f0b) " (a rate on an optimiser bound)" else ""),
        "", sprintf(paste("(1) Directional tests (own S->D rate inside the interval, 1 df): windows [E, E - %g] Ma, epochs [B, 0] Ma; z = sign(log ratio)",
                          "x sqrt(LR), Holm within family; p_boot per interval, p_fw = max-z within family; effect = Holm p < 0.05 and p_fw < 0.05;",
                          "robust = Holm p < 0.05 in >= 80 %% of the ROBUST trees; an attributed or separate effect (1a) is localised only if",
                          "det_loc >= 0.8 at loc_mult:"), WIDTH),
        pr(dc[, .(family, event, LR, z, p_holm, p_boot, p_fw, ratio, at_bound, frac_robust, loc_mult, loc_power, call)]),
        if (!is.null(bt)) sprintf("bootstrap: %d labellings fitted (deep-count tolerance %.2f, acceptance %.3g)", bt$n, bt$tol, bt$acc),
        "", paste("(1a) Attribution: per family the effect interval with the largest z; every other effect interval refitted conditional on it",
                  "(own S->D rate on its part outside the attributed interval vs the background, or, if nested inside it, on itself vs the rest",
                  "of the attributed interval; 1 df, one-sided asymptotic p, Holm within family); separate effect = conditional Holm p < 0.05:"),
        if (nrow(cd)) pr(cd[, .(family, attr_event, event, nested, free_bl, ratio, LR, z, p1, p_holm, at_bound)]) else
          "no family with two or more effect intervals",
        "", "(1b) Secondary (ED): 4-parameter era model (own q_SD and q_DS in the window), asymptotic LRT df = 4 - npar, Holm over the windows:",
        if (is.null(e4)) "not fitted" else pr(e4[, .(event, LR, df, p, p_holm, ratio_SD, ratio_DS, at_bound)]),
        "", sprintf(paste("(1c) Power and localisation (%d labellings per cell, background S->D rescaled to the observed deep count): power = rule (1)",
                          "for the planted interval, loc = its z is the largest of the family, det_loc = both:"), NPOW),
        if (is.null(ps)) "power module not run" else pr(ps),
        "", sprintf("(2) %d stochastic maps under the homogeneous fit; whole-tree S->D / expected = %.3f, D->S = %.3f (1 expected at the ML rates):",
                    NMAP, ss[window == "all", ratio_SD], ss[window == "all", ratio_DS]),
        pr(ss),
        "", if (is.null(ow)) "(3) no deep clades" else c(sprintf(
          "(3) %d origins (maximal all-D clades), %d origin-age draws; R = origins / expected (constant colonisation per shallow lineage-time), p_rand = %d random WIDTH windows per draw in [%g, %g] Ma:",
          og$N, NDRAW, nrow(EVT), WIDTH, og$Hm), pr(ow),
          "predate: deep clades whose crown is older than each epoch boundary (supported = P(D) at the crown >= 0.95 in the maps):", pr(pd)),
        "", sprintf("(4) robustness: %d ROBUST trees read (%s), %d fitted; %d POINT trees fitted", length(RB), ROB,
                    if (nrow(dr)) uniqueN(dr$tree) else 0L, if (nrow(pt)) nrow(unique(pt[, .(file, tree)])) else 0L),
        if (nrow(dr)) pr(dr[, .(n = .N, frac_holm = mean((p_holm < 0.05) %in% TRUE), ratio_med = median(ratio, na.rm = TRUE),
                                n_crown_older_med = median(n_crown_older), n_crown_older_min = min(n_crown_older),
                                n_crown_older_max = max(n_crown_older)), by = .(family, event)]),
        if (nrow(pt)) pr(pt[, .(file, tree, family, event, LR, z, p_holm, ratio, root_age, n_crown_older)]),
        "", paste("Caveat: from extant tips a raised recent S->D rate cannot be separated from extinction of old deep lineages",
                  "(state-dependent extinction); both read as recent recruitment of the extant deep fauna."))
writeLines(rp, file.path(OUT, "report.txt")); lg("done:", OUT)
