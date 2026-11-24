#!/usr/bin/env Rscript
# A1 (PREREG_addon.md): annotation-route, genome-GC and proteome-quality confounds of the physicochemistry PGLS
# (32 / 44), with simulation-calibrated p values.
# Usage: Rscript 53_route_confound.R [TREE] [OUTDIR]     (default work/final_v3/date/lsd.timetree.nwk, work/addon/A1)
# env LABEL_COL (lab), DATED (1), NCORES (16), NSIM_SP (1000 Mk labellings, species null), NSIM (the first 100 of them
#   reused for the gene null), NBOOT (200), SEED (20261025), PART (all | species | gene | null | decide; all = in order),
#   REF44 (44 controls_species.tsv) / REF32 (32 pgls_gene.tsv): reproduction checks of S_route3 / G_ref (check_ref.txt)
# Species (traits P6 + GC set + busco_* analogues of the whole_* P6 traits; SD units of all labelled species, as 44):
#   S_route3 deep + route (= 44 ref_all), S_noroute, S_route2 (annot2), S_noTSA, S_mp, S_an (one route only),
#   S_int (no TSA, deep x other, other = non-annot), S_qc (+ protC + protD), S_hq (protC >= 90), S_ggc (no TSA, + gGC);
#   subsets need >= 6 D and >= 6 S (44's clade rule). GC decomposition: G1 gGC ~ deep + route; G3 (busco_GC3 - gGC),
#   (busco_cds_GC - gGC) ~ deep + route; G4 P6 ~ deep + route + busco_GC3. T5: Wald CI (t, n - k df), phylolm
#   parametric bootstrap (P6, S_route3), p_sim of S_route3 over Mk ARD labellings conditioned on the deep count.
# Gene level (32's core traits + f_T, f_F; locus >= 20 species, >= 5 D, >= 5 S; y scaled within locus, as 32):
#   G_ref (= 32), G_route, G_mp, G_an; calibrated sign test of G_route on the first NSIM labellings.
# Outputs (plain tsv; fread of .gz would need R.utils): species_specs.tsv, gc_decomp.tsv, mk_fit.tsv, sim_labels.tsv,
#   sim_t_species.tsv, species_null.tsv, gene_specs.tsv, gene_summary.tsv, gene_null/sim_NNN.tsv, gene_null_summary.tsv,
#   decision.tsv, report.txt,
#   check_ref.txt
fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(Sys.getenv("ADDON_LIB", file.path(dirname(normalizePath(sub("^--file=", "", fa[1]))), "addon_lib.R")))
suppressPackageStartupMessages({ library(phylolm); library(parallel) })
a <- commandArgs(TRUE)
TREE <- if (length(a) >= 1) a[1] else file.path(W, "date/lsd.timetree.nwk")
OUT <- if (length(a) >= 2) a[2] else file.path(R, "work/addon/A1")
dir.create(file.path(OUT, "gene_null"), showWarnings = FALSE, recursive = TRUE)
if (!nzchar(Sys.getenv("DATED"))) Sys.setenv(DATED = "1")
NC <- as.integer(Sys.getenv("NCORES", "16")); NSIM_SP <- as.integer(Sys.getenv("NSIM_SP", "1000"))
NSIM <- as.integer(Sys.getenv("NSIM", "100")); NBOOT <- as.integer(Sys.getenv("NBOOT", "200"))
SEED <- as.integer(Sys.getenv("SEED", "20261025")); PART <- Sys.getenv("PART", "all")
lg <- function(...) { cat(format(Sys.time(), "%F %T"), ..., "\n"); flush(stdout()) }
L <- load_species(read.tree(TREE)); sp <- L$sp; tr <- L$tr; rm(L)
lg("labels", Sys.getenv("LABEL_COL", "lab"), ": D", sum(sp$deep), "S", sum(1 - sp$deep), "; tips", Ntip(tr), "; route",
   paste(names(table(sp$route)), table(sp$route), collapse = " "), "; gGC", sum(is.finite(sp$gGC)), "; protC",
   sum(is.finite(sp$protC)))
AN <- intersect(sub("^whole_", "busco_", grep("^whole_", P6, value = TRUE)), names(sp))   # busco_* analogues
TR <- intersect(unique(c(P6, GCSET, AN)), names(sp)); TRN <- intersect(c(P6, GCSET), TR)
if (length(setdiff(c(P6, GCSET), TR))) lg("WARNING missing traits:", setdiff(c(P6, GCSET), TR))
CORE <- c("acid_base", "net_charge7_per_res", "pI", "R_RK", "f_C", "f_P", "f_G", "gravy", "GC3", "GC4", "IVYWREL")  # 32
GT <- c(CORE, "f_T", "f_F"); GN <- c("pI", "R_RK", "acid_base", "net_charge7_per_res", "f_T", "f_F", "GC3")
cov_ggc <- sp[route != "tsa", mean(is.finite(gGC))]
yof <- function(tv) { d <- sp[is.finite(get(tv))]; d[, y := zs(get(tv))]; d }   # SD units of all labelled species
fw <- function(...) suppressWarnings(fit1(...))
wci <- function(r) { k <- nrow(r); r[, `:=`(ci_lo = est - qt(0.975, n - k) * se, ci_hi = est + qt(0.975, n - k) * se)] }
ok_n <- function(d, m = 6) sum(d$deep) >= m && sum(1 - d$deep) >= m
SPECS <- list(   # subset of the trait table, right-hand side, route variable (added when it varies)
  S_route3 = list(function(d) d, "deep", "route"), S_noroute = list(function(d) d, "deep", NULL),
  S_route2 = list(function(d) d, "deep", "annot2"), S_noTSA = list(function(d) d[route != "tsa"], "deep", "route"),
  S_mp = list(function(d) d[route == "miniprot"], "deep", NULL), S_an = list(function(d) d[route == "annot"], "deep", NULL),
  S_qc = list(function(d) d[is.finite(protC) & is.finite(protD)][, `:=`(zC = zs(protC), zD = zs(protD))],
              "deep + zC + zD", "route"),
  S_hq = list(function(d) d[protC >= 90], "deep", "route"),
  S_ggc = list(function(d) d[route != "tsa" & is.finite(gGC)][, zG := zs(gGC)], "deep + zG", "route"))
fit_int <- function(d) {   # S_int: no TSA; deep x other (other = non-annot route); deep effect within each stratum
  df <- as.data.frame(d[route != "tsa"]); rownames(df) <- df$code; df$other <- as.integer(df$annot2 == "other")
  m <- tryCatch(suppressWarnings(phylolm(y ~ deep * other, data = df, phy = keep.tip(tr, df$code), model = "lambda")),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- summary(m)$coefficients; V <- if (is.null(m$vcov)) vcov(m) else m$vcov; n <- nrow(df); k <- nrow(cf)
  b <- c("deep", "deep:other"); e <- c(cf[b, "Estimate"], sum(cf[b, "Estimate"])); s <- c(cf[b, "StdErr"], sqrt(sum(V[b, b])))
  data.table(term = c("deep", "deep:other", "deep|other"), est = e, se = s, t = e / s, p = 2 * pt(-abs(e / s), n - k),
             lambda = m$optpar, n = n, n_deep = sum(df$deep), aic = m$aic, rhs = "deep * other",
             n_deep_stratum = c(sum(df$deep[df$other == 0]), NA, sum(df$deep[df$other == 1])),
             ci_lo = e - qt(0.975, n - k) * s, ci_hi = e + qt(0.975, n - k) * s)
}
run_species <- function() {
  lg("species:", length(TR), "traits x", length(SPECS) + 1, "specs; bootstrap", NBOOT, "; gGC cover (no TSA)", round(cov_ggc, 3))
  res <- mclapply(TR, function(tv) {
    d0 <- yof(tv); o <- list()
    for (nm in names(SPECS)) {
      s <- SPECS[[nm]]; d <- s[[1]](copy(d0)); if (!ok_n(d)) next
      r <- fw(d, s[[2]], tr, rv = s[[3]], boot = if (nm == "S_route3" && tv %in% P6) NBOOT else 0L)
      if (!is.null(r)) o[[length(o) + 1]] <- wci(r)[term != "(Intercept)"][, `:=`(trait = tv, spec = nm)]
    }
    r <- fit_int(d0); if (!is.null(r)) o[[length(o) + 1]] <- r[, `:=`(trait = tv, spec = "S_int")]
    rbindlist(o, fill = TRUE)
  }, mc.cores = min(NC, length(TR)))
  bad <- !vapply(res, is.data.table, TRUE); if (any(bad)) lg("WARNING failed traits:", TR[bad])
  ss <- rbindlist(res[!bad], fill = TRUE); ss[, q := p.adjust(p, "BH"), by = .(spec, term)]
  fwrite(ss, file.path(OUT, "species_specs.tsv"), sep = "\t")
  go <- list(); dg <- sp[is.finite(gGC)]
  addg <- function(r, nm, tv) if (!is.null(r)) go[[length(go) + 1]] <<- wci(r)[term != "(Intercept)"][, `:=`(trait = tv, spec = nm)]
  if (ok_n(dg)) {
    addg(fw(copy(dg)[, y := zs(gGC)], "deep", tr), "G1", "genome_GC")
    for (tv in intersect(c("busco_GC3", "busco_cds_GC"), names(sp)))
      addg(fw(dg[is.finite(get(tv))][, y := zs(get(tv) - gGC)], "deep", tr), "G3", paste0(tv, "-gGC"))
  }
  for (tv in intersect(P6, names(sp)))
    addg(fw(yof(tv)[is.finite(busco_GC3)][, zGC3 := zs(busco_GC3)], "deep + zGC3", tr), "G4", tv)
  fwrite(rbindlist(go, fill = TRUE), file.path(OUT, "gc_decomp.tsv"), sep = "\t")
  lg("species specs", nrow(ss), "rows; GC decomposition", length(go), "fits")
  ss
}
run_species_null <- function(ss) {   # Mk ARD on the observed labels; labellings with the deep count within 10 % (20 %)
  x <- setNames(sp$deep, sp$code); mk <- mk2_fit(tr, x, "ARD")
  sc <- mk2_sim_cond(tr, mk$q_SD, mk$q_DS, NSIM_SP, sum(x), seed = SEED); S <- sc$sims
  colnames(S) <- sprintf("s%04d", seq_len(ncol(S)))
  mf <- data.table(model = mk$model, q_SD = mk$q_SD, q_DS = mk$q_DS, loglik = mk$loglik, conv = mk$conv,
                   parsimony = mk$k_pars, n_D = sum(x), n_S = sum(1 - x), n_sim = ncol(S), tol = sc$tol,
                   acc_rate = sc$acc_rate, n_draw = sc$n_draw, sim_nD_min = min(colSums(S)), sim_nD_max = max(colSums(S)))
  fwrite(mf, file.path(OUT, "mk_fit.tsv"), sep = "\t"); print(mf)
  fwrite(data.table(code = rownames(S), S), file.path(OUT, "sim_labels.tsv"), sep = "\t")
  D0 <- lapply(setNames(TRN, TRN), yof)
  tn <- mclapply(seq_len(ncol(S)), function(j) {
    xs <- S[, j]
    vapply(TRN, function(tv) { d <- copy(D0[[tv]])[, deep := as.integer(xs[code])]
      r <- fw(d, "deep", tr); v <- if (is.null(r)) numeric(0) else r[term == "deep", t]
      if (length(v) == 1) v else NA_real_ }, 0)
  }, mc.cores = NC)
  ok <- vapply(tn, is.numeric, TRUE); if (any(!ok)) lg("WARNING null workers failed:", sum(!ok))
  TN <- do.call(rbind, tn[ok]); fwrite(data.table(sim = colnames(S)[ok], TN), file.path(OUT, "sim_t_species.tsv"), sep = "\t")
  tob <- ss[spec == "S_route3" & term == "deep"]
  nl <- rbindlist(lapply(TRN, function(tv) {
    to <- tob[trait == tv, t]; ts <- TN[, tv]; ts <- ts[is.finite(ts)]; if (length(to) != 1 || !length(ts)) return(NULL)
    data.table(trait = tv, t_obs = to, p_wald = tob[trait == tv, p], n_sim = length(ts),
               p_sim = (1 + sum(abs(ts) >= abs(to))) / (1 + length(ts)), abs_t_sim_q95 = quantile(abs(ts), 0.95),
               mean_t_sim = mean(ts), sd_t_sim = sd(ts))
  }))
  fwrite(nl, file.path(OUT, "species_null.tsv"), sep = "\t"); print(nl)
}
load_genes <- function() {   # as 32: gene_metrics of the labelled species, deep and route attached, one table per locus
  g <- fread(cmd = paste("zcat", shQuote(file.path(R, "work/physchem/gene_metrics.tsv.gz"))))
  split(merge(g[code %in% sp$code], sp[, .(code, deep, route)], by = "code"), by = "busco")
}
gfit <- function(gb, traits, specs, lab = NULL) {   # one locus; lab = named 0/1 simulated labelling (NULL = observed)
  if (!is.null(lab)) gb <- copy(gb)[, deep := as.integer(lab[code])]
  if (nrow(gb) < 20 || sum(gb$deep) < 5 || sum(1 - gb$deep) < 5) return(NULL)
  o <- list()
  for (tv in intersect(traits, names(gb))) {
    d <- gb[is.finite(get(tv))]; if (nrow(d) < 20) next
    d[, y := zs(get(tv))]; if (!is.finite(sd(d$y))) next
    for (nm in specs) {
      dd <- switch(nm, G_mp = d[route == "miniprot"], G_an = d[route == "annot"], d)
      if (nm %in% c("G_mp", "G_an") && (nrow(dd) < 20 || sum(dd$deep) < 5 || sum(1 - dd$deep) < 5)) next
      r <- fw(dd, "deep", tr, rv = if (nm == "G_route") "route" else NULL)
      if (!is.null(r)) o[[length(o) + 1]] <- r[term == "deep", .(busco = gb$busco[1], trait = tv, spec = nm, est, se, t, p,
                                                                 lambda, n, n_deep)]
    }
  }
  rbindlist(o)
}
keep_dt <- function(l) rbindlist(Filter(function(z) is.data.table(z) && nrow(z), l))
run_gene <- function() {
  G <- load_genes(); lg("gene level:", length(G), "loci; traits", paste(intersect(GT, names(G[[1]])), collapse = " "))
  gl <- mclapply(G, gfit, traits = GT, specs = c("G_ref", "G_route", "G_mp", "G_an"), mc.cores = NC)
  ne <- sum(vapply(gl, inherits, TRUE, "try-error")); gl <- keep_dt(gl); gl[, q := p.adjust(p, "BH"), by = .(trait, spec)]
  fwrite(gl, file.path(OUT, "gene_specs.tsv"), sep = "\t")
  gs <- gl[, .(n_loci = .N, frac_pos = mean(est > 0), sign_p = binom.test(sum(est > 0), .N)$p.value,
               n_q05 = sum(q < 0.05), median_est = median(est)), by = .(trait, spec)]
  fwrite(gs, file.path(OUT, "gene_summary.tsv"), sep = "\t"); lg("gene fits", nrow(gl), "; worker errors", ne); print(gs)
}
run_null <- function() {   # calibrated sign test: G_route on the first NSIM simulated labellings (resumable per labelling)
  sl <- fread(file.path(OUT, "sim_labels.tsv")); S <- as.matrix(sl[, -1]); rownames(S) <- sl$code
  G <- load_genes(); ns <- min(NSIM, ncol(S)); lg("gene null:", ns, "labellings x", length(G), "loci x", length(GN), "traits")
  for (j in seq_len(ns)) {
    f <- file.path(OUT, "gene_null", sprintf("sim_%03d.tsv", j)); if (file.exists(f)) next
    gl <- mclapply(G, gfit, traits = GN, specs = "G_route", lab = S[, j], mc.cores = NC)
    ne <- sum(vapply(gl, inherits, TRUE, "try-error")); if (ne) lg("WARNING sim", j, "worker errors", ne)
    tf <- sub("\\.tsv$", ".part", f); fwrite(keep_dt(gl), tf, sep = "\t"); file.rename(tf, f)
    if (j %% 10 == 0) lg("gene null", j, "/", ns)
  }
  ob <- fread(file.path(OUT, "gene_specs.tsv"))[spec == "G_route", .(f_obs = mean(est > 0), n_loci_obs = .N), by = trait]
  fs <- list.files(file.path(OUT, "gene_null"), "^sim_[0-9]+\\.tsv$", full.names = TRUE)
  nu <- rbindlist(lapply(fs, function(ff) { x <- fread(ff); if (!nrow(x)) return(NULL); sid <- basename(ff)
    x[, .(f = mean(est > 0), n = .N), by = trait][, sim := sid] }))
  gn <- merge(nu, ob, by = "trait")[, .(f_obs = f_obs[1], n_loci_obs = n_loci_obs[1], n_sim = .N, f_sim_mean = mean(f),
    f_sim_lo = quantile(f, 0.025), f_sim_hi = quantile(f, 0.975), n_loci_sim = median(n),
    p_cal = (1 + sum(abs(f - 0.5) >= abs(f_obs[1] - 0.5))) / (1 + .N)), by = trait]
  fwrite(gn, file.path(OUT, "gene_null_summary.tsv"), sep = "\t"); print(gn)
}
run_decide <- function() {   # decision rules (a)-(g) per P6 trait; GC3 withdrawal rule; reproduction checks
  ss <- fread(file.path(OUT, "species_specs.tsv")); nl <- fread(file.path(OUT, "species_null.tsv"))
  gnf <- file.path(OUT, "gene_null_summary.tsv"); gn <- if (file.exists(gnf)) fread(gnf) else data.table(trait = "", p_cal = NA_real_)
  v <- function(tv, s, tm = "deep", col = "est") { z <- ss[trait == tv & spec == s & term == tm][[col]]; if (length(z)) z[1] else NA_real_ }
  ref <- ss[spec == "S_route3" & term == "deep" & trait %in% P6, .(trait, holm = p.adjust(p, "holm"))]
  dec <- rbindlist(lapply(ref$trait, function(tv) {
    E <- function(s, tm = "deep") v(tv, s, tm); Pv <- function(s, tm = "deep") v(tv, s, tm, "p")
    s0 <- sign(E("S_route3")); sg <- function(z) isTRUE(sign(z) == s0); an <- sub("^whole_", "busco_", tv)
    has_an <- tv != an && an %in% ss$trait
    data.table(trait = tv, est = E("S_route3"), p = Pv("S_route3"), holm = ref[trait == tv, holm],
      mp_est = E("S_mp"), mp_p = Pv("S_mp"), an_est = E("S_an"), an_p = Pv("S_an"), int_p = Pv("S_int", "deep:other"),
      within_annot = E("S_int", "deep"), within_other = E("S_int", "deep|other"), ggc_est = E("S_ggc"), ggc_p = Pv("S_ggc"),
      qc_est = E("S_qc"), qc_p = Pv("S_qc"), p_sim = nl[trait == tv, p_sim][1], analogue = if (has_an) an else NA_character_,
      analogue_est = if (has_an) v(an, "S_route3") else NA_real_, analogue_p = if (has_an) v(an, "S_route3", col = "p") else NA_real_,
      crit_a = ref[trait == tv, holm] < 0.05,
      crit_b = isTRUE(Pv("S_mp") < 0.05) && sg(E("S_mp")) && sg(E("S_an")),
      crit_c = isTRUE(Pv("S_int", "deep:other") >= 0.05) || (sg(E("S_int", "deep")) && sg(E("S_int", "deep|other"))),
      crit_d = if (cov_ggc >= 0.9) isTRUE(Pv("S_ggc") < 0.05) && sg(E("S_ggc")) else NA,
      crit_e = isTRUE(Pv("S_qc") < 0.05) && sg(E("S_qc")), crit_f = isTRUE(nl[trait == tv, p_sim][1] < 0.05),
      crit_g = if (has_an) sg(v(an, "S_route3")) else NA)
  }))
  cr <- grep("^crit_", names(dec), value = TRUE)
  dec[, robust := apply(as.matrix(.SD), 1, all, na.rm = TRUE), .SDcols = cr]
  dec[, failed := apply(as.matrix(.SD), 1, function(z) paste(sub("crit_", "", cr[!is.na(z) & !z]), collapse = ",")), .SDcols = cr]
  dec[, call := fifelse(robust, "route-robust", "route- or quality-sensitive")]
  fwrite(dec, file.path(OUT, "decision.tsv"), sep = "\t")
  gc3_cal <- gn[trait == "GC3", p_cal][1]; gc3_mp <- v("busco_GC3", "S_mp", col = "p")
  gc3_keep <- isTRUE(gc3_cal < 0.05) && isTRUE(gc3_mp < 0.05)
  rep <- c(sprintf("A1 decision (%s; gGC cover of non-TSA species %.3f, criterion d %s)", format(Sys.time(), "%F %T"),
                   cov_ggc, if (cov_ggc >= 0.9) "evaluated" else "not evaluated (< 0.90)"),
           "criteria (b), (c), (d), (e) and (g) require the sign of the S_route3 deep effect (stricter reading of PREREG)",
           sprintf("%-28s est %7.3f holm %.2g p_sim %.3g -> %s%s", dec$trait, dec$est, dec$holm, dec$p_sim, dec$call,
                   ifelse(nzchar(dec$failed), paste0(" (fails ", dec$failed, ")"), "")),
           sprintf("gene-level GC3: calibrated G_route sign p %s, S_mp busco_GC3 p %s -> GC3 statement %s",
                   signif(gc3_cal, 3), signif(gc3_mp, 3), if (gc3_keep) "kept" else "WITHDRAWN"))
  writeLines(rep, file.path(OUT, "report.txt")); cat(rep, sep = "\n")
  ck <- character(0); tc <- function(expr) tryCatch(expr, error = function(e) paste("check failed:", conditionMessage(e)))
  f44 <- Sys.getenv("REF44"); if (nzchar(f44) && file.exists(f44)) ck <- c(ck, tc({
    m <- merge(ss[spec == "S_route3" & term == "deep", .(trait, est, p)],
               fread(f44)[scope == "Animalia" & model == "ref_all" & term == "deep", .(trait, e44 = est, p44 = p)], by = "trait")
    sprintf("S_route3 vs 44 ref_all (%s): %d traits, max |d est| %.3g, max |d log10 p| %.3g", f44, nrow(m),
            max(abs(m$est - m$e44)), max(abs(log10(m$p) - log10(m$p44))))
  }))
  f32 <- Sys.getenv("REF32"); gf <- file.path(OUT, "gene_specs.tsv")
  if (nzchar(f32) && file.exists(f32) && file.exists(gf)) ck <- c(ck, tc({
    m <- merge(fread(gf)[spec == "G_ref", .(busco, trait, est, p)], fread(f32)[, .(busco, trait, e32 = est, p32 = p)],
               by = c("busco", "trait"))
    sprintf("G_ref vs 32 pgls_gene (%s): %d locus x trait, max |d est| %.3g, share |d est| < 1e-6 %.4f", f32,
            nrow(m), max(abs(m$est - m$e32)), mean(abs(m$est - m$e32) < 1e-6))
  }))
  if (length(ck)) { writeLines(ck, file.path(OUT, "check_ref.txt")); cat(ck, sep = "\n") }
}
if (PART %in% c("all", "species")) run_species_null(run_species())
if (PART %in% c("all", "gene")) run_gene()
if (PART %in% c("all", "null")) run_null()
if (PART %in% c("all", "decide")) run_decide()
lg("53 done PART", PART)
