#!/usr/bin/env Rscript
# R1 (PREREG_nature.md): does the RELAX relaxation excess survive lineage-wide confounders?
# Inputs: work/addon/A4/sel_features.tsv (54a), work/addon/R1/row_cov.tsv + tip_cov.tsv (62a), meta/labels_v2.tsv,
#   work/final_v3/date/lsd.timetree.nwk, work/addon/A2/{origins,contrasts}.tsv. Output dir work/addon/R1.
# Models as 54b (binomial, HOG-clustered CR1): M1, M3 = M1 + d_omega, M5 = M1 + d_omega + d_gc + d_buscoC + d_tbl,
#   M6 = M5 + d_temp (exploratory); HQ subset (fg_nonannot = bg_nonannot = 0 in both rows, min BUSCO C >= 80) M0/M1;
#   calibration (null-row model -> expected real rate); dose (real rows, fg_depth); per-clade PGLS of log omega_tip ~ D;
#   threonine origin contrasts vs dGC3/dGC. Scenario A/B/C rule on relax_relaxed is prestated in the PREREG.
suppressPackageStartupMessages({library(data.table); library(ape); library(nlme)})
R <- path.expand("~/deepsea_conv"); O <- file.path(R, "work/addon/R1")
REP <- character(0); lg <- function(...) { s <- paste0(...); cat(s, "\n"); REP <<- c(REP, s) }
d <- fread(file.path(R, "work/addon/A4/sel_features.tsv"), na.strings = c("NA", "nan", "NaN", ""))
rc <- fread(file.path(O, "row_cov.tsv"), na.strings = c("NA", "nan", "NaN", ""))
d <- merge(d[set %in% c("real", "null1") & in_nullhy == 1], rc, by = c("clade", "hog", "set"))
d[, is_real := as.integer(set == "real")]; d[, hid := paste(clade, hog, sep = ":")]
lt <- function(x) fifelse(is.na(x), NA_integer_, as.integer(x < 0.05))
d[, y_relax := lt(relax_p)]
d[, y_relax_relaxed := fifelse(is.na(y_relax) | is.na(relax_K), NA_integer_, as.integer(y_relax == 1L & relax_K < 1))]
cl_vc <- function(m, g) {
  X <- model.matrix(m)[, !is.na(coef(m)), drop = FALSE]
  S <- rowsum(X * (m$y - fitted(m)), g); B <- summary(m)$cov.unscaled
  G <- nrow(S); N <- nrow(X); K <- ncol(X); B %*% crossprod(S) %*% B * (G / (G - 1)) * ((N - 1) / (N - K))
}
fitm <- function(x, rhs, terms) {
  v <- setdiff(all.vars(as.formula(paste("~", rhs))), "clade")
  x <- x[complete.cases(x[, c("y", v), with = FALSE])]
  x <- x[clade %in% x[, .(e = sum(y)), by = clade][e > 0]$clade]
  if (nrow(x) < 20 || sum(x$y) == 0) return(NULL)
  m <- suppressWarnings(glm(as.formula(paste("y ~", rhs)), binomial, data = x)); cf <- coef(summary(m)); V <- cl_vc(m, x$hid)
  rbindlist(lapply(intersect(terms, rownames(cf)), function(t) data.table(term = t, est = cf[t, 1], se_cl = sqrt(V[t, t]),
    p_cl = 2 * pnorm(-abs(cf[t, 1] / sqrt(V[t, t]))), n_rows = nrow(x), n_hog = uniqueN(x$hid), n_events = sum(x$y),
    rate_real = mean(x$y[x$is_real == 1]), rate_null = mean(x$y[x$is_real == 0]))))
}
paired <- function(x) x[hid %in% x[, .(k = any(is_real == 1) & any(is_real == 0)), by = hid][k == TRUE]$hid]
F1 <- "fg_nonannot + log(n_tips) + gapdiff + clade"
SPEC <- list(M0 = "is_real + clade", M1 = paste("is_real +", F1), M3 = paste("is_real + d_omega +", F1),
             M5 = paste("is_real + d_omega + d_gc + d_buscoC + d_tbl +", F1),
             M5b = paste("is_real + d_omega + d_gc3 + d_buscoC + d_tbl +", F1),
             M6 = paste("is_real + d_omega + d_gc + d_buscoC + d_tbl + d_temp +", F1))
TERMS <- c("is_real", "d_omega", "d_gc", "d_gc3", "d_buscoC", "d_tbl", "d_temp", "fg_depth")
MOD <- list(); CAL <- list(); DOSE <- list()
for (ind in c("relax_relaxed", "relax")) {
  x0 <- d[!is.na(get(paste0("y_", ind)))][, y := get(paste0("y_", ind))]
  # common sample: rows with every M5/M6 covariate, paired
  cc <- paired(x0[complete.cases(x0[, .(d_omega, d_gc, d_gc3, d_buscoC, d_tbl, d_temp, fg_nonannot, gapdiff)])])
  for (fr in c("all", "common")) {
    x <- if (fr == "all") paired(x0) else cc
    for (s in names(SPEC)) { r <- fitm(x, SPEC[[s]], TERMS); if (!is.null(r)) MOD[[length(MOD) + 1]] <- r[, `:=`(ind = ind, frame = fr, model = s)] }
  }
  hq <- cc[, .(ok = all(fg_nonannot == 0 & bg_nonannot == 0) & all(min_buscoC >= 80)), by = hid][ok == TRUE]$hid
  xh <- cc[hid %in% hq]; lg(ind, ": common frame ", uniqueN(cc$hid), " HOGs; HQ subset ", length(hq), " HOGs")
  for (s in c("M0", "M1")) { r <- fitm(xh, SPEC[[s]], TERMS); if (!is.null(r)) MOD[[length(MOD) + 1]] <- r[, `:=`(ind = ind, frame = "HQ", model = s)] }
  # calibration: null-row model, expected real rate
  xn <- cc[is_real == 0]; xr <- cc[is_real == 1]
  mn <- suppressWarnings(glm(as.formula(paste("y ~ d_omega + d_gc + d_buscoC + d_tbl +", F1)), binomial, data = xn))
  xr <- xr[clade %in% xn$clade]; e <- predict(mn, xr, type = "response")
  bs <- replicate(500, { i <- sample(nrow(xr), replace = TRUE); mean(xr$y[i]) - mean(e[i]) })
  CAL[[length(CAL) + 1]] <- data.table(ind = ind, n_real = nrow(xr), obs = mean(xr$y), exp = mean(e), obs_minus_exp = mean(xr$y) - mean(e),
    lo = quantile(bs, .025), hi = quantile(bs, .975), ratio = mean(xr$y) / mean(e))
  r <- fitm(cc[is_real == 1], paste("fg_depth + d_omega +", F1), c("fg_depth", "d_omega"))
  if (!is.null(r)) DOSE[[length(DOSE) + 1]] <- r[, ind := ind]
}
MOD <- rbindlist(MOD, fill = TRUE); fwrite(MOD, file.path(O, "models_R1.tsv"), sep = "\t")
fwrite(rbindlist(CAL), file.path(O, "calibration_R1.tsv"), sep = "\t"); fwrite(rbindlist(DOSE), file.path(O, "dose_R1.tsv"), sep = "\t")
# prestated decision, relax_relaxed, common frame
g <- function(fr, m) MOD[ind == "relax_relaxed" & frame == fr & model == m & term == "is_real"]
m1 <- g("common", "M1"); m5 <- g("common", "M5"); h1 <- g("HQ", "M1")
ratio <- m5$est / m1$est
hq_ok <- nrow(h1) == 1 && h1$est > 0 && (h1$p_cl < 0.05 || h1$n_hog < 100)
scen <- if (ratio >= 0.5 && m5$p_cl < 0.05 && hq_ok) "A" else if (ratio < 0.25 || m5$p_cl >= 0.05) "C" else "B"
lg("relax_relaxed is_real logOR: M1 ", signif(m1$est, 3), " (p ", signif(m1$p_cl, 3), "), M5 ", signif(m5$est, 3), " (p ",
   signif(m5$p_cl, 3), "), ratio ", signif(ratio, 3), "; HQ M1 ", if (nrow(h1)) paste0(signif(h1$est, 3), " (p ", signif(h1$p_cl, 3),
   ", ", h1$n_hog, " HOGs)") else "NA", " -> scenario ", scen)
writeLines(scen, file.path(O, "SCENARIO.txt"))
# lineage-level PGLS: log omega_tip ~ D, per clade, Brownian on the pruned LSD tree; inverse-variance meta
tc <- fread(file.path(O, "tip_cov.tsv"), na.strings = c("NA", "nan")); lab <- fread(file.path(R, "meta/labels_v2.tsv"))
tr <- read.tree(file.path(R, "work/final_v3/date/lsd.timetree.nwk"))
tc <- merge(tc, lab[, .(code, lab, clade_lab = clade)], by = "code")[lab %in% c("D", "S") & is.finite(omega_tip) & omega_tip > 0 & n_branch >= 100]
PG <- list()
for (cd in unique(tc$clade_dir)) {
  x <- as.data.frame(tc[clade_dir == cd]); x <- x[x$code %in% tr$tip.label, ]
  if (sum(x$lab == "D") < 3 || sum(x$lab == "S") < 3) next
  t2 <- keep.tip(tr, x$code); rownames(x) <- x$code; x$deep <- as.integer(x$lab == "D"); x$lw <- log(x$omega_tip)
  for (v in c("lw", "gc3", "buscoC")) {
    f <- tryCatch(gls(as.formula(paste(v, "~ deep")), data = x, correlation = corBrownian(1, t2, form = ~code), method = "ML"),
                  error = function(e) NULL)
    if (is.null(f)) next; s <- summary(f)$tTable
    PG[[length(PG) + 1]] <- data.table(clade_dir = cd, var = v, n = nrow(x), nD = sum(x$deep), est = s["deep", 1], se = s["deep", 2], p = s["deep", 4],
      raw_diff = mean(x[[v]][x$deep == 1]) - mean(x[[v]][x$deep == 0]))
  }
}
PG <- rbindlist(PG)
if (nrow(PG)) {
  MA <- PG[, .(k = .N, est = sum(est / se^2) / sum(1 / se^2), se = sqrt(1 / sum(1 / se^2))), by = var][, p := 2 * pnorm(-abs(est / se))]
  fwrite(PG, file.path(O, "pgls_tip.tsv"), sep = "\t"); fwrite(MA, file.path(O, "pgls_meta.tsv"), sep = "\t")
  for (i in seq_len(nrow(MA))) lg("PGLS meta ", MA$var[i], " deep effect ", signif(MA$est[i], 3), " (SE ", signif(MA$se[i], 3), ", p ", signif(MA$p[i], 3), ", k ", MA$k[i], ")")
}
# threonine origin contrasts vs GC3 / GC
og <- fread(file.path(R, "work/addon/A2/origins.tsv")); ct <- fread(file.path(R, "work/addon/A2/contrasts.tsv"))
ct <- ct[trait == "busco_f_T" & adj == "raw"]
tv <- tc[, .(code, gc3, genome_gc, buscoC)]; setkey(tv, code)
gm <- function(s, v) { s <- strsplit(s, ",")[[1]]; mean(tv[J(s)][[v]], na.rm = TRUE) }
og[, `:=`(dgc3 = mapply(function(a, b) gm(a, "gc3") - gm(b, "gc3"), D_tips, S_tips),
          dgc = mapply(function(a, b) gm(a, "genome_gc") - gm(b, "genome_gc"), D_tips, S_tips),
          minC = mapply(function(a, b) min(tv[J(c(strsplit(a, ",")[[1]], strsplit(b, ",")[[1]]))]$buscoC, na.rm = TRUE), D_tips, S_tips))]
z <- merge(ct, og[, .(set, node, clade, dgc3, dgc, minC)], by = c("set", "node", "clade"))
TH <- list()
for (s in unique(z$set)) for (sub in c("all", "HQ")) {
  x <- z[set == s & is.finite(dgc3) & is.finite(dgc)]; if (sub == "HQ") x <- x[minC >= 80]
  if (nrow(x) < 8) next
  for (f in c("contrast ~ dgc3", "contrast ~ dgc3 + dgc")) {
    m <- lm(as.formula(f), data = x); s1 <- summary(m)$coefficients; adjc <- x$contrast - (fitted(m) - coef(m)[1])
    TH[[length(TH) + 1]] <- data.table(set = s, subset = sub, model = f, n = nrow(x), raw_median = median(x$contrast),
      n_pos = sum(x$contrast > 0), intercept = s1[1, 1], int_p = s1[1, 4], slope_gc3 = s1["dgc3", 1], slope_gc3_p = s1["dgc3", 4],
      adj_median = median(adjc), adj_n_pos = sum(adjc > 0), adj_wilcox_p = wilcox.test(adjc)$p.value)
  }
}
TH <- rbindlist(TH); fwrite(TH, file.path(O, "threonine_gc_R1.tsv"), sep = "\t")
for (i in seq_len(nrow(TH))) lg("Thr ", TH$set[i], "/", TH$subset[i], " ", TH$model[i], ": n ", TH$n[i], " intercept ", signif(TH$intercept[i], 3),
  " (p ", signif(TH$int_p[i], 3), "), adj pos ", TH$adj_n_pos[i], "/", TH$n[i], " Wilcoxon p ", signif(TH$adj_wilcox_p[i], 3))
writeLines(REP, file.path(O, "report_R1.txt")); file.create(file.path(O, "R1B.DONE"))
