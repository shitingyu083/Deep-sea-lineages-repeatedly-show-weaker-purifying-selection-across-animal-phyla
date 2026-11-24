#!/usr/bin/env Rscript
# EXPLORATORY (not preregistered): R1 HQ-subset sensitivity. The prestated HQ rule in 62b returned 0 HOGs because
#   bg_nonannot is NA for most clades and min BUSCO C >= 80 holds for 16.8% of rows; the confirmatory scenario stays B.
#   Here: relaxed HQ definitions, M0/M1/M5 on relax_relaxed and relax, common frame. Output work/addon/R1/explore_hq.
#   Model code copied verbatim from 62b (lines 9-40); 62b itself is frozen and untouched.
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
OX <- file.path(O, "explore_hq"); dir.create(OX, showWarnings = FALSE)
lg("EXPLORATORY HQ sensitivity (not preregistered); prestated HQ rule gave 0 HOGs -> confirmatory scenario remains B")
MOD <- list(); NH <- list()
for (ind in c("relax_relaxed", "relax")) {
  x0 <- d[!is.na(get(paste0("y_", ind)))][, y := get(paste0("y_", ind))]
  cc <- paired(x0[complete.cases(x0[, .(d_omega, d_gc, d_gc3, d_buscoC, d_tbl, d_temp, fg_nonannot, gapdiff)])])
  hm <- cc[, .(mb = min(min_buscoC), fa = all(fg_nonannot == 0), ba = all(is.na(bg_nonannot) | bg_nonannot == 0), clade = clade[1]), by = hid]
  hm[, cmed := median(mb, na.rm = TRUE), by = clade]
  DEF <- list(H_busco80 = hm[mb >= 80]$hid,
              H_busco80_fg0 = hm[mb >= 80 & fa]$hid,
              H_fg0_bgNAok_busco80 = hm[mb >= 80 & fa & ba]$hid,
              H_busco70 = hm[mb >= 70]$hid,
              H_cladeTopHalf = hm[mb >= cmed]$hid,
              H_cladeTopHalf_fg0 = hm[mb >= cmed & fa]$hid)
  for (h in names(DEF)) {
    xh <- cc[hid %in% DEF[[h]]]
    NH[[length(NH) + 1]] <- data.table(ind = ind, def = h, n_hog = uniqueN(xh$hid), n_clade = uniqueN(xh$clade),
      clades = paste(xh[, .N, by = clade][order(clade), paste0(clade, ":", N / 2)], collapse = ","))
    for (s in c("M0", "M1", "M5")) { r <- fitm(xh, SPEC[[s]], TERMS); if (!is.null(r)) MOD[[length(MOD) + 1]] <- r[, `:=`(ind = ind, frame = h, model = s)] }
  }
}
MOD <- rbindlist(MOD, fill = TRUE); NH <- rbindlist(NH)
fwrite(MOD, file.path(OX, "models_hq_explore.tsv"), sep = "\t"); fwrite(NH, file.path(OX, "hq_sizes.tsv"), sep = "\t")
for (i in seq_len(nrow(NH))) lg(NH$ind[i], " ", NH$def[i], ": ", NH$n_hog[i], " HOGs, ", NH$n_clade[i], " clades (", NH$clades[i], ")")
for (i in which(MOD$term == "is_real")) lg(MOD$ind[i], " ", MOD$frame[i], " ", MOD$model[i], " is_real logOR ", signif(MOD$est[i], 3),
  " (cl SE ", signif(MOD$se_cl[i], 3), ", p ", signif(MOD$p_cl[i], 3), "; ", MOD$n_hog[i], " HOGs, ", MOD$n_events[i], " events)")
writeLines(REP, file.path(OX, "report_hq_explore.txt")); file.create(file.path(OX, "HQX.DONE"))
