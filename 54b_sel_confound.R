#!/usr/bin/env Rscript
# A4 (PREREG_addon.md): selection-test confounds and false-positive rates, on 54a's sel_features.tsv.
# Usage: Rscript 54b_sel_confound.R [FEATURES] [OUTDIR] [FINAL_MAIN]
#   defaults work/addon/A4/sel_features.tsv, work/addon/A4, work/combine25c/final_main.combined.tsv (25c primary)
# Indicators: relax (RELAX p < 0.05), relax_relaxed / relax_intens (and K < 1 / K > 1), bph (BUSTED-PH FG p < 0.05),
#   bph_strict (and difference p < 0.05; NA when either is NA, as 25b), cs_perm (CSUBST permutation p < 0.05),
#   cs_qual (>= 1 qualified pair; the indicator form of 25c "CSUBST qualified pairs").
# Frames = the paired frames of 25b / 25c: HyPhy, nullhy HOGs with a real and a null1 value (rows real, null1);
#   CSUBST, hylist HOGs with a real and >= 1 null value (rows real and every null). HOG id = clade + hog.
# Models per indicator (y = 0/1): M1 (primary) glm y ~ is_real + fg_nonannot + log(n_tips) + gapdiff + clade,
#   binomial, HOG-clustered sandwich SE (CR1: G/(G-1) (N-1)/(N-K)); M0 y ~ is_real + clade; M2 = M1 + bg_nonannot +
#   fg_tsa + log(aln_len); C1 conditional logistic y ~ is_real + fg_nonannot + gapdiff + strata(HOG) (survival::clogit,
#   within-HOG contrast); per clade M1 without clade (>= 5 events); N1 route effect within the null rows
#   (y ~ fg_nonannot + log(n_tips) + gapdiff + clade, clustered SE). Clades without an event in a frame leave the
#   pooled fits (their intercept diverges; they carry no information on is_real given clade).
# Route-matched pairs: per HOG the real row and its null row with the smallest |fg_nonannot difference| (HyPhy: null1;
#   ties: lowest null index), kept when <= MATCH_TOL (0.2); exact McNemar (two-sided sign test on the discordant
#   pairs); the same on all pairs (CSUBST: real vs null1) for reference.
# Gap tertiles: per HOG d = y_real - mean(y_null), tertiles of the real row's gapdiff (pooled); mean d, sign test.
# Null false-positive rate per clade x indicator: the null rows of the frame, Wilson 95 % CI, flag when > 0.10.
# Decision (each 25c primary test with p_holm < 0.05): keep when M1 is_real > 0 with clustered p < 0.05, else the
#   main-text wording of that test is downgraded.
# Outputs: feature_summary.tsv, models.tsv, perclade.tsv, nullroute.tsv, matched.tsv, gap_tertiles.tsv,
#   null_fpr.tsv, decision.tsv, report.txt
suppressPackageStartupMessages(library(data.table))
HAS_SURV <- suppressWarnings(suppressPackageStartupMessages(require(survival, quietly = TRUE)))
a <- commandArgs(TRUE)
R <- Sys.getenv("DSC_ROOT", path.expand("~/deepsea_conv"))
FEAT <- if (length(a) >= 1) a[1] else file.path(R, "work/addon/A4/sel_features.tsv")
OUT <- if (length(a) >= 2) a[2] else file.path(R, "work/addon/A4")
FM <- if (length(a) >= 3) a[3] else file.path(R, "work/combine25c/final_main.combined.tsv")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TOL <- as.numeric(Sys.getenv("MATCH_TOL", "0.2")); FLAG <- 0.10
REP <- character(0)
lg <- function(...) { s <- paste0(...); cat(format(Sys.time(), "%F %T"), s, "\n"); flush(stdout()); REP <<- c(REP, s) }
d <- fread(FEAT, sep = "\t", na.strings = c("NA", "nan", "NaN", ""))
d[, is_real := as.integer(set == "real")]
d[, nidx := fifelse(set == "real", 0L, suppressWarnings(as.integer(sub("^null", "", set))))]
d[, hid := paste(clade, hog, sep = ":")]
lt <- function(x) fifelse(is.na(x), NA_integer_, as.integer(x < 0.05))
d[, `:=`(y_relax = lt(relax_p), y_bph = lt(bph_p), y_cs_perm = lt(cs_perm_p),
         y_cs_qual = fifelse(is.na(cs_nqual), NA_integer_, as.integer(cs_nqual >= 1)))]
d[, y_bph_strict := fifelse(is.na(y_bph) | is.na(lt(bph_diff_p)), NA_integer_, y_bph * lt(bph_diff_p))]
d[, y_relax_relaxed := fifelse(is.na(y_relax) | is.na(relax_K), NA_integer_, as.integer(y_relax == 1L & relax_K < 1))]
d[, y_relax_intens := fifelse(is.na(y_relax) | is.na(relax_K), NA_integer_, as.integer(y_relax == 1L & relax_K > 1))]
IND <- c(relax = "hy", relax_relaxed = "hy", relax_intens = "hy", bph = "hy", bph_strict = "hy", cs_perm = "cs",
         cs_qual = "cs")
T25 <- c(relax = "RELAX p<0.05", bph = "BUSTED-PH FG p<0.05", bph_strict = "BUSTED-PH strict (FG p<0.05 and diff p<0.05)",
         cs_perm = "CSUBST perm p<0.05", cs_qual = "CSUBST qualified pairs")
lg("A4 54b: ", nrow(d), " rows, ", uniqueN(d$hid), " HOGs, clades ", paste(sort(unique(d$clade)), collapse = " "),
   "; fg_match = 0 in ", sum(d$fg_match == 0, na.rm = TRUE), " rows; survival ", HAS_SURV)
frame_of <- function(ind) {   # the 25b / 25c paired frame of an indicator; y = its 0/1 value
  x <- d[!is.na(get(paste0("y_", ind)))][, y := get(paste0("y_", ind))]
  x <- if (IND[[ind]] == "hy") x[in_nullhy == 1 & set %in% c("real", "null1")] else x[in_hylist == 1]
  ok <- x[, .(r = any(is_real == 1L), n = any(is_real == 0L)), by = hid][r & n]$hid
  x[hid %in% ok]
}
cl_vc <- function(m, g) {     # HOG-clustered sandwich covariance (CR1) of a binomial-logit glm
  X <- model.matrix(m)[, !is.na(coef(m)), drop = FALSE]
  S <- rowsum(X * (m$y - fitted(m)), g); B <- summary(m)$cov.unscaled
  G <- nrow(S); N <- nrow(X); K <- ncol(X)
  B %*% crossprod(S) %*% B * (G / (G - 1)) * ((N - 1) / (N - K))
}
fit_glm <- function(x, rhs, terms = "is_real") {
  v <- setdiff(all.vars(as.formula(paste("~", rhs))), "clade")
  x <- x[complete.cases(x[, c("y", v), with = FALSE])]
  if ("aln_len" %in% v) x <- x[aln_len > 0]
  if (grepl("clade", rhs)) {
    x <- x[clade %in% x[, .(e = sum(y)), by = clade][e > 0]$clade]
    if (uniqueN(x$clade) < 2) rhs <- sub("\\s*\\+\\s*clade", "", rhs)
  }
  if (nrow(x) < 10 || sum(x$y) == 0 || sum(x$y) == nrow(x)) return(NULL)
  w <- character(0)
  m <- withCallingHandlers(glm(as.formula(paste("y ~", rhs)), binomial, data = x),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  cf <- coef(summary(m)); V <- cl_vc(m, x$hid)
  rbindlist(lapply(intersect(terms, rownames(cf)), function(t) {
    se <- sqrt(V[t, t]); data.table(term = t, est = cf[t, 1], se = cf[t, 2], p = cf[t, 4], se_cl = se,
      p_cl = 2 * pnorm(-abs(cf[t, 1] / se)), n_rows = nrow(x), n_hog = uniqueN(x$hid), n_events = sum(x$y),
      rate_real = mean(x$y[x$is_real == 1L]), rate_null = mean(x$y[x$is_real == 0L]), n_clades = uniqueN(x$clade),
      converged = m$converged, warn = if (length(w)) substr(paste(unique(w), collapse = "; "), 1, 120) else "")
  }))
}
fit_clogit <- function(x) {   # within-HOG contrast; strata need both set types and variation in y
  if (!HAS_SURV) return(NULL)
  x <- x[complete.cases(x[, .(y, fg_nonannot, gapdiff)])]
  x <- x[hid %in% x[, .(k = uniqueN(y) > 1 && any(is_real == 1L) && any(is_real == 0L)), by = hid][k == TRUE]$hid]
  if (uniqueN(x$hid) < 10) return(NULL)
  m <- tryCatch(suppressWarnings(clogit(y ~ is_real + fg_nonannot + gapdiff + strata(hid), data = x)),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- coef(summary(m))
  rbindlist(lapply(intersect(c("is_real", "fg_nonannot", "gapdiff"), rownames(cf)), function(t)
    data.table(term = t, est = cf[t, "coef"], se = cf[t, "se(coef)"], p = cf[t, "Pr(>|z|)"], se_cl = NA_real_,
               p_cl = NA_real_, n_rows = nrow(x), n_hog = uniqueN(x$hid), n_events = sum(x$y),
               rate_real = mean(x$y[x$is_real == 1L]), rate_null = mean(x$y[x$is_real == 0L]),
               n_clades = uniqueN(x$clade), converged = TRUE, warn = "")))
}
wilson <- function(k, n, z = qnorm(0.975)) {
  p <- k / n; den <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / den
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  list(lo = c0 - h, hi = c0 + h)
}
sign_p <- function(np, nn) if (np + nn > 0) binom.test(np, np + nn)$p.value else NA_real_
F1 <- "is_real + fg_nonannot + log(n_tips) + gapdiff + clade"
SPEC <- list(M1 = F1, M0 = "is_real + clade", M2 = paste(F1, "+ bg_nonannot + fg_tsa + log(aln_len)"))
TERMS <- c("is_real", "fg_nonannot", "gapdiff", "log(n_tips)", "bg_nonannot", "fg_tsa")
MOD <- list(); PC <- list(); NR <- list(); MP <- list(); GT <- list(); FP <- list()
for (ind in names(IND)) {
  x <- frame_of(ind)
  lg(ind, " (", IND[[ind]], " frame): ", uniqueN(x$hid), " HOGs, ", nrow(x), " rows; real rate ",
     signif(mean(x$y[x$is_real == 1L]), 3), ", null rate ", signif(mean(x$y[x$is_real == 0L]), 3))
  if (!nrow(x)) next
  for (s in names(SPEC)) { r <- fit_glm(x, SPEC[[s]], TERMS); if (!is.null(r)) MOD[[length(MOD) + 1]] <- r[, `:=`(ind = ind, model = s)] }
  r <- fit_clogit(x); if (!is.null(r)) MOD[[length(MOD) + 1]] <- r[, `:=`(ind = ind, model = "C1")]
  for (cl in sort(unique(x$clade))) {
    xc <- x[clade == cl]
    if (sum(xc$y) < 5) { PC[[length(PC) + 1]] <- data.table(ind = ind, clade = cl, term = "is_real", n_events = sum(xc$y), note = "< 5 events"); next }
    r <- fit_glm(xc, sub(" + clade", "", F1, fixed = TRUE), c("is_real", "fg_nonannot"))
    if (!is.null(r)) PC[[length(PC) + 1]] <- r[, `:=`(ind = ind, clade = cl, note = "")]
  }
  r <- fit_glm(x[is_real == 0L], "fg_nonannot + log(n_tips) + gapdiff + clade", c("fg_nonannot", "gapdiff", "log(n_tips)"))
  if (!is.null(r)) NR[[length(NR) + 1]] <- r[, ind := ind]
  pr <- merge(x[is_real == 1L, .(hid, clade, yr = y, fr = fg_nonannot, gr = gapdiff)],
              x[is_real == 0L, .(hid, nidx, yn = y, fn = fg_nonannot)], by = "hid")
  pr[, dlt := abs(fr - fn)]
  pm <- pr[!is.na(dlt)][order(hid, dlt, nidx)][, .SD[1], by = hid][dlt <= TOL]
  for (k in c("matched", "all")) {
    q <- if (k == "matched") pm else pr[order(hid, nidx)][, .SD[1], by = hid]
    b <- sum(q$yr == 1L & q$yn == 0L); cc <- sum(q$yr == 0L & q$yn == 1L)
    MP[[length(MP) + 1]] <- data.table(ind = ind, pairs = k, n_pairs = nrow(q), real_only = b, null_only = cc,
      p_mcnemar = sign_p(b, cc), mean_delta_fg_nonannot = mean(q$fr - q$fn, na.rm = TRUE), tol = TOL)
  }
  h <- x[, .(dd = mean(y[is_real == 1L]) - mean(y[is_real == 0L]), g = gapdiff[is_real == 1L][1]), by = hid][!is.na(g)]
  if (nrow(h) >= 3) {
    br <- unique(quantile(h$g, c(0, 1 / 3, 2 / 3, 1)))
    h[, tert := if (length(br) > 1) cut(g, br, include.lowest = TRUE, labels = FALSE) else 1L]
    GT[[length(GT) + 1]] <- h[, .(n_hog = .N, gap_lo = min(g), gap_hi = max(g), mean_d = mean(dd), n_pos = sum(dd > 0),
                                  n_neg = sum(dd < 0), p_sign = sign_p(sum(dd > 0), sum(dd < 0))), by = tert][order(tert)][, ind := ind]
  }
  z <- rbind(x[is_real == 0L, .(n = .N, k = sum(y), real_rate = mean(x$y[x$is_real == 1L & x$clade == clade[1]])), by = clade],
             x[is_real == 0L, .(clade = "all", n = .N, k = sum(y), real_rate = mean(x$y[x$is_real == 1L]))])
  z[, `:=`(rate = k / n, lo = wilson(k, n)$lo, hi = wilson(k, n)$hi)][, flag := rate > FLAG]
  FP[[length(FP) + 1]] <- z[, ind := ind]
}
wr <- function(x, f) {        # always writes the file (empty when there is nothing to report)
  x <- if (is.list(x) && !is.data.table(x)) rbindlist(x, fill = TRUE) else x
  if (!ncol(x)) writeLines(character(0), file.path(OUT, f)) else fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")
}
st <- d[, .(n = .N, fg_nonannot = mean(fg_nonannot, na.rm = TRUE), bg_nonannot = mean(bg_nonannot, na.rm = TRUE),
            fg_tsa = mean(fg_tsa, na.rm = TRUE), fg_gap = mean(fg_gap, na.rm = TRUE), bg_gap = mean(bg_gap, na.rm = TRUE),
            gapdiff = mean(gapdiff, na.rm = TRUE), fg_n = mean(fg_n, na.rm = TRUE)),
        by = .(clade, set_type = fifelse(set == "real", "real", fifelse(set == "null1", "null1", "null2+")))]
wr(st[order(clade, set_type)], "feature_summary.tsv")
MD <- rbindlist(MOD, fill = TRUE)
if (!nrow(MD)) MD <- data.table(ind = character(), model = character(), term = character(), est = numeric(),
                                p = numeric(), p_cl = numeric())
setcolorder(MD, c("ind", "model", "term"))
wr(MD, "models.tsv"); wr(PC, "perclade.tsv"); wr(NR, "nullroute.tsv"); wr(MP, "matched.tsv")
wr(GT, "gap_tertiles.tsv"); wr(FP, "null_fpr.tsv")
fm <- if (file.exists(FM)) fread(FM, sep = "\t", na.strings = c("NA", "nan", "NaN", "")) else NULL
gv <- function(i, m, col) { v <- MD[ind == i & model == m & term == "is_real"][[col]]; if (length(v)) v[1] else NA_real_ }
DEC <- rbindlist(lapply(names(T25), function(i) {
  ph <- if (is.null(fm)) NA_real_ else suppressWarnings(as.numeric(fm[test == T25[[i]]]$p_holm[1]))
  e1 <- gv(i, "M1", "est"); p1 <- gv(i, "M1", "p_cl")
  data.table(test = T25[[i]], ind = i, p_holm_25c = ph, M0_logOR = gv(i, "M0", "est"), M0_p_cl = gv(i, "M0", "p_cl"),
             M1_logOR = e1, M1_OR = exp(e1), M1_p_cl = p1, M2_logOR = gv(i, "M2", "est"), M2_p_cl = gv(i, "M2", "p_cl"),
             C1_logOR = gv(i, "C1", "est"), C1_p = gv(i, "C1", "p"),
             decision = if (is.na(ph)) "25c result missing" else if (ph >= 0.05) "not significant in 25c"
                        else if (!is.na(e1) && !is.na(p1) && e1 > 0 && p1 < 0.05) "keep" else "downgrade")
}))
wr(DEC, "decision.tsv")
lg("25c final_main: ", if (is.null(fm)) paste("missing", FM) else FM)
for (i in seq_len(nrow(DEC))) with(DEC[i], lg(sprintf("%-46s 25c Holm %s | M0 logOR %.3f (p %.3g) -> M1 %.3f (OR %.2f, clustered p %.3g); C1 %.3f (p %.3g) => %s",
  test, format(signif(p_holm_25c, 3)), M0_logOR, M0_p_cl, M1_logOR, M1_OR, M1_p_cl, C1_logOR, C1_p, decision)))
NRt <- if (length(NR)) rbindlist(NR)[term == "fg_nonannot"] else data.table()
for (r in seq_len(nrow(NRt))) with(NRt[r],
  lg(sprintf("null rows, %s: fg_nonannot logOR %.3f (clustered p %.3g), %d rows", ind, est, p_cl, n_rows)))
FPt <- if (length(FP)) rbindlist(FP)[flag == TRUE] else data.table()
lg("null false-positive rate > ", FLAG, ": ", if (nrow(FPt)) paste(FPt[, sprintf("%s/%s %.3f [%.3f, %.3f]", ind, clade, rate, lo, hi)], collapse = "; ") else "none")
writeLines(REP, file.path(OUT, "report.txt"))
cat("A4_DONE\n")
