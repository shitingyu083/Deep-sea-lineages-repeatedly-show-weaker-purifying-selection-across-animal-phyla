#!/usr/bin/env Rscript
# A3 (PREREG_addon.md): leave-one-clade-out (LOCO) prediction of deep-sea habitat from proteome traits.
# Usage: Rscript 56_loco_predict.R [TREE] [OUTDIR]     (default work/final_v3/date/lsd.timetree.nwk, work/addon/A3)
# env LABEL_COL (lab), DATED (1), MA (1), NSIM (1000), SEED (20261025), MIN_DS (5)
# Species = 44's (addon_lib load_species) with every M_aa feature finite. Features: M_aa = 44's traits (busco_* /
#   whole_*, no *_len) without the GC / codon family; M_all = all 44 traits (GC / codon NA -> training medians);
#   M_aa_res = M_aa residualised on route by least squares in the training set; M_route = one-hot route. Route coding
#   per fold: route (annot / miniprot / tsa), annot2 if the held-out clade has a route level absent from training.
# Each clade with >= MIN_DS D and >= MIN_DS S held out in turn: ridge logistic regression (glmnet alpha 0,
#   standardised), lambda.min of the inner cross-validation with folds = training clades (binomial deviance).
# Per held-out clade x model: AUC with the DeLong 95 % CI; calibrated p = (1 + #AUC_sim >= AUC_obs) / (1 + N_valid)
#   over NSIM labellings simulated on the held-out subtree under its Mk fit (addon_lib mk2_fit ARD, ER if ARD fails;
#   mk2_sim_cond: deep count within 10 %, 20 % if fewer than 1 % are kept), predictions fixed, the same labellings for
#   every model. Combined per model: weighted Stouffer (w = sqrt n species; p capped at NSIM / (NSIM + 1) for z);
#   pooled AUC of all held-out predictions (descriptive). Decision (PREREG A3): predictable = Stouffer p(M_aa_res)
#   < 0.05 and M_route AUC >= M_aa_res AUC in at most half of the clades.
# Fold ALL: the same fits on all species (descriptive coefficient signature, no evaluation).
# Outputs (plain tsv): features.tsv, predictions.tsv, folds.tsv, coef.tsv, combined.tsv, decision.tsv, report.txt
fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(Sys.getenv("ADDON_LIB", file.path(dirname(normalizePath(sub("^--file=", "", fa[1]))), "addon_lib.R")))
suppressPackageStartupMessages(library(glmnet))
a <- commandArgs(TRUE)
TREE <- if (length(a) >= 1) a[1] else file.path(W, "date/lsd.timetree.nwk")
OUT <- if (length(a) >= 2) a[2] else file.path(R, "work/addon/A3")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (!nzchar(Sys.getenv("DATED"))) Sys.setenv(DATED = "1")
NSIM <- as.integer(Sys.getenv("NSIM", "1000")); SEED <- as.integer(Sys.getenv("SEED", "20261025"))
MIN_DS <- as.integer(Sys.getenv("MIN_DS", "5"))
lg <- function(...) { cat(format(Sys.time(), "%F %T"), ..., "\n"); flush(stdout()) }
tt <- if (Sys.getenv("MA", "1") == "1") load_tree(TREE, ma = TRUE) else read.tree(TREE)
L <- load_species(tt); sp <- L$sp; tr <- L$tr; rm(L)
TR44 <- setdiff(grep("^(busco|whole)_", names(sp), value = TRUE), c("busco_len", "whole_len"))
F_AA <- setdiff(TR44, GCFAM); F_GC <- intersect(TR44, GCFAM)
for (v in TR44) set(sp, j = v, value = suppressWarnings(as.numeric(sp[[v]])))
ok <- Reduce(`&`, lapply(F_AA, function(v) is.finite(sp[[v]]))) & !is.na(sp$route)
lg("labels", Sys.getenv("LABEL_COL", "lab"), ": species", nrow(sp), "with every M_aa feature", sum(ok), "( D", sum(sp$deep[ok]),
   "S", sum(1 - sp$deep[ok]), "); features M_aa", length(F_AA), "GC / codon", length(F_GC))
sp <- sp[ok]; tr <- keep.tip(tr, sp$code)
fwrite(data.table(feature = TR44, family = fifelse(TR44 %in% F_AA, "aa", "gc_codon"),
                  n_missing = vapply(TR44, function(v) sum(!is.finite(sp[[v]])), 0L)), file.path(OUT, "features.tsv"), sep = "\t")
cl <- sp[, .(n = .N, n_D = sum(deep), n_S = sum(1 - deep)), by = clade][order(-n)]
HO <- cl[n_D >= MIN_DS & n_S >= MIN_DS, clade]
lg("held-out clades:", paste(cl[clade %in% HO, sprintf("%s (D %d, S %d)", clade, n_D, n_S)], collapse = ", "))
if (length(HO) < 2) stop("fewer than 2 clades qualify for leave-one-clade-out")
MODELS <- c("M_aa", "M_aa_res", "M_all", "M_route")

# design matrices of one model, training rows itr, held-out rows iho (route coding per fold, see header)
mk_x <- function(model, itr, iho) {
  rv <- if (all(sp$route[iho] %in% sp$route[itr])) "route" else "annot2"
  lev <- sort(unique(sp[[rv]][itr]))
  oh <- function(i) { f <- factor(sp[[rv]][i], levels = lev); stopifnot(!anyNA(f))
    m <- matrix(0, length(i), length(lev), dimnames = list(NULL, lev)); m[cbind(seq_along(i), as.integer(f))] <- 1; m }
  if (model == "M_route") return(list(tr = oh(itr), ho = oh(iho), rv = rv))
  fs <- if (model == "M_all") TR44 else F_AA
  X <- as.matrix(sp[, ..fs]); Xt <- X[itr, , drop = FALSE]; Xh <- X[iho, , drop = FALSE]
  if (model == "M_all") {
    keep <- colSums(is.finite(Xt)) > 0; Xt <- Xt[, keep, drop = FALSE]; Xh <- Xh[, keep, drop = FALSE]
    for (j in intersect(F_GC, colnames(Xt))) { md <- median(Xt[, j], na.rm = TRUE)
      Xt[!is.finite(Xt[, j]), j] <- md; Xh[!is.finite(Xh[, j]), j] <- md } }
  if (model == "M_aa_res") { Dt <- cbind(1, oh(itr)[, -1, drop = FALSE]); Dh <- cbind(1, oh(iho)[, -1, drop = FALSE])
    B <- qr.coef(qr(Dt), Xt); B[is.na(B)] <- 0; Xt <- Xt - Dt %*% B; Xh <- Xh - Dh %*% B }
  list(tr = Xt, ho = Xh, rv = rv)
}
ridge <- function(X, y, fold) {                      # glmnet needs >= 2 columns: a constant column is ignored by it
  if (ncol(X$tr) < 2) { X$tr <- cbind(X$tr, zero_ = 0); X$ho <- cbind(X$ho, zero_ = 0) }
  cv <- suppressWarnings(cv.glmnet(X$tr, y, family = "binomial", alpha = 0, foldid = fold, type.measure = "deviance"))
  b <- as.matrix(coef(cv, s = "lambda.min"))[, 1]; sdv <- c(NA, apply(X$tr, 2, sd))
  list(p = if (nrow(X$ho)) as.numeric(predict(cv, newx = X$ho, s = "lambda.min", type = "response")) else numeric(0),
       lambda = cv$lambda.min, coef = data.table(feature = names(b), coef = b, coef_std = b * sdv))
}
# AUC (Mann-Whitney, ties 1/2) with the DeLong standard error and 95 % CI (clipped to [0, 1])
delong <- function(s, y) {
  x1 <- s[y == 1]; x0 <- s[y == 0]; n1 <- length(x1); n0 <- length(x0)
  psi <- outer(x1, x0, ">") + 0.5 * outer(x1, x0, "==")
  v10 <- rowMeans(psi); v01 <- colMeans(psi); auc <- mean(v10); se <- sqrt(var(v10) / n1 + var(v01) / n0)
  c(auc = auc, se = se, lo = max(0, auc - qnorm(0.975) * se), hi = min(1, auc + qnorm(0.975) * se))
}
auc_sims <- function(s, Y) {                         # AUC of fixed scores for every simulated labelling (columns of Y)
  r <- rank(s); n1 <- colSums(Y); n0 <- length(s) - n1
  v <- (colSums(r * Y) - n1 * (n1 + 1) / 2) / (n1 * n0); v[n1 == 0 | n0 == 0] <- NA; v
}
stouffer <- function(p, w) { z <- qnorm(pmin(p, NSIM / (NSIM + 1)), lower.tail = FALSE); Z <- sum(w * z) / sqrt(sum(w^2))
  c(Z = Z, p = pnorm(Z, lower.tail = FALSE)) }

PR <- list(); FO <- list(); CO <- list()
for (k in seq_along(HO)) {
  h <- HO[k]; iho <- which(sp$clade == h); itr <- which(sp$clade != h); y <- sp$deep[iho]
  fold <- as.integer(factor(sp$clade[itr]))
  st <- keep.tip(tr, sp$code[iho]); x <- setNames(y, sp$code[iho])
  mk <- mk2_fit(st, x, "ARD"); sm <- mk2_sim_cond(st, mk$q_SD, mk$q_DS, NSIM, sum(y), seed = SEED + k)
  Y <- sm$sims[sp$code[iho], , drop = FALSE]
  lg(h, ": n", length(iho), "D", sum(y), "| Mk", mk$model, "q_SD", signif(mk$q_SD, 3), "q_DS", signif(mk$q_DS, 3),
     "| null labellings", ncol(Y), "tol", sm$tol, "acceptance", signif(sm$acc_rate, 3))
  pr <- data.table(code = sp$code[iho], clade = h, deep = y, route = sp$route[iho])
  for (m in MODELS) {
    X <- mk_x(m, itr, iho); f <- ridge(X, sp$deep[itr], fold); pr[, (m) := f$p]
    dl <- delong(f$p, y); asim <- auc_sims(f$p, Y); nv <- sum(!is.na(asim))
    FO[[length(FO) + 1]] <- data.table(clade = h, model = m, n = length(iho), n_D = sum(y), n_S = sum(1 - y), route_coding = X$rv,
      n_features = ncol(X$tr), lambda = f$lambda, auc = dl[["auc"]], auc_se = dl[["se"]], auc_lo = dl[["lo"]], auc_hi = dl[["hi"]],
      null_mean_auc = mean(asim, na.rm = TRUE), n_null = nv, p_sim = (1 + sum(asim >= dl[["auc"]] - 1e-9, na.rm = TRUE)) / (1 + nv),
      mk_model = mk$model, q_SD = mk$q_SD, q_DS = mk$q_DS, sim_tol = sm$tol, sim_acc = sm$acc_rate)
    CO[[length(CO) + 1]] <- f$coef[, `:=`(fold = h, model = m)]
  }
  PR[[k]] <- pr
}
iall <- seq_len(nrow(sp)); fall <- as.integer(factor(sp$clade))
for (m in MODELS) CO[[length(CO) + 1]] <- ridge(mk_x(m, iall, iall), sp$deep, fall)$coef[, `:=`(fold = "ALL", model = m)]
PR <- rbindlist(PR); FO <- rbindlist(FO); CO <- rbindlist(CO)
fwrite(PR, file.path(OUT, "predictions.tsv"), sep = "\t"); fwrite(FO, file.path(OUT, "folds.tsv"), sep = "\t")
fwrite(CO[, .(fold, model, feature, coef, coef_std)], file.path(OUT, "coef.tsv"), sep = "\t")
rte <- FO[model == "M_route", .(clade, auc_route = auc)]
CB <- rbindlist(lapply(MODELS, function(m) {
  d <- merge(FO[model == m], rte, by = "clade"); s <- stouffer(d$p_sim, sqrt(d$n)); pl <- delong(PR[[m]], PR$deep)
  data.table(model = m, n_clades = nrow(d), stouffer_Z = s[["Z"]], stouffer_p = s[["p"]], median_auc = median(d$auc),
             min_auc = min(d$auc), max_auc = max(d$auc), n_p05 = sum(d$p_sim < 0.05), n_route_ge = sum(d$auc_route >= d$auc),
             pooled_auc = pl[["auc"]], pooled_lo = pl[["lo"]], pooled_hi = pl[["hi"]], n_pred = nrow(PR))
}))
fwrite(CB, file.path(OUT, "combined.tsv"), sep = "\t")
r <- CB[model == "M_aa_res"]
DEC <- data.table(stouffer_p_M_aa_res = r$stouffer_p, n_clades = r$n_clades, n_route_ge_M_aa_res = r$n_route_ge,
                  route_majority = r$n_route_ge > r$n_clades / 2, predictable = r$stouffer_p < 0.05 & !(r$n_route_ge > r$n_clades / 2))
fwrite(DEC, file.path(OUT, "decision.tsv"), sep = "\t")
rp <- c(sprintf("A3 leave-one-clade-out prediction, labels %s, %d species (D %d, S %d), %d held-out clades, NSIM %d",
                Sys.getenv("LABEL_COL", "lab"), nrow(sp), sum(sp$deep), sum(1 - sp$deep), length(HO), NSIM),
        FO[, sprintf("  %-15s %-8s n %3d (D %3d) AUC %.3f [%.3f, %.3f] null %.3f p_sim %.4f", clade, model, n, n_D, auc, auc_lo,
                     auc_hi, null_mean_auc, p_sim)],
        CB[, sprintf("combined %-8s Stouffer Z %.2f p %.2g; median AUC %.3f (%.3f-%.3f); %d/%d clades p_sim < 0.05; M_route >= in %d; pooled %.3f [%.3f, %.3f]",
                     model, stouffer_Z, stouffer_p, median_auc, min_auc, max_auc, n_p05, n_clades, n_route_ge, pooled_auc, pooled_lo, pooled_hi)],
        sprintf("decision: predictable = %s (Stouffer p M_aa_res %.2g; M_route AUC >= M_aa_res AUC in %d of %d clades)",
                DEC$predictable, DEC$stouffer_p_M_aa_res, DEC$n_route_ge_M_aa_res, DEC$n_clades))
writeLines(rp, file.path(OUT, "report.txt")); cat(rp, sep = "\n")
lg("done")
