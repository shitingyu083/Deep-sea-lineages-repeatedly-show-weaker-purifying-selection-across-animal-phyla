#!/usr/bin/env Rscript
# R1 depth-threshold sensitivity (PREREG_nature.md L42-44), descriptive. Same data assembly and model code as 62b
# (62b untouched). For threshold T in 200/500/1000 m, keep HOGs whose real FG tips all have depth_median >= T
# (FG unchanged when tips between 200 m and T are recoded to U; FG tips without a depth value -> HOG dropped),
# with their null1 rows; fit M0 and M1 for relax_relaxed and relax_intensified. Output work/addon/R1/depth_thresh.
suppressPackageStartupMessages({library(data.table)})
R <- path.expand("~/deepsea_conv"); O <- file.path(R, "work/addon/R1/depth_thresh")
d <- fread(file.path(R, "work/addon/A4/sel_features.tsv"), na.strings = c("NA", "nan", "NaN", ""))
rc <- fread(file.path(R, "work/addon/R1/row_cov.tsv"), na.strings = c("NA", "nan", "NaN", ""))
fd <- fread(file.path(O, "fg_depth_rows.tsv"), na.strings = c("NA", "nan", "NaN", ""))
d <- merge(d[set %in% c("real", "null1") & in_nullhy == 1], rc, by = c("clade", "hog", "set"))
d[, is_real := as.integer(set == "real")]; d[, hid := paste(clade, hog, sep = ":")]
lt <- function(x) fifelse(is.na(x), NA_integer_, as.integer(x < 0.05))
d[, y_relax := lt(relax_p)]
d[, y_relax_relaxed := fifelse(is.na(y_relax) | is.na(relax_K), NA_integer_, as.integer(y_relax == 1L & relax_K < 1))]
d[, y_relax_intensified := fifelse(is.na(y_relax) | is.na(relax_K), NA_integer_, as.integer(y_relax == 1L & relax_K > 1))]
cl_vc <- function(m, g) {
  X <- model.matrix(m)[, !is.na(coef(m)), drop = FALSE]
  S <- rowsum(X * (m$y - fitted(m)), g); B <- summary(m)$cov.unscaled
  G <- nrow(S); N <- nrow(X); K <- ncol(X); B %*% crossprod(S) %*% B * (G / (G - 1)) * ((N - 1) / (N - K))
}
fitm <- function(x, rhs) {
  v <- setdiff(all.vars(as.formula(paste("~", rhs))), "clade")
  x <- x[complete.cases(x[, c("y", v), with = FALSE])]
  x <- x[clade %in% x[, .(e = sum(y)), by = clade][e > 0]$clade]
  if (nrow(x) < 20 || sum(x$y) == 0) return(NULL)
  m <- suppressWarnings(glm(as.formula(paste("y ~", rhs)), binomial, data = x)); cf <- coef(summary(m)); V <- cl_vc(m, x$hid)
  data.table(est = cf["is_real", 1], se_cl = sqrt(V["is_real", "is_real"]),
    p_cl = 2 * pnorm(-abs(cf["is_real", 1] / sqrt(V["is_real", "is_real"]))), n_rows = nrow(x), n_hog = uniqueN(x$hid),
    n_clades = uniqueN(x$clade), rate_real = mean(x$y[x$is_real == 1]), rate_null = mean(x$y[x$is_real == 0]))
}
paired <- function(x) x[hid %in% x[, .(k = any(is_real == 1) & any(is_real == 0)), by = hid][k == TRUE]$hid]
F1 <- "fg_nonannot + log(n_tips) + gapdiff + clade"
SPEC <- list(M0 = "is_real + clade", M1 = paste("is_real +", F1))
fd[, hid := paste(clade, hog, sep = ":")]
OUTL <- list()
for (thr in c(200, 500, 1000)) {
  keep <- fd[n_fg_na == 0 & !is.na(fg_min_depth) & fg_min_depth >= thr]$hid
  for (ind in c("relax_relaxed", "relax_intensified")) {
    x <- paired(d[hid %in% keep & !is.na(get(paste0("y_", ind)))][, y := get(paste0("y_", ind))])
    for (s in names(SPEC)) {
      r <- fitm(x, SPEC[[s]])
      if (!is.null(r)) OUTL[[length(OUTL) + 1]] <- r[, `:=`(threshold_m = thr, ind = ind, model = s)]
    }
  }
}
OUT <- rbindlist(OUTL); setcolorder(OUT, c("threshold_m", "ind", "model"))
fwrite(OUT, file.path(O, "models_depth_thresh.tsv"), sep = "\t")
print(OUT[, .(threshold_m, ind, model, OR = round(exp(est), 2), p_cl = signif(p_cl, 3), n_hog, n_clades,
              rate_real = round(rate_real, 3), rate_null = round(rate_null, 3))])
file.create(file.path(O, "DT.DONE"))
