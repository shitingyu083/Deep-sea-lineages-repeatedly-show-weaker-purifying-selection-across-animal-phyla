#!/usr/bin/env Rscript
# R5 D3 (PREREG_r5.md): route-pure pruned subsets. HOG-clustered CR1 binomial GLM (cl_vc/fitm copied from 62b):
#   R0 = is_real + clade, R1 = is_real + log(n_tips) + gapdiff + clade on the pruned rows; R0 on the unpruned rows of
#   the same HOGs (reference). Indicators relax_relaxed (primary) and relax (P < 0.05). Input R5/D3/rows.tsv (69d).
#   A single-clade fit drops the clade term (it has one level). Output R5/D3/models.tsv + verdict.tsv.
suppressPackageStartupMessages(library(data.table))
R <- path.expand("~/deepsea_conv"); O <- file.path(R, "work/addon/R5/D3")
d <- fread(file.path(O, "rows.tsv"), na.strings = c("NA", "nan", "NaN", ""))
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
  if (uniqueN(x$clade) < 2) rhs <- sub(" \\+ clade$", "", rhs)
  m <- suppressWarnings(glm(as.formula(paste("y ~", rhs)), binomial, data = x)); cf <- coef(summary(m)); V <- cl_vc(m, x$hid)
  rbindlist(lapply(intersect(terms, rownames(cf)), function(t) data.table(term = t, est = cf[t, 1], se_cl = sqrt(V[t, t]),
    p_cl = 2 * pnorm(-abs(cf[t, 1] / sqrt(V[t, t]))), n_rows = nrow(x), n_hog = uniqueN(x$hid), n_clade = uniqueN(x$clade),
    n_events = sum(x$y), rate_real = mean(x$y[x$is_real == 1]), rate_null = mean(x$y[x$is_real == 0]))))
}
paired <- function(x) x[hid %in% x[, .(k = any(is_real == 1) & any(is_real == 0)), by = hid][k == TRUE]$hid]
SPEC <- list(R0 = "is_real + clade", R1 = "is_real + log(n_tips) + gapdiff + clade")
MOD <- list()
for (rt in c("A", "M")) for (ind in c("relax_relaxed", "relax")) {
  x0 <- copy(d[route == rt]); x0[, y := get(paste0("y_", ind))]
  xp <- paired(x0[frame == "pruned" & !is.na(y)]); xu <- paired(x0[frame == "unpruned" & !is.na(y) & hid %in% xp$hid])
  xp <- xp[hid %in% xu$hid]   # same HOGs in both frames
  for (s in names(SPEC)) { r <- fitm(xp, SPEC[[s]], "is_real"); if (!is.null(r)) MOD[[length(MOD) + 1]] <- r[, `:=`(route = rt, ind = ind, frame = "pruned", model = s)] }
  r <- fitm(xu, SPEC$R0, "is_real"); if (!is.null(r)) MOD[[length(MOD) + 1]] <- r[, `:=`(route = rt, ind = ind, frame = "unpruned", model = "R0")]
}
MOD <- rbindlist(MOD); fwrite(MOD, file.path(O, "models.tsv"), sep = "\t")
print(MOD[, .(route, ind, frame, model, est = round(est, 3), se_cl = round(se_cl, 3), p_cl = signif(p_cl, 3), n_hog, n_clade, n_events,
              rate_real = round(rate_real, 3), rate_null = round(rate_null, 3))])
V <- character(0); surv <- c()
for (rt in c("A", "M")) {
  a <- MOD[route == rt & ind == "relax_relaxed" & frame == "pruned" & model == "R1"]
  u <- MOD[route == rt & ind == "relax_relaxed" & frame == "unpruned" & model == "R0"]
  p0 <- MOD[route == rt & ind == "relax_relaxed" & frame == "pruned" & model == "R0"]
  info <- nrow(a) == 1 && a$n_clade >= 2 && a$n_hog >= 50
  sv <- info && a$est > 0 && a$p_cl < 0.05; surv[rt] <- if (info) sv else NA
  V <- c(V, sprintf("D3\t%s\tinformative=%s\tR1_logOR=%s\tR1_p=%s\tn_hog=%s\tn_clade=%s\tsurvives=%s\tratio_prunedR0_unprunedR0=%s",
    rt, info, if (nrow(a)) round(a$est, 3) else NA, if (nrow(a)) signif(a$p_cl, 3) else NA, if (nrow(a)) a$n_hog else 0,
    if (nrow(a)) a$n_clade else 0, if (info) sv else NA, if (nrow(p0) && nrow(u)) round(p0$est / u$est, 3) else NA))
}
s <- surv[!is.na(surv)]
ov <- if (length(s) < 2) { paste0("only_", paste(names(s), collapse = ""), "_informative:", if (length(s)) ifelse(s, "survives", "fails") else "none")
      } else if (all(s)) "route_robust" else if (!any(s)) "route_dependent" else paste0("single_route_", names(s)[s])
V <- c(V, paste0("D3\toverall\t", ov)); writeLines(V, file.path(O, "verdict.tsv")); cat(V, sep = "\n")
