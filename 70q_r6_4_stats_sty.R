# R6-4 dose-response / Ne-proxy regression (ADDENDUM_R6_4.md part 2).
# Usage: Rscript 70q_r6_4_stats_sty.R <covariates.tsv> <OUT dir> <tag>=<stats dir> [<tag>=<stats dir> ...]
# Each stats dir holds r6_2_lineages.tsv + boot_lineage_eff.tsv written by 70o. Tags: A, B, AB.
suppressPackageStartupMessages(library(parallel))
a <- commandArgs(TRUE); COV <- a[1]; OUT <- a[2]; SETS <- a[-(1:2)]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NC <- 16; NPERM <- 10000; FOCAL <- c("depth", "age", "enc_ratio")
cv <- read.delim(COV, stringsAsFactors = FALSE)
cv$key <- paste(cv$clade, cv$canon, sep = ":")
z <- function(x) (x - mean(x)) / sd(x)
log <- file(file.path(OUT, "r6_4_report.txt"), "w")
say <- function(...) { s <- sprintf(...); cat(s, "\n"); writeLines(s, log) }
res <- list()

fit <- function(d, form, w = NULL) {
  m <- if (is.null(w)) lm(form, d) else lm(form, d, weights = w)
  cf <- summary(m)$coefficients; cf
}
tstat <- function(d, form, v, w = NULL) { cf <- fit(d, form, w); if (v %in% rownames(cf)) cf[v, 3] else NA }
perm_p <- function(d, form, v, seed, w = NULL) {
  t0 <- tstat(d, form, v, w); idx <- split(seq_len(nrow(d)), d$clade)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  tp <- unlist(mclapply(seq_len(NPERM), function(i) {
    dd <- d; for (g in idx) dd[[v]][g] <- if (length(g) > 1) sample(d[[v]][g]) else d[[v]][g]
    tstat(dd, form, v, w) }, mc.cores = NC, mc.set.seed = TRUE))
  c(t = t0, p = (1 + sum(abs(tp) >= abs(t0), na.rm = TRUE)) / (1 + sum(!is.na(tp))))
}
prep <- function(dir, extra = character(0)) {
  L <- read.delim(file.path(dir, "r6_2_lineages.tsv"), stringsAsFactors = FALSE)
  L$key <- paste(L$clade, L$canon, sep = ":")
  d <- merge(L[, c("key", "eff", "n_hog")], cv, by = "key")
  need <- c("eff", FOCAL, "chemo", extra)
  ok <- complete.cases(d[, need]); d <- d[ok, ]
  for (v in c(FOCAL, extra)) d[[v]] <- z(d[[v]])
  d$clade <- factor(d$clade); attr(d, "ndrop") <- sum(!ok); attr(d, "nlin") <- nrow(L); d
}
boot_ci <- function(dir, d, form, vars) {
  B <- read.delim(file.path(dir, "boot_lineage_eff.tsv"), check.names = FALSE)
  B <- as.matrix(B[, d$key, drop = FALSE])
  bb <- t(apply(B, 1, function(e) { dd <- d; dd$eff <- e; cf <- coef(lm(form, dd)); cf[vars] }))
  apply(bb, 2, quantile, c(.025, .975), na.rm = TRUE)
}
PRIM <- eff ~ clade + depth + age + enc_ratio + chemo

for (s in SETS) {
  tag <- sub("=.*", "", s); dir <- sub("^[^=]*=", "", s)
  d <- prep(dir)
  say("== set %s (%s): lineages %d, used %d, dropped (NA covariate) %d, clades %d",
      tag, dir, attr(d, "nlin"), nrow(d), attr(d, "ndrop"), nlevels(d$clade))
  cf <- fit(d, PRIM); ci <- boot_ci(dir, d, PRIM, FOCAL)
  pp <- t(sapply(seq_along(FOCAL), function(i) perm_p(d, PRIM, FOCAL[i], 20261101 + i)))
  holm <- p.adjust(pp[, "p"], "holm")
  for (i in seq_along(FOCAL)) {
    v <- FOCAL[i]
    say("PRIMARY %s %-9s beta %+.4f boot95 [%+.4f, %+.4f] t %+.3f perm_p %.4f holm %.4f",
        tag, v, cf[v, 1], ci[1, v], ci[2, v], pp[i, "t"], pp[i, "p"], holm[i])
    res[[length(res) + 1]] <- data.frame(set = tag, model = "primary", term = v, beta = cf[v, 1],
      lo = ci[1, v], hi = ci[2, v], t = pp[i, "t"], perm_p = pp[i, "p"], holm = holm[i], n = nrow(d))
  }
  if ("chemo" %in% rownames(cf)) say("PRIMARY %s chemo beta %+.4f (OLS p %.4f; not tested)", tag, cf["chemo", 1], cf["chemo", 4])
  # secondary (descriptive; OLS p only)
  sec <- list(weighted = list(PRIM, "n_hog"),
              gsize = list(update(PRIM, . ~ . + log_gsize), NULL, "log_gsize"),
              T_hab = list(update(PRIM, . ~ . + T_hab), NULL, "T_hab"),
              only_depth = list(eff ~ clade + depth, NULL), only_age = list(eff ~ clade + age, NULL),
              only_enc = list(eff ~ clade + enc_ratio, NULL))
  for (nm in names(sec)) {
    sp <- sec[[nm]]; ex <- if (length(sp) > 2) sp[[3]] else character(0)
    ds <- prep(dir, ex); w <- if (!is.null(sp[[2]])) ds[[sp[[2]]]] else NULL
    cf2 <- fit(ds, sp[[1]], w)
    for (v in intersect(c(FOCAL, ex), rownames(cf2))) {
      say("SECONDARY %s %-10s %-9s beta %+.4f OLS_p %.4f n %d", tag, nm, v, cf2[v, 1], cf2[v, 4], nrow(ds))
      res[[length(res) + 1]] <- data.frame(set = tag, model = nm, term = v, beta = cf2[v, 1], lo = NA, hi = NA,
        t = cf2[v, 3], perm_p = NA, holm = NA, n = nrow(ds))
    }
  }
}
R <- do.call(rbind, res)
write.table(R, file.path(OUT, "r6_4_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
P <- subset(R, model == "primary")
g <- function(st, v, col) { x <- P[P$set == st & P$term == v, col]; if (length(x)) x else NA }
if (all(c("A", "B") %in% P$set)) {
  rule <- function(v) isTRUE(g("B", v, "holm") < .05 && sign(g("B", v, "beta")) == sign(g("A", v, "beta")))
  dr <- rule("depth") || rule("age"); ne <- rule("enc_ratio")
  say("RULE dose-response (depth or age: B holm<.05 & same sign in A): %s [depth %s, age %s]",
      if (dr) "SUPPORTED" else "NOT SUPPORTED", rule("depth"), rule("age"))
  say("RULE Ne-proxy (enc_ratio: B holm<.05 & same sign in A): %s -> wording: %s",
      if (ne) "SUPPORTED" else "NOT SUPPORTED",
      if (ne) "consistent with reduced selection efficiency" else "mechanism undetermined")
  for (v in FOCAL) say("SIGN %s: expected %s, A %+.4f, B %+.4f", v,
                       c(depth = "<0", age = "none", enc_ratio = ">0")[v], g("A", v, "beta"), g("B", v, "beta"))
}
close(log); quit(status = 0)
