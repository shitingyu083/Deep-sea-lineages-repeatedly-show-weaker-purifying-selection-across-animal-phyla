#!/usr/bin/env Rscript
## 70k2_r6_5_power.R  -- corrected pre-declared power computation for ADDENDUM_R6_5.
## The first run of 70k_r6_5_heterogeneity.R added N(0, sqrt(mean(Vm))) at gene level, where Vm is
## the variance of the SYNTHETIC GENE MEANS (0.099); that sd (0.314) is the between-gene spread of the
## finished synthetic means, not the sampling noise of one gene's synthetic mean. Power therefore came
## out 1.00 even at sigma2_bio = 0 and the table is discarded. Corrected here: the added biological
## component is drawn at the gene level with variance sigma2_bio, on top of the ACTUAL synthetic
## construction (residual resampling at the observed group counts and tip strata). Logged in
## DEVIATIONS_r6.md; it does not touch any E1-E8 estimand or the verdict.

suppressPackageStartupMessages(library(data.table))
setDTthreads(min(32L, parallel::detectCores()))
R <- path.expand("~/deepsea_conv"); R6 <- file.path(R, "work/addon/R6")
OUT <- file.path(R6, "R6_5"); SEED <- 20261102
KLO <- 0.01; KHI <- 49; NB <- 2000
P <- function(...) cat(sprintf(...), "\n", sep = "")
H <- function(x) cat("\n", strrep("=",78), "\n", x, "\n", strrep("=",78), "\n", sep = "")

d <- fread(file.path(R6, "R6_2b", "r6_2_fits.tsv"), sep = "\t", data.table = TRUE)
d[, `:=`(clade = as.character(clade), hog = as.character(hog), set = as.character(set),
         K = suppressWarnings(as.numeric(K)), n_tips = suppressWarnings(as.integer(n_tips)))]
d <- d[fit == "ok" & is.finite(K) & K > 0]
d[, l := log10(K)]
d <- d[l >= log10(KLO) & l <= log10(KHI)]
d[, `:=`(gid = paste(clade, hog, sep = "|"), bucket = pmin(n_tips, 4L))]

pg <- d[, .(m = mean(l), n = .N), by = .(clade, hog, set)]
w  <- dcast(pg, clade + hog ~ set, value.var = c("m", "n"))
w  <- w[!is.na(n_real) & !is.na(n_null1) & n_real >= 1 & n_null1 >= 1]
w[, gid := paste(clade, hog, sep = "|")]
mnull <- setNames(w$m_null1, w$gid)

nm  <- d[set == "null1", .(m = mean(l), n = .N), by = .(clade, gid)][n >= 2]
res <- merge(d[set == "null1", .(clade, gid, bucket, l)], nm[, .(clade, gid, m)], by = c("clade","gid"))
res[, r := l - m]
pdt <- res[, .(r = list(r)), by = .(clade, bucket)]
cd  <- res[, .(r = list(r)), by = clade]
poolv <- setNames(pdt$r, paste(pdt$clade, pdt$bucket, sep = "|"))
clw   <- setNames(cd$r,  cd$clade)

req <- d[set == "real" & gid %in% w$gid, .(n = .N), by = .(gid, clade, bucket)]
cb  <- unique(req[, .(clade, bucket)]); cb[, k := paste(clade, bucket, sep = "|")]
rq  <- split(req, paste(req$clade, req$bucket, sep = "|"))

## one draw of the synthetic gene means, exactly as in E4
draw_mu <- function() {
  parts <- vector("list", 0)
  for (i in seq_len(nrow(cb))) {
    s <- rq[[cb$k[i]]]; if (is.null(s) || !nrow(s)) next
    pv <- poolv[[cb$k[i]]]; if (is.null(pv) || length(pv) < 30) pv <- clw[[cb$clade[i]]]
    rid <- rep(seq_len(nrow(s)), s$n)
    parts[[length(parts) + 1]] <-
      data.table(gid = s$gid[rid], v = pv[sample.int(length(pv), length(rid), TRUE)])
  }
  rbindlist(parts)[, .(mu = mean(v)), by = gid]
}

set.seed(SEED)
nullVm <- numeric(NB); nullmu <- vector("list", NB)
for (b in seq_len(NB)) {
  ag <- draw_mu(); nullmu[[b]] <- setNames(ag$mu, ag$gid)
  nullVm[b] <- var(mnull[ag$gid] + ag$mu)
}
thr <- quantile(nullVm, .95, names = FALSE)
P("E4 null distribution: mean %.5f sd %.5f ; 95%% threshold %.5f ; observed %.5f",
  mean(nullVm), sd(nullVm), thr, var(w$m_real))

## sd of one gene's synthetic mean (the piece the first run mis-scaled)
P("sd of a gene's synthetic mean (mu): %.4f ; between-gene sd of null means: %.4f",
  mean(vapply(nullmu[1:200], sd, numeric(1))), sd(w$m_null1))

grid <- c(0, 0.0005, 0.001, 0.002, 0.005, 0.010, 0.020, 0.040)
base_v <- var(w$m_null1)
pw <- vapply(grid, function(s2) {
  set.seed(SEED + 7)
  hit <- 0L
  for (r in 1:400) {
    mu <- nullmu[[((r - 1) %% NB) + 1]]
    syn <- mnull[names(mu)] + mu + rnorm(length(mu), 0, sqrt(s2))
    if (var(syn) > thr) hit <- hit + 1L
  }
  hit / 400
}, numeric(1))
H("POWER  corrected (set B structure)")
for (i in seq_along(grid))
  P("  sigma2_bio %.4f (%5.1f %% of null V_mean)  sd_bio %.3f  power %.2f",
    grid[i], 100 * grid[i] / base_v, sqrt(grid[i]), pw[i])
md <- grid[which(pw >= 0.80)[1]]
P("minimum detectable sigma2_bio at power >= 0.80: %s (sd_bio %s)",
  ifelse(is.na(md), "> 0.040", sprintf("%.4f", md)),
  ifelse(is.na(md), "> 0.20", sprintf("%.3f", sqrt(md))))
P("interpretation: the deep side would have to differ from what pseudo-lineages produce by this much")
P("  gene-level sd (in log10 K units) before E4 would reject at alpha 0.05 with power 0.80.")

writeLines(c(
  sprintf("sigma2_bio\tpower\tsd_bio\tpct_of_null_Vmean"),
  sprintf("%.4f\t%.2f\t%.3f\t%.1f", grid, pw, sqrt(grid), 100*grid/base_v)),
  file.path(OUT, "power_B_corrected.tsv"))
sink(file.path(OUT, "power_B_corrected.txt"))
cat("E4 null mean\t", mean(nullVm), "\n", sep = "")
cat("E4 null sd\t", sd(nullVm), "\n", sep = "")
cat("E4 95pct threshold\t", thr, "\n", sep = "")
cat("observed V_mean real\t", var(w$m_real), "\n", sep = "")
cat("min detectable sigma2_bio (power .80)\t", md, "\n", sep = "")
sink()
cat("\nDONE\n")
