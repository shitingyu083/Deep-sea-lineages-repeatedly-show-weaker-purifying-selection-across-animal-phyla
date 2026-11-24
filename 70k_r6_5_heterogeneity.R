#!/usr/bin/env Rscript
## 70k_r6_5_heterogeneity.R
## ADDENDUM_R6_5 (work/addon/R6/ADDENDUM_R6_5.md, md5 0973c9e50ed97b27583df99aa456c421,
##               frozen 2026-11-02T15:43:11Z, before any set-B statistic)
##
## Question: is the deep-branch relaxation gene-specific (heterogeneous across genes,
## possibly functionally structured), or a near-uniform shift of the whole gene distribution?
## Internal noise reference: the exchangeable-label (null1) side of the SAME genes -- same
## orthogroups, same gene trees, same RELAX settings, near-identical tip-count distribution,
## no deep/shallow biology.
##
## Unit = (clade, hog). Groups are NOT paired across sides (deep and pseudo-lineage tip sets
## are disjoint by construction), so no group-level real-minus-null difference is formed here.

suppressPackageStartupMessages(library(data.table))
setDTthreads(min(32L, parallel::detectCores()))

R     <- path.expand("~/deepsea_conv")
R6    <- file.path(R, "work/addon/R6")
OUT   <- file.path(R6, "R6_5"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SEED  <- 20261102
KLO   <- 0.01; KHI <- 49          # PREREG_r6 boundary range (primary)
NB    <- 2000                     # bootstrap / synthetic-null replicates (pre-declared)
NPERM <- 10000                    # GO permutation replicates (pre-declared)

H <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
P <- function(...) cat(sprintf(...), "\n", sep = "")

## ---------------------------------------------------------------- load
read_set <- function(sub) {
  f <- file.path(R6, sub, "r6_2_fits.tsv")
  if (!file.exists(f)) stop("missing ", f)
  d <- fread(f, sep = "\t", data.table = TRUE)
  need <- c("clade", "hog", "set", "group", "n_tips", "K", "fit")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("missing columns in ", f, ": ", paste(miss, collapse = ", "))
  d[, `:=`(clade = as.character(clade), hog = as.character(hog),
           set = as.character(set), group = as.character(group),
           K = suppressWarnings(as.numeric(K)),
           n_tips = suppressWarnings(as.integer(n_tips)))]
  d
}

prepare <- function(d, tag, lo = KLO, hi = KHI) {
  n0 <- nrow(d)
  d <- d[fit == "ok"];            n1 <- nrow(d)
  d <- d[is.finite(K) & K > 0];   n2 <- nrow(d)
  d[, l := log10(K)]
  nb <- d[l < log10(lo) | l > log10(hi), .N]
  P("[%s] rows %d -> fit ok %d -> K finite >0 %d ; dropped at K bound %d (%.1f %%)",
    tag, n0, n1, n2, nb, 100 * nb / max(n2, 1))
  d[l >= log10(lo) & l <= log10(hi)]
}

pergene <- function(d) {
  d[, .(m = mean(l), v = if (.N >= 2) var(l) else NA_real_, n = .N, mt = mean(n_tips)),
    by = .(clade, hog, set)]
}

wide <- function(pg) {
  w <- dcast(pg, clade + hog ~ set, value.var = c("m", "v", "n", "mt"))
  for (cc in c("m_real", "m_null1", "v_real", "v_null1", "n_real", "n_null1",
               "mt_real", "mt_null1"))
    if (!cc %in% names(w)) w[, (cc) := NA_real_]
  w[, D := m_real - m_null1]
  w[]
}

boot_ratio <- function(xr, xn, B = NB) {
  ok <- is.finite(xr) & is.finite(xn)
  xr <- xr[ok]; xn <- xn[ok]
  n  <- length(xr)
  r_obs <- var(xr) / var(xn)
  if (n < 5) return(c(r_obs, NA_real_, NA_real_, n))
  r <- vapply(seq_len(B), function(b) {
    i <- sample.int(n, n, TRUE)
    a <- var(xr[i]); bb <- var(xn[i])
    if (!is.finite(a) || !is.finite(bb) || bb == 0) NA_real_ else a / bb
  }, numeric(1))
  r <- r[is.finite(r)]
  c(r_obs, quantile(r, .025, names = FALSE), quantile(r, .975, names = FALSE), n)
}

## ------------------------------------------------- synthetic deep side
## Residual = log10 K minus that gene's own null-side mean, so the pool is centred and carries
## only within-gene pseudo-lineage noise at matched tip depth. Pooling per clade x depth stratum.
resid_pool <- function(d, minnull = 2L) {
  d <- copy(d); d[, gid := paste(clade, hog, sep = "|")]
  nm <- d[set == "null1", .(m = mean(l), n = .N), by = .(clade, gid)][n >= minnull]
  res <- merge(d[set == "null1", .(clade, gid, n_tips, l)], nm[, .(clade, gid, m)],
               by = c("clade", "gid"))
  res[, `:=`(r = l - m, bucket = pmin(n_tips, 4L))]
  pdt <- res[, .(r = list(r)), by = .(clade, bucket)]
  cd  <- res[, .(r = list(r)), by = clade]
  list(n = nrow(res), genes = uniqueN(res$gid),
       poolv = setNames(pdt$r, paste(pdt$clade, pdt$bucket, sep = "|")),
       clw   = setNames(cd$r, cd$clade))
}

synth_perm <- function(w, d, B = NB, minnull = 2L) {
  d <- copy(d)
  d[, `:=`(gid = paste(clade, hog, sep = "|"), bucket = pmin(n_tips, 4L))]
  w <- copy(w); w[, gid := paste(clade, hog, sep = "|")]
  rp <- resid_pool(d, minnull)
  mnull <- setNames(w$m_null1, w$gid)

  ## group counts at each tip-depth stratum, per gene -- exactly as the true deep side has
  req <- d[set == "real" & gid %in% w$gid, .(n = .N), by = .(gid, clade, bucket)]
  cb  <- unique(req[, .(clade, bucket)])
  cb[, k := paste(clade, bucket, sep = "|")]
  rq  <- split(req, paste(req$clade, req$bucket, sep = "|"))

  Vm <- numeric(B); Vi <- numeric(B)
  for (b in seq_len(B)) {
    parts <- vector("list", 0)
    for (i in seq_len(nrow(cb))) {
      s  <- rq[[cb$k[i]]]
      if (is.null(s) || !nrow(s)) next
      pv <- rp$poolv[[cb$k[i]]]
      if (is.null(pv) || length(pv) < 30) pv <- rp$clw[[cb$clade[i]]]
      if (is.null(pv) || !length(pv)) next
      rid <- rep(seq_len(nrow(s)), s$n)
      parts[[length(parts) + 1]] <-
        data.table(gid = s$gid[rid], v = pv[sample.int(length(pv), length(rid), TRUE)])
    }
    dd  <- rbindlist(parts)
    ag  <- dd[, .(mu = mean(v), vv = if (.N >= 2) var(v) else NA_real_), by = gid]
    syn <- mnull[ag$gid] + ag$mu
    Vm[b] <- var(syn[is.finite(syn)])
    Vi[b] <- mean(ag$vv, na.rm = TRUE)
  }
  list(Vm = Vm, Vi = Vi, pool_n = rp$n, pool_genes = rp$genes)
}

## ------------------------------------------------------------- analysis
analyse <- function(tag, sub, nboot = NB) {
  H(sprintf("SET %s   (%s)", tag, sub))
  d  <- prepare(read_set(sub), tag)
  pg <- pergene(d); w <- wide(pg)

  both <- w[!is.na(n_real) & !is.na(n_null1) & n_real >= 1 & n_null1 >= 1]
  P("genes: total %d ; >=1 valid group on BOTH sides %d ; >=2 on both %d",
    nrow(w), nrow(both),
    nrow(w[!is.na(n_real) & !is.na(n_null1) & n_real >= 2 & n_null1 >= 2]))
  P("clades %d ; valid groups real %d null1 %d", uniqueN(w$clade),
    sum(w$n_real, na.rm = TRUE), sum(w$n_null1, na.rm = TRUE))

  ## ---- E1 gene-mean dispersion
  H("E1  var of per-gene mean log10 K, real vs null1 (same genes, both sides)")
  b1 <- boot_ratio(both$m_real, both$m_null1)
  P("n genes %d | V_mean real %.5f  null1 %.5f | R = %.3f  95%%CI [%.3f, %.3f]",
    b1[4], var(both$m_real), var(both$m_null1), b1[1], b1[2], b1[3])
  P("mean  m: real %.4f null1 %.4f | median m: real %.4f null1 %.4f",
    mean(both$m_real), mean(both$m_null1), median(both$m_real), median(both$m_null1))

  eq <- both[n_real == n_null1]
  eb <- boot_ratio(eq$m_real, eq$m_null1)
  P("E1b equal valid group counts: n %d | R %.3f 95%%CI [%.3f, %.3f]", eb[4], eb[1], eb[2], eb[3])

  tc <- both[!is.na(mt_real) & !is.na(mt_null1) & abs(mt_real - mt_null1) <= 1]
  tb <- boot_ratio(tc$m_real, tc$m_null1)
  P("E1c |mean n_tips diff| <= 1: n %d | R %.3f 95%%CI [%.3f, %.3f]", tb[4], tb[1], tb[2], tb[3])

  ## ---- E2 within-gene dispersion
  H("E2  within-gene dispersion of log10 K (>=2 valid groups on that side)")
  wr <- both[!is.na(v_real)]
  wn <- both[!is.na(v_null1)]
  vb <- boot_ratio(both$v_real, both$v_null1)
  P("mean   V_in  real %.5f (n=%d)  null1 %.5f (n=%d)  ratio %.3f",
    mean(wr$v_real), nrow(wr), mean(wn$v_null1), nrow(wn), mean(wr$v_real) / mean(wn$v_null1))
  P("median V_in  real %.5f  null1 %.5f  ratio %.3f",
    median(wr$v_real), median(wn$v_null1), median(wr$v_real) / median(wn$v_null1))
  P("implied sd(log10 K): real %.3f  null1 %.3f",
    sqrt(mean(wr$v_real)), sqrt(mean(wn$v_null1)))
  P("paired bootstrap, genes with v on both sides: R_in %.3f 95%%CI [%.3f, %.3f] n %d",
    vb[1], vb[2], vb[3], vb[4])

  ## ---- E3 shift
  H("E3  pooled shift")
  P("mean D %.4f  median D %.4f  share D<0 %.3f  n %d",
    mean(both$D), median(both$D), mean(both$D < 0), nrow(both))
  cen <- copy(both)
  cen[, `:=`(mrc = m_real - mean(m_real), mnc = m_null1 - mean(m_null1)), by = clade]
  cen[, Dc := mrc - mnc]
  P("after centring each side on its own clade mean: mean %.4f median %.4f",
    mean(cen$Dc), median(cen$Dc))
  pc <- cen[, .(md = median(D), n = .N), by = clade][order(md)]
  P("per-clade median D: %s", paste(sprintf("%s %.3f(n=%d)", pc$clade, pc$md, pc$n), collapse = "  "))

  ## ---- E4 / E5 self-calibrated synthetic null
  H("E4  synthetic deep side resampled from the pseudo-lineage residual pool (PRIMARY)")
  set.seed(SEED)
  sn  <- synth_perm(both, d, B = nboot, minnull = 2L)
  obs <- var(both$m_real)
  pE4 <- (1 + sum(sn$Vm >= obs)) / (1 + length(sn$Vm))
  zE4 <- (obs - mean(sn$Vm)) / sd(sn$Vm)
  P("residual pool: %d residuals from %d genes with >=2 null groups", sn$pool_n, sn$pool_genes)
  P("observed V_mean(real) %.5f | synthetic mean %.5f sd %.5f | z %.3f | upper-tail p %.4f",
    obs, mean(sn$Vm), sd(sn$Vm), zE4, pE4)
  P("synthetic null quantiles: 50%% %.5f  95%% %.5f  99%% %.5f",
    quantile(sn$Vm, .5), quantile(sn$Vm, .95), quantile(sn$Vm, .99))

  obsi <- mean(both$v_real, na.rm = TRUE)
  pE5  <- (1 + sum(sn$Vi >= obsi)) / (1 + length(sn$Vi))
  P("E5 observed V_in(real) %.5f | synthetic mean %.5f | upper-tail p %.4f",
    obsi, mean(sn$Vi), pE5)

  set.seed(SEED)
  sn1  <- synth_perm(both, d, B = nboot, minnull = 1L)
  pE4s <- (1 + sum(sn1$Vm >= obs)) / (1 + length(sn1$Vm))
  P("SENSITIVITY pool>=1 null group: synthetic mean %.5f sd %.5f | p %.4f",
    mean(sn1$Vm), sd(sn1$Vm), pE4s)

  list(tag = tag, sub = sub, d = d, w = w, both = both,
       R1 = b1[1], R1ci = b1[2:3], Rin = vb[1], Rin_ci = vb[2:3],
       pE4 = pE4, zE4 = zE4, pE5 = pE5, pE4s = pE4s, snVm = sn$Vm)
}

## ------------------------------------------------------------------ run
res <- list()
res$A <- analyse("A", "R6_2")
res$B <- analyse("B", "R6_2b")

## --------------------------------------------------------------- E6 GO
H("E6  functional structure: between-term weighted variance of D, permutation within clade (SET B)")
load_go <- function(clades) {
  go <- rbindlist(lapply(clades, function(cc) {
    for (cand in unique(c(cc, sub("_m$", "", cc)))) {
      p <- file.path(R, "work/func_v2", cand, "hog_annot.tsv")
      if (file.exists(p)) {
        a <- fread(p, sep = "\t", select = c("hog", "go"), data.table = TRUE)
        return(data.table(clade = cc, hog = as.character(a$hog), go = as.character(a$go)))
      }
    }
    P("  [GO] no annotation file for %s", cc); NULL
  }), fill = TRUE)
  go <- go[!is.na(go) & nzchar(go)]
  go[, .(go = unlist(strsplit(go, ","))), by = .(clade, hog)][nzchar(go)]
}

go_test <- function(w, nperm = NPERM, minN = 15) {
  g <- load_go(unique(w$clade))
  g <- merge(g, w[, .(clade, hog, D)], by = c("clade", "hog"))
  P("annotated gene-terms in analysis set: %d (genes %d)", nrow(g),
    uniqueN(paste(g$clade, g$hog)))
  set.seed(SEED + 1)
  W <- 0; Wn <- numeric(nperm); terms_used <- 0L; genes_used <- 0L
  for (cc in unique(w$clade)) {
    s <- g[clade == cc]
    if (!nrow(s)) next
    X <- as.matrix(table(s$hog, s$go)); X[X > 0] <- 1
    keep <- colSums(X) >= minN
    if (!sum(keep)) next
    X  <- X[, keep, drop = FALSE]
    wc <- w[clade == cc]
    Dv <- wc$D[match(rownames(X), wc$hog)]
    if (all(is.na(Dv))) next
    Dv[is.na(Dv)] <- 0
    nn <- colSums(X); terms_used <- terms_used + ncol(X); genes_used <- genes_used + nrow(X)
    gr <- mean(Dv)
    W  <- W + sum(nn * (colSums(X * Dv) / nn - gr)^2)
    for (b0 in seq(1, nperm, by = 1000)) {
      cols <- b0:min(b0 + 999, nperm)
      DD <- vapply(cols, function(i) sample(Dv, length(Dv), FALSE), numeric(length(Dv)))
      S  <- crossprod(X, DD)
      Wn[cols] <- Wn[cols] + colSums(S^2 / nn) - 2 * gr * colSums(S) + gr^2 * sum(nn)
    }
  }
  P("genes entering %d ; terms with >= %d genes: %d ; global W %.1f",
    genes_used, minN, terms_used, W)
  p <- (1 + sum(Wn >= W)) / (1 + nperm)
  P("E6 global permutation p = %.4f  (nperm %d)", p, nperm)
  p
}

pE6 <- go_test(res$B$both)

## ------------------------------------------------------------ E7 sets
H("E7  a-priori functional sets on D (SET B)")
pE7 <- NA_real_
f7  <- file.path(R6, "R6_3", "r6_3_fits.tsv")
if (file.exists(f7)) {
  s <- fread(f7, sep = "\t", data.table = TRUE)
  s[, `:=`(clade = as.character(clade), hog = as.character(hog))]
  s <- s[, .(LIGHT = max(LIGHT), CIRCADIAN = max(CIRCADIAN), PIGMENT = max(PIGMENT)),
         by = .(clade, hog)]
  m <- merge(res$B$both, s, by = c("clade", "hog"))
  P("genes of set B covered by the a-priori set table: %d", nrow(m))
  if (nrow(m) >= 20) {
    pv <- c()
    for (st in c("LIGHT", "CIRCADIAN", "PIGMENT")) {
      idx <- m[[st]] == 1
      if (sum(idx) < 5) { P("%-10s n %d -- too few to test", st, sum(idx)); pv <- c(pv, NA); next }
      tt <- suppressWarnings(wilcox.test(m$D[idx], m$D[!idx]))
      P("%-10s n %d  median D %.3f vs %.3f  Wilcoxon p %.3f",
        st, sum(idx), median(m$D[idx]), median(m$D[!idx]), tt$p.value)
      pv <- c(pv, tt$p.value)
    }
    pE7 <- suppressWarnings(min(p.adjust(pv, "holm"), na.rm = TRUE))
    P("Holm over 3 sets: min adjusted p = %.3f", pE7)
  }
} else P("a-priori set table not found: %s", f7)

## --------------------------------------------------------- E8 / power
H("E8  set A alongside (already inspected before freezing -- disclosure in the addendum)")
P("A: R_mean %.3f [%.3f, %.3f]  E4 p %.4f  |  B: R_mean %.3f [%.3f, %.3f]  E4 p %.4f",
  res$A$R1, res$A$R1ci[1], res$A$R1ci[2], res$A$pE4,
  res$B$R1, res$B$R1ci[1], res$B$R1ci[2], res$B$pE4)
P("A: V_in ratio %.3f  E5 p %.4f  |  B: V_in ratio %.3f  E5 p %.4f",
  res$A$Rin, res$A$pE5, res$B$Rin, res$B$pE5)

H("POWER  minimum detectable gene-level variance (set B structure)")
thr    <- quantile(res$B$snVm, .95, names = FALSE)
grid   <- c(0, 0.0025, 0.005, 0.0075, 0.010, 0.015, 0.020, 0.030)
base_v <- var(res$B$both$m_null1)
pw <- vapply(grid, function(s2) {
  set.seed(SEED + 7)
  hit <- 0L
  for (r in 1:200) {
    b   <- rnorm(length(res$B$both$m_null1), 0, sqrt(s2))
    syn <- res$B$both$m_null1 + b + rnorm(length(b), 0, sqrt(mean(res$B$snVm)))
    if (var(syn) > thr) hit <- hit + 1L
  }
  hit / 200
}, numeric(1))
for (i in seq_along(grid))
  P("  sigma2_bio %.4f (%.1f %% of observed null V_mean)  power %.2f",
    grid[i], 100 * grid[i] / base_v, pw[i])
P("observed null V_mean(set B) %.5f ; E4 95%% threshold %.5f", base_v, thr)

## ----------------------------------------------------------- decision
H("DECISION  (pre-declared, judged on SET B only)")
C1 <- (res$B$R1ci[1] > 1) || (res$B$pE4 < 0.05)
C2 <- is.finite(pE6) && pE6 < 0.05
P("C1 gene-specific variance beyond pseudo-lineage noise : %s", ifelse(C1, "HOLDS", "FAILS"))
P("     R_mean %.3f CI [%.3f, %.3f] ; E4 p %.4f ; z %.2f",
  res$B$R1, res$B$R1ci[1], res$B$R1ci[2], res$B$pE4, res$B$zE4)
P("C2 functional structure of D (GO)                    : %s", ifelse(C2, "HOLDS", "FAILS"))
P("     E6 global p %.4f ; a-priori sets Holm min p %s",
  pE6, ifelse(is.na(pE7), "NA", sprintf("%.3f", pE7)))
verdict <- if (C1 && C2) "REORGANIZATION SUPPORTED" else
           if (!C1) "NOISE-LIMITED HETEROGENEITY" else "UNSTRUCTURED HETEROGENEITY"
P("VERDICT: %s", verdict)

## -------------------------------------------------------------- write
fwrite(res$A$both, file.path(OUT, "pergene_A.tsv"), sep = "\t")
fwrite(res$B$both, file.path(OUT, "pergene_B.tsv"), sep = "\t")
sink(file.path(OUT, "r6_5_summary.txt"))
cat("verdict\t", verdict, "\n", sep = "")
cat("B n_genes\t", nrow(res$B$both), "\n", sep = "")
cat("B R_mean\t", res$B$R1, "\tlo\t", res$B$R1ci[1], "\thi\t", res$B$R1ci[2], "\n", sep = "")
cat("B E4 p\t", res$B$pE4, "\tz\t", res$B$zE4, "\n", sep = "")
cat("B E4 p sensitivity pool>=1\t", res$B$pE4s, "\n", sep = "")
cat("B E5 p\t", res$B$pE5, "\n", sep = "")
cat("B E6 p\t", pE6, "\n", sep = "")
cat("B E7 holm\t", pE7, "\n", sep = "")
cat("A n_genes\t", nrow(res$A$both), "\n", sep = "")
cat("A R_mean\t", res$A$R1, "\tlo\t", res$A$R1ci[1], "\thi\t", res$A$R1ci[2], "\n", sep = "")
cat("A E4 p\t", res$A$pE4, "\n", sep = "")
sink()
cat("\nWROTE", OUT, "\n")
