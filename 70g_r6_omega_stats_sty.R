# R6-3 statistics (PREREG_r6.md R6-3a / R6-3b) on 70f output. Usage: Rscript 70g_r6_omega_stats_sty.R DIR
args <- commandArgs(TRUE); D <- if (length(args)) args[1] else "~/deepsea_conv/work/addon/R6/R6_3"
x <- read.delim(file.path(D, "r6_3_fits.tsv"), na.strings = "NA", stringsAsFactors = FALSE)
NB <- 2000; set.seed(20261031)
x$relaxed <- as.integer(!is.na(x$relax_p) & !is.na(x$relax_K) & x$relax_p < 0.05 & x$relax_K < 1)
x$Dpur <- log10(x$wpur_T + 0.001) - log10(x$wpur_B + 0.001)
x$Dpos <- x$ppos_T - x$ppos_B
CL <- unique(x$clade); out <- list()
put <- function(...) out[[length(out) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
ci <- function(b) quantile(b, c(.025, .975), na.rm = TRUE, names = FALSE)
bp <- function(b) min(1, 2 * min(mean(b <= 0, na.rm = TRUE), mean(b >= 0, na.rm = TRUE)))

## R6-3a paired real - null1 (HOGs with both BUSTED-PH fits)
ok <- x[x$bph_ok == 1, ]
pr <- merge(ok[ok$label == "real", ], ok[ok$label == "null1", ], by = c("clade", "hog"), suffixes = c(".r", ".n"))
for (v in c("Dpur", "Dpos")) {
  d <- pr[[paste0(v, ".r")]] - pr[[paste0(v, ".n")]]; keep <- !is.na(d)
  dd <- split(d[keep], pr$clade[keep])
  cm <- sapply(dd, mean)
  est <- mean(cm)
  b <- replicate(NB, mean(sapply(dd, function(z) mean(sample(z, replace = TRUE)))))
  put(analysis = "R6-3a_paired", metric = v, clade = "ALL", n = sum(keep), estimate = est,
      lo = ci(b)[1], hi = ci(b)[2], p_boot = bp(b), n_clades_pos = sum(cm > 0), n_clades = length(cm))
  for (k in names(dd)) {
    bk <- replicate(NB, mean(sample(dd[[k]], replace = TRUE)))
    put(analysis = "R6-3a_paired", metric = v, clade = k, n = length(dd[[k]]), estimate = mean(dd[[k]]),
        lo = ci(bk)[1], hi = ci(bk)[2], p_boot = bp(bk), n_clades_pos = NA, n_clades = NA)
  }
}
## R6-3a secondary: real-relaxed vs null1-relaxed (unpaired, clade-equal weight)
for (v in c("Dpur", "Dpos")) {
  r <- ok[ok$relaxed == 1 & !is.na(ok[[v]]), ]
  g <- split(r, r$clade); g <- g[sapply(g, function(z) all(c("real", "null1") %in% z$label))]
  f <- function(z) mean(z[[v]][z$label == "real"]) - mean(z[[v]][z$label == "null1"])
  est <- mean(sapply(g, f))
  b <- replicate(NB, mean(sapply(g, function(z) {
    zr <- z[z$label == "real", ]; zn <- z[z$label == "null1", ]
    mean(sample(zr[[v]], replace = TRUE)) - mean(sample(zn[[v]], replace = TRUE)) })))
  put(analysis = "R6-3a_relaxed_only", metric = v, clade = "ALL", n = nrow(r), estimate = est,
      lo = ci(b)[1], hi = ci(b)[2], p_boot = bp(b), n_clades_pos = sum(sapply(g, f) > 0), n_clades = length(g))
}
## R6-3b functional concentration: relaxed ~ label * inSet + clade, HOG-cluster bootstrap within clade
y <- x[x$has_go == 1 & !is.na(x$relax_p) & !is.na(x$relax_K), ]
y$real <- as.integer(y$label == "real")
hogs <- split(seq_len(nrow(y)), paste(y$clade, y$hog))
hc <- sapply(strsplit(names(hogs), " "), `[`, 1)
fitI <- function(d, s) {
  m <- suppressWarnings(glm(d$relaxed ~ d$real * d[[s]] + factor(d$clade), family = binomial))
  unname(coef(m)[paste0("d$real:d[[s]]")])
}
res3b <- list()
for (s in c("LIGHT", "CIRCADIAN", "PIGMENT")) {
  est <- fitI(y, s)
  b <- replicate(NB, {
    idx <- unlist(lapply(split(seq_along(hogs), hc), function(i) sample(i, replace = TRUE)))
    fitI(y[unlist(hogs[idx]), ], s) })
  # excess shares: sum_c n_real,c,S * (rate_real - rate_null) / total excess
  ex <- function(m) sum(sapply(split(y[m, ], y$clade[m]), function(z) {
    zr <- z[z$real == 1, ]; zn <- z[z$real == 0, ]
    nrow(zr) * (mean(zr$relaxed) - if (nrow(zn)) mean(zn$relaxed) else NA) }), na.rm = TRUE)
  inS <- y[[s]] == 1
  res3b[[s]] <- data.frame(analysis = "R6-3b_set_interaction", metric = s, clade = "ALL",
      n = sum(inS & y$real == 1), estimate = est, lo = ci(b)[1], hi = ci(b)[2], p_boot = bp(b),
      n_clades_pos = NA, n_clades = NA, excess_share = ex(inS) / ex(rep(TRUE, nrow(y))),
      hog_share = mean(inS[y$real == 1]), stringsAsFactors = FALSE)
}
r3 <- do.call(rbind, res3b); r3$p_holm <- p.adjust(r3$p_boot, "holm")
o <- do.call(rbind, out); o$excess_share <- NA; o$hog_share <- NA; o$p_holm <- NA
o <- rbind(o, r3)
write.table(o, file.path(D, "r6_3_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# genome-wide label effect for context (same data as R6-3b)
m0 <- glm(relaxed ~ real + factor(clade), family = binomial, data = y)
cat("genome-wide logOR real vs null1 (has_go fits):", coef(m0)["real"], "\n")
cat("rates:\n"); print(aggregate(relaxed ~ clade + label, data = x, FUN = mean))
print(o, digits = 3)
