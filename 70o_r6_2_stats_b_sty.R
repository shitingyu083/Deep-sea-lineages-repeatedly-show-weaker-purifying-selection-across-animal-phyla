# R6-2b / combined statistics (ADDENDUM_R6_4.md part 1). Copy of 70k_r6_2_stats_sty.R (statistics unchanged) with
# input = one or more fit directories (each with r6_2_fits.tsv, own CUT / cut100.txt filter) and output directory;
# additionally writes the per-replicate bootstrap lineage effects (boot_lineage_eff.tsv) for R6-4.
# Usage: Rscript 70o_r6_2_stats_b_sty.R <OUT dir> <tag> <fit dir> [<fit dir> ...]
a <- commandArgs(TRUE); OUTD <- a[1]; TAG <- a[2]; IN <- a[-(1:2)]
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
nall <- 0; att <- NULL
d <- do.call(rbind, lapply(IN, function(G) {
  x <- read.delim(file.path(G, "r6_2_fits.tsv"), colClasses = "character"); nall <<- nall + nrow(x)
  if (file.exists(file.path(G, "CUT"))) {
    cu <- read.table(file.path(G, "cut100.txt"), colClasses = "character")
    x <- x[paste(x$clade, x$hog) %in% paste(cu$V1, cu$V2), ]
  }
  x
}))
stopifnot(!any(duplicated(d[, c("clade", "hog", "set", "group")])))
att <- as.data.frame(table(set = d$set, fit = d$fit))
d <- d[d$fit == "ok" & d$canon != "merged", ]
d$logK <- log10(pmin(pmax(as.numeric(d$K), .001), 50))
d$p <- as.numeric(d$p_joint)

eff <- function(x, L) {
  r <- x[x$set == "real", ]; n <- x[x$set == "null1", ]
  nm <- tapply(n$logK, n$clade, median)
  lm <- tapply(r$logK, paste(r$clade, r$canon), median)
  e <- lm[paste(L$clade, L$canon)] - nm[L$clade]
  data.frame(clade = L$clade, canon = L$canon, eff = as.numeric(e))
}
pool <- function(e, cl = NULL) {
  if (!is.null(cl)) e <- e[e$clade %in% cl, ]
  e <- e[!is.na(e$eff), ]
  k <- table(e$clade); ok <- names(k)[k >= 2]
  if (!length(ok)) return(NA)
  mean(tapply(e$eff, e$clade, mean)[ok])
}
r <- d[d$set == "real", ]
nh <- tapply(r$hog, paste(r$clade, r$canon), function(v) length(unique(v)))
L <- unique(r[, c("clade", "canon")]); L <- L[nh[paste(L$clade, L$canon)] >= 5, ]
L <- L[order(L$clade, L$canon), ]
L$n_hog <- as.integer(nh[paste(L$clade, L$canon)])
cn <- read.delim(file.path(IN[1], "canon.tsv"), colClasses = "character")
L$size <- as.numeric(cn$n_species[match(paste(L$clade, L$canon), paste(cn$clade, cn$canon))])
E <- eff(d, L); E$n_hog <- L$n_hog; E$size <- L$size
nl <- table(E$clade); cl2 <- names(nl)[nl >= 2]
top2 <- names(sort(nl[cl2], decreasing = TRUE))[seq_len(min(2, length(cl2)))]

drop_largest <- function(e) { e <- e[!is.na(e$eff), ]
  do.call(rbind, lapply(split(e, e$clade), function(z) z[-which.max(z$size), ])) }
V <- list(pooled = function(e) pool(e),
          lineage_weighted = function(e) mean(e$eff, na.rm = TRUE),
          share_lt0 = function(e) mean(e$eff < 0, na.rm = TRUE),
          without_top2_clades = function(e) pool(e, setdiff(cl2, top2)),
          without_largest_lineage = function(e) pool(drop_largest(e)))
for (c in cl2) V[[paste0("loco_", c)]] <- local({ cc <- c; function(e) pool(e, setdiff(cl2, cc)) })
for (c in cl2) V[[paste0("clade_", c)]] <- local({ cc <- c; function(e) mean(e$eff[e$clade == cc], na.rm = TRUE) })
obs <- sapply(V, function(f) f(E))

set.seed(20261031)
byc <- split(unique(d[, c("clade", "hog")]), unique(d[, c("clade", "hog")])$clade)
idx <- split(seq_len(nrow(d)), paste(d$clade, d$hog))
nv <- length(V)
B <- t(replicate(2000, {
  s <- unlist(lapply(byc, function(z) { h <- sample(z$hog, nrow(z), TRUE); unlist(idx[paste(z$clade[1], h)], use.names = FALSE) }))
  e <- eff(d[s, ], L); e$size <- L$size
  c(sapply(V, function(f) f(e)), e$eff)
}))
BL <- B[, -(1:nv), drop = FALSE]; B <- B[, 1:nv, drop = FALSE]; colnames(B) <- names(V)
colnames(BL) <- paste(L$clade, L$canon, sep = ":")
ci <- t(apply(B, 2, quantile, c(.025, .975), na.rm = TRUE))
res <- data.frame(stat = names(obs), est = obs, lo = ci[, 1], hi = ci[, 2], row.names = NULL)

pc <- obs["pooled"]; claim <- !is.na(pc) && pc < 0 && ci["pooled", 2] < 0 && obs["share_lt0"] > .5 &&
  !is.na(obs["without_top2_clades"]) && obs["without_top2_clades"] < 0
rob <- function(s) !is.na(obs[s]) && ci[s, 2] < 0
jp <- tapply(d$p < .05, d$set, mean)
write.table(res, file.path(OUTD, "r6_2_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(E, file.path(OUTD, "r6_2_lineages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(BL, file.path(OUTD, "boot_lineage_eff.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sink(file.path(OUTD, "r6_2_report.txt"))
cat("R6-2", TAG, "per-origin RELAX (ADDENDUM_R6_4 part 1; statistics = ADDENDUM 0e9d4fe6)\ninputs:", IN,
    "\nrows in fits:", nall, " CUT:", sapply(IN, function(G) file.exists(file.path(G, "CUT"))),
    " STOP:", sapply(IN, function(G) file.exists(file.path(G, "STOP"))), "\nfit status after cut filter:\n")
print(att[att$Freq > 0, ])
cat("\nsummarised lineages:", nrow(E), " clades with >=2 lineages:", paste(cl2, collapse = ","),
    "\ntop2 clades by lineage count:", paste(top2, collapse = ","), "\n\n")
print(res, digits = 3); cat("\nper-lineage effects:\n"); print(E, digits = 3)
cat("\njoint p<0.05 share by set (reported, not used):\n"); print(round(jp, 3))
cat("\nRULE (pooled<0 & CI excl 0 & share>.5 & without top2 <0):", if (claim) "MET" else "NOT MET", "\n")
cat("ROBUST without_top2_clades (CI excl 0):", rob("without_top2_clades"),
    " ROBUST without_largest_lineage (CI excl 0):", rob("without_largest_lineage"), "\n")
sink()
cat("ok\n")
