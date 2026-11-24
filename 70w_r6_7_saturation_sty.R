# R6-7 models (ADDENDUM_R6_7.md md5 3250794e, POST HOC). Run from ~/deepsea_conv/work after 70w_r6_7_saturation_sty.py.
set.seed(20261105)
d <- read.delim("addon/R6/R6_7/r6_7_fits.tsv")
d <- d[is.finite(d$p) & is.finite(d$K), ]
d$relaxed <- as.integer(d$p < .05 & d$K < 1); d$real <- as.integer(d$label == "real")
g <- d[d$label == "real", c("clade","hog","D_gene")]
d$Dg <- g$D_gene[match(paste(d$clade, d$hog), paste(g$clade, g$hog))]; d$Dg[is.na(d$Dg)] <- d$D_gene[is.na(d$Dg)]
cat("fits", nrow(d), " real", sum(d$real), " null1", sum(1 - d$real), "\n")
rr <- tapply(d$relaxed, d$label, mean); cm <- tapply(d$relaxed, list(d$clade, d$label), mean)
cat(sprintf("GATE pooled relaxed: real %.4f null1 %.4f (E1 0.326 / 0.149, tol 0.005); clade-mean real %.4f null1 %.4f\n",
            rr["real"], rr["null1"], mean(cm[, "real"]), mean(cm[, "null1"])))
gate <- abs(rr["real"] - .326) <= .005 & abs(rr["null1"] - .149) <= .005
cat("GATE:", if (gate) "PASS" else "FAIL -> stop per addendum", "\n")
if (!gate) quit(save = "no")
## descriptive
cat("\nper-clade divergence (real fits): median D_gene, median max branch, share genes any branch > 1\n")
R <- d[d$real == 1, ]
print(round(do.call(rbind, lapply(split(R, R$clade), function(x) c(n = nrow(x), med_D = median(x$Dg),
  med_max = median(x$max_br), any_gt1 = mean(x$max_br > 1)))), 4))
key <- paste(d$clade, d$hog); idxl <- split(seq_len(nrow(d)), key); uk <- names(idxl); uc <- sub(" .*", "", uk)
boot <- function(f, B = 2000) replicate(B, { k <- unlist(lapply(unique(uc), function(c) { i <- uk[uc == c]; sample(i, length(i), TRUE) }))
  f(d[unlist(idxl[k], use.names = FALSE), ]) })
lor <- function(x, fo = relaxed ~ real + factor(clade)) coef(glm(fo, binomial, x))["real"]
ci <- function(b) quantile(b, c(.025, .975), na.rm = TRUE)
## P1
q4 <- quantile(d$Dg[d$real == 1], .75); d$top <- d$Dg > q4
cat(sprintf("\nP1 pooled top-quartile threshold D_gene > %.4f; genes in top quartile per clade:\n", q4))
print(table(d$clade[d$top & d$real == 1]))
all0 <- lor(d); top0 <- lor(d[d$top, ]); rest0 <- lor(d[!d$top, ])
bt <- boot(function(x) c(lor(x), lor(x[x$top, ]), lor(x[!x$top, ])))
cat(sprintf("real log OR all %.3f (%.3f-%.3f); top quartile %.3f (%.3f-%.3f); excluding top %.3f (%.3f-%.3f)\n",
            all0, ci(bt[1, ])[1], ci(bt[1, ])[2], top0, ci(bt[2, ])[1], ci(bt[2, ])[2], rest0, ci(bt[3, ])[1], ci(bt[3, ])[2]))
## P2
fo2 <- relaxed ~ real + log10(D_fg) + log10(Dg) + factor(clade)
a0 <- lor(d, fo2); b2 <- boot(function(x) lor(x, fo2))
cat(sprintf("P2 adjusted real log OR %.3f (%.3f-%.3f); unadjusted same rows %.3f; ratio %.3f\n",
            a0, ci(b2)[1], ci(b2)[2], all0, a0 / all0))
print(round(summary(glm(fo2, binomial, d))$coefficients[2:4, ], 4))
rob <- ci(bt[2, ])[1] > 0 & ci(b2)[1] > 0 & a0 / all0 >= .5
cat("DECISION (ADDENDUM_R6_7):", if (rob) "ROBUST" else "QUALIFIED", "\n")
## secondary
d$Q <- cut(d$Dg, quantile(d$Dg[d$real == 1], 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
m0 <- glm(relaxed ~ real + Q + factor(clade), binomial, d); m1 <- glm(relaxed ~ real * Q + factor(clade), binomial, d)
cat(sprintf("\ns1 quartile x real LRT P = %.4f; real log OR by quartile:\n", anova(m0, m1, test = "LRT")[2, "Pr(>Chi)"]))
print(round(sapply(levels(d$Q), function(q) lor(d[d$Q %in% q, ])), 3))
print(round(tapply(d$relaxed, list(d$Q, d$label), mean), 4))
cat("s2 per-clade real log OR in the within-clade top tertile of D_gene:\n")
print(round(sapply(split(d, d$clade), function(x) { t3 <- quantile(x$Dg[x$real == 1], 2 / 3)
  y <- x[x$Dg > t3, ]; coef(glm(relaxed ~ real, binomial, y))["real"] }), 3))
pr <- merge(d[d$real == 1, c("clade","hog","D_fg")], d[d$real == 0, c("clade","hog","D_fg")], by = c("clade","hog"))
cat(sprintf("s3 paired D_fg real/null1: n %d, median ratio %.3f, Wilcoxon P = %.3g\n", nrow(pr),
            median(pr$D_fg.x / pr$D_fg.y), wilcox.test(log(pr$D_fg.x), log(pr$D_fg.y), paired = TRUE)$p.value))
print(round(sapply(split(pr, pr$clade), function(x) median(x$D_fg.x / x$D_fg.y)), 3))
d$lK <- log10(pmin(pmax(d$K, .001), 50))
cat("s4 Spearman log10 K vs log10 D_fg, by clade x arm:\n")
print(round(sapply(split(d, d$clade), function(x) sapply(c(real = 1, null1 = 0), function(a) {
  y <- x[x$real == a, ]; cor(y$lK, log10(y$D_fg), method = "spearman") })), 3))
cat("s5 per-clade median D_gene and real log OR:\n")
print(round(t(sapply(split(d, d$clade), function(x) c(med_D = median(x$Dg[x$real == 1]),
  lor = unname(coef(glm(relaxed ~ real, binomial, x))["real"])))), 3))
