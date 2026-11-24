# DEVIATIONS_r6 (e) Delta: R6-1 U-arm log OR (real vs null1, clade fixed effects; same model as 70g m0)
# on the random-stratum HOGs with replicate calls, replicate calls minus original calls.
# Failed replicate fits: reported; Delta uses HOGs where both U-arm replicate fits exist (paired sample).
x <- read.delim("~/deepsea_conv/work/addon/R6/refit/repl844_compare.tsv", stringsAsFactors = FALSE)
u <- x[x$stratum == "random" & x$arm == "U", ]
both <- tapply(u$repl == "ok", u$hog, all); keep <- names(both)[both]
u <- u[u$hog %in% keep, ]; u$real <- as.integer(u$set == "real")
lor <- function(y) unname(coef(glm(y ~ u$real + factor(u$clade), family = binomial))[2])
o <- lor(u$rel_o); r <- lor(u$rel_r)
cat(sprintf("U arm: %d HOGs with both replicate fits (of %d random-stratum HOGs)\n", length(keep), length(unique(x$hog[x$stratum == "random"]))))
cat(sprintf("log OR original %.3f  replicate %.3f  Delta %.3f\n", o, r, r - o))
