# R6-1 statistics (PREREG_r6.md) on 70e output. Usage: Rscript 70h_r6_route_stats_sty.R [DIR]
args <- commandArgs(TRUE); D <- if (length(args)) args[1] else "~/deepsea_conv/work/addon/R6"
x <- read.delim(file.path(D, "r6_1_fits.tsv"), na.strings = "NA", stringsAsFactors = FALSE)
NB <- 2000; set.seed(20261031)
x$real <- as.integer(x$set == "real")
x$relaxed <- as.integer(x$fit == "ok" & x$relax_p < 0.05 & x$relax_K < 1)
x$intens <- as.integer(x$fit == "ok" & x$relax_p < 0.05 & x$relax_K > 1)
x$Kbound <- as.integer(x$fit == "ok" & (x$relax_K < 0.01 | x$relax_K > 49))
ou <- x[x$arm %in% c("O", "U"), ]
full <- tapply(ou$fit == "ok", paste(ou$clade, ou$hog), all)
x$key <- paste(x$clade, x$hog); A <- x[x$key %in% names(full)[full], ]
out <- list(); put <- function(...) out[[length(out) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
ci <- function(b) quantile(b, c(.025, .975), na.rm = TRUE, names = FALSE)
lor <- function(d, y = "relaxed") {
  if (length(unique(d$clade)) > 1) unname(coef(suppressWarnings(glm(d[[y]] ~ d$real + factor(d$clade), family = binomial)))[2])
  else { t <- table(factor(d$real, 0:1), factor(d[[y]], 0:1)) + 0.5 * any(table(factor(d$real, 0:1), factor(d[[y]], 0:1)) == 0)
         log(t[2, 2] * t[1, 1] / (t[2, 1] * t[1, 2])) }
}
inter <- function(d, y = "relaxed") {
  d <- d[d$arm %in% c("O", "U"), ]; u <- as.integer(d$arm == "U")
  unname(coef(suppressWarnings(glm(d[[y]] ~ d$real * u + factor(d$clade), family = binomial)))["d$real:u"])
}
stat <- function(d, y) {
  lu <- lor(d[d$arm == "U", ], y); lo <- lor(d[d$arm == "O", ], y)
  c(U = lu, O = lo, P = lor(d[d$arm == "P" & d$fit == "ok", ], y), inter = inter(d, y), R = lu / lo)
}
idx <- split(seq_len(nrow(A)), A$key); kc <- sub(" .*", "", names(idx))
boot <- function(f) replicate(NB, { k <- unlist(lapply(split(seq_along(idx), kc), function(i) sample(i, replace = TRUE)))
  f(A[unlist(idx[k]), ]) })
for (y in c("relaxed", "intens")) {
  est <- stat(A, y); b <- boot(function(d) stat(d, y))
  for (k in names(est)) put(outcome = y, estimand = k, clade = "ALL", n_hog = length(idx), estimate = est[[k]],
                            lo = ci(b[k, ])[1], hi = ci(b[k, ])[2])
}
for (cl in unique(A$clade)) {
  d <- A[A$clade == cl, ]; ii <- split(seq_len(nrow(d)), d$key)
  e <- lor(d[d$arm == "U", ]); eo <- lor(d[d$arm == "O", ])
  b <- replicate(NB, { dd <- d[unlist(ii[sample(length(ii), replace = TRUE)]), ]; lor(dd[dd$arm == "U", ]) })
  put(outcome = "relaxed", estimand = "U_clade", clade = cl, n_hog = length(ii), estimate = e, lo = ci(b)[1], hi = ci(b)[2])
  put(outcome = "relaxed", estimand = "O_clade", clade = cl, n_hog = length(ii), estimate = eo, lo = NA, hi = NA)
  dl <- A[A$clade != cl & A$arm == "U", ]
  put(outcome = "relaxed", estimand = "U_leave_one_clade_out", clade = cl, n_hog = length(unique(dl$key)),
      estimate = lor(dl), lo = NA, hi = NA)
}
o <- do.call(rbind, out)
write.table(o, file.path(D, "r6_1_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# paired log K (U - O), agreement, diagnostics
w <- reshape(A[A$arm %in% c("O", "U"), c("key", "clade", "set", "arm", "relax_K", "relaxed")],
             idvar = c("key", "clade", "set"), timevar = "arm", direction = "wide")
sink(file.path(D, "r6_1_report.txt"))
cat("analysis set HOGs:", length(idx), "of built", length(full), "\n")
print(o, digits = 3)
for (s in c("real", "null1")) {
  z <- w[w$set == s, ]; dk <- log(z$relax_K.U) - log(z$relax_K.O)
  t <- table(O = z$relaxed.O, U = z$relaxed.U); po <- sum(diag(t)) / sum(t)
  pe <- sum(rowSums(t) * colSums(t)) / sum(t)^2
  cat("\n", s, ": median logK(U)-logK(O) =", median(dk), " Wilcoxon p =", wilcox.test(dk)$p.value,
      "\n relaxed agreement O vs U:\n"); print(t); cat(" kappa =", (po - pe) / (1 - pe), "\n")
}
cat("\nDiagnostics by set x arm (all built HOGs):\n")
print(aggregate(cbind(fit_fail = fit != "ok", Kbound, L, gap_frac) ~ set + arm, data = x[x$arm != "P", ],
                FUN = function(v) round(mean(v, na.rm = TRUE), 4), na.action = na.pass))
cat("\nRates by clade x set x arm (analysis set):\n")
print(aggregate(cbind(relaxed, intens) ~ clade + set + arm, data = A, FUN = function(v) round(mean(v), 3)))
g <- o[o$outcome == "relaxed" & o$clade == "ALL", ]; U <- g[g$estimand == "U", ]; Rr <- g[g$estimand == "R", ]
dec <- if (U$estimate > 0 && U$lo > 0 && Rr$estimate >= 0.5) "SUCCESS" else if (U$estimate > 0 && U$lo > 0) "ATTENUATED" else "FAIL"
cat("\nPRE-DECLARED DECISION:", dec, " logOR_U =", U$estimate, "[", U$lo, ",", U$hi, "]  R =", Rr$estimate, "\n")
sink()
cat(readLines(file.path(D, "r6_1_report.txt")), sep = "\n")
