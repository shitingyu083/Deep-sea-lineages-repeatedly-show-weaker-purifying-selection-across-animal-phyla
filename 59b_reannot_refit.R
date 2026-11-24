#!/usr/bin/env Rscript
# A1 re-annotation control, corrected refit (PREREG_addon.md A1; input deltas.tsv of 59_route_reannot.py --collect).
# delta = annot - miniprot (raw units) of a re-annotated annot-route species; shift_L = mean delta of the re-annotated
# species with label L (D / S). Every miniprot-route species of label L gets y + shift_L / sd (y in SD units of all
# labelled species of the original trait, as 44 / 53: the corrected estimate keeps the original scale), then
# S_route3 (deep + route, = 44 ref_all) and S_noroute (deep) are refitted (fit1, Pagel's lambda). Bootstrap: the
# re-annotated species resampled with replacement within D and within S (NBOOT draws, the same draw for every trait),
# shifts recomputed, refitted; percentile 95 % interval of the corrected deep effect. (h) "re-annotation consistent":
# the corrected S_route3 deep effect has the sign of the original, p < 0.05 and a bootstrap interval excluding 0;
# not evaluable with < 5 re-annotated species of either label.
# Usage: Rscript 59b_reannot_refit.R [A1_DIR] [TREE]  (default work/addon/A1_reannot, work/final_v3/date/lsd.timetree.nwk)
# env LABEL_COL (lab), DATED (1), NCORES (8), NBOOT (200), SEED (20261025)
# Outputs in A1_DIR: refit.tsv, refit_boot.tsv, refit_report.txt, REFIT.DONE
fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SELF <- normalizePath(sub("^--file=", "", fa[1]))
source(Sys.getenv("ADDON_LIB", file.path(dirname(SELF), "addon_lib.R")))
suppressPackageStartupMessages({ library(phylolm); library(parallel) })
a <- commandArgs(TRUE)
A1 <- if (length(a) >= 1) a[1] else file.path(R, "work/addon/A1_reannot")
TREE <- if (length(a) >= 2) a[2] else file.path(W, "date/lsd.timetree.nwk")
if (!nzchar(Sys.getenv("DATED"))) Sys.setenv(DATED = "1")
NC <- as.integer(Sys.getenv("NCORES", "8")); NBOOT <- as.integer(Sys.getenv("NBOOT", "200"))
SEED <- as.integer(Sys.getenv("SEED", "20261025")); NMIN <- 5L
lg <- function(...) { cat(format(Sys.time(), "%F %T"), ..., "\n"); flush(stdout()) }
fw <- function(...) suppressWarnings(fit1(...))
head0 <- sprintf("59b_reannot_refit.R %s (md5 %s); A1 %s; tree %s; label %s; NBOOT %d; SEED %d", format(Sys.time(), "%F %T"),
                 unname(tools::md5sum(SELF)), A1, TREE, Sys.getenv("LABEL_COL", "lab"), NBOOT, SEED)
dl <- fread(file.path(A1, "deltas.tsv"), sep = "\t", colClasses = "character")
if (!all(c("code", "lab", "trait", "delta_raw") %in% names(dl)) || !nrow(dl)) {
  writeLines(c(head0, "no completed re-annotation: nothing refitted"), file.path(A1, "refit_report.txt"))
  writeLines(format(Sys.time(), "%F %T no data"), file.path(A1, "REFIT.DONE")); quit(save = "no", status = 0)
}
L <- load_species(read.tree(TREE)); sp <- L$sp; tr <- L$tr; rm(L)
dl <- dl[, .(code, lab, trait, delta = num(delta_raw))][lab %in% c("D", "S") & is.finite(delta)]
TR <- intersect(c(P6, GCSET), intersect(names(sp), unique(dl$trait)))
RC <- list(D = sort(unique(dl[lab == "D", code])), S = sort(unique(dl[lab == "S", code])))
cd <- sort(unique(dl$code)); M <- matrix(NA_real_, length(cd), length(TR), dimnames = list(cd, TR))
x <- dl[trait %in% TR]; M[cbind(match(x$code, cd), match(x$trait, TR))] <- x$delta
mshift <- function(dc, sc) {   # label-specific mean delta per trait (raw units); codes with multiplicity
  f <- function(cc) { v <- colMeans(M[cc, , drop = FALSE], na.rm = TRUE); v[!is.finite(v)] <- NA_real_; v }
  rbind(D = f(dc), S = f(sc))
}
nre <- rbind(D = colSums(is.finite(M[RC$D, , drop = FALSE])), S = colSums(is.finite(M[RC$S, , drop = FALSE])))
lg("labels: D", sum(sp$deep), "S", sum(1 - sp$deep), "; route", paste(names(table(sp$route)), table(sp$route), collapse = " "),
   "; re-annotated D", length(RC$D), "S", length(RC$S), "; traits", length(TR))
D0 <- lapply(setNames(TR, TR), function(tv) { d <- sp[is.finite(get(tv))]; v <- d[[tv]]
  d[, `:=`(y0 = (v - mean(v)) / sd(v), xsd = sd(v))]; d })
SPEC <- list(S_route3 = "route", S_noroute = NULL)
fit_tv <- function(tv, sh, ver) {   # deep effect of both specs; sh = c(D, S) raw shifts, NULL = original
  d <- copy(D0[[tv]])
  if (is.null(sh)) d[, y := y0] else {
    if (!all(is.finite(sh))) return(NULL)
    d[, y := y0 + fifelse(route == "miniprot", fifelse(deep == 1L, sh[["D"]], sh[["S"]]) / xsd, 0)]
  }
  rbindlist(lapply(names(SPEC), function(nm) {
    r <- fw(d, "deep", tr, rv = SPEC[[nm]])
    if (is.null(r)) NULL else r[term == "deep"][, `:=`(trait = tv, spec = nm, version = ver)]
  }), fill = TRUE)
}
if (!length(TR)) stop("59b: no P6 / GC-set trait in both deltas.tsv and species_metrics")
S0 <- mshift(RC$D, RC$S)
res <- mclapply(TR, function(tv) rbindlist(list(fit_tv(tv, NULL, "original"), fit_tv(tv, S0[, tv], "corrected")), fill = TRUE),
                mc.cores = min(NC, length(TR)))
bad <- !vapply(res, is.data.table, TRUE); if (any(bad)) lg("WARNING failed traits:", TR[bad])
res <- rbindlist(res[!bad], fill = TRUE)
set.seed(SEED)
BD <- lapply(seq_len(NBOOT), function(b) lapply(RC, function(v) v[sample.int(length(v), length(v), replace = TRUE)]))
bt <- mclapply(seq_len(NBOOT), function(b) {
  sh <- mshift(BD[[b]]$D, BD[[b]]$S)
  rbindlist(lapply(TR, function(tv) {
    r <- fit_tv(tv, sh[, tv], "boot")
    if (is.null(r) || !nrow(r)) NULL else r[, .(b = b, trait, spec, shift_D = sh["D", tv], shift_S = sh["S", tv], est, p)]
  }))
}, mc.cores = NC)
okb <- vapply(bt, is.data.table, TRUE); if (any(!okb)) lg("WARNING bootstrap workers failed:", sum(!okb))
BT <- rbindlist(bt[okb], fill = TRUE)
fwrite(BT, file.path(A1, "refit_boot.tsv"), sep = "\t")
if (nrow(BT)) BT <- BT[is.finite(est)]
sdt <- data.table(trait = TR, sd_trait = vapply(TR, function(tv) D0[[tv]]$xsd[1], 0), shift_D = S0["D", TR],
                  shift_S = S0["S", TR], n_re_D = nre["D", TR], n_re_S = nre["S", TR])
sdt[, `:=`(shift_D_sd = shift_D / sd_trait, shift_S_sd = shift_S / sd_trait)]
out <- merge(res, sdt, by = "trait")
if (nrow(BT)) {
  bq <- BT[, .(boot_n = .N, boot_lo = quantile(est, 0.025), boot_med = median(est), boot_hi = quantile(est, 0.975),
               boot_frac_p05 = mean(p < 0.05)), by = .(trait, spec)][, version := "corrected"]
  out <- merge(out, bq, by = c("trait", "spec", "version"), all.x = TRUE)
} else out[, `:=`(boot_n = 0L, boot_lo = NA_real_, boot_med = NA_real_, boot_hi = NA_real_, boot_frac_p05 = NA_real_)]
out <- merge(out, out[version == "original", .(trait, spec, est_orig = est, p_orig = p)], by = c("trait", "spec"), all.x = TRUE)
out[, h := NA_character_]
out[spec == "S_route3" & version == "corrected", h := fifelse(n_re_D < NMIN | n_re_S < NMIN, "not_evaluable",
  fifelse(sign(est) == sign(est_orig) & p < 0.05 & is.finite(boot_lo) &
            ((est_orig > 0 & boot_lo > 0) | (est_orig < 0 & boot_hi < 0)), "consistent", "not_consistent"))]
setorder(out, trait, spec, -version)
fwrite(out, file.path(A1, "refit.tsv"), sep = "\t")
f3 <- function(v) ifelse(is.finite(v), formatC(v, digits = 3, format = "g"), "NA")
cs <- file.path(W, "pgls_controls/controls_summary.tsv"); chk <- "44 ref_all reproduction: controls_summary.tsv not found"
if (file.exists(cs)) {   # original S_route3 should reproduce 44 ref_all (column all_est) on the final tree
  c44 <- fread(cs, sep = "\t", colClasses = "character")
  if ("all_est" %in% names(c44)) {
    z <- merge(out[spec == "S_route3" & version == "original", .(trait, est)], c44[, .(trait, e44 = num(all_est))], by = "trait")
    chk <- sprintf("44 ref_all reproduction (original S_route3): %d traits, max |est - all_est| = %s", nrow(z),
                   f3(if (nrow(z)) max(abs(z$est - z$e44), na.rm = TRUE) else NA_real_))
  }
}
rp <- c(head0, chk, sprintf("re-annotated species with a delta: D %d, S %d; miniprot-route species shifted: D %d, S %d",
                            length(RC$D), length(RC$S), sp[route == "miniprot", sum(deep)], sp[route == "miniprot", sum(1 - deep)]),
        "trait\tspec\tn_re D/S\tshift D/S (SD)\toriginal est (p)\tcorrected est (p)\tboot 95% [lo, hi] (n)\t(h)")
for (i in seq_len(nrow(out[version == "corrected"]))) {
  o <- out[version == "corrected"][i]
  rp <- c(rp, sprintf("%s\t%s\t%d/%d\t%s/%s\t%s (%s)\t%s (%s)\t[%s, %s] (%d)\t%s", o$trait, o$spec, o$n_re_D, o$n_re_S,
                      f3(o$shift_D_sd), f3(o$shift_S_sd), f3(o$est_orig), f3(o$p_orig), f3(o$est), f3(o$p),
                      f3(o$boot_lo), f3(o$boot_hi), as.integer(o$boot_n), ifelse(is.na(o$h), "", o$h)))
}
writeLines(rp, file.path(A1, "refit_report.txt")); cat(rp, sep = "\n")
writeLines(format(Sys.time(), "%F %T"), file.path(A1, "REFIT.DONE"))
lg("done")
