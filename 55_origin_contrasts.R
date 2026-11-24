#!/usr/bin/env Rscript
# A2 (PREREG_addon.md): physicochemical contrasts across independent deep-sea origins.
# Usage: Rscript 55_origin_contrasts.R [TREE] [OUTDIR]     (default work/final_v3/date/lsd.timetree.nwk, work/addon/A2)
# env LABEL_COL (lab), ASR_COL (asr), DATED (1), NCORES (8), NBOOT (2000), SEED (20261025), PART (all | contrast | influence)
# Origins = maximal all-D clades (addon_lib deep_clades) of the labelled tree in Ma, comparator = the S tips of the
#   sister subtree; contrast = mean z(D) - mean z(S) per trait, z = SD units of all labelled species (as 44);
#   route_adj: z minus the S_route3 route effects (y ~ deep + route, all species). Sets: disjoint (primary; greedy,
#   youngest stem first, the comparator loses S tips used by an accepted origin, skipped if none are left), strict
#   (skipped on any overlap), all (full comparators, reused), asr_disjoint (maximal deep clades of the ASR_COL labels
#   as in 31_asr_depth.R, tips restricted to the PGLS species with the same state; disjoint rule).
# Per set x adj x trait (P6 + GC set): Wilcoxon signed-rank (Holm over P6 within set x adj), sign test, weighted mean
#   (w = sqrt(nD nS / (nD + nS)), NBOOT origin resamples for the 95 % CI), route-matched origins (every tip one route;
#   descriptive if < 15), Spearman contrast vs stem age. Decision (P6): repeated = disjoint raw Holm p < 0.05 and the
#   median has the S_route3 sign; route_robust = repeated and the same holds for disjoint route_adj.
# Influence: S_route3 deep effect leaving out each clade with >= 6 D and >= 6 S and the D tips of each of the 10
#   largest origins (z kept on the all-species scale).
# Outputs (plain tsv): origins.tsv, contrasts.tsv, summary.tsv, decision.tsv, influence.tsv, report.txt, report_influence.txt
fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(Sys.getenv("ADDON_LIB", file.path(dirname(normalizePath(sub("^--file=", "", fa[1]))), "addon_lib.R")))
suppressPackageStartupMessages({ library(phylolm); library(parallel) })
a <- commandArgs(TRUE)
TREE <- if (length(a) >= 1) a[1] else file.path(W, "date/lsd.timetree.nwk")
OUT <- if (length(a) >= 2) a[2] else file.path(R, "work/addon/A2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (!nzchar(Sys.getenv("DATED"))) Sys.setenv(DATED = "1")
NC <- as.integer(Sys.getenv("NCORES", "8")); NBOOT <- as.integer(Sys.getenv("NBOOT", "2000"))
SEED <- as.integer(Sys.getenv("SEED", "20261025")); PART <- Sys.getenv("PART", "all"); ASR_COL <- Sys.getenv("ASR_COL", "asr")
lg <- function(...) { cat(format(Sys.time(), "%F %T"), ..., "\n"); flush(stdout()) }
# MA = 1: LSD2 tree rescaled to Ma before pruning (root age from lsd.timetree.lsd); MA = 0: TREE is already in Ma
tt <- if (Sys.getenv("MA", "1") == "1") load_tree(TREE, ma = TRUE) else read.tree(TREE)
L <- load_species(tt); sp <- L$sp; tr <- L$tr; rm(L)
TRT <- intersect(c(P6, GCSET), names(sp)); PR <- intersect(P6, TRT)
lg("labels", Sys.getenv("LABEL_COL", "lab"), ": D", sum(sp$deep), "S", sum(1 - sp$deep), "; tree height",
   round(max(node.depth.edgelength(tr)), 1), "Ma; traits", length(TRT))
yof <- function(tv) { d <- sp[is.finite(get(tv))]; d[, y := zs(get(tv))]; d }
ref <- rbindlist(lapply(TRT, function(tv) { r <- suppressWarnings(fit1(yof(tv), "deep", tr)); if (!is.null(r)) r[, trait := tv]; r }))
rd <- ref[term == "deep", .(trait, ref_est = est, ref_p = p)]
zval <- function(tv, adj) { d <- yof(tv); y <- setNames(d$y, d$code)
  if (adj) for (b in split(ref[trait == tv & startsWith(term, "route")], by = "term")) {
    k <- d$route == sub("^route", "", b$term); y[k] <- y[k] - b$est }
  y }
ZV <- list(raw = sapply(TRT, zval, adj = FALSE, simplify = FALSE), route_adj = sapply(TRT, zval, adj = TRUE, simplify = FALSE))

tips <- function(s) strsplit(s, ",", fixed = TRUE)[[1]]
# greedy, youngest stem first. trim: the comparator loses S tips already used and the origin is skipped if none are
# left (all S tips of the sister subtree share the same MRCA with the origin); strict: skipped on any overlap. Either
# way every species enters at most one contrast and the contrasts are uncorrelated under Brownian motion.
disjoint <- function(og, trim = TRUE) {
  og <- og[order(stem_age, -n_D, node)][, n_S_full := n_S]; used <- character(0); k <- logical(nrow(og)); S2 <- og$S_tips
  for (i in seq_len(nrow(og))) { d <- tips(og$D_tips[i]); s <- tips(og$S_tips[i])
    s <- if (trim) setdiff(s, used) else if (any(s %in% used)) character(0) else s
    if (length(s) && !any(d %in% used)) { k[i] <- TRUE; used <- c(used, d, s); S2[i] <- paste(s, collapse = ",") } }
  og[, S_tips := S2][k][, n_S := lengths(strsplit(S_tips, ",", fixed = TRUE))] }
og_all <- deep_clades(tr, setNames(sp$deep, sp$code))
SETS <- list(disjoint = disjoint(og_all), strict = disjoint(og_all, trim = FALSE), all = copy(og_all)[, n_S_full := n_S])
lb <- fread(file.path(R, "meta/labels_v2.tsv"), sep = "\t", quote = "")
if (ASR_COL %in% names(lb)) {
  st <- setNames(lb[[ASR_COL]], lb$code); st <- st[st %in% c("D", "S") & names(st) %in% tt$tip.label]
  oa <- deep_clades(keep.tip(tt, names(st)), setNames(as.integer(st == "D"), names(st)))
  rs <- function(v, keep) vapply(v, function(s) paste(intersect(tips(s), keep), collapse = ","), "", USE.NAMES = FALSE)
  kD <- sp[deep == 1, code]; kS <- sp[deep == 0, code]
  oa[, `:=`(D_tips = rs(D_tips, kD), S_tips = rs(S_tips, kS))]
  oa <- oa[nzchar(D_tips) & nzchar(S_tips)][, `:=`(n_D = lengths(strsplit(D_tips, ",")), n_S = lengths(strsplit(S_tips, ",")))]
  SETS$asr_disjoint <- disjoint(oa)
} else lg("labels column", ASR_COL, "absent: asr_disjoint skipped")
for (s in names(SETS)) if (!NROW(SETS[[s]])) { lg("set", s, "has no origins: skipped"); SETS[[s]] <- NULL }
rt_of <- setNames(sp$route, sp$code); cl_of <- setNames(sp$clade, sp$code)
SETS <- lapply(SETS, function(og) { og <- copy(og)
  og[, route_match := vapply(paste(D_tips, S_tips, sep = ","), function(s) { r <- unique(rt_of[tips(s)])
    if (length(r) == 1L) r else NA_character_ }, "", USE.NAMES = FALSE)]
  og[, clade := vapply(D_tips, function(s) names(which.max(table(cl_of[tips(s)]))), "", USE.NAMES = FALSE)]
  og })
lg("origins: all", nrow(SETS$all), "disjoint", NROW(SETS$disjoint), "strict", NROW(SETS$strict), "asr_disjoint", NROW(SETS$asr_disjoint))

con1 <- function(og, z) rbindlist(lapply(seq_len(nrow(og)), function(i) {
  d <- z[tips(og$D_tips[i])]; s <- z[tips(og$S_tips[i])]; d <- d[is.finite(d)]; s <- s[is.finite(s)]
  data.table(node = og$node[i], stem_age = og$stem_age[i], route_match = og$route_match[i], clade = og$clade[i],
             nd = length(d), ns = length(s), contrast = if (length(d) && length(s)) mean(d) - mean(s) else NA_real_) }))
summ <- function(d) {
  d <- d[is.finite(contrast)]; n <- nrow(d); if (n < 3) return(data.table(n = n))
  w <- sqrt(d$nd * d$ns / (d$nd + d$ns)); x <- d$contrast; np <- sum(x > 0); nz <- sum(x != 0)
  set.seed(SEED); bs <- replicate(NBOOT, { i <- sample.int(n, n, TRUE); sum(w[i] * x[i]) / sum(w[i]) })
  wx <- suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE))
  sr <- suppressWarnings(cor.test(x, d$stem_age, method = "spearman", exact = FALSE)); xr <- x[!is.na(d$route_match)]
  data.table(n = n, median = median(x), mean = mean(x), wmean = sum(w * x) / sum(w), wmean_lo = unname(quantile(bs, 0.025)),
             wmean_hi = unname(quantile(bs, 0.975)), n_pos = np, sign_p = if (nz) binom.test(np, nz)$p.value else NA_real_,
             wilcox_V = unname(wx$statistic), wilcox_p = wx$p.value, spearman_rho = unname(sr$estimate), spearman_p = sr$p.value,
             rm_n = length(xr), rm_median = if (length(xr)) median(xr) else NA_real_,
             rm_wilcox_p = if (length(xr) >= 3) suppressWarnings(wilcox.test(xr, exact = FALSE))$p.value else NA_real_,
             rm_descriptive = length(xr) < 15)
}
if (PART %in% c("all", "contrast")) {
  fwrite(rbindlist(SETS, idcol = "set", use.names = TRUE, fill = TRUE), file.path(OUT, "origins.tsv"), sep = "\t", quote = FALSE)
  CON <- rbindlist(lapply(names(SETS), function(s) rbindlist(lapply(names(ZV), function(ad) rbindlist(lapply(TRT, function(tv)
    con1(SETS[[s]], ZV[[ad]][[tv]])[, `:=`(set = s, adj = ad, trait = tv)]))))))
  fwrite(CON, file.path(OUT, "contrasts.tsv"), sep = "\t", quote = FALSE, na = "NA")
  SM <- rbindlist(lapply(split(CON, by = c("set", "adj", "trait")), function(d) cbind(d[1, .(set, adj, trait)], summ(d))), fill = TRUE)
  SM[trait %in% PR, wilcox_holm := p.adjust(wilcox_p, "holm"), by = .(set, adj)]
  SM <- merge(SM, rd, by = "trait", all.x = TRUE, sort = FALSE)[, agree := sign(median) == sign(ref_est)]
  fwrite(SM, file.path(OUT, "summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
  DEC <- rd[trait %in% PR]
  for (k in list(c("disjoint", "raw"), c("disjoint", "route_adj"), c("strict", "raw"), c("all", "raw"), c("asr_disjoint", "raw"))) {
    x <- SM[set == k[1] & adj == k[2]]; i <- match(DEC$trait, x$trait)
    DEC[, (paste(k, collapse = "_")) := x$wilcox_holm[i] < 0.05 & x$agree[i]] }
  DEC[, `:=`(repeated = disjoint_raw %in% TRUE, route_robust = disjoint_raw %in% TRUE & disjoint_route_adj %in% TRUE)]
  fwrite(DEC, file.path(OUT, "decision.tsv"), sep = "\t", quote = FALSE, na = "NA")
  nd <- sum(lengths(strsplit(SETS$disjoint$D_tips, ",")))
  rp <- c(sprintf("A2 origin contrasts %s; labels %s; tree %s", format(Sys.time(), "%F %T"), Sys.getenv("LABEL_COL", "lab"), TREE),
          sprintf("origins: all %d, disjoint %d (D tips %d of %d, S tips %d of %d), strict %d, asr_disjoint %s", nrow(SETS$all),
                  nrow(SETS$disjoint), nd, sum(sp$deep), sum(SETS$disjoint$n_S), sum(1 - sp$deep), NROW(SETS$strict),
                  if (is.null(SETS$asr_disjoint)) "skipped" else nrow(SETS$asr_disjoint)))
  for (tv in PR) { r <- SM[trait == tv & set == "disjoint" & adj == "raw"]; q <- SM[trait == tv & set == "disjoint" & adj == "route_adj"]
    if (!nrow(r) || !is.finite(r$n) || r$n < 3) { rp <- c(rp, sprintf("%-26s fewer than 3 origins", tv)); next }
    rp <- c(rp, sprintf("%-26s ref %+.3f (p %.2g) | disjoint n %d median %+.3f wmean %+.3f [%+.3f, %+.3f] pos %d Wilcoxon p %.2g Holm %.2g | route_adj median %+.3f Holm %.2g | rm n %d | repeated %s route_robust %s",
      tv, r$ref_est, r$ref_p, r$n, r$median, r$wmean, r$wmean_lo, r$wmean_hi, r$n_pos, r$wilcox_p, r$wilcox_holm, q$median,
      q$wilcox_holm, r$rm_n, DEC[trait == tv, repeated], DEC[trait == tv, route_robust])) }
  writeLines(rp, file.path(OUT, "report.txt")); cat(rp, sep = "\n")
}
if (PART %in% c("all", "influence")) {
  cls <- sp[, .(nD = sum(deep), nS = sum(1 - deep)), by = clade][nD >= 6 & nS >= 6][order(clade), clade]
  big <- SETS$all[order(-n_D, node)][seq_len(min(10L, nrow(SETS$all)))]
  jb <- c(lapply(cls, function(c) list(type = "clade", id = c, drop = sp[clade == c, code])),
          lapply(seq_len(nrow(big)), function(i) list(type = "origin", drop = tips(big$D_tips[i]),
                 id = sprintf("%s_nD%d_%s", big$clade[i], big$n_D[i], tips(big$D_tips[i])[1]))))
  gr <- CJ(j = seq_along(jb), tv = PR)
  res <- mclapply(seq_len(nrow(gr)), function(k) { b <- jb[[gr$j[k]]]; tv <- gr$tv[k]
    r <- suppressWarnings(fit1(yof(tv)[!code %in% b$drop], "deep", tr)); if (is.null(r)) return(NULL)
    r[term == "deep", .(trait = tv, left_out = b$type, id = b$id, n_drop = length(b$drop), est, se, p, n, n_deep)] },
    mc.cores = NC, mc.preschedule = FALSE)
  IN <- rbindlist(Filter(is.data.table, res)); lg("influence fits", nrow(IN), "of", nrow(gr))
  IN <- merge(IN, rd, by = "trait", sort = FALSE)[, same_sign_sig := sign(est) == sign(ref_est) & p < 0.05]
  fwrite(IN, file.path(OUT, "influence.tsv"), sep = "\t", quote = FALSE, na = "NA")
  ri <- IN[, .(n_fits = .N, kept = sum(same_sign_sig), est_min = min(est), est_max = max(est), p_max = max(p),
               worst = id[which.max(p)]), by = .(trait, left_out)]
  rp <- c(sprintf("A2 influence %s: clades %s; largest origins %d", format(Sys.time(), "%F %T"), paste(cls, collapse = ","), nrow(big)),
          sprintf("%-26s %-6s kept p<0.05 & sign %d/%d, est %+.3f..%+.3f, max p %.2g (%s)", ri$trait, ri$left_out, ri$kept, ri$n_fits,
                  ri$est_min, ri$est_max, ri$p_max, ri$worst))
  writeLines(rp, file.path(OUT, "report_influence.txt")); cat(rp, sep = "\n")
}
lg("done")
