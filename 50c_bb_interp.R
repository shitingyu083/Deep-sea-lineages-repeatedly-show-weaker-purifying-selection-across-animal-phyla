#!/usr/bin/env Rscript
# Ages of all nodes of the final tree from the backbone MCMCTree posterior (2026-10-25; PREREG_addon A6,
# DEVIATIONS_fig4 2026-10-25b). Full tree = the ML tree of the same topology (work/final_v3/bl.treefile, substitutions
# per site), pruned to the time-tree tips, rooted on the time tree's root split (root branch halved, branches
# < 1e-8 set to 1e-8). For each of the 100 posterior draws of 50b_mcmctree_collect.py (pooled clock-2 chains) the
# backbone node ages are fixed (age.min = age.max) and the other node ages estimated by penalized likelihood
# (ape::chronos, correlated rates, lambda 1). Point tree = the same with the posterior-mean ages of the clock-2 chains;
# sensitivity: chronos relaxed and clock models, and the clock-3 posterior means (correlated).
# Writes OUT/interp/{draws.nwk (chronograms, Ma), point_<variant>.nwk, node_ages.tsv, report.txt}; node_ages.tsv
# = per full-tree node: two defining tips, n_tips, backbone flag, LSD2 age, point age, draw median / 2.5 / 97.5 %,
# sensitivity ages (Ma).
# Usage: Rscript 50c_bb_interp.R [OUT = work/mcmctree_bb] [NCORES = 8] [NDRAW = all]
suppressPackageStartupMessages({ library(ape); library(data.table); library(parallel) })
source(file.path(Sys.getenv("DSC_ROOT", path.expand("~/deepsea_conv")), "scripts/addon_lib.R"))
a <- commandArgs(TRUE)
OUT <- normalizePath(if (length(a) >= 1) a[1] else file.path(R, "work/mcmctree_bb"))
NC <- if (length(a) >= 2) as.integer(a[2]) else 8L
ND <- if (length(a) >= 3) as.integer(a[3]) else NA_integer_
D <- file.path(OUT, "use"); O <- file.path(OUT, "interp"); dir.create(O, showWarnings = FALSE)
u <- strsplit(readLines(file.path(D, "USE"), warn = FALSE)[1], " ")[[1]]; sfx <- if (u[1] == "none") "" else u[1]
rep_ <- c(sprintf("backbone run %s, chains suffix '%s', convergence %s", basename(normalizePath(D)), sfx, u[2]))
# ---- trees ----------------------------------------------------------------------------------------------------
tt <- load_tree(ma = TRUE); n <- Ntip(tt)                       # LSD2, Ma
kids_of <- function(tr, v) tr$edge[tr$edge[, 1] == v, 2]
tips_under <- function(tr, v) if (v <= Ntip(tr)) tr$tip.label[v] else extract.clade(tr, v)$tip.label
side <- tips_under(tt, kids_of(tt, n + 1L)[1])
ml <- read.tree(file.path(W, "bl.treefile")); ml <- keep.tip(ml, tt$tip.label)
if (length(side) > 1) ml <- root(ml, outgroup = setdiff(tt$tip.label, side)[1], resolve.root = TRUE)  # side = a clade
ml <- root(ml, outgroup = side, resolve.root = TRUE)
rc <- kids_of(ml, n + 1L)
if (length(rc) != 2 || !any(vapply(rc, function(v) setequal(tips_under(ml, v), side), TRUE)))
  stop("ML tree cannot be rooted on the time tree's root split")
re <- which(ml$edge[, 1] == n + 1L); ml$edge.length[re] <- sum(ml$edge.length[re]) / 2
ml$edge.length[ml$edge.length < 1e-8] <- 1e-8
topo_ok <- isTRUE(all.equal.phylo(ml, tt, use.edge.length = FALSE))
rep_ <- c(rep_, sprintf("ML tree %d tips, rooted on the time-tree root split; topology identical to the LSD2 time tree: %s",
                        n, topo_ok))
if (!topo_ok) { writeLines(rep_, file.path(O, "report.txt")); stop("topology differs from the time tree") }
# full-tree node definitions: first tip of the first and of the second child
first_tip <- function(tr, v) { while (v > Ntip(tr)) v <- kids_of(tr, v)[1]; tr$tip.label[v] }
inn <- n + seq_len(ml$Nnode)
def <- rbindlist(lapply(inn, function(v) { k <- kids_of(ml, v)
  data.table(node = v, tipA = first_tip(ml, k[1]), tipB = first_tip(ml, k[2]), n_tips = length(tips_under(ml, v))) }))
ages_at <- function(tr, A, B) {                                  # node ages (time before present) at MRCA(A, B)
  d <- node.depth.edgelength(tr); H <- max(d[seq_len(Ntip(tr))]); M <- mrca(tr)
  H - d[M[cbind(match(A, rownames(M)), match(B, colnames(M)))]]
}
def[, lsd2_Ma := ages_at(tt, tipA, tipB)]
# ---- backbone ages ------------------------------------------------------------------------------------------------
na <- fread(file.path(D, sprintf("node_ages%s.tsv", sfx)), na.strings = "NA")
mc <- function(ck) { cs <- intersect(sprintf("%s_r%d%s_mean", ck, 1:2, sfx), names(na)); if (!length(cs)) return(NULL)
  m <- rowMeans(as.matrix(na[, ..cs]), na.rm = TRUE); if (all(!is.finite(m))) NULL else m }
bbn <- vapply(seq_len(nrow(na)), function(i) as.integer(getMRCA(ml, c(na$tipA[i], na$tipB[i]))), 1L)
stopifnot(!anyDuplicated(bbn), (n + 1L) %in% bbn)
def[, backbone := as.integer(node %in% bbn)]
dr <- fread(file.path(D, sprintf("draws%s.tsv", sfx)))
dr[, bb := bbn[match(node, na$node)]]; stopifnot(!anyNA(dr$bb))
ids <- sort(unique(dr$draw)); if (is.finite(ND)) ids <- head(ids, ND)
# ---- chronos ------------------------------------------------------------------------------------------------------
chr <- function(ages, model = "correlated") {
  cal <- data.frame(node = bbn, age.min = ages, age.max = ages, soft.bounds = FALSE)
  r <- tryCatch(suppressWarnings(chronos(ml, lambda = 1, model = model, quiet = TRUE, calibration = cal,
                                         control = chronos.control(iter.max = 1e4, eval.max = 1e4))),
                error = function(e) conditionMessage(e))
  if (is.character(r)) return(list(tree = NULL, msg = r))
  msg <- paste(attr(r, "message"), collapse = " "); pl <- attr(r, "ploglik")
  r <- list(edge = r$edge, edge.length = r$edge.length, tip.label = r$tip.label, Nnode = r$Nnode); class(r) <- "phylo"
  list(tree = r, msg = msg, ploglik = pl)
}
t0 <- Sys.time()
fits <- mclapply(ids, function(k) { x <- dr[draw == k]; chr(x$age_Ma[match(bbn, x$bb)]) }, mc.cores = NC, mc.preschedule = FALSE)
m2 <- mc("c2"); m3 <- mc("c3")
spec <- list(point_correlated = list(m2, "correlated"), point_relaxed = list(m2, "relaxed"), point_clock = list(m2, "clock"),
             c3_correlated = list(m3, "correlated"))
sens <- mclapply(names(spec), function(v) if (is.null(spec[[v]][[1]])) list(tree = NULL, msg = "no posterior means") else
  chr(spec[[v]][[1]], spec[[v]][[2]]), mc.cores = min(NC, length(spec)), mc.preschedule = FALSE)
names(sens) <- names(spec)
good <- function(f) is.list(f) && !is.null(f$tree)               # a killed worker returns NULL / try-error
msg_of <- function(f) if (is.list(f) && length(f$msg)) as.character(f$msg)[1] else "worker failed"
for (v in names(sens)) if (!good(sens[[v]])) sens[[v]] <- list(tree = NULL, msg = msg_of(sens[[v]]))
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
# ---- summaries ----------------------------------------------------------------------------------------------------
ok <- Filter(good, fits)
A <- vapply(ok, function(f) ages_at(f$tree, def$tipA, def$tipB), numeric(nrow(def)))
A <- matrix(A, nrow(def))
q_ <- function(p) if (ncol(A)) apply(A, 1, quantile, p, names = FALSE) else NA_real_
def[, `:=`(draw_med = q_(0.5), draw_lo = q_(0.025), draw_hi = q_(0.975))]
for (v in names(sens)) def[, (paste0(v, "_Ma")) := if (is.null(sens[[v]]$tree)) NA_real_ else ages_at(sens[[v]]$tree, tipA, tipB)]
fwrite(def, file.path(O, "node_ages.tsv"), sep = "\t", na = "NA", quote = FALSE)
if (length(ok)) { trs <- lapply(ok, `[[`, "tree"); class(trs) <- "multiPhylo"; write.tree(trs, file.path(O, "draws.nwk")) }
for (v in names(sens)) if (!is.null(sens[[v]]$tree)) write.tree(sens[[v]]$tree, file.path(O, paste0(v, ".nwk")))
cr <- function(x, y, k = rep(TRUE, length(x))) { k <- k & is.finite(x) & is.finite(y)
  if (sum(k) < 3) return("n < 3")
  sprintf("Pearson r %.4f, median |diff| %.2f Ma, n %d", cor(x[k], y[k]), median(abs(x[k] - y[k])), sum(k)) }
nb <- def$backbone == 0; bi <- match(bbn, def$node)
r_bb <- if (!is.null(m2)) cor(m2, def$lsd2_Ma[bi]) else NA_real_
rep_ <- c(rep_,
  sprintf("backbone nodes %d of %d (root fixed); draws %d, chronos fits OK %d (%.1f min, %d cores)", length(bbn), nrow(def),
          length(ids), length(ok), el, NC),
  if (length(ok) < length(fits)) paste("failed draws:", paste(head(unique(vapply(Filter(Negate(good), fits), msg_of, "")), 3),
                                                              collapse = " | ")),
  sprintf("PREREG_fig4 C (LSD2 main unless r < 0.95), backbone nodes = MCMCTree estimates: r(clock-2 posterior mean, LSD2) %.4f -> %s",
          r_bb, if (!is.finite(r_bb)) "NA" else if (r_bb >= 0.95) "LSD2 stays main" else "r < 0.95: MCMCTree ages main"),
  paste("point (correlated) vs LSD2, all nodes:", cr(def$point_correlated_Ma, def$lsd2_Ma)),
  paste("point (correlated) vs LSD2, interpolated nodes:", cr(def$point_correlated_Ma, def$lsd2_Ma, nb)),
  sprintf("root: point %.1f Ma, LSD2 %.1f Ma, draws %.1f (%.1f-%.1f)", def$point_correlated_Ma[1], def$lsd2_Ma[1],
          def$draw_med[1], def$draw_lo[1], def$draw_hi[1]),
  if (!is.null(m2) && !is.null(sens$point_correlated$tree))
    sprintf("fixed backbone ages reproduced by the point tree: max |diff| %.2e Ma", max(abs(def$point_correlated_Ma[bi] - m2))),
  sprintf("draw 95 %% interval width, median: interpolated nodes %.1f Ma, backbone %.1f Ma; LSD2 age inside: all %.3f, interpolated %.3f",
          median((def$draw_hi - def$draw_lo)[nb]), median((def$draw_hi - def$draw_lo)[!nb]),
          mean(def$lsd2_Ma >= def$draw_lo & def$lsd2_Ma <= def$draw_hi), mean((def$lsd2_Ma >= def$draw_lo & def$lsd2_Ma <= def$draw_hi)[nb])),
  vapply(setdiff(names(sens), "point_correlated"), function(v) paste0(v, " vs point (all nodes): ",
         cr(def[[paste0(v, "_Ma")]], def$point_correlated_Ma)), ""),
  vapply(names(sens), function(v) sprintf("%s chronos: %s", v, substr(msg_of(sens[[v]]), 1, 120)), ""))
writeLines(rep_, file.path(O, "report.txt")); cat(rep_, sep = "\n")
