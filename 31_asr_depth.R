#!/usr/bin/env Rscript
# Ancestral depth-habitat reconstruction and count of independent shallow->deep transitions.
# Usage: Rscript 31_asr_depth.R TREE OUTDIR [NSIM=1000]
#  TREE   rooted (or rooted here on Ctenophora) species tree, tips = species codes
#  states from meta/labels_v2.tsv (34a_labels.py), column ASR_COL (default 'asr'): D = deep,
#  S = shallow-marine background + non-marine species; anything else (eurybathic, no depth data,
#  anchialine, parasitic) -> pruned. ASR_COL=asr_strict / asr_1000 for the sensitivity runs.
# Outputs: OUTDIR/asr_fit.tsv (ER/ARD/irreversible AIC), asr_nodes.tsv (P(deep) per node),
#          transitions.tsv (simmap counts S->D, D->S per simulation), deep_origins.tsv
#          (parsimony-independent origins: maximal deep clades with their crown/stem node).
suppressMessages({library(ape); library(phytools); library(data.table)})
a <- commandArgs(TRUE)
tr <- read.tree(a[1]); out <- a[2]; nsim <- if (length(a) > 2) as.integer(a[3]) else 1000
dir.create(out, showWarnings = FALSE, recursive = TRUE)
lb <- fread("~/deepsea_conv/meta/labels_v2.tsv", sep = "\t", quote = "")
col <- Sys.getenv("ASR_COL", "asr")
st <- lb[[col]]; st[!st %in% c("D", "S")] <- NA; names(st) <- lb$code
cat("labels: meta/labels_v2.tsv column", col, "\n")
if (!is.rooted(tr)) {   # root before pruning so outgroups without depth labels can still root
  og <- intersect(c("Mnelei"), tr$tip.label)   # ctenophore-sister rooting (see Methods)
  og_env <- Sys.getenv("OUTGROUP"); if (nzchar(og_env)) og <- intersect(strsplit(og_env, ",")[[1]], tr$tip.label)
  tr <- if (length(og)) root(tr, outgroup = og, resolve.root = TRUE) else midpoint.root(tr)
}
keep <- intersect(tr$tip.label, names(st)[!is.na(st)])
tr <- drop.tip(tr, setdiff(tr$tip.label, keep))
if (is.null(tr$edge.length)) tr <- compute.brlen(tr)
tr$edge.length[tr$edge.length <= 0] <- 1e-6
x <- setNames(st[tr$tip.label], tr$tip.label)
cat("tips", length(x), "deep", sum(x == "D"), "shallow", sum(x == "S"), "\n")

# phytools' internal ARD optimiser can land below the nested ER optimum (seen on Annelida: ARD -10.85 < ER -10.41).
# Own ML: optimise log-rates of the free Q cells, scoring with fitMk(fixedQ=), from several starts; keep the best.
fitMk <- function(tree, x, model, ...) {
  f0 <- phytools::fitMk(tree, x, model = model, ...)
  idx <- f0$index.matrix; k <- max(idx, na.rm = TRUE)
  mkQ <- function(lr) { Q <- matrix(0, nrow(idx), ncol(idx), dimnames = dimnames(idx))
    Q[!is.na(idx) & idx > 0] <- exp(lr)[idx[!is.na(idx) & idx > 0]]; diag(Q) <- -rowSums(Q); Q }
  nll <- function(lr) -as.numeric(logLik(phytools::fitMk(tree, x, model = model, fixedQ = mkQ(lr), ...)))
  er <- phytools::fitMk(tree, x, model = "ER", ...)$rates
  starts <- list(log(pmax(f0$rates, 1e-4)), rep(log(er), k), rep(log(er) - 1, k), rep(log(er) + 1, k))
  best <- f0; bl <- as.numeric(logLik(f0))
  for (s in starts) {
    o <- tryCatch(optim(s, nll, method = if (k == 1) "BFGS" else "Nelder-Mead"), error = function(e) NULL)
    if (!is.null(o) && -o$value > bl + 1e-6) {
      best <- phytools::fitMk(tree, x, model = model, fixedQ = mkQ(o$par), ...); bl <- -o$value }
  }
  best$k_free <- k; best
}
set.seed(7)
fits <- list(
  ER  = fitMk(tr, x, model = "ER"),
  ARD = fitMk(tr, x, model = "ARD"),
  # rows = from, cols = to (states sorted D, S): only shallow->deep, only deep->shallow
  IRR_SD = fitMk(tr, x, model = matrix(c(0, 0, 1, 0), 2, byrow = TRUE, dimnames = list(c("D", "S"), c("D", "S")))),
  IRR_DS = fitMk(tr, x, model = matrix(c(0, 1, 0, 0), 2, byrow = TRUE, dimnames = list(c("D", "S"), c("D", "S"))))
)
aic <- sapply(fits, function(f) 2 * f$k_free - 2 * as.numeric(logLik(f)))  # df = free rates (fixedQ refits report wrong df)
fwrite(data.table(model = names(fits), logLik = sapply(fits, logLik), AIC = aic,
                  dAIC = aic - min(aic)), file.path(out, "asr_fit.tsv"), sep = "\t")
best <- names(which.min(aic)); cat("best model", best, "\n")
Q <- if (startsWith(best, "IRR")) "ARD" else best  # simmap on ER/ARD; irreversible fits reported only
anc <- ancr(fits[[best]])                         # marginal ancestral states (phytools >= 2.0)
fwrite(data.table(node = rownames(anc$ace), P_deep = anc$ace[, "D"]), file.path(out, "asr_nodes.tsv"), sep = "\t")

set.seed(1)
sm <- make.simmap(tr, x, model = Q, nsim = nsim, pi = "fitzjohn", message = FALSE)
cs <- countSimmap(sm)$Tr
cs <- as.data.table(cs); fwrite(cs, file.path(out, "transitions.tsv"), sep = "\t")
sd_col <- grep("S,D", names(cs), value = TRUE); ds_col <- grep("D,S", names(cs), value = TRUE)
cat(sprintf("S->D transitions: median %g (95%% %g-%g); D->S: median %g (95%% %g-%g)\n",
            median(cs[[sd_col]]), quantile(cs[[sd_col]], .025), quantile(cs[[sd_col]], .975),
            median(cs[[ds_col]]), quantile(cs[[ds_col]], .025), quantile(cs[[ds_col]], .975)))

# maximal all-deep clades = independent deep lineages under parsimony (as in 23_label_trees.py)
nt <- Ntip(tr); desc <- lapply(seq_len(nt + tr$Nnode), function(n) if (n <= nt) n else getDescendants(tr, n)[getDescendants(tr, n) <= nt])
alld <- sapply(desc, function(d) all(x[d] == "D"))
par <- setNames(tr$edge[, 1], tr$edge[, 2])
maxd <- which(alld & !sapply(seq_along(alld), function(n) { p <- par[as.character(n)]; !is.na(p) && alld[p] }))
fwrite(data.table(node = maxd, n_tips = sapply(desc[maxd], length),
                  tips = sapply(desc[maxd], function(d) paste(tr$tip.label[d], collapse = ","))),
       file.path(out, "deep_origins.tsv"), sep = "\t")
cat("maximal deep clades", length(maxd), "\n")
saveRDS(list(tree = tr, x = x, fits = fits, simmap_counts = cs), file.path(out, "asr.rds"))
