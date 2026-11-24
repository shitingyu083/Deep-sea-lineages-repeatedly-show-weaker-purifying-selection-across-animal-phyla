#!/usr/bin/env Rscript
# Relative evolutionary rates (RERconverge 0.3.0) of the species-tree loci and their association with deep-sea habitat.
# Usage: Rscript 37_rerconverge.R TREES MASTER OUTDIR [LABEL_COL=lab] [NPERM=1000] [NCORES=16]
#  TREES   "locus<TAB>newick" gene trees, branch lengths on the fixed species topology (36_rer_trees.sh)
#  MASTER  rooted species tree with branch lengths (the dated tree), tips = species codes
#  Species: D (deep) and S (shallow-marine background) of meta/labels_v2.tsv column LABEL_COL (lab / lab_strict /
#   lab_1000); all other species are pruned from the gene trees and the master tree.
#  Gene trees in the master's form, as estimatePhangornTreeAll builds them: the master pruned to the gene's species with
#   the gene's branch lengths, transferred by bipartition (any unmatched bipartition = topology differs: error); read
#   with the unrooted master (readTrees; the master's own branch lengths are kept: minSpecs unset),
#   RER = getAllResiduals(transform = "sqrt", min.sp = 10; package defaults otherwise), run step by step without its
#   diagnostic plot (identical values; see below).
#  Binary: foreground = D species, foreground2Paths(clade = "all"), Kendall tau per gene (correlateWithBinaryPhenotype,
#   min.pos = 2); permulation p-values (permpvalcor) from NPERM categorical permulations (getPermsBinary permmode = "cc":
#   ER Mk model fitted on the master tree, tip-state counts preserved). clade = "terminal": observed-only sensitivity.
#  Continuous: log10(max(depth_median, 1 m)) of the same species; char2Paths + winsorised Pearson (correlateWithContinuousPhenotype);
#   NPERM phylogenetic permulations (getPermsContinuous type = "simperm") on the rooted master tree.
#  Genome-wide: median Rho over genes and the numbers of genes with P < 0.05 and Rho > 0 / Rho < 0, observed vs. the
#   same statistics of every permulated phenotype (empirical p, add-one).
#  Permulations run in NCORES chunks (L'Ecuyer-CMRG streams, chunk seeds 20261019 + i).
#  env EXCLUDE_CHEMO=1 (sensitivity, 46_rer_nochemo.sh): the chemo_broad species of meta/chemo_habitat.tsv
#   (45_chemo_flag.py: vent / seep / fall endemics, chemosymbiotic hosts, vent-collected specimens) are pruned too.
# Outputs in OUTDIR: binary.tsv, binary_terminal.tsv, continuous.tsv, genomewide.tsv, gene_trees.txt, perms_binary.rds,
#   perms_continuous.rds, rer.rds, trees.rds (readTrees object), rer_diagnostics.pdf, run.txt                DONE
suppressMessages({library(RERconverge); library(ape); library(data.table); library(parallel)})
a <- commandArgs(TRUE)
tf <- a[1]; mf <- a[2]; out <- a[3]
lc <- if (length(a) > 3) a[4] else "lab"
nperm <- if (length(a) > 4) as.integer(a[5]) else 1000L
nc <- if (length(a) > 5) as.integer(a[6]) else 16L
dir.create(out, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
say <- function(...) { cat(format(Sys.time(), "%F %T"), "[37]", ..., "\n"); flush.console() }
run <- file.path(out, "run.txt")
note <- function(...) cat(..., "\n", file = run, append = TRUE, sep = "")
cat("", file = run)
note("RERconverge ", as.character(packageVersion("RERconverge")), "; R ", R.version$major, ".", R.version$minor)
note("trees\t", normalizePath(tf), "\nmaster\t", normalizePath(mf), "\nlabels\t", lc, "\nnperm\t", nperm)

lb <- fread("~/deepsea_conv/meta/labels_v2.tsv", sep = "\t", quote = "")
mt <- read.tree(mf)
if (!is.rooted(mt)) stop("MASTER must be rooted")
lb <- lb[code %in% mt$tip.label & get(lc) %in% c("D", "S")]
if (Sys.getenv("EXCLUDE_CHEMO") == "1") {
  chf <- path.expand("~/deepsea_conv/meta/chemo_habitat.tsv")
  ch <- fread(chf, sep = "\t", quote = "")
  dropc <- intersect(lb$code, ch[chemo_broad == 1, code])
  say("EXCLUDE_CHEMO: pruning", length(dropc), "chemo_broad species (D", lb[code %in% dropc & get(lc) == "D", .N], ")")
  note("exclude_chemo\t", length(dropc), "\t", paste(sort(dropc), collapse = ","), "\nchemo_md5\t", unname(tools::md5sum(chf)))
  lb <- lb[!code %in% dropc]
}
gl <- strsplit(readLines(tf), "\t", fixed = TRUE)
gt <- setNames(lapply(gl, function(x) read.tree(text = x[2])), vapply(gl, `[`, "", 1))
gsp <- unique(unlist(lapply(gt, `[[`, "tip.label")))
if (length(nog <- setdiff(lb$code, gsp))) {      # e.g. taxa excluded from the gene trees (meta/topo_fix.tsv drop)
  say("species in no gene tree, dropped:", paste(nog, collapse = ","))
  note("dropped_no_gene_tree\t", paste(nog, collapse = ","))
  lb <- lb[code %in% gsp]
}
sp <- lb$code
fg <- lb[get(lc) == "D", code]
mr <- multi2di(keep.tip(mt, sp))                 # rooted, binary: permulation simulations
mr$edge.length[mr$edge.length <= 0] <- 1e-6
mu <- unroot(mr)                                 # readTrees master (gene trees are read unrooted)
say("labels", lc, ": D", length(fg), "S", length(sp) - length(fg), "; master tips", Ntip(mu))
note("species\t", length(sp), " (D ", length(fg), ", S ", length(sp) - length(fg), ")")

# Gene trees in the master's form. readTrees maps the written ancestor-descendant node pairs of each gene tree onto the
# master (allPathsMasterRelativeTT); a tree written with another root (IQ-TREE roots at the first sequence) is
# "discordant" and gets an all-zero path row even when its unrooted topology is the master's.
spc <- setNames(seq_along(sp), sp)
bkey <- function(phy) {              # bipartition of every edge: the side without the lowest-coded species, sorted codes
  tc <- spc[phy$tip.label]; r0 <- which.min(tc); n <- length(tc)
  vapply(phangorn::Descendants(phy, phy$edge[, 2], "tips"),
         function(s) { if (r0 %in% s) s <- seq_len(n)[-s]; paste(sort(tc[s]), collapse = ",") }, "")
}
al <- lapply(gt, function(g) {
  s <- intersect(g$tip.label, sp)
  if (length(s) < 3) return(NULL)
  g <- unroot(keep.tip(g, s)); p <- unroot(keep.tip(mu, s))
  i <- match(bkey(p), bkey(g))
  if (anyNA(i) || anyDuplicated(i) || nrow(p$edge) != nrow(g$edge)) return(NA)   # topology differs from the master
  p$edge.length <- g$edge.length[i]
  p
})
bad <- names(al)[vapply(al, function(x) identical(x, NA), TRUE)]
if (length(bad)) stop(length(bad), " gene trees differ in topology from the master, e.g. ", bad[1])
al <- al[!vapply(al, is.null, TRUE)]
gtf <- file.path(out, "gene_trees.txt")
writeLines(paste(names(al), vapply(al, write.tree, "", digits = 10), sep = "\t"), gtf)
# report recomputed over all master tips (readTrees builds it from the species of the trees that set a new maximum only)
tr <- suppressMessages(readTrees(gtf, masterTree = mu))
mtl <- tr$masterTree$tip.label
tr$report <- do.call(rbind, lapply(tr$trees, function(x) as.integer(mtl %in% x$tip.label)))
dimnames(tr$report) <- list(names(tr$trees), mtl)
tr$maxSp <- length(mtl)
stopifnot(isTRUE(all.equal(sum(tr$masterTree$edge.length), sum(mu$edge.length))))
nsp <- rowSums(tr$report)
say("readTrees:", tr$numTrees, "gene trees; species per tree median", median(nsp), "range", min(nsp), "-", max(nsp))
z <- rowSums(tr$paths != 0, na.rm = TRUE) == 0     # readTrees' "discordant topology" rows (all zero, not NA)
if (any(z)) stop(sum(z), " gene trees discordant with the master in readTrees, e.g. ", rownames(tr$report)[which(z)[1]])
note("gene_trees\t", tr$numTrees, "\nspecies_per_tree_median\t", median(nsp))

saveRDS(tr, file.path(out, "trees.rds"))

# RER = getAllResiduals(tr, transform = "sqrt", min.sp = 10), run as its own steps (RERconverge 0.3.0: transformPaths,
# coreGetResiduals, getRMat with their defaults) with computeWeightsAllVar(plot = FALSE): the plot's "After" panel refits
# every gene on all of its ~1e4 ancestor-descendant paths with the dense O(n^2) weighted solver (~10 s per gene here,
# hours per run) and leaves the weights unchanged. transformPaths' kNN-imputed paths are not computed: they are used only
# with n.pcs > 0 (n.pcs = 0).
ns <- asNamespace("RERconverge")
tp <- tr
nv <- sqrt(ns$getColMeansNV(tp$paths))
tp$paths <- sqrt(tp$paths); tp$transform <- "sqrt"; tp$transformFunc <- sqrt; tp$transformInv <- function(x) x^2
tp$weights <- ns$computeWeightsAllVar(tp$paths, nv = nv, transform = "none", plot = FALSE)
rer <- ns$getRMat(ns$coreGetResiduals(tp, min.sp = 10))
saveRDS(rer, file.path(out, "rer.rds"))
say("RER matrix", nrow(rer), "genes x", ncol(rer), "paths;", sum(is.finite(rer)), "values")
# diagnostic: the weighted residuals should show no trend of variance with (sqrt) path length
ok <- which(is.finite(rer) & rer != 0 & is.finite(tp$paths))
pdf(file.path(out, "rer_diagnostics.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
ns$plotAsBox(tp$paths[ok], log(rer[ok]^2), n = 50, main = "log RER^2 by sqrt path length (bins)", xlab = "")
hist(rer[ok], 200, main = "RER", xlab = "RER")
invisible(dev.off())
rm(tp)

chunks <- split(seq_len(nperm), cut(seq_len(nperm), min(nc, nperm), labels = FALSE))
RNGkind("L'Ecuyer-CMRG")
perms <- function(f, seed0) {
  res <- mclapply(seq_along(chunks), function(i) {
    set.seed(seed0 + i)
    r <- NULL
    invisible(capture.output(r <- suppressMessages(f(length(chunks[[i]])))))
    lapply(r[c("corP", "corRho")], function(m) { m <- as.matrix(m); colnames(m) <- paste0("perm", chunks[[i]]); m })
  }, mc.cores = nc, mc.preschedule = FALSE)
  bad <- which(vapply(res, function(x) inherits(x, "try-error") || is.null(x), TRUE))
  if (length(bad)) stop("permulation chunk(s) failed: ", paste(bad, collapse = ","), "\n", as.character(res[[bad[1]]]))
  list(corP = do.call(cbind, lapply(res, `[[`, "corP")), corRho = do.call(cbind, lapply(res, `[[`, "corRho")))
}
gw <- function(obs, pm, test) {        # genome-wide statistics: observed vs. permulated
  st <- function(rho, p) c(median_rho = median(rho, na.rm = TRUE),
                           n_pos_p05 = sum(p < 0.05 & rho > 0, na.rm = TRUE),
                           n_neg_p05 = sum(p < 0.05 & rho < 0, na.rm = TRUE))
  o <- st(obs$Rho, obs$P)
  nl <- sapply(seq_len(ncol(pm$corRho)), function(j) st(pm$corRho[, j], pm$corP[, j]))
  rbindlist(lapply(names(o), function(k) {
    x <- nl[k, ]; n <- sum(is.finite(x))
    data.table(test = test, statistic = k, observed = o[[k]], null_mean = mean(x, na.rm = TRUE), null_sd = sd(x, na.rm = TRUE),
               null_q025 = quantile(x, 0.025, na.rm = TRUE), null_q975 = quantile(x, 0.975, na.rm = TRUE),
               p_upper = (1 + sum(x >= o[[k]], na.rm = TRUE)) / (1 + n),
               p_lower = (1 + sum(x <= o[[k]], na.rm = TRUE)) / (1 + n), n_perm = n)
  }))
}
tab <- function(cr, pp = NULL) {
  d <- data.table(gene = rownames(cr), Rho = cr$Rho, N = cr$N, P = cr$P, p.adj = cr$p.adj)
  if (!is.null(pp)) d[, `:=`(permP = pp$permpval[match(gene, rownames(pp))], permStat = pp$permstats[match(gene, rownames(pp))])][
    , permQ := p.adjust(permP, "BH")]
  d
}

# binary: deep foreground
phen <- foreground2Paths(fg, tr, clade = "all")
cb <- correlateWithBinaryPhenotype(rer, phen, min.sp = 10, min.pos = 2)
cbt <- correlateWithBinaryPhenotype(rer, foreground2Paths(fg, tr, clade = "terminal"), min.sp = 10, min.pos = 2)
fwrite(tab(cbt), file.path(out, "binary_terminal.tsv"), sep = "\t")
say("binary observed:", sum(is.finite(cb$Rho)), "genes tested; median Rho", signif(median(cb$Rho, na.rm = TRUE), 3))
pb <- perms(function(k) getPermsBinary(k, fg, NA, fg[1], rer, tr, mr, permmode = "cc", calculateenrich = FALSE), 20261019L)
saveRDS(pb, file.path(out, "perms_binary.rds"))
fwrite(tab(cb, permpvalcor(cb, pb)), file.path(out, "binary.tsv"), sep = "\t")
g1 <- gw(cb, pb, "binary_deep")
say("binary permulations done:", ncol(pb$corRho))

# continuous: log10 median depth (metres, floored at 1 m as in 32_pgls_physchem.R)
dm <- setNames(log10(pmax(as.numeric(lb$depth_median), 1)), lb$code)
dm <- dm[is.finite(dm)]
mrc <- keep.tip(mr, names(dm))
cp <- char2Paths(dm, tr)
cc <- correlateWithContinuousPhenotype(rer, cp, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
say("continuous observed:", length(dm), "species,", sum(is.finite(cc$Rho)), "genes tested")
note("continuous_species\t", length(dm))
pc <- perms(function(k) getPermsContinuous(k, dm, rer, NULL, tr, mrc, calculateenrich = FALSE, type = "simperm"),
            20261019L + 100000L)
saveRDS(pc, file.path(out, "perms_continuous.rds"))
fwrite(tab(cc, permpvalcor(cc, pc)), file.path(out, "continuous.tsv"), sep = "\t")
g2 <- gw(cc, pc, "continuous_log10depth")
fwrite(rbind(g1, g2), file.path(out, "genomewide.tsv"), sep = "\t")
print(rbind(g1, g2))
note("runtime_min\t", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
file.create(file.path(out, "DONE"))
say("DONE")
