#!/usr/bin/env Rscript
# Phylogenetic regression of proteome / CDS physicochemistry on deep-sea habitat.
# Usage: Rscript 32_pgls_physchem.R TREE OUTDIR
# Species level: trait ~ deep (+ log10 depth sensitivity) with phylolm (Pagel's lambda), across all
# Animalia and per clade (>= 6 deep and >= 6 shallow). Annotation route (source) as covariate.
# Gene level: per BUSCO locus present in >= 20 species, phylolm trait ~ deep; BH across loci;
# reports the fraction of loci shifting in the species-level direction (sign test).
suppressMessages({library(ape); library(phylolm); library(data.table)})
a <- commandArgs(TRUE); tr <- read.tree(a[1]); out <- a[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
R <- "~/deepsea_conv"
sp <- fread(file.path(R, "work/physchem/species_metrics.tsv"))
lb <- fread(file.path(R, "meta/labels_v2.tsv"), sep = "\t", quote = "")
lc <- Sys.getenv("LABEL_COL", "lab")   # lab (default) / lab_strict / lab_1000
sp <- merge(sp, lb[, .(code, L = get(lc))], by = "code")
sp <- sp[L %in% c("D", "S")]
sp[, deep := as.integer(L == "D")]
cat("labels:", lc, "D", sum(sp$deep), "S", sum(1 - sp$deep), "\n")
sp[, ldepth := log10(pmax(as.numeric(depth_median), 1))]
sp[, route := ifelse(source %in% c("ncbi_annot", "local", "extra_local"), "annot", source)]
sp[, clade := fcase(class == "Actinopteri", "Actinopteri",
                    class %in% c("Chondrichthyes", "Myxini", "Hyperoartia"), "Chondrichthyes",
                    phylum %in% c("Annelida", "Nemertea"), "Annelida",
                    default = phylum)]
sp <- sp[code %in% tr$tip.label]
tr <- keep.tip(tr, sp$code)
if (!is.rooted(tr)) tr <- phytools::midpoint.root(tr)
tr$edge.length[tr$edge.length <= 0] <- 1e-6
dated <- Sys.getenv("DATED") == "1" || is.ultrametric(tr, tol = 1e-4, option = 1)
if (!dated) {   # relative-time ultrametric when no dated tree is given
  tr <- chronos(tr, lambda = 1, quiet = TRUE); class(tr) <- "phylo"
}
cat("tree:", if (dated) "dated tree as given" else "chronos lambda = 1 (relative time)", "; tips", Ntip(tr), "\n")
traits <- grep("^(busco|whole)_", names(sp), value = TRUE)
traits <- setdiff(traits, c("busco_len", "whole_len"))
fit1 <- function(d, tt, f) {
  df <- as.data.frame(d); rownames(df) <- df$code
  t2 <- keep.tip(tt, df$code)
  m <- tryCatch(phylolm(f, data = df, phy = t2, model = "lambda"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- summary(m)$coefficients
  data.table(term = rownames(cf), est = cf[, 1], se = cf[, 2], p = cf[, 4], lambda = m$optpar, n = nrow(df))
}
res <- list()
for (tv in traits) {
  d <- sp[is.finite(get(tv))]
  d[, y := scale(get(tv))[, 1]]   # standardised effect (SD units)
  r <- fit1(d, tr, y ~ deep + route); if (!is.null(r)) res[[length(res) + 1]] <- r[term == "deep"][, `:=`(trait = tv, scope = "Animalia", model = "deep+route")]
  r <- fit1(d[is.finite(ldepth)], tr, y ~ ldepth + route); if (!is.null(r)) res[[length(res) + 1]] <- r[term == "ldepth"][, `:=`(trait = tv, scope = "Animalia", model = "log10depth+route")]
  for (cl in unique(d$clade)) {
    dc <- d[clade == cl]
    if (sum(dc$deep) < 6 || sum(1 - dc$deep) < 6) next
    f <- if (length(unique(dc$route)) > 1) y ~ deep + route else y ~ deep
    r <- fit1(dc, tr, f); if (!is.null(r)) res[[length(res) + 1]] <- r[term == "deep"][, `:=`(trait = tv, scope = cl, model = "deep")]
  }
}
res <- rbindlist(res)
res[, q := p.adjust(p, "BH"), by = .(scope, model)]
fwrite(res[order(scope, model, p)], file.path(out, "pgls_species.tsv"), sep = "\t")
print(res[scope == "Animalia" & model == "deep+route"][order(p)][, .(trait, est, se, p, q, lambda, n)])

# gene level for the traits significant at species level (or the core set)
core <- c("acid_base", "net_charge7_per_res", "pI", "R_RK", "f_C", "f_P", "f_G", "gravy", "GC3", "GC4", "IVYWREL")
g <- fread(cmd = paste("zcat", file.path(R, "work/physchem/gene_metrics.tsv.gz")))
g <- g[code %in% sp$code]
g <- merge(g, sp[, .(code, deep, route)], by = "code")
gene_fit <- function(b) {
  gb <- g[busco == b]; o <- list()
  if (nrow(gb) < 20 || sum(gb$deep) < 5 || sum(1 - gb$deep) < 5) return(NULL)
  for (tv in intersect(core, names(gb))) {
    d <- gb[is.finite(get(tv))]; if (nrow(d) < 20) next
    d[, y := scale(get(tv))[, 1]]
    if (!is.finite(sd(d$y))) next
    r <- fit1(d, tr, y ~ deep); if (!is.null(r)) o[[length(o) + 1]] <- r[term == "deep"][, `:=`(busco = b, trait = tv)]
  }
  rbindlist(o)
}
gl <- parallel::mclapply(unique(g$busco), gene_fit, mc.cores = as.integer(Sys.getenv("NCORES", "16")))
if (!length(Filter(function(z) is.data.table(z) && nrow(z), gl))) { cat("no loci pass gene-level filters
"); quit(save = "no") }
gl <- rbindlist(Filter(is.data.table, gl)); gl[, q := p.adjust(p, "BH"), by = trait]
fwrite(gl, file.path(out, "pgls_gene.tsv"), sep = "\t")
sm <- gl[, .(n_loci = .N, frac_pos = mean(est > 0), sign_p = binom.test(sum(est > 0), .N)$p.value,
             n_q05 = sum(q < 0.05), median_est = median(est)), by = trait]
fwrite(sm, file.path(out, "pgls_gene_summary.tsv"), sep = "\t"); print(sm)
