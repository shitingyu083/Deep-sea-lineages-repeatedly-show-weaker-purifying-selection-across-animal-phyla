#!/usr/bin/env Rscript
# Robustness for 66c (not preregistered; DEVIATIONS D-R-4): permulation-calibrated GO rank test for the fish RER
# (Chikina et al. 2016 / Saputra et al. 2021 logic). For each of the 1,000 binary permulations in perms_binary.rds the
# gene statistic s = sign(Rho) * -log10(P) is recomputed, term AUCs (Mann-Whitney U / (k(N-k)), genes with GO only)
# recomputed, and the empirical two-sided p = (1 + #|AUC_perm - .5| >= |AUC_obs - .5|) / 1001. BH within namespace.
# Read-only on inputs. Usage: Rscript 66f_fishrer_go_perm_sty.R FUNC_DIR
a <- commandArgs(TRUE); FUNC <- a[1]
B <- path.expand("~/deepsea_conv"); O <- file.path(B, "work/addon/R3/fishrer/go")
ann <- file.path(FUNC, "Actinopteri_m/hog_annot.tsv")
if (!file.exists(ann)) ann <- file.path(FUNC, "Actinopteri/hog_annot.tsv")
an <- read.delim(ann, colClasses = "character")
obs <- read.delim(file.path(B, "work/addon/R3/fishrer/lab/binary.tsv"), colClasses = c(gene = "character"))
obs <- obs[!is.na(obs$P) & !is.na(obs$Rho), ]
go <- strsplit(an$go, ","); names(go) <- an$hog
genes <- intersect(obs$gene, an$hog[nchar(an$go) > 0])
pr <- readRDS(file.path(B, "work/addon/R3/fishrer/lab/perms_binary.rds"))
genes <- intersect(genes, rownames(pr$corP))
cat("genes with GO and perms", length(genes), "\n")
terms0 <- read.delim(file.path(O, "fishrer_go_rank.tsv"), colClasses = c(term = "character"))
mem <- stack(lapply(go[genes], function(x) x[nchar(x) > 0])); colnames(mem) <- c("term", "gene")
mem <- mem[mem$term %in% terms0$term, ]
tl <- terms0$term; gi <- match(as.character(mem$gene), genes); ti <- match(mem$term, tl)
M <- matrix(0, length(tl), length(genes)); M[cbind(ti, gi)] <- 1
sfun <- function(p, r) sign(r) * -log10(pmax(p, 1e-300))
auc_mat <- function(S) {  # S genes x reps, NA allowed
  ok <- !is.na(S); R <- apply(S, 2, function(x) { r <- rank(x, na.last = "keep"); r[is.na(r)] <- 0; r })
  N <- colSums(ok); K <- M %*% ok; RS <- M %*% R
  U <- RS - K * (K + 1) / 2; U / (K * (rep(N, each = nrow(K)) - K))
}
so <- sfun(obs$P[match(genes, obs$gene)], obs$Rho[match(genes, obs$gene)])
A0 <- as.vector(auc_mat(matrix(so, ncol = 1)))
SP <- sfun(pr$corP[genes, , drop = FALSE], pr$corRho[genes, , drop = FALSE]); rm(pr); gc()
AP <- auc_mat(SP)
ex <- rowSums(abs(AP - 0.5) >= abs(A0 - 0.5), na.rm = TRUE); np <- rowSums(!is.na(AP))
pe <- (1 + ex) / (1 + np)
z <- (A0 - rowMeans(AP, na.rm = TRUE)) / apply(AP, 1, sd, na.rm = TRUE)
res <- data.frame(term = tl, namespace = terms0$namespace, name = terms0$name, k = terms0$k,
                  AUC = A0, AUC_param = terms0$AUC, p_param = terms0$p, q_param = terms0$q,
                  perm_mean = rowMeans(AP, na.rm = TRUE), perm_sd = apply(AP, 1, sd, na.rm = TRUE),
                  z_perm = z, n_perm = np, p_perm = pe, q_perm = NA)
for (ns in unique(res$namespace)) { i <- res$namespace == ns; res$q_perm[i] <- p.adjust(res$p_perm[i], "BH") }
res <- res[order(res$p_perm, -abs(res$z_perm)), ]
write.table(res, file.path(O, "fishrer_go_perm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("max |AUC - AUC_param|", max(abs(res$AUC - res$AUC_param), na.rm = TRUE), "\n")
cat("terms", nrow(res), "q_param<.05", sum(res$q_param < .05), "q_perm<.05", sum(res$q_perm < .05),
    "q_perm<.05 & AUC>.5", sum(res$q_perm < .05 & res$AUC > .5), "\n")
print(head(res[, c("term", "namespace", "name", "k", "AUC", "z_perm", "p_perm", "q_perm", "q_param")], 25), row.names = FALSE)
cat("Wang 2019 terms:\n")
print(res[res$term %in% c("GO:0031214", "GO:0001503", "GO:0006643", "GO:0006457"),
          c("term", "name", "k", "AUC", "z_perm", "p_perm", "q_perm")], row.names = FALSE)
file.create(file.path(O, "GOPERM.DONE"))
