# Extra R packages for the comparative/figure layer (2026-10-20). Installs only missing packages,
# never upgrades installed ones (RERconverge stack must stay as is).
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
lib <- .libPaths()[1]
cran <- c("ggrepel", "corHMM", "ggnewscale", "metafor", "ggridges", "phylopath", "BiocManager")
miss <- cran[!cran %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss, lib = lib, dependencies = c("Depends", "Imports", "LinkingTo"))
if (!"fgsea" %in% rownames(installed.packages()))
  BiocManager::install("fgsea", update = FALSE, ask = FALSE, lib = lib)
all <- c(cran, "fgsea")
ok <- sapply(all, function(p) suppressWarnings(requireNamespace(p, quietly = TRUE)))
print(ok)
cat(if (all(ok)) "RPKGS_OK\n" else "RPKGS_INCOMPLETE\n")
