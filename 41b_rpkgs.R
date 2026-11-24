# Second pass of the extra R packages (2026-10-20), fixing the failures of 41_rpkgs.R:
#   GenSA >= the current CRAN release uses Rf_allocLang (R >= 4.4.1): newest archived release that builds on R 4.3.3;
#   Rmpfr needs the MPFR headers: conda-forge mpfr is added to envs/rphy before this script runs;
#   fgsea 1.28 compiles with CXX_STD = CXX11 but BH 1.90 needs C++14: CXX11STD overridden through
#     R_MAKEVARS_USER for this process only (no ~/.R/Makevars change);
#   ggm (phylopath) needs Bioconductor 'graph'.
# Installs only missing packages and never upgrades installed ones (RERconverge stack must stay as is).
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
lib <- .libPaths()[1]
have <- function(p) p %in% rownames(installed.packages(lib.loc = lib))
mv <- file.path(getwd(), "work", "rpkgs", "Makevars_cxx14")
dir.create(dirname(mv), showWarnings = FALSE, recursive = TRUE)
writeLines(c("CXX11STD = -std=gnu++14"), mv)
Sys.setenv(R_MAKEVARS_USER = mv)

bioc <- c("graph", "qvalue", "fgsea")
for (p in bioc) if (!have(p)) BiocManager::install(p, update = FALSE, ask = FALSE, lib = lib)

if (!have("GenSA")) {
  arch <- "https://cloud.r-project.org/src/contrib/Archive/GenSA/"
  idx <- readLines(arch, warn = FALSE)
  m <- regmatches(idx, regexpr("GenSA_[0-9][0-9.-]*\\.tar\\.gz", idx))
  m <- unique(m)
  v <- sub("^GenSA_(.*)\\.tar\\.gz$", "\\1", m)
  m <- m[order(numeric_version(v), decreasing = TRUE)]
  for (f in m) {
    cat("trying", f, "\n")
    try(install.packages(paste0(arch, f), repos = NULL, type = "source", lib = lib))
    if (have("GenSA")) { cat("GenSA installed from", f, "\n"); break }
  }
}
pre <- normalizePath(file.path(R.home(), "..", ".."))   # envs/rphy (not $CONDA_PREFIX: may be the base env)
if (!have("Rmpfr"))
  install.packages("Rmpfr", lib = lib, configure.args = c(
    paste0("--with-mpfr-include=", file.path(pre, "include")), paste0("--with-mpfr-lib=", file.path(pre, "lib"))))
for (p in c("corHMM", "ggm", "phylopath"))
  if (!have(p)) install.packages(p, lib = lib, dependencies = c("Depends", "Imports", "LinkingTo"))

all <- c("ggrepel", "corHMM", "ggnewscale", "metafor", "ggridges", "phylopath", "BiocManager", "fgsea", "qvalue",
         "graph", "GenSA", "Rmpfr")
ok <- sapply(all, function(p) suppressWarnings(requireNamespace(p, quietly = TRUE)))
print(ok)
cat(if (all(ok)) "RPKGS_OK\n" else "RPKGS_INCOMPLETE\n")
