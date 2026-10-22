# Install all packages needed for the 16S DADA2 -> phyloseq walkthrough.
# Run once with: source("install-packages.R")
# Takes ~10-20 minutes (compiles Bioconductor packages); then you're set.
#
# Prompt-proofing: when a script is source()d, any interactive question
# (install from source? update old packages?) silently consumes the NEXT
# LINE OF THE SCRIPT as its "answer" and dies with
#   Error in askYesNo(msg) : Unrecognized response ...
# The option below plus ask/update/upgrade = never settings everywhere
# keep this script from ever asking a question mid-run.
options(install.packages.compile.from.source = "never")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c("dada2", "phyloseq", "Biostrings", "DESeq2", "DECIPHER", "microbiome",
    "ALDEx2", "ANCOMBC"),
  ask = FALSE, update = FALSE
)

install.packages(c("ggplot2", "corrr", "patchwork", "UpSetR", "metacoder"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# MaAsLin 3 (one of the four differential-abundance methods) REQUIRES
# R >= 4.4 -- every release, and the GitHub version, declares that floor.
# On R < 4.4 (e.g. R 4.3.3 / Bioconductor 3.18) it cannot be installed at
# all, so skip it with a clear message instead of erroring out and
# stopping the rest of this script. To use MaAsLin 3, upgrade R to
# 4.4 or newer first, then re-run this script.
if (getRversion() < "4.4") {
  message("Skipping maaslin3: it requires R >= 4.4 (you have R ",
          as.character(getRversion()), "). Everything else will install.")
} else if (!requireNamespace("maaslin3", quietly = TRUE)) {
  try(suppressWarnings(BiocManager::install("maaslin3", ask = FALSE,
                                             update = FALSE)),
      silent = TRUE)
  if (!requireNamespace("maaslin3", quietly = TRUE)) {
    try(remotes::install_github("biobakery/maaslin3", upgrade = "never"),
        silent = TRUE)
  }
}

# Optional: chkMocks powers the quantitative mock-composition check.
# GitHub-only package. Its Bioconductor dependencies (DECIPHER, microbiome)
# are installed above — remotes cannot fetch Bioconductor packages itself.
remotes::install_github("microsud/chkMocks", upgrade = "never")

if (requireNamespace("maaslin3", quietly = TRUE)) {
  message("All packages installed. You can now knit 16S-DADA2-phyloseq-walkthrough.Rmd.")
} else {
  message("Done -- except maaslin3 (needs R >= 4.4, see note above).")
}
