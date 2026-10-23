# Taxonomic classification with DADA2's naive Bayesian classifier
# against the SILVA v138 reference.
library(dada2)

seqtab.nochim <- readRDS(snakemake@input[["seqtab"]])

silva_file <- "data/silva_nr99_v138.1_train_set.fa.gz"
# ~131 MB file: raise R's 60-second default download timeout, and download to a
# temp file so a failed download leaves no partial file behind.
if (!file.exists(silva_file) || file.size(silva_file) < 1.3e8) {
  dir.create(dirname(silva_file), showWarnings = FALSE, recursive = TRUE)
  options(timeout = 1200)
  tmp <- tempfile(fileext = ".fa.gz")
  download.file(snakemake@params[["silva_url"]], destfile = tmp,
                method = "libcurl", mode = "wb")
  file.rename(tmp, silva_file)
}

taxa <- assignTaxonomy(seqtab.nochim, silva_file, multithread = TRUE)
saveRDS(taxa, snakemake@output[["taxa"]])
