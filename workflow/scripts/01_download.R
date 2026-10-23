# Download the example dataset (mothur MiSeq SOP). Cached after the first run.
url <- snakemake@params[["url"]]
outdir <- snakemake@output[["fastq_dir"]]

dir.create(dirname(outdir), showWarnings = FALSE, recursive = TRUE)
tmpzip <- file.path(dirname(outdir), "miseqsopdata.zip")
download.file(url, destfile = tmpzip)
# the zip contains a top-level MiSeq_SOP/ folder, which lands at data/MiSeq_SOP
unzip(tmpzip, exdir = dirname(outdir))
