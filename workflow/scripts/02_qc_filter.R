# Quality control: filter and trim based on the quality profiles.
# Truncation lengths come from config (chosen from the data, not defaults).
library(dada2)

fastq_dir <- snakemake@input[["fastq_dir"]]
filtFs_dir <- snakemake@output[["filt_fs"]]
filtRs_dir <- snakemake@output[["filt_rs"]]
dir.create(filtFs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(filtRs_dir, showWarnings = FALSE, recursive = TRUE)

fnFs <- sort(list.files(fastq_dir, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(fastq_dir, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(filtFs_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtRs_dir, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = snakemake@params[["trunc_len"]],
                     maxN = 0,
                     maxEE = snakemake@params[["max_ee"]],
                     truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)

dir.create(dirname(snakemake@output[["stats"]]), showWarnings = FALSE, recursive = TRUE)
write.table(data.frame(sample = rownames(out), out, row.names = NULL),
            file = snakemake@output[["stats"]],
            sep = "\t", row.names = FALSE, quote = FALSE)
