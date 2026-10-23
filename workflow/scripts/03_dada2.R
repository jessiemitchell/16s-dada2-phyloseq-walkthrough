# DADA2 denoising: learn error rates, infer ASVs, merge pairs,
# build the sequence table, remove chimeras.
library(dada2)

filtFs <- sort(list.files(snakemake@input[["filt_fs"]], full.names = TRUE))
filtRs <- sort(list.files(snakemake@input[["filt_rs"]], full.names = TRUE))
# filtered files are named SAMPLENAME_F_filt.fastq.gz
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(filtFs) <- sample.names
names(filtRs) <- sample.names

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)

saveRDS(seqtab.nochim, snakemake@output[["seqtab"]])

# read-tracking table: no step should lose a majority of reads
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
write.table(data.frame(sample = rownames(track), track, row.names = NULL),
            file = snakemake@output[["track"]],
            sep = "\t", row.names = FALSE, quote = FALSE)
