# Consensus differential abundance: DESeq2 + ALDEx2 + ANCOM-BC2 + MaAsLin 3.
# Mirrors the "Differential abundance" section of the R Markdown walkthrough.
# Rationale (Nearing et al. 2022, Nat Commun 13:342): DA methods disagree so
# much on real data that no single method's calls should be trusted alone;
# a taxon counts as robust when at least 2 of the 4 methods call it.
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(ANCOMBC)
library(maaslin3)
set.seed(1)

ps <- readRDS(snakemake@input[["ps"]])

# 10% prevalence filter (the benchmark's filtered analysis), applied before
# testing so all four methods vote on the same taxa
ps.da <- filter_taxa(ps, function(x) sum(x > 0) > 0.10 * length(x), prune = TRUE)

# DESeq2: negative binomial model on raw counts; padj < 0.05
dds <- phyloseq_to_deseq2(ps.da, ~ When)
dds <- DESeq(dds, test = "Wald", fitType = "parametric")
deseq.res <- as(results(dds, cooksCutoff = FALSE), "data.frame")
sig.deseq <- rownames(deseq.res)[which(deseq.res$padj < 0.05)]

# ALDEx2: Dirichlet Monte Carlo + CLR; expected BH-adjusted Wilcoxon p < 0.05
aldex.res <- aldex(t(as(otu_table(ps.da), "matrix")),
                   as.character(sample_data(ps.da)$When),
                   test = "t", effect = TRUE, denom = "all")
sig.aldex <- rownames(aldex.res)[aldex.res$we.eBH < 0.05]

# ANCOM-BC2: bias-corrected log-fold changes; the package's robust call
# (significant at q < 0.05, BH, and stable in the pseudocount sensitivity
# analysis)
out.bc <- ancombc2(data = ps.da, fix_formula = "When", group = "When",
                   p_adj_method = "BH", alpha = 0.05)
res.bc <- out.bc$res
q.col <- grep("^q_", colnames(res.bc), value = TRUE)
q.col <- q.col[!grepl("Intercept", q.col)][1]  # the tested coefficient, e.g. q_WhenLate
suffix <- sub("^q_", "", q.col)
robust.col <- paste0("diff_robust_", suffix)
if (!robust.col %in% colnames(res.bc)) robust.col <- paste0("diff_", suffix)
sig.bc <- res.bc$taxon[res.bc[[robust.col]] %in% TRUE]

# MaAsLin 3: prevalence (logistic) + abundance (linear) hurdle model;
# joint q < 0.1 is the authors' significance convention
meta.m3 <- as(sample_data(ps.da), "data.frame")
meta.m3$When <- factor(meta.m3$When, levels = c("Early", "Late"))
fit.m3 <- maaslin3(input_data = as(otu_table(ps.da), "matrix"),
                   input_metadata = meta.m3,
                   output = snakemake@params[["maaslin_dir"]],
                   formula = ~ When,
                   plot_summary_plot = FALSE, plot_associations = FALSE)
m3.res <- rbind(fit.m3$fit_data_abundance$results,
                fit.m3$fit_data_prevalence$results)
m3.res <- m3.res[m3.res$metadata == "When", ]
sig.m3 <- unique(m3.res$feature[m3.res$qval_joint < 0.1])

# the vote: one row per ASV, one column per method
votes <- data.frame(DESeq2    = taxa_names(ps.da) %in% sig.deseq,
                    ALDEx2    = taxa_names(ps.da) %in% sig.aldex,
                    ANCOM_BC2 = taxa_names(ps.da) %in% sig.bc,
                    MaAsLin3  = taxa_names(ps.da) %in% sig.m3,
                    row.names = taxa_names(ps.da))
votes$n_methods <- rowSums(votes)

dir.create(dirname(snakemake@output[["votes"]]), showWarnings = FALSE, recursive = TRUE)
votes.out <- cbind(data.frame(ASV = rownames(votes), row.names = NULL),
                   votes,
                   DESeq2_log2FC = deseq.res[rownames(votes), "log2FoldChange"],
                   DESeq2_padj = deseq.res[rownames(votes), "padj"],
                   as(tax_table(ps.da), "matrix"))
write.table(votes.out, file = snakemake@output[["votes"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

consensus <- votes.out[votes.out$n_methods >= 2, ]
consensus <- consensus[order(-consensus$n_methods, -abs(consensus$DESeq2_log2FC)), ]
write.table(consensus, file = snakemake@output[["consensus"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

votes01 <- as.data.frame(lapply(votes[, 1:4], as.integer),
                         row.names = rownames(votes))
png(snakemake@output[["upset"]], width = 1400, height = 800, res = 150)
UpSetR::upset(votes01, sets = c("DESeq2", "ALDEx2", "ANCOM_BC2", "MaAsLin3"),
              order.by = "freq",
              mainbar.y.label = "ASVs called by exactly this combination",
              sets.x.label = "Significant ASVs per method")
dev.off()
