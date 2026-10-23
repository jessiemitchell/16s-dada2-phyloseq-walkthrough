# Ecological analysis in phyloseq: import DADA2 outputs, alpha/beta diversity,
# composition plots, and a metacoder heat tree. (Differential abundance is
# 07_differential_abundance.R.)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(metacoder)
theme_set(theme_bw())
set.seed(1)

seqtab.nochim <- readRDS(snakemake@input[["seqtab"]])
taxa <- readRDS(snakemake@input[["taxa"]])

# sample metadata is encoded in the filenames (mouse subject, day post-weaning)
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject, 1, 1)
subject <- substr(subject, 2, 999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject = subject, Gender = gender, Day = day)
samdf$When <- "Early"
samdf$When[samdf$Day > 100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(samdf),
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)  # drop the mock control

# short ASV names for plots; full sequences kept in the refseq slot
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
saveRDS(ps, snakemake@output[["ps"]])

dir.create(dirname(snakemake@output[["alpha"]]), showWarnings = FALSE, recursive = TRUE)

p_alpha <- plot_richness(ps, x = "Day",
                         measures = c("Shannon", "Simpson"), color = "When")
ggsave(snakemake@output[["alpha"]], p_alpha, width = 7, height = 4)

ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
ord <- ordinate(ps.prop, method = "NMDS", distance = "bray")
p_ord <- plot_ordination(ps.prop, ord, color = "When", title = "Bray-Curtis NMDS")
ggsave(snakemake@output[["ordination"]], p_ord, width = 6, height = 5)

top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps.top20 <- prune_taxa(top20,
                       transform_sample_counts(ps, function(otu) otu / sum(otu)))
p_comp <- plot_bar(ps.top20, x = "Day", fill = "Family") +
  facet_wrap(~When, scales = "free_x")
ggsave(snakemake@output[["composition"]], p_comp, width = 8, height = 5)

# heat tree (metacoder): node size = ASV count, colour = mean abundance,
# taxonomy cut at genus level
obj <- parse_phyloseq(ps)
obj$data$taxon_abund <- calc_taxon_abund(obj, data = "otu_table")
obj$data$taxon_abund$mean_abund <-
  rowMeans(obj$data$taxon_abund[, -1, drop = FALSE])
obj <- filter_taxa(obj, tolower(taxon_ranks) == "genus", supertaxa = TRUE)
p_tree <- heat_tree(obj,
                    node_label = taxon_names,
                    node_size = n_obs,
                    node_color = mean_abund,
                    node_size_axis_label = "ASV count",
                    node_color_axis_label = "Mean abundance",
                    layout = "davidson-harel",
                    initial_layout = "reingold-tilford")
ggsave(snakemake@output[["heat_tree"]], p_tree, width = 8, height = 8)
