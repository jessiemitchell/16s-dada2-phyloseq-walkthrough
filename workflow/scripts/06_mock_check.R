# Mock-community composition check with chkMocks (custom-mock route):
# are the observed strain proportions quantitatively accurate, or just present?
# Mirrors the walkthrough's chkMocks section. Outputs results/tables/mock_check.tsv.
library(Biostrings)
library(phyloseq)
library(chkMocks)

seqtab.nochim <- readRDS(snakemake@input[["seqtab"]])
mock.ref <- getSequences(snakemake@params[["mock_ref"]])

unqs.mock <- seqtab.nochim["Mock", ]
unqs.mock <- sort(unqs.mock[unqs.mock > 0], decreasing = TRUE)

# exact-match each ASV to a strain via the reference headers (abbreviated
# species codes; several strains have multiple distinct operon sequences)
strain.of <- sapply(names(unqs.mock), function(x) {
  hit <- which(sapply(mock.ref, function(r) grepl(x, r)))
  if (length(hit) >= 1) sub("\\.[0-9]+$", "", names(mock.ref)[hit[1]]) else "Other"
})
obs <- tapply(unqs.mock, strain.of, sum)

# the 20 recovered strains (H. pylori is in the pool and the reference file
# but yields no reads; see walkthrough text), representative genome sizes
# (Mb) and 16S operon counts
hmp.strains <- data.frame(
  Code = c("A.baumannii", "A.odontolyticus", "B.cereus", "B.vulgatus",
           "C.beijerinckii", "D.radiodurans", "E.faecalis", "E.coli",
           "L.gasseri", "L.monocytogenes", "N.meningitidis", "P.gingivalis",
           "P.acnes", "P.aeruginosa", "R.sphaeroides", "S.aureus",
           "S.epidermidis", "S.agalactiae", "S.mutans", "S.pneumoniae"),
  GenomeMb = c(3.98, 2.80, 5.43, 5.16, 6.00, 3.28, 3.22, 4.90, 1.89, 2.94,
               2.20, 2.35, 2.56, 6.26, 4.60, 2.80, 2.60, 2.13, 2.03, 2.10),
  rrn      = c(6, 2, 13, 7, 14, 8, 4, 7, 6, 6,
               4, 4, 3, 4, 3, 6, 5, 7, 5, 4))
rownames(hmp.strains) <- hmp.strains$Code

# two bracketing expectations: equal DNA mass per strain (5% each), and
# template-weighted (operon copies per genome / genome size)
hmp.strains$EvenDNA  <- 100 / nrow(hmp.strains)
hmp.strains$Template <- with(hmp.strains,
                             100 * (rrn / GenomeMb) / sum(rrn / GenomeMb))

obs.pct <- 100 * obs / sum(obs)
hmp.strains$Observed <- 0
shared <- intersect(names(obs.pct), rownames(hmp.strains))
hmp.strains[shared, "Observed"] <- obs.pct[shared]

# chkMocks compares samples held side by side in one phyloseq OTU table
comp <- hmp.strains[, c("EvenDNA", "Template", "Observed")]
colnames(comp) <- c("EvenDNA", "TemplateWeighted", "Mock")
ps.chk <- phyloseq(otu_table(as.matrix(comp), taxa_are_rows = TRUE))

chk.out <- data.frame(
  Strain = hmp.strains$Code,
  GenomeMb = hmp.strains$GenomeMb,
  rrn = hmp.strains$rrn,
  EvenDNA_pct = hmp.strains$EvenDNA,
  TemplateWeighted_pct = hmp.strains$Template,
  Observed_pct = hmp.strains$Observed,
  Spearman_vs_template = compare2theorectical(ps.chk,
                                              theoretical_id = "TemplateWeighted")$TemplateWeighted[1],
  Spearman_vs_evenDNA = compare2theorectical(ps.chk,
                                             theoretical_id = "EvenDNA")$EvenDNA[1])

dir.create(dirname(snakemake@output[["mock"]]), showWarnings = FALSE,
           recursive = TRUE)
write.table(chk.out, file = snakemake@output[["mock"]],
            sep = "\t", row.names = FALSE, quote = FALSE)
