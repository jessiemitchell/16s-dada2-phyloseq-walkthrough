# End-to-end 16S rRNA analysis: DADA2 denoising to phyloseq ecology

A fully reproducible walkthrough of a standard 16S rRNA amplicon analysis — raw paired-end fastq files through QC, DADA2 denoising to amplicon sequence variants (ASVs), SILVA taxonomy assignment, and ecological analysis in phyloseq (alpha/beta diversity, composition, consensus differential abundance). Rendered notebook: `16S-DADA2-phyloseq-walkthrough.Rmd` (knit to HTML in RStudio).

## Data

The mothur MiSeq SOP example dataset: 2x250 bp Illumina MiSeq sequencing of the V4 region of the 16S rRNA gene, from fecal samples collected longitudinally from mice post-weaning, plus the HMP mock community used to validate the pipeline: genomic DNA from 21 bacterial strains, of which the 20 with reference 16S sequences in `HMP_MOCK.v35.fasta` are recovered exactly (the 21st, *Helicobacter pylori*, yields no reads — see the walkthrough). Public data, downloaded automatically on first run from the mothur project. Taxonomy via the SILVA v138 reference (downloaded automatically from Zenodo on first run).

## How to run

Two ways, same analysis:

**A. The walkthrough (R Markdown).** Open `16S-DADA2-phyloseq-walkthrough.Rmd` in RStudio and knit. Best for reading — full narrative, QC plots inline, and the reasoning behind each parameter choice.

**B. The pipeline (Snakemake).** For rerunning on new data:

```bash
cd workflow
snakemake --use-conda --cores 4
```

The workflow (`workflow/Snakefile`) runs the same steps as modular R scripts (`workflow/scripts/`): download → QC/filter → DADA2 denoising → SILVA taxonomy → phyloseq analysis, plus a mock-community composition check (`06_mock_check.R`) and a consensus differential-abundance analysis (`07_differential_abundance.R`). Parameters (URLs, truncation lengths, error thresholds) live in `workflow/config/config.yaml`; the conda environment (R + Bioconductor) is pinned in `workflow/envs/r-env.yaml`. Outputs land in `results/`: filtered reads, the ASV table, taxonomy, read-tracking tables, the phyloseq object, figures, the per-method DA vote table and consensus table, and the mock-check table.

Requires R with the following packages (Bioconductor for `dada2`, `phyloseq`, `Biostrings`, `DESeq2`; CRAN for `ggplot2`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "DESeq2"))
install.packages("ggplot2")
```

Simplest is to run `install-packages.R` from the repo root instead — it installs everything above plus the rest of the walkthrough's packages (DECIPHER, microbiome, the four-method DA panel, chkMocks, MaAsLin 3).

Then open `16S-DADA2-phyloseq-walkthrough.Rmd` in RStudio and knit (or run the Snakemake pipeline above). Expected runtime: ~15–30 minutes on a laptop (most of it in error-rate learning and taxonomy assignment); data and reference downloads happen once and are cached.

## Contents

- QC-driven parameter choices (truncation lengths chosen from quality profiles, not defaults)
- DADA2 error modeling, denoising, read merging, chimera removal
- A read-tracking table through every pipeline step
- Mock-community validation, two layers: exact sequence recovery (20 strains, 20 exact matches), then a chkMocks composition check — are the recovered proportions quantitatively right, not just present?
- phyloseq ecology: alpha diversity, Bray-Curtis NMDS ordination, composition bar plots, consensus differential abundance — DESeq2, ALDEx2, ANCOM-BC2 and MaAsLin 3 vote on every taxon, because no single DA method is trustworthy alone (Nearing et al. 2022)

## Author

Jessica Mitchell — PhD computational biologist (Harvard). Microbiome, metagenomics, and NGS analysis.

## License

MIT — see `LICENSE`.
