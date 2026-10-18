 End-to-end 16S rRNA analysis: DADA2 denoising to phyloseq ecology

A fully reproducible walkthrough of a standard 16S rRNA amplicon analysis — raw paired-end fastq files through QC, DADA2 denoising to amplicon sequence variants (ASVs), SILVA taxonomy assignment, and ecological analysis in phyloseq (alpha/beta diversity, composition, DESeq2 differential abundance). Rendered notebook: `16S-DADA2-phyloseq-walkthrough.Rmd` (knit to HTML in RStudio).

## Data

The mothur MiSeq SOP example dataset: 2x250 bp Illumina MiSeq sequencing of the V4 region of the 16S rRNA gene, from fecal samples collected longitudinally from mice post-weaning, plus a mock community of 20 known strains used to validate the pipeline. Public data, downloaded automatically on first run from the mothur project. Taxonomy via the SILVA v138 reference (downloaded automatically from Zenodo on first run).

## How to run

Two ways, same analysis:

**A. The walkthrough (R Markdown).** Open `16S-DADA2-phyloseq-walkthrough.Rmd` in RStudio and knit. Best for reading — full narrative, QC plots inline, and the reasoning behind each parameter choice.

**B. The pipeline (Snakemake).** For rerunning on new data:

```bash
cd workflow
snakemake --use-conda --cores 4
```

The workflow (`workflow/Snakefile`) runs the same steps as modular R scripts (`workflow/scripts/`): download → QC/filter → DADA2 denoising → SILVA taxonomy → phyloseq analysis. Parameters (URLs, truncation lengths, error thresholds) live in `workflow/config/config.yaml`; the conda environment (R + Bioconductor) is pinned in `workflow/envs/r-env.yaml`. Outputs land in `results/`: filtered reads, the ASV table, taxonomy, read-tracking tables, the phyloseq object, figures, and the differential-abundance table.

Requires R with the following packages (Bioconductor for `dada2`, `phyloseq`, `Biostrings`, `DESeq2`; CRAN for `ggplot2`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "DESeq2"))
install.packages("ggplot2")
```

Then open `16S-DADA2-phyloseq-walkthrough.Rmd` in RStudio and knit (or run the Snakemake pipeline above). Expected runtime: ~15–30 minutes on a laptop (most of it in error-rate learning and taxonomy assignment); data and reference downloads happen once and are cached.

## What it shows

- QC-driven parameter choices (truncation lengths chosen from quality profiles, not defaults)
- DADA2 error modeling, denoising, read merging, chimera removal
- A read-tracking table through every pipeline step — the first thing I check in any amplicon project
- Mock-community validation: 20 strains in, 20 exact sequence matches out
- phyloseq ecology: alpha diversity, Bray-Curtis NMDS ordination, composition bar plots, DESeq2 differential abundance

## Author

Jessica Mitchell — PhD computational biologist (Harvard). Microbiome, metagenomics, and NGS analysis.

## License

MIT — see `LICENSE`.
