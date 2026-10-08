# CHARM-Vis

<img width="973" alt="image" src="https://github.com/user-attachments/assets/e24016db-167e-4fcb-a1e8-0964b9723941" />


## Description

This repository contains code for reproducing results presented in the manuscript: *Charting the single cell transcriptional landscape governing visual imprinting*, available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.23.660422v1).

The repository includes the results of the differential expression and enrichment analyses used in the manuscript, so the figures can be generated without re-running the analyses.

| Folder / file | Content |
|---|---|
| `1_differential_analysis.R` | Pseudo-bulk differential expression analysis (DESeq2) |
| `1_differential_analysis/` | Results of `1_differential_analysis.R` used in the manuscript |
| `2_enrichment_analysis.R` | Gene set enrichment analysis (GSEA) on GO biological processes |
| `2_enrichment_analysis/` | Results of `2_enrichment_analysis.R` used in the manuscript |
| `Figure_1/`, `Figure_2/`, `Figure_3/` | One script per figure panel; panels are written as PDF and PNG (600 dpi) at their final size |
| `Figure_4/` | `Figure_4.R`, producing all panels of Figure 4 |
| `Figure_5/` | `Panel_f.R`, producing Figure 5f (in situ hybridization, IMM) and the same analysis for NeoS |
| `ancillary/` | Functions and package loading shared by the scripts |
| `data/` | Annotation and sample information; the Seurat object must be downloaded (see below) |
| `data/figure_4/` | Behavioural data and validation measurements used for Figure 4 |
| `data/figure_5/` | In situ hybridization measurements (one value per cell) used for Figure 5f |
| `renv.lock` | R package versions of the main environment |
| `Figure_3/renv_panel_a/renv.lock` | R package versions used only for Figure 3a |

## Instructions

1. Clone the repository on your machine.
2. Download `GSE299793_combined_sets.rds` (2.1 GB) from [GSE299793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299793) and save it as `data/combined_sets.rds`. It is needed by `1_differential_analysis.R` and by the Figure 1 and Figure 2a scripts.
3. Install R 4.5.2 (see [r-project.org](https://www.r-project.org/)) and set up the package environments (see [Reproducibility](#reproducibility)).
4. Run any of the scripts in the `Figure_1` to `Figure_5` folders, from inside that folder. For example:
   ```
   cd Figure_1
   Rscript Panel_b.R
   ```
   Each script writes its panel in a subfolder with the same name (e.g. `Figure_1/Panel_b/`). `Figure_4/Figure_4.R` analyses all genes, sides and brain regions and writes plots and statistics in `Figure_4/Results/<gene>/<side>_<region>/`, and the panels of Figure 4 in `Figure_4/Panel_a/` to `Figure_4/Panel_g/` (listed at the top of the script). `Figure_5/Panel_f.R` writes the dotplot, the plotted values and the mixed model results (`results.txt`) for IMM (Figure 5f) and NeoS in `Figure_5/Panel_f/<area>/`.
5. Optionally, re-run the analyses from the repository root with `Rscript 1_differential_analysis.R` and `Rscript 2_enrichment_analysis.R`. This overwrites the results in `1_differential_analysis/` and `2_enrichment_analysis/` (see the notes below on what to expect).

## Reproducibility

The package versions are recorded with [renv](https://rstudio.github.io/renv/). The project uses R 4.5.2 and Bioconductor 3.22.

**Main environment.** From the repository root, run:

```
Rscript -e 'renv::restore()'
```

renv installs itself if needed (through the `.Rprofile` file) and then installs the packages listed in `renv.lock` in a project library. The `.Rprofile` files in the `Figure_*` folders activate the same environment when a script is started from those folders.

**Figure 3a environment.** Figure 3a was made with simplifyEnrichment 1.14.1, which needs older versions of some packages (e.g. GOSemSim 2.30.2) than the rest of the analysis. Its packages are in a separate lockfile, used automatically by `Figure_3/Panel_a.R`. Install them once with:

```
cd Figure_3/renv_panel_a
Rscript -e 'renv::restore()'
```

**Notes**

- Some packages are installed from source. On macOS this requires the Xcode command line tools and the GNU Fortran compiler from [mac.r-project.org/tools](https://mac.r-project.org/tools/); on Linux, the usual system libraries for R packages (e.g. libxml2, libcurl, openssl, fontconfig, harfbuzz, fribidi, freetype, libpng).
- Run R with a UTF-8 locale (e.g. `LANG=en_US.UTF-8`). `data/phenoData.csv` starts with a byte order mark, which R only removes in a UTF-8 locale.
- **Differential expression**: re-running `1_differential_analysis.R` reproduces the results in `1_differential_analysis/` up to numerical precision.
- **Enrichment**: `2_enrichment_analysis.R` uses the GO gene sets of the original analysis (`data/GO_BP_gene_sets_mouse_ensembl_GO_3.21.rds`, from org.Mm.eg.db and GO.db 3.21.0), so the results do not depend on the installed GO annotation, and a fixed seed, so that repeated runs give identical results. GSEA p-values are estimated with random permutations and the original analysis did not fix the seed: a re-run gives the same gene sets and essentially the same enrichment scores (correlation > 0.9999), but terms close to the significance thresholds can change. For this reason, figures made from re-run enrichment results can differ slightly from those in the manuscript.
- **Figure 4 data**: in `data/figure_4/GLUBK89/behavioural.xls` and `data/figure_4/LUC7L/behavioural.xls`, one trained chick was originally labelled `Trained ` (with a trailing space) and was not counted in the trained group; the label is corrected in these files. In `data/figure_4/ROBO1/` and `data/figure_4/RORA/`, the sample identifier of the right IMM of chick a28 is corrected from 254 to 154. `Figure_4.R` removes leading and trailing spaces from the labels and checks the input files (expected labels, same samples in the same order in both files).
- `2_enrichment_analysis.R` takes several hours on a laptop (mostly `simplify()`); the number of parallel workers is set by `ncores`.
- With these environments, the figure scripts give identical images across runs (checked on macOS). Small rendering differences can occur with other operating systems or graphics devices.
