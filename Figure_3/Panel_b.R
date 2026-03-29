# Script for plotting single pathways

#### set up ####

# memory and libraries
rm(list = ls())
source('../ancillary/libraries.R')

# control panel
gse_results_file <- '../2_enrichment_analysis/Cell_type-Glutamatergic_neurons/Cell_type-Glutamatergic_neurons_gseGO_res.rds'
go_ids <- c('GO:0007612', 'GO:0008306', 'GO:0007613')
res_folder <- 'Panel_b'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# loading the results
gse_res <- readRDS(gse_results_file)

# selecting the GOs 
ids <- which(gse_res@result$ID %in% go_ids)

# gsea plot
p <- gseaplot2(gse_res, geneSetID = ids, base_size = 30, rel_heights = c(2, 0.5, 1))
png(filename = file.path(res_folder, 'Panel_b.png'),
    width = 5200, height = 3000, res = 300)
p
dev.off()

