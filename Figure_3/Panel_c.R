# Script for plotting single pathways

#### set up ####

# memory and libraries
rm(list = ls())
source('../libraries.R')

# control panel
gse_results_file <- '../../../analysis_snRNA-seq/2024_02_24/analysis_cell_bender_last/5.4_enrichment_analysis_mouse/integration_approach_1/Glutamatergic_neurons/Glutamatergic_neurons_gseGO_res.rds'
go_ids <- c('GO:0007612', 'GO:0008306', 'GO:0007613')
res_folder <- 'Panel_c'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# loading the results
gse_res <- readRDS(gse_results_file)

# selecting the GOs 
ids <- which(gse_res@result$ID %in% go_ids)

# gsea plot
p <- gseaplot2(gse_res, geneSetID = ids, base_size = 30, rel_heights = c(2, 0.5, 1))
png(filename = file.path(res_folder, 'Panel_c.png'),
    width = 5200, height = 3000, res = 300)
p
dev.off()
  
