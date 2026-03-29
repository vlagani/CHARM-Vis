# Script for Figure 3c

# set up
rm(list = ls())
source('../ancillary/libraries.R')

# control panel
enr_res_file <- '../2_enrichment_analysis/Cell_type-Glutamatergic_neurons/Cell_type-Glutamatergic_neurons_gseGO_res_filtered.rds' # enrichment file to load
significance_threshold <- 0.1
res_folder <- 'Panel_c'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# loading the enrichment results
enr_res <- readRDS(enr_res_file)

# focusing on up regulated, significant pathways for plotting
gse_res_up <- enr_res
gse_res_up@result <- gse_res_up@result[gse_res_up@result$NES > 0, ]
gse_res_up@result <- gse_res_up@result[gse_res_up@result$p.adjust <= significance_threshold, ]

# dotplot
p <- dotplot(gse_res_up)
png(filename = file.path(res_folder, 'Panel_c.png'),
    width = 2100, height = 3000, res = 300)
plot(p)
dev.off()
