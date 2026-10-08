
# Script Figure 1b, UMAP of the data with cell types

# memory and library
rm(list = ls())
source('../ancillary/libraries.R')
source('../ancillary/figure_settings.R')

# control panel 
data_file <- '../data/combined_sets.rds'
res_folder <- './Panel_b'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
panel_width <- half_width # mm
panel_height <- 80 # mm
universe <- list(OPC = c(20, 25),
                 Astrocytes = c(6, 11),
                 Ependymal = 23,
                 Other = c(32, 35),
                 `Gabaergic neurons`= c(26, 22, 34, 27, 15, 
                                        28, 14, 31, 19, 29))
universe$`Glutamatergic neurons` <- setdiff(0:35, unlist(universe))

# loading the combined sets
combined_sets <- readRDS(data_file)

# creating the cell identity
combined_sets$`Cell identity` <- ''
for(i in 1:length(universe)){
  idx <- combined_sets$seurat_clusters %in% as.character(universe[[i]])
  combined_sets$`Cell identity`[idx] <- names(universe)[i]
}

# umap
DefaultAssay(combined_sets) <- 'integrated'
combined_sets$`Cell identity` <- factor(combined_sets$`Cell identity`, 
                                        levels = names(cell_type_colors))
Idents(combined_sets) <- 'Cell identity'
p <- DimPlot(combined_sets, pt.size = 3, label = FALSE, alpha = 0.9, 
             cols = cell_type_colors, raster = TRUE, 
             raster.dpi = round(600 * c(panel_width, panel_height) / 25.4)) + 
  theme_figure() + NoAxes() + 
  theme(legend.position = 'bottom', 
        legend.key.size = unit(3, 'mm'), 
        legend.key.spacing.y = unit(0.5, 'mm'),
        legend.margin = margin(0, 0, 0, 0)) + 
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 1.5)))
save_panel(file.path(res_folder, 'panel_b'), function() plot(p), 
           panel_width, panel_height)
