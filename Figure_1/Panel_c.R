
# Script for Panel c, heatmap with main markers

# memory and library
rm(list = ls())
source('../ancillary/libraries.R')
source('../ancillary/figure_settings.R')

# control panel 
data_file <- '../data/combined_sets.rds'
res_folder <- 'Panel_c'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
panel_width <- half_width # mm
panel_height <- 85 # mm
markers <- list(gabaergic = c('GAD1', 'GAD2', 'SLC32A1'), 
                glutamatergic = c('ARPP21', 'SV2B', 'SLC17A6', 'SATB2'), 
                astrocite = c('GLI3', 'AQP4'), 
                ependymal = c('ENSGALG00010011199'), # 'LRRIQ1'
                oligodendrocite = c('SOX10', 'PDGFRA'))#c('OLIG2', 'SOX10', 'PDGFRA'))

# main cell types
cell_types <- list(`Gabaergic neurons` = sort(c(26, 22, 34, 27, 15, 28, 14, 31, 19, 29)), 
                   OPC = c(20, 25),
                   Astrocytes = c(6, 11), 
                   Ependymal = 23,
                   Other = c(32, 35))
cell_types$`Glutamatergic neurons` <- setdiff(0:35, unlist(cell_types))

# loading the combined sets
combined_sets <- readRDS(data_file)

# setting the default assay
DefaultAssay(combined_sets) <- "RNA"

# retaining only the markers in our data
for(i in 1:length(markers)){
  markers[[i]] <- intersect(markers[[i]], rownames(combined_sets))
}

# cosmetic change
combined_sets[['Clusters']] <- paste0('C', combined_sets$seurat_clusters)

# custom data aggregation
plain_markers <- intersect(rownames(combined_sets), unlist(markers))
tmp <- as.matrix(combined_sets@assays$RNA@data[plain_markers, ])
tmp <- data.frame(t(tmp), cluster = combined_sets[['Clusters']])
tmp <- tmp %>% group_by(Clusters) %>% 
  summarise_all(mean)
tmp <- as.data.frame(tmp)
Clusters <- tmp$Clusters
tmp$Clusters <- NULL
rownames(tmp) <- Clusters
tmp <- t(tmp)

# renaming ENSGALG00010011199 to LRRIQ1
rownames(tmp)[rownames(tmp) == 'ENSGALG00010011199'] <- 'LRRIQ1'
markers$ependymal <- 'LRRIQ1'

# building the meta data
meta_data <- data.frame(Clusters = colnames(tmp))
rownames(meta_data) <- meta_data$Clusters
meta_data$Cell_type <- '' 
for(i in 1:length(cell_types)){
  idx <- paste0('C', cell_types[[i]])
  meta_data[idx, 'Cell_type'] <- names(cell_types)[i]
}

# SummarizedExperiment for dittoHeatmap
aggregated_sets <- SummarizedExperiment(assays = list(RNA = tmp), 
                                        colData = DataFrame(meta_data))

# plotting main types (the colors of the cell types are explained by the
# legend of panel b)
p <- dittoSeq::dittoHeatmap(aggregated_sets, 
                            genes = c(markers$gabaergic, 
                                      markers$glutamatergic, 
                                      markers$oligodendrocite, 
                                      markers$astrocite, 
                                      markers$ependymal),
                            annot.by = 'Cell_type', 
                            annotation_colors = list(Cell_type = cell_type_colors), 
                            annotation_legend = FALSE, 
                            annotation_names_col = FALSE, 
                            show_colnames = TRUE, 
                            cluster_cols = TRUE, 
                            fontsize = font_text, 
                            fontsize_row = font_small, 
                            fontsize_col = font_small, 
                            treeheight_row = 15, 
                            treeheight_col = 15, 
                            silent = TRUE)
# smaller color scale (the legend uses at most the height of its viewport)
legend_idx <- which(p$gtable$layout$name == 'legend')
p$gtable$grobs[[legend_idx]] <- 
  grid::gTree(children = grid::gList(p$gtable$grobs[[legend_idx]]), 
              vp = grid::viewport(y = 1, just = 'top', 
                                  height = grid::unit(25, 'mm')))

draw_heatmap <- function(){
  # small margin on the right, so that the labels of the color scale are not cut
  grid::pushViewport(grid::viewport(x = 0, just = 'left', 
                                    width = grid::unit(1, 'npc') - grid::unit(2, 'mm')))
  grid::grid.draw(p$gtable)
  grid::popViewport()
}
save_panel(file.path(res_folder, 'panel_c'), draw_heatmap, panel_width, panel_height)
