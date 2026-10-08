
# Script Figure 1d, UMAP of the data with detailed cell types

# memory and library
rm(list = ls())
source('../ancillary/libraries.R')
source('../ancillary/figure_settings.R')
set.seed(12345)

# control panel 
data_file <- '../data/combined_sets.rds'
classification_file <- '../data/cell_classification.csv'
res_folder <- './Panel_d'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
panel_width <- half_width # mm
panel_height <- 85 # mm

# loading the plot information
cluster_annotation <- read.csv(classification_file, stringsAsFactors = FALSE)
cluster_annotation$Cell_subtype <- case_match(cluster_annotation$Cell_subtype, 
                                                         'Blood_and_immune_system_cells' ~ 'Other',
                                                         'Microglia_Endothelial_Vascular' ~ 'Other', 
                                                         'Oligodendrocytes_PC' ~ 'OPC',
                                                         .default = cluster_annotation$Cell_subtype)
cluster_annotation$Cell_subtype <- gsub('_', ' ', cluster_annotation$Cell_subtype)

# loading the combined sets
combined_sets <- readRDS(data_file)

# re-coding the seurat_clusters
mapping <- setNames(cluster_annotation$Cell_subtype,
                    cluster_annotation$Cluster)
combined_sets@meta.data[['Cell identity']] <- 
  mapping[as.character(combined_sets@meta.data$seurat_clusters)]

# setting assay and identity
DefaultAssay(combined_sets) <- 'integrated'
Idents(combined_sets) <- 'Cell identity'

# Number of clusters / identities
idents <- levels(Idents(combined_sets))
idents <- sort(idents)
n_idents <- length(idents)

# palette
final_palette <- dittoColors()[seq_len(n_idents)]
names(final_palette) <- idents
shared <- intersect(names(final_palette), names(cell_type_colors))
final_palette[shared] <- cell_type_colors[shared]
final_palette[names(final_palette) == 'EXC GLU-7'] <- '#AD0000'

# umap
p <- DimPlot(combined_sets, pt.size = 3, label = FALSE, alpha = 0.6, 
             cols = final_palette, raster = TRUE, 
             raster.dpi = round(600 * c(panel_width, panel_height) / 25.4)) + 
  theme_figure() + NoAxes() + NoLegend()
p <- LabelClusters(p, id = 'ident', fontface = 'bold', color = 'black', 
                   repel = TRUE, size = font_text / ggplot2::.pt, 
                   family = font_family)
save_panel(file.path(res_folder, 'panel_d'), function() plot(p), 
           panel_width, panel_height)
