
# Script Figure 1d, UMAP of the data with detailed cell types

# memory and library
rm(list = ls())
source('../ancillary/libraries.R')
set.seed(12345)

# control panel 
data_file <- '../data/combined_sets.rds'
classification_file <- '../data/cell_classification.csv'
res_folder <- './Panel_d'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

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
final_palette[names(final_palette) == 'Astrocytes'] <- '#00BA38'
final_palette[names(final_palette) == 'OPC'] <- '#619CFF'
final_palette[names(final_palette) == 'Ependymal'] <- '#01BFC4'
final_palette[names(final_palette) == 'Other'] <- '#F564E3'
final_palette[names(final_palette) == 'EXC GLU-7'] <- '#AD0000'

# umap
p <- DimPlot(combined_sets, pt.size = 0.1, label.size = 6, repel = TRUE,
             label = FALSE, alpha = 0.6, cols = final_palette) + NoLegend()
p <- LabelClusters(p, id = 'ident', fontface = "bold", 
                   color = "black", repel = TRUE, size = 6)
png(filename = file.path(res_folder, 'panel_d.png'), 
    width = 7500, height = 6600, res = 600)
plot(p)
dev.off()
