
# Script for Panel c, heatmap with main markers

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
data_file <- '../data/combined_sets.rds'
res_folder <- 'Panel_c'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
markers <- list(gabaergic = c('GAD1', 'GAD2', 'SLC32A1', 'SLC6A1'), 
                glutamatergic = c('ARPP21', 'SV2B', 'SLC17A6', 'SATB2'), 
                astrocite = c('GLI3', 'SLC39A12', 'AQP4'), 
                oligodendrocite = c('OLIG2', 'SOX10'))

# main cell types
cell_types <- list(`Gabaergic neurons` = sort(c(26, 22, 34, 27, 15, 28, 14, 31, 19, 29)), 
                   Oligodendrocytes = c(20, 25),
                   Astrocytes = c(6, 11, 23), 
                   Other = c(32, 35))
cell_types$`Glutamatergic neurons` <- setdiff(0:35, unlist(cell_types))
cell_types_color <- c(Astrocytes = '#00BF7C',
                      `Gabaergic neurons` = '#A3A400',
                      `Glutamatergic neurons` = '#F8766C', 
                      Oligodendrocytes = '#00AFF5',
                      Other = '#E76BF2')

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

# plotting main types
png(filename = file.path(res_folder, 'panel_c.png'), 
    width = 5400, height = 4200, res = 600)
dittoSeq::dittoHeatmap(aggregated_sets, 
                       annot.colors = cell_types_color, 
                       genes = c(markers$gabaergic, 
                                 markers$glutamatergic, 
                                 markers$oligodendrocite, 
                                 markers$astrocite),
                       annot.by = 'Cell_type', 
                       show_colnames = TRUE, 
                       cluster_cols = TRUE)
dev.off()

