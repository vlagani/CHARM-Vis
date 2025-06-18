
# Script Figure 1b, UMAP of the data with cell types

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
data_file <- '../data/combined_sets.rds'
res_folder <- './Panel_b'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
universe <- list(Oligodendrocytes = c(20, 25),
                 Astrocytes = c(6, 11, 23),
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
Idents(combined_sets) <- 'Cell identity'
p <- DimPlot(combined_sets, pt.size = 0.1, 
             label = FALSE, alpha = 0.9)
png(filename = file.path(res_folder, 'panel_b.png'), 
    width = 7500, height = 6600, res = 600)
plot(p)
dev.off()
j