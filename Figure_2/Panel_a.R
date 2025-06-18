
# Script Figure 2a, UMAP of the good learners vs untrained

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
data_file <- '../data/combined_sets.rds'
res_folder <- './Panel_a'
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

# extracing the meta data
meta_data <- combined_sets@meta.data

# adding the umap information
umap <- combined_sets@reductions$umap@cell.embeddings
meta_data <- cbind(meta_data, 
                   umap_1 = umap[rownames(meta_data), 'umap_1'],
                   umap_2 = umap[rownames(meta_data), 'umap_2'])

# plotting the umap for contrasting good leaners and untrained chicks distribution
p <- ggplot(data = meta_data, 
            mapping = aes(x = umap_1, y = umap_2, 
                          color = group)) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  scale_x_continuous(name = NULL, breaks = NULL) + 
  scale_y_continuous(name = NULL, breaks = NULL) +
  scale_color_manual(name = 'Group', 
                     labels = c('Good learners', 'Untrained chicks'),
                     values = c(Untrained = 'blue', 
                                Good_learner = 'red')) + 
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) + 
  guides(color = guide_legend(override.aes = 
                                list(size = 15, shape = 'square'))) 
png(filename = file.path(res_folder, 'umap_good_learners_vs_untrained.png'), 
    width = 6000, height = 6300, res = 600)
plot(p)
dev.off()

