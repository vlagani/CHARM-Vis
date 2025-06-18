
# Script Figure 1d, UMAP with relevant markers

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
data_file <- '../data/combined_sets.rds'
res_folder <- './Panel_d'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# loading the plot information
plot_info <- read.csv('panel_d_notation.csv', stringsAsFactors = FALSE)

# set to grey clusters where we do not characterize their markers
plot_info$Color[plot_info$Notation == ''] <- 'grey90'

# loading the combined sets
combined_sets <- readRDS(data_file)

# extracting the meta data
meta_data <- combined_sets@meta.data

# merging meta data and plot information
meta_data$seurat_clusters <- as.numeric(as.character(meta_data$seurat_clusters))
meta_data$cell_id <- rownames(meta_data)
colnames(plot_info)[1] <- 'seurat_clusters'
meta_data <- merge(meta_data, plot_info)
rownames(meta_data) <- meta_data$cell_id

# merging meta data and umap information
umap <- as.data.frame(combined_sets@reductions$umap@cell.embeddings)
umap$cell_id <- rownames(umap)
meta_data <- merge(meta_data, umap)
rownames(meta_data) <- meta_data$cell_id

# making seurat_clusters a character
meta_data$seurat_clusters <- as.character(meta_data$seurat_clusters)

# merging cluster 0 and 1
meta_data$seurat_clusters[meta_data$seurat_clusters == '0' | 
                            meta_data$seurat_clusters == '1'] <- '0_1'

# computing the centroid and other info for each cluster
cluster_info <- meta_data %>% group_by(seurat_clusters) %>% 
  summarise(notation = unique(Notation), 
            color = unique(Color), 
            mean_umap_1 = mean(umap_1, trim = 0.1),
            mean_umap_2 = mean(umap_2, trim = 0.1)) %>%
  as.data.frame()

# plotting
color_values <- cluster_info$color
names(color_values) <- cluster_info$seurat_clusters
p <- ggplot() + 
  geom_point(mapping = aes(x = umap_1, y = umap_2, 
                           color = seurat_clusters), 
             data = meta_data, size = 0.1, alpha = 0.9) + 
  scale_x_continuous(name = NULL, labels = NULL) + 
  scale_y_continuous(name = NULL, labels = NULL) + 
  scale_color_manual(values = color_values) + 
  geom_text_repel(mapping = aes(x = mean_umap_1, y = mean_umap_2, 
                                label = notation), data = cluster_info,
            color = 'black', size = 7) + 
  theme_classic() + 
  theme(legend.position = 'none') 
png(filename = file.path(res_folder, 'panel_d.png'), 
    width = 6600, height = 6600, res = 600)
plot(p)
dev.off()



