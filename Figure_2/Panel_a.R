
# Script Figure 2a, UMAP of the good learners vs untrained

# memory and library
rm(list = ls())
source('../ancillary/libraries.R')
source('../ancillary/figure_settings.R')

# control panel 
data_file <- '../data/combined_sets.rds'
res_folder <- './Panel_a'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
panel_width <- half_width # mm
panel_height <- 80 # mm

# loading the combined sets
combined_sets <- readRDS(data_file)

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
  scattermore::geom_scattermore(pointsize = 3, alpha = 0.5, 
                                pixels = round(600 * c(panel_width, panel_height) / 25.4)) + 
  scale_color_manual(name = 'Group', 
                     labels = c('Good learners', 'Untrained chicks'),
                     values = c(Untrained = 'blue', 
                                Good_learner = 'red')) + 
  theme_void() + theme_figure() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        legend.position = 'bottom', 
        legend.key.size = unit(3, 'mm'), 
        legend.margin = margin(0, 0, 0, 0)) + 
  guides(color = guide_legend(override.aes = list(size = 1.8, shape = 16, alpha = 1))) 
save_panel(file.path(res_folder, 'umap_good_learners_vs_untrained'), 
           function() plot(p), panel_width, panel_height)
