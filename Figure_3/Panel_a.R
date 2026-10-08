
#### Script for creating Figure 3 a ####

# Figure 3a has its own package environment (renv_panel_a/renv.lock):
# it needs simplifyEnrichment 1.14.1, which is not compatible with the
# packages of the main environment (see README.md)
Sys.setenv(RENV_PROJECT = normalizePath('renv_panel_a'))
source('renv_panel_a/renv/activate.R')

# memory and library
rm(list = ls())
suppressPackageStartupMessages({
  library(simplifyEnrichment)
  library(ComplexHeatmap)
  library(circlize)
})
source('../ancillary/simplifyGOFromMultipleLists.R')
source('../ancillary/figure_settings.R')

# control panel
enr_folder <- '../2_enrichment_analysis'
minClusterSize <- 20
seed <- 12345 # the colors of the keywords are chosen at random
res_folder <- 'Panel_a'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)
panel_width <- full_width # mm
panel_height <- 90 # mm

#### cell subtypes ####

# choosing the analysis
en_res_names <- dir(enr_folder)
en_res_names <- en_res_names[grep('EXC|INH', en_res_names)]
num_res <- length(en_res_names)

# loading the enrichment results
enr_res <- vector('list', length(en_res_names))
names(enr_res) <- en_res_names
for(i in 1:num_res){
  
  # files to read
  file_to_read <- file.path(enr_folder, en_res_names[i], 
                            #paste0(en_res_names[i], '_gseGO_res_filtered.rds'))
                            paste0(en_res_names[i], '_gseGO_res.rds'))
  
  # read them
  enr_res[[i]] <- readRDS(file_to_read)
  
}

# formatting names
names(enr_res) <- gsub('Cell_type-|Cell_subtype-', '', names(enr_res))
names(enr_res) <- gsub('_', ' ', names(enr_res))
names(enr_res) <- gsub('PC', '(PC)', names(enr_res))

# selecting the up regulated
up_regulated <- lapply(enr_res, function(x){
  tmp <- x@result # same as as.data.frame(x), without loading clusterProfiler
  tmp <- tmp[tmp$NES > 0, ]
})

# font sizes of the heatmaps
ht_opt(legend_title_gp = gpar(fontsize = font_text, fontface = 'bold'), 
       legend_labels_gp = gpar(fontsize = font_text), 
       heatmap_column_names_gp = gpar(fontsize = font_text))

# simplify enrichment (the seed is set before each drawing, so that the pdf and
# png files have the same colors)
draw_panel <- function(){
  set.seed(seed)
  simplifyGOFromMultipleLists(lt = up_regulated, 
                              method = "dynamicTreeCut", 
                              control = list(minClusterSize = minClusterSize), 
                              column_title = character(0), 
                              show_bar_labels = FALSE, 
                              fontsize_label = font_text, 
                              fontsize_axis = font_small, 
                              column_width = unit(2.5, 'mm'), 
                              fontsize_range = c(font_small, 8), 
                              word_cloud_grob_param = list(max_width = 45))
}
save_panel(file.path(res_folder, 'Panel_a'), draw_panel, panel_width, panel_height)


