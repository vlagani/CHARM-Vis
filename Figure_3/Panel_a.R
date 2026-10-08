
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

# control panel
enr_folder <- '../2_enrichment_analysis'
minClusterSize <- 20
seed <- 12345 # the colors of the keywords are chosen at random
res_folder <- 'Panel_a'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

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

# simplify enrichment
set.seed(seed)
png(filename = file.path(res_folder, 'Panel_a.png'), 
    width = 4200, height = 2100, res = 300)
cell_type_up_res <- simplifyGOFromMultipleLists(lt = up_regulated, 
                                                method = "dynamicTreeCut", 
                                                control = list(minClusterSize = minClusterSize))
dev.off()


