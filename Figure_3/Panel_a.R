
# Script for creating Figure 3 a

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
enr_folder <- '../2_enrichment_analysis'
res_folder <- './Panel_a'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# choosing the analysis
cell_type_list <- list(Gabaergic_neurons = c(14, 15, 19, 22, 26, 27, 28, 29, 31, 34), 
                       Glutamatergic_neurons = c(0:5, 7:10, 12, 13, 16:18, 21, 24, 30, 33),
                       Astrocytes = c(6, 11, 23), Oligodendrocytes = c(20, 25))

# loading the enrichment results
enr_res <- vector('list', length(cell_type_list))
names(enr_res) <- names(cell_type_list)
for(i in 1:length(enr_res)){

  # file to read
  file_to_read <- file.path(enr_folder, names(cell_type_list[i]), 
                            paste0(names(cell_type_list[i]), '_gseGO_res.rds'))
  
  # id the file exist, read it
  if(file.exists(file_to_read)){
    enr_res[[i]] <- readRDS(file_to_read)
  }
  
}


# formatting names
names(enr_res) <- gsub('_', ' ', names(enr_res))

# selecting the up regulated
up_regulated <- lapply(enr_res, function(x){
  tmp <- as.data.frame(x)
  tmp <- tmp[tmp$NES > 0, ]
})

# simplify enrichment
png(filename = file.path(res_folder, 'Panel_a.png'), 
    width = 3600, height = 2100, res = 300)
cell_type_up_res <- simplifyGOFromMultipleLists(lt = up_regulated, 
                                                method = "louvain", 
                                                control = list(resolution = 1.1))
dev.off()

