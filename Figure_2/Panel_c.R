
# Script Figure 2c

# memory and library
rm(list = ls())
source('../ancillary/libraries.R')

# control panel 
fc_threshold <- 1
padj_threshold <- 0.05
min_base_mean <- 100 
heatmap_limits <- c(-1, 4)
min_num_cells <- 50
data_file <- '../data/combined_sets.rds'
pheno_file <- '../data/phenoData.csv'
classification_file <- '../data/cell_classification.csv'
de_folder <- '../1_differential_analysis'
res_folder <- './Panel_c'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)


# loading classification
classification <- read.csv(classification_file)

# choosing the analysis
cell_types <- unique(classification$Cell_type)
cell_subtypes <- unique(classification$Cell_subtype)

# colors!
cell_types_color <- c(Astrocytes = '#00BA38',
                      `Gabaergic neurons` = '#B79F00',
                      `Glutamatergic neurons` = '#F8766D',
                      OPC = '#619CFF',
                      Ependymal = '#01BFC4',
                      Other = '#F564E3',
                      `Cell types` = 'white'
                      #XXX , IMM = 'grey'
)

#### cell subtype ####

# loading the DE results 
cell_subtypes_de_res <- vector('list', length(cell_subtypes))
names(cell_subtypes_de_res) <- paste0('Cell_subtype-', cell_subtypes)
for(i in 1:length(cell_subtypes_de_res)){
  
  # reading the number of cells 
  num_cells <- read.table(header = TRUE, file.path(de_folder, names(cell_subtypes_de_res)[i], 
                                                   'num_cells_per_sample.txt'))
  
  # proceeding only if enough cell are measured for each sample
  if(min(num_cells) < min_num_cells){
    next()
  }
  
  # reading the actual results
  cell_subtypes_de_res[[i]] <- read.csv(file.path(de_folder, 
                                                  names(cell_subtypes_de_res)[i], 
                                                  paste0(names(cell_subtypes_de_res)[i], '_results.csv')))
}

# which de results have actual results?
usable_de_res <- cell_subtypes_de_res[which(!sapply(cell_subtypes_de_res, is.null))]
num_res <- length(usable_de_res)

#### logFC correlation matrix ####

# initializng the matrix
corr_m <- matrix(NA, num_res, num_res)
rownames(corr_m) <- colnames(corr_m) <- names(usable_de_res)
symmetric_corr_m <- corr_m

# computing the correlation
for(i in 1:(num_res - 1)){
  
  # first set of lgFC values
  idx <- usable_de_res[[i]]$baseMean >= min_base_mean & 
    !is.na(usable_de_res[[i]]$log2FoldChange)
  de_i <- usable_de_res[[i]]$gene[idx]
  logFC_i <- usable_de_res[[i]]$log2FoldChange[idx]
  names(logFC_i) <- de_i
  
  # looping over the other results
  for(j in (i+1):num_res){
    
    # second set of lgFC values
    idx <- usable_de_res[[j]]$baseMean >= min_base_mean & 
      !is.na(usable_de_res[[j]]$log2FoldChange)
    de_j <- usable_de_res[[j]]$gene[idx]
    logFC_j <- usable_de_res[[j]]$log2FoldChange[idx]
    names(logFC_j) <- de_j
    
    # common logFC values
    common_de <- intersect(de_i, de_j)
    tmp <- cor.test(logFC_i[common_de], logFC_j[common_de], method = 'spearman')
    corr_m[i, j] <- symmetric_corr_m[i, j] <- symmetric_corr_m[j, i] <- tmp$estimate
    corr_m[j, i] <- tmp$p.value
    
  }
  
}

# self correlation to 1
for(i in 1:num_res){
  symmetric_corr_m[i, i] <- 1
}

#### plotting ####

# vector for indicating the cell types for each cluster
splitting_vector <- c()
for(i in 1:num_res){
  idx <- which(paste0('Cell_subtype-', classification$Cell_subtype) == names(usable_de_res)[i])[1]
  if(grepl('EXC', names(usable_de_res)[i]) | grepl('INH', names(usable_de_res)[i])){
    splitting_vector <- c(splitting_vector, 
                          classification$Cell_type[idx])
  }else{
    splitting_vector <- c(splitting_vector, 
                          classification$Cell_subtype[idx])
  }
  
}

# beautifying
splitting_vector <- gsub('_', ' ', splitting_vector)
splitting_vector <- gsub('Oligodendrocytes PC', 'OPC', splitting_vector)
splitting_vector <- factor(splitting_vector, levels = unique(splitting_vector))
rownames(symmetric_corr_m) <- colnames(symmetric_corr_m) <- 
  rownames(corr_m) <- colnames(corr_m) <- gsub('Cell_subtype-', '', colnames(corr_m))
rownames(symmetric_corr_m) <- colnames(symmetric_corr_m) <- 
  rownames(corr_m) <- colnames(corr_m) <- gsub('_', ' ', colnames(corr_m))
rownames(symmetric_corr_m) <- colnames(symmetric_corr_m) <- 
  rownames(corr_m) <- colnames(corr_m) <- gsub('Oligodendrocytes PC', 'OPC', colnames(corr_m))

# annotations
column_ha = HeatmapAnnotation(`Cell type` = splitting_vector, 
                              show_annotation_name = FALSE, 
                              show_legend = FALSE, 
                              col = list(`Cell type` = cell_types_color))

# which matrix to plot?
to_plot <- symmetric_corr_m # corr_m

# removing correlations less than 0.5
#to_plot[to_plot < 0.65] <- 0

# logFC heatmap between cell_subtypes
p <- Heatmap(to_plot, name = "Correlation", 
             #col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c('blue3', 'lightblue', "white", "pink", 'red3')), 
             col = colorRamp2(c(0, 0.5, 1), c("white", "pink", 'red3')), 
             top_annotation = column_ha, 
             show_column_names = FALSE,
             show_row_dend = FALSE,
             clustering_distance_rows = 'spearman',
             clustering_distance_columns = 'spearman', 
             column_dend_height = unit(3, "cm"),
             heatmap_legend_param = list(direction = 'horizontal')
)
png(filename = file.path(res_folder, 'logFC_correlation_between_cell_subtypes.png'), 
    width = 6000, height = 5700, res = 600)
plot(p)
dev.off()

