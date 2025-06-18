
# Script Figure 2c

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
min_base_mean <- 100 # we want to focus on enough expressed genes
min_num_cells <- 50
data_file <- '../data/combined_sets.rds'
pheno_file <- '../data/phenoData.csv'
de_folder <- '../1_differential_analysis'
res_folder <- './Panel_c'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# choosing the analysis
cell_type_list <- list(Gabaergic_neurons = c(14, 15, 19, 22, 26, 27, 28, 29, 31, 34), 
                      Glutamatergic_neurons = c(0:5, 7:10, 12, 13, 16:18, 21, 24, 30, 33),
                      Astrocytes = c(6, 11, 23), Oligodendrocytes = c(20, 25))
clusters <- unlist(cell_type_list)

# colors!
cell_type_colors <- c(`Gabaergic neurons` = '#A3A400', 
                      `Glutamatergic neurons` = '#F8766C',
                      Astrocytes = '#00BF7C', 
                      Oligodendrocytes = '#00AFF5', 
                      `Cell types` = 'white')

# loading the combined sets and meta data
combined_sets <- readRDS(data_file)
DefaultAssay(combined_sets) <- 'RNA'
Idents(combined_sets) <- 'seurat_clusters'

# loading the meta data
pheno_data <- read.csv(pheno_file)
pheno_data <- pheno_data[pheno_data$X10X_Serial_ID %in% combined_sets$orig.ident, ]
pheno_data <- pheno_data[order(pheno_data$Group), ]

# loading the DE results (clusters)
clusters_de_res <- vector('list', length(clusters))
names(clusters_de_res) <- paste0('Cluster_', clusters)
for(i in 1:length(clusters_de_res)){
  
  # reading the number of cells 
  num_cells <- read.table(header = TRUE, file.path(de_folder, names(clusters_de_res)[i], 
                                                   'num_cells_per_sample.txt'))
  
  # proceeding only if enough cell are measured for each sample
  if(min(num_cells) < min_num_cells){
    next()
  }
  
  # reading the actual results
  clusters_de_res[[i]] <- read.csv(file.path(de_folder, names(clusters_de_res)[i], 
                                             paste0(names(clusters_de_res)[i], '_results.csv')))
}

# renaming due to legacy changes
de_res <- clusters_de_res
num_res <- length(de_res)

# which de results have actual results?
usable_de_res <- which(!sapply(clusters_de_res, is.null))

#### logFC correlation matrix ####

# initializng the matrix
corr_m <- matrix(NA, num_res, num_res)
rownames(corr_m) <- colnames(corr_m) <- gsub('_', ' ', names(de_res))
symmetric_corr_m <- corr_m

# computing the correlation
for(i in 1:(num_res - 1)){
  
  # first set of lgFC values
  if(!is.null(de_res[[i]])){
    idx <- de_res[[i]]$baseMean >= min_base_mean & 
      !is.na(de_res[[i]]$log2FoldChange)
    de_i <- de_res[[i]]$gene[idx]
    logFC_i <- de_res[[i]]$log2FoldChange[idx]
    names(logFC_i) <- de_i
  }else{
    next()
  }
  
  # looping over the other results
  for(j in (i+1):num_res){
    
    # if no j results, skip?
    if(is.null(de_res[[j]])){
      next()
    }
    
    # second set of lgFC values
    idx <- de_res[[j]]$baseMean >= min_base_mean & 
      !is.na(de_res[[j]]$log2FoldChange)
    de_j <- de_res[[j]]$gene[idx]
    logFC_j <- de_res[[j]]$log2FoldChange[idx]
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
for(i in 1:length(cell_type_list)){
  splitting_vector <- c(splitting_vector, 
                        rep(names(cell_type_list)[i], 
                            length(cell_type_list[[i]])))
}
splitting_vector <- gsub('_', ' ', splitting_vector)
splitting_vector <- factor(splitting_vector, levels = unique(splitting_vector))

# keeping only the usable part of the splitting vector
splitting_vector <- splitting_vector[usable_de_res]

# annotations
column_ha = HeatmapAnnotation(`Cell type` = splitting_vector, 
                              show_annotation_name = FALSE, 
                              show_legend = FALSE, 
                              col = list(`Cell type` = cell_type_colors))

# keeping only usable correlations
to_plot <- symmetric_corr_m[usable_de_res, usable_de_res]

# logFC heatmap between clusters
p <- Heatmap(to_plot, name = "Correlation", 
             col = colorRamp2(c(0, 0.5, 1), c("white", "pink", 'red3')), 
             top_annotation = column_ha, 
             show_column_names = FALSE,
             show_row_dend = FALSE,
             clustering_distance_rows = 'spearman',
             clustering_distance_columns = 'spearman', 
             column_dend_height = unit(3, "cm"),
             heatmap_legend_param = list(direction = 'horizontal')
             )
png(filename = file.path(res_folder, 'logFC_correlation_between_clusters.png'), 
    width = 6000, height = 5700, res = 600)
plot(p)
dev.off()

