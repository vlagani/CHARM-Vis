
# Script for Figure 2 d

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
fc_threshold <- 1
padj_threshold <- 0.05
min_base_mean <- 100
num_genes_to_plot <- 5
min_num_cells <- 50
data_file <- '../data/combined_sets.rds'
pheno_file <- '../data/phenoData.csv'
de_folder <- '../1_differential_analysis'
res_folder <- './Panel_d'
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
                      `.` = 'white')

# loading the combined sets and meta data
combined_sets <- readRDS(data_file)
DefaultAssay(combined_sets) <- 'RNA'
Idents(combined_sets) <- 'seurat_clusters'

# loading the meta data
pheno_data <- read.csv(pheno_file)
pheno_data <- pheno_data[pheno_data$X10X_Serial_ID %in% combined_sets$orig.ident, ]
pheno_data <- pheno_data[order(pheno_data$Group), ]

# loading the DE results (cell type)
de_res <- vector('list', length(cell_type_list))
names(de_res) <- names(cell_type_list)
for(i in 1:length(cell_type_list)){
  de_res[[i]] <- read.csv(file.path(de_folder, names(cell_type_list)[i], 
                                    paste0(names(cell_type_list)[i], '_results.csv')))
}

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

# joining the DE results
de_res <- c(de_res, clusters_de_res)

# which de results have actual results?
usable_de_res <- which(!sapply(clusters_de_res, is.null))

# compiling the list of num_genes_to_plot DE genes with highest logFC
de_genes <- data.frame(cell_type = c(), gene = c())
for(i in 1:length(cell_type_list)){
  tmp <- de_res[[i]]
  tmp <- tmp[tmp$log2FoldChange >= fc_threshold & tmp$padj <= padj_threshold & 
               tmp$baseMean >= min_base_mean, ]
  tmp <- tmp[!duplicated(tmp$gene), ]
  #tmp <- tmp[!grepl('ENSGALG', tmp$gene), ]
  tmp <- tmp[order(abs(tmp$log2FoldChange), decreasing = TRUE), ]
  genes <- tmp$gene[1:num_genes_to_plot]
  de_genes <- rbind(de_genes, data.frame(cell_type = names(de_res)[i], gene = genes))
}

# matrix to plot
to_plot <- matrix(NA, dim(de_genes)[1], length(clusters_de_res) + 1)
colnames(to_plot) <- c('Whole cell type', names(clusters_de_res))
rownames(to_plot) <- de_genes$gene

# filling the matrix
for(i in 1:dim(de_genes)[1]){
  
  # retrieving the logFC for the cell type
  de_gene <- de_genes$gene[i]
  de_gene_cell_type <- de_genes$cell_type[i]
  idx <- which(de_res[[de_gene_cell_type]]$gene == de_gene)[1]
  logFC_de_gene <- de_res[[de_gene_cell_type]]$log2FoldChange[idx]
  to_plot[i, 1] <- logFC_de_gene
  
  # retrieving the logFC for each cluster
  for(j in 1:(length(clusters_de_res))){
    
    # finding the logFC
    current_res <- clusters_de_res[[j]]
    if(is.null(current_res)){ #usable results?
      idx <- NULL
    }else{
      idx <- which(current_res$gene == de_gene)[1]
    }
    current_logFC <- 0
    if(length(idx) > 0){
      if(!is.na(current_res$padj[idx]) &&
         current_res$padj[idx] <= padj_threshold){
        current_logFC <- current_res$log2FoldChange[idx]
      }
    }
    
    # storing in the matrix
    to_plot[i, j + 1] <- current_logFC
    
  }
}

# NAs to 0
to_plot[is.na(to_plot)] <- 0

# removing underscores from the names of the heatmap
rownames(to_plot) <- gsub('_', ' ', rownames(to_plot), fixed = TRUE)
colnames(to_plot) <- gsub('_', ' ', colnames(to_plot), fixed = TRUE)

# keeping only the cluster results that are usable
# " + 1" because we have the whole cells results in the first column
to_plot <- to_plot[, c(1, usable_de_res + 1)] 

# splitting vectors
row_splitting_vector <- gsub('_', ' ', de_genes$cell_type, fixed = TRUE)
col_splitting_vector <- c('.')
for(i in 1:length(cell_type_list)){
  col_splitting_vector <- c(col_splitting_vector, 
                            rep(names(cell_type_list)[i], 
                                length(cell_type_list[[i]])))
}
col_splitting_vector <- gsub('_', ' ', col_splitting_vector)

# keeping only the usable splitting vector elements 
col_splitting_vector <- col_splitting_vector[c(1, usable_de_res + 1)] 

# annotations
row_ha = rowAnnotation(`Cell type` = row_splitting_vector, 
                       show_annotation_name = FALSE, 
                       show_legend = TRUE, 
                       col = list(`Cell type` = cell_type_colors))
column_ha = HeatmapAnnotation(`Cell type` = col_splitting_vector, 
                              show_annotation_name = FALSE, 
                              show_legend = FALSE,
                              col = list(`Cell type` = cell_type_colors))

# selecting cell types
to_plot <- to_plot[row_splitting_vector == 'Glutamatergic neurons', 
                   col_splitting_vector == 'Glutamatergic neurons' | 
                     colnames(to_plot) ==  'Whole cell type']
colnames(to_plot)[1] <- 'Glut. cells'

# plotting
p <- Heatmap(t(to_plot), name = "log2(FC)", border = TRUE,
             col = colorRamp2(c(0, max(max(to_plot))), c("white", "red")), 
             #show_row_dend = FALSE, show_column_dend = FALSE,  
             cluster_rows = FALSE, cluster_columns = FALSE,
             column_names_rot = 50, row_names_side = "left",
             heatmap_legend_param = list(direction = 'horizontal'))
png(filename = file.path(res_folder, 'panel_d.png'), 
    width = 3000, height = 3600, res = 600)
plot(p)
dev.off()

