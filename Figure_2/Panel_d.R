
# Script for Figure 2 d

# memory and library
rm(list = ls())
source('../ancillary/libraries.R')

# control panel 
fc_threshold <- 1
padj_threshold <- 0.05
min_base_mean <- 0
num_genes_to_plot <- 25
min_num_cells <- 50
de_folder <- '../1_differential_analysis'
res_folder <- './Panel_d'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# choosing the analysis
de_res_names <- c('Cell_type-Glutamatergic_neurons', 
                  paste0('Cell_subtype-', 
                         c(paste0('EXC_GLU-', 1:7), 'EXC_IMN')))
num_res <- length(de_res_names)

# loading the DE results (cell type)
de_res <- vector('list', num_res)
names(de_res) <- de_res_names
for(i in 1:num_res){
  
  # enough samples?
  num_cells <-  t(read.table(file = file.path(de_folder, de_res_names[i], 
                                              'num_cells_per_sample.txt')))
  if(min(as.numeric(num_cells[, 2])) < min_num_cells){
    next()
  }
  
  # storing results
  de_res[[i]] <- read.csv(file.path(de_folder, de_res_names[i], 
                                    paste0(de_res_names[i], '_results.csv')))
  
}

# keeping only the successful de res
to_keep <- !sapply(de_res, is.null)
de_res <- de_res[to_keep]
num_res <- length(de_res)

# compiling the list of num_genes_to_plot DE genes with highest logFC
tmp <- de_res$`Cell_type-Glutamatergic_neurons`
tmp <- tmp[tmp$log2FoldChange >= fc_threshold & tmp$padj <= padj_threshold & 
             tmp$baseMean >= min_base_mean, ]
tmp <- tmp[!duplicated(tmp$gene), ]
tmp <- tmp[order(abs(tmp$log2FoldChange), decreasing = TRUE), ]
de_genes <- tmp$gene[1:num_genes_to_plot]

# matrix to plot
to_plot <- matrix(NA, num_genes_to_plot, num_res)
colnames(to_plot) <- gsub('Cell_subtype-', '', names(de_res))
colnames(to_plot) <- gsub('_', ' ', colnames(to_plot))
colnames(to_plot)[1] <- 'Glut. cells'
rownames(to_plot) <- de_genes

# retrieving the logFC for the cell type
current_res <- de_res$`Cell_type-Glutamatergic_neurons`
idx <- which(current_res$gene %in% de_genes)
logFC_de_gene <- current_res$log2FoldChange[idx]
names(logFC_de_gene) <- de_genes
logFC_de_gene <- logFC_de_gene[de_genes]
to_plot[de_genes, 1] <- logFC_de_gene

# retrieving the logFC for each cluster
for(j in 2:num_res){
  
  current_res <- de_res[[j]]
  idx <- which(current_res$gene %in% de_genes & 
                 current_res$log2FoldChange >= fc_threshold & 
                 current_res$padj <= padj_threshold & 
                 current_res$baseMean >= min_base_mean)
  logFC_de_gene <- current_res$log2FoldChange[idx]
  current_de_genes <- current_res$gene[idx]
  names(logFC_de_gene) <- current_de_genes
  logFC_de_gene <- logFC_de_gene[current_de_genes]
  to_plot[current_de_genes, j] <- logFC_de_gene
  
  
}

# NAs to 0
to_plot[is.na(to_plot)] <- 0

# plotting
p <- Heatmap(to_plot, name = "log2(FC)", border = TRUE,
             col = colorRamp2(c(0, max(max(to_plot))), c("white", "red")), 
             cluster_rows = FALSE, cluster_columns = FALSE,
             column_names_rot = 50, row_names_side = "left",
             heatmap_legend_param = list(direction = 'vertical')) 

png(filename = file.path(res_folder, 'panel_d.png'), 
    width = 3000, height = 3600, res = 600)
plot(p)
dev.off()

