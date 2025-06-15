# Script for performing DE analysis

#### set up ####

# memory and libraries
rm(list = ls())
source('libraries.R')
options(future.globals.maxSize = 2000 * 1024^2)

# control panel
min_baseMean <- 10
num_processors <- 8
num_genes_pca <- 500
factor_to_investigate <- c('Group', 'Sex')
res_folder <- '1_differential_analysis'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# choosing the analysis
clusters_list <- list(Gabaergic_neurons = c(26, 22, 34, 27, 15, 28, 14, 31, 19, 29),
                      Glutamatergic_neurons = c(0:5, 7:10, 12, 13, 16:18, 21, 24, 30, 33),
                      Astrocytes = c(6, 11, 23), Oligodendrocytes = c(20, 25))
for(i in 0:35){ # each single cluster
  clusters_list[[paste0('Cluster_', i)]] <- i
}
analysis_names <- names(clusters_list)
n_analyses <- length(analysis_names)

# creating the parallel environment 
plan("multisession", workers = num_processors)

# loading the combined sets
combined_sets <- readRDS(file.path('data', 'combined_sets.rds'))

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(combined_sets) <- "RNA"

# loading phenotype data
phenoData <- read.csv(file.path('data', 'phenoData.csv'), check.names = FALSE)
rownames(phenoData) <- phenoData$`10X_Serial_ID`
phenoData <- phenoData[phenoData$`10X_Serial_ID` %in% unique(combined_sets$orig.ident), ]

# loading gene information translation
gene_ids_table <- read.csv(file.path('data', 'gene_ids_biomart_2024_02_25.csv'))
names(gene_ids_table) <- c('ensembl', 'symbol', 'entrez')

# sex genes
sex_genes <- readRDS(file = file.path('data', 'sex_genes_GRCg7b_108.rds'))

# iterating over analyses
for(iter in 1:n_analyses){
  
  # specific analysis
  name_analysis <- analysis_names[iter]
  clusters <- clusters_list[[iter]]
  
  # creating the results folder
  current_res_folder <- file.path(res_folder, name_analysis)
  dir.create(current_res_folder, showWarnings = FALSE, recursive = TRUE)
  
  # extracting the desired cell type
  Idents(combined_sets) <- 'seurat_clusters'
  selected_cells <- subset(combined_sets, idents = clusters)
  
  # writing the number of cells for each sample
  sink(file.path(current_res_folder, 'num_cells_per_sample.txt'))
  print(table(selected_cells$orig.ident))
  sink()
  
  #### Pseudo bulk ####
  
  # pseudo counts!
  tmp <- selected_cells@assays$RNA@counts
  meta_data <- selected_cells@meta.data
  pseudo <- matrix(NA, nrow = nrow(tmp), ncol = nrow(phenoData))
  colnames(pseudo) <- phenoData$`10X_Serial_ID`
  rownames(pseudo) <- rownames(tmp)
  for(i in 1:dim(phenoData)[1]){
    pseudo[ , phenoData$`10X_Serial_ID`[i]] <- 
      apply(tmp[ , meta_data$orig.ident == phenoData$`10X_Serial_ID`[i], 
                 drop = FALSE], 1, sum)
  }
  
  # DESeq2 object
  phenoData$Group <- factor(phenoData$Group)
  phenoData$Group <- relevel(phenoData$Group, 'Untrained')
  phenoData$Sex <- factor(phenoData$Sex)
  phenoData$Sex <- relevel(phenoData$Sex, 'F')
  dds <- DESeqDataSetFromMatrix(pseudo, phenoData, ~ Group + Sex)
  
  # DESeq2 analysis
  dds <- tryCatch(DESeq(dds), error = function(e){
    dds <- DESeqDataSetFromMatrix(pseudo + 1, phenoData, ~ Group + Sex)
    DESeq(dds)
  })
  
  # saving pseudo counts
  saveRDS(dds, file = file.path(current_res_folder, 'pseudoBulk_dds.RData'))
  
  # log fold change shrinking
  res <- results(dds)
  
  # results
  resTable <- as.data.frame(res)
  
  # removing genes with NA p-values
  resTable <- resTable[!is.na(resTable$pvalue), ]
  
  # removing sex genes
  resTable <- resTable[!(rownames(resTable) %in% sex_genes), ]
  
  # removing genes with low average expression
  resTable <- resTable[resTable$baseMean >= min_baseMean,  ]
  
  # adjusting the p-values again
  resTable$padj <- p.adjust(resTable$pvalue, method = 'fdr')
  
  # attaching gene names
  resTable <- cbind(gene = rownames(resTable), resTable)

  # add ensembl ids for all the genes...
  resTable <- merge(resTable, gene_ids_table, 
                    by.x = 'gene', by.y = 'symbol', 
                    all.x = TRUE)
  resTable$entrez <- NULL
  idx <- which(is.na(resTable$ensembl) & grepl('ENSGALG', resTable$gene))
  resTable$ensembl[idx] <- resTable$gene[idx]
  
  # add symbol when available...
  resTable <- merge(resTable, gene_ids_table, 
                    by.x = 'ensembl', by.y = 'ensembl', 
                    all.x = TRUE)
  idx <- which(grepl('ENSGALG', resTable$gene) & (resTable$symbol != '' & !is.na(resTable$symbol)))
  resTable$gene[idx] <- resTable$symbol[idx]
  resTable$symbol <- NULL
  resTable$entrez <- NULL

  # ordering
  resTable <- resTable[order(resTable$pvalue), ]
  
  # only unique values
  resTable <- unique(resTable)
  
  # printing
  fileName <- file.path(current_res_folder, paste0(name_analysis, '_results.csv'))
  write.csv(resTable, file = fileName, row.names = FALSE)
  
  # extracting log counts
  tmp <- assay(rlog(dds))
  
  # PCA plot
  sds <- apply(tmp, 1, sd)
  toKeep <- which(rank(-sds) <= num_genes_pca)
  tmp <- tmp[toKeep, ]
  tmp <- t(tmp)
  pcres <- prcomp(tmp)
  evar <- pcres$sdev^2 / sum(pcres$sdev^2)
  toPlot <- data.frame(PC1 = pcres$x[, 1], 
                       PC2 = pcres$x[, 2], 
                       Group = phenoData$Group, 
                       Sex = phenoData$Sex, 
                       sample = phenoData$`10X_Serial_ID`)
  p <- ggplot(toPlot, aes(x = PC1, y = PC2, 
                          color = Group, 
                          shape = Sex,
                          label = sample)) + 
    geom_point() + 
    geom_text_repel() + 
    scale_x_continuous(name = paste0('PC1, var: ', round(evar[1], 3) * 100, '%')) + 
    scale_y_continuous(name = paste0('PC2, var: ', round(evar[2], 3) * 100, '%')) +
    theme_bw()
  png(filename = file.path(current_res_folder, 'PCA.png'), 
      width = 2100, height = 1800, res = 300)
  plot(p)
  dev.off()
  
}

# closing cluster
plan(sequential)
