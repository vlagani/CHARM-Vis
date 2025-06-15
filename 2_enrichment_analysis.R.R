# Script for performing enrichemnt analysis 

#### set up ####

# memory and libraries
rm(list = ls())
source('libraries.R')

# control panel
de_folder <- '1_differential_analysis'
min_num_cells <- 50
min_num_genes <- 20
max_num_genes <- 200
significance_threshold <- 0.1
ncores <- 8
res_folder <- '2_enrichment_analysis'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# starting the parallel environment
registerDoFuture()
plan("multisession", workers = ncores)
options(future.globals.maxSize = 20000 * 1024^2)

# how many differential analyses?
analyses <- c('Gabaergic_neurons', 'Glutamatergic_neurons', 
              'Astrocytes', 'Oligodendrocytes')
for(i in 0:35){
  analyses <- c(analyses, paste0('Cluster_', i))
}
n_analyses <- length(analyses)

# loading information for gene translation
source('gene_conversion_function.R')
load(file.path('data', 'annot_table_mouse_chicken_2024_02_09.RData'))

#### enrichment analysis ####

# iterating over analyses
sinking <- foreach(iter = 1:n_analyses, .packages = c('clusterProfiler', 
                                           'enrichplot', 
                                           'org.Gg.eg.db')) %dopar% {

  # specific analysis
  analysis <- analyses[iter]
  print(iter)
  
  # reading the number of cells 
  num_cells <- read.table(header = TRUE, file.path(de_folder, 
                                    analysis, 'num_cells_per_sample.txt'))
  
  # proceeding only if enough cell are measured for each sample
  if(min(num_cells) < min_num_cells){
    return(NULL)
  }
  
  # creating the results folder
  current_res_folder <- file.path(res_folder, analysis)
  dir.create(current_res_folder, showWarnings = FALSE, recursive = TRUE)
  
  # loading the de results
  de_results <- read.csv(file.path(de_folder, analysis, 
                                   paste0(analysis, '_results.csv')))
  
  # translating
  transl_res <- name_conversion(names_to_convert = de_results$gene, 
                                organism_from = 'chicken', 
                                organism_to = 'mouse', 
                                type_from = 'mixed', 
                                type_to = 'ensembl', 
                                annot_table = annot_table)
  de_results$ensembl <- transl_res$converted_names
  
  # preparing the gene list for clusterProfiler
  tmp <- de_results[, c('ensembl', 'stat')] 
  tmp <- unique(tmp)
  tmp <- tmp[!is.na(tmp$ensembl), ]
  tmp <- tmp %>% group_by(ensembl) %>% summarise(stat = mean(stat))
  geneList <- tmp$stat
  names(geneList) <- tmp$ensembl
  geneList <- sort(geneList, decreasing = TRUE)
  
  # # enrichment analysis
  gse_res <- gseGO(geneList = geneList, ont = 'BP', pvalueCutoff = 1,
                   OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL',
                   minGSSize = min_num_genes, maxGSSize = max_num_genes)
  
  # readable for plots
  gse_res <- setReadable(gse_res, 'org.Mm.eg.db', 'ENSEMBL')
  
  # saving...
  saveRDS(gse_res, file = file.path(current_res_folder, paste0(analysis, '_gseGO_res.rds')))
  gse_res_df <- gse_res@result
  write.csv(gse_res_df, file = file.path(current_res_folder,
                                         paste0(analysis, '_gseGO_res_df.csv')))
  
  # simplifying
  gse_res <- simplify(gse_res)

  # saving the filtered results
  saveRDS(gse_res, file = file.path(current_res_folder, 
                                    paste0(analysis, 
                                           '_gseGO_res_filtered.rds')))
  gse_res_df <- gse_res@result
  write.csv(gse_res_df, file = file.path(current_res_folder,
                                         paste0(analysis, 
                                                '_gseGO_res_filtered_df.csv')))
  
  # focusing on up regulated, significant pathways for plotting
  gse_res_up <- gse_res
  gse_res_up@result <- gse_res_up@result[gse_res_up@result$NES > 0, ]
  gse_res_up@result <- gse_res_up@result[gse_res_up@result$p.adjust <= significance_threshold, ]
  
  # any significant pathway?
  if(dim(gse_res_up@result)[1] == 0){
    return(NULL)
  }
  
  # dotplot
  p <- dotplot(gse_res_up)
  png(filename = file.path(current_res_folder, 'gse_res_up_dotplot.png'),
      width = 2100, height = 3000, res = 300)
  plot(p)
  dev.off()

}

# closing cluster
plan(sequential)
