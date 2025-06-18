
# Script Figure 2b

# memory and library
rm(list = ls())
source('../libraries.R')

# control panel 
fc_threshold <- 1
padj_threshold <- 0.05
num_genes_to_plot <- 5
data_file <- '../data/combined_sets.rds'
pheno_file <- '../data/phenoData.csv'
de_folder <- '../1_differential_analysis/integration_approach_1'
res_folder <- './Panel_b'
dir.create(res_folder, showWarnings = FALSE, recursive = TRUE)

# choosing the analysis
cell_type_list <- list(Gabaergic_neurons = c(26, 22, 34, 27, 15, 28, 14, 31, 19, 29), 
                      Glutamatergic_neurons = c(0:5, 7:10, 12, 13, 16:18, 21, 24, 30, 33),
                      Astrocytes = c(6, 11, 23), Oligodendrocytes = c(20, 25))

# colors!
cell_type_colors <- c(`Gabaergic neurons` = '#A3A400', 
                      `Glutamatergic neurons` = '#F8766C',
                      Astrocytes = '#00BF7C', 
                      Oligodendrocytes = '#00AFF5')

# loading the combined sets and meta data
combined_sets <- readRDS(data_file)
DefaultAssay(combined_sets) <- 'RNA'
Idents(combined_sets) <- 'seurat_clusters'

# loading the meta data
pheno_data <- read.csv(pheno_file)
pheno_data <- pheno_data[pheno_data$X10X_Serial_ID %in% combined_sets$orig.ident, ]
pheno_data <- pheno_data[order(pheno_data$Group), ]

# loading the DE results
de_res <- vector('list', length(cell_type_list))
names(de_res) <- names(cell_type_list)
for(i in 1:length(cell_type_list)){
  de_res[[i]] <- read.csv(file.path(de_folder, names(cell_type_list)[i], 
                                    paste0(names(cell_type_list)[i], '_results.csv')))
}

# compiling the list of de genes
i <- 1
de_genes <- data.frame(cell_type = c(), gene = c())
all_genes <- data.frame(cell_type = c(), gene = c())
for(i in 1:length(cell_type_list)){
  tmp <- de_res[[i]]
  genes <- tmp$gene[abs(tmp$log2FoldChange) >= fc_threshold & tmp$padj <= padj_threshold]
  de_genes <- rbind(de_genes, data.frame(cell_type = names(de_res)[i], gene = genes))
  all_genes <- rbind(all_genes, data.frame(cell_type = names(de_res)[i], gene = tmp$gene))
}

# beautifying
de_genes$cell_type <- gsub('_', ' ', de_genes$cell_type)

# panel b, barplot stacked by lncRNA

# loading gene type information
gtf_selected <- readRDS(file = '../data/gene_biotype_GRCg7b_108.rds')

# merge gene type information with de_genes
de_genes <- merge(de_genes, gtf_selected, 
                  by.x = 'gene', by.y = 'gene_id', all.x = TRUE)
de_genes$gene_biotype[is.na(de_genes$gene_biotype)] <- 'protein_coding'

# summarize information to plot
to_plot <- de_genes %>% group_by(cell_type) %>% 
  summarise(num_total = n(), num_lncRNA = sum(gene_biotype == 'lncRNA'))

# computing p-values vs background information
i <- 1
background <- to_plot
background$pvalue <- NA
for(i in 1:dim(to_plot)[1]){
  
  # computing the total number of lncRNA and transcripts
  current_all_genes <- all_genes[all_genes$cell_type == 
                                   gsub(' ', '_', to_plot$cell_type[i]), ]
  num_transcripts <- dim(current_all_genes)[1]
  idx1 <- which(current_all_genes$gene %in% gtf$gene_id)
  tmp1 <- sum(gtf$gene_id %in% current_all_genes$gene[idx1] & 
                gtf$gene_biotype == 'lncRNA')
  idx2 <- which(current_all_genes$gene %in% gtf$gene_name)
  tmp2 <- sum(gtf$gene_name %in% current_all_genes$gene[idx2] & 
                gtf$gene_biotype == 'lncRNA')
  num_lncRNA <- tmp1 + tmp2
  background$num_total[i] <- num_transcripts
  background$num_lncRNA[i] <- num_lncRNA
  tmp <- matrix(as.numeric(c(to_plot[i, 2:3], num_transcripts, num_lncRNA)), 
                nrow = 2)
  background$pvalue[i] <- chisq.test(tmp)$p.value

}

# beautifying the background
tmp <- to_plot
colnames(tmp) <- c("cell_type",  "num_de",  "num_lncRNA_de")
background <- merge(tmp, background)
background$de_ratio <- background$num_lncRNA_de / background$num_de
background$bg_ratio <- background$num_lncRNA / background$num_total
write.csv(background, file = file.path(res_folder, 'lncRNA_de_statistics.csv'))

# plotting!
p <- ggplot(data = to_plot, 
            mapping = aes(x = reorder(cell_type, -num_total), 
                          y = num_total, fill = cell_type)) + 
  geom_col() + 
  geom_text(mapping = aes(x = reorder(cell_type, -num_total), 
                          y = num_total, label = num_total), 
            color = 'white', size = 3.5, nudge_y = -7) + 
  scale_x_discrete(name = NULL) + 
  scale_y_continuous(name = 'Number of differentially expressed genes') + 
  scale_fill_manual(name = NULL, values = cell_type_colors) + 
  theme_bw() + 
  theme(legend.position = 'none') +
  geom_col(mapping = aes(x = reorder(cell_type, -num_total), 
                         y = num_lncRNA), 
           fill = 'black', alpha = 0.5) + 
  geom_text(mapping = aes(x = reorder(cell_type, -num_total), 
                          y = num_lncRNA, label = num_lncRNA), 
            color = 'white', size = 3.5, nudge_y = -7)

png(filename = file.path(res_folder, 'number_DE_genes_by_cell_type.png'), 
    width = 3300, height = 2100, res = 600)
plot(p)
dev.off()
