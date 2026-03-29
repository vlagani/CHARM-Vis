# loading / installing libraries
if(!require(pak)){
  install.packages('pak', quiet = TRUE)
  library(pak)
}
if(!require(BiocManager)){
  install.packages('BiocManager', quiet = TRUE)
  library(BiocManager)
}
if(!require(Seurat)){
  pak('Seurat', ask = FALSE)
  library(Seurat)
}
if(!require(tidyverse)){
  pak('tidyverse', ask = FALSE)
  library(tidyverse)
}
if(!require(ggrepel)){
  pak('ggrepel', ask = FALSE)
  library(ggrepel)
}
if(!require(openxlsx)){
  pak('openxlsx', ask = FALSE)
  library(openxlsx)
}
if(!require(rtracklayer)){
  pak('rtracklayer', ask = FALSE)
  library(rtracklayer)
}
if(!require(DESeq2)){
  pak('DESeq2', ask = FALSE)
  library(DESeq2)
}
if(!require(future)){
  pak('future', ask = FALSE)
  library(future)
}
if(!require(doFuture)){
  pak('doFuture', ask = FALSE)
  library(doFuture)
}
if(!require(fgsea)){
  BiocManager::install('fgsea', ask = FALSE, update = FALSE)
  library(fgsea)
}
if(!require(DOSE)){
  BiocManager::install('DOSE', ask = FALSE, update = FALSE)
  library(DOSE)
}
if(!require(enrichplot)){
  BiocManager::install('enrichplot', ask = FALSE, update = FALSE)
  library(enrichplot)
}
if(!require(clusterProfiler)){
  pak('clusterProfiler', ask = FALSE)
  library(clusterProfiler)
}
if(!require(enrichplot)){
  pak('enrichplot', ask = FALSE)
  library(enrichplot)
}
if(!require(org.Mm.eg.db)){
  BiocManager::install('org.Mm.eg.db', ask = FALSE, update = FALSE)
  library(org.Mm.eg.db)
}
if(!require(dittoSeq)){
  pak('dittoSeq', ask = FALSE)
  library(dittoSeq)
}
if(!require(SummarizedExperiment)){
  pak('SummarizedExperiment', ask = FALSE)
  library(SummarizedExperiment)
}
if(!require(ComplexHeatmap)){
  pak('ComplexHeatmap', ask = FALSE)
  library(ComplexHeatmap)
}
if(!require(circlize)){
  pak('circlize', ask = FALSE)
  library(circlize)
}
if(!require(simplifyEnrichment)){
  pak('simplifyEnrichment', ask = FALSE)
  library(simplifyEnrichment)
}
