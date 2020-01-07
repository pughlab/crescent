### Shared by Samah El Ghamrasni <samah.elghamrasni@uhnresearch.ca>
### Jan 7th, 2020 to Javier and Suluxan
 

bc.data.p2 <- Read10X(data.dir = "/Volumes/Samwise/projects/ATAC-SEQ/RNA-SEQ/Pugh_Samah__TNBC_xeno_human-mouse/raw_gene_bc_matrices/hg19/")

bc.TNBC.p2 <- CreateSeuratObject(raw.data = bc.data.p2, min.cells = 0, min.genes = 0, 
                              project = "10X_bc.TNBC.p2")
metadata.hg19 <- bc.TNBC.p2@meta.data
metadata.hg19 <- metadata.hg19[order(metadata.hg19$nUMI
                                     , decreasing = TRUE),]
library(dropbead)
library(data.table)
plotCumulativeFractionOfReads(metadata.hg19,cutoff = 20000)

x <- metadata.hg19[1:1139,]
b<-rownames(x)

subset.matrix <- bc.TNBC.p2@raw.data[,b] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
bc.TNBC.p2 <- CreateSeuratObject(subset.matrix)
