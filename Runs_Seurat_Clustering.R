####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on http://satijalab.org/seurat/pbmc3k_tutorial.html
### Things missing from this tutorial:
### 1) Further subdivisions within cell types (i.e. granularity of clusters)
### 2) Assigning cell type identity to clusters (needs supervised annotations, maybe based on Gene Set Enrichment analysis)
### 3) Using saveRDS
### 4) Use knitr() to produce better html plot layout (https://yihui.name/knitr/demo/stitch/)
####################################

####################################
### Required libraries
####################################
### 'optparse'   to handle one-line-commands
### 'data.table' to read tables quicker than read.table - only needed is using '-t Dropseq'
### 'Seurat'     to run QC, differential gene expression and clustering analyses
### 'dplyr'      needed by Seurat for data manupulation
####################################

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )


####################################
### Get inputs from command line argumets
####################################

option_list <- list(
  make_option(c("-i", "--input"), default="NA",
              help="Either the path/name to a 10X *directory* with barcodes.tsv, genes.tsv and matrix.mtx files;
                or path/name of a  <tab> delimited Dropseq *file* with cell barcodes in columns and genes in rows"),

  make_option(c("-t", "--input_type"), default="NA",
              help="Indicates if input is either a '10X' directory or a 'Dropseq' matrix file"),
  
  make_option(c("-r", "--resolution"), default="0.8",
              help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain
                a larger (smaller) number of communities"),
  
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh,  e.g. '0.01'"),

  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),

  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID")
)
opt <- parse_args(OptionParser(option_list=option_list))

Input          <- opt$input
InputType      <- opt$input_type
Outdir         <- opt$outdir
ThreshReturn   <- opt$return_threshold
PrefixOutfiles <- c(paste(opt$prefix_outfiles,"_res",opt$resolution,sep=""))
Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

####################################
### Define tailored parameters
####################################
### Some of these parameters are the defaults provided by Seurat developers, others are tailored according to clusters/t-SNE granularity
###
### Parameters for Seurat filters

MinCells<-3
MinGenes<-200
LowThresholds<-c(200,-Inf)
HighThresholds<-c(2500,0.05)
### Parameters for Seurat normalization
ScaleFactor<-10000
### Parameters for Seurat variable gene detection
XLowCutoff<-0.0125
XHighCutoff<-3
YCutoff<-0.5
### Parameters for PCA
PrintPCA.PcsPrint<-1:5
PrintPCA.GenesPrint<-5
VizPCA.PcsUse<-1:6
PCHeatmapCellsUse<-500
PCHeatmapComponentsToPlot<-18
JackStrawNumReplicate<-100
JackStrawPlotPcs<-1:18 ### These are the number of PCs to plot to see their influence data, but won't influence the plots themselves
### Parameters for Clustering
FindClusters.DimsUse<-1:10 ### These should be tailored based on the number of PCs found to influence data (e.g. from JackStrawPlotPcs and/or PCElbowPlot plots)
### Parameters for Cluster Biomarkers
FindAllMarkers.MinPct    <- 0.25
FindAllMarkers.ThreshUse <- 0.25
FindAllMarkers.PrintTopN <- 10
FindMarkers.Pseudocount  <- 1e-99 ### Default is 1, which sounds high for a Log level correction. Also see https://goo.gl/3VzQ3L
NumberOfGenesToPlotFeatures <-16
NumberOfGenesPerClusterToPlotTsne <-2
NumberOfGenesPerClusterToPlotHeatmap <-10

StartTimeOverall<-Sys.time()


####################################
### Create outdirs
####################################

dir.create(file.path(Outdir, "SEURAT"), recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Load data
####################################
if(regexpr("^10X$", InputType, ignore.case = T)[1] == 1) {
  print("Loading 10X infiles")
  input.matrix <- Read10X(data.dir = Input)
}else if (regexpr("^Dropseq$", InputType, ignore.case = T)[1] == 1) {
  print("Loading Drop-seq matrix")
  library(data.table)
  input.matrix <- data.frame(fread(Input),row.names=1)
}else{
  stop(paste("Unexpected type of infile: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}
dim(input.matrix)

####################################
### Create a Seurat object
####################################

seurat.object  <- CreateSeuratObject(raw.data = input.matrix, min.cells = MinCells, min.genes = MinGenes, project = PrefixOutfiles)
seurat.object

####################################
### Get  mitochondrial genes
####################################
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat.object@data), value = T)
percent.mito <- Matrix::colSums(seurat.object@raw.data[mito.genes, ])/Matrix::colSums(seurat.object@raw.data)
seurat.object <- AddMetaData(object = seurat.object, metadata = percent.mito, col.name = "percent.mito")

####################################
### Voilin plots for data UNfiltered by Seurat
####################################
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.pdf", sep=""))
VlnPlot(object = seurat.object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()
### For HTML
VlnPlot(object = seurat.object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

####################################
### Gene scatter plots for data UNfiltered by Seurat
####################################
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_GenePlot.pdf", sep=""),width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "nGene")
dev.off()
### For HTML
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "nGene")

####################################
### Filter cells based gene counts and mitochondrial representation
####################################
seurat.object<-FilterCells(object = seurat.object, subset.names = c("nGene", "percent.mito"), low.thresholds = LowThresholds, high.thresholds = HighThresholds)
seurat.object

####################################
### Voilin plots for data filtered by Seurat
####################################
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.seurat_filtered.pdf", sep=""))
VlnPlot(object = seurat.object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()
### For HTML
VlnPlot(object = seurat.object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

####################################
### Gene scatter plots for data filtered by Seurat
####################################
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_GenePlot.seurat_filtered.pdf", sep=""), width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "nGene")
dev.off()
### For HTML
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "nGene")

####################################
### Normalize data
####################################
seurat.object <- NormalizeData(object = seurat.object, normalization.method = "LogNormalize", scale.factor = ScaleFactor)

####################################
### Detect, save list and plot variable genes
####################################
### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object <- FindVariableGenes(object = seurat.object, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff, do.plot=F)
length(x = seurat.object@var.genes)
write(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VariableGenes.txt", sep=""), x=seurat.object@var.genes)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VariableGenes.pdf", sep=""))
VariableGenePlot(object = seurat.object,  x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff)
dev.off()
### For HTML
VariableGenePlot(object = seurat.object,  x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff)

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
seurat.object <- ScaleData(object = seurat.object, vars.to.regress = c("nUMI", "percent.mito"), display.progress=F)

####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
####################################
seurat.object <- RunPCA(object = seurat.object, pc.genes = seurat.object@var.genes, do.print = T, pcs.print = PrintPCA.PcsPrint, genes.print = PrintPCA.GenesPrint)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
VizPCA(object = seurat.object, pcs.use = VizPCA.PcsUse)
dev.off()
### For HTML
VizPCA(object = seurat.object, pcs.use = VizPCA.PcsUse)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCAPlot.pdf", sep=""))
PCAPlot(object = seurat.object, dim.1 = 1, dim.2 = 2, no.legend=T)
dev.off()
### For HTML
PCAPlot(object = seurat.object, dim.1 = 1, dim.2 = 2, no.legend=T)
#
seurat.object <- ProjectPCA(object = seurat.object, do.print = F)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCHeatmap.C1.pdf", sep=""))
PCHeatmap(object = seurat.object, pc.use = 1, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F)
dev.off()
### For HTML
PCHeatmap(object = seurat.object, pc.use = 1, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCHeatmap.C1to",PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
PCHeatmap(object = seurat.object, pc.use = 1:PCHeatmapComponentsToPlot, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F, use.full = F)
dev.off()
### For HTML
PCHeatmap(object = seurat.object, pc.use = 1:PCHeatmapComponentsToPlot, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F, use.full = F)

####################################
### Determine statistically significant principal components
####################################
### NOTE: This process can take a long time for big datasets, comment out for
### expediency.  More approximate techniques such as those implemented in
### PCElbowPlot() can be used to reduce computation time
seurat.object <- JackStraw(object = seurat.object, num.replicate = JackStrawNumReplicate, display.progress = F)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_JackStraw.C1toC12.pdf", sep=""))
JackStrawPlot(object = seurat.object, PCs = JackStrawPlotPcs)
dev.off()
### For HTML
JackStrawPlot(object = seurat.object, PCs = JackStrawPlotPcs)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCElbowPlot.pdf", sep=""))
PCElbowPlot(object = seurat.object)
dev.off()
### For HTML
PCElbowPlot(object = seurat.object)

####################################
### Cluster the cells
####################################
StartTimeClustering<-Sys.time()
seurat.object <- FindClusters(object = seurat.object, reduction.type = "pca", dims.use = FindClusters.DimsUse, resolution = opt$resolution, print.output = 0, save.SNN = T)
EndTimeClustering<-Sys.time()

CellNames<-seurat.object@cell.names
ClusterIdent<-seurat.object@ident
Headers<-paste("CLUSTERS","seurat_clusters",sep="\t")
clusters_data<-paste(CellNames,ClusterIdent,sep="\t")
OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)


####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
seurat.object <- RunTSNE(object = seurat.object, dims.use = FindClusters.DimsUse, do.fast = T)
### Note that you can set do.label=T to help label individual clusters
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot.pdf", sep=""))
TSNEPlot(object = seurat.object, do.label = T,label.size=10)
dev.off()
### For HTML
TSNEPlot(object = seurat.object, do.label = T,label.size=10)

####################################
### Finding differentially expressed genes (cluster biomarkers)
####################################
### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
#########
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object, ident.1 = 1, min.pct = FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object, ident.1 = 5, ident.2 = c(0, 3), min.pct = FindAllMarkers.MinPct)
### print(x = head(x = cluster5.markers, n = 5))
########
### NOTE: FindAllMarkers() uses return.thresh = 0.01 as defaults, but FindMarkers() displays all genes passing previous filters.
###       Thus to make the outputs between these two commands identical to each other use return.thresh = 1
#########

StartTimeFindAllMarkers<-Sys.time()
seurat.object.markers <- FindAllMarkers(object = seurat.object, only.pos = T, min.pct = FindAllMarkers.MinPct, return.thresh = ThreshReturn, thresh.use = FindAllMarkers.ThreshUse, pseudocount.use=FindMarkers.Pseudocount)
EndTimeFindAllMarkers<-Sys.time()

write.table(data.frame("GENE"=rownames(seurat.object.markers),seurat.object.markers),paste(Tempdir,"/",PrefixOutfiles,".SEURAT_MarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)
### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne<-(seurat.object.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(summary(top_genes_by_cluster_for_tsne[,"cluster"]))

####################################
### Violin plots for top genes are not shown by now
####################################
NumberOfPanesForFeaturesPlot<-(NumberOfClusters*NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot_AfterClusters.pdf", sep=""))
VlnPlot(object = seurat.object, features.plot = c(top_genes_by_cluster_for_tsne.list), use.raw = T, y.log = T, adjust.use=1,point.size.use = 0.5)
dev.off()
### For HTML
VlnPlot(object = seurat.object, features.plot = c(top_genes_by_cluster_for_tsne.list), use.raw = T, y.log = T, adjust.use=1,point.size.use = 0.5)

####################################
### t-SNE plots showing each cluster top genes
####################################
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""))
FeaturePlot(object = seurat.object, features.plot = c(top_genes_by_cluster_for_tsne.list), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
### For HTML
FeaturePlot(object = seurat.object, features.plot = c(top_genes_by_cluster_for_tsne.list), cols.use = c("grey", "blue"), reduction.use = "tsne")

####################################
### Heatmaps
####################################
top_genes_by_cluster_for_heatmap <- seurat.object.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotHeatmap, avg_logFC)
# setting slim.col.label to T will print just the cluster IDS instead of every cell name
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_Heatmap.pdf", sep=""))
DoHeatmap(object = seurat.object, genes.use = top_genes_by_cluster_for_heatmap$gene, slim.col.label = T, remove.key = T, title = PrefixOutfiles, cex.row = 6)
dev.off()
### For HTML
DoHeatmap(object = seurat.object, genes.use = top_genes_by_cluster_for_heatmap$gene, slim.col.label = T, remove.key = T, title = PrefixOutfiles, cex.row = 6)

####################################
### Report used options
####################################
OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

countOptions<-0
for (optionInput in opt) {
  countOptions = countOptions + 1
  write(file = OutfileOptionsUsed, x=paste(names(opt[countOptions]), optionInput, sep="\t"), append = T)
}

####################################
### Report time used
####################################
EndTimeOverall<-Sys.time()

TookTimeClustering     <-format(difftime(EndTimeClustering,     StartTimeClustering,     units = "min"))
TookTimeFindAllMarkers <-format(difftime(EndTimeFindAllMarkers, StartTimeFindAllMarkers, units = "min"))
TookTimeOverall        <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "min"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CPUusage.tsv", sep="")
ReportTime<-c(
  paste("clustering",TookTimeClustering,collapse = "\t"),
  paste("FindAllMarkers",TookTimeFindAllMarkers,collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into ourdir
####################################
outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".SEURAT_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){ 
  file.rename(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/SEURAT/",eachFile,sep=""))
})

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################

print("END - All done!!! Took time:")
print(ReportTime)

quit()

