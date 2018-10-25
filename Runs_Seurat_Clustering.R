####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on http://satijalab.org/seurat/pbmc3k_tutorial.html
### Things missing from this tutorial:
### 1) Further subdivisions within cell types (i.e. granularity of clusters)
### 2) Assigning cell type identity to clusters (needs supervised annotations, maybe based on Gene Set Enrichment analysis)
### 3) Using saveRDS
### 4) Use knitr() to produce better html plot layout (https://yihui.name/knitr/demo/stitch/)
###
### Other things missing
### 5) Add a lists of ENSEMBL Ids for mitochondrial genes instead of just MT- and mt-
### 6) To get full path to the User's home to replace '~/' in input paths (if provided in such way) to make the *summary_plots.pdf outfile
###    See 'Create summary plots outfile' part
### 7) Tried to make an option to make optional to -generate_plots but it was conflicting with creating some of the plots. Tried with both:
###    if(regexpr("^y$", GeneratePlots, ignore.case = T)[1] == 1) {
###    if(GeneratePlots == "y") {
####################################

####################################
### Required libraries
####################################
### 'optparse'   to handle one-line-commands
### 'data.table' to read tables quicker than read.table - only needed is using '-t Dropseq'
### 'Seurat'     to run QC, differential gene expression and clustering analyses
### 'dplyr'      needed by Seurat for data manupulation
### 'staplr'     only if using option '-s y', note it needs pdftk
### 'fmsb'       to calculate the percentages of extra properties to be t-SNE plotted
####################################

####################################
### Required external packages
####################################
### 'pdftk'   to merge selected *pdf files into *summary_plots.pdf using libary(staplr)
###           and to add statistics to violin plots overlapping pdf files
###           in Mac you can install it with something like:
###           'sudo port install pdftk'
####################################

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fmsb))

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

####################################
### Get inputs from command line argumets
####################################
#
option_list <- list(
  make_option(c("-i", "--input"), default="NA",
              help="Either the path/name to a 10X *directory* with barcodes.tsv, genes.tsv and matrix.mtx files;
                or path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns"),
#
  make_option(c("-t", "--input_type"), default="NA",
              help="Indicates if input is either a '10X' directory or a 'DGE' file"),
#
  make_option(c("-r", "--resolution"), default="1",
              help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain
                a larger (smaller) number of communities"),
#
  make_option(c("-o", "--outdir"), default="SEURAT_OUTPUTS",
              help="A path/name for the results directory"),
#
  make_option(c("-p", "--prefix_outfiles"), default="your_sample",
              help="A prefix for outfile names, e.g. your project ID"),
#
  make_option(c("-s", "--summary_plots"), default="y",
              help="Indicates if a *summary_plots.pdf file should be generated [use 'y'] or not [use 'n']
                Note this needs 'pdftk' and R library(staplr)"),
#
  make_option(c("-c", "--infile_colour_tsne_discrete"), default="NA",
              help="A <tab> delimited table of barcodes and discrete properties to colour the t-SNE, like:
                Barcode              CellClass    InOtherDatasets
                AAACCTGAGCGGCTTC-1   1            yes
                AAACCTGAGTCGAGTG-1   1            no
                AAACCTGCAAAGGAAG-1   2            yes
                AAACCTGGTCTCATCC-1   2            no"),
#
  make_option(c("-g", "--list_genes"), default="NA",
              help="A <comma> delimited list of gene identifiers whose expression will be mapped into the t-SNE plots"),
#
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check these two files:
                *JackStraw*pdf, use the number of PC's where the solid curve shows a plateau along the dotted line, and
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers"),
#
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh,  e.g. '0.01'")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input          <- opt$input
InputType      <- opt$input_type
Resolution     <- as.numeric(opt$resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir         <- opt$outdir
PrefixOutfiles <- opt$prefix_outfiles
SummaryPlots   <- opt$summary_plots
ListGenes      <- opt$list_genes
ColourTsne     <- opt$infile_colour_tsne
PcaDimsUse     <- c(1:as.numeric(opt$pca_dimensions))
ThreshReturn   <- as.numeric(opt$return_threshold)

PrefixOutfiles <- c(paste(PrefixOutfiles,"_res",Resolution,sep=""))
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
PCHeatmapCellsUse<-300
PCHeatmapComponentsToPlot<-18
JackStrawNumReplicate<-100
JackStrawPlotPcs<-1:18 ### These are the number of PCs to plot to see their influence data, but won't influence the plots themselves
### Parameters for Cluster Biomarkers
FindAllMarkers.MinPct    <- 0.25
FindAllMarkers.ThreshUse <- 0.25
FindAllMarkers.PrintTopN <- 10
NumberOfGenesToPlotFeatures <- 16
NumberOfGenesPerClusterToPlotTsne <- 2
NumberOfGenesPerClusterToPlotHeatmap <- 10
### Parameters for Violin plots of top biomarkers
VlnPlotSizeTitle <- 10 # 10 is good for ENSEMBL ID's in a 4 column matrix-style plot
### Parameters for t-SNE plots
BasePlotSizeTsneSelectedGenes<-14
BasePlotSizeTsneExtraProperties<-14
MaxNumberOfPlotsPerRowInOutfileTsneSelectedGenes<-4
MaxNumberOfPlotsPerRowInOutfileTsneExtraProperties<-2

StartTimeOverall<-Sys.time()

####################################
### Create outdirs
####################################

CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
#
dir.create(file.path(Outdir, "SEURAT"), recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Load scRNA-seq data
####################################
if(regexpr("^10X$", InputType, ignore.case = T)[1] == 1) {
  print("Loading 10X infiles")
  input.matrix <- Read10X(data.dir = Input)
}else if (regexpr("^DGE$", InputType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  library(data.table)
  input.matrix <- data.frame(fread(Input),row.names=1)
}else{
  stop(paste("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
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

mitoRegExpressions<- paste(c("^MT", "^mt"),collapse = "|")
mito.genes <- grep(pattern = mitoRegExpressions, x = rownames(x = seurat.object@data), value = T)
percent.mito <- Matrix::colSums(seurat.object@raw.data[mito.genes, ])/Matrix::colSums(seurat.object@raw.data)
seurat.object <- AddMetaData(object = seurat.object, metadata = percent.mito, col.name = "percent.mito")

####################################
### Voilin plots for data UNfiltered by Seurat
####################################
VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

### Get mean and median
nGene_mean<-(mean(seurat.object@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object@meta.data[,"percent.mito"]))

### Make a file with mean and median
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.B.pdf", sep="")
pdf(file=VlnPlotPdfB, width = 7, height = 7)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                  sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,2),"\n",
                   "median=",round(percent_mito_median,2)),
                 sep = "", collapse="")
mtext(nGeneStats,at = 0.15, cex = 1, col = "blue")
mtext(nUMIStats,at = 0.6, cex = 1, col = "blue")
mtext(mitoStats,at = 1,    cex = 1, col = "blue")
dev.off()

### Merge VlnPlot, and mean and median, pdf files
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)

####################################
### Gene scatter plots for data UNfiltered by Seurat
####################################
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_GenePlot.pdf", sep=""),width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "nGene")
dev.off()

####################################
### Filter cells based gene counts and mitochondrial representation
####################################
seurat.object<-FilterCells(object = seurat.object, subset.names = c("nGene", "percent.mito"), low.thresholds = LowThresholds, high.thresholds = HighThresholds)
seurat.object

####################################
### Voilin plots for data filtered by Seurat
####################################
VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.seurat_filtered.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

### Get mean and median
nGene_mean<-(mean(seurat.object@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object@meta.data[,"percent.mito"]))

### Make a file with mean and median
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.seurat_filtered.B.pdf", sep="")
pdf(file=VlnPlotPdfB, width = 7, height = 7)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                 sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,2),"\n",
                   "median=",round(percent_mito_median,2)),
                 sep = "", collapse="")
mtext(nGeneStats,at = 0.15, cex = 1, col = "blue")
mtext(nUMIStats,at = 0.6, cex = 1, col = "blue")
mtext(mitoStats,at = 1,    cex = 1, col = "blue")
dev.off()

### Merge VlnPlot, and mean and median, pdf files
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.seurat_filtered.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)

####################################
### Gene scatter plots for data filtered by Seurat
####################################
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_GenePlot.seurat_filtered.pdf", sep=""), width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object, gene1 = "nUMI", gene2 = "nGene")
dev.off()

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

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VariableGenes.pdf", sep=""))
VariableGenePlot(object = seurat.object,  x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff)
dev.off()

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

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
VizPCA(object = seurat.object, pcs.use = VizPCA.PcsUse)
dev.off()

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCAPlot.pdf", sep=""))
PCAPlot(object = seurat.object, dim.1 = 1, dim.2 = 2, no.legend=T)
dev.off()

seurat.object <- ProjectPCA(object = seurat.object, do.print = F)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCHeatmap.C1.pdf", sep=""))
PCHeatmap(object = seurat.object, pc.use = 1, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F)
dev.off()

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCHeatmap.C1to",PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
PCHeatmap(object = seurat.object, pc.use = 1:PCHeatmapComponentsToPlot, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F, use.full = F)
dev.off()

####################################
### Determine statistically significant principal components
####################################
### NOTE: This process can take a long time for big datasets, comment out for
### expediency.  More approximate techniques such as those implemented in
### PCElbowPlot() can be used to reduce computation time

# seurat.object <- JackStraw(object = seurat.object, num.replicate = JackStrawNumReplicate, display.progress = F)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_JackStraw.C1toC12.pdf", sep=""))
# JackStrawPlot(object = seurat.object, PCs = JackStrawPlotPcs)
dev.off()

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCElbowPlot.pdf", sep=""))
# PCElbowPlot(object = seurat.object)
dev.off()

####################################
### Cluster the cells
####################################
StartTimeClustering<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object <- FindClusters(object = seurat.object, reduction.type = "pca", dims.use = PcaDimsUse, resolution = Resolution, print.output = 0, save.SNN = T)
EndTimeClustering<-Sys.time()

CellNames<-seurat.object@cell.names
ClusterIdent<-seurat.object@ident
Headers<-paste("CLUSTERS","seurat_clusters",sep="\t")
clusters_data<-paste(CellNames,ClusterIdent,sep="\t")
#
OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)

####################################
### Get average expression for each cluster for each gene
####################################
cluster.averages<-AverageExpression(object = seurat.object, use.raw = T)
OutfileClusterAverages<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
seurat.object <- RunTSNE(object = seurat.object, dims.use = PcaDimsUse, do.fast = T)

### Note that you can set do.label=T to help label individual clusters
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot.pdf", sep=""))
TSNEPlot(object = seurat.object, do.label = T,label.size=10)
dev.off()

####################################
### Colour t-SNE by nGene, nUMI, and percent.mito
####################################

CellPropertiesToTsne<-c("nGene", "nUMI", "percent.mito")
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 15, height = 5)
FeaturePlot(object = seurat.object, features.plot = CellPropertiesToTsne, cols.use = c("grey", "blue"), reduction.use = "tsne", nCol = 3, pt.size = 1.5)
dev.off()

####################################
### Colour t-SNE by -infile_colour_tsne_discrete
####################################

if (ColourTsne == "NA") {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  seurat.object.meta.data<-seurat.object@meta.data
  ExtraCellProperties <- data.frame(read.table(ColourTsne, header = T, row.names = 1))
  
  # This is because Seurat removes the last dash-digit from barcode ID's
  # and we need to match those ID's against the inputted -infile_colour_tsne_discrete
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object <- AddMetaData(object = seurat.object, metadata = ExtraCellProperties)

  # Generating outfile
  # Note TSNEPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...))
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""))
  for (property in colnames(ExtraCellProperties)) {
    TSNEPlot(object = seurat.object, group.by = property, plot.title = property)
  }
  dev.off()
}

####################################
### t-SNE plots showing each requested gene
####################################

if (ListGenes == "NA") {
  print("No selected genes for t-SNE plots")
}else{
  ListOfGenesForTsnes<-unlist(strsplit(ListGenes, ","))
  ListOfGenesForTsnes
  if (length(ListOfGenesForTsnes) < MaxNumberOfPlotsPerRowInOutfileTsneSelectedGenes) {
    pdfWidth<-(BasePlotSizeTsneSelectedGenes / MaxNumberOfPlotsPerRowInOutfileTsneSelectedGenes) * length(ListOfGenesForTsnes)
    pdfHeight<-BasePlotSizeTsneSelectedGenes / MaxNumberOfPlotsPerRowInOutfileTsneSelectedGenes
  }else{
    NumberOfRowsInPlot<-as.integer(length(ListOfGenesForTsnes) / MaxNumberOfPlotsPerRowInOutfileTsneSelectedGenes)
    pdfWidth<-BasePlotSizeTsneSelectedGenes
    pdfHeight<-NumberOfRowsInPlot * BasePlotSizeTsneSelectedGenes
  }
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  FeaturePlot(object = seurat.object, features.plot = c(ListOfGenesForTsnes), cols.use = c(rgb(red = 0.9, green = 0.9, blue = 0.9, alpha = 0.1), rgb(red = 0, green = 0, blue = 1, alpha = 0.01)), reduction.use = "tsne")
  dev.off()
}

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
### NOTES:
### 1) FindAllMarkers() uses return.thresh = 0.01 as defaults, but FindMarkers() displays all genes passing previous filters.
###    Thus to make the outputs between these two commands identical to each other use return.thresh = 1
###
### 2) Default pseudocount.use=1, which sounds high for a Log level correction
###    An earlier version of this script was using 1e-99, but it was probably too small
###    Now using the inverse of the number of cells in the data.
###    This is sufficiently small as to not compress logGER magnitudes,
###    while keeping comparisons with zero reasonably close to the range of potential logGER values (Innes and Bader, 2018, F1000 Research)
#########

StartTimeFindAllMarkers<-Sys.time()

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object@meta.data))
seurat.object.markers <- FindAllMarkers(object = seurat.object, only.pos = T, min.pct = FindAllMarkers.MinPct, return.thresh = ThreshReturn, thresh.use = FindAllMarkers.ThreshUse, pseudocount.use=FindMarkers.Pseudocount)
EndTimeFindAllMarkers<-Sys.time()

write.table(data.frame("GENE"=rownames(seurat.object.markers),seurat.object.markers),paste(Tempdir,"/",PrefixOutfiles,".SEURAT_MarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)
### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne<-(seurat.object.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers[["cluster"]]))

####################################
### Violin plots for top genes
####################################
NumberOfPanesForFeaturesPlot<-(NumberOfClusters*NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfWidth<-7
pdfHeight<-NumberOfClusters
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
VlnPlot(object = seurat.object, features.plot = c(top_genes_by_cluster_for_tsne.list), use.raw = T, y.log = T, adjust.use=1,point.size.use = 0.5,size.title.use=VlnPlotSizeTitle)
dev.off()

####################################
### t-SNE plots showing each cluster top genes
####################################
pdfWidth<-7
pdfHeight<-NumberOfClusters
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
FeaturePlot(object = seurat.object, features.plot = c(top_genes_by_cluster_for_tsne.list), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()

####################################
### Heatmaps
####################################
top_genes_by_cluster_for_heatmap <- seurat.object.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotHeatmap, avg_logFC)
# setting slim.col.label to T will print just the cluster IDS instead of every cell name
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_Heatmap.pdf", sep=""))
DoHeatmap(object = seurat.object, genes.use = top_genes_by_cluster_for_heatmap$gene, slim.col.label = T, remove.key = T, title = PrefixOutfiles, cex.row = 6)
dev.off()

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

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CPUusage.txt", sep="")
ReportTime<-c(
  paste("clustering",TookTimeClustering,collapse = "\t"),
  paste("FindAllMarkers",TookTimeFindAllMarkers,collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Create summary plots outfile
####################################
if (SummaryPlots == "y") {
suppressPackageStartupMessages(library(staplr))
SummaryPlotsPdf<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_summary_plots.pdf", sep = "", collapse = "")
ListOfPdfFilesToMerge<-c(paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_VlnPlot.pdf", sep = "", collapse = ""),
                         paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_Heatmap.pdf", sep = "", collapse = ""),
                         paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
                         paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
                         paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
                         )
  if (ListGenes == "NA") {
    print("Create summary file")
  }else{
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    print("Create summary file including t-SNE's for selected genes")
  }
staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge, output_filepath = SummaryPlotsPdf)
}

####################################
### Moving outfiles into ourdir
####################################
outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".SEURAT_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/SEURAT/",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
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

