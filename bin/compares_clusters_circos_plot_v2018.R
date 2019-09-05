####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://xgaoo.github.io/ClusterMap/ClusterMap.html
###
### IMPORTANT: this script needs library(Seurat) version 2
###            because library(ClusterMap) needs to be updated to use Seurat v3
###            If needed, one can install the two versions of Seurat (v2 and v3) in different paths,
###            as shown in script '~/r_programs/installs_two_versions_of_seurat.R'
###            and load the v2 like:
###            PathForV2Libs<-paste(.libPaths(), "/Seurat_v2", sep = "", collapse = "")
###            library("Seurat",lib.loc = PathForV2Libs)
###
### NOTE: If the MTX inputs come from Cell Ranger v2 (genes.tsv, matrix.tsv, barcodes.tsv) they will enter directly to this script
###       If the MTX inputs come from Cell Ranger v3 (features.tsv.gv, matrix.tsv.gz, barcodes.tsv.gz) they will first be converted to v2 format
###
### THINGS MISSING / TO DO
### 1) In "Load scRNA-seq data" we use Seurat's Read10X() function. In this command, when all barcodes come from the same sample (i.e. finish with the same digit), like:
###    CTCTACGCAAGAGGCT-1
###    CTGAAACCAAGAGGCT-1
###    CTGAAACCAAGAGGCT-1
###    ... etc
###    Read10X will remove the '-digit'
###
###    Hence we need to implement code to remove the '-digit' from --input_clusters barcode ID's as well WHEN all barcodes come from the same sample
###    For now, this script is removing the digit always. And user must provide the inputs like:
###    1-CTCTACGCAAGAGGCT
###    2-CTCGAAAAGCTAACAA
###    3-CTGCCTAGTGCAGGTA
###    When the same barcode appears in multiple samples, instead of:
###    CTCTACGCAAGAGGCT-1
###    CTCGAAAAGCTAACAA-2
###    CTGCCTAGTGCAGGTA-3
###
### 2) To implement loading preexisting Seurat object *RDS files
###   
####################################

####################################
### Required libraries
####################################
PathForV2Libs<-paste(.libPaths(), "/Seurat_v2", sep = "", collapse = "")
suppressPackageStartupMessages(library("Seurat",lib.loc = PathForV2Libs))  # to run QC, differential gene expression and clustering analyses, change path as needed for Seurat v2
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(ClusterMap))   # to draw circos and other plots. Install as: library(devtools), then install_github('xgaoo/ClusterMap')
suppressPackageStartupMessages(library(VennDiagram))  # to get gene and cell Venn diagrams
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table - only needed if using '-t DGE'
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

####################################
### Get inputs from command line argumets
####################################

option_list <- list(
  make_option(c("-i", "--input_1"), default="NA",
              help="Path/name to the first datset to compare
              It can be wither a MTX *directory* with barcodes.tsv, genes.tsv and matrix.mtx files
              or path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns"),
  #
  make_option(c("-j", "--input_2"), default="NA",
              help="Same as -i but indicating the second datset to compare
              NOTE: the two datasets to compare can have different genes, but this program will work on their intersection"),
  #
  make_option(c("-t", "--input_type_1"), default="NA",
              help="Indicates if --input_1 is either a 'MTX' *directory* or a 'DGE' file"),
  #
  make_option(c("-u", "--input_type_2"), default="NA",
              help="Same as -t but for --input_2"),
  #
  make_option(c("-k", "--prefix_input_1"), default="NA",
              help="A prefix for labels referring to input_1, e.g. sample_1_ID"),
  #
  make_option(c("-l", "--prefix_input_2"), default="NA",
              help="A prefix for labels referring to input_2, e.g. sample_2_ID"),
  #
  make_option(c("-q", "--input_1_resolution"), default="1",
              help="Value of the resolution parameter to use for --input_1
              Use a value above (below) 1.0 if you want to obtain
              a larger (smaller) number of communities"),
  #
  make_option(c("-r", "--input_2_resolution"), default="1",
              help="Same as -q but for --input_2"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID"),
  #
  make_option(c("-s", "--summary_plots"), default="y",
              help="Indicates if a *summary_plots.pdf file should be generated [use 'y'] or not [use 'n']
              Note this needs 'pdftk' and R library(staplr)"),
  #
  make_option(c("-x", "--infile_colour_tsne_extra_1"), default="NA",
              help="A <tab> delimited table of barcodes and extra properties to colour the t-SNE of input_1, like:
              Barcode              Extra_property_1   Extra_property_2
              AAACCTGAGCGGCTTC-1   1                  yes
              AAACCTGAGTCGAGTG-1   1                  no
              AAACCTGCAAAGGAAG-1   2                  yes
              AAACCTGGTCTCATCC-1   2                  no"),
  #
  make_option(c("-y", "--infile_colour_tsne_extra_2"), default="NA",
              help="A <tab> delimited table of barcodes and extra properties to colour the t-SNE of input_2,
              similar format than --infile_colour_tsne_extra_1 but for barcodes from input_2"),
  #
  make_option(c("-z", "--infile_colour_tsne_extra_merged"), default="NA",
              help="A <tab> delimited table of barcodes and extra properties to colour the t-SNE of the dataset merging (input_1 and input_2),
              similar format than --infile_colour_tsne_extra_1 but for barcodes IDs should look like:
              Dataset1-AAACCCAAGCTGAAGC
              Dataset2-AAACCCAGTACTCCCC"),
  #
  make_option(c("-g", "--list_genes"), default="NA",
              help="A <comma> delimited list of gene identifiers whose expression will be mapped into the t-SNE plots"),
  #
  make_option(c("-a", "--opacity"), default="0.1",
              help="If using a --list_genes, this parameter provides a value for the minimal opacity of gene expression. Use a value between 0 and 1"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
              FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
              Typically '10' is enough, if unsure use '10' and afterwards check these two files:
              *JackStraw*pdf, use the number of PC's where the solid curve shows a plateau along the dotted line, and
              *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers"),
  #
  make_option(c("-m", "--percent_mito"), default="Inf,0.05",
              help="<comma> delimited min,max number of percentage of mitochondrial gene counts in a cell to be included in normalization and clustering analyses
              Use 'Inf' as min if no minumum limit should be used, e.g. 'Inf,0.05'"),
  #
  make_option(c("-n", "--n_genes"), default="200,8000",
              help="<comma> delimited min,max number of unique gene counts in a cell to be included in normalization and clustering analyses. E.g '200,8000'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh,  e.g. '0.01'"),
  
  make_option(c("-f", "--edge_cutoff"), default="0.1",
              help="The edge length cutoff to decide the sub-nodes to merge or not between two the two datasets
              Also controls the cutoff to draw lines between clusters in the circos plot
              Default is 0.1. If lines are too weak you can lower it e.g. 0.05")
  )

opt <- parse_args(OptionParser(option_list=option_list))

Input1         <- opt$input_1
Input2         <- opt$input_2
InputType1     <- opt$input_type_1
InputType2     <- opt$input_type_2
PrefixInput1   <- opt$prefix_input_1
PrefixInput2   <- opt$prefix_input_2
Resolution1    <- as.numeric(opt$input_1_resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Resolution2    <- as.numeric(opt$input_2_resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir         <- opt$outdir
PrefixOutfiles <- opt$prefix_outfiles
SummaryPlots   <- opt$summary_plots
ListGenes      <- opt$list_genes
ColourTsne1    <- opt$infile_colour_tsne_extra_1
ColourTsne2    <- opt$infile_colour_tsne_extra_2
ColourTsneM    <- opt$infile_colour_tsne_extra_merged
Opacity        <- as.numeric(opt$opacity)
PcaDimsUse     <- c(1:as.numeric(opt$pca_dimensions))
StrNGenes      <- opt$n_genes
StrPmito       <- opt$percent_mito
ThreshReturn   <- as.numeric(opt$return_threshold)
EdgeCutoff     <- as.numeric(opt$edge_cutoff)

Resolution12   <-max(Resolution1,Resolution2)

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

####################################
### Define tailored parameters
####################################
### Some of these parameters are the defaults provided by Seurat developers, others are tailored according to clusters/t-SNE granularity
###
### Parameters for Seurat filters
ListNGenes<-unlist(strsplit(StrNGenes, ","))
MinCells<-3
MinGenes<-as.numeric(ListNGenes[1])
MaxGenes<-as.numeric(ListNGenes[2])
ListPmito<-unlist(strsplit(StrPmito, ","))
MinPercentMito<-"Inf"
MinPercentMito<-gsub(pattern = "Inf",replacement = "-Inf", x=MinPercentMito)
MinPercentMito<-as.numeric(MinPercentMito)
MaxPercentMito<-as.numeric(ListPmito[2])
LowThresholds<-c(MinGenes,MinPercentMito)
HighThresholds<-c(MaxGenes,MaxPercentMito)
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
# Parameters for Colouring datasets
ColorDataset1Tsne<-rgb(0.90,0.60,0,0.3) # orange with strong transparency
ColorDataset2Tsne<-rgb(0.35,0.70,0.90,0.3) # blue with strong transparency
ColorDataset1Venn<-rgb(0.90,0.60,0,0.7) # orange with weak transparency
ColorDataset2Venn<-rgb(0.35,0.70,0.90,0.7) # blue with weak transparency

StartTimeOverall <-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input_1", "input_2", "input_type_1", "input_type_2", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs
####################################
writeLines("\n*** Create outdirs ***\n")
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
Outdir<-gsub("/$", "", Outdir)
Tempdir<-gsub("/$", "", Tempdir)
#
dir.create(file.path(Outdir, "SEURAT"), recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Check if MTX files need to be downgraded to v2 format
####################################
writeLines("\n*** Check scRNA-seq data 1 format ***\n")
if (regexpr("^MTX$", InputType1, ignore.case = T)[1] == 1) {

  ### To look for *gz source files
  OriginFeaturesFileGz <- paste(Input1, "/features.tsv.gz",   sep = "", collapse = "")
  OriginBarcodesFileGz <- paste(Input1, "/barcodes.tsv.gz",   sep = "", collapse = "")
  OriginMatrixFileGz   <- paste(Input1, "/matrix.mtx.gz",     sep = "", collapse = "")

  if (file.exists(OriginFeaturesFileGz) == T) {
    NewInput1 <- paste(Outdir, "/SEURAT/MTX_DATASET_1", sep = "", collapse = "")
    dir.create(file.path(NewInput1), recursive = T)
    
    ### Destination files
    DestBarcodesFileGz   <- paste(NewInput1, "/barcodes.tsv.gz",   sep = "", collapse = "")
    DestMatrixFileGz     <- paste(NewInput1, "/matrix.mtx.gz",     sep = "", collapse = "")
    DestGenesFileGz      <- paste(NewInput1, "/genes.tsv",         sep = "", collapse = "")
    
    file.copy(from=OriginBarcodesFileGz, to=DestBarcodesFileGz,  overwrite=T)
    file.copy(from=OriginMatrixFileGz,   to=DestMatrixFileGz,    overwrite=T)
    system(command = paste("gunzip -f ", DestBarcodesFileGz, sep = "", collapse = ""), wait = T)
    system(command = paste("gunzip -f ", DestMatrixFileGz,   sep = "", collapse = ""), wait = T)
    system(command = paste("zmore ", OriginFeaturesFileGz, " | cut -f 1,2 > ",  DestGenesFileGz, sep = "", collapse = ""))
    
    ### New infiles are downgraded v2 format
    Input1<-NewInput1
    ToRemoveDowngraded1<-1
  }
}

writeLines("\n*** Check scRNA-seq data 2 format ***\n")
if (regexpr("^MTX$", InputType2, ignore.case = T)[1] == 1) {
  
  ### To look for *gz source files
  OriginFeaturesFileGz <- paste(Input2, "/features.tsv.gz",   sep = "", collapse = "")
  OriginBarcodesFileGz <- paste(Input2, "/barcodes.tsv.gz",   sep = "", collapse = "")
  OriginMatrixFileGz   <- paste(Input2, "/matrix.mtx.gz",     sep = "", collapse = "")
  
  if (file.exists(OriginFeaturesFileGz) == T) {
    NewInput2 <- paste(Outdir, "/SEURAT/MTX_DATASET_2", sep = "", collapse = "")
    dir.create(file.path(NewInput2), recursive = T)
    
    ### Destination files
    DestBarcodesFileGz   <- paste(NewInput2, "/barcodes.tsv.gz",   sep = "", collapse = "")
    DestMatrixFileGz     <- paste(NewInput2, "/matrix.mtx.gz",     sep = "", collapse = "")
    DestGenesFileGz      <- paste(NewInput2, "/genes.tsv",         sep = "", collapse = "")
    
    file.copy(from=OriginBarcodesFileGz, to=DestBarcodesFileGz,  overwrite=T)
    file.copy(from=OriginMatrixFileGz,   to=DestMatrixFileGz,    overwrite=T)
    system(command = paste("gunzip -f ", DestBarcodesFileGz, sep = "", collapse = ""), wait = T)
    system(command = paste("gunzip -f ", DestMatrixFileGz,   sep = "", collapse = ""), wait = T)
    system(command = paste("zmore ", OriginFeaturesFileGz, " | cut -f 1,2 > ",  DestGenesFileGz, sep = "", collapse = ""))
    
    ### New infiles are downgraded v2 format
    Input2<-NewInput2
    ToRemoveDowngraded2<-1
  }
}

####################################
### Load scRNA-seq data and Create Seurat objects
####################################
writeLines("\n*** Load scRNA-seq data 1 ***\n")
if (regexpr("^MTX$", InputType1, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  input.matrix_1 <- Read10X(data.dir = Input1)
  seurat.object_1_full  <- CreateSeuratObject(raw.data = input.matrix_1, min.cells = MinCells, min.genes = MinGenes, project = PrefixOutfiles)
}else if (regexpr("^DGE$", InputType1, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  input.matrix_1 <- data.frame(fread(Input1),row.names=1)
  seurat.object_1_full  <- CreateSeuratObject(raw.data = input.matrix_1, min.cells = MinCells, min.genes = MinGenes, project = PrefixOutfiles)
  stop(paste("Unexpected type of input: ", InputType1, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

writeLines("\n*** Load scRNA-seq data 2 ***\n")
if (regexpr("^MTX$", InputType2, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  input.matrix_2 <- Read10X(data.dir = Input2)
  seurat.object_2_full <- CreateSeuratObject(raw.data = input.matrix_2, min.cells = MinCells, min.genes = MinGenes, project = PrefixOutfiles)
}else if (regexpr("^DGE$", InputType2, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  library(data.table)
  input.matrix_2 <- data.frame(fread(Input2),row.names=1)
  seurat.object_2_full <- CreateSeuratObject(raw.data = input.matrix_2, min.cells = MinCells, min.genes = MinGenes, project = PrefixOutfiles)
}else{
  stop(paste("Unexpected type of input: ", InputType2, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

####################################
### Get intersection of genes from the two datasets
####################################
writeLines("\n*** Get intersection of genes from the two datasets ***\n")

listGenes_Union         <- union(row.names(seurat.object_1_full@raw.data),row.names(seurat.object_2_full@raw.data))
listGenes_Intersection  <- intersect(row.names(seurat.object_1_full@raw.data),row.names(seurat.object_2_full@raw.data))
listGenes_Complement1   <- setdiff(row.names(seurat.object_1_full@raw.data),row.names(seurat.object_2_full@raw.data))
listGenes_Complement2   <- setdiff(row.names(seurat.object_2_full@raw.data),row.names(seurat.object_1_full@raw.data))

subset.matrix_1      <- seurat.object_1_full@raw.data[listGenes_Intersection, ]
seurat.object_1      <- CreateSeuratObject(subset.matrix_1)
subset.matrix_2      <- seurat.object_2_full@raw.data[listGenes_Intersection, ]
seurat.object_2      <- CreateSeuratObject(subset.matrix_2)

####################################
### Get intersection of cell barcodes from the two datasets
####################################
writeLines("\n*** Get intersection of cell barcodes from the two datasets ***\n")

listCells_Union         <- union(colnames(seurat.object_1_full@raw.data),colnames(seurat.object_2_full@raw.data))
listCells_Intersection  <- intersect(colnames(seurat.object_1_full@raw.data),colnames(seurat.object_2_full@raw.data))
listCells_Complement1   <- setdiff(colnames(seurat.object_1_full@raw.data),colnames(seurat.object_2_full@raw.data))
listCells_Complement2   <- setdiff(colnames(seurat.object_2_full@raw.data),colnames(seurat.object_1_full@raw.data))

####################################
### Get gene and cell Venn diagrams
####################################
writeLines("\n*** Get gene and cell Venn diagrams ***\n")

### Gene Venn diagram
GeneVennDiagram<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_GenesVennDiagram.pdf", sep="")
pdf(file=GeneVennDiagram, width = 7, height = 7)

plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,bty="n",xlab="",ylab="", main = "Genes in datasets")
title(sub=paste("Union=", length(listGenes_Union), sep = "", collapse = NULL), adj=1, line=3, font=2)
par(new=T)

coloursVenn<-c(ColorDataset1Venn,ColorDataset2Venn)
venn.plot <- draw.pairwise.venn(
  area1 = length(row.names(seurat.object_1_full@raw.data)),
  area2 = length(row.names(seurat.object_2_full@raw.data)),
  cross.area = length(listGenes_Intersection),
  category = c(PrefixInput1, PrefixInput2),
  fill = coloursVenn,
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  cat.col = coloursVenn,
  margin=0.1
)
dev.off()

### Cells Venn diagram
CellsVennDiagram<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_CellsVennDiagram.pdf", sep="")
pdf(file=CellsVennDiagram, width = 7, height = 7)

plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,bty="n",xlab="",ylab="", main = "Cells in datasets")
title(sub=paste("Union=", length(listCells_Union), sep = "", collapse = NULL), adj=1, line=3, font=2)
par(new=T)

coloursVenn<-c(ColorDataset1Venn,ColorDataset2Venn)
venn.plot <- draw.pairwise.venn(
  area1 = length(colnames(seurat.object_1_full@raw.data)),
  area2 = length(colnames(seurat.object_2_full@raw.data)),
  cross.area = length(listCells_Intersection),
  category = c(PrefixInput1, PrefixInput2),
  fill = coloursVenn,
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  cat.col = coloursVenn,
  margin=0.1
)
dev.off()

####################################
### Create tables of genes and cells per dataset
####################################

writeLines("\n*** Create tables of genes and cells per dataset ***\n")

OutfileGenesPerDataset<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_GenesPerDataset.tsv", sep="")
write.table(paste("GeneID", PrefixInput1, PrefixInput2, sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste(listGenes_Intersection,"yes","yes", sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement1,"yes","no"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement2,"no","yes"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)

OutfileGenesPerDataset<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CellsPerDataset.tsv", sep="")
write.table(paste("CellBarcode", PrefixInput1, PrefixInput2, sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste(listGenes_Intersection,"yes","yes", sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement1,"yes","no"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement2,"no","yes"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)

####################################
### Create merged Seurat object
####################################

writeLines("\n*** Get intersection of genes from the two datasets ***\n")

print("Original datasets:")
print("Dataset 1:")
seurat.object_1_full
print("Dataset 2:")
seurat.object_2_full

print("Gene intersection datasets:")
print("Dataset 1:")
seurat.object_1
print("Dataset 2:")
seurat.object_2
nCellsInOriginalMatrix_1<-length(colnames(seurat.object_1@raw.data))
nCellsInOriginalMatrix_2<-length(colnames(seurat.object_2@raw.data))

writeLines("\n*** Create merged Seurat object ***\n")

colnames(input.matrix_1)=paste0('Dataset1-', colnames(input.matrix_1))
colnames(input.matrix_2)=paste0('Dataset2-', colnames(input.matrix_2))
CombDatasets.data = cbind(input.matrix_1, input.matrix_2)
seurat.object_c12  <- CreateSeuratObject(raw.data = CombDatasets.data, min.cells = MinCells, min.genes = MinGenes, project = PrefixOutfiles)
seurat.object_c12
nCellsInOriginalMatrix_c12<-length(colnames(seurat.object_c12@raw.data))

####################################
### Define mitochindrial genes
####################################
writeLines("\n*** Define mitochindrial genes ***\n")

mitoRegExpressions<- paste(c("^MT-", "^mt-"),collapse = "|")


####################################
####################################
####################################
### Process dataset 1
####################################
####################################
####################################
writeLines("\n*** Process dataset 1 ***\n")

writeLines("\n*** Get  mitochondrial genes dataset 1 ***\n")

mito.genes <- grep(pattern = mitoRegExpressions, x = rownames(x = seurat.object_1@data), value = T)
percent.mito <- Matrix::colSums(seurat.object_1@raw.data[mito.genes, ])/Matrix::colSums(seurat.object_1@raw.data)
seurat.object_1 <- AddMetaData(object = seurat.object_1, metadata = percent.mito, col.name = "percent.mito")

####################################
### Violin plots for data UNfiltered by Seurat
####################################
### Note: files  VlnPlotPdfA (violin plots) and VlnPlotPdfB (details on means/medians and nCells)
### are generated independently because I haven't found a way to make VlnPlot() to include extra data

writeLines("\n*** Violin plots for data UNfiltered by Seurat dataset 1 ***\n")

### Get violin plots data_1
VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object_1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.x.use = 0)
dev.off()

### Get mean and median data_1
nGene_mean<-(mean(seurat.object_1@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object_1@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object_1@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object_1@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object_1@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object_1@meta.data[,"percent.mito"]))
#
### Make a file with mean and median data_1
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VlnPlot.B.pdf", sep="")
#
par()$oma
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
pdf(file=VlnPlotPdfB, width = 7, height = 7) ## keep this after redefining par(oma=...)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
#
FigureTitle<-paste(c("nCells_before_filters=",nCellsInOriginalMatrix_1), sep = "", collapse="")
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                 sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,4),"\n",
                   "median=",round(percent_mito_median,4)),
                 sep = "", collapse="")
#
mtext(text=FigureTitle,  side=1, line=4, outer=F, at = 0.5,  cex = 1, col = "blue")
mtext(text=nGeneStats,   side=3, line=0.6, outer=F, at = 0.15, cex = 1, col = "blue")
mtext(text=nUMIStats,    side=3, line=0.6, outer=F, at = 0.6,  cex = 1, col = "blue")
mtext(text=mitoStats,    side=3, line=0.6, outer=F, at = 1,    cex = 1, col = "blue")
dev.off()
#
### Merge VlnPlot, and mean and median, pdf files data_1
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VlnPlot.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)

####################################
### Gene scatter plots for data UNfiltered by Seurat
####################################
writeLines("\n*** Gene scatter plots for data UNfiltered by Seurat dataset 1 ***\n")

pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_GenePlot.pdf", sep=""),width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object_1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object_1, gene1 = "nUMI", gene2 = "nGene")
dev.off()

####################################
### Filter cells based gene counts and mitochondrial representation
####################################
writeLines("\n*** Filter cells based gene counts and mitochondrial representation dataset 1 ***\n")

seurat.object_1<-FilterCells(object = seurat.object_1, subset.names = c("nGene", "percent.mito"), low.thresholds = LowThresholds, high.thresholds = HighThresholds)
seurat.object_1
nCellsInFilteredMatrix_1<-length(seurat.object_1@meta.data$percent.mito)


####################################
### Violin plots for data filtered by Seurat
####################################
### Note: files  VlnPlotPdfA (violin plots) and VlnPlotPdfB (details on means/medians, nCells and filters)
### are generated independently because I haven't found a way to make VlnPlot() to include extra data

writeLines("\n*** Violin plots for data filtered by Seurat dataset 1 ***\n")

VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VlnPlot.seurat_filtered.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object_1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.x.use = 0, size.title.use = 16)
dev.off()
#
### Get mean and median
nGene_mean<-(mean(seurat.object_1@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object_1@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object_1@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object_1@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object_1@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object_1@meta.data[,"percent.mito"]))
#
### Here making a pdf file that has details on means/medians, nCells and filters
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VlnPlot.seurat_filtered.B.pdf", sep="")
#
par()$oma
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
pdf(file=VlnPlotPdfB, width = 7, height = 7) ## keep this after redefining par(oma=...)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
#
FigureTitle<-paste(c("nCells_after_filters=",nCellsInFilteredMatrix_1, " ",
                     "nGene(min=", LowThresholds[1], ",max=" ,HighThresholds[1], ")
                     %mito(min=", LowThresholds[2], ",max=" ,HighThresholds[2], ")"
), sep = "", collapse="")
FigureTitle<-gsub("\n| +", "  ", FigureTitle, perl=T)
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                 sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,4),"\n",
                   "median=",round(percent_mito_median,4)),
                 sep = "", collapse="")
#
mtext(text=FigureTitle,  side=1, line=4,   outer=F, at = 0.5,  cex = 1, col = "blue")
mtext(text=nGeneStats,   side=3, line=0.6, outer=F, at = 0.15, cex = 1, col = "blue")
mtext(text=nUMIStats,    side=3, line=0.6, outer=F, at = 0.6,  cex = 1, col = "blue")
mtext(text=mitoStats,    side=3, line=0.6, outer=F, at = 1,    cex = 1, col = "blue")
dev.off()
#
### Merge VlnPlot, and mean and median, pdf files
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VlnPlot.seurat_filtered.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)

####################################
### Gene scatter plots for data filtered by Seurat
####################################
writeLines("\n*** Gene scatter plots for data filtered by Seurat dataset 1 ***\n")

pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_GenePlot.seurat_filtered.pdf", sep=""), width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object_1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object_1, gene1 = "nUMI", gene2 = "nGene")
dev.off()

####################################
### Normalize data
####################################
writeLines("\n*** Normalize data dataset 1 ***\n")

seurat.object_1 <- NormalizeData(object = seurat.object_1, normalization.method = "LogNormalize", scale.factor = ScaleFactor)

####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes dataset 1 ***\n")

### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.

seurat.object_1 <- FindVariableGenes(object = seurat.object_1, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff, do.plot=F)
length(x = seurat.object_1@var.genes)
write(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VariableGenes.txt", sep=""), x=seurat.object_1@var.genes)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VariableGenes.pdf", sep=""))
VariableGenePlot(object = seurat.object_1,  x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff)
dev.off()

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation dataset 1 ***\n")

seurat.object_1 <- ScaleData(object = seurat.object_1, vars.to.regress = c("nUMI", "percent.mito"), display.progress=F)

####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA dataset 1 ***\n")

seurat.object_1 <- RunPCA(object = seurat.object_1, pc.genes = seurat.object_1@var.genes, do.print = T, pcs.print = PrintPCA.PcsPrint, genes.print = PrintPCA.GenesPrint)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
VizPCA(object = seurat.object_1, pcs.use = VizPCA.PcsUse)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_PCAPlot.pdf", sep=""))
PCAPlot(object = seurat.object_1, dim.1 = 1, dim.2 = 2, no.legend=T)
dev.off()
#
seurat.object_1 <- ProjectPCA(object = seurat.object_1, do.print = F)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_PCHeatmap.C1.pdf", sep=""))
PCHeatmap(object = seurat.object_1, pc.use = 1, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_PCHeatmap.C1to",PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
PCHeatmap(object = seurat.object_1, pc.use = 1:PCHeatmapComponentsToPlot, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F, use.full = F)
dev.off()

####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components dataset 1 ***\n")

### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time

pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_PCElbowPlot.pdf", sep=""))
PCElbowPlot(object = seurat.object_1)
dev.off()

####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells dataset 1 ***\n")

StartTimeClustering_1<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object_1 <- FindClusters(object = seurat.object_1, reduction.type = "pca", dims.use = PcaDimsUse, resolution = Resolution1, print.output = 0, save.SNN = T)
EndTimeClustering_1<-Sys.time()
#
CellNames<-seurat.object_1@cell.names
ClusterIdent<-seurat.object_1@ident
Headers<-paste("Cell_barcode","seurat_clusters",sep="\t")
clusters_data<-paste(CellNames,ClusterIdent,sep="\t")
#
OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution1,  "_", PrefixInput1, ".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters_1<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution1,  "_", PrefixInput1, ".SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters_1,file = OutfileNumbClusters)

####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene dataset 1 ***\n")

cluster.averages<-AverageExpression(object = seurat.object_1, use.raw = T)
OutfileClusterAverages<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
writeLines("\n*** Run Non-linear dimensional reduction (tSNE) dataset 1 ***\n")

### NOTE: if the datasets is small you may get
### "Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large."
### One can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10

seurat.object_1 <- RunTSNE(object = seurat.object_1, dims.use = PcaDimsUse, do.fast = T)
#
### Note that you can set do.label=T to help label individual clusters
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot.pdf", sep=""))
TSNEPlot(object = seurat.object_1, do.label = T,label.size=10, plot.title='Dataset1')
dev.off()
outfileRDs_1<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".Dataset.RDS", sep="")
saveRDS(seurat.object_1, file = outfileRDs_1)

####################################
### Colour t-SNE by nGene, nUMI, and percent.mito
####################################
writeLines("\n*** Colour t-SNE by nGene, nUMI, and percent.mito dataset 1 ***\n")

CellPropertiesToTsne<-c("nGene", "nUMI", "percent.mito")
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 15, height = 5)
FeaturePlot(object = seurat.object_1, features.plot = CellPropertiesToTsne, cols.use = c("grey", "blue"), reduction.use = "tsne", nCol = 3, pt.size = 1.5)
dev.off()

####################################
### Colour t-SNE by -infile_colour_tsne_extra_1
####################################
writeLines("\n*** Colour t-SNE by -infile_colour_tsne_extra_1 ***\n")

if (regexpr("^NA$", ColourTsne1, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  
  seurat.object_1.meta.data<-seurat.object_1@meta.data
  ExtraCellProperties <- data.frame(read.table(ColourTsne1, header = T, row.names = 1))
  #
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_tsne_extra_1 and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object_1 <- AddMetaData(object = seurat.object_1, metadata = ExtraCellProperties)
  #
  # Generating outfile
  # Note TSNEPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...))
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""))
  for (property in colnames(ExtraCellProperties)) {
    TSNEPlot(object = seurat.object_1, group.by = property, plot.title = property)
  }
  dev.off()
  
}

####################################
### Colour t-SNE plots showing each requested gene
####################################
writeLines("\n*** Colour t-SNE plots showing each requested gene dataset 1 ***\n")

if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
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
  
  pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  FeaturePlot(object = seurat.object_1, features.plot = c(ListOfGenesForTsnes), cols.use = c(rgb(red = 0.9, green = 0.9, blue = 0.9, alpha = 0.1), rgb(red = 0, green = 0, blue = 1, alpha = Opacity)), reduction.use = "tsne")
  dev.off()
  
}

####################################
### Finding differentially expressed genes (cluster biomarkers)
####################################
writeLines("\n*** Finding differentially expressed genes (cluster biomarkers) dataset 1 ***\n")

StartTimeFindAllMarkers_1<-Sys.time()
FindMarkers.Pseudocount_1 <- 1/length(rownames(seurat.object_1@meta.data))
seurat.object_1.markers <- FindAllMarkers(object = seurat.object_1, only.pos = T, min.pct = FindAllMarkers.MinPct, return.thresh = ThreshReturn, thresh.use = FindAllMarkers.ThreshUse, pseudocount.use=FindMarkers.Pseudocount_1)
EndTimeFindAllMarkers_1<-Sys.time()
#
OufileTsvMarkers_1<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_MarkersPerCluster.tsv",sep="")
OufileCsvMarkers_1<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_MarkersPerCluster.csv",sep="")
write.table(data.frame("GENE"=rownames(seurat.object_1.markers),seurat.object_1.markers), file = OufileTsvMarkers_1, row.names = F,sep="\t",quote = F)
write.csv(seurat.object_1.markers, file = OufileCsvMarkers_1)
#
### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne_1<-(seurat.object_1.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters_1<-length(unique(seurat.object_1.markers[["cluster"]]))

####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes dataset 1 ***\n")

pdfWidth<-7

NumberOfPanesForFeaturesPlot<-(NumberOfClusters_1*NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne_1.list<-top_genes_by_cluster_for_tsne_1[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfHeight<-NumberOfClusters_1
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
VlnPlot(object = seurat.object_1, features.plot = c(top_genes_by_cluster_for_tsne_1.list), use.raw = T, y.log = T, adjust.use=1,point.size.use = 0.5,size.title.use=VlnPlotSizeTitle)
dev.off()

####################################
### t-SNE plots showing each cluster top genes
####################################
writeLines("\n*** t-SNE plots showing each cluster top genes dataset 1 ***\n")

pdfWidth<-7

pdfHeight<-NumberOfClusters_1
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
FeaturePlot(object = seurat.object_1, features.plot = c(top_genes_by_cluster_for_tsne_1.list), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()

####################################
### Heatmaps
####################################
writeLines("\n*** Heatmaps dataset 1 ***\n")

top_genes_by_cluster_for_heatmap <- seurat.object_1.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotHeatmap, avg_logFC)
### setting slim.col.label to T will print just the cluster IDS instead of every cell name
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_Heatmap.pdf", sep=""))
DoHeatmap(object = seurat.object_1, genes.use = top_genes_by_cluster_for_heatmap$gene, slim.col.label = T, remove.key = T, title = PrefixInput1, cex.row = 6)
dev.off()

####################################
####################################
####################################
### Process dataset 2
####################################
####################################
####################################
writeLines("\n*** Process dataset 2 ***\n")
#
writeLines("\n*** Get  mitochondrial genes ***\n")
#
mito.genes <- grep(pattern = mitoRegExpressions, x = rownames(x = seurat.object_2@data), value = T)
percent.mito <- Matrix::colSums(seurat.object_2@raw.data[mito.genes, ])/Matrix::colSums(seurat.object_2@raw.data)
seurat.object_2 <- AddMetaData(object = seurat.object_2, metadata = percent.mito, col.name = "percent.mito")
#
####################################
### Violin plots for data UNfiltered by Seurat
####################################
### Note: files  VlnPlotPdfA (violin plots) and VlnPlotPdfB (details on means/medians and nCells)
### are generated independently because I haven't found a way to make VlnPlot() to include extra data
#
writeLines("\n*** Violin plots for data UNfiltered by Seurat dataset 2 ***\n")
#
### Get violin plots data_2
VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VlnPlot.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object_2, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.x.use = 0)
dev.off()
#
### Get mean and median data_2
nGene_mean<-(mean(seurat.object_2@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object_2@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object_2@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object_2@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object_2@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object_2@meta.data[,"percent.mito"]))
#
### Make a file with mean and median data_2
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VlnPlot.B.pdf", sep="")
#
par()$oma
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
pdf(file=VlnPlotPdfB, width = 7, height = 7) ## keep this after redefining par(oma=...)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
#
FigureTitle<-paste(c("nCells_before_filters=",nCellsInOriginalMatrix_2), sep = "", collapse="")
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                 sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,4),"\n",
                   "median=",round(percent_mito_median,4)),
                 sep = "", collapse="")
#
mtext(text=FigureTitle,  side=1, line=4, outer=F, at = 0.5,  cex = 1, col = "blue")
mtext(text=nGeneStats,   side=3, line=0.6, outer=F, at = 0.15, cex = 1, col = "blue")
mtext(text=nUMIStats,    side=3, line=0.6, outer=F, at = 0.6,  cex = 1, col = "blue")
mtext(text=mitoStats,    side=3, line=0.6, outer=F, at = 1,    cex = 1, col = "blue")
dev.off()
#
### Merge VlnPlot, and mean and median, pdf files data_2
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VlnPlot.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)
#
####################################
### Gene scatter plots for data UNfiltered by Seurat
####################################
writeLines("\n*** Gene scatter plots for data UNfiltered by Seurat dataset 2 ***\n")
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_GenePlot.pdf", sep=""),width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object_2, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object_2, gene1 = "nUMI", gene2 = "nGene")
dev.off()
#
####################################
### Filter cells based gene counts and mitochondrial representation
####################################
writeLines("\n*** Filter cells based gene counts and mitochondrial representation dataset 2 ***\n")
#
seurat.object_2<-FilterCells(object = seurat.object_2, subset.names = c("nGene", "percent.mito"), low.thresholds = LowThresholds, high.thresholds = HighThresholds)
seurat.object_2
nCellsInFilteredMatrix_2<-length(seurat.object_2@meta.data$percent.mito)
#
####################################
### Violin plots for data filtered by Seurat
####################################
### Note: files  VlnPlotPdfA (violin plots) and VlnPlotPdfB (details on means/medians, nCells and filters)
### are generated independently because I haven't found a way to make VlnPlot() to include extra data
#
writeLines("\n*** Violin plots for data filtered by Seurat dataset 2 ***\n")
#
VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VlnPlot.seurat_filtered.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object_2, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.x.use = 0, size.title.use = 16)
dev.off()
#
### Get mean and median
nGene_mean<-(mean(seurat.object_2@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object_2@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object_2@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object_2@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object_2@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object_2@meta.data[,"percent.mito"]))
#
### Here making a pdf file that has details on means/medians, nCells and filters
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VlnPlot.seurat_filtered.B.pdf", sep="")
#
par()$oma
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
pdf(file=VlnPlotPdfB, width = 7, height = 7) ## keep this after redefining par(oma=...)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
#
FigureTitle<-paste(c("nCells_after_filters=",nCellsInFilteredMatrix_2, " ",
                     "nGene(min=", LowThresholds[1], ",max=" ,HighThresholds[1], ")
                     %mito(min=", LowThresholds[2], ",max=" ,HighThresholds[2], ")"
), sep = "", collapse="")
FigureTitle<-gsub("\n| +", "  ", FigureTitle, perl=T)
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                 sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,4),"\n",
                   "median=",round(percent_mito_median,4)),
                 sep = "", collapse="")
#
mtext(text=FigureTitle,  side=1, line=4, outer=F,   at = 0.5,  cex = 1, col = "blue")
mtext(text=nGeneStats,   side=3, line=0.6, outer=F, at = 0.15, cex = 1, col = "blue")
mtext(text=nUMIStats,    side=3, line=0.6, outer=F, at = 0.6,  cex = 1, col = "blue")
mtext(text=mitoStats,    side=3, line=0.6, outer=F, at = 1,    cex = 1, col = "blue")
dev.off()
#
### Merge VlnPlot, and mean and median, pdf files
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VlnPlot.seurat_filtered.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)
#
####################################
### Gene scatter plots for data filtered by Seurat
####################################
writeLines("\n*** Gene scatter plots for data filtered by Seurat dataset 2 ***\n")
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_GenePlot.seurat_filtered.pdf", sep=""), width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object_2, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object_2, gene1 = "nUMI", gene2 = "nGene")
dev.off()
#
####################################
### Normalize data
####################################
writeLines("\n*** Normalize data dataset 2 ***\n")
#
seurat.object_2 <- NormalizeData(object = seurat.object_2, normalization.method = "LogNormalize", scale.factor = ScaleFactor)
#
####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes dataset 2 ***\n")
#
### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
#
seurat.object_2 <- FindVariableGenes(object = seurat.object_2, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff, do.plot=F)
length(x = seurat.object_2@var.genes)
write(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VariableGenes.txt", sep=""), x=seurat.object_2@var.genes)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VariableGenes.pdf", sep=""))
VariableGenePlot(object = seurat.object_2,  x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff)
dev.off()
#
####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation dataset 2 ***\n")
#
seurat.object_2 <- ScaleData(object = seurat.object_2, vars.to.regress = c("nUMI", "percent.mito"), display.progress=F)
#
####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA dataset 2 ***\n")
#
seurat.object_2 <- RunPCA(object = seurat.object_2, pc.genes = seurat.object_2@var.genes, do.print = T, pcs.print = PrintPCA.PcsPrint, genes.print = PrintPCA.GenesPrint)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
VizPCA(object = seurat.object_2, pcs.use = VizPCA.PcsUse)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_PCAPlot.pdf", sep=""))
PCAPlot(object = seurat.object_2, dim.1 = 1, dim.2 = 2, no.legend=T)
dev.off()
#
seurat.object_2 <- ProjectPCA(object = seurat.object_2, do.print = F)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_PCHeatmap.C1.pdf", sep=""))
PCHeatmap(object = seurat.object_2, pc.use = 1, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_PCHeatmap.C1to",PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
PCHeatmap(object = seurat.object_2, pc.use = 1:PCHeatmapComponentsToPlot, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F, use.full = F)
dev.off()
#
####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components dataset 2 ***\n")
#
### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_PCElbowPlot.pdf", sep=""))
PCElbowPlot(object = seurat.object_2)
dev.off()
#
####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells dataset 2 ***\n")
#
StartTimeClustering_2<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object_2 <- FindClusters(object = seurat.object_2, reduction.type = "pca", dims.use = PcaDimsUse, resolution = Resolution2, print.output = 0, save.SNN = T)
EndTimeClustering_2<-Sys.time()
#
CellNames<-seurat.object_2@cell.names
ClusterIdent<-seurat.object_2@ident
Headers<-paste("CLUSTERS","seurat_clusters",sep="\t")
clusters_data<-paste(CellNames,ClusterIdent,sep="\t")
#
OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution2,  "_", PrefixInput2, ".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters_2<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution2,  "_", PrefixInput2, ".SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters_2,file = OutfileNumbClusters)
#
####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene dataset 2 ***\n")
#
cluster.averages<-AverageExpression(object = seurat.object_2, use.raw = T)
OutfileClusterAverages<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
writeLines("\n*** Run Non-linear dimensional reduction (tSNE) dataset 2 ***\n")
#
### NOTE: if the datasets is small you may get
### "Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large."
### One can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
#
seurat.object_2 <- RunTSNE(object = seurat.object_2, dims.use = PcaDimsUse, do.fast = T)
#
### Note that you can set do.label=T to help label individual clusters
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot.pdf", sep=""))
TSNEPlot(object = seurat.object_2, do.label = T,label.size=10, plot.title='Dataset2')
dev.off()
outfileRDs_2<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".Dataset.RDS", sep="")
saveRDS(seurat.object_2, file = outfileRDs_2)
#
#
####################################
### Colour t-SNE by nGene, nUMI, and percent.mito
####################################
writeLines("\n*** Colour t-SNE by nGene, nUMI, and percent.mito dataset 2 ***\n")
#
CellPropertiesToTsne<-c("nGene", "nUMI", "percent.mito")
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 15, height = 5)
FeaturePlot(object = seurat.object_2, features.plot = CellPropertiesToTsne, cols.use = c("grey", "blue"), reduction.use = "tsne", nCol = 3, pt.size = 1.5)
dev.off()
#
####################################
### Colour t-SNE by -infile_colour_tsne_extra_2
####################################
writeLines("\n*** Colour t-SNE by -infile_colour_tsne_extra_2 ***\n")
#
if (regexpr("^NA$", ColourTsne2, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  seurat.object_2.meta.data<-seurat.object_2@meta.data
  ExtraCellProperties <- data.frame(read.table(ColourTsne2, header = T, row.names = 1))
  #
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_tsne_extra_2 and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object_2 <- AddMetaData(object = seurat.object_2, metadata = ExtraCellProperties)
  #
  # Generating outfile
  # Note TSNEPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...))
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""))
  for (property in colnames(ExtraCellProperties)) {
    TSNEPlot(object = seurat.object_2, group.by = property, plot.title = property)
  }
  dev.off()
}

####################################
### Colour t-SNE plots showing each requested gene
####################################
writeLines("\n*** Colour t-SNE plots showing each requested gene dataset 2 ***\n")

if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
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
  pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  FeaturePlot(object = seurat.object_2, features.plot = c(ListOfGenesForTsnes), cols.use = c(rgb(red = 0.9, green = 0.9, blue = 0.9, alpha = 0.1), rgb(red = 0, green = 0, blue = 1, alpha = Opacity)), reduction.use = "tsne")
  dev.off()
}

####################################
### Finding differentially expressed genes (cluster biomarkers)
####################################
writeLines("\n*** Finding differentially expressed genes (cluster biomarkers) dataset 2 ***\n")
#
StartTimeFindAllMarkers_2<-Sys.time()
FindMarkers.Pseudocount_2 <- 1/length(rownames(seurat.object_2@meta.data))
seurat.object_2.markers <- FindAllMarkers(object = seurat.object_2, only.pos = T, min.pct = FindAllMarkers.MinPct, return.thresh = ThreshReturn, thresh.use = FindAllMarkers.ThreshUse, pseudocount.use=FindMarkers.Pseudocount_2)
EndTimeFindAllMarkers_2<-Sys.time()
#
OufileTsvMarkers_2<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_MarkersPerCluster.tsv",sep="")
OufileCsvMarkers_2<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_MarkersPerCluster.csv",sep="")
write.table(data.frame("GENE"=rownames(seurat.object_2.markers),seurat.object_2.markers), file = OufileTsvMarkers_2, row.names = F,sep="\t",quote = F)
write.csv(seurat.object_2.markers, file = OufileCsvMarkers_2)
#
### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne_2<-(seurat.object_2.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters_2<-length(unique(seurat.object_2.markers[["cluster"]]))
#
####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes dataset 2 ***\n")
#
pdfWidth<-7
#
NumberOfPanesForFeaturesPlot<-(NumberOfClusters_2*NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne_2.list<-top_genes_by_cluster_for_tsne_2[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfHeight<-NumberOfClusters_2
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
VlnPlot(object = seurat.object_2, features.plot = c(top_genes_by_cluster_for_tsne_2.list), use.raw = T, y.log = T, adjust.use=1,point.size.use = 0.5,size.title.use=VlnPlotSizeTitle)
dev.off()
#
####################################
### t-SNE plots showing each cluster top genes
####################################
writeLines("\n*** t-SNE plots showing each cluster top genes dataset 2 ***\n")
#
pdfWidth<-7
#
pdfHeight<-NumberOfClusters_2
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
FeaturePlot(object = seurat.object_2, features.plot = c(top_genes_by_cluster_for_tsne_2.list), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
#
####################################
### Heatmaps
####################################
writeLines("\n*** Heatmaps dataset 2 ***\n")
#
top_genes_by_cluster_for_heatmap <- seurat.object_2.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotHeatmap, avg_logFC)
### setting slim.col.label to T will print just the cluster IDS instead of every cell name
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_Heatmap.pdf", sep=""))
DoHeatmap(object = seurat.object_2, genes.use = top_genes_by_cluster_for_heatmap$gene, slim.col.label = T, remove.key = T, title = PrefixInput2, cex.row = 6)
dev.off()






####################################
####################################
####################################
### Process merged dataset
####################################
####################################
####################################
writeLines("\n*** Process merged dataset ***\n")
#
writeLines("\n*** Get  mitochondrial genes merged dataset ***\n")
#
mito.genes <- grep(pattern = mitoRegExpressions, x = rownames(x = seurat.object_c12@data), value = T)
percent.mito <- Matrix::colSums(seurat.object_c12@raw.data[mito.genes, ])/Matrix::colSums(seurat.object_c12@raw.data)
seurat.object_c12 <- AddMetaData(object = seurat.object_c12, metadata = percent.mito, col.name = "percent.mito")
#
####################################
### Violin plots for data UNfiltered by Seurat
####################################
### Note: files  VlnPlotPdfA (violin plots) and VlnPlotPdfB (details on means/medians and nCells)
### are generated independently because I haven't found a way to make VlnPlot() to include extra data
#
writeLines("\n*** Violin plots for data UNfiltered by Seurat merged dataset ***\n")
#
OufileTsvMarkers_2<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_MarkersPerCluster.tsv",sep="")
#
#
### Get violin plots data_c12
VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VlnPlot.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object_c12, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.x.use = 0)
dev.off()
#
### Get mean and median data_c12
nGene_mean<-(mean(seurat.object_c12@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object_c12@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object_c12@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object_c12@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object_c12@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object_c12@meta.data[,"percent.mito"]))
#
### Make a file with mean and median data_c12
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VlnPlot.B.pdf", sep="")
#
par()$oma
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
pdf(file=VlnPlotPdfB, width = 7, height = 7) ## keep this after redefining par(oma=...)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
#
FigureTitle<-paste(c("nCells_before_filters=",nCellsInOriginalMatrix_c12), sep = "", collapse="")
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                 sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,4),"\n",
                   "median=",round(percent_mito_median,4)),
                 sep = "", collapse="")
#
mtext(text=FigureTitle,  side=1, line=4, outer=F, at = 0.5,  cex = 1, col = "blue")
mtext(text=nGeneStats,   side=3, line=0.6, outer=F, at = 0.15, cex = 1, col = "blue")
mtext(text=nUMIStats,    side=3, line=0.6, outer=F, at = 0.6,  cex = 1, col = "blue")
mtext(text=mitoStats,    side=3, line=0.6, outer=F, at = 1,    cex = 1, col = "blue")
dev.off()
#
### Merge VlnPlot, and mean and median, pdf files data_c12
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VlnPlot.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)
#
####################################
### Gene scatter plots for data UNfiltered by Seurat
####################################
writeLines("\n*** Gene scatter plots for data UNfiltered by Seurat merged dataset ***\n")
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_GenePlot.pdf", sep=""),width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object_c12, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object_c12, gene1 = "nUMI", gene2 = "nGene")
dev.off()
#
####################################
### Filter cells based gene counts and mitochondrial representation
####################################
writeLines("\n*** Filter cells based gene counts and mitochondrial representation merged dataset ***\n")
#
seurat.object_c12<-FilterCells(object = seurat.object_c12, subset.names = c("nGene", "percent.mito"), low.thresholds = LowThresholds, high.thresholds = HighThresholds)
seurat.object_c12
nCellsInFilteredMatrix_c12<-length(seurat.object_c12@meta.data$percent.mito)
#
####################################
### Violin plots for data filtered by Seurat
####################################
### Note: files  VlnPlotPdfA (violin plots) and VlnPlotPdfB (details on means/medians, nCells and filters)
### are generated independently because I haven't found a way to make VlnPlot() to include extra data
#
writeLines("\n*** Violin plots for data filtered by Seurat merged dataset ***\n")
#
VlnPlotPdfA<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VlnPlot.seurat_filtered.A.pdf", sep="")
pdf(file=VlnPlotPdfA, width = 7, height = 7)
VlnPlot(object = seurat.object_c12, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, size.x.use = 0, size.title.use = 16)
dev.off()
#
### Get mean and median
nGene_mean<-(mean(seurat.object_c12@meta.data[,"nGene"]))
nGene_median<-(median(seurat.object_c12@meta.data[,"nGene"]))
nUMI_mean<-(mean(seurat.object_c12@meta.data[,"nUMI"]))
nUMI_median<-(median(seurat.object_c12@meta.data[,"nUMI"]))
percent_mito_mean<-(mean(seurat.object_c12@meta.data[,"percent.mito"]))
percent_mito_median<-(median(seurat.object_c12@meta.data[,"percent.mito"]))
#
### Here making a pdf file that has details on means/medians, nCells and filters
VlnPlotPdfB<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VlnPlot.seurat_filtered.B.pdf", sep="")
#
par()$oma
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
pdf(file=VlnPlotPdfB, width = 7, height = 7) ## keep this after redefining par(oma=...)
plot(x=NA,y=NA,xlim=c(0,1),ylim=c(0,1),axes=F,bty="n",xlab="",ylab="")
#
FigureTitle<-paste(c("nCells_after_filters=",nCellsInFilteredMatrix_c12, " ",
                     "nGene(min=", LowThresholds[1], ",max=" ,HighThresholds[1], ")
                     %mito(min=", LowThresholds[2], ",max=" ,HighThresholds[2], ")"
), sep = "", collapse="")
FigureTitle<-gsub("\n| +", "  ", FigureTitle, perl=T)
nGeneStats<-paste(c("mean=",round(nGene_mean,0),"\n",
                    "median=",round(nGene_median,0)),
                  sep = "", collapse="")
nUMIStats<-paste(c("mean=",round(nUMI_mean,0),"\n",
                   "median=",round(nUMI_median,0)),
                 sep = "", collapse="")
mitoStats<-paste(c("mean=",round(percent_mito_mean,4),"\n",
                   "median=",round(percent_mito_median,4)),
                 sep = "", collapse="")
#
mtext(text=FigureTitle,  side=1, line=4, outer=F,   at = 0.5,  cex = 1, col = "blue")
mtext(text=nGeneStats,   side=3, line=0.6, outer=F, at = 0.15, cex = 1, col = "blue")
mtext(text=nUMIStats,    side=3, line=0.6, outer=F, at = 0.6,  cex = 1, col = "blue")
mtext(text=mitoStats,    side=3, line=0.6, outer=F, at = 1,    cex = 1, col = "blue")
dev.off()
#
### Merge VlnPlot, and mean and median, pdf files
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VlnPlot.seurat_filtered.pdf", sep="")
CommandToOverlapPdfs<-paste(c("pdftk ", VlnPlotPdfA, " background ", VlnPlotPdfB, " output ", VlnPlotPdf), sep = "", collapse = "")
system(CommandToOverlapPdfs, wait = T)
file.remove(VlnPlotPdfA,VlnPlotPdfB)
#
####################################
### Gene scatter plots for data filtered by Seurat
####################################
writeLines("\n*** Gene scatter plots for data filtered by Seurat merged dataset ***\n")
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_GenePlot.seurat_filtered.pdf", sep=""), width=14, height=7)
par(mfrow = c(1, 2))
GenePlot(object = seurat.object_c12, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seurat.object_c12, gene1 = "nUMI", gene2 = "nGene")
dev.off()
#
####################################
### Normalize data
####################################
writeLines("\n*** Normalize data merged dataset ***\n")
#
seurat.object_c12 <- NormalizeData(object = seurat.object_c12, normalization.method = "LogNormalize", scale.factor = ScaleFactor)
#
####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes merged dataset ***\n")
#
### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
#
seurat.object_c12 <- FindVariableGenes(object = seurat.object_c12, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff, do.plot=F)
length(x = seurat.object_c12@var.genes)
write(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VariableGenes.txt", sep=""), x=seurat.object_c12@var.genes)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VariableGenes.pdf", sep=""))
VariableGenePlot(object = seurat.object_c12,  x.low.cutoff = XLowCutoff, x.high.cutoff = XHighCutoff, y.cutoff = YCutoff)
dev.off()
#
####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation merged dataset ***\n")
#
seurat.object_c12 <- ScaleData(object = seurat.object_c12, vars.to.regress = c("nUMI", "percent.mito"), display.progress=F)
#
####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA merged dataset ***\n")
#
seurat.object_c12 <- RunPCA(object = seurat.object_c12, pc.genes = seurat.object_c12@var.genes, do.print = T, pcs.print = PrintPCA.PcsPrint, genes.print = PrintPCA.GenesPrint)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
VizPCA(object = seurat.object_c12, pcs.use = VizPCA.PcsUse)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_PCAPlot.pdf", sep=""))
PCAPlot(object = seurat.object_c12, dim.1 = 1, dim.2 = 2, no.legend=T)
dev.off()
#
seurat.object_c12 <- ProjectPCA(object = seurat.object_c12, do.print = F)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_PCHeatmap.C1.pdf", sep=""))
PCHeatmap(object = seurat.object_c12, pc.use = 1, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_PCHeatmap.C1to",PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
PCHeatmap(object = seurat.object_c12, pc.use = 1:PCHeatmapComponentsToPlot, cells.use = PCHeatmapCellsUse, do.balanced = T, label.columns = F, use.full = F)
dev.off()
#
####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components merged dataset ***\n")
#
### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_PCElbowPlot.pdf", sep=""))
PCElbowPlot(object = seurat.object_c12)
dev.off()
#
####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells merged dataset ***\n")
#
StartTimeClustering_c12<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object_c12 <- FindClusters(object = seurat.object_c12, reduction.type = "pca", dims.use = PcaDimsUse, resolution = Resolution12, print.output = 0, save.SNN = T)
EndTimeClustering_c12<-Sys.time()
#
CellNames<-seurat.object_c12@cell.names
ClusterIdent<-seurat.object_c12@ident
Headers<-paste("CLUSTERS","seurat_clusters",sep="\t")
clusters_data<-paste(CellNames,ClusterIdent,sep="\t")
#
OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution12, "_merged.SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters_c12<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution12, "_merged.SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters_c12,file = OutfileNumbClusters)
#
####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene merged dataset ***\n")
#
cluster.averages<-AverageExpression(object = seurat.object_c12, use.raw = T)
OutfileClusterAverages<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
writeLines("\n*** Run Non-linear dimensional reduction (tSNE) merged dataset ***\n")
#
### NOTE: if the datasets is small you may get
### "Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large."
### One can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
#
seurat.object_c12 <- RunTSNE(object = seurat.object_c12, dims.use = PcaDimsUse, do.fast = T)
#
### Note that you can set do.label=T to help label individual clusters
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_TSNEPlot.pdf", sep=""))
TSNEPlot(object = seurat.object_c12, do.label = T,label.size=10, plot.title=paste("Merged ", PrefixInput1, " and ", PrefixInput2, sep = "", collapse = ""))
dev.off()
outfileRDs_c12<-paste(Tempdir,"/",PrefixOutfiles,"_merged.Dataset.RDS", sep="")
saveRDS(seurat.object_c12, file = outfileRDs_c12)
#

####################################
### Colour t-SNE by nGene, nUMI, and percent.mito
####################################
writeLines("\n*** Colour t-SNE by nGene, nUMI, and percent.mito each dataset ***\n")
#
CellPropertiesToTsne<-c("nGene", "nUMI", "percent.mito")
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_TSNEPlot_QC.pdf", sep=""), width = 15, height = 5)
FeaturePlot(object = seurat.object_c12, features.plot = CellPropertiesToTsne, cols.use = c("grey", "blue"), reduction.use = "tsne", nCol = 3, pt.size = 1.5)
dev.off()

####################################
### Colour t-SNE by -infile_colour_tsne_extra_merged
####################################
writeLines("\n*** Colour t-SNE by -infile_colour_tsne_extra_merged ***\n")
#
if (regexpr("^NA$", ColourTsneM, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  #
  seurat.object_c12.meta.data<-seurat.object_c12@meta.data
  ExtraCellProperties <- data.frame(read.table(ColourTsneM, header = T, row.names = 1))
  #
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_tsne_extra_merged and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object_c12 <- AddMetaData(object = seurat.object_c12, metadata = ExtraCellProperties)
  #
  # Generating outfile
  # Note TSNEPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...))
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""))
  for (property in colnames(ExtraCellProperties)) {
    TSNEPlot(object = seurat.object_c12, group.by = property, plot.title = property)
  }
  dev.off()
  #
}
#
####################################
### Colour t-SNE plots showing each requested gene
####################################
writeLines("\n*** Colour t-SNE plots showing each requested gene merged dataset ***\n")
#
if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
  print("No selected genes for t-SNE plots")
}else{
  #
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
  #
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  FeaturePlot(object = seurat.object_c12, features.plot = c(ListOfGenesForTsnes), cols.use = c(rgb(red = 0.9, green = 0.9, blue = 0.9, alpha = 0.1), rgb(red = 0, green = 0, blue = 1, alpha = Opacity)), reduction.use = "tsne")
  dev.off()
  #
}
#
####################################
### Finding differentially expressed genes (cluster biomarkers)
####################################
writeLines("\n*** Finding differentially expressed genes (cluster biomarkers) merged dataset ***\n")
#
StartTimeFindAllMarkers_c12<-Sys.time()
FindMarkers.Pseudocount_c12 <- 1/length(rownames(seurat.object_c12@meta.data))
seurat.object_c12.markers <- FindAllMarkers(object = seurat.object_c12, only.pos = T, min.pct = FindAllMarkers.MinPct, return.thresh = ThreshReturn, thresh.use = FindAllMarkers.ThreshUse, pseudocount.use=FindMarkers.Pseudocount_c12)
EndTimeFindAllMarkers_c12<-Sys.time()
#
OufileTsvMarkers_c12<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_MarkersPerCluster.tsv",sep="")
OufileCsvMarkers_c12<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_MarkersPerCluster.csv",sep="")
write.table(data.frame("GENE"=rownames(seurat.object_c12.markers),seurat.object_c12.markers), file = OufileTsvMarkers_c12, row.names = F,sep="\t",quote = F)
write.csv(seurat.object_c12.markers, file = OufileCsvMarkers_c12)
#
### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne_c12<-(seurat.object_c12.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters_c12<-length(unique(seurat.object_c12.markers[["cluster"]]))
#
####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes merged dataset ***\n")
#
pdfWidth<-7
#
NumberOfPanesForFeaturesPlot<-(NumberOfClusters_c12*NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne_c12.list<-top_genes_by_cluster_for_tsne_c12[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfHeight<-NumberOfClusters_c12
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
VlnPlot(object = seurat.object_c12, features.plot = c(top_genes_by_cluster_for_tsne_c12.list), use.raw = T, y.log = T, adjust.use=1,point.size.use = 0.5,size.title.use=VlnPlotSizeTitle)
dev.off()
#
####################################
### t-SNE plots showing each cluster top genes
####################################
writeLines("\n*** t-SNE plots showing each cluster top genes merged dataset ***\n")
#
pdfWidth<-7
#
pdfHeight<-NumberOfClusters_c12
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
FeaturePlot(object = seurat.object_c12, features.plot = c(top_genes_by_cluster_for_tsne_c12.list), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
#
####################################
### Heatmaps
####################################
writeLines("\n*** Heatmaps merged dataset ***\n")
#
top_genes_by_cluster_for_heatmap <- seurat.object_c12.markers %>% group_by(cluster) %>% top_n(NumberOfGenesPerClusterToPlotHeatmap, avg_logFC)
### setting slim.col.label to T will print just the cluster IDS instead of every cell name
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_Heatmap.pdf", sep=""))
DoHeatmap(object = seurat.object_c12, genes.use = top_genes_by_cluster_for_heatmap$gene, slim.col.label = T, remove.key = T, title = paste("Merged ", PrefixInput1, " and ", PrefixInput2, sep = "", collapse = ""), cex.row = 6)
dev.off()
#
####################################
### t-SNE plots coloured by dataset
####################################
writeLines("\n*** t-SNE plots coloured by dataset ***\n")
#
### Save a table with barcode_ids and datasets
OutfileCellDatasets<-paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_cell_ids.tsv", sep="")
ColNames_1<-colnames(seurat.object_c12@raw.data)[grep(pattern = "Dataset1", x = colnames(seurat.object_c12@raw.data))]
ColNames_2<-colnames(seurat.object_c12@raw.data)[grep(pattern = "Dataset2", x = colnames(seurat.object_c12@raw.data))]
Headers<-paste("CELL_BARCODE", paste("Dataset", sep = "", collapse = "") ,sep="\t")
write.table(Headers,file = OutfileCellDatasets, row.names = F, col.names = F, sep="\t", quote = F)
write.table(x=paste(ColNames_1, PrefixInput1, sep = "\t"), file=OutfileCellDatasets, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(x=paste(ColNames_2, PrefixInput2, sep = "\t"), file=OutfileCellDatasets, row.names = F, col.names = F, sep="\t", quote = F, append = T)
DatasetPerBarcode <- data.frame(read.table(OutfileCellDatasets, header = T, row.names = 1))
seurat.object_c12 <- AddMetaData(object = seurat.object_c12, metadata = DatasetPerBarcode)
#
#################
## t-SNE showing both datasets overlapped, coloured per dataset
pdf(file=paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_TSNEPlot_ColourByDataset.pdf", sep=""))
TSNEPlot(object = seurat.object_c12, group.by = "Dataset",
         plot.title = "Colour by dataset", colors.use = c(ColorDataset1Tsne,ColorDataset2Tsne), pt.size = 0.5
)
dev.off()
#
#################
## t-SNE's showing each dataset, coordinates and colours based on clustering of seurat.object_c12
seurat.object_c12_1<-SubsetData(object = seurat.object_c12, cells.use = ColNames_1)
seurat.object_c12_2<-SubsetData(object = seurat.object_c12, cells.use = ColNames_2)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_CoordsAndColoursFromMerged.pdf", sep=""))
TSNEPlot(object = seurat.object_c12_1, plot.title = paste(PrefixInput1, " - coordinates and colours from merged", sep = "", collapse = "") , do.label = T, label.size = 10)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_CoordsAndColoursFromMerged.pdf", sep=""))
TSNEPlot(object = seurat.object_c12_2, plot.title = paste(PrefixInput2, " - coordinates and colours from merged", sep = "", collapse = "") , do.label = T, label.size = 10)
dev.off()
#
#################
## t-SNE's showing each dataset, coordinates based on clustering of seurat.object_c12, colours based on each dataset clustering
InfileClusters1<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution1,  "_", PrefixInput1, ".SEURAT_CellClusters.tsv", sep="")
InfileClusters2<-paste(Tempdir,"/",PrefixOutfiles, "_res", Resolution2,  "_", PrefixInput2, ".SEURAT_CellClusters.tsv", sep="")
Clusters1 <- data.frame(read.table(InfileClusters1, header = T, row.names = 1))
Clusters2 <- data.frame(read.table(InfileClusters2, header = T, row.names = 1))
rownames(Clusters1)<-paste("Dataset1",rownames(Clusters1),sep = "-")
rownames(Clusters2)<-paste("Dataset2",rownames(Clusters2),sep = "-")
seurat.object_c12_1 <- AddMetaData(object = seurat.object_c12_1, metadata = Clusters1)
seurat.object_c12_2 <- AddMetaData(object = seurat.object_c12_2, metadata = Clusters2)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_CoordsFromMerged_ColoursFromItself.pdf", sep=""))
TSNEPlot(object = seurat.object_c12_1, group.by = "seurat_clusters", plot.title = paste(PrefixInput1, " - coordinates from merged, colours from ", PrefixInput1, sep = "", collapse = "") , do.label = T, label.size = 10)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_CoordsFromMerged_ColoursFromItself.pdf", sep=""))
TSNEPlot(object = seurat.object_c12_2, group.by = "seurat_clusters", plot.title = paste(PrefixInput2, " - coordinates from merged, colours from ", PrefixInput2, sep = "", collapse = "") , do.label = T, label.size = 10)
dev.off()
#
#################
## t-SNE's showing each dataset, coordinates based on clustering of seurat.object_c12, colours based on QC nGene, nUMI, percent.mito for each dataset
CellPropertiesToTsne<-c("nGene", "nUMI", "percent.mito")
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_QC_CoordsFromMerge.pdf", sep=""), width = 15, height = 5)
FeaturePlot(object = seurat.object_c12_1, features.plot = CellPropertiesToTsne, cols.use = c("grey", "blue"), reduction.use = "tsne", nCol = 3, pt.size = 1.5)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_QC_CoordsFromMerge.pdf", sep=""), width = 15, height = 5)
FeaturePlot(object = seurat.object_c12_2, features.plot = CellPropertiesToTsne, cols.use = c("grey", "blue"), reduction.use = "tsne", nCol = 3, pt.size = 1.5)
dev.off()





####################################
####################################
####################################
##### Generate Circos plot
####################################
####################################
####################################
PrefixCombOutfiles<-paste(Tempdir,"/",PrefixOutfiles,"_CombDatasets", sep="")
#
marker_file_list <- c(Dataset1 = OufileCsvMarkers_1, Dataset2 = OufileCsvMarkers_2)
fList <- c(Dataset1 = outfileRDs_1, Dataset2 = outfileRDs_2, CombDatasets = outfileRDs_c12)
objList <- lapply(fList, readRDS)
objList
#
single_obj_list <- c(Dataset1 = objList$Dataset1, Dataset2 = objList$Dataset2)
single_obj_list
#
### edge_cutoff = 0.1 ## default
res <- cluster_map(marker_file_list, edge_cutoff = EdgeCutoff, output = PrefixCombOutfiles, single_obj_list = single_obj_list, comb_obj = objList$CombDatasets)
res ### This command is needed to make the circos plot
#
#################
##### Generate log2 and log10 scatter plot of number of genes
#################
### adds correl
#
rowMeansNoLogDataset_1<-rowMeans(input.matrix_1)
rowMeansNoLogDatasetWo0_1<-rowMeansNoLogDataset_1
#
rowMeansNoLogDataset_2<-rowMeans(input.matrix_2)
rowMeansNoLogDatasetWo0_2<-rowMeansNoLogDataset_2
#
ValueToReplace0sInLog<-min(c(rowMeansNoLogDataset_1[which(rowMeansNoLogDataset_1 > 0)], rowMeansNoLogDataset_2[which(rowMeansNoLogDataset_2 > 0)]))
#
rowMeansNoLogDatasetWo0_2[rowMeansNoLogDatasetWo0_2==0]<-ValueToReplace0sInLog
rowMeansLog10DatasetWo0_2<-log10(rowMeansNoLogDatasetWo0_2)
rowMeansNoLogDatasetWo0_1[rowMeansNoLogDatasetWo0_1==0]<-ValueToReplace0sInLog
rowMeansLog10DatasetWo0_1<-log10(rowMeansNoLogDatasetWo0_1)
#
r.cor.nolog <- round(cor(x=rowMeansNoLogDataset_1, y=rowMeansNoLogDataset_2,  method = "pearson",  use = "pairwise.complete.obs"), digits = 3)
s.cor.nolog <- round(cor(x=rowMeansNoLogDataset_1, y=rowMeansNoLogDataset_2,  method = "spearman", use = "pairwise.complete.obs"), digits = 3)
r.cor.log10 <- round(cor(x=rowMeansLog10DatasetWo0_1, y=rowMeansLog10DatasetWo0_2,  method = "pearson",  use = "pairwise.complete.obs"), digits = 3)
s.cor.log10 <- round(cor(x=rowMeansLog10DatasetWo0_1, y=rowMeansLog10DatasetWo0_2,  method = "spearman", use = "pairwise.complete.obs"), digits = 3)
#
##### log10 scatter plot of number of genes
pdf(file=paste(PrefixCombOutfiles, ".reads_per_gene_scatter.log10.pdf", sep=""))
plot(x=rowMeansLog10DatasetWo0_1, y=rowMeansLog10DatasetWo0_2,
     pch=16,col=rgb(0.35,0.70,0.90,0.25),
     xlim=c(min(rowMeansLog10DatasetWo0_1,rowMeansLog10DatasetWo0_2), max(rowMeansLog10DatasetWo0_1,rowMeansLog10DatasetWo0_2)), ylim=c(min(rowMeansLog10DatasetWo0_1,rowMeansLog10DatasetWo0_2), max(rowMeansLog10DatasetWo0_1,rowMeansLog10DatasetWo0_2)),
     xlab=PrefixInput1,
     ylab=PrefixInput2
)
#
### highlight MT- and mt-
points(x = rowMeansLog10DatasetWo0_1[mito.genes], y=rowMeansLog10DatasetWo0_2[mito.genes], col="red")
MaxOfPlot    <-max(rowMeansLog10DatasetWo0_1,rowMeansLog10DatasetWo0_2)
MinOfPlot    <-min(rowMeansLog10DatasetWo0_1,rowMeansLog10DatasetWo0_2)
MiddleOfPlot <-mean(c(MinOfPlot,MaxOfPlot))
MitoXposition<-MiddleOfPlot + ((MaxOfPlot - MiddleOfPlot) * 0.7)
MitoYposition<-MiddleOfPlot + ((MiddleOfPlot - MinOfPlot) * 0.75)
text(x = MitoXposition, y= MitoYposition, labels = "mito.genes", col = "red")
text(x = MinOfPlot+0.5, y= MaxOfPlot-0.5, pos=3, labels = paste("r=",r.cor.log10,sep = "", collapse = ""))
text(x = MinOfPlot+0.5, y= MaxOfPlot-0.5, pos=1,labels = paste("s=",s.cor.log10,sep = "", collapse = ""))
dev.off()
#
#################
##### Generate table log10(mean gene expression)
#################
label_av_1    <- paste(PrefixInput1, "_av", collapse = "", sep ="")
label_log10_1 <- paste(PrefixInput1, "_log10", collapse = "", sep ="")
label_av_2    <- paste(PrefixInput2, "_av", collapse = "", sep ="")
label_log10_2 <- paste(PrefixInput2, "_log10", collapse = "", sep ="")
#
mat<- data.frame("Gene" = rownames(input.matrix_1),
                 V1 = rowMeans(input.matrix_1),
                 V2 = log10(rowMeans(input.matrix_1)),
                 V3 = rowMeans(input.matrix_2),
                 V4 = log10(rowMeans(input.matrix_2))
)
colnames(mat) <- c("Gene", label_av_1, label_log10_1, label_av_2, label_log10_2)
write.table(file = paste(PrefixCombOutfiles, "_raw_counts.tsv", sep = ""), x=mat[rev(order(mat[1])),], row.names = F, sep = "\t", quote = F)
#
#################
##### Generate overlapping histograms
#################
h1 <- hist(rowMeansLog10DatasetWo0_1,breaks=30, plot = F)
h2 <- hist(rowMeansLog10DatasetWo0_2,breaks=30, plot = F)
pdf(file=paste(PrefixCombOutfiles, ".hist.pdf", sep=""))
plot( h1, col=rgb(0,0,1,1/4), border= rgb(0,0,1,1/4), xlim=c(MinOfPlot,MaxOfPlot), main = "",xlab = "log10(mean gene expression)")
plot( h2, col=rgb(1,0,0,1/4), border= rgb(1,0,0,1/4), xlim=c(MinOfPlot,MaxOfPlot), add=T)
legend("topright",legend = c(PrefixInput1, PrefixInput2), pch=15, col=c(rgb(0,0,1,1/3),rgb(1,0,0,1/3)), bty ="n", cex=1.3)
dev.off()



####################################
####################################
####################################
### Create summary plots outfiles
####################################
####################################
####################################
#
writeLines("\n*** Create summary plots outfile ***\n")
#
if (regexpr("^y$", SummaryPlots, ignore.case = T)[1] == 1) {
  suppressPackageStartupMessages(library(staplr))
  #
  SummaryPlotsPdf_1<-paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput1, ".SEURAT_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge_1<-c(
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput1, ".SEURAT_VlnPlot.seurat_filtered.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput1, ".SEURAT_Heatmap.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput1, ".SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput1, ".SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput1, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
  )
  #
  SummaryPlotsPdf_2<-paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput2, ".SEURAT_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge_2<-c(
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput2, ".SEURAT_VlnPlot.seurat_filtered.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput2, ".SEURAT_Heatmap.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput2, ".SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput2, ".SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles,  "_", PrefixInput2, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
  )
  #
  SummaryPlotsPdf_c12<-paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge_c12<-c(
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_VlnPlot.seurat_filtered.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_Heatmap.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
  )
  #
  SummaryPlotsPdf_1and2_c12<-paste(PrefixCombOutfiles, ".SEURAT_AND_CIRCOS.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge_1and2_c12<-c(
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_TSNEPlot_ColourByDataset.pdf", sep = "", collapse = ""),
    paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_CoordsFromMerged_ColoursFromItself.pdf", sep=""),
    paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_CoordsFromMerged_ColoursFromItself.pdf", sep=""),
    paste(Tempdir,"/",PrefixOutfiles, "_CombDatasets.circos.pdf", sep=""),
    paste(Tempdir,"/",PrefixOutfiles, "_merged.SEURAT_TSNEPlot_QC.pdf", sep=""),
    paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_QC_CoordsFromMerge.pdf", sep=""),
    paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_QC_CoordsFromMerge.pdf", sep=""),
    paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_CellsVennDiagram.pdf", sep=""),
    paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_GenesVennDiagram.pdf", sep="")
  )
  #
  #### Adding t-SNEs of requested genes (if applicable)
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("Create summary files")
  }else{
    ListOfPdfFilesToMerge_1         <-c(ListOfPdfFilesToMerge_1,         paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput1, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    ListOfPdfFilesToMerge_2         <-c(ListOfPdfFilesToMerge_2,         paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput2, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    ListOfPdfFilesToMerge_c12       <-c(ListOfPdfFilesToMerge_c12,       paste(Tempdir,"/",PrefixOutfiles,"_merged.SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    print("Create summary files including t-SNE's for selected genes")
  }
  
  ### Stappling *pdf's
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge_1,         output_filepath = SummaryPlotsPdf_1)
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge_2,         output_filepath = SummaryPlotsPdf_2)
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge_c12,       output_filepath = SummaryPlotsPdf_c12)
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge_1and2_c12, output_filepath = SummaryPlotsPdf_1and2_c12)
}

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")
OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))
#
for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}
#
####################################
### Report time used
####################################
writeLines("\n*** Report time used ***\n")
#
EndTimeOverall<-Sys.time()
#
TookTimeClustering_1       <-format(difftime(EndTimeClustering_1,     StartTimeClustering_1,     units = "min"))
TookTimeFindAllMarkers_1   <-format(difftime(EndTimeFindAllMarkers_1, StartTimeFindAllMarkers_1, units = "min"))
#
TookTimeClustering_2       <-format(difftime(EndTimeClustering_2,     StartTimeClustering_2,     units = "min"))
TookTimeFindAllMarkers_2   <-format(difftime(EndTimeFindAllMarkers_2, StartTimeFindAllMarkers_2, units = "min"))
#
TookTimeClustering_c12     <-format(difftime(EndTimeClustering_c12,     StartTimeClustering_c12,     units = "min"))
TookTimeFindAllMarkers_c12 <-format(difftime(EndTimeFindAllMarkers_c12, StartTimeFindAllMarkers_c12, units = "min"))
#
TookTimeOverall            <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "min"))
#
OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CPUusage.txt", sep="")
ReportTime<-c(
  paste("clustering_dataset_1",TookTimeClustering_1,collapse = "\t"),
  paste("FindAllMarkers_dataset_1",TookTimeFindAllMarkers_1,collapse = "\t"),
  paste("clustering_dataset_2",TookTimeClustering_2,collapse = "\t"),
  paste("FindAllMarkers_dataset_2",TookTimeFindAllMarkers_2,collapse = "\t"),
  paste("clustering_dataset_c12",TookTimeClustering_c12,collapse = "\t"),
  paste("FindAllMarkers_dataset_c12",TookTimeFindAllMarkers_c12,collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)
#
write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"/SEURAT/",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir,pattern = PrefixOutfiles, full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/SEURAT/",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

####################################
### Remove temporary files
####################################

if (ToRemoveDowngraded1 == 1) {
  system(command = paste("rm -r ", NewInput1, sep = "", collapse = ""))
}

if (ToRemoveDowngraded2 == 1) {
  system(command = paste("rm -r ", NewInput2, sep = "", collapse = ""))
}

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

