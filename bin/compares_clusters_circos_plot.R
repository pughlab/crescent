####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://xgaoo.github.io/ClusterMap/ClusterMap.html
###
### NEW IMPLEMENTATIONS SINCE Seurat v2:
### 1) Rewritten with Seurat v3 commands (including ability to read output from Cell Ranger v3)
###    Main differences vs. Seurat v2 include:
###    a) new function names
###    b) a new function DimPlot() is used for 'uma', 'tsne' and 'pca' plots, instead of PCAPlot() and TSNEPlot()
###       Hence, now we define 'dimensions' instead of principal components
###    c) since Cell Ranger v3 allows now to have multiple features (not only genes),
###       in general all references to 'genes' in v2 are now called 'features'
### 2) Default parameters are indexed using list(DefaultParameters) instead of variable names. This allows list(DefaultParameters)
###    to be reported in a *log file
### 3) Implemented QC violin plots directly in ggplots instead of Seurat's VlnPlot() to have more control of layout, legends, titles, etc.
### 4) Using R base plot() to create feature-vs-feature scatter plots instead of FeatureScatter()
### 5) Implemented all ggplots and Seurat plotting functions (which are based on ggplots) with print() function, like:
###    `print(FeaturePlot(...))` instead of `FeaturePlot(...)` alone
###    Otherwise using ggplots and Seurat plots inside if/else loops cause errors
###    https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-do-lattice_002ftrellis-graphics-not-work_003f
### 6) Violin plots and t-SNE QC plots show now percentage of ribosomal protein genes in cells
###
### THINGS TO DO:
### 1) Pick the right number of dimension components
###    E.g. try to automatically get the inflection point from the PCElbowPlot() function output
### 2) Pick the right resolution from FindClusters()
###    E.g. implement Iness and Bader paper https://f1000research.com/articles/7-1522/v1 approach
###    By picking the number of clusters based on differentially expressed genes
### 3) Add a lists of ENSEMBL Ids for mitochondrial genes instead of just MT- and mt- (at gene names)
###    Need to do it for both Human and Mouse
### 4) In Seurat v2, Suluxan reported that the -c example Javier provided called a duplicated row names error
###    Need to see if it's still happening
###
### THINGS NICE TO HAVE:
### 1) Assigning cell type identity to clusters (needs supervised annotations, maybe based on GSVA)
### 2) Use knitr() to produce html plot layouts (https://yihui.name/knitr/demo/stitch/)
### 3) In "Load data" we use Seurat(Read10X). In this command, when all barcodes come from the same sample (i.e. finish with the same digit), like:
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
###    Instead of:
###    CTCTACGCAAGAGGCT-1
###    CTCGAAAAGCTAACAA-2
###    CTGCCTAGTGCAGGTA-3
###   
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
### Seurat v3 can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table - only needed is using '-t DGE'
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # To run parallel processes
suppressPackageStartupMessages(library(ClusterMap))   # to draw circos and other plots. Install as: library(devtools), then install_github('xgaoo/ClusterMap')
suppressPackageStartupMessages(library(VennDiagram))  # to get gene and cell Venn diagrams
### library(staplr)     only if using option '-s y', note it needs pdftk
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
              help="Path/name to:
                a) a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files; or
                b) a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns; or
                c) a *RDS Seurat object with pre-computed cell clusters and differentially expressed genes for each cell cluster

                Notes:
                The 'MTX' files can be for example the output from Cell Ranger `/path_to/outs/filtered_feature_bc_matrix/`
                The 'DGE' file can be for example the output from Dropseq tools
                The 'RDS' file can be for example the *RDS output from Runs_Seurat_v3.R

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-j", "--input_2"), default="NA",
              help="Same as -i but indicating the second datset to compare
                NOTE: the two datasets to compare can have different genes, but this program will work on their intersection

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--input_type_1"), default="NA",
              help="Indicates if input is either a 'MTX' directory or a 'DGE' file or a 'RDS' file

                Default = 'No default. It's mandatory to specify this parameter'"),
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
  make_option(c("-m", "--percent_mito"), default="0,0.05",
              help="<comma> delimited min,max number of percentage of mitochondrial gene counts in a cell to be included in normalization and clustering analyses
                Use 'Inf' as min if no minumum limit should be used, e.g. '0,0.05'"),
  #
  make_option(c("-n", "--n_genes"), default="50,8000",
              help="<comma> delimited min,max number of unique gene counts in a cell to be included in normalization and clustering analyses. E.g '50,8000'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh,  e.g. '0.01'"),
  
  make_option(c("-f", "--edge_cutoff"), default="0.1",
              help="The edge length cutoff to decide the sub-nodes to merge or not between two the two datasets
                Also controls the cutoff to draw lines between clusters in the circos plot
                Default is 0.1. If lines are too weak you can lower it e.g. 0.05"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system
                Default = 'MAX'")
  
)

opt <- parse_args(OptionParser(option_list=option_list))

Input_d1         <- opt$input_1
Input_d2         <- opt$input_2
InputType_d1     <- opt$input_type_1
InputType_d2     <- opt$input_type_2
Resolution_d1    <- as.numeric(opt$resolution_1) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Resolution_d2    <- as.numeric(opt$resolution_2) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
PrefixInput_d1   <- opt$prefix_input_1
PrefixInput_d2   <- opt$prefix_input_2
Outdir           <- opt$outdir
PrefixOutfiles   <- opt$prefix_outfiles
SummaryPlots     <- opt$summary_plots
ListGenes        <- opt$list_genes
InfileColourTsne <- opt$infile_colour_tsne
Opacity          <- as.numeric(opt$opacity)
PcaDimsUse       <- c(1:as.numeric(opt$pca_dimensions))
ListPMito        <- opt$percent_mito
ListNGenes       <- opt$n_genes
ThreshReturn     <- as.numeric(opt$return_threshold)
EdgeCutoff       <- as.numeric(opt$edge_cutoff)
NumbCores        <- opt$number_cores

Input_d1         <- "~/SINGLE_CELL/10X/KETELA_DATA/190430_NB501085_0333_AHCY7JBGXB_Pugh_Laura_Javier/Pugh_Laura_Javier__GBM1070r/filtered_feature_bc_matrix/"
Input_d2         <- "~/SINGLE_CELL/10X/KETELA_DATA/190430_NB501085_0333_AHCY7JBGXB_Pugh_Laura_Javier/Pugh_Laura_Javier__GBM1070r_4hr_test/filtered_feature_bc_matrix"
InputType_d1     <- "MTX"
InputType_d2     <- "MTX"
Resolution_d1    <- 0.8
Resolution_d2    <- 0.8
PrefixInput_d1   <- "0hr"
PrefixInput_d2   <- "4hr"
Outdir           <- "~/SINGLE_CELL/10X/KETELA_DATA/190430_NB501085_0333_AHCY7JBGXB_Pugh_Laura_Javier/CIRCOS/"
PrefixOutfiles   <- "GBM1070r_0hr_vs_4hr"
SummaryPlots     <- "y"
ListGenes        <- "MALAT1,GAPDH"
InfileColourTsne <- "NA"
Opacity          <- 0.3
PcaDimsUse       <- c(1:10)
ListPMito        <- "0,0.05"
ListNGenes       <- "50,8000"
ThreshReturn     <- 0.01
EdgeCutoff       <- 0.1
NumbCores        <- 4

Resolution_c12   <- max(Resolution_d1,Resolution_d2)
PrefixInput_c12  <- "merged"

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

####################################
### Define the number of cores for parallelization
####################################

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (is.numeric(NumbCores) == T) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

cat("Using ", NumbCoresToUse, "cores")

plan(strategy = "multicore", workers = NumbCoresToUse)

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 850mb of RAM (850*1024^2 = 891289600), use:
### https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r
options(future.globals.maxSize= 891289600)

####################################
### Define default parameters
####################################
### Some of these default parameters are provided by Seurat developers,
### others are tailored according to clusters/t-SNE granularity

ListNGenes = unlist(strsplit(ListNGenes, ","))
MinGenes   = as.numeric(ListNGenes[1])
MaxGenes   = as.numeric(ListNGenes[2])
#
ListPMito  = unlist(strsplit(ListPMito,  ","))
MinPMito   = as.numeric(ListPMito[1])
MaxPMito   = as.numeric(ListPMito[2])

DefaultParameters <- list(
  
  ### Parameters for Seurat filters
  MinCells = 3,
  MinGenes = MinGenes,
  MaxGenes = MaxGenes,
  MinPMito = MinPMito,
  MaxPMito = MaxPMito,
  
  ### Parameters for Seurat normalization
  ScaleFactor = 10000,
  
  ### Parameters for Seurat variable gene detection
  XLowCutoff = 0.0125,
  XHighCutoff = 3,
  YCutoff = 0.5,
  
  ### Parameters for PCA
  PrintPCA.PcsPrint = 1:5,
  PrintPCA.GenesPrint = 5,
  VizPCA.PcsUse = 1:6,
  VizPCA.nGenesToPlot = 20,
  PCHeatmapCellsUse = 300,
  PCHeatmapComponentsToPlot = 18,
  
  ### Parameters for Cluster Biomarkers
  FindAllMarkers.MinPct     =  0.25,
  FindAllMarkers.ThreshUse  =  0.25,
  NumberOfGenesPerClusterToPlotTsne  =  2,
  NumberOfGenesPerClusterToPlotHeatmap  =  10,
  
  ### Parameters for t-SNE plots
  BaseSizeSingleTnePlot  = 7,
  BaseSizeMultipleWidth  = 3.7,
  BaseSizeMultipleHeight = 3
)

### Colour definitions
ColourDefinitions<-list("orange"        = "#E69F00",
                        "bluishgreen"   = "#009E73",
                        "reddishpurple" = "#CC79A7",
                        "dodgerblue"    = "#1E90FF",
                        "vermillion"    = "#D55E00",
                        "snow4"         = "#8B8989",
                        "yellow"        = "#FFD700",
                        "seagreen3"     = "#43CD80",
                        "pink"          = "#FFC0CB",
                        "skyblue"       = "#56B4E9",
                        "orchid4"       = "#8B4789",
                        "blue"          = "#0072B2",
                        "black"         = "#000000"
)
ColoursQCViolinPlots <- c(ColourDefinitions[["skyblue"]][[1]], ColourDefinitions[["orange"]][[1]])

# Parameters for Colouring datasets
ColorDataset1Tsne<-rgb(0.90,0.60,0,0.3) # orange with strong transparency
ColorDataset2Tsne<-rgb(0.35,0.70,0.90,0.3) # blue with strong transparency
ColorDataset1Venn<-rgb(0.90,0.60,0,0.7) # orange with weak transparency
ColorDataset2Venn<-rgb(0.35,0.70,0.90,0.7) # blue with weak transparency

####################################
### Start stopwatch
####################################

StartTimeOverall <-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input_d1", "input_d2", "input_type_d1", "input_type_d2", "outdir", "prefix_outfiles")
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
### Load scRNA-seq data and create UNfiltered Seurat objects
####################################

writeLines("\n*** Load scRNA-seq data and create UNfiltered Seurat objects ***\n")

writeLines("\n*** Load scRNA-seq data 1 ***\n")
if (regexpr("^MTX$", InputType_d1, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  input.matrix_d1 <- Read10X(data.dir = Input_d1)
}else if (regexpr("^DGE$", InputType_d1, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  input.matrix_d1 <- data.frame(fread(Input_d1),row.names=1)
}else{
  stop(paste("Unexpected type of input: ", InputType_d1, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

writeLines("\n*** Load scRNA-seq data 2 ***\n")
if (regexpr("^MTX$", InputType_d2, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  input.matrix_d2 <- Read10X(data.dir = Input_d2)
}else if (regexpr("^DGE$", InputType_d2, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  library(data.table)
  input.matrix_d2 <- data.frame(fread(Input_d2),row.names=1)
}else{
  stop(paste("Unexpected type of input: ", InputType_d2, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

####################################
### Create separate dataset Seurat objects
####################################
writeLines("\n*** Create separate dataset Seurat objects ***\n")

seurat.object_d1.u <- CreateSeuratObject(counts = input.matrix_d1, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = PrefixInput_d1)
seurat.object_d2.u <- CreateSeuratObject(counts = input.matrix_d2, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = PrefixInput_d2)

nCellsInOriginalMatrix_d1<-length(seurat.object_d1.u@meta.data$orig.ident)
nCellsInOriginalMatrix_d2<-length(seurat.object_d2.u@meta.data$orig.ident)

####################################
### Get intersection of genes from the two datasets
####################################
writeLines("\n*** Get intersection of genes from the two datasets ***\n")

listGenes_Union         <- union(rownames(seurat.object_d1.u),rownames(seurat.object_d2.u))
listGenes_Intersection  <- intersect(rownames(seurat.object_d1.u),rownames(seurat.object_d2.u))
listGenes_Complement1   <- setdiff(rownames(seurat.object_d1.u),rownames(seurat.object_d2.u))
listGenes_Complement2   <- setdiff(rownames(seurat.object_d2.u),rownames(seurat.object_d1.u))

subset.matrix_d1     <- seurat.object_d1.u@assays$RNA@counts[listGenes_Intersection, ]
seurat.object_d1.ig  <- CreateSeuratObject(counts = subset.matrix_d1, project = PrefixInput_d1)
subset.matrix_d2     <- seurat.object_d2.u@assays$RNA@counts[listGenes_Intersection, ]
seurat.object_d2.ig  <- CreateSeuratObject(counts = subset.matrix_d2, project = PrefixInput_d2)

####################################
### Get intersection of cell barcodes from the two datasets
####################################
writeLines("\n*** Get intersection of cell barcodes from the two datasets ***\n")

listCells_Union         <- union(colnames(seurat.object_d1.u),colnames(seurat.object_d2.u))
listCells_Intersection  <- intersect(colnames(seurat.object_d1.u),colnames(seurat.object_d2.u))
listCells_Complement1   <- setdiff(colnames(seurat.object_d1.u),colnames(seurat.object_d2.u))
listCells_Complement2   <- setdiff(colnames(seurat.object_d2.u),colnames(seurat.object_d1.u))

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
  area1 = length(row.names(seurat.object_d1.u@assays$RNA@counts)),
  area2 = length(row.names(seurat.object_d2.u@assays$RNA@counts)),
  cross.area = length(listGenes_Intersection),
  category = c(PrefixInput_d1, PrefixInput_d2),
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
  area1 = length(colnames(seurat.object_d1.u@assays$RNA@counts)),
  area2 = length(colnames(seurat.object_d2.u@assays$RNA@counts)),
  cross.area = length(listCells_Intersection),
  category = c(PrefixInput_d1, PrefixInput_d2),
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
write.table(paste("GeneID", PrefixInput_d1, PrefixInput_d2, sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste(listGenes_Intersection,"yes","yes", sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement1,"yes","no"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement2,"no","yes"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)

OutfileGenesPerDataset<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CellsPerDataset.tsv", sep="")
write.table(paste("CellBarcode", PrefixInput_d1, PrefixInput_d2, sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste(listGenes_Intersection,"yes","yes", sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement1,"yes","no"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(paste(listGenes_Complement2,"no","yes"  , sep="\t", collapse = "\n"), file = OutfileGenesPerDataset, row.names = F, col.names = F, sep="\t", quote = F, append = T)

####################################
### Create merged Seurat object, still UNfiltered by mito and ribo contents (that will be done below)
####################################

writeLines("\n*** Create merged Seurat object ***\n")

print("Original datasets:")
print("Dataset 1:")
seurat.object_d1.u
print("Dataset 2:")
seurat.object_d2.u

print("Gene intersection datasets:")
print("Dataset 1:")
seurat.object_d1.ig
print("Dataset 2:")
seurat.object_d2.ig

input.matrix_d1<-seurat.object_d1.u@assays$RNA@counts
input.matrix_d2<-seurat.object_d2.u@assays$RNA@counts

colnames(input.matrix_d1)=paste0('Dataset1-', colnames(input.matrix_d1))
colnames(input.matrix_d2)=paste0('Dataset2-', colnames(input.matrix_d2))
CombDatasets.data = cbind(input.matrix_d1[listGenes_Intersection,], input.matrix_d2[listGenes_Intersection,])

### Note the unfiltered *u and filtered *ig seurat objects for the merged dataset are identical
### This is because to create it we are restrincting to intersected genes
seurat.object_c12.u   <- CreateSeuratObject(counts = CombDatasets.data, project = PrefixOutfiles)
seurat.object_c12.ig  <- seurat.object_c12.u

print("Merged dataset:")
seurat.object_c12.u
seurat.object_c12.ig
nCellsInOriginalMatrix_c12<-length(colnames(seurat.object_c12.ig@assays$RNA@counts))

####################################
### Get mitochondrial genes
####################################
writeLines("\n*** Get  mitochondrial genes ***\n")

mitoRegExpressions<- paste(c("^MT-"),collapse = "|")
mito.features <- grep(pattern = mitoRegExpressions, ignore.case = T, x = unique(c(rownames(seurat.object_d1.u),rownames(seurat.object_d2.u))), value = T)

if (length(mito.features)[[1]] > 0) {
  percent.mito_d1 <-   Matrix::colSums(x = GetAssayData(object = seurat.object_d1.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object_d1.u, slot = 'counts'))
  seurat.object_d1.u[['percent.mito']] <-   percent.mito_d1
  percent.mito_d2 <-   Matrix::colSums(x = GetAssayData(object = seurat.object_d2.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object_d2.u, slot = 'counts'))
  seurat.object_d2.u[['percent.mito']] <-   percent.mito_d2
  percent.mito_c12 <- Matrix::colSums(x = GetAssayData(object = seurat.object_c12.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object_c12.u, slot = 'counts'))
  seurat.object_c12.u[['percent.mito']] <- percent.mito_c12
  }else{
  percent.mito_d1 <-   0 / Matrix::colSums(x = GetAssayData(object = seurat.object_d1.u, slot = 'counts'))
  seurat.object_d1.u[['percent.mito']] <- percent.mito_d1
  percent.mito_d2 <-   0 / Matrix::colSums(x = GetAssayData(object = seurat.object_d2.u, slot = 'counts'))
  seurat.object_d2.u[['percent.mito']] <- percent.mito_d2
  percent.mito_c12 <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object_c12.u, slot = 'counts'))
  seurat.object_c12.u[['percent.mito']] <- percent.mito_c12
}

####################################
### Get ribosomal protein genes
####################################
writeLines("\n*** Get ribosomal protein genes ***\n")

riboRegExpressions<- paste(c("^MRPL", "^MRPS", "^RPL", "^RPS"),collapse = "|")
ribo.features <- grep(pattern = riboRegExpressions, ignore.case = T, x = rownames(x = seurat.object_c12.ig), value = T)

if (length(ribo.features)[[1]] > 0) {
  percent.ribo_d1  <- Matrix::colSums(x = GetAssayData(object = seurat.object_d1.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object_d1.u, slot = 'counts'))
  seurat.object_d1.u[['percent.ribo']] <- percent.ribo_d1
  percent.ribo_d2  <- Matrix::colSums(x = GetAssayData(object = seurat.object_d2.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object_d2.u, slot = 'counts'))
  seurat.object_d2.u[['percent.ribo']] <- percent.ribo_d2
  percent.ribo_c12 <- Matrix::colSums(x = GetAssayData(object = seurat.object_c12.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object_c12.u, slot = 'counts'))
  seurat.object_c12.u[['percent.ribo']] <- percent.ribo_c12
}else{
  percent.ribo_d1  <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object_d1.u, slot = 'counts'))
  seurat.object_d1.u[['percent.ribo']] <- percent.ribo_d1
  percent.ribo_d2  <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object_d2.u, slot = 'counts'))
  seurat.object_d2.u[['percent.ribo']] <- percent.ribo_d2
  percent.ribo_c12 <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object_c12.u, slot = 'counts'))
  seurat.object_c12.u[['percent.ribo']] <- percent.ribo_c12
}

####################################
### Filter cells based gene counts and mitochondrial representation
####################################
writeLines("\n*** Filter cells based gene counts and mitochondrial representation ***\n")

if (length(mito.features)[[1]] > 0) {
  seurat.object_d1.f  <-subset(x = seurat.object_d1.ig,  subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes & percent.mito > DefaultParameters$MinPMito & percent.mito < DefaultParameters$MaxPMito)
  seurat.object_d2.f  <-subset(x = seurat.object_d2.ig,  subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes & percent.mito > DefaultParameters$MinPMito & percent.mito < DefaultParameters$MaxPMito)
  seurat.object_c12.f <-subset(x = seurat.object_c12.ig, subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes & percent.mito > DefaultParameters$MinPMito & percent.mito < DefaultParameters$MaxPMito)
}else{
  seurat.object_d1.f  <-subset(x = seurat.object_d1.ig,  subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes)
  seurat.object_d2.f  <-subset(x = seurat.object_d2.ig,  subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes)
  seurat.object_c12.f <-subset(x = seurat.object_c12.ig, subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes)
}





####################################
####################################
####################################
### Process dataset 1
####################################
####################################
####################################
writeLines("\n*** Process dataset 1 ***\n")

####################################
### QC EDA violin plots dataset 1
####################################
writeLines("\n*** QC EDA violin plots dataset 1 ***\n")

### Get unfiltered data QC statistics
nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object_d1.u@meta.data$nFeature_RNA, nGenes = 1)
nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object_d1.u@meta.data$nCount_RNA,   nCount_RNA = 1)
percent.mito.u.df  <-data.frame(Expression_level = seurat.object_d1.u@meta.data$percent.mito, percent.mito = 1)
percent.ribo.u.df  <-data.frame(Expression_level = seurat.object_d1.u@meta.data$percent.ribo, percent.ribo = 1)
#
nFeature_RNAStats.u<-paste(c(" mean = ",round(mean(seurat.object_d1.u@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object_d1.u@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.u  <-paste(c( "mean = ",round(mean(seurat.object_d1.u@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object_d1.u@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.u     <-paste(c(" mean = ",round(mean(seurat.object_d1.u@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object_d1.u@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
percent.ribo.u     <-paste(c(" mean = ",round(mean(seurat.object_d1.u@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object_d1.u@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")

### Get filtered data QC statistics
nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object_d1.f@meta.data$nFeature_RNA, nGenes = 2)
nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object_d1.f@meta.data$nCount_RNA,   nCount_RNA = 2)
percent.mito.f.df  <-data.frame(Expression_level = seurat.object_d1.f@meta.data$percent.mito, percent.mito = 2)
percent.ribo.f.df  <-data.frame(Expression_level = seurat.object_d1.f@meta.data$percent.ribo, percent.ribo = 2)
#
nFeature_RNAStats.f<-paste(c(" mean = ",round(mean(seurat.object_d1.f@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object_d1.f@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.f  <-paste(c(" mean = ",round(mean(seurat.object_d1.f@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object_d1.f@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.f     <-paste(c(" mean = ",round(mean(seurat.object_d1.f@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object_d1.f@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
percent.ribo.f     <-paste(c(" mean = ",round(mean(seurat.object_d1.f@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object_d1.f@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")

### Put QC statistics together
nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
percent.mito.m.df  <-data.frame(rbind(percent.mito.u.df,percent.mito.f.df))
percent.ribo.m.df  <-data.frame(rbind(percent.ribo.u.df,percent.ribo.f.df))
LabelUnfiltered    <-paste("Before filters: No. of cells = ", nrow(seurat.object_d1.u@meta.data), sep ="", collapse = "")
LabelFiltered      <-paste("After filters:  No. of cells = ", nrow(seurat.object_d1.f@meta.data), sep ="", collapse = "")
#LabelFilters       <-paste("Filters: No. of genes per cell = ", ListNGenes, "Fraction mitochondrial protein genes per cell = ", ListPMito), sep ="", collapse = "")

### Commands for violin ggplot's
DataForHeader.df<-data.frame(forx = c(0.4,0.4), fory = c(0.09,0.03), label = c(LabelFiltered,LabelUnfiltered))
Headers.plot<-ggplot(data=DataForHeader.df, aes(x = forx, y = fory)) + theme_void() + 
  geom_point(colour = ColoursQCViolinPlots, size = 7) + xlim(0,1) + ylim(0,0.12) +
  geom_text(hjust = -0.08, label = c(LabelUnfiltered,LabelFiltered), size = 4)

nFeature_RNA.plot<-ggplot(data=nFeature_RNA.m.df, aes(x = factor(nGenes), y = Expression_level)) +
  geom_violin(aes(fill = factor(nGenes))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="No. of genes") +
  annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])

nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
  geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="No. of reads") +
  annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])

percent.mito.plot<-ggplot(data=percent.mito.m.df, aes(x = factor(percent.mito), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.mito))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Mitochondrial genes (fraction)") +
  annotate("text", x = 1 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.f, col = ColoursQCViolinPlots[[2]])

percent.ribo.plot<-ggplot(data=percent.ribo.m.df, aes(x = factor(percent.ribo), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.ribo))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Ribosomal protein genes (fraction)") +
  annotate("text", x = 1 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.f, col = ColoursQCViolinPlots[[2]])

bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, percent.mito.plot, percent.ribo.plot, ncol = 4)

### Create a *pdf file with the violin ggplot's

VlnPlotPdf<-paste(Tempdir, "/" , PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_QC_VlnPlot.pdf", sep="")
pdf(file=VlnPlotPdf, width = 12, height = 7)
print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
dev.off()

####################################
### Feature-vs-feature scatter plots
####################################
writeLines("\n*** Feature-vs-feature scatter plots dataset 1 ***\n")

UnfilteredData.df<-data.frame(nCount_RNA = seurat.object_d1.u@meta.data$nCount_RNA,
                              nGene = seurat.object_d1.u@meta.data$nFeature_RNA,
                              percent.mito = seurat.object_d1.u@meta.data$percent.mito,
                              filtered_out = colnames(seurat.object_d1.u) %in% colnames(seurat.object_d1.f))
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = ColoursQCViolinPlots[[1]])
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = FALSE, replacement = ColoursQCViolinPlots[[2]])
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = 4)
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = FALSE, replacement = 16)

FeatureVsFeaturePlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_NumbReadsVsNumbGenesAndMito_VlnPlot.pdf", sep="")
pdf(file=FeatureVsFeaturePlotPdf, width = 10, height = 5)
par(mfrow=c(1,2))
## No. of reads vs. Mitochond. %
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$percent.mito, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "Mitochond. %")
legend("topright", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))

## No. of reads vs. No. of genes
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$nGene, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "No. of genes")
legend("topleft", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))
dev.off()

####################################
### Remove the Unfiltered seurat object
####################################
writeLines("\n*** Remove the Unfiltered seurat object dataset 1 ***\n")

rm(seurat.object_d1.u)
rm(UnfilteredData.df)

####################################
### Normalize data
####################################
writeLines("\n*** Normalize data dataset 1 ***\n")

seurat.object_d1.f <- NormalizeData(object = seurat.object_d1.f, normalization.method = "LogNormalize", scale.factor = DefaultParameters$ScaleFactor)

####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes dataset 1 ***\n")

### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object_d1.f <- FindVariableFeatures(object = seurat.object_d1.f, selection.method = 'mean.var.plot', mean.cutoff = c(DefaultParameters$XLowCutoff, DefaultParameters$XHighCutoff), dispersion.cutoff = c(DefaultParameters$YCutoff, Inf))
VariableGenes<-VariableFeatures(object = seurat.object_d1.f)
length(VariableGenes)

write(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_VariableGenes.txt", sep=""), x=VariableGenes)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_VariableGenes.pdf", sep=""))
print(VariableFeaturePlot(object = seurat.object_d1.f, cols = c("blue", "red")))
dev.off()

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation dataset 1 ***\n")

seurat.object_d1.f <- ScaleData(object = seurat.object_d1.f, vars.to.regress = c("nCount_RNA", "percent.mito"), features = rownames(x = seurat.object_d1.f), display.progress=F)

####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
### Note: DimPlot can now handle 'umap' and 'tsne' in addition to 'pca', but for 'umap'
### you must first install the umap-learn python package (e.g. via pip install umap-learn)
### https://github.com/satijalab/seurat/blob/master/R/dimensional_reduction.R
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA dataset 1 ***\n")

seurat.object_d1.f <- RunPCA(object = seurat.object_d1.f, features = VariableGenes, verbose = T, do.print = T, ndims.print = DefaultParameters$PrintPCA.PcsPrint, nfeatures.print = DefaultParameters$PrintPCA.GenesPrint)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
print(VizDimLoadings(object = seurat.object_d1.f, reduction = "pca", dims = DefaultParameters$VizPCA.PcsUse, nfeatures = DefaultParameters$VizPCA.nGenesToPlot))
dev.off()

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_PCAPlot.pdf", sep=""))
print(DimPlot(object = seurat.object_d1.f, dims = c(1,2), reduction = "pca") + theme(legend.position="none"))
dev.off()

seurat.object_d1.f <- ProjectDim(object = seurat.object_d1.f, overwrite = T, verbose = T, nfeatures.print = 10)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_PCHeatmap.C1to",DefaultParameters$PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
print(DimHeatmap(object = seurat.object_d1.f, dims = 1:DefaultParameters$PCHeatmapComponentsToPlot, cells = DefaultParameters$PCHeatmapCellsUse, balanced = T))
dev.off()

####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components dataset 1 ***\n")

### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_PCElbowPlot.pdf", sep=""))
print(ElbowPlot(object = seurat.object_d1.f))
dev.off()

####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells dataset 1 ***\n")

StartTimeClustering<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object_d1.f <- FindNeighbors(object = seurat.object_d1.f, dims = PcaDimsUse) ## This step was part of FindClusters() in Seurat v2
seurat.object_d1.f <- FindClusters(object = seurat.object_d1.f, reduction.type = "pca", resolution = Resolution_d1, print.output = 0, save.SNN = T)
EndTimeClustering<-Sys.time()

CellNames<-rownames(seurat.object_d1.f@meta.data)
ClusterIdent<-seurat.object_d1.f@active.ident
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution_d1, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)

####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene dataset 1 ***\n")

cluster.averages<-AverageExpression(object = seurat.object_d1.f, use.scale = F, use.counts = F)
OutfileClusterAverages<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("r", Resolution_d1, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
writeLines("\n*** Run Non-linear dimensional reduction (tSNE) dataset 1 ***\n")

### NOTE: if the datasets is small you may get
### "Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large."
### One can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10

seurat.object_d1.f <- RunTSNE(object = seurat.object_d1.f, dims.use = PcaDimsUse, do.fast = T)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot.pdf", sep=""), width = 7, height = 7)
print(DimPlot(object = seurat.object_d1.f, reduction = 'tsne', group.by = 'ident', label = T, label.size=10))
dev.off()

####################################
### Saving the R object
####################################
writeLines("\n*** Saving the R object dataset 1 ***\n")

OutfileRDS_d1<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_object.rds", sep="")
saveRDS(seurat.object_d1.f, file = OutfileRDS_d1)

####################################
### Write out t-SNE coordinates
####################################
writeLines("\n*** Write out t-SNE coordinates dataset 1 ***\n")

OutfileTsneCoordinates<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNECoordinates.tsv", sep="")

Headers<-paste("Barcode",paste(colnames(seurat.object_d1.f@reductions$tsne@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileTsneCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
write.table(seurat.object_d1.f@reductions$tsne@cell.embeddings, file = OutfileTsneCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Colour t-SNE by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo
####################################
writeLines("\n*** Colour t-SNE by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo dataset 1 ***\n")

CellPropertiesToTsne<-c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 21, height = 5)
print(FeaturePlot(object = seurat.object_d1.f, features = CellPropertiesToTsne, cols = c("lightgrey", "blue"), reduction = "tsne", ncol = 4, pt.size = 1.5))
dev.off()


####################################
### Colour t-SNE by -infile_colour_tsne
####################################
writeLines("\n*** Colour t-SNE by -infile_colour_tsne dataset 1 ***\n")

if (regexpr("^NA$", InfileColourTsne, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  seurat.object.meta.data<-seurat.object_d1.f@meta.data
  ExtraCellProperties <- data.frame(read.table(InfileColourTsne, header = T, row.names = 1))
  head(ExtraCellProperties)
  
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_tsne and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object_d1.f <- AddMetaData(object = seurat.object_d1.f, metadata = ExtraCellProperties)
  
  # Generating outfile
  # Note DimPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...)) or  par(mfrow=())
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""), height = DefaultParameters$BaseSizeSingleTnePlot, width = DefaultParameters$BaseSizeSingleTnePlot)
  for (property in colnames(ExtraCellProperties)) {
    print(DimPlot(object = seurat.object_d1.f, reduction = 'tsne', group.by = property, combine = T, legend = "none") + ggtitle(property))
  }
  dev.off()
}

####################################
### Colour t-SNE plots showing each requested gene
####################################
writeLines("\n*** Colour t-SNE plots showing each requested gene dataset 1 ***\n")

### To program layout() for more than 3 genes in multiple rows
if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
  print("No selected genes for t-SNE plots")
}else{
  ListOfGenesForTsnes<-unlist(strsplit(ListGenes, ","))
  if (length(ListOfGenesForTsnes) <= 4) {
    pdfWidth  <- (length(ListOfGenesForTsnes) * DefaultParameters$BaseSizeMultipleWidth)
    pdfHeight <- DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- length(ListOfGenesForTsnes)
  }else{
    pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
    pdfHeight <- (as.integer(length(ListOfGenesForTsnes) / 4) + 1) * DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- 4
  }
  
  pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  print(FeaturePlot(object = seurat.object_d1.f, ncol = nColFeaturePlot, features = c(ListOfGenesForTsnes), cols = c("lightgrey", "blue"), reduction = "tsne"))
  dev.off()
}

####################################
### Finding differentially expressed genes for each cell cluster
####################################
writeLines("\n*** Finding differentially expressed genes for each cell cluster dataset 1 ***\n")

### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
#########
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object_d1.f, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object_d1.f, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
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

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object_d1.f@meta.data))

StartTimeFindAllMarkers<-Sys.time()
seurat.object.markers_d1 <- FindAllMarkers(object = seurat.object_d1.f, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
EndTimeFindAllMarkers<-Sys.time()
#
OufileTsvMarkers_d1<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_MarkersPerCluster.tsv",sep="")
OufileCsvMarkers_d1<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_MarkersPerCluster.csv",sep="")
write.table(data.frame("GENE"=rownames(seurat.object.markers_d1),seurat.object.markers_d1), file = OufileTsvMarkers_d1, row.names = F,sep="\t",quote = F)
write.csv(seurat.object.markers_d1, file = OufileCsvMarkers_d1)

### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne<-(seurat.object.markers_d1 %>% group_by(cluster) %>% top_n(DefaultParameters$NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers_d1[["cluster"]]))

####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes dataset 1 ***\n")

NumberOfPanesForFeaturesPlot<-(NumberOfClusters*DefaultParameters$NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5

### Here need to program a better way to control the y-axis labels
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object_d1.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), slot = 'counts', log = T, adjust = 1, pt.size = 0.5))
dev.off()

####################################
### t-SNE plots showing each cluster top genes
####################################
writeLines("\n*** t-SNE plots showing each cluster top genes dataset 1 ***\n")

pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
pdfHeight <- NumberOfClusters * DefaultParameters$BaseSizeMultipleHeight / 2
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(FeaturePlot(object = seurat.object_d1.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), cols = c("lightgrey", "blue"), reduction = "tsne"))
dev.off()

####################################
### Cell clusters heatmap
####################################
writeLines("\n*** Cell clusters heatmap dataset 1 ***\n")

top_genes_by_cluster_for_heatmap <- seurat.object.markers_d1 %>% group_by(cluster) %>% top_n(n = DefaultParameters$NumberOfGenesPerClusterToPlotHeatmap, wt = avg_logFC)
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_Heatmap.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(DoHeatmap(object = seurat.object_d1.f, features = top_genes_by_cluster_for_heatmap$gene, label = T, group.bar = T, raster = F, angle = 0) + NoLegend() + ggtitle("Cell clusters"))
dev.off()

####################################
### Create summary plots outfile
####################################
writeLines("\n*** Create summary plots outfile dataset 1 ***\n")

if (regexpr("^y$", SummaryPlots, ignore.case = T)[1] == 1) {
  suppressPackageStartupMessages(library(staplr))
  SummaryPlotsPdf<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge<-c(paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_QC_VlnPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_Heatmap.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
  )
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("Create summary file")
  }else{
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    print("Create summary file including t-SNE's for selected genes")
  }
  
  ### Stappling *pdf's
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge, output_filepath = SummaryPlotsPdf)
}





####################################
####################################
####################################
### Process dataset 2
####################################
####################################
####################################
writeLines("\n*** Process dataset 2 ***\n")

####################################
### QC EDA violin plots dataset 2
####################################
writeLines("\n*** QC EDA violin plots dataset 2 ***\n")

### Get unfiltered data QC statistics
nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object_d2.u@meta.data$nFeature_RNA, nGenes = 1)
nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object_d2.u@meta.data$nCount_RNA,   nCount_RNA = 1)
percent.mito.u.df  <-data.frame(Expression_level = seurat.object_d2.u@meta.data$percent.mito, percent.mito = 1)
percent.ribo.u.df  <-data.frame(Expression_level = seurat.object_d2.u@meta.data$percent.ribo, percent.ribo = 1)
#
nFeature_RNAStats.u<-paste(c(" mean = ",round(mean(seurat.object_d2.u@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object_d2.u@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.u  <-paste(c( "mean = ",round(mean(seurat.object_d2.u@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object_d2.u@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.u     <-paste(c(" mean = ",round(mean(seurat.object_d2.u@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object_d2.u@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
percent.ribo.u     <-paste(c(" mean = ",round(mean(seurat.object_d2.u@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object_d2.u@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")

### Get filtered data QC statistics
nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object_d2.f@meta.data$nFeature_RNA, nGenes = 2)
nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object_d2.f@meta.data$nCount_RNA,   nCount_RNA = 2)
percent.mito.f.df  <-data.frame(Expression_level = seurat.object_d2.f@meta.data$percent.mito, percent.mito = 2)
percent.ribo.f.df  <-data.frame(Expression_level = seurat.object_d2.f@meta.data$percent.ribo, percent.ribo = 2)
#
nFeature_RNAStats.f<-paste(c(" mean = ",round(mean(seurat.object_d2.f@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object_d2.f@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.f  <-paste(c(" mean = ",round(mean(seurat.object_d2.f@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object_d2.f@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.f     <-paste(c(" mean = ",round(mean(seurat.object_d2.f@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object_d2.f@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
percent.ribo.f     <-paste(c(" mean = ",round(mean(seurat.object_d2.f@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object_d2.f@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")

### Put QC statistics together
nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
percent.mito.m.df  <-data.frame(rbind(percent.mito.u.df,percent.mito.f.df))
percent.ribo.m.df  <-data.frame(rbind(percent.ribo.u.df,percent.ribo.f.df))
LabelUnfiltered    <-paste("Before filters: No. of cells = ", nrow(seurat.object_d2.u@meta.data), sep ="", collapse = "")
LabelFiltered      <-paste("After filters:  No. of cells = ", nrow(seurat.object_d2.f@meta.data), sep ="", collapse = "")
#LabelFilters       <-paste("Filters: No. of genes per cell = ", ListNGenes, "Fraction mitochondrial protein genes per cell = ", ListPMito), sep ="", collapse = "")

### Commands for violin ggplot's
DataForHeader.df<-data.frame(forx = c(0.4,0.4), fory = c(0.09,0.03), label = c(LabelFiltered,LabelUnfiltered))
Headers.plot<-ggplot(data=DataForHeader.df, aes(x = forx, y = fory)) + theme_void() + 
  geom_point(colour = ColoursQCViolinPlots, size = 7) + xlim(0,1) + ylim(0,0.12) +
  geom_text(hjust = -0.08, label = c(LabelUnfiltered,LabelFiltered), size = 4)

nFeature_RNA.plot<-ggplot(data=nFeature_RNA.m.df, aes(x = factor(nGenes), y = Expression_level)) +
  geom_violin(aes(fill = factor(nGenes))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="No. of genes") +
  annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])

nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
  geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="No. of reads") +
  annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])

percent.mito.plot<-ggplot(data=percent.mito.m.df, aes(x = factor(percent.mito), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.mito))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Mitochondrial genes (fraction)") +
  annotate("text", x = 1 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.f, col = ColoursQCViolinPlots[[2]])

percent.ribo.plot<-ggplot(data=percent.ribo.m.df, aes(x = factor(percent.ribo), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.ribo))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Ribosomal protein genes (fraction)") +
  annotate("text", x = 1 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.f, col = ColoursQCViolinPlots[[2]])

bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, percent.mito.plot, percent.ribo.plot, ncol = 4)

### Create a *pdf file with the violin ggplot's

VlnPlotPdf<-paste(Tempdir, "/" , PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_QC_VlnPlot.pdf", sep="")
pdf(file=VlnPlotPdf, width = 12, height = 7)
print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
dev.off()

####################################
### Feature-vs-feature scatter plots
####################################
writeLines("\n*** Feature-vs-feature scatter plots dataset 2 ***\n")

UnfilteredData.df<-data.frame(nCount_RNA = seurat.object_d2.u@meta.data$nCount_RNA,
                              nGene = seurat.object_d2.u@meta.data$nFeature_RNA,
                              percent.mito = seurat.object_d2.u@meta.data$percent.mito,
                              filtered_out = colnames(seurat.object_d2.u) %in% colnames(seurat.object_d2.f))
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = ColoursQCViolinPlots[[1]])
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = FALSE, replacement = ColoursQCViolinPlots[[2]])
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = 4)
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = FALSE, replacement = 16)

FeatureVsFeaturePlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_NumbReadsVsNumbGenesAndMito_VlnPlot.pdf", sep="")
pdf(file=FeatureVsFeaturePlotPdf, width = 10, height = 5)
par(mfrow=c(1,2))
## No. of reads vs. Mitochond. %
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$percent.mito, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "Mitochond. %")
legend("topright", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))

## No. of reads vs. No. of genes
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$nGene, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "No. of genes")
legend("topleft", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))
dev.off()

####################################
### Remove the Unfiltered seurat object
####################################
writeLines("\n*** Remove the Unfiltered seurat object dataset 2 ***\n")

rm(seurat.object_d2.u)
rm(UnfilteredData.df)

####################################
### Normalize data
####################################
writeLines("\n*** Normalize data dataset 2 ***\n")

seurat.object_d2.f <- NormalizeData(object = seurat.object_d2.f, normalization.method = "LogNormalize", scale.factor = DefaultParameters$ScaleFactor)

####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes dataset 2 ***\n")

### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object_d2.f <- FindVariableFeatures(object = seurat.object_d2.f, selection.method = 'mean.var.plot', mean.cutoff = c(DefaultParameters$XLowCutoff, DefaultParameters$XHighCutoff), dispersion.cutoff = c(DefaultParameters$YCutoff, Inf))
VariableGenes<-VariableFeatures(object = seurat.object_d2.f)
length(VariableGenes)

write(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_VariableGenes.txt", sep=""), x=VariableGenes)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_VariableGenes.pdf", sep=""))
print(VariableFeaturePlot(object = seurat.object_d2.f, cols = c("blue", "red")))
dev.off()

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation dataset 2 ***\n")

seurat.object_d2.f <- ScaleData(object = seurat.object_d2.f, vars.to.regress = c("nCount_RNA", "percent.mito"), features = rownames(x = seurat.object_d2.f), display.progress=F)

####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
### Note: DimPlot can now handle 'umap' and 'tsne' in addition to 'pca', but for 'umap'
### you must first install the umap-learn python package (e.g. via pip install umap-learn)
### https://github.com/satijalab/seurat/blob/master/R/dimensional_reduction.R
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA dataset 2 ***\n")

seurat.object_d2.f <- RunPCA(object = seurat.object_d2.f, features = VariableGenes, verbose = T, do.print = T, ndims.print = DefaultParameters$PrintPCA.PcsPrint, nfeatures.print = DefaultParameters$PrintPCA.GenesPrint)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
print(VizDimLoadings(object = seurat.object_d2.f, reduction = "pca", dims = DefaultParameters$VizPCA.PcsUse, nfeatures = DefaultParameters$VizPCA.nGenesToPlot))
dev.off()

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_PCAPlot.pdf", sep=""))
print(DimPlot(object = seurat.object_d2.f, dims = c(1,2), reduction = "pca") + theme(legend.position="none"))
dev.off()

seurat.object_d2.f <- ProjectDim(object = seurat.object_d2.f, overwrite = T, verbose = T, nfeatures.print = 10)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_PCHeatmap.C1to",DefaultParameters$PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
print(DimHeatmap(object = seurat.object_d2.f, dims = 1:DefaultParameters$PCHeatmapComponentsToPlot, cells = DefaultParameters$PCHeatmapCellsUse, balanced = T))
dev.off()

####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components dataset 2 ***\n")

### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_PCElbowPlot.pdf", sep=""))
print(ElbowPlot(object = seurat.object_d2.f))
dev.off()

####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells dataset 2 ***\n")

StartTimeClustering<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object_d2.f <- FindNeighbors(object = seurat.object_d2.f, dims = PcaDimsUse) ## This step was part of FindClusters() in Seurat v2
seurat.object_d2.f <- FindClusters(object = seurat.object_d2.f, reduction.type = "pca", resolution = Resolution_d2, print.output = 0, save.SNN = T)
EndTimeClustering<-Sys.time()

CellNames<-rownames(seurat.object_d2.f@meta.data)
ClusterIdent<-seurat.object_d2.f@active.ident
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution_d2, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)

####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene dataset 2 ***\n")

cluster.averages<-AverageExpression(object = seurat.object_d2.f, use.scale = F, use.counts = F)
OutfileClusterAverages<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("r", Resolution_d2, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
writeLines("\n*** Run Non-linear dimensional reduction (tSNE) dataset 2 ***\n")

### NOTE: if the datasets is small you may get
### "Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large."
### One can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10

seurat.object_d2.f <- RunTSNE(object = seurat.object_d2.f, dims.use = PcaDimsUse, do.fast = T)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot.pdf", sep=""), width = 7, height = 7)
print(DimPlot(object = seurat.object_d2.f, reduction = 'tsne', group.by = 'ident', label = T, label.size=10))
dev.off()

####################################
### Saving the R object
####################################
writeLines("\n*** Saving the R object dataset 2 ***\n")

OutfileRDS_d2<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_object.rds", sep="")
saveRDS(seurat.object_d2.f, file = OutfileRDS_d2)

####################################
### Write out t-SNE coordinates
####################################
writeLines("\n*** Write out t-SNE coordinates dataset 2 ***\n")

OutfileTsneCoordinates<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNECoordinates.tsv", sep="")

Headers<-paste("Barcode",paste(colnames(seurat.object_d2.f@reductions$tsne@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileTsneCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
write.table(seurat.object_d2.f@reductions$tsne@cell.embeddings, file = OutfileTsneCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Colour t-SNE by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo
####################################
writeLines("\n*** Colour t-SNE by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo dataset 2 ***\n")

CellPropertiesToTsne<-c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 21, height = 5)
print(FeaturePlot(object = seurat.object_d2.f, features = CellPropertiesToTsne, cols = c("lightgrey", "blue"), reduction = "tsne", ncol = 4, pt.size = 1.5))
dev.off()


####################################
### Colour t-SNE by -infile_colour_tsne
####################################
writeLines("\n*** Colour t-SNE by -infile_colour_tsne dataset 2 ***\n")

if (regexpr("^NA$", InfileColourTsne, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  seurat.object.meta.data<-seurat.object_d2.f@meta.data
  ExtraCellProperties <- data.frame(read.table(InfileColourTsne, header = T, row.names = 1))
  head(ExtraCellProperties)
  
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_tsne and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object_d2.f <- AddMetaData(object = seurat.object_d2.f, metadata = ExtraCellProperties)
  
  # Generating outfile
  # Note DimPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...)) or  par(mfrow=())
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""), height = DefaultParameters$BaseSizeSingleTnePlot, width = DefaultParameters$BaseSizeSingleTnePlot)
  for (property in colnames(ExtraCellProperties)) {
    print(DimPlot(object = seurat.object_d2.f, reduction = 'tsne', group.by = property, combine = T, legend = "none") + ggtitle(property))
  }
  dev.off()
}

####################################
### Colour t-SNE plots showing each requested gene
####################################
writeLines("\n*** Colour t-SNE plots showing each requested gene dataset 2 ***\n")

### To program layout() for more than 3 genes in multiple rows
if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
  print("No selected genes for t-SNE plots")
}else{
  ListOfGenesForTsnes<-unlist(strsplit(ListGenes, ","))
  if (length(ListOfGenesForTsnes) <= 4) {
    pdfWidth  <- (length(ListOfGenesForTsnes) * DefaultParameters$BaseSizeMultipleWidth)
    pdfHeight <- DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- length(ListOfGenesForTsnes)
  }else{
    pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
    pdfHeight <- (as.integer(length(ListOfGenesForTsnes) / 4) + 1) * DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- 4
  }
  
  pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  print(FeaturePlot(object = seurat.object_d2.f, ncol = nColFeaturePlot, features = c(ListOfGenesForTsnes), cols = c("lightgrey", "blue"), reduction = "tsne"))
  dev.off()
}

####################################
### Finding differentially expressed genes for each cell cluster
####################################
writeLines("\n*** Finding differentially expressed genes for each cell cluster dataset 2 ***\n")

### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
#########
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object_d2.f, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object_d2.f, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
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

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object_d2.f@meta.data))

StartTimeFindAllMarkers<-Sys.time()
seurat.object.markers_d2 <- FindAllMarkers(object = seurat.object_d2.f, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
EndTimeFindAllMarkers<-Sys.time()
#
OufileTsvMarkers_d2<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_MarkersPerCluster.tsv",sep="")
OufileCsvMarkers_d2<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_MarkersPerCluster.csv",sep="")
write.table(data.frame("GENE"=rownames(seurat.object.markers_d2),seurat.object.markers_d2), file = OufileTsvMarkers_d2, row.names = F,sep="\t",quote = F)
write.csv(seurat.object.markers_d2, file = OufileCsvMarkers_d2)

### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne<-(seurat.object.markers_d2 %>% group_by(cluster) %>% top_n(DefaultParameters$NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers_d2[["cluster"]]))

####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes dataset 2 ***\n")

NumberOfPanesForFeaturesPlot<-(NumberOfClusters*DefaultParameters$NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5

### Here need to program a better way to control the y-axis labels
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object_d2.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), slot = 'counts', log = T, adjust = 1, pt.size = 0.5))
dev.off()

####################################
### t-SNE plots showing each cluster top genes
####################################
writeLines("\n*** t-SNE plots showing each cluster top genes dataset 2 ***\n")

pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
pdfHeight <- NumberOfClusters * DefaultParameters$BaseSizeMultipleHeight / 2
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(FeaturePlot(object = seurat.object_d2.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), cols = c("lightgrey", "blue"), reduction = "tsne"))
dev.off()

####################################
### Cell clusters heatmap
####################################
writeLines("\n*** Cell clusters heatmap dataset 2 ***\n")

top_genes_by_cluster_for_heatmap <- seurat.object.markers_d2 %>% group_by(cluster) %>% top_n(n = DefaultParameters$NumberOfGenesPerClusterToPlotHeatmap, wt = avg_logFC)
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_Heatmap.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(DoHeatmap(object = seurat.object_d2.f, features = top_genes_by_cluster_for_heatmap$gene, label = T, group.bar = T, raster = F, angle = 0) + NoLegend() + ggtitle("Cell clusters"))
dev.off()

####################################
### Create summary plots outfile
####################################
writeLines("\n*** Create summary plots outfile dataset 2 ***\n")

if (regexpr("^y$", SummaryPlots, ignore.case = T)[1] == 1) {
  suppressPackageStartupMessages(library(staplr))
  SummaryPlotsPdf<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge<-c(paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_QC_VlnPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_Heatmap.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
  )
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("Create summary file")
  }else{
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    print("Create summary file including t-SNE's for selected genes")
  }
  
  ### Stappling *pdf's
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge, output_filepath = SummaryPlotsPdf)
}





####################################
####################################
####################################
### Process merged dataset
####################################
####################################
####################################
writeLines("\n*** Process merged dataset ***\n")

####################################
### QC EDA violin plots merged dataset
####################################
writeLines("\n*** QC EDA violin plots merged dataset ***\n")

### Get unfiltered data QC statistics
nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object_c12.u@meta.data$nFeature_RNA, nGenes = 1)
nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object_c12.u@meta.data$nCount_RNA,   nCount_RNA = 1)
percent.mito.u.df  <-data.frame(Expression_level = seurat.object_c12.u@meta.data$percent.mito, percent.mito = 1)
percent.ribo.u.df  <-data.frame(Expression_level = seurat.object_c12.u@meta.data$percent.ribo, percent.ribo = 1)
#
nFeature_RNAStats.u<-paste(c(" mean = ",round(mean(seurat.object_c12.u@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object_c12.u@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.u  <-paste(c( "mean = ",round(mean(seurat.object_c12.u@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object_c12.u@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.u     <-paste(c(" mean = ",round(mean(seurat.object_c12.u@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object_c12.u@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
percent.ribo.u     <-paste(c(" mean = ",round(mean(seurat.object_c12.u@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object_c12.u@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")

### Get filtered data QC statistics
nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object_c12.f@meta.data$nFeature_RNA, nGenes = 2)
nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object_c12.f@meta.data$nCount_RNA,   nCount_RNA = 2)
percent.mito.f.df  <-data.frame(Expression_level = seurat.object_c12.f@meta.data$percent.mito, percent.mito = 2)
percent.ribo.f.df  <-data.frame(Expression_level = seurat.object_c12.f@meta.data$percent.ribo, percent.ribo = 2)
#
nFeature_RNAStats.f<-paste(c(" mean = ",round(mean(seurat.object_c12.f@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object_c12.f@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.f  <-paste(c(" mean = ",round(mean(seurat.object_c12.f@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object_c12.f@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.f     <-paste(c(" mean = ",round(mean(seurat.object_c12.f@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object_c12.f@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
percent.ribo.f     <-paste(c(" mean = ",round(mean(seurat.object_c12.f@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object_c12.f@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")

### Put QC statistics together
nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
percent.mito.m.df  <-data.frame(rbind(percent.mito.u.df,percent.mito.f.df))
percent.ribo.m.df  <-data.frame(rbind(percent.ribo.u.df,percent.ribo.f.df))
LabelUnfiltered    <-paste("Before filters: No. of cells = ", nrow(seurat.object_c12.u@meta.data), sep ="", collapse = "")
LabelFiltered      <-paste("After filters:  No. of cells = ", nrow(seurat.object_c12.f@meta.data), sep ="", collapse = "")
#LabelFilters       <-paste("Filters: No. of genes per cell = ", ListNGenes, "Fraction mitochondrial protein genes per cell = ", ListPMito), sep ="", collapse = "")

### Commands for violin ggplot's
DataForHeader.df<-data.frame(forx = c(0.4,0.4), fory = c(0.09,0.03), label = c(LabelFiltered,LabelUnfiltered))
Headers.plot<-ggplot(data=DataForHeader.df, aes(x = forx, y = fory)) + theme_void() + 
  geom_point(colour = ColoursQCViolinPlots, size = 7) + xlim(0,1) + ylim(0,0.12) +
  geom_text(hjust = -0.08, label = c(LabelUnfiltered,LabelFiltered), size = 4)

nFeature_RNA.plot<-ggplot(data=nFeature_RNA.m.df, aes(x = factor(nGenes), y = Expression_level)) +
  geom_violin(aes(fill = factor(nGenes))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="No. of genes") +
  annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])

nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
  geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="No. of reads") +
  annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])

percent.mito.plot<-ggplot(data=percent.mito.m.df, aes(x = factor(percent.mito), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.mito))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Mitochondrial genes (fraction)") +
  annotate("text", x = 1 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.f, col = ColoursQCViolinPlots[[2]])

percent.ribo.plot<-ggplot(data=percent.ribo.m.df, aes(x = factor(percent.ribo), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.ribo))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Ribosomal protein genes (fraction)") +
  annotate("text", x = 1 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.f, col = ColoursQCViolinPlots[[2]])

bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, percent.mito.plot, percent.ribo.plot, ncol = 4)

### Create a *pdf file with the violin ggplot's

VlnPlotPdf<-paste(Tempdir, "/" , PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_QC_VlnPlot.pdf", sep="")
pdf(file=VlnPlotPdf, width = 12, height = 7)
print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
dev.off()

####################################
### Feature-vs-feature scatter plots
####################################
writeLines("\n*** Feature-vs-feature scatter plots merged dataset ***\n")

UnfilteredData.df<-data.frame(nCount_RNA = seurat.object_c12.u@meta.data$nCount_RNA,
                              nGene = seurat.object_c12.u@meta.data$nFeature_RNA,
                              percent.mito = seurat.object_c12.u@meta.data$percent.mito,
                              filtered_out = colnames(seurat.object_c12.u) %in% colnames(seurat.object_c12.f))
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = ColoursQCViolinPlots[[1]])
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = FALSE, replacement = ColoursQCViolinPlots[[2]])
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = 4)
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = FALSE, replacement = 16)

FeatureVsFeaturePlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_NumbReadsVsNumbGenesAndMito_VlnPlot.pdf", sep="")
pdf(file=FeatureVsFeaturePlotPdf, width = 10, height = 5)
par(mfrow=c(1,2))
## No. of reads vs. Mitochond. %
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$percent.mito, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "Mitochond. %")
legend("topright", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))

## No. of reads vs. No. of genes
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$nGene, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "No. of genes")
legend("topleft", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))
dev.off()

####################################
### Remove the Unfiltered seurat object
####################################
writeLines("\n*** Remove the Unfiltered seurat object merged dataset ***\n")

rm(seurat.object_c12.u)
rm(UnfilteredData.df)

####################################
### Normalize data
####################################
writeLines("\n*** Normalize data merged dataset ***\n")

seurat.object_c12.f <- NormalizeData(object = seurat.object_c12.f, normalization.method = "LogNormalize", scale.factor = DefaultParameters$ScaleFactor)

####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes merged dataset ***\n")

### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object_c12.f <- FindVariableFeatures(object = seurat.object_c12.f, selection.method = 'mean.var.plot', mean.cutoff = c(DefaultParameters$XLowCutoff, DefaultParameters$XHighCutoff), dispersion.cutoff = c(DefaultParameters$YCutoff, Inf))
VariableGenes<-VariableFeatures(object = seurat.object_c12.f)

write(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_VariableGenes.txt", sep=""), x=VariableGenes)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_VariableGenes.pdf", sep=""))
print(VariableFeaturePlot(object = seurat.object_c12.f, cols = c("blue", "red")))
dev.off()

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation merged dataset ***\n")

seurat.object_c12.f <- ScaleData(object = seurat.object_c12.f, vars.to.regress = c("nCount_RNA", "percent.mito"), features = rownames(x = seurat.object_c12.f), display.progress=F)

####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
### Note: DimPlot can now handle 'umap' and 'tsne' in addition to 'pca', but for 'umap'
### you must first install the umap-learn python package (e.g. via pip install umap-learn)
### https://github.com/satijalab/seurat/blob/master/R/dimensional_reduction.R
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA merged dataset ***\n")

seurat.object_c12.f <- RunPCA(object = seurat.object_c12.f, features = VariableGenes, verbose = T, do.print = T, ndims.print = DefaultParameters$PrintPCA.PcsPrint, nfeatures.print = DefaultParameters$PrintPCA.GenesPrint)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
print(VizDimLoadings(object = seurat.object_c12.f, reduction = "pca", dims = DefaultParameters$VizPCA.PcsUse, nfeatures = DefaultParameters$VizPCA.nGenesToPlot))
dev.off()

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_PCAPlot.pdf", sep=""))
print(DimPlot(object = seurat.object_c12.f, dims = c(1,2), reduction = "pca") + theme(legend.position="none"))
dev.off()

seurat.object_c12.f <- ProjectDim(object = seurat.object_c12.f, overwrite = T, verbose = T, nfeatures.print = 10)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_PCHeatmap.C1to",DefaultParameters$PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
print(DimHeatmap(object = seurat.object_c12.f, dims = 1:DefaultParameters$PCHeatmapComponentsToPlot, cells = DefaultParameters$PCHeatmapCellsUse, balanced = T))
dev.off()

####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components merged dataset ***\n")

### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_PCElbowPlot.pdf", sep=""))
print(ElbowPlot(object = seurat.object_c12.f))
dev.off()

####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells merged dataset ***\n")

StartTimeClustering<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object_c12.f <- FindNeighbors(object = seurat.object_c12.f, dims = PcaDimsUse) ## This step was part of FindClusters() in Seurat v2
seurat.object_c12.f <- FindClusters(object = seurat.object_c12.f, reduction.type = "pca", resolution = Resolution_c12, print.output = 0, save.SNN = T)
EndTimeClustering<-Sys.time()

CellNames<-rownames(seurat.object_c12.f@meta.data)
ClusterIdent<-seurat.object_c12.f@active.ident
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution_c12, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)

####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene merged dataset ***\n")

cluster.averages<-AverageExpression(object = seurat.object_c12.f, use.scale = F, use.counts = F)
OutfileClusterAverages<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("r", Resolution_c12, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
writeLines("\n*** Run Non-linear dimensional reduction (tSNE) merged dataset ***\n")

### NOTE: if the datasets is small you may get
### "Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large."
### One can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10

seurat.object_c12.f <- RunTSNE(object = seurat.object_c12.f, dims.use = PcaDimsUse, do.fast = T)

pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot.pdf", sep=""), width = 7, height = 7)
print(DimPlot(object = seurat.object_c12.f, reduction = 'tsne', group.by = 'ident', label = T, label.size=10))
dev.off()

####################################
### Saving the R object
####################################
writeLines("\n*** Saving the R object merged dataset ***\n")

OutfileRDS_c12<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_object.rds", sep="")
saveRDS(seurat.object_c12.f, file = OutfileRDS_c12)

####################################
### Write out t-SNE coordinates
####################################
writeLines("\n*** Write out t-SNE coordinates merged dataset ***\n")

OutfileTsneCoordinates<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNECoordinates.tsv", sep="")

Headers<-paste("Barcode",paste(colnames(seurat.object_c12.f@reductions$tsne@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileTsneCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
write.table(seurat.object_c12.f@reductions$tsne@cell.embeddings, file = OutfileTsneCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Colour t-SNE by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo
####################################
writeLines("\n*** Colour t-SNE by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo merged dataset ***\n")

CellPropertiesToTsne<-c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 21, height = 5)
print(FeaturePlot(object = seurat.object_c12.f, features = CellPropertiesToTsne, cols = c("lightgrey", "blue"), reduction = "tsne", ncol = 4, pt.size = 1.5))
dev.off()

####################################
### Colour t-SNE by -infile_colour_tsne
####################################
writeLines("\n*** Colour t-SNE by -infile_colour_tsne merged dataset ***\n")

if (regexpr("^NA$", InfileColourTsne, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  seurat.object.meta.data<-seurat.object_c12.f@meta.data
  ExtraCellProperties <- data.frame(read.table(InfileColourTsne, header = T, row.names = 1))
  head(ExtraCellProperties)
  
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_tsne and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object_c12.f <- AddMetaData(object = seurat.object_c12.f, metadata = ExtraCellProperties)
  
  # Generating outfile
  # Note DimPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...)) or  par(mfrow=())
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""), height = DefaultParameters$BaseSizeSingleTnePlot, width = DefaultParameters$BaseSizeSingleTnePlot)
  for (property in colnames(ExtraCellProperties)) {
    print(DimPlot(object = seurat.object_c12.f, reduction = 'tsne', group.by = property, combine = T, legend = "none") + ggtitle(property))
  }
  dev.off()
}

####################################
### Colour t-SNE plots showing each requested gene
####################################
writeLines("\n*** Colour t-SNE plots showing each requested gene merged dataset ***\n")

### To program layout() for more than 3 genes in multiple rows
if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
  print("No selected genes for t-SNE plots")
}else{
  ListOfGenesForTsnes<-unlist(strsplit(ListGenes, ","))
  if (length(ListOfGenesForTsnes) <= 4) {
    pdfWidth  <- (length(ListOfGenesForTsnes) * DefaultParameters$BaseSizeMultipleWidth)
    pdfHeight <- DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- length(ListOfGenesForTsnes)
  }else{
    pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
    pdfHeight <- (as.integer(length(ListOfGenesForTsnes) / 4) + 1) * DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- 4
  }
  
  pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  print(FeaturePlot(object = seurat.object_c12.f, ncol = nColFeaturePlot, features = c(ListOfGenesForTsnes), cols = c("lightgrey", "blue"), reduction = "tsne"))
  dev.off()
}

####################################
### Finding differentially expressed genes for each cell cluster
####################################
writeLines("\n*** Finding differentially expressed genes for each cell cluster merged dataset ***\n")

### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
#########
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object_c12.f, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object_c12.f, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
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

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object_c12.f@meta.data))

StartTimeFindAllMarkers<-Sys.time()
seurat.object.markers_c12 <- FindAllMarkers(object = seurat.object_c12.f, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
EndTimeFindAllMarkers<-Sys.time()
#
OufileTsvMarkers_c12<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_MarkersPerCluster.tsv",sep="")
OufileCsvMarkers_c12<-paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_MarkersPerCluster.csv",sep="")
write.table(data.frame("GENE"=rownames(seurat.object.markers_c12),seurat.object.markers_c12), file = OufileTsvMarkers_c12, row.names = F,sep="\t",quote = F)
write.csv(seurat.object.markers_c12, file = OufileCsvMarkers_c12)

### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne<-(seurat.object.markers_c12 %>% group_by(cluster) %>% top_n(DefaultParameters$NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers_c12[["cluster"]]))

####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes merged dataset ***\n")

NumberOfPanesForFeaturesPlot<-(NumberOfClusters*DefaultParameters$NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5

### Here need to program a better way to control the y-axis labels
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object_c12.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), slot = 'counts', log = T, adjust = 1, pt.size = 0.5))
dev.off()

####################################
### t-SNE plots showing each cluster top genes
####################################
writeLines("\n*** t-SNE plots showing each cluster top genes merged dataset ***\n")

pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
pdfHeight <- NumberOfClusters * DefaultParameters$BaseSizeMultipleHeight / 2
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(FeaturePlot(object = seurat.object_c12.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), cols = c("lightgrey", "blue"), reduction = "tsne"))
dev.off()

####################################
### Cell clusters heatmap
####################################
writeLines("\n*** Cell clusters heatmap merged dataset ***\n")

top_genes_by_cluster_for_heatmap <- seurat.object.markers_c12 %>% group_by(cluster) %>% top_n(n = DefaultParameters$NumberOfGenesPerClusterToPlotHeatmap, wt = avg_logFC)
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5
pdf(file=paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_Heatmap.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(DoHeatmap(object = seurat.object_c12.f, features = top_genes_by_cluster_for_heatmap$gene, label = T, group.bar = T, raster = F, angle = 0) + NoLegend() + ggtitle("Cell clusters"))
dev.off()

####################################
### Create summary plots outfile
####################################
writeLines("\n*** Create summary plots outfile merged dataset ***\n")

if (regexpr("^y$", SummaryPlots, ignore.case = T)[1] == 1) {
  suppressPackageStartupMessages(library(staplr))
  SummaryPlotsPdf<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge<-c(paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_QC_VlnPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_Heatmap.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
  )
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("Create summary file")
  }else{
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_c12, ".SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    print("Create summary file including t-SNE's for selected genes")
  }
  
  ### Stappling *pdf's
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge, output_filepath = SummaryPlotsPdf)
}

####################################
####################################
####################################
##### t-SNE plots coloured by dataset
####################################
####################################
####################################

####################################
### t-SNE plots coloured by dataset
####################################
writeLines("\n*** t-SNE plots coloured by dataset ***\n")
#
### Save a table with barcode_ids and datasets
OutfileCellDatasets<-paste(Tempdir,"/",PrefixOutfiles, "_" , PrefixInput_c12, ".SEURAT_cell_ids.tsv", sep="")
ColNames_d1<-colnames(seurat.object_c12.f@assays$RNA@counts)[grep(pattern = "Dataset1", x = colnames(seurat.object_c12.f@assays$RNA@counts))]
ColNames_d2<-colnames(seurat.object_c12.f@assays$RNA@counts)[grep(pattern = "Dataset2", x = colnames(seurat.object_c12.f@assays$RNA@counts))]
Headers<-paste("CELL_BARCODE", paste("Dataset", sep = "", collapse = "") ,sep="\t")
write.table(Headers,file = OutfileCellDatasets, row.names = F, col.names = F, sep="\t", quote = F)
write.table(x=paste(ColNames_d1, PrefixInput_d1, sep = "\t"), file=OutfileCellDatasets, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(x=paste(ColNames_d2, PrefixInput_d2, sep = "\t"), file=OutfileCellDatasets, row.names = F, col.names = F, sep="\t", quote = F, append = T)
DatasetPerBarcode <- data.frame(read.table(OutfileCellDatasets, header = T, row.names = 1))
seurat.object_c12.f <- AddMetaData(object = seurat.object_c12.f, metadata = DatasetPerBarcode)
#
## t-SNE showing both datasets overlapped, coloured per dataset
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_" , PrefixInput_c12, ".SEURAT_TSNEPlot_ColourByDataset.pdf", sep=""))
DimPlot(object = seurat.object_c12.f, reduction = 'tsne', group.by = "Dataset", plot.title = "Colour by dataset", cols = c(ColorDataset1Tsne,ColorDataset2Tsne), pt.size = 0.5)
dev.off()
#
## t-SNE's showing each dataset, coordinates and colours based on clustering of seurat.object_c12.f
seurat.object_c12.f_d1<-SubsetData(object = seurat.object_c12.f, cells.use = ColNames_d1)
seurat.object_c12.f_d2<-SubsetData(object = seurat.object_c12.f, cells.use = ColNames_d2)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_CoordsAndColoursFromMerged.pdf", sep=""))
DimPlot(object = seurat.object_c12.f_d1, reduction = 'tsne', plot.title = paste(PrefixInput_d1, " - coords and colours from merged", sep = "", collapse = ""), do.label = T,  label.size = 10)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_CoordsAndColoursFromMerged.pdf", sep=""))
DimPlot(object = seurat.object_c12.f_d2, reduction = 'tsne', plot.title = paste(PrefixInput_d2, " - coords and colours from merged", sep = "", collapse = ""), do.label = T,  label.size = 10)
dev.off()
#
## t-SNE's showing each dataset, coordinates based on clustering of seurat.object_c12.f, colours based on each dataset clustering
InfileClusters1<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_CellClusters.tsv", sep="")
InfileClusters2<-paste(Tempdir,"/", PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_CellClusters.tsv", sep="")
Clusters1 <- data.frame(read.table(InfileClusters1, header = T, row.names = 1))
Clusters2 <- data.frame(read.table(InfileClusters2, header = T, row.names = 1))
rownames(Clusters1)<-paste("Dataset1",rownames(Clusters1),sep = "-")
rownames(Clusters2)<-paste("Dataset2",rownames(Clusters2),sep = "-")
seurat.object_c12.f_d1 <- AddMetaData(object = seurat.object_c12.f_d1, metadata = Clusters1)
seurat.object_c12.f_d2 <- AddMetaData(object = seurat.object_c12.f_d2, metadata = Clusters2)
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d1, ".SEURAT_TSNEPlot_CoordsFromMerged_ColoursFromItself.pdf", sep=""))
DimPlot(object =  seurat.object_c12.f_d1, reduction = 'tsne', plot.title = paste(PrefixInput_d1, " - coords from merged, colours from", PrefixInput_d1, sep = "", collapse = ""), do.label = T,  label.size = 10)
dev.off()
#
pdf(file=paste(Tempdir,"/",PrefixOutfiles, "_", PrefixInput_d2, ".SEURAT_TSNEPlot_CoordsFromMerged_ColoursFromItself.pdf", sep=""))
DimPlot(object =  seurat.object_c12.f_d2, reduction = 'tsne', plot.title = paste(PrefixInput_d2, " - coords from merged, colours from", PrefixInput_d2, sep = "", collapse = ""), do.label = T,  label.size = 10)
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
marker_file_list <- c(Dataset1 = OufileCsvMarkers_d1, Dataset2 = OufileCsvMarkers_d2)
fList <- c(Dataset1 = OutfileRDS_d1, Dataset2 = OutfileRDS_d2, CombDatasets = OutfileRDS_c12)
objList <- lapply(fList, readRDS)
objList
#
single_obj_list <- c(Dataset1 = objList$Dataset1, Dataset2 = objList$Dataset2)
single_obj_list
#


HERE APPARENTLY NEED TO MODIFY library(ClusterMap)
to be able to use Seurat v3


### edge_cutoff = 0.1 ## default
res <- cluster_map(marker_file_list, edge_cutoff = EdgeCutoff, output = PrefixCombOutfiles, single_obj_list = single_obj_list, comb_obj = objList$CombDatasets)
res ### This command is needed to make the circos plot
#
#################
##### Generate log2 and log10 scatter plot of number of genes
#################
### adds correl
#
rowMeansNoLogDataset_d1<-rowMeans(input.matrix_d1)
rowMeansNoLogDatasetWo0_d1<-rowMeansNoLogDataset_d1
#
rowMeansNoLogDataset_d2<-rowMeans(input.matrix_d2)
rowMeansNoLogDatasetWo0_d2<-rowMeansNoLogDataset_d2
#
ValueToReplace0sInLog<-min(c(rowMeansNoLogDataset_d1[which(rowMeansNoLogDataset_d1 > 0)], rowMeansNoLogDataset_d2[which(rowMeansNoLogDataset_d2 > 0)]))
#
rowMeansNoLogDatasetWo0_d2[rowMeansNoLogDatasetWo0_d2==0]<-ValueToReplace0sInLog
rowMeansLog10DatasetWo0_d2<-log10(rowMeansNoLogDatasetWo0_d2)
rowMeansNoLogDatasetWo0_d1[rowMeansNoLogDatasetWo0_d1==0]<-ValueToReplace0sInLog
rowMeansLog10DatasetWo0_d1<-log10(rowMeansNoLogDatasetWo0_d1)
#
r.cor.nolog <- round(cor(x=rowMeansNoLogDataset_d1, y=rowMeansNoLogDataset_d2,  method = "pearson",  use = "pairwise.complete.obs"), digits = 3)
s.cor.nolog <- round(cor(x=rowMeansNoLogDataset_d1, y=rowMeansNoLogDataset_d2,  method = "spearman", use = "pairwise.complete.obs"), digits = 3)
r.cor.log10 <- round(cor(x=rowMeansLog10DatasetWo0_d1, y=rowMeansLog10DatasetWo0_d2,  method = "pearson",  use = "pairwise.complete.obs"), digits = 3)
s.cor.log10 <- round(cor(x=rowMeansLog10DatasetWo0_d1, y=rowMeansLog10DatasetWo0_d2,  method = "spearman", use = "pairwise.complete.obs"), digits = 3)
#
##### log10 scatter plot of number of genes
pdf(file=paste(PrefixCombOutfiles, ".reads_per_gene_scatter.log10.pdf", sep=""))
plot(x=rowMeansLog10DatasetWo0_d1, y=rowMeansLog10DatasetWo0_d2,
     pch=16,col=rgb(0.35,0.70,0.90,0.25),
     xlim=c(min(rowMeansLog10DatasetWo0_d1,rowMeansLog10DatasetWo0_d2), max(rowMeansLog10DatasetWo0_d1,rowMeansLog10DatasetWo0_d2)), ylim=c(min(rowMeansLog10DatasetWo0_d1,rowMeansLog10DatasetWo0_d2), max(rowMeansLog10DatasetWo0_d1,rowMeansLog10DatasetWo0_d2)),
     xlab=PrefixInput_d1,
     ylab=PrefixInput_d2
)
#
### highlight mito.features
points(x = rowMeansLog10DatasetWo0_d1[mito.features], y=rowMeansLog10DatasetWo0_d2[mito.features], col="red")
MaxOfPlot    <-max(rowMeansLog10DatasetWo0_d1,rowMeansLog10DatasetWo0_d2)
MinOfPlot    <-min(rowMeansLog10DatasetWo0_d1,rowMeansLog10DatasetWo0_d2)
MiddleOfPlot <-mean(c(MinOfPlot,MaxOfPlot))
MitoXposition<-MiddleOfPlot + ((MaxOfPlot - MiddleOfPlot) * 0.7)
MitoYposition<-MiddleOfPlot + ((MiddleOfPlot - MinOfPlot) * 0.75)
text(x = MitoXposition, y= MitoYposition, labels = "mito.features", col = "red")
text(x = MinOfPlot+0.5, y= MaxOfPlot-0.5, pos=3, labels = paste("r=",r.cor.log10,sep = "", collapse = ""))
text(x = MinOfPlot+0.5, y= MaxOfPlot-0.5, pos=1,labels = paste("s=",s.cor.log10,sep = "", collapse = ""))
dev.off()
#
#################
##### Generate table log10(mean gene expression)
#################
label_av_d1    <- paste(PrefixInput_d1, "_av", collapse = "", sep ="")
label_log10_d1 <- paste(PrefixInput_d1, "_log10", collapse = "", sep ="")
label_av_d2    <- paste(PrefixInput_d2, "_av", collapse = "", sep ="")
label_log10_d2 <- paste(PrefixInput_d2, "_log10", collapse = "", sep ="")
#
mat<- data.frame("Gene" = rownames(input.matrix_d1),
                 V1 = rowMeans(input.matrix_d1),
                 V2 = log10(rowMeans(input.matrix_d1)),
                 V3 = rowMeans(input.matrix_d2),
                 V4 = log10(rowMeans(input.matrix_d2))
)
colnames(mat) <- c("Gene", label_av_d1, label_log10_d1, label_av_d2, label_log10_d2)
write.table(file = paste(PrefixCombOutfiles, "_raw_counts.tsv", sep = ""), x=mat[rev(order(mat[1])),], row.names = F, sep = "\t", quote = F)
#
#################
##### Generate overlapping histograms
#################
h1 <- hist(rowMeansLog10DatasetWo0_d1,breaks=30, plot = F)
h2 <- hist(rowMeansLog10DatasetWo0_d2,breaks=30, plot = F)
pdf(file=paste(PrefixCombOutfiles, ".hist.pdf", sep=""))
plot( h1, col=rgb(0,0,1,1/4), border= rgb(0,0,1,1/4), xlim=c(MinOfPlot,MaxOfPlot), main = "",xlab = "log10(mean gene expression)")
plot( h2, col=rgb(1,0,0,1/4), border= rgb(1,0,0,1/4), xlim=c(MinOfPlot,MaxOfPlot), add=T)
legend("topright",legend = c(PrefixInput_d1, PrefixInput_d2), pch=15, col=c(rgb(0,0,1,1/3),rgb(1,0,0,1/3)), bty ="n", cex=1.3)
dev.off()


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
TookTimeClustering_d1       <-format(difftime(EndTimeClustering_d1,     StartTimeClustering_d1,     units = "min"))
TookTimeFindAllMarkers_d1   <-format(difftime(EndTimeFindAllMarkers_d1, StartTimeFindAllMarkers_d1, units = "min"))
#
TookTimeClustering_d2       <-format(difftime(EndTimeClustering_d2,     StartTimeClustering_d2,     units = "min"))
TookTimeFindAllMarkers_d2   <-format(difftime(EndTimeFindAllMarkers_d2, StartTimeFindAllMarkers_d2, units = "min"))
#
TookTimeClustering_c12     <-format(difftime(EndTimeClustering_c12,     StartTimeClustering_c12,     units = "min"))
TookTimeFindAllMarkers_c12 <-format(difftime(EndTimeFindAllMarkers_c12, StartTimeFindAllMarkers_c12, units = "min"))
#
TookTimeOverall            <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "min"))
#
OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CPUusage.txt", sep="")
ReportTime<-c(
  paste("clustering_dataset_d1",TookTimeClustering_d1,collapse = "\t"),
  paste("FindAllMarkers_dataset_d1",TookTimeFindAllMarkers_d1,collapse = "\t"),
  paste("clustering_dataset_d2",TookTimeClustering_d2,collapse = "\t"),
  paste("FindAllMarkers_dataset_d2",TookTimeFindAllMarkers_d2,collapse = "\t"),
  paste("clustering_dataset_c12",TookTimeClustering_c12,collapse = "\t"),
  paste("FindAllMarkers_dataset_c12",TookTimeFindAllMarkers_c12,collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)
#
write(file = OutfileCPUusage, x=c(ReportTime))
#
####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"/SEURAT/",sep="",collapse = ""))
#
outfiles_to_move <- list.files(Tempdir,pattern = PrefixOutfiles, full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/SEURAT/",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})
#
####################################
### Turning warnings on
####################################
options(warn = oldw)
#
####################################
### Finish
####################################
print("END - All done!!! Took time:")
print(ReportTime)

quit()
