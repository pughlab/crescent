####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
### 
### Requires Seurat v3. Main differences from versions using Seurat v2 include:
### 1) Updated function names
### 2) Reads MTX inputs from Cell Ranger v3 (incl. features.tsv.gz instead of genes.tsv)
### 3) A new function DimPlot() is used for 'uma', 'tsne' and 'pca' plots, instead of PCAPlot() and TSNEPlot()
###    Hence, now we define 'dimensions' instead of principal components
### 4) Default parameters are indexed using list(DefaultParameters) instead of variable names. This allows list(DefaultParameters)
###    to be reported in a *log file
### 5) Implemented QC violin plots directly in ggplots instead of Seurat's VlnPlot() to have more control of layout, legends, titles, etc.
### 6) Using R base plot() to create feature-vs-feature scatter plots instead of FeatureScatter()
### 7) Implemented all ggplots and Seurat plotting functions (which are based on ggplots) with print() function, like:
###    `print(FeaturePlot(...))` instead of `FeaturePlot(...)` alone
###    Otherwise using ggplots and Seurat plots inside if/else loops cause errors
###    https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-do-lattice_002ftrellis-graphics-not-work_003f
### 8) Violin plots and t-SNE QC plots show now percentage of ribosomal protein genes in cells
### 9) Implemented library(future) for parallelization
### 10) Added a stopwatch to every step
### 11) Added a step to define outdirs and CWL parameters, which can be tuned with option `-w`, variable 'RunsCwl'
### 12) Added a step to transform select *pdf files into *png
### 
### THINGS TO DO:
### 1) Parallelize FindAllMarkers() to send each cluster to a separate core
###    Although it seems that FindMarkers() is already parallelized by Seurat:
###    https://satijalab.org/seurat/v3.0/future_vignette.html
### 3) Pick the right number of dimension components
###    E.g. try to automatically get the inflection point from the PCElbowPlot() function output
### 4) Pick the right resolution from FindClusters()
###    E.g. implement Iness and Bader paper https://f1000research.com/articles/7-1522/v1 approach
###    By picking the number of clusters based on differentially expressed genes
### 5) Add a lists of ENSEMBL Ids for mitochondrial genes instead of just MT- and mt- (at gene names)
###    Need to do it for both Human and Mouse
### 6) In Seurat v2, Suluxan reported that the -c example Javier provided called a duplicated row names error
###    Need to see if it's still happening
### 7) To implement a flag to see if inputted dataset matches expected percent.mito
###    In particular, if imputting a whole cell sample and using `-m 0,0.05` likely will filterout most cells
###    and the script may crash at step 'Perform linear dimensional reduction by PCA'
###    because the filtered matrix will be too sparse and small to get the gene PC's
### 8) To compare SCTransform() vs. [NormalizeData(), ScaleData(), and FindVariableFeatures()]
###    Do they produce similar results? Maybe use Circos plot to evaluate
### 9) Since we are creating dimension reduction plots (UMAP and tSNE) in a for() loop,
###    need to record list(StopWatchStart) and list(StopWatchEnd) for these steps
###
### THINGS NICE TO HAVE:
### 1) Assigning cell type identity to clusters (needs supervised annotations, maybe based on GSVA)
### 2) Use knitr() to produce html plot layouts (https://yihui.name/knitr/demo/stitch/)
### 3) In "Load data" we use Seurat's Read10X() function. In this command, when all barcodes come from the same sample (i.e. finish with the same digit), like:
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
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
### Requires Seurat v3 (tested on v3.0.3.9023), which can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # To run parallel processes
### library(staplr)     only if using option '-s y', note it needs pdftk
####################################

####################################
### Required external packages
####################################
### 'pdftk'   to merge selected *pdf files into *summary_plots.pdf using libary(staplr)
###           and to add statistics to violin plots overlapping pdf files
###           in Mac you can install it with something like:
###           'sudo port install pdftk'
###           https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/
### 'convert' to transofrm *pdf files into *png
###           it's part of ImageMagic (https://imagemagick.org/index.php)
####################################

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
              help="Either the path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files;
                or path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns
                Notes:
                The 'MTX' files can be for example the output from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--input_type"), default="NA",
              help="Indicates if input is either a 'MTX' directory or a 'DGE' file
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-b", "--normalize_and_scale_sample"), default="Y",
              help="Indicates if input is either raw counts and hence needs to be normalized and scaled [use 'y']
                or it's already normalized (e.g. if inputting transcripts per million counts (TPMs)) [use 'n']
                Default = 'Y'"),
  #
  make_option(c("-r", "--resolution"), default="1",
              help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of cell clusters
                Default = '1'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-s", "--summary_plots"), default="y",
              help="Indicates if a *summary_plots.pdf file should be generated [use 'y'] or not [use 'n']
                Note this needs 'pdftk' and R library(staplr)
                Default = 'y'"),
  #
  make_option(c("-c", "--infile_colour_dim_red_plots"), default="NA",
              help="A <tab> delimited table of barcodes and discrete properties to colour the t-SNE, like:
                Barcode              CellClass    InOtherDatasets
                AAACCTGAGCGGCTTC-1   1            yes
                AAACCTGAGTCGAGTG-1   1            no
                AAACCTGCAAAGGAAG-1   2            yes
                AAACCTGGTCTCATCC-1   2            no
                Default = 'NA' (i.e. no --infile_colour_dim_red_plots is provided)"),
  #
  make_option(c("-g", "--list_genes"), default="NA",
              help="A <comma> delimited list of gene identifiers whose expression will be mapped into the t-SNE plots
                Default = 'NA' (no --list_genes provided)"),
  #
  make_option(c("-a", "--opacity"), default="0.2",
              help="If using a --list_genes, this parameter provides a value for the minimal opacity of gene expression. Use a value between 0 and 1
                Default = '0.2'"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Number of PCA dimensions to use for clustering and dimension reduction functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file *PCElbowPlot.pdf,
                use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                Default = '10'"),
  #
  make_option(c("-m", "--percent_mito"), default="0,0.05",
              help="<comma> delimited min,max fraction of gene counts of mitochondrial originin a cell to be included in normalization and clustering analyses
                For example, for whole cell scRNA-seq use '0,0.2', or for Nuc-seq use '0,0.05'
                For negative values (e.g. if using TPM in log scale refer negative values with an 'n', like this 'n1,0.5')
                Default = '0,0.05'"),
  #
  make_option(c("-n", "--n_genes"), default="50,8000",
              help="<comma> delimited min,max number of unique genes measured in a cell to be included in normalization and clustering analyses
                Default = '50,8000'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh
                Default = '0.01'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system
                Default = 'MAX'"),
  #
  make_option(c("-w", "--run_cwl"), default="N",
              help="Indicates if this script is running inside a virtual machine container, such that outfiles are written directly into the 'HOME' . Type 'y/Y' or 'n/N'.
                Note, if using 'y/Y' this supersedes option -o
                Default = 'N'")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input                   <- opt$input
InputType               <- opt$input_type
NormalizeAnsScale       <- opt$normalize_and_scale_sample
Resolution              <- as.numeric(opt$resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
SummaryPlots            <- opt$summary_plots
InfileColourDimRedPlots <- opt$infile_colour_dim_red_plots
ListGenes               <- opt$list_genes
Opacity                 <- as.numeric(opt$opacity)
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
ListPMito               <- opt$percent_mito
ListNGenes              <- opt$n_genes
ThreshReturn            <- as.numeric(opt$return_threshold)
NumbCores               <- opt$number_cores
RunsCwl                 <- opt$run_cwl

### To be able to take negative values for ListPMito entered using optparse (e.g. -m m1,1 means mito.percentage from -1 to 1)
ListPMito <- gsub("m","-",ListPMito, ignore.case = T)

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

ProgramOutdir <- "SEURAT"

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Outfiles will be stored into `ProgramOutdir` directory
  #PrefixOutfiles <- "cwl_run" 
  PrefixOutfiles  <- opt$prefix_outfiles
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) ## Note Tempdir will be the final out-directory as well
}else{
  PrefixOutfiles <- c(paste(PrefixOutfiles,"_res",Resolution,sep=""))
  ## Using `Tempdir` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  Tempdir        <- "~/temp" 
  #
  CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
  UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
  #
  Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
  Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
  Outdir<-gsub("/$", "", Outdir)
  Tempdir<-gsub("/$", "", Tempdir)
  #
  dir.create(file.path(Outdir, ProgramOutdir), recursive = T)
  dir.create(file.path(Tempdir), showWarnings = F, recursive = T)
}

####################################
### Define number of cores and RAM for parallelization
####################################

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (is.numeric(NumbCores) == T) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript Runs_Seurat_v3.R -h\n\n", sep=""))
}

cat("Using ", NumbCoresToUse, "cores")

plan(strategy = "multicore", workers = NumbCoresToUse)

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4Gb of RAM (4000*1024^2), use:
options(future.globals.maxSize = 4000 * 1024^2)

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
  ScaleFactor = 1000000, ### Using 1000000 to set scale.factor as counts per million (CPM)
  NormalizationMethod = "LogNormalize",
  
  ### Parameters for Seurat variable gene detection
  VarGeneDetectSelectMethod = "mean.var.plot",
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
  
  ### Parameters for dimmension reduction plots
  BaseSizeSinglePlotPdf  = 7,
  BaseSizeSinglePlotPng  = 480,
  BaseSizeMultiplePlotPdfWidth  = 3.7,
  BaseSizeMultiplePlotPdfHeight = 3
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

### Dimension reduction methods
DimensionReductionMethods<-list()
DimensionReductionMethods$umap$name <-"UMAP"
DimensionReductionMethods$tsne$name <-"TSNE"
DimensionReductionMethods$umap$run  <-as.function(RunUMAP)
DimensionReductionMethods$tsne$run  <-as.function(RunTSNE)

####################################
### Start stopwatches
####################################

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("input", "input_type", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

StopWatchStart$LoadScRNAseqData  <- Sys.time()

if (regexpr("^MTX$", InputType, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  input.matrix <- Read10X(data.dir = Input)
}else if (regexpr("^DGE$", InputType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  ## Note `check.names = F` is needed for both `fread` and `data.frame`
  input.matrix <- as.matrix(data.frame(fread(Input, check.names = F), row.names=1, check.names = F))
}else{
  stop(paste("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}
dim(input.matrix)

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Create a Seurat object
####################################
writeLines("\n*** Create a Seurat object ***\n")

StopWatchStart$CreateSeuratObject  <- Sys.time()

seurat.object.u  <- CreateSeuratObject(counts = input.matrix, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = PrefixOutfiles)
nCellsInOriginalMatrix<-length(seurat.object.u@meta.data$orig.ident)

StopWatchEnd$CreateSeuratObject  <- Sys.time()

####################################
### Get mitochondrial genes
####################################
writeLines("\n*** Get  mitochondrial genes ***\n")

StopWatchStart$GetMitoGenes  <- Sys.time()

mitoRegExpressions<- paste(c("^MT-"),collapse = "|")
mito.features <- grep(pattern = mitoRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)

if (length(mito.features)[[1]] > 0) {
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['percent.mito']] <- percent.mito
}else{
  percent.mito <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['percent.mito']] <- percent.mito
}

StopWatchEnd$GetMitoGenes  <- Sys.time()

####################################
### Get ribosomal protein genes
####################################
writeLines("\n*** Get ribosomal protein genes ***\n")

StopWatchStart$GetRiboGenes  <- Sys.time()

riboRegExpressions<- paste(c("^MRPL", "^MRPS", "^RPL", "^RPS"),collapse = "|")
ribo.features <- grep(pattern = riboRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)

if (length(ribo.features)[[1]] > 0) {
  percent.ribo <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['percent.ribo']] <- percent.ribo
}else{
  percent.ribo <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['percent.ribo']] <- percent.ribo
}

StopWatchEnd$GetRiboGenes  <- Sys.time()

####################################
### Filter cells based gene counts and mitochondrial representation
####################################
writeLines("\n*** Filter cells based gene counts and mitochondrial representation ***\n")

StopWatchStart$FilterCells  <- Sys.time()

rm(seurat.object.f)

if (length(mito.features)[[1]] > 0) {
  seurat.object.f<-subset(x = seurat.object.u, subset = nFeature_RNA >= DefaultParameters$MinGenes & nFeature_RNA <= DefaultParameters$MaxGenes & percent.mito >= DefaultParameters$MinPMito & percent.mito <= DefaultParameters$MaxPMito)
  
  ### Get list of barcodes excluded by mito and/or nFeature_RNA
  BarcodesExcludedByMito            <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = percent.mito >= DefaultParameters$MinPMito & percent.mito <= DefaultParameters$MaxPMito)))
  BarcodesExcludedByNFeature        <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA >= DefaultParameters$MinGenes & nFeature_RNA <= DefaultParameters$MaxGenes)))
  BarcodesExcludedByMitoAndNFeature <- intersect(BarcodesExcludedByMito, BarcodesExcludedByNFeature)

}else{
  seurat.object.f<-subset(x = seurat.object.u, subset = nFeature_RNA >= DefaultParameters$MinGenes & nFeature_RNA <= DefaultParameters$MaxGenes)

  ### Get list of barcodes excluded by nFeature_RNA
  BarcodesExcludedByMito            <- NULL
  BarcodesExcludedByNFeature        <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA >= DefaultParameters$MinGenes & nFeature_RNA <= DefaultParameters$MaxGenes)))
  BarcodesExcludedByMitoAndNFeature <- NULL

}
NumberOfBarcodesExcludedByMito            <- length(BarcodesExcludedByMito)
NumberOfBarcodesExcludedByNFeature        <- length(BarcodesExcludedByNFeature)
NumberOfBarcodesExcludedByMitoAndNFeature <- length(BarcodesExcludedByMitoAndNFeature)

StopWatchEnd$FilterCells  <- Sys.time()

### Just reporting the summary of the UNfiltered and filtered objects
seurat.object.u
seurat.object.f

####################################
### QC EDA violin plots
####################################
writeLines("\n*** QC EDA violin plots ***\n")

StopWatchStart$QCviolinplots  <- Sys.time()

QCStats <- list()
QCStats$unfiltered$mean$nFeature_RNA<-round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$unfiltered$median$nFeature_RNA<-round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$unfiltered$mean$nCount_RNA<-round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$unfiltered$median$nCount_RNA<-round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$unfiltered$mean$percent.mito<-round(mean(seurat.object.u@meta.data[,"percent.mito"]),3)
QCStats$unfiltered$median$percent.mito<-round(median(seurat.object.u@meta.data[,"percent.mito"]),3)
QCStats$unfiltered$mean$percent.ribo<-round(mean(seurat.object.u@meta.data[,"percent.ribo"]),3)
QCStats$unfiltered$median$percent.ribo<-round(median(seurat.object.u@meta.data[,"percent.ribo"]),3)
#
QCStats$filtered$mean$nFeature_RNA<-round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$filtered$median$nFeature_RNA<-round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$filtered$mean$nCount_RNA<-round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$filtered$median$nCount_RNA<-round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$filtered$mean$percent.mito<-round(mean(seurat.object.u@meta.data[,"percent.mito"]),3)
QCStats$filtered$median$percent.mito<-round(median(seurat.object.u@meta.data[,"percent.mito"]),3)
QCStats$filtered$mean$percent.ribo<-round(mean(seurat.object.u@meta.data[,"percent.ribo"]),3)
QCStats$filtered$median$percent.ribo<-round(median(seurat.object.u@meta.data[,"percent.ribo"]),3)

### Get unfiltered data QC statistics
nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
percent.mito.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$percent.mito, percent.mito = 1)
percent.ribo.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$percent.ribo, percent.ribo = 1)
#
nFeature_RNAStats.u<-paste(c(" mean = ", QCStats$unfiltered$mean$nFeature_RNA,"\n", "median = ", QCStats$unfiltered$median$nFeature_RNA), sep = "", collapse="")
nCount_RNAStats.u  <-paste(c(" mean = ", QCStats$unfiltered$mean$nCount_RNA,  "\n", "median = ", QCStats$unfiltered$median$nCount_RNA),   sep = "", collapse="")
percent.mito.u     <-paste(c(" mean = ", QCStats$unfiltered$mean$percent.mito,"\n", "median = ", QCStats$unfiltered$median$percent.mito), sep = "", collapse="")
percent.ribo.u     <-paste(c(" mean = ", QCStats$unfiltered$mean$percent.ribo,"\n", "median = ", QCStats$unfiltered$median$percent.ribo), sep = "", collapse="")

### Get filtered data QC statistics
nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$nFeature_RNA, nGenes = 2)
nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object.f@meta.data$nCount_RNA,   nCount_RNA = 2)
percent.mito.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$percent.mito, percent.mito = 2)
percent.ribo.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$percent.ribo, percent.ribo = 2)
#
nFeature_RNAStats.f<-paste(c(" mean = ", QCStats$filtered$mean$nFeature_RNA,"\n", "median = ", QCStats$filtered$median$nFeature_RNA), sep = "", collapse="")
nCount_RNAStats.f  <-paste(c(" mean = ", QCStats$filtered$mean$nCount_RNA,  "\n", "median = ", QCStats$filtered$median$nCount_RNA),   sep = "", collapse="")
percent.mito.f     <-paste(c(" mean = ", QCStats$filtered$mean$percent.mito,"\n", "median = ", QCStats$filtered$median$percent.mito), sep = "", collapse="")
percent.ribo.f     <-paste(c(" mean = ", QCStats$filtered$mean$percent.ribo,"\n", "median = ", QCStats$filtered$median$percent.ribo), sep = "", collapse="")

### Put QC statistics together
nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
percent.mito.m.df  <-data.frame(rbind(percent.mito.u.df,percent.mito.f.df))
percent.ribo.m.df  <-data.frame(rbind(percent.ribo.u.df,percent.ribo.f.df))
LabelUnfiltered    <-paste("Before filters: No. of cells = ", nrow(seurat.object.u@meta.data), sep ="", collapse = "")
LabelFiltered      <-paste("After filters:  No. of cells = ", nrow(seurat.object.f@meta.data), sep ="", collapse = "")
NumberOfCells<- list()
NumberOfCells[["unfiltered"]] <- nrow(seurat.object.u@meta.data)
NumberOfCells[["filtered"]]   <- nrow(seurat.object.f@meta.data)

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
  # labs(x="No. of genes") +
  labs(x=paste("No. of genes", "\n", "Filter: ", "min=", ListNGenes[[1]], " max=", ListNGenes[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNFeature, sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])

nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
  geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x=paste("No. of reads", "\n", "no filter was used", "\n", sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])

percent.mito.plot<-ggplot(data=percent.mito.m.df, aes(x = factor(percent.mito), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.mito))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Mitochondrial genes (fraction)") +
  labs(x=paste("Mitochondrial genes", "\n", "Filter: ", "min=", ListPMito[[1]], " max=", ListPMito[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByMito, sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.f, col = ColoursQCViolinPlots[[2]])

percent.ribo.plot<-ggplot(data=percent.ribo.m.df, aes(x = factor(percent.ribo), y = Expression_level)) +
  geom_violin(aes(fill = factor(percent.ribo))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x=paste("Ribosomal protein genes (fraction)", "\n", "no filter was used", "\n", sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.f, col = ColoursQCViolinPlots[[2]])

bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, percent.mito.plot, percent.ribo.plot, ncol = 4)

### Create a *pdf file with the violin ggplot's
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_QC_VlnPlot.pdf", sep="")
pdf(file=VlnPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 1.7, height = DefaultParameters$BaseSizeSinglePlotPdf)
print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
dev.off()

StopWatchEnd$QCviolinplots  <- Sys.time()

####################################
### QC summary table
####################################
writeLines("\n*** QC summary table ***\n")

StopWatchStart$QCsummarytable  <- Sys.time()

Headers<-paste("Run_prefix","Filter_status", "Number_of_cells", "Metric", "Measurement", "Value", sep = "\t", collapse = "")
write.table(Headers,paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_QC_summary.tsv",sep=""), row.names = F, sep="", quote = F, col.names = F)

for (filter_status in names(QCStats)) {
  for (metric in names(QCStats[[filter_status]])) {
    for (measurement in names(QCStats[[filter_status]][[metric]])) {
      write.table(paste(PrefixOutfiles, filter_status, NumberOfCells[[filter_status]], metric, measurement, QCStats[[filter_status]][[metric]][[measurement]],  sep = "\t", collapse = ""), paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_QC_summary.tsv",sep=""), row.names = F, sep="\t", quote = F, col.names = F, append = T)
    }
  }
}

StopWatchEnd$QCsummarytable  <- Sys.time()


####################################
### Feature-vs-feature scatter plot
####################################
writeLines("\n*** Feature-vs-feature scatter plot ***\n")

StopWatchStart$FeatureVsFeatureplot  <- Sys.time()

UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                              nGene = seurat.object.u@meta.data$nFeature_RNA,
                              percent.mito = seurat.object.u@meta.data$percent.mito,
                              filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.f))
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = ColoursQCViolinPlots[[1]])
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = FALSE, replacement = ColoursQCViolinPlots[[2]])
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = 4)
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = FALSE, replacement = 16)

FeatureVsFeaturePlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_NumbReadsVsNumbGenesAndMito_VlnPlot.pdf", sep="")
pdf(file=FeatureVsFeaturePlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 2, height = DefaultParameters$BaseSizeSinglePlotPdf)
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

StopWatchEnd$FeatureVsFeatureplot  <- Sys.time()

####################################
### Remove the Unfiltered seurat object
####################################
writeLines("\n*** Remove the Unfiltered seurat object ***\n")

rm(seurat.object.u)
rm(UnfilteredData.df)

####################################
### Normalize data (if applicable)
####################################
if (regexpr("^Y$", NormalizeAnsScale, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Normalize data ***\n")
  
  StopWatchStart$NormalizeData  <- Sys.time()
  
  seurat.object.f <- NormalizeData(object = seurat.object.f, normalization.method = NormalizationMethod, scale.factor = DefaultParameters$ScaleFactor)
  
  StopWatchEnd$NormalizeData  <- Sys.time()
  
}else{
  writeLines("\n*** Data assumed to be already normalized ***\n")
}

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  
  # output the normalized count matrix for violin plots
  writeLines("\n*** Outputting normalized count matrix as tsv ***\n")
  
  normalized_count_matrix <- as.matrix(seurat.object.f@assays[["RNA"]]@data)
  write.table(normalized_count_matrix, file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_normalized_count_matrix.tsv", sep=""), sep="\t", row.names=TRUE, col.names=TRUE)
} else {
  writeLines("\n*** Skipping normalized count matrix tsv output ***\n")
}

####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes ***\n")

StopWatchStart$FindVariableGenes  <- Sys.time()

### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object.f <- FindVariableFeatures(object = seurat.object.f, selection.method = VarGeneDetectSelectMethod, mean.cutoff = c(DefaultParameters$XLowCutoff, DefaultParameters$XHighCutoff), dispersion.cutoff = c(DefaultParameters$YCutoff, Inf))
VariableGenes<-VariableFeatures(object = seurat.object.f)
length(VariableGenes)

write(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_VariableGenes.txt", sep=""), x=VariableGenes)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_VariableGenes.pdf", sep=""))
print(VariableFeaturePlot(object = seurat.object.f, cols = c("blue", "red")))
dev.off()

StopWatchEnd$FindVariableGenes  <- Sys.time()

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation ***\n")

StopWatchStart$ScaleData  <- Sys.time()

seurat.object.f <- ScaleData(object = seurat.object.f, vars.to.regress = c("nCount_RNA", "percent.mito"), features = rownames(x = seurat.object.f), display.progress=F)

StopWatchEnd$ScaleData  <- Sys.time()

####################################
### Perform linear dimensional reduction by PCA
####################################
### Examine and visualize PCA results a few different ways
### Note: DimPlot can now handle 'umap' and 'tsne' in addition to 'pca', but for 'umap'
### you must first install the umap-learn python package (e.g. via pip install umap-learn)
### https://github.com/satijalab/seurat/blob/master/R/dimensional_reduction.R
###
### Note on Error:
### `Error in irlba(A = t(x = object), nv = npcs, ...) :
### max(nu, nv) must be strictly less than min(nrow(A), ncol(A))
### Calls: RunPCA ... RunPCA -> RunPCA.Assay -> RunPCA -> RunPCA.default -> irlba
### Execution halted`
### This error has to do with `-m` parameter
### For whole cell scRNA-seq use '0,0.2', or for Nuc-seq use '0,0.05'
### If you use '0,0.05' for whole cell you will filter out most of the cells and won't have enough data to get the requested -d (PC's)
### Also see here: https://github.com/satijalab/seurat/issues/127
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA ***\n")

StopWatchStart$RunPCA  <- Sys.time()

seurat.object.f <- RunPCA(object = seurat.object.f, features = VariableGenes, verbose = T, do.print = T, ndims.print = DefaultParameters$PrintPCA.PcsPrint, nfeatures.print = DefaultParameters$PrintPCA.GenesPrint)

StopWatchEnd$RunPCA  <- Sys.time()

StopWatchStart$PCAPlots  <- Sys.time()

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_VizPCA.pdf", sep=""), width=7, height=10)
print(VizDimLoadings(object = seurat.object.f, reduction = "pca", dims = DefaultParameters$VizPCA.PcsUse, nfeatures = DefaultParameters$VizPCA.nGenesToPlot))
dev.off()

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_PCAPlot.pdf", sep=""))
print(DimPlot(object = seurat.object.f, dims = c(1,2), reduction = "pca") + theme(legend.position="none"))
dev.off()

seurat.object.f <- ProjectDim(object = seurat.object.f, overwrite = T, verbose = T, nfeatures.print = 10)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_PCHeatmap.C1to",DefaultParameters$PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
print(DimHeatmap(object = seurat.object.f, dims = 1:DefaultParameters$PCHeatmapComponentsToPlot, cells = DefaultParameters$PCHeatmapCellsUse, balanced = T))
dev.off()

StopWatchEnd$PCAPlots  <- Sys.time()

####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components ***\n")

StopWatchStart$GetSignificantPCs  <- Sys.time()

### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time

ForElbowPlot<-ElbowPlot(object = seurat.object.f, ndims = 50, reduction = "pca")
MaxYAxis<-as.integer(max(ForElbowPlot$data$stdev)+1)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_PCElbowPlot.pdf", sep=""))
print(ForElbowPlot
      + scale_x_continuous(breaks =  seq(from = 0, to = 50, by=5))
      + geom_vline(xintercept = seq(from = 0, to = 50, by=5), linetype='dotted', col="red")
      + scale_y_continuous(breaks =  seq(from = 0, to = MaxYAxis, by=0.5))
      + geom_hline(yintercept = seq(from = 0, to = MaxYAxis, by=0.5), linetype='dotted', col="red")
)
dev.off()

StopWatchEnd$GetSignificantPCs  <- Sys.time()

####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells ***\n")

StopWatchStart$ClusterCells  <- Sys.time()

?FindNeighbors

options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.f <- FindNeighbors(object = seurat.object.f, dims = PcaDimsUse) ## This step was part of FindClusters() in Seurat v2
seurat.object.f <- FindClusters(object = seurat.object.f, reduction.type = "pca", resolution = Resolution, print.output = 0, save.SNN = T)
EndTimeClustering<-Sys.time()

StopWatchEnd$ClusterCells  <- Sys.time()

StopWatchStart$CellClusterTables  <- Sys.time()

CellNames<-rownames(seurat.object.f@meta.data)
ClusterIdent<-seurat.object.f@meta.data$seurat_clusters
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)

StopWatchEnd$CellClusterTables  <- Sys.time()

####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene ***\n")

StopWatchStart$AverageGeneExpression  <- Sys.time()

cluster.averages<-AverageExpression(object = seurat.object.f, use.scale = F, use.counts = F)
OutfileClusterAverages<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("r", Resolution, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$AverageGeneExpression  <- Sys.time()

####################################
### Run and plot dimension reductions
####################################
writeLines("\n*** Run and plot dimension reductions ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {

  ####################################
  ### Run non-linear dimensional reductions
  ####################################
  writeLines(paste("\n*** Run ", DimensionReductionMethods[[dim_red_method]][["name"]], " ***\n", sep = "", collapse = ""))

  StopWatchStart$DimensionReduction$dim_red_method  <- Sys.time()
  
  ### NOTES:
  ### In RunTSNE: if the datasets is small user may get error:
  ### `Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large.`
  ### User can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
  
  seurat.object.f <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.f, dims = PcaDimsUse, do.fast = T)

  StopWatchEnd$DimensionReduction$dim_red_method  <- Sys.time()

  StopWatchStart$DimensionReductionPlot$dim_red_method  <- Sys.time()
    
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot.pdf", sep="", collapse = ""))
  print(DimPlot(object = seurat.object.f, reduction = dim_red_method, group.by = 'ident', label = T, label.size=10))
  dev.off()

  StopWatchEnd$DimensionReductionPlot$dim_red_method  <- Sys.time()
  
  ####################################
  ### Write out coordinates
  ####################################
  writeLines(paste("\n*** Write out ", DimensionReductionMethods[[dim_red_method]][["name"]], " coordinates ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  OutfileCoordinates<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Coordinates.tsv", sep="", collapse = "")
  Headers<-paste("Barcode",paste(colnames(seurat.object.f@reductions[[dim_red_method]]@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
  write.table(Headers,file = OutfileCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(seurat.object.f@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

  StopWatchEnd$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo ***\n", sep = "", collapse = ""))
  
  StopWatchStart$QCDimRedPlots$dim_red_method  <- Sys.time()
  
  CellPropertiesToColour<-c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_QC.pdf", sep=""), width = DefaultParameters$BaseSizeSinglePlotPdf * 1.5, height = DefaultParameters$BaseSizeSinglePlotPdf * 1.5)
  print(FeaturePlot(object = seurat.object.f, label = T, order = T, features = CellPropertiesToColour, cols = c("lightgrey", "blue"), reduction = dim_red_method, ncol = 2, pt.size = 1.5))
  dev.off()

  StopWatchEnd$QCDimRedPlots$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by -infile_colour_dim_red_plots
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by -infile_colour_dim_red_plots ***\n", sep = "", collapse = ""))

  if (regexpr("^NA$", InfileColourDimRedPlots, ignore.case = T)[1] == 1) {
    print("No extra barcode-attributes will be used for dimension reduction plots")
  }else{
    
    StopWatchStart$DimRedPlotsColuredByOptC$dim_red_method  <- Sys.time()
    
    seurat.object.meta.data<-seurat.object.f@meta.data
    ExtraCellProperties <- data.frame(read.table(InfileColourDimRedPlots, header = T, row.names = 1))
    print(head(ExtraCellProperties))
    
    # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
    # so that barcodes from --infile_colour_dim_red_plots and --input can match each other
    rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
    seurat.object.f <- AddMetaData(object = seurat.object.f, metadata = ExtraCellProperties)
    
    # Generating outfile
    # Note DimPlot() takes the entire current device (pdf)
    # even if using layout(matrix(...)) or  par(mfrow=())
    # Thus each property plot is written to a separate page of a single *pdf outfile
    pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ExtraProperties.pdf", sep=""), width = DefaultParameters$BaseSizeSinglePlotPdf, height = DefaultParameters$BaseSizeSinglePlotPdf)
    for (property in colnames(ExtraCellProperties)) {
      print(DimPlot(object = seurat.object.f, reduction = dim_red_method, group.by = property, combine = T, legend = "none") + ggtitle(property))
    }
    dev.off()

    StopWatchEnd$DimRedPlotsColuredByOptC$dim_red_method  <- Sys.time()
    
  }
  
  ####################################
  ### Colour dimension reduction plots showing each requested gene
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot showing each requested gene ***\n", sep = "", collapse = ""))
  
  ### To program layout() for more than 3 genes in multiple rows
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("No selected genes for dimension reduction plots")
  }else{

    StopWatchStart$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()

    ListOfGenesForDimRedPlots<-unlist(strsplit(ListGenes, ","))
    if (length(ListOfGenesForDimRedPlots) <= 4) {
      pdfWidth  <- (length(ListOfGenesForDimRedPlots) * DefaultParameters$BaseSizeMultiplePlotPdfWidth)
      pdfHeight <- DefaultParameters$BaseSizeMultiplePlotPdfHeight
      nColFeaturePlot <- length(ListOfGenesForDimRedPlots)
    }else{
      pdfWidth  <- 4 * DefaultParameters$BaseSizeMultiplePlotPdfWidth
      pdfHeight <- (as.integer(length(ListOfGenesForDimRedPlots) / 4) + 1) * DefaultParameters$BaseSizeMultiplePlotPdfHeight
      nColFeaturePlot <- 4
    }
    
    pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
    print(FeaturePlot(object = seurat.object.f, ncol = nColFeaturePlot, features = c(ListOfGenesForDimRedPlots), cols = c("lightgrey", "blue"), reduction = dim_red_method))
    dev.off()

    StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()

  }
}

####################################
### Finding differentially expressed genes for each cell cluster
####################################
writeLines("\n*** Finding differentially expressed genes for each cell cluster ***\n")

### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object.f, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object.f, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster5.markers, n = 5))
###
### NOTES:
### 1) FindAllMarkers() uses return.thresh = 0.01 as defaults, but FindMarkers() displays all genes passing previous filters.
###    Thus to make the outputs between these two commands identical to each other use return.thresh = 1
###
### 2) Default pseudocount.use=1, which sounds high for a Log level correction
###    An earlier version of this script was using 1e-99, but it was probably too small
###    Now using the inverse of the number of cells in the data.
###    This is sufficiently small as to not compress logGER magnitudes,
###    while keeping comparisons with zero reasonably close to the range of potential logGER values (Innes and Bader, 2018, F1000 Research)

StopWatchStart$FindDiffMarkers  <- Sys.time()

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.f@meta.data))

seurat.object.markers <- FindAllMarkers(object = seurat.object.f, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)

write.table(data.frame("GENE"=rownames(seurat.object.markers),seurat.object.markers),paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_MarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)
### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne<-(seurat.object.markers %>% group_by(cluster) %>% top_n(DefaultParameters$NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers[["cluster"]]))

StopWatchEnd$FindDiffMarkers  <- Sys.time()

####################################
### Saving the R object
####################################
writeLines("\n*** Saving the R object ***\n")

StopWatchStart$SaveRDS  <- Sys.time()

OutfileRDS<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_object.rds", sep="")
saveRDS(seurat.object.f, file = OutfileRDS)

StopWatchEnd$SaveRDS  <- Sys.time()

###################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes ***\n")

StopWatchStart$DiffMarkerViolinPlots  <- Sys.time()

NumberOfPanesForFeaturesPlot<-(NumberOfClusters*DefaultParameters$NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5

### Here need to program a better way to control the y-axis labels
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_VlnPlot_CountsLog10_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), slot = 'counts', log = T, adjust = 1, pt.size = 0.5))
dev.off()

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_VlnPlot_Norm_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), adjust = 1, pt.size = 0.5))
dev.off()

StopWatchEnd$DiffMarkerViolinPlots  <- Sys.time()

####################################
### Plot dimension reductions showing each cluster top genes
####################################
writeLines("\n*** Plot dimension reductions showing each cluster top genes ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {
  
  ####################################
  ### Dimension reduction plots showing each cluster top genes
  ####################################
  writeLines("\n*** Dimension reduction plots showing each cluster top genes ***\n")
  
  pdfWidth  <- 4 * DefaultParameters$BaseSizeMultiplePlotPdfWidth
  pdfHeight <- NumberOfClusters * DefaultParameters$BaseSizeMultiplePlotPdfHeight / 2
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  print(FeaturePlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), cols = c("lightgrey", "blue"), reduction = dim_red_method))
  dev.off()
  
}

####################################
### Cell clusters heatmap
####################################
writeLines("\n*** Cell clusters heatmap ***\n")

StopWatchStart$CellClustersHeatmap  <- Sys.time()

top_genes_by_cluster_for_heatmap <- seurat.object.markers %>% group_by(cluster) %>% top_n(n = DefaultParameters$NumberOfGenesPerClusterToPlotHeatmap, wt = avg_logFC)
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_Heatmap.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(DoHeatmap(object = seurat.object.f, features = top_genes_by_cluster_for_heatmap$gene, label = T, group.bar = T, raster = F, angle = 0) + NoLegend() + ggtitle("Cell clusters"))
dev.off()

StopWatchEnd$CellClustersHeatmap  <- Sys.time()

####################################
### Create summary plots outfile
####################################
writeLines("\n*** Create summary plots outfile ***\n")

StopWatchStart$SummaryPlots  <- Sys.time()

if (regexpr("^y$", SummaryPlots, ignore.case = T)[1] == 1) {
  suppressPackageStartupMessages(library(staplr))
  SummaryPlotsPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge<-c(paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_QC_VlnPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_Heatmap.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_UMAPPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_TSNEPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_VlnPlot_CountsLog10_AfterClusters.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_UMAPPlot_EachTopGene.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
  )
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("Create summary file")
  }else{
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_UMAPPlot_SelectedGenes.pdf", sep="", collapse = ""))
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    print("Create summary file including dimension reduction plots for selected genes")
  }
  
  ### Stappling *pdf's
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge, output_filepath = SummaryPlotsPdf)
}

StopWatchEnd$SummaryPlots  <- Sys.time()

####################################
### Transform select *pdf files into *png
####################################
writeLines("\n*** Transform select *pdf files into *png ***\n")

StopWatchStart$TransformPdfToPng  <- Sys.time()

ListOfPdfFilesToPng<-c(paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_PCElbowPlot.pdf", sep = "", collapse = ""),
                        paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_QC_VlnPlot.pdf",  sep = "", collapse = ""),
                        paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_TSNEPlot_EachTopGene.pdf",  sep = "", collapse = "")
)

sapply(ListOfPdfFilesToPng,FUN=function(eachFile) {
  inpdf  <- eachFile
  outpng <- gsub(".pdf$",".png", x= inpdf, ignore.case = T, perl = T)
  CommandConvert <- paste("convert -density 150 ", inpdf, " -quality 300 ", outpng, sep = "", collapse = "")
  system(command = CommandConvert, input = NULL, wait = T)
})

StopWatchEnd$TransformPdfToPng  <- Sys.time()

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_CPUusage.txt", sep="")
write(file = OutfileCPUusage, x = paste("Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
Headers<-paste("Step", "Time(minutes)", sep="\t")
write.table(Headers,file = OutfileCPUusage, row.names = F, col.names = F, sep="\t", quote = F, append = T)

for (stepToClock in names(StopWatchStart)) {
  TimeStart <- StopWatchStart[[stepToClock]]
  TimeEnd   <- StopWatchEnd[[stepToClock]]
  TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
  ReportTime<-c(paste(stepToClock, TimeDiff, sep = "\t", collapse = ""))
  write(file = OutfileCPUusage, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
}

####################################
### Moving outfiles into outdir or keeping them at tempdir (if using CWL)
####################################

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(paste(Tempdir, sep="", collapse = ""))
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  writeLines(paste(Outdir,"/", ProgramOutdir, "/", sep="", collapse = ""))
  
  outfiles_to_move <- list.files(Tempdir, pattern = paste(PrefixOutfiles, ".", ProgramOutdir, "_", sep=""), full.names = F)
  sapply(outfiles_to_move,FUN=function(eachFile) {
    ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
    file.copy(from=paste(Tempdir,"/",eachFile,sep=""), to=paste(Outdir, "/", ProgramOutdir, "/", eachFile,sep=""), overwrite=T)
    file.remove(paste(Tempdir,"/",eachFile,sep=""))
  })
}

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
writeLines(paste("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report", sep = "", collapse = ""))

quit()
