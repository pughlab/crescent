####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
### and https://satijalab.org/seurat/v3.0/sctransform_vignette.html
### 
### THINGS TO DO:
### 1) Parallelize FindAllMarkers() to send each cluster to a separate core
###    Although it seems that FindMarkers() is already parallelized by Seurat:
###    https://satijalab.org/seurat/v3.0/future_vignette.html
### 2) Determine automatically the inflection point from the ElbowPlot() and use it as `dims` for FindNeighbors(), RunUMAP() and RunTSNE()
### 3) Determine automatically the optimal resolution from FindClusters() following Iness and Bader (F1000Research, 2019)
###    By picking the number of clusters based on differentially expressed genes. See https://f1000research.com/articles/7-1522/v1
### 4) Add a lists of ENSEMBL Ids for mitochondrial genes instead of just MT- and mt- (at gene names)
###    Need to do it for both Human and Mouse
### 5) To implement a flag to see if inputted dataset matches expected mito.fraction
###    In particular, if imputting a whole cell sample and using `-m 0,0.05` likely will filterout most cells
###    and the script may crash at step 'Perform linear dimensional reduction by PCA'
###    because the filtered matrix will be too sparse and small to get the gene PC's
### 6) To compare SCTransform() vs. [NormalizeData(), ScaleData(), and FindVariableFeatures()]
###    Do they produce similar results? Maybe use Circos plot to evaluate
###    Implement variants with option -b
### 7) Make *MarkersPerCluster.tsv and *TopTwoMarkersPerCluster.tsv a single file
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
### COMMENTS ON FINDING DIFFERENTIALLY EXPRESSED GENES
### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object.each_sample, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object.each_sample, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
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
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(Seurat))       # (CRAN) to run QC, differential gene expression and clustering analyses
### Package 'Seurat' version 3 is needed - tested using v3.1.1
### Seurat's latest stable version can be installed like:
### install.packages('Seurat')
###
### Or the development version can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
###
### Or a specific older version (e.g. v3.0.3.9023) can be installed from GitHub like:
### `PathForV3_0_3_9023Libs<-paste(.libPaths(), "/Seurat_V3_0_3_9023", sep = "", collapse = "")`
### `dir.create(file.path(PathForV3_0_3_9023Libs), showWarnings = F, recursive = T)``
### `devtools::install_github("satijalab/seurat", ref = "556e598", lib = PathForV3_0_3_9023Libs)`
### ### To check that v3.0.3.9023 is installed
### `library("Seurat",lib.loc = PathForV3_0_3_9023Libs)``
### `packageVersion("Seurat")`
### ### should return:
### [1] ‘3.0.3.9023’

suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # (CRAN) to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(staplr))       # (CRAN) to merge pdf files. Note it needs pdftk available. If not available use `SummaryPlots <- "N"`
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
  make_option(c("-j", "--inputs_remove_barcodes"), default="NA",
              help="Path/name to a <tab> delimited list of barcodes to be removed from analysis, like:
                dataset1_id  AAACCTGAGCTCCCAG
                dataset2_id  AAACCTGTCACCATAG
                dataset3_id  AAACCTGTCAGCTTAG
                Or type 'NA' to include all barcodes
                Default = 'NA'"),
  #
  make_option(c("-b", "--normalize_and_scale_sample"), default="2",
              help="Indicates one of three options:
                '1' --input contains raw counts (e.g. from cellranger) and this script should use Seurat's NormalizeData() function
                '2' --input contains raw counts (e.g. from cellranger) and this script should use Seurat's SCTransform() function
                '3' --input contains pre-normalized counts (e.g. transcripts per million counts [TPMs]) and the normalization step should be skipped
                Default = '2'"),
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
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Number of PCA dimensions to use for clustering and dimension reduction functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file *PCElbowPlot.pdf,
                use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                Default = '10'"),
  #
  make_option(c("-m", "--percent_mito"), default="0,0.05",
              help="<comma> delimited min,max fraction of gene counts of mitochondrial origin a cell to be included in normalization and clustering analyses
                For example, for whole cell scRNA-seq use '0,0.2', or for Nuc-seq use '0,0.05'
                For negative values (e.g. if using TPM in log scale refer negative values with an 'n', like this 'n1,0.5')
                Default = '0,0.05'"),
  #
  make_option(c("-q", "--percent_ribo"), default="0,0.75",
              help="<comma> delimited min,max fraction of gene counts of ribosomal proteins to be included in normalization and clustering analyses
                For negative values (e.g. if using TPM in log scale refer negative values with an 'n', like this 'n1,0.5')
                Default = '0,0.75'"),
  #
  make_option(c("-n", "--n_genes"), default="50,8000",
              help="<comma> delimited min,max number of unique genes measured in a cell to be included in normalization and clustering analyses
                Default = '50,8000'"),
  #
  make_option(c("-v", "--n_reads"), default="1,80000",
              help="<comma> delimited min,max number of reads measured in a cell to be included in normalization and clustering analyses
                Default = '1,80000'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh
                Default = '0.01'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system
                Default = 'MAX'"),
  #
  make_option(c("-s", "--save_r_object"), default="N",
              help="Indicates if a R object with the data and analyzes from the run should be saved
                Note that this may be time consuming. Type 'y/Y' or 'n/N'
                Default = 'N'"),
  #
  make_option(c("-w", "--run_cwl"), default="N",
              help="Indicates if this script is running inside a virtual machine container, such that outfiles are written directly into the 'HOME' . Type 'y/Y' or 'n/N'.
                Note, if using 'y/Y' this supersedes option -o
                Default = 'N'")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input                   <- opt$input
InputType               <- opt$input_type
InfileRemoveBarcodes    <- opt$inputs_remove_barcodes
NormalizeAndScale       <- as.numeric(opt$normalize_and_scale_sample)
Resolution              <- as.numeric(opt$resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
InfileColourDimRedPlots <- opt$infile_colour_dim_red_plots
ListGenes               <- opt$list_genes
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
ListPMito               <- opt$percent_mito
ListNGenes              <- opt$n_genes
ListPRibo               <- opt$percent_ribo
ListNReads              <- opt$n_reads
ThreshReturn            <- as.numeric(opt$return_threshold)
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl

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
  dir.create(file.path(paste0(Tempdir,"/normalized")), showWarnings = F) ## Note Tempdir will be the final out-directory as well
  dir.create(file.path(paste0(Tempdir,"/coordinates")), showWarnings = F) ## Note Tempdir will be the final out-directory as well
  dir.create(file.path(paste0(Tempdir,"/raw")), showWarnings = F) ## Note Tempdir will be the final out-directory as well
  dir.create(file.path(paste0(Tempdir,"/metadata")), showWarnings = F) ## Note Tempdir will be the final out-directory as well
  dir.create(file.path(paste0(Tempdir,"/markers")), showWarnings = F) ## Note Tempdir will be the final out-directory as well
  dir.create(file.path(paste0(Tempdir,"/qc")), showWarnings = F) ## Note Tempdir will be the final out-directory as well
  
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
  Outdir<-gsub("/+", "/", Outdir, perl = T)
  Tempdir<-gsub("/+", "/", Tempdir, perl = T)
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

### To be able to take negative values for ListPMito entered using optparse (e.g. -m m1,1 means mito.percentage from -1 to 1)
ListPMito <- gsub("m","-",ListPMito, ignore.case = T)
ListPRibo <- gsub("m","-",ListPRibo, ignore.case = T)

ListNGenes = unlist(strsplit(ListNGenes, ","))
MinGenes   = as.numeric(ListNGenes[1])
MaxGenes   = as.numeric(ListNGenes[2])
#
ListPMito  = unlist(strsplit(ListPMito,  ","))
MinPMito   = as.numeric(ListPMito[1])
MaxPMito   = as.numeric(ListPMito[2])
#
ListNReads = unlist(strsplit(ListNReads, ","))
MinReads   = as.numeric(ListNReads[1])
MaxReads   = as.numeric(ListNReads[2])
#
ListPRibo  = unlist(strsplit(ListPRibo,  ","))
MinPRibo   = as.numeric(ListPRibo[1])
MaxPRibo   = as.numeric(ListPRibo[2])

DefaultParameters <- list(
  
  ### Parameters for QC plots
  CellPropertiesToQC = c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction"),

  ### Parameters for Seurat filters
  MinCells = 3,
  MinReads = MinReads,
  MaxGenes = MaxGenes,
  MinGenes = MinGenes,
  MaxReads = MaxReads,
  MinPMito = MinPMito,
  MaxPMito = MaxPMito,
  MinPRibo = MinPRibo,
  MaxPRibo = MaxPRibo,
  
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
  MinNumberOfCellsToReducePerplexity = 150,
  ReducedPerplexity = 7,
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

seurat.object.u   <- CreateSeuratObject(counts = input.matrix, project = PrefixOutfiles)

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
  mito.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['mito.fraction']] <- mito.fraction
}else{
  mito.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['mito.fraction']] <- mito.fraction
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
  ribo.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['ribo.fraction']] <- ribo.fraction
}else{
  ribo.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['ribo.fraction']] <- ribo.fraction
}

StopWatchEnd$GetRiboGenes  <- Sys.time()

####################################
### Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation
####################################
writeLines("\n*** Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation ***\n")

StopWatchStart$FilterCells  <- Sys.time()

if (length(mito.features)[[1]] > 0) {
  seurat.object.f<-subset(x = seurat.object.u, subset = 
                            nFeature_RNA >= DefaultParameters$MinGenes
                          & nFeature_RNA <= DefaultParameters$MaxGenes 
                          & nCount_RNA   >= DefaultParameters$MinReads
                          & nCount_RNA   <= DefaultParameters$MaxReads 
                          & mito.fraction >= DefaultParameters$MinPMito
                          & mito.fraction <= DefaultParameters$MaxPMito
                          & ribo.fraction >= DefaultParameters$MinPRibo
                          & ribo.fraction <= DefaultParameters$MaxPRibo)

  ### Get list of barcodes excluded by nFeature_RNA, nCount_RNA, mito.fraction or ribo.fraction
  BarcodesExcludedByNFeature <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA >= DefaultParameters$MinGenes & nFeature_RNA <= DefaultParameters$MaxGenes)))
  BarcodesExcludedByNReads   <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nCount_RNA   >= DefaultParameters$MinReads & nCount_RNA   <= DefaultParameters$MaxReads)))
  BarcodesExcludedByMito     <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = mito.fraction >= DefaultParameters$MinPMito & mito.fraction <= DefaultParameters$MaxPMito)))
  BarcodesExcludedByRibo     <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = ribo.fraction >= DefaultParameters$MinPRibo & ribo.fraction <= DefaultParameters$MaxPRibo)))

}else{
  seurat.object.f<-subset(x = seurat.object.u, subset = 
                            nFeature_RNA >= DefaultParameters$MinGenes
                          & nFeature_RNA <= DefaultParameters$MaxGenes 
                          & nCount_RNA   >= DefaultParameters$MinReads
                          & nCount_RNA   <= DefaultParameters$MaxReads 
                          & ribo.fraction >= DefaultParameters$MinPRibo
                          & ribo.fraction <= DefaultParameters$MaxPRibo)
  
  ### Get list of barcodes excluded by nFeature_RNA, nCount_RNA or ribo.fraction
  BarcodesExcludedByNFeature <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA >= DefaultParameters$MinGenes & nFeature_RNA <= DefaultParameters$MaxGenes)))
  BarcodesExcludedByNReads   <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nCount_RNA   >= DefaultParameters$MinReads & nCount_RNA   <= DefaultParameters$MaxReads)))
  BarcodesExcludedByMito     <- 0
  BarcodesExcludedByRibo     <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = ribo.fraction >= DefaultParameters$MinPRibo & ribo.fraction <= DefaultParameters$MaxPRibo)))
}
NumberOfBarcodesExcludedByNFeature <- length(BarcodesExcludedByNFeature)
NumberOfBarcodesExcludedByNReads   <- length(BarcodesExcludedByNReads)
NumberOfBarcodesExcludedByMito     <- length(BarcodesExcludedByMito)
NumberOfBarcodesExcludedByRibo     <- length(BarcodesExcludedByRibo)

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
QCStats$unfiltered$mean$mito.fraction<-round(mean(seurat.object.u@meta.data[,"mito.fraction"]),3)
QCStats$unfiltered$median$mito.fraction<-round(median(seurat.object.u@meta.data[,"mito.fraction"]),3)
QCStats$unfiltered$mean$ribo.fraction<-round(mean(seurat.object.u@meta.data[,"ribo.fraction"]),3)
QCStats$unfiltered$median$ribo.fraction<-round(median(seurat.object.u@meta.data[,"ribo.fraction"]),3)
#
QCStats$filtered$mean$nFeature_RNA<-round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$filtered$median$nFeature_RNA<-round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$filtered$mean$nCount_RNA<-round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$filtered$median$nCount_RNA<-round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$filtered$mean$mito.fraction<-round(mean(seurat.object.u@meta.data[,"mito.fraction"]),3)
QCStats$filtered$median$mito.fraction<-round(median(seurat.object.u@meta.data[,"mito.fraction"]),3)
QCStats$filtered$mean$ribo.fraction<-round(mean(seurat.object.u@meta.data[,"ribo.fraction"]),3)
QCStats$filtered$median$ribo.fraction<-round(median(seurat.object.u@meta.data[,"ribo.fraction"]),3)

### Get unfiltered data QC statistics
nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
mito.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$mito.fraction, mito.fraction = 1)
ribo.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$ribo.fraction, ribo.fraction = 1)
#
nFeature_RNAStats.u<-paste(c(" mean = ", QCStats$unfiltered$mean$nFeature_RNA,"\n", "median = ", QCStats$unfiltered$median$nFeature_RNA), sep = "", collapse="")
nCount_RNAStats.u  <-paste(c(" mean = ", QCStats$unfiltered$mean$nCount_RNA,  "\n", "median = ", QCStats$unfiltered$median$nCount_RNA),   sep = "", collapse="")
mito.fraction.u     <-paste(c(" mean = ", QCStats$unfiltered$mean$mito.fraction,"\n", "median = ", QCStats$unfiltered$median$mito.fraction), sep = "", collapse="")
ribo.fraction.u     <-paste(c(" mean = ", QCStats$unfiltered$mean$ribo.fraction,"\n", "median = ", QCStats$unfiltered$median$ribo.fraction), sep = "", collapse="")

### Get filtered data QC statistics
nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$nFeature_RNA, nGenes = 2)
nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object.f@meta.data$nCount_RNA,   nCount_RNA = 2)
mito.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$mito.fraction, mito.fraction = 2)
ribo.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$ribo.fraction, ribo.fraction = 2)
#
nFeature_RNAStats.f<-paste(c(" mean = ", QCStats$filtered$mean$nFeature_RNA,"\n", "median = ", QCStats$filtered$median$nFeature_RNA), sep = "", collapse="")
nCount_RNAStats.f  <-paste(c(" mean = ", QCStats$filtered$mean$nCount_RNA,  "\n", "median = ", QCStats$filtered$median$nCount_RNA),   sep = "", collapse="")
mito.fraction.f     <-paste(c(" mean = ", QCStats$filtered$mean$mito.fraction,"\n", "median = ", QCStats$filtered$median$mito.fraction), sep = "", collapse="")
ribo.fraction.f     <-paste(c(" mean = ", QCStats$filtered$mean$ribo.fraction,"\n", "median = ", QCStats$filtered$median$ribo.fraction), sep = "", collapse="")

### Put QC statistics together
nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
mito.fraction.m.df  <-data.frame(rbind(mito.fraction.u.df,mito.fraction.f.df))
ribo.fraction.m.df  <-data.frame(rbind(ribo.fraction.u.df,ribo.fraction.f.df))
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
  labs(x=paste("No. of genes", "\n", "Filter: ", "min=", ListNGenes[[1]], " max=", ListNGenes[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNFeature, sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])

nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
  geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x=paste("No. of reads", "\n", "Filter: ", "min=", ListNReads[[1]], " max=", ListNReads[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNReads, sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])

mito.fraction.plot<-ggplot(data=mito.fraction.m.df, aes(x = factor(mito.fraction), y = Expression_level)) +
  geom_violin(aes(fill = factor(mito.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Mitochondrial genes (fraction)") +
  labs(x=paste("Mitochondrial genes", "\n", "Filter: ", "min=", ListPMito[[1]], " max=", ListPMito[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByMito, sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.f, col = ColoursQCViolinPlots[[2]])

ribo.fraction.plot<-ggplot(data=ribo.fraction.m.df, aes(x = factor(ribo.fraction), y = Expression_level)) +
  geom_violin(aes(fill = factor(ribo.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x=paste("Ribosomal protein genes (fraction)", "\n", "Filter: ", "min=", ListPRibo[[1]], " max=", ListPRibo[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByRibo, sep = "", collapse = "")) +
  annotate("text", x = 1 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.f, col = ColoursQCViolinPlots[[2]])

bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, mito.fraction.plot, ribo.fraction.plot, ncol = 4)

### Create a *pdf file with the violin ggplot's
VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_QC_VlnPlot.pdf", sep="")
pdf(file=VlnPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 1.7, height = DefaultParameters$BaseSizeSinglePlotPdf)
print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
dev.off()

StopWatchEnd$QCviolinplots  <- Sys.time()

####################################
### CWL interactive qc plots
####################################
writeLines("\n*** CWL interactive qc plots ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  # unfiltered
  interactive_qc_plot_u  <-data.frame(Barcodes = row.names(seurat.object.u@meta.data), Number_of_Genes = seurat.object.u@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.u@meta.data$nCount_RNA, Mitochondrial_Genes_Fraction = seurat.object.u@meta.data$percent.mito, Ribosomal_Protein_Genes_Fraction = seurat.object.u@meta.data$percent.ribo)
  colnames(interactive_qc_plot_u) <- c("Barcodes","Number of Genes","Number of Reads","Mitochondrial Genes Fraction","Ribosomal Protein Genes Fraction")
  write.table(interactive_qc_plot_u, paste(Tempdir,"/","qc/","BeforeFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
  
  # filtered
  interactive_qc_plot_f  <-data.frame(Barcodes = row.names(seurat.object.f@meta.data), Number_of_Genes = seurat.object.f@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.f@meta.data$nCount_RNA, Mitochondrial_Genes_Fraction = seurat.object.f@meta.data$percent.mito, Ribosomal_Protein_Genes_Fraction = seurat.object.f@meta.data$percent.ribo )
  colnames(interactive_qc_plot_f) <- c("Barcodes","Number of Genes","Number of Reads","Mitochondrial Genes Fraction","Ribosomal Protein Genes Fraction")
  write.table(interactive_qc_plot_f, paste(Tempdir,"/","qc/","AfterFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
} 

####################################
### Feature-vs-feature scatter plot
####################################
writeLines("\n*** Feature-vs-feature scatter plot ***\n")

StopWatchStart$FeatureVsFeatureplot  <- Sys.time()

UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                              nGene = seurat.object.u@meta.data$nFeature_RNA,
                              mito.fraction = seurat.object.u@meta.data$mito.fraction,
                              filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.f))
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = ColoursQCViolinPlots[[1]])
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = F, replacement = ColoursQCViolinPlots[[2]])
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = 4)
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = F, replacement = 16)

FeatureVsFeaturePlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_NumbReadsVsNumbGenesAndMito_ScatterPlot.pdf", sep="")
pdf(file=FeatureVsFeaturePlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 2, height = DefaultParameters$BaseSizeSinglePlotPdf)
par(mfrow=c(1,2))
## No. of reads vs. Mitochond. %
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$mito.fraction, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "Mitochond. %")
legend("topright", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))

## No. of reads vs. No. of genes
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$nGene, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "No. of genes")
legend("topleft", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))
dev.off()

StopWatchEnd$FeatureVsFeatureplot  <- Sys.time()

####################################
### Write out filter details and number of filtered cells
####################################
writeLines(paste("\n*** Write out filter details and number of filtered cells ***\n", sep = "", collapse = ""))

StopWatchStart$OutTablesFilterDetailsAndFilteredCells  <- Sys.time()

OutTableFilterDetails<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_FilterDetails.tsv", sep="", collapse = "")
Headers<-paste("Step", "Filter_min", "Filter_max", "Mean_before_filter", "Median_before_filter", "Mean_after_filter", "Median_after_filter", "Excluded_cells", sep = "\t", collapse = "")
write.table(Headers, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F)

FilterDetails.nFeature_RNA  <- paste("nFeature_RNA", ListNGenes[[1]], ListNGenes[[2]],
                                     QCStats$unfiltered$mean$nFeature_RNA,  QCStats$unfiltered$median$nFeature_RNA,  QCStats$filtered$mean$nFeature_RNA,  QCStats$filtered$median$nFeature_RNA, 
                                     NumberOfBarcodesExcludedByNFeature, sep = "\t", collapse = "")
FilterDetails.nCount_RNA    <- paste("nCount_RNA", ListNReads[[1]], ListNReads[[2]], 
                                     QCStats$unfiltered$mean$nCount_RNA,    QCStats$unfiltered$median$nCount_RNA,    QCStats$filtered$mean$nCount_RNA,    QCStats$filtered$median$nCount_RNA, 
                                     NumberOfBarcodesExcludedByNReads, sep = "\t", collapse = "")
FilterDetails.mito.fraction <- paste("mito.fraction", ListPMito[[1]], ListPMito[[2]],
                                     QCStats$unfiltered$mean$mito.fraction, QCStats$unfiltered$median$mito.fraction, QCStats$filtered$mean$mito.fraction, QCStats$filtered$median$mito.fraction, 
                                     NumberOfBarcodesExcludedByMito, sep = "\t", collapse = "")
FilterDetails.ribo.fraction <- paste("ribo.fraction", ListPRibo[[1]], ListPRibo[[2]],
                                     QCStats$unfiltered$mean$ribo.fraction, QCStats$unfiltered$median$ribo.fraction, QCStats$filtered$mean$ribo.fraction, QCStats$filtered$median$ribo.fraction, 
                                     NumberOfBarcodesExcludedByRibo, sep = "\t", collapse = "")

write.table(FilterDetails.nFeature_RNA,  file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(FilterDetails.nCount_RNA,    file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(FilterDetails.mito.fraction, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(FilterDetails.ribo.fraction, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)

OutTableFilteredCells<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_NumberOfFilteredCells.tsv", sep="", collapse = "")
write.table(paste("Number_of_cells_before_filters", NumberOfCells[["unfiltered"]], sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste("Number_of_cells_after_filters", NumberOfCells[["filtered"]],    sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$OutTablesFilterDetailsAndFilteredCells  <- Sys.time()

####################################
### Write out QC data
####################################
writeLines("\n*** Write out QC data ***\n")

StopWatchStart$WriteOutQCData  <- Sys.time()

Headers<-paste("Cell_barcode", paste(DefaultParameters$CellPropertiesToQC, sep = "", collapse = "\t") ,sep="\t")

BarcodeIdsBeforeFilters <- unlist(x = colnames(seurat.object.u))
OutfileQCMetadataBeforeFilters<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_", "Before_filters_QC_metadata.tsv", sep = "", collapse = "")
write.table(Headers, file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(BarcodeIdsBeforeFilters, seurat.object.u@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)

BarcodeIdsAfterFilters <- unlist(x = colnames(seurat.object.f))
OutfileQCMetadataAfterFilters<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_", "After_filters_QC_metadata.tsv", sep = "", collapse = "")
write.table(Headers, file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(BarcodeIdsAfterFilters, seurat.object.f@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$WriteOutQCData  <- Sys.time()

####################################
### Remove barcodes by parameter -j (if applicable)
####################################

if (regexpr("^NA$", InfileRemoveBarcodes , ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Ignoring option -j ***\n")

}else{
  
  StopWatchStart$RemoveBarcodes <- Sys.time()
  
  writeLines("\n*** Removing barcodes by parameter -j ***\n")
  
  InfileAllBarcodesToRemove<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), InfileRemoveBarcodes)
  AllBarcodesToRemove.tab<-read.table(InfileAllBarcodesToRemove, header = F, row.names = NULL, stringsAsFactors = FALSE)
  colnames(AllBarcodesToRemove.tab) <- c("Barcode")
  
  seurat.object.full    <- seurat.object.f
  seurat.object.subset  <- subset(seurat.object.full, cells = colnames(seurat.object.full[,AllBarcodesToRemove.tab[,"Barcode"]]))
  print(paste(paste("Before:", ncol(seurat.object.full), sep = "", collapse =""), paste("After:", ncol(seurat.object.subset), sep = "", collapse =""), sep = "  ", collapse = "\n"))
  seurat.object.f       <- seurat.object.subset

  StopWatchEnd$RemoveBarcodes <- Sys.time()
  
}

####################################
### Remove the temporal Seurat objects
####################################
writeLines("\n*** Remove the temporal Seurat objects ***\n")

rm(seurat.object.u)
rm(UnfilteredData.df)
rm(seurat.object.full)
rm(seurat.object.subset)

####################################
### Normalize data (if applicable)
####################################
if (NormalizeAndScale == 1) {
  
  writeLines("\n*** Normalize data using NormalizeData() ***\n")
  
  StopWatchStart$NormalizeData  <- Sys.time()
  
  seurat.object.f <- NormalizeData(object = seurat.object.f, normalization.method = DefaultParameters$NormalizationMethod, scale.factor = DefaultParameters$ScaleFactor)
  
  StopWatchEnd$NormalizeData  <- Sys.time()
  
}else if (NormalizeAndScale == 2) {
    
    writeLines("\n*** Normalize data using NormalizeData() ***\n")
    
    StopWatchStart$SCTransform  <- Sys.time()
    
    seurat.object.f <- SCTransform(object = seurat.object.f, verbose = T)
    
    StopWatchEnd$SCTransform  <- Sys.time()
    
}else if (NormalizeAndScale == 3) {
  
  writeLines("\n*** Data assumed to be already normalized ***\n")
  
} else {
  stop(paste("Unexpected option -b", NormalizeAndScale, ". Only options '1', '2', or '3' are allowed. \n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  suppressPackageStartupMessages(library(loomR)) 
  
  # output the normalized count matrix for violin plots
  writeLines("\n*** Outputting normalized count matrix as loom ***\n")
  
  normalized_count_matrix <- as.matrix(seurat.object.f@assays[["RNA"]]@data)
  
  features_tsv <- as.data.frame(rownames(normalized_count_matrix))
  write.table(features_tsv, file=paste(Tempdir,"/","raw/","features.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  loom_file <- paste(Tempdir,"/","normalized/","normalized_counts.loom", sep="")
  create(loom_file, normalized_count_matrix)
  
}


####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes ***\n")

StopWatchStart$FindVariableGenes  <- Sys.time()

### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object.f <- FindVariableFeatures(object = seurat.object.f, selection.method = DefaultParameters$VarGeneDetectSelectMethod, mean.cutoff = c(DefaultParameters$XLowCutoff, DefaultParameters$XHighCutoff), dispersion.cutoff = c(DefaultParameters$YCutoff, Inf))
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

seurat.object.f <- ScaleData(object = seurat.object.f, vars.to.regress = c("nCount_RNA", "mito.fraction"), features = rownames(x = seurat.object.f), verbose = T)

StopWatchEnd$ScaleData  <- Sys.time()

####################################
### Perform linear dimensional reduction by PCA
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

options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.f <- FindNeighbors(object = seurat.object.f, dims = PcaDimsUse) ## This step was part of FindClusters() in Seurat v2
seurat.object.f <- FindClusters(object = seurat.object.f, resolution = Resolution)

StopWatchEnd$ClusterCells  <- Sys.time()

StopWatchStart$CellClusterTables  <- Sys.time()

CellNames<-rownames(seurat.object.f@meta.data)
ClusterIdent<-seurat.object.f@meta.data$seurat_clusters
Headers<-paste("Cell_barcode", paste("seurat_cluster_resolution", Resolution, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  OutfileClusters<-paste(Tempdir,"/","groups.tsv", sep="")
} else {
  OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_CellClusters.tsv", sep="")
}
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
  ### `Error in .check_tsne_params(nrow(X), dims = dims, perplexity = perplexity,  : 
  ### perplexity is too large for the number of samples`
  ### User can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
  
  if (("tsne" %in% dim_red_method) & (length(colnames(seurat.object.f)) < DefaultParameters$MinNumberOfCellsToReducePerplexity)) {
    writeLines(paste("\n*** Using reduced perplexity = ", DefaultParameters$ReducedPerplexity, " because found ",  length(colnames(seurat.object.f)), " cells", " ***\n", sep = "", collapse = ""))
    seurat.object.f <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.f, dims = PcaDimsUse, perplexity = DefaultParameters$ReducedPerplexity)
  }else{
    seurat.object.f <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.f, dims = PcaDimsUse)
  }

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
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    OutfileCoordinates<-paste(Tempdir,"/","coordinates/",DimensionReductionMethods[[dim_red_method]][["name"]], "Coordinates.tsv", sep="", collapse = "")
  } else {
    OutfileCoordinates<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Coordinates.tsv", sep="", collapse = "")
  }
  
  Headers<-paste("Barcode",paste(colnames(seurat.object.f@reductions[[dim_red_method]]@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
  write.table(Headers,file = OutfileCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(seurat.object.f@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

  StopWatchEnd$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction ***\n", sep = "", collapse = ""))
  
  StopWatchStart$QCDimRedPlots$dim_red_method  <- Sys.time()
  
  CellPropertiesToColour<-c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction")
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
    print(FeaturePlot(object = seurat.object.f, ncol = nColFeaturePlot, features = c(ListOfGenesForDimRedPlots), cols = c("lightgrey", "blue"), reduction = dim_red_method, order = T))
    dev.off()

    StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()

  }
}

####################################
### Finding differentially expressed genes for each cell cluster
####################################
writeLines("\n*** Finding differentially expressed genes for each cell cluster ***\n")

print(paste("NumberOfClusters=", NumberOfClusters, sep = "", collapse = ""))

StopWatchStart$FindDiffMarkers  <- Sys.time()

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.f@meta.data))

seurat.object.markers <- FindAllMarkers(object = seurat.object.f, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)

SimplifiedDiffExprGenes.df <- seurat.object.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC")]
write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_MarkersPerCluster.tsv",sep=""), row.names = F, sep="\t", quote = F)

### Get top-2 genes sorted by cluster, then by avg_logFC
top_genes_by_cluster_for_tsne<-(seurat.object.markers %>% group_by(cluster) %>% top_n(DefaultParameters$NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers[["cluster"]]))

StopWatchEnd$FindDiffMarkers  <- Sys.time()

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  markers_file <- top_genes_by_cluster_for_tsne[,c("gene","cluster","p_val","avg_logFC")]
  write.table(markers_file, paste(Tempdir,"/","markers/","TopTwoMarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)
} 

####################################
### Saving the R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving the R object ***\n")
  
  StopWatchStart$SaveRDS  <- Sys.time()
  
  OutfileRDS<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_object.rds", sep="")
  saveRDS(seurat.object.f, file = OutfileRDS)
  
  StopWatchEnd$SaveRDS  <- Sys.time()
  
}else{
  
  writeLines("\n*** Not saving R object ***\n")
  
}

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
  print(FeaturePlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), cols = c("lightgrey", "blue"), reduction = dim_red_method, order = T))
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

### Make this option "N" if library(staplr) or pdftk are not available
SummaryPlots <- "Y"

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

ListOfPdfFilesToPng<-c(paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_PCElbowPlot.pdf", sep = "", collapse = "")
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
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo.txt", sep="")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt", sep="")
write(file = OutfileCPUusage, x = paste("Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
Headers<-paste("Step", "Time(minutes)", sep="\t")
write.table(Headers,file = OutfileCPUusage, row.names = F, col.names = F, sep="\t", quote = F, append = T)

for (stepToClock in names(StopWatchStart)) {
  if (regexpr("POSIXct", class(StopWatchStart[[stepToClock]]), ignore.case = T)[1] == 1) {
    TimeStart <- StopWatchStart[[stepToClock]]
    TimeEnd   <- StopWatchEnd[[stepToClock]]
    TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
    ReportTime<-c(paste(stepToClock, TimeDiff, sep = "\t", collapse = ""))
    write(file = OutfileCPUusage, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
  }else{
    ### This is NOT being printed out
    ### Because the StopWatch[Start|End]$STEP$SUB_STEP
    ### need to be split and programmed are using the word provided by SUB_STEP itself as key
    ### need to pass the SUB_STEP value instead
    for (substep in rownames(DimensionReductionMethods)) {
      if (regexpr("POSIXct", class(StopWatchStart[[stepToClock]][[substep]]), ignore.case = T)[1] == 1) {
        TimeStart <- StopWatchStart[[stepToClock]][[substep]]
        TimeEnd   <- StopWatchEnd[[stepToClock]][[substep]]
        TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
        ReportTime<-c(paste(stepToClock, TimeDiff, substep, sep = "\t", collapse = ""))
        write(file = OutfileCPUusage, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
      }
    }
  }
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
