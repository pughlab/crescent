####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/v3.0/merge_vignette.html
### and https://satijalab.org/seurat/v3.0/sctransform_vignette.html
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
suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
### Seurat v3 can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
suppressPackageStartupMessages(library(earlycross))   # to handle reading and writing mtx files
### Which can be installed like:
### install.packages('devtools')
### devtools::install_github("daskelly/earlycross")
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
### 'UMAP' can be installed using `pip install umap-learn`
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
  make_option(c("-i", "--inputs_list"), default="NA",
              help="Path/name to a <tab> delimited file with the list of dataset IDs, infile path/names, and dataset type (e.g. batch) like:
                dataset1_id  /path_to/dataset1_mtx  dataset1_type  dataset_format
                dataset2_id  /path_to/dataset2_mtx  dataset2_type  dataset_format
                dataset3_id  /path_to/dataset3_mtx  dataset3_type  dataset_format
                ...etc
                
                Important note:
                The order of the list of datasets in --inputs_list influences the results,
                including number of clusters, t-SNE/UMAP and differentially expressed genes
                
                Each dataset input must be in either format 'MTX' or 'DGE':
                a) the path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files; or
                b) the path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns
                Notes:
                      The 'MTX' files can be the outputs from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                      Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz

                Datasets will be normalized using SC transform and three levels of integration will be used for clustering:
                1) cluster cells from each dataset
                2) integrate cells by 'dataset_type' (column 3) and cluster them
                3) integrate cells from all datasets and cluster them
                
                Default = 'No default. It's mandatory to specify this parameter'"),
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
  make_option(c("-a", "--opacity"), default="0.1",
              help="If using a --list_genes, this parameter provides a value for the minimal opacity of gene expression. Use a value between 0 and 1
                Default = 'y'"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file:
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                Default = '10'"),
  #
  make_option(c("-m", "--percent_mito"), default="0,0.05",
              help="<comma> delimited min,max number of percentage of mitochondrial gene counts in a cell to be included in normalization and clustering analyses
               For example, for whole cell scRNA-seq use '0,0.2', or for Nuc-seq use '0,0.05'
               Default = '0,0.05'"),
  #
  make_option(c("-n", "--n_genes"), default="50,8000",
              help="<comma> delimited min,max number of unique gene counts in a cell to be included in normalization and clustering analyses
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

InputsList              <- opt$inputs_list
Resolution              <- as.numeric(opt$resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
InfileColourDimRedPlots <- opt$infile_colour_dim_red_plots
ListGenes               <- opt$list_genes
Opacity                 <- as.numeric(opt$opacity)
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
ListPMito               <- opt$percent_mito
ListNGenes              <- opt$n_genes
ThreshReturn            <- as.numeric(opt$return_threshold)
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

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
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
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
  
  ### Parameters for Cluster Biomarkers
  FindAllMarkers.MinPct     =  0.25,
  FindAllMarkers.ThreshUse  =  0.25,

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

ListMandatory<-list("infiles_list", "inputs_type", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

### Need to program this to start from an R object
# seurat.object.integrated <- readRDS(file = "~/SINGLE_CELL/10X/CHATURVEDI_LAB/Seurat_Integration/SEURAT_MERGED_LV_EXCLUDING_AO/R_0.4/SEURAT/Rajiv_samples_res0.4.SEURAT_integrated_object.rds")

InputsList<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), InputsList)
InputsTable<-read.table(InputsList, header = F, row.names = 1, stringsAsFactors = FALSE)
colnames(InputsTable)<-c("PathToDataset","DatasetType","DatasetFormat")

SeuratObjects        <-list()
DatasetIds           <-list()
list_DatasetToType   <-list()
list_TypeToDatasets  <-list()
list_DatasetToFormat <-list()

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  Dataset.SO <-paste(dataset, ".so",    sep = "", collapse = "")
  
  PathToDataset <- InputsTable[dataset,"PathToDataset"]
  PathToDataset <- gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), PathToDataset)
  DatasetType   <- InputsTable[dataset,"DatasetType"]
  DatasetFormat <- InputsTable[dataset,"DatasetFormat"]
  list_DatasetToType[[dataset]] <- DatasetType
  list_DatasetToFormat[[dataset]] <- DatasetFormat
  list_TypeToDatasets[[DatasetType]] <- append(list_TypeToDatasets[[DatasetType]], dataset)

  if (regexpr("^MTX$|^DGE$", DatasetFormat, ignore.case = T, perl = T)[1] == 1) {

    ####################################
    ### Loading MTX or DGE infiles
    ####################################
  
    StopWatchStart$LoadScRNAseqData$dataset <- Sys.time()
  
    if (regexpr("^MTX$", DatasetFormat, ignore.case = T)[1] == 1) {
      writeLines(paste("\n*** Loading MTX infiles for ", dataset, " from: ", PathToDataset, " ***\n", sep = "", collapse = ""))
      expression_matrix <- Read10X(data.dir = PathToDataset)
    }else if (regexpr("^DGE$", DatasetFormat, ignore.case = T)[1] == 1) {
      print("Loading Digital Gene Expression matrix")
      writeLines(paste("\n*** Loading Digital Gene Expression matrix for ", dataset, " from: ", PathToDataset, " ***\n", sep = "", collapse = ""))
      ## Note `check.names = F` is needed for both `fread` and `data.frame`
      expression_matrix <- as.matrix(data.frame(fread(PathToDataset, check.names = F), row.names=1, check.names = F))
    }
    dim(expression_matrix)
  
    StopWatchEnd$LoadScRNAseqData$dataset  <- Sys.time()

    ####################################
    ### Create seurat object
    ####################################
    
    StopWatchStart$CreateSeuratObject$dataset  <- Sys.time()
    
    writeLines(paste("\n*** Create seurat object for ", dataset, " ***\n", sep = "", collapse = ""))
    seurat.object.u  <- CreateSeuratObject(counts = expression_matrix, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = paste(PrefixOutfiles, "_", dataset, sep = "", collapse = ""))
    
    StopWatchEnd$CreateSeuratObject$dataset  <- Sys.time()
    
    ####################################
    ### Add dataset label
    ####################################
    writeLines(paste("\n*** Add dataset label for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$AddDatasetLabel$dataset  <- Sys.time()
    
    seurat.object.u[['dataset.label']] <- dataset
    
    StopWatchEnd$AddDatasetLabel$dataset  <- Sys.time()

    ####################################
    ### Add dataset type label
    ####################################
    writeLines(paste("\n*** Add dataset type label for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$AddDatasetLabel$dataset  <- Sys.time()
    
    seurat.object.u[['dataset_type.label']] <- list_DatasetToType[[dataset]]
    
    StopWatchEnd$AddDatasetLabel$dataset  <- Sys.time()

    ####################################
    ### Get mitochondrial genes
    ####################################
    writeLines(paste("\n*** Get  mitochondrial genes for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$GetMitoGenes$dataset  <- Sys.time()
    
    mitoRegExpressions<- paste(c("^MT-"), collapse = "|")
    mito.features <- grep(pattern = mitoRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)
    
    if (length(mito.features)[[1]] > 0) {
      percent.mito <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.mito']] <- percent.mito
    }else{
      percent.mito <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.mito']] <- percent.mito
    }
    
    StopWatchEnd$GetMitoGenes$dataset  <- Sys.time()
    
    ####################################
    ### Get ribosomal protein genes
    ####################################
    writeLines(paste("\n*** Get ribosomal protein genes for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$GetRiboGenes$dataset  <- Sys.time()
    
    riboRegExpressions<- paste(c("^MRPL", "^MRPS", "^RPL", "^RPS"),collapse = "|")
    ribo.features <- grep(pattern = riboRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)
    
    if (length(ribo.features)[[1]] > 0) {
      percent.ribo <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.ribo']] <- percent.ribo
    }else{
      percent.ribo <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.ribo']] <- percent.ribo
    }
    
    StopWatchEnd$GetRiboGenes$dataset  <- Sys.time()
    
    ####################################
    ### Filter cells based gene counts and mitochondrial representation
    ####################################
    writeLines(paste("\n*** Filter cells based gene counts and mitochondrial representation for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$FilterCells$dataset  <- Sys.time()
    
    if (length(mito.features)[[1]] > 0) {
      seurat.object.integrated<-subset(x = seurat.object.u, subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes & percent.mito > DefaultParameters$MinPMito & percent.mito < DefaultParameters$MaxPMito)
    }else{
      seurat.object.integrated<-subset(x = seurat.object.u, subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes)
    }
    
    StopWatchEnd$FilterCells$dataset  <- Sys.time()
    
    ### Just reporting the summary of the UNfiltered and filtered objects
    seurat.object.u
    seurat.object.integrated
    
    ####################################
    ### QC EDA violin plots
    ####################################
    writeLines(paste("\n*** QC EDA violin plots for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$QCviolinplots$dataset  <- Sys.time()
    
    ### Get unfiltered data QC statistics
    nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
    nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
    percent.mito.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$percent.mito, percent.mito = 1)
    percent.ribo.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$percent.ribo, percent.ribo = 1)
    #
    nFeature_RNAStats.u<-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
    nCount_RNAStats.u  <-paste(c( "mean = ",round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
    percent.mito.u     <-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object.u@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
    percent.ribo.u     <-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object.u@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")
    
    ### Get filtered data QC statistics
    nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object.integrated@meta.data$nFeature_RNA, nGenes = 2)
    nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object.integrated@meta.data$nCount_RNA,   nCount_RNA = 2)
    percent.mito.f.df  <-data.frame(Expression_level = seurat.object.integrated@meta.data$percent.mito, percent.mito = 2)
    percent.ribo.f.df  <-data.frame(Expression_level = seurat.object.integrated@meta.data$percent.ribo, percent.ribo = 2)
    #
    nFeature_RNAStats.f<-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object.integrated@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
    nCount_RNAStats.f  <-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object.integrated@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
    percent.mito.f     <-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object.integrated@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
    percent.ribo.f     <-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object.integrated@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")
    
    ### Put QC statistics together
    nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
    nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
    percent.mito.m.df  <-data.frame(rbind(percent.mito.u.df,percent.mito.f.df))
    percent.ribo.m.df  <-data.frame(rbind(percent.ribo.u.df,percent.ribo.f.df))
    LabelUnfiltered    <-paste("Before filters: No. of cells = ", nrow(seurat.object.u@meta.data), sep ="", collapse = "")
    LabelFiltered      <-paste("After filters:  No. of cells = ", nrow(seurat.object.integrated@meta.data), sep ="", collapse = "")
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
    
    VlnPlotPdf<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_QC_VlnPlot.pdf", sep="", collapse = )
    pdf(file=VlnPlotPdf, width = 12, height = 7)
    print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
    dev.off()
    
    StopWatchEnd$QCviolinplots$dataset  <- Sys.time()
    
    ####################################
    ### Feature-vs-feature scatter plot
    ####################################
    writeLines(paste("\n*** Feature-vs-feature scatter plot for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$FeatureVsFeatureplot$dataset  <- Sys.time()
    
    UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                                  nGene = seurat.object.u@meta.data$nFeature_RNA,
                                  percent.mito = seurat.object.u@meta.data$percent.mito,
                                  filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.integrated))
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = ColoursQCViolinPlots[[1]])
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = FALSE, replacement = ColoursQCViolinPlots[[2]])
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = 4)
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = FALSE, replacement = 16)
    
    FeatureVsFeaturePlotPdf<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_NumbReadsVsNumbGenesAndMito_VlnPlot.pdf", sep="", collapse = "")
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
    
    StopWatchEnd$FeatureVsFeatureplot$dataset  <- Sys.time()
    
    ####################################
    ### Remove the Unfiltered seurat object
    ####################################
    writeLines(paste("\n*** Remove the Unfiltered seurat object for ", dataset, " ***\n", sep = "", collapse = ""))
    
    rm(seurat.object.u)
    rm(UnfilteredData.df)
    
    ####################################
    ### Assign data to Datasets lists
    ####################################
    writeLines(paste("\n*** Assign data to Datasets lists: ", dataset, " ***\n", sep = "", collapse = ""))
    
    SeuratObjects[[as.character(NumberOfDatasets)]]  <- seurat.object.integrated
    DatasetIds[[as.character(NumberOfDatasets)]]     <- dataset
    
  }else{
    stop(paste("Unexpected type of input: ", DatasetType, "\n\nFor help type:\n\nintegrates_datasets_with_seurat.R -h\n\n", sep=""))
  }
}

####################################
### Merge Seurat objects
####################################
writeLines("\n*** Merge Seurat objects ***\n")

StopWatchStart$MergeSeuratObjects  <- Sys.time()

FirstSeuratObject <- SeuratObjects[[1]]
FirstSampleId     <- DatasetIds[[1]]
RestOfSeuratObjects <- SeuratObjects[c(2:NumberOfDatasets)]
RestOfSamplesIds    <- unlist(DatasetIds[c(2:NumberOfDatasets)])

seurat.object.merged <- merge(FirstSeuratObject, y = RestOfSeuratObjects, 
                              add.cell.ids = DatasetIds,
                              project = PrefixOutfiles)

StopWatchEnd$MergeSeuratObjects  <- Sys.time()

dataset.label.metadata<-cbind.data.frame(sample=seurat.object.merged@meta.data$dataset.label)
rownames(dataset.label.metadata)<-colnames(seurat.object.merged)

seurat.object.merged.withmetadata <- CreateSeuratObject(seurat.object.merged@assays$RNA@counts, meta.data = dataset.label.metadata)
seurat.object.list <- SplitObject(seurat.object.merged.withmetadata, split.by = "sample")

####################################
### Running SCTransform
####################################
writeLines("\n*** Running SCTransform ***\n")

StopWatchStart$SCTransform  <- Sys.time()

for (i in 1:length(seurat.object.list)) {
  seurat.object.list[[i]] <- SCTransform(seurat.object.list[[i]], verbose = T)
}

StopWatchEnd$SCTransform  <- Sys.time()

####################################
### Integrating datasets
####################################
writeLines("\n*** Integrating datasets ***\n")

StopWatchStart$Integration  <- Sys.time()

seurat.object.integratedfeatures <- SelectIntegrationFeatures(object.list = seurat.object.list, nfeatures = 3000)
seurat.object.list <- PrepSCTIntegration(object.list = seurat.object.list, anchor.features = seurat.object.integratedfeatures, verbose = FALSE)
seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, normalization.method = "SCT", anchor.features = seurat.object.integratedfeatures, verbose = FALSE)
seurat.object.integrated <- IntegrateData(anchorset = seurat.object.anchors, normalization.method = "SCT", verbose = FALSE)

StopWatchEnd$Integration  <- Sys.time()


seurat.object.integrated


####################################
### Reducing dimensions
####################################
writeLines("\n*** Reducing dimensions ***\n")

StopWatchStart$RunPCA  <- Sys.time()

seurat.object.integrated <- RunPCA(seurat.object.integrated, verbose = FALSE)

StopWatchEnd$RunPCA  <- Sys.time()

StopWatchStart$PCAPlots  <- Sys.time()

ForElbowPlot<-ElbowPlot(object = seurat.object.integrated, ndims = 50, reduction = "pca")
MaxYAxis<-as.integer(max(ForElbowPlot$data$stdev)+1)

pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_PCElbowPlot_Integrated.pdf", sep=""))
print(ForElbowPlot
      + scale_x_continuous(breaks =  seq(from = 0, to = 50, by=5))
      + geom_vline(xintercept = seq(from = 0, to = 50, by=5), linetype='dotted', col="red")
      + scale_y_continuous(breaks =  seq(from = 0, to = MaxYAxis, by=0.5))
      + geom_hline(yintercept = seq(from = 0, to = MaxYAxis, by=0.5), linetype='dotted', col="red")
)
dev.off()

StopWatchEnd$PCAPlots  <- Sys.time()

####################################
### Globally cluster cells using integrated data
####################################
writeLines("\n*** Globally cluster cells using integrated data ***\n")

StopWatchStart$ClusterAllCells  <- Sys.time()

options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.integrated <- FindNeighbors(object = seurat.object.integrated, dims = PcaDimsUse)
seurat.object.integrated <- FindClusters(object = seurat.object.integrated, reduction.type = "pca", resolution = Resolution, print.output = 0, save.SNN = T)

StopWatchEnd$ClusterAllCells  <- Sys.time()

StopWatchStart$AllCellClusterTables  <- Sys.time()

CellNames<-rownames(seurat.object.integrated@meta.data)
ClusterIdent <-seurat.object.integrated@meta.data$seurat_clusters
NumberOfClusters<-length(unique(ClusterIdent))

Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)

StopWatchEnd$AllCellClusterTables  <- Sys.time()

####################################
### Saving global cell cluster identities
####################################
writeLines("\n*** Saving global cell cluster identities ***\n")

seurat.object.integrated$GlobalCellClusterIdentities <- Idents(object = seurat.object.integrated)

####################################
### Get average gene expression for each global cluster
####################################
writeLines("\n*** Get average gene expression for each global cluster ***\n")

StopWatchStart$AverageGeneExpression  <- Sys.time()

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)
#
OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllSamples_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllSamples_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllSamples_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("r", Resolution, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$AverageGeneExpression  <- Sys.time()

####################################
### Finding differentially expressed genes for each global cell cluster
####################################
writeLines("\n*** Finding differentially expressed genes for each global cell cluster ***\n")

StopWatchStart$FindDiffMarkers  <- Sys.time()

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))

seurat.object.integrated.markers <- FindAllMarkers(object = seurat.object.integrated, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)

write.table(data.frame("GENE"=rownames(seurat.object.integrated.markers),seurat.object.integrated.markers),paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "MarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)

StopWatchEnd$FindDiffMarkers  <- Sys.time()

####################################
### FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data
####################################
writeLines("\n*** FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {
  ####################################
  ### Run non-linear dimension reductions using integrated data
  ####################################
  writeLines(paste("\n*** Run ", DimensionReductionMethods[[dim_red_method]][["name"]], " ***\n", sep = "", collapse = ""))

  StopWatchStart$DimensionReduction$dim_red_method  <- Sys.time()

  ### NOTES:
  ### In RunTSNE: if the datasets is small user may get error:
  ### `Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large.`
  ### User can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
  
  seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, do.fast = T)
  
  StopWatchEnd$DimensionReduction$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plot by global cell clusters using integrated data
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by global cell clusters ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimRedOPlotColourByCellCluster$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("seurat_clusters"), combine = FALSE, reduction = dim_red_method, label = T, label.size = 10)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
  IntegratedDimRedPlotPdf<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir,  "_GlobalClustering_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$DimRedOPlotColourByCellCluster$dim_red_method  <- Sys.time()

  ####################################
  ### Write out dimension reduction plot coordinates using integrated data
  ####################################
  writeLines(paste("\n*** Write out ", DimensionReductionMethods[[dim_red_method]][["name"]], " coordinates ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  OutfileCoordinates<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_Coordinates.tsv", sep="", collapse = "")
  Headers<-paste("Barcode",paste(colnames(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
  write.table(Headers,file = OutfileCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by dataset using integrated data
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by dataset ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimRedPlotsByDataset$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("sample"), combine = FALSE, reduction = "umap")
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
  IntegratedDimRedPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByDataset.pdf", sep="")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$DimRedPlotsByDataset$dim_red_method  <- Sys.time()
  
  
  ####################################
  ### Colour dimension reduction plots by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo using integrated data
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo ***\n", sep = "", collapse = ""))
  
  StopWatchStart$QCDimRedPlots$dim_red_method  <- Sys.time()
  
  CellPropertiesToColour<-c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByQC.pdf", sep=""), width = DefaultParameters$BaseSizeSinglePlotPdf * 1.5, height = DefaultParameters$BaseSizeSinglePlotPdf * 1.5)
  print(FeaturePlot(object = seurat.object.integrated, label = T, order = T, features = CellPropertiesToColour, cols = c("lightgrey", "blue"), reduction = dim_red_method, pt.size = 1.5, ncol = 2))
  dev.off()

  StopWatchEnd$QCDimRedPlots$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by -infile_colour_dim_red_plots using integrated data
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by -infile_colour_dim_red_plots ***\n", sep = "", collapse = ""))
  
  if (regexpr("^NA$", InfileColourDimRedPlots, ignore.case = T)[1] == 1) {
    print("No extra barcode-attributes will be used for dimension reduction plots")
  }else{
    
    StopWatchStart$DimRedPlotsColuredByMetadata$dim_red_method  <- Sys.time()
    
    seurat.object.meta.data<-seurat.object.integrated@meta.data
    ExtraCellProperties <- data.frame(read.table(InfileColourDimRedPlots, header = T, row.names = 1))
    print(head(ExtraCellProperties))
    
    # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
    # so that barcodes from --infile_colour_dim_red_plots and --input can match each other
    rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
    seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = ExtraCellProperties)
    
    # Generating outfile
    # Note DimPlot() takes the entire current device (pdf)
    # even if using layout(matrix(...)) or  par(mfrow=())
    # Thus each property plot is written to a separate page of a single *pdf outfile
    pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByExtraProperties.pdf", sep=""), width = DefaultParameters$BaseSizeSinglePlotPdf, height = DefaultParameters$BaseSizeSinglePlotPdf)
    for (property in colnames(ExtraCellProperties)) {
      if ( (sum(ExtraCellProperties[,property] %in% 0:1 == T)) == (nrow(ExtraCellProperties)) ) { ## is binary
        CellsToHighlight <- rownames(ExtraCellProperties)[ExtraCellProperties[,property]==1]
        print(DimPlot(object = seurat.object.integrated, reduction = dim_red_method, group.by = property, combine = T, legend = "none", cells.highlight = CellsToHighlight) + ggtitle(property))
      }else{
        print(DimPlot(object = seurat.object.integrated, reduction = dim_red_method, group.by = property, combine = T, legend = "none") + ggtitle(property))
      }
    }
    dev.off()
    
    StopWatchEnd$DimRedPlotsColuredByMetadata$dim_red_method  <- Sys.time()
    
  }
  
  ####################################
  ### Colour dimension reduction plots showing each requested gene using integrated data
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
    
    ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
    ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
    seurat.object.integrated.sa <- CreateSeuratObject(seurat.object.integrated@assays$RNA@data)
    seurat.object.integrated.sa@reductions$umap <- seurat.object.integrated@reductions$umap
    seurat.object.integrated.sa@reductions$tsne <- seurat.object.integrated@reductions$tsne
    
    pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
    print(FeaturePlot(object = seurat.object.integrated.sa, ncol = nColFeaturePlot, features = ListOfGenesForDimRedPlots, cols = c("lightgrey", "blue"),
                      reduction = dim_red_method, order = T, slot = "data", pt.size = 1.5, min.cutoff = "q0.1", max.cutoff = "q90"))
    dev.off()
    
    rm(seurat.object.integrated.sa)
    
    StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
    
  }
}

####################################
### FOR EACH SAMPLE merge sample_id and GLOBAL CLUSTER identities to make sample-specific identities
### Then get DGE and colour dimension reduction plots
####################################
writeLines("\n*** FOR EACH SAMPLE merge sample_id and GLOBAL CLUSTER to make sample-specific identities. Then get DGE and colour dimension reduction plots  ***\n")

EachSampleGlobalClusteredCellClusters <- unlist(x = strsplit(x = paste(seurat.object.integrated$sample, seurat.object.integrated$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleGlobalClusteredCellClusters, col.name = "EachSampleGlobalClusteredCellClusters")

# switch the identity class of all cells to reflect sample-specific identities
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleGlobalClusteredCellClusters

####################################
### Get average gene expression for each sample based on global clusters
####################################
writeLines("\n*** Get average gene expression for each sample based on global clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSample_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSample_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSample_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Colour dimension reduction plots for each sample based on global clusters
####################################

writeLines("\n*** Colour dimension reduction plots for each sample based on global clusters ***\n")

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_sample") == T) {
    rm(seurat.object.each_sample)
  }
  seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_sample)

  for (dim_red_method in names(DimensionReductionMethods)) {
    
    StopWatchStart$DimRedPlotsByDatasetIntegratedCellClusters$dataset$dim_red_method  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_sample, group.by = c("seurat_clusters"), combine = FALSE, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
    
    StopWatchEnd$DimRedPlotsByDatasetIntegratedCellClusters$dataset$dim_red_method  <- Sys.time()
    
    ### To program layout() for more than 3 genes in multiple rows
    if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
      print("No selected genes for dimension reduction plots")
    }else{

      StopWatchStart$DimRedPlotsByDatasetColuredByGenes$dim_red_method  <- Sys.time()
      
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
      
      ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
      ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
      seurat.object.each_sample.sa <- CreateSeuratObject(seurat.object.each_sample@assays$RNA@data)
      seurat.object.each_sample.sa@reductions$umap <- seurat.object.each_sample@reductions$umap
      seurat.object.each_sample.sa@reductions$tsne <- seurat.object.each_sample@reductions$tsne
      
      pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
      print(FeaturePlot(object = seurat.object.each_sample.sa, ncol = nColFeaturePlot, features = ListOfGenesForDimRedPlots, cols = c("lightgrey", "blue"),
                        reduction = dim_red_method, order = T, slot = "data", pt.size = 1.5, min.cutoff = "q0.1", max.cutoff = "q90"))
      dev.off()
      
      StopWatchEnd$DimRedPlotsByDatasetColuredByGenes$dim_red_method  <- Sys.time()
  
      rm(seurat.object.each_sample.sa)
      
    }
  }
}

####################################
### Colour dimension reduction plots for each sample by -infile_colour_dim_red_plots using integrated data
####################################

if (regexpr("^NA$", InfileColourDimRedPlots, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for dimension reduction plots")
}else{
  writeLines("\n*** Colour dimension reduction plots for each sample by -infile_colour_dim_red_plots using integrated data ***\n")
    
  seurat.object.meta.data<-seurat.object.integrated@meta.data
  ExtraCellProperties <- data.frame(read.table(InfileColourDimRedPlots, header = T, row.names = 1))
  print(head(ExtraCellProperties))
  
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_dim_red_plots and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = ExtraCellProperties)

  NumberOfDatasets <- 0
  for (dataset in rownames(InputsTable)) {
    NumberOfDatasets <- NumberOfDatasets + 1
    print(NumberOfDatasets)
    
    if (exists(x = "seurat.object.each_sample") == T) {
      rm(seurat.object.each_sample)
    }
    seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
    print(seurat.object.each_sample)
    
    for (dim_red_method in names(DimensionReductionMethods)) {
      
      StopWatchStart$DimRedPlotsByMetadata$dataset$dim_red_method  <- Sys.time()
  
      # Generating outfile
      # Note DimPlot() takes the entire current device (pdf)
      # even if using layout(matrix(...)) or  par(mfrow=())
      # Thus each property plot is written to a separate page of a single *pdf outfile
      pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByExtraProperties.pdf", sep=""), width = DefaultParameters$BaseSizeSinglePlotPdf, height = DefaultParameters$BaseSizeSinglePlotPdf)
      for (property in colnames(ExtraCellProperties)) {
        
        if ( (sum(ExtraCellProperties[,property] %in% 0:1 == T)) == (nrow(ExtraCellProperties)) ) { ## is binary
          CellsToHighlight <- rownames(ExtraCellProperties)[ExtraCellProperties[,property]==1]
          print(DimPlot(object = seurat.object.each_sample, reduction = dim_red_method, group.by = property, combine = T, legend = "none", cells.highlight = CellsToHighlight) + ggtitle(property))
        }else{
          print(DimPlot(object = seurat.object.each_sample, reduction = dim_red_method, group.by = property, combine = T, legend = "none") + ggtitle(property))
        }
      }
      dev.off()
      
      StopWatchEnd$DimRedPlotsByMetadata$dataset$dim_red_method  <- Sys.time()
        
    }
  }
}


####################################
### FOR EACH SAMPLE uses integrated counts and RE-CLUSTER cells
### Then get DGE and colour dimension reduction plots 
####################################
writeLines("\n*** FOR EACH SAMPLE uses integrated counts and RE-CLUSTER cells. Then get DGE and colour dimension reduction plots ***\n")

# switch the identity class of all cells to reflect sample
Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample

# Headers for sample-specific cell clusters
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = ""), sep="\t")
OutfileEachSampleClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "CellClusters.tsv", sep="")
write.table(Headers,file = OutfileEachSampleClusters, row.names = F, col.names = F, sep="\t", quote = F)

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_sample") == T) {
    rm(seurat.object.each_sample)
  }
  seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_sample)

  StopWatchStart$ClusterEachSampleCellsFromInteg$dataset  <- Sys.time()
  
  ####################################
  ### Re-cluster each sample cells using integrated data
  ####################################
  
  options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
  seurat.object.each_sample <- FindNeighbors(object = seurat.object.each_sample, dims = PcaDimsUse)
  seurat.object.each_sample <- FindClusters(object = seurat.object.each_sample, reduction.type = "pca", resolution = Resolution, print.output = 0, save.SNN = T)
  ClustersThisSample <- unlist(x = strsplit(x = paste(dataset, seurat.object.each_sample@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
  
  StopWatchEnd$ClusterEachSampleCellsFromInteg$dataset  <- Sys.time()
  
  ####################################
  ### Write out each sample cell clusters
  ####################################
  writeLines("\n*** Write out each sample cell clusters ***\n")
  
  StopWatchStart$WriteClustersEachSampleCellClusterTables$dataset <- Sys.time()
  
  write.table(paste(colnames(seurat.object.each_sample),ClustersThisSample, sep = "\t", collapse = "\n"),
              file = OutfileEachSampleClusters, row.names = F, col.names = F, quote = F, append = T)

  CellNames<-rownames(seurat.object.each_sample@meta.data)
  ClusterIdent <-seurat.object.each_sample@meta.data$seurat_clusters
  NumberOfClusters<-length(unique(ClusterIdent))
  OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", dataset,"_NumbCellClusters", ".tsv", sep="")
  write(x=NumberOfClusters,file = OutfileNumbClusters)
  
  StopWatchEnd$WriteClustersEachSampleCellClusterTables$dataset <- Sys.time()
  
  ####################################
  ### Finding differentially expressed genes for each sample re-clustered cell clusters
  ####################################
  writeLines("\n*** Finding differentially expressed genes for each sample re-clustered cell clusters ***\n")
  
  StopWatchStart$FindDiffMarkers$dataset  <- Sys.time()
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_sample@meta.data))
  
  seurat.object.each_sample.markers <- FindAllMarkers(object = seurat.object.each_sample, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
  
  write.table(data.frame("GENE"=rownames(seurat.object.each_sample.markers),seurat.object.each_sample.markers), paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", dataset, "_MarkersPerCluster.tsv", sep=""),row.names = F,sep="\t",quote = F)
  
  StopWatchEnd$FindDiffMarkers$dataset  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots for each sample re-clustered cells
  ####################################
  writeLines("\n*** Colour dimension reduction plots for each sample re-clustered cells ***\n")
  
  for (dim_red_method in names(DimensionReductionMethods)) {
    
    StopWatchStart$DimRedPlotsByEachSampleReclusteredCellCluster$dataset$dim_red_method  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_sample, group.by = c("seurat_clusters"), combine = FALSE, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
    
    StopWatchEnd$DimRedPlotsByEachSampleReclusteredCellCluster$dataset$dim_red_method <- Sys.time()
    
  }
}

####################################
### Load each sample re-cluster assignments
####################################
writeLines("\n*** Load each sample re-cluster assignments ***\n")

EachSampleCellReClusters <- data.frame(read.table(OutfileEachSampleClusters, header = T, row.names = 1, sep = "\t"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleCellReClusters, col.name = "EachSampleCellReClusters")

# switch the identity class of all cells to reflect each sample clusters
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleCellReClusters

####################################
### Get average gene expression for each sample re-clustered clusters
####################################
writeLines("\n*** Get average gene expression for each sample re-clustered clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "PerSample_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "PerSample_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "PerSample_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)


####################################
### FOR EACH SAMPLE TYPE uses integrated/normalized counts and generates dimension reduction plots for each sample type using:
### a) sample type
### b) global clusters
### c) re-clustered cells of each sample - also gets DGE
####################################
writeLines("\n*** FOR EACH SAMPLE TYPE uses integrated/normalized counts and generates dimension reduction plots for each sample type using:
a) sample type
b) global clusters
c) re-clustered cells of each sample - also gets DGE  ***\n")

# Get sample_type assignments replacing sample ids using list_DatasetToType()
seurat.object.integrated$sample_type<-seurat.object.integrated$sample
for (dataset in rownames(InputsTable)) {
seurat.object.integrated$sample_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$sample_type)
}
SampleTypes <- unique(seurat.object.integrated$sample_type)

####################################
### Colour dimension reduction plots by sample type
####################################
writeLines("\n*** Colour dimension reduction plots by sample type ***\n")

# switch identity to sample type
Idents(object = seurat.object.integrated) <- seurat.object.integrated$sample_type

for (dim_red_method in names(DimensionReductionMethods)) {
  
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " coloured by sample_type ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimRedOPlotColourBySampleType$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("sample_type"), combine = FALSE, reduction = dim_red_method, label = F, label.size = 10)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
  IntegratedDimRedPlotPdf<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySampleType.pdf", sep="")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  print(IntegratedDimRedPlotPdf)
  
  StopWatchEnd$DimRedOPlotColourBySampleType$dim_red_method  <- Sys.time()
  
}

####################################
### Write out each sample type global cluster assignments
####################################
writeLines("\n*** Write out each sample type global cluster assignments ***\n")

# Headers for sample-type-specific cell clusters
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = ""), sep="\t")
OutfileEachSampleTypeGlobalClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "EachSampleType_CellClusters.tsv", sep="")
GlobalClustersBySampleType<- unlist(x = strsplit(x = paste(seurat.object.integrated@meta.data$sample_type, seurat.object.integrated@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
write.table(Headers,file = OutfileEachSampleTypeGlobalClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste(colnames(seurat.object.integrated),GlobalClustersBySampleType, sep = "\t", collapse = "\n"),
            file = OutfileEachSampleTypeGlobalClusters, row.names = F, col.names = F, quote = F, append = T)

####################################
### Load each sample type global cluster assignments
####################################
writeLines("\n*** Load each sample type global cluster assignments ***\n")

EachSampleTypeGlobalCellClusters <- data.frame(read.table(OutfileEachSampleTypeGlobalClusters, header = T, row.names = 1, sep = "\t"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleTypeGlobalCellClusters, col.name = "EachSampleTypeGlobalCellClusters")

# switch the identity class of all cells to reflect each sample clusters
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleTypeGlobalCellClusters

####################################
### Get average gene expression for each sample type using global clusters
####################################
writeLines("\n*** Get average gene expression for each sample type using global clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSampleType_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSampleType_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSampleType_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Colour dimension reduction plots for each sample type by glabal clusters and gets DGE
####################################
writeLines("\n*** Colour dimension reduction plots for each sample type by glabal clusters and gets DGE ***\n")

# switch identity to sample type global cluster
Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample_type

NumberOfSampleTypes <- 0
for (sample_type in SampleTypes) {
  NumberOfSampleTypes <- NumberOfSampleTypes + 1
  print(NumberOfSampleTypes)
  
  if (exists(x = "seurat.object.each_sample_type") == T) {
    rm(seurat.object.each_sample_type)
  }
  seurat.object.each_sample_type <- subset(x = seurat.object.integrated, idents = sample_type)

  ####################################
  ### Colour dimension reduction plots by each sample type using global cell clusters
  ####################################

  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for ", sample_type, " using global cell clusters ***\n", sep = "", collapse = ""))

  for (dim_red_method in names(DimensionReductionMethods)) {

    StopWatchStart$DimRedOPlotSampleTypeColourByCellCluster$sample_type$dim_red_method  <- Sys.time()

    plots <- DimPlot(seurat.object.each_sample_type, group.by = c("seurat_clusters"), combine = FALSE, reduction = dim_red_method, label = T, label.size = 10)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", sample_type,"_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()

    StopWatchEnd$DimRedOPlotSampleTypeColourByCellCluster$sample_type$dim_red_method  <- Sys.time()

  }

  ####################################
  ### Finding differentially expressed genes for each sample type using global cell clusters
  ####################################
  writeLines("\n*** Finding differentially expressed genes for each sample type using global cell clusters ***\n")
  
  Idents(object = seurat.object.each_sample_type) <- seurat.object.each_sample_type@meta.data$EachSampleTypeGlobalCellClusters

  StopWatchStart$FindDiffMarkers$sample_type  <- Sys.time()
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_sample_type@meta.data))

  seurat.object.each_sample_type.markers <- FindAllMarkers(object = seurat.object.each_sample_type, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)

  write.table(data.frame("GENE"=rownames(seurat.object.each_sample_type.markers),seurat.object.each_sample_type.markers), paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_" , sample_type, "_MarkersPerCluster.tsv", sep=""),row.names = F,sep="\t",quote = F)
  
  StopWatchEnd$FindDiffMarkers$sample_type  <- Sys.time()
}

####################################
### Reclusters each sample type cells and gets DGE
####################################

### Write out headers for sample type cell re-cluster outfile
OutfileEachSampleTypeReClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "CellClusters.tsv", sep="")
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = ""), sep="\t")
write.table(Headers,file = OutfileEachSampleTypeReClusters, row.names = F, col.names = F, sep="\t", quote = F)

NumberOfSampleTypes <- 0
for (sample_type in SampleTypes) {
  NumberOfSampleTypes <- NumberOfSampleTypes + 1
  print(NumberOfSampleTypes)
  
  if (exists(x = "seurat.object.each_sample_type") == T) {
    rm(seurat.object.each_sample_type)
  }
  seurat.object.each_sample_type <- subset(x = seurat.object.integrated, idents = sample_type)
  print(seurat.object.each_sample_type)
  
  ####################################
  ### Re-cluster each sample type cells using integrated data
  ####################################
  writeLines("\n*** Re-cluster each sample type cells using integrated data ***\n")
  
  StopWatchStart$ClusterEachSampleTypeCellsFromInteg$sample_type  <- Sys.time()

  options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
  seurat.object.each_sample_type <- FindNeighbors(object = seurat.object.each_sample_type, dims = PcaDimsUse)
  seurat.object.each_sample_type <- FindClusters(object = seurat.object.each_sample_type, reduction.type = "pca", resolution = Resolution, print.output = 0, save.SNN = T)
  ClustersThisSampleType <- unlist(x = strsplit(x = paste(sample_type, seurat.object.each_sample_type@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
  
  StopWatchEnd$ClusterEachSampleTypeCellsFromInteg$sample_type  <- Sys.time()
  
  ####################################
  ### Write out each sample type re-clustered cell clusters
  ####################################
  writeLines("\n*** Write out each sample type re-clustered cell clusters ***\n")

  StopWatchStart$WriteClustersEachSampleTypeCellClusterTables$sample_type <- Sys.time()
  
  write.table(paste(colnames(seurat.object.each_sample_type), ClustersThisSampleType, sep = "\t", collapse = "\n"),
              file = OutfileEachSampleTypeReClusters, row.names = F, col.names = F, quote = F, append = T)
  
  CellNames<-rownames(seurat.object.each_sample_type@meta.data)
  ClusterIdent <-seurat.object.each_sample_type@meta.data$seurat_clusters
  NumberOfClusters<-length(unique(ClusterIdent))
  OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", sample_type, "_NumbCellClusters", ".tsv", sep="")
  write(x=NumberOfClusters,file = OutfileNumbClusters)
  
  StopWatchEnd$WriteClustersEachSampleTypeCellClusterTables$sample_type <- Sys.time()
  
  ####################################
  ### Finding differentially expressed genes for each sample type re-clustered cell clusters
  ####################################
  writeLines("\n*** Finding differentially expressed genes for each sample type re-clustered cell clusters ***\n")
  
  StopWatchStart$FindDiffMarkers$sample_type  <- Sys.time()
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_sample_type@meta.data))
  
  seurat.object.each_sample_type.markers <- FindAllMarkers(object = seurat.object.each_sample_type, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
  
  write.table(data.frame("GENE"=rownames(seurat.object.each_sample_type.markers),seurat.object.each_sample_type.markers),paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", sample_type, "_MarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)
  
  StopWatchEnd$FindDiffMarkers$sample_type  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by each sample type re-clustered cells
  ####################################

  writeLines("\n*** Colour dimension reduction plots by each sample type re-clustered cells ***\n")
  
  for (dim_red_method in names(DimensionReductionMethods)) {

    StopWatchStart$UmapAndTsneColourByEachSampleTypeCellCluster$sample_type  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_sample_type, group.by = c("seurat_clusters"), combine = FALSE, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
    IntegratedUMAPPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", sample_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedUMAPPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
  
    StopWatchEnd$UmapAndTsneColourByEachSampleTypeCellCluster$sample_type <- Sys.time()
    
    ### To program layout() for more than 3 genes in multiple rows
    if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
      print("No selected genes for dimension reduction plots")
    }else{
      
      StopWatchEnd$DimRedPlotsBySampleTypeColuredByGenes$dim_red_method  <- Sys.time()
  
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
      
      ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
      ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
      seurat.object.each_sample_type.sa <- CreateSeuratObject(seurat.object.each_sample_type@assays$RNA@data)
      seurat.object.each_sample_type.sa@reductions$umap <- seurat.object.each_sample_type@reductions$umap
      seurat.object.each_sample_type.sa@reductions$tsne <- seurat.object.each_sample_type@reductions$tsne
      
      pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_", sample_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
      print(FeaturePlot(object = seurat.object.each_sample_type.sa, ncol = nColFeaturePlot, features = ListOfGenesForDimRedPlots, cols = c("lightgrey", "blue"),
                        reduction = dim_red_method, order = T, slot = "data", pt.size = 1.5, min.cutoff = "q0.1", max.cutoff = "q90"))
      dev.off()
      
      rm(seurat.object.each_sample_type.sa)
      
      StopWatchEnd$DimRedPlotsBySampleTypeColuredByGenes$dim_red_method  <- Sys.time()
    
    }
  }
}

####################################
### Load each sample type re-cluster assignments
####################################
writeLines("\n*** Load each sample type re-cluster assignments ***\n")

EachSampleTypeReclusteredCellClusters <- data.frame(read.table(OutfileEachSampleTypeReClusters, header = T, row.names = 1, sep = "\t"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleTypeReclusteredCellClusters, col.name = "EachSampleTypeReclusteredCellClusters")

# switch the identity class of all cells to reflect each sample_type clusters
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleTypeReclusteredCellClusters

####################################
### Get average expression for each sample type re-clustered clusters
####################################
writeLines("\n*** Get average expression for each sample type re-clustered clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "PerSampleType_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "PerSampleType_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "PerSampleType_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)


####################################
### Saving the R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {

  writeLines("\n*** Saving the R object ***\n")
  
  StopWatchStart$SaveRDS  <- Sys.time()
  
  OutfileRDS<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_integrated_object.rds", sep="")
  saveRDS(seurat.object.integrated, file = OutfileRDS)
  
  StopWatchEnd$SaveRDS  <- Sys.time()

}else{
  
  writeLines("\n*** Not saving R object ***\n")
  
}

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_UsedOptions.txt", sep="")
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
    ### Because the StopWatch[Start|End]$LoadScRNAseqData$dataset
    ### are using the word `dataset` itself as key instead of the dataset 'ID'
    for (dataset in rownames(InputsTable)) {
      if (regexpr("POSIXct", class(StopWatchStart[[stepToClock]][[dataset]]), ignore.case = T)[1] == 1) {
        TimeStart <- StopWatchStart[[stepToClock]][[dataset]]
        TimeEnd   <- StopWatchEnd[[stepToClock]][[dataset]]
        TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
        ReportTime<-c(paste(stepToClock, TimeDiff, dataset, sep = "\t", collapse = ""))
        write(file = OutfileCPUusage, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
      }
    }
  }
}

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")

writeLines(paste(Outdir,"/", ProgramOutdir, "/", sep="", collapse = ""))

outfiles_to_move <- list.files(Tempdir, pattern = PrefixOutfiles, full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""), to=paste(Outdir,"/", ProgramOutdir, "/", eachFile,sep=""), overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
writeLines(paste("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report", sep = "", collapse = ""))

quit()
