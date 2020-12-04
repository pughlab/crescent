####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/v3.2/sctransform_vignette.html (SCtransform normalization)
### https://satijalab.org/seurat/v3.2/integration.html (general integration)
### https://satijalab.org/seurat/v3.2/immune_alignment.html (control vs. treatment)
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads integrated datasets R object,  presumably from script `Runs_Seurat_v3_MultiDatasets_Integration.R'
### 2) Process integrated datasets as a whole, including:
###    - dimension reduction
###    - 'global' cell clustering
###    - UMAP/tSNE plots by global cell clusters, sample, sample type, requested genes and metadata
###    - average gene expression
### 3) For each dataset on its own uses global clusters and repeats analysis from step 2
### 4) For each dataset on its own re-clusters cells and repeats analysis from step 2
### 5) For each dataset type on its own uses global clusters and repeats analysis from step 2
### 6) For each dataset type on its own re-clusters cells and repeats analysis from step 2
### 7) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_Seurat_v3_MultiDatasets_PCA_Clustering_DimReduction.R -h'
### for help
####################################

####################################
### Tested with R v4.0.2
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(Seurat))       # (CRAN) tested with v3.2.1. To run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # (CRAN) to calculate the percentages of extra properties to be plotted
suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(gtools))       # (CRAN) to do alphanumeric sorting. Only needed if using `-w Y`.
suppressPackageStartupMessages(library(loomR))        # (GitHub mojaveazure/loomR) needed for fron-end display of data. Only needed if using `-w Y`.
suppressPackageStartupMessages(library(tidyr))        # (CRAN) to handle tibbles and data.frames
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_Seurat_v3_MultiDatasets_PCA_Clustering_DimReduction.R"
ProgramOutdir  <- "SEURAT"

####################################
### Get inputs from command line argumets
####################################
#####

option_list <- list(
  make_option(c("-i", "--inputs_list"), default="NA",
              help="Path/name to a <tab> delimited file with one dataset per row and the following 12 columns specifying filters/details for each dataset:
                1) unique dataset ID (e.g. 'd1')
                2) /path_to/dataset
                3) dataset type (e.g. 'control' or 'treatment') or use 'type' to skip using dataset types
                4) dataset format (either 'MTX', 'TSV' or 'HDF5')
                5) dataset minimum mitochondrial fraction (e.g. '0')
                6) dataset maximum mitochondrial fraction (e.g. '0.05' for NucSeq or '0.2' for whole cell scRNA-seq)
                7) dataset minimum ribosomal fraction (e.g. '0')
                8) dataset maximum ribosomal fraction (e.g. '0.75')
                9) dataset minimum number of genes (e.g. '50')
                10) dataset maximum number of genes (e.g. '8000')
                11) dataset minimum number of reads (e.g. '1')
                12) dataset maximum number of reads (e.g. '80000')

                Notes:
                (a) The order of the list of datasets in --inputs_list may influence the results, including number of clusters,
                t-SNE/UMAP and differentially expressed genes. List datasets better measured first.
                
                (b) middle dashes '-' are not allowed in columns 1 or 3, if ocurring they will be replaced by low dashes '_'
                
                (c) Column 4 indicates the dataset format. It must be in either:
                'MTX'  and column 2 must be the path/name to an MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
                       'MTX' files can be the outputs from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                       Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                'TSV'  and column 2 must be the path/name to a <tab> delimited *file* with genes in rows vs. barcodes in columns
                'HDF5' and column 2 must be the path/name of a *file* in hdf5 format (e.g. from Cell Ranger)

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-j", "--infile_r_object"), default="NA",
              help="Path/name to the Integration R/Seurat object

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-r", "--resolution"), default="1",
              help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of cell clusters
              
                Default = '1'"),
  #
  make_option(c("-v", "--clustering_inputs"), default="1",
              help="One or more of the following options. If more than one, pass the choice <comma> delimited, e.g. '1,2,3'
                '1' = globaly clusters the integrated datasets (mandatory in this script)
                '2' = re-clusters each datasets split from the integrated datasets. Needed for `-f 6`
                '3' = re-clusters each dataset type (column 3 from -inputs_list) from the integrated datasets. Needed for `-f 7`
                
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
  make_option(c("-c", "--infile_metadata"), default="NA",
              help="A <tab> delimited table of barcodes and discrete properties to colour the dimension reduction plots, like:
                Barcode                CellClass    InOtherDatasets
                d1_AAACCTGAGCGGCTTC-1   1            yes
                d2_AAACCTGAGTCGAGTG-1   1            no
                d3_AAACCTGCAAAGGAAG-1   2            yes
                Note: the barcode id must include the dataset ID
                
                Default = 'NA' (i.e. no metadata to be used for plots)"),
  #
  make_option(c("-g", "--infile_selected_genes"), default="NA",
              help="A path/name to a file with a list of gene identifiers (one-per-row) whose expression will be mapped into the dimension reduction plots
              
                Default = 'NA' (i.e. no genes to be plotted)"),
  #
  make_option(c("-m", "--apply_list_genes"), default="0",
              help="One or more of the following options. If more than one, pass the choice <comma> delimited, e.g. '1,2,3'
                '0' = no genes to be mapped into dimension reduction plots
                '1' = using cells from all integrated datasets to map genes requested by `-g`
                '2' = using cells from each dataset to map genes requested by `-g`
                '3' = using cells from each dataset type to map genes requested by `-g`
                
                Default = '0'"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file:
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                
                Default = '10'"),
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
                
                Default = 'N'"),
  #
  make_option(c("-x", "--minio_path"), default="NA",
              help="Used with the 'run_cwl' option above for mounting input data files in 'inputs_list' to CWL containers
              
                Default = 'NA'"),
  #
  make_option(c("-a", "--max_global_variables"), default="10000",
              help="Indicates maximum allowed total size (in bytes) of global variables identified. Used by library(future) to prevent too large exports
              
                Default = '10000' for 10000 MiB")
)

#####

opt <- parse_args(OptionParser(option_list=option_list))

InputsList              <- opt$inputs_list
InfileRobject           <- opt$infile_r_object
Resolution              <- as.numeric(opt$resolution)
ClusteringInputs        <- opt$clustering_inputs
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
InfileMetadata          <- opt$infile_metadata
InfileSelectedGenes     <- opt$infile_selected_genes
ApplySelectedGenes      <- opt$apply_list_genes
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl
MinioPath               <- opt$minio_path
MaxGlobalVariables      <- as.numeric(opt$max_global_variables)

#####
OneLineCommands <- paste0("\n", "One-line-commands used:", "\n", "`Rscript /path_to/", ThisScriptName)
for (optionInput in option_list) {
  OneLineCommands <- paste0(OneLineCommands, paste0(" ", optionInput@short_flag, " ", opt[optionInput@dest]))
}
OneLineCommands <- paste0(OneLineCommands, paste0("`\n"))
writeLines(OneLineCommands)

####################################
### Start stopwatches
####################################

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Using `-w Y` will make Tempdir, which takes the value of ProgramOutdir, and it will be the final out-directory
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) 
  
  FILE_TYPE_OUT_DIRECTORIES = c(
    "CRESCENT_CLOUD",
    "CRESCENT_CLOUD/frontend_normalized",
    "CRESCENT_CLOUD/frontend_coordinates",
    "CRESCENT_CLOUD/frontend_raw",
    "CRESCENT_CLOUD/frontend_groups",
    "AVERAGE_GENE_EXPRESSION_TABLES",
    "CELL_CLUSTER_IDENTITIES", 
    "CELL_FRACTIONS", 
    "DIMENSION_REDUCTION_COORDINATE_TABLES", 
    "DIMENSION_REDUCTION_PLOTS",
    "LOG_FILES",
    "QC_PLOTS",
    "R_OBJECTS",
    "PSEUDO_BULK",
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS",
    "UNFILTERED_DATA_MATRICES"
  )
  
}else{
  ## Using `Tempdir/DIRECTORY` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## 'DIRECTORY' is one of the directories specified at FILE_TYPE_OUT_DIRECTORIES
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  Tempdir        <- "~/temp" 
  #
  UserHomeDirectory <- Sys.getenv("HOME")[[1]]
  #
  Outdir<-gsub("^~/",paste0(UserHomeDirectory,"/"), Outdir)
  Tempdir<-gsub("^~/",paste0(UserHomeDirectory,"/"), Tempdir)
  Outdir<-gsub("/+", "/", Outdir, perl = T)
  Tempdir<-gsub("/+", "/", Tempdir, perl = T)
  Outdir<-gsub("/$", "", Outdir)
  Tempdir<-gsub("/$", "", Tempdir)
  #
  dir.create(file.path(Outdir, ProgramOutdir), recursive = T)
  dir.create(file.path(Tempdir), showWarnings = F, recursive = T)
  
  FILE_TYPE_OUT_DIRECTORIES = c(
    "AVERAGE_GENE_EXPRESSION_TABLES",
    "CELL_CLUSTER_IDENTITIES", 
    "CELL_FRACTIONS", 
    "DIMENSION_REDUCTION_COORDINATE_TABLES", 
    "DIMENSION_REDUCTION_PLOTS",
    "LOG_FILES",
    "QC_PLOTS",
    "R_OBJECTS",
    "PSEUDO_BULK",
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS",
    "UNFILTERED_DATA_MATRICES"
  )
}

sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(eachdir) {
  dir.create(file.path(paste0(Tempdir, "/", eachdir)), showWarnings = F, recursive = T)
})

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

StopWatchStart$ReportUsedOptions  <- Sys.time()

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_LogFiles_", "UsedOptions", ".txt")

TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=paste0("Run started: ", TimeOfRun, "\n"))

write(file = OutfileOptionsUsed, x=c("Commands used:"), append = T)
for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

write(file = OutfileOptionsUsed, x = OneLineCommands)

StopWatchEnd$ReportUsedOptions  <- Sys.time()

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_LogFiles_", "RSessionInfo", ".txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

####################################
### Define number of cores and RAM for parallelization
####################################

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (regexpr("^[0-9]+$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste0("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript ", ThisScriptName, " -h\n\n"))
}

if (NumbCoresToUse == 1) {
  plan(strategy = "sequential", workers = NumbCoresToUse)
  writeLines(paste0("\n", "*** Running: ", ThisScriptName, " in 'sequential' mode with ", NumbCoresToUse, " core ***", "\n"))
}else if (NumbCoresToUse > 1) {
  plan(strategy = "multicore", workers = NumbCoresToUse)
  writeLines(paste0("\n", "*** Running: ", ThisScriptName, " in 'multicore' mode with ", NumbCoresToUse, " cores ***", "\n"))
}else{
  stop(paste0("Unexpected --number_cores = ", NumbCoresToUse))
}

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4GiB of global variables identified (4000*1024^2), use: `options(future.globals.maxSize = 4000 * 1024^2)`
options(future.globals.maxSize = MaxGlobalVariables * 1024^2)

####################################
### Define default parameters
####################################
### Some of these default parameters are provided by Seurat developers, others are tailored empirically

RequestedClusteringInputs        = unlist(strsplit(ClusteringInputs, ","))
RequestedApplySelectedGenes      = unlist(strsplit(ApplySelectedGenes, ","))

DefaultParameters <- list(
  ### Parameters for QC plots
  CellPropertiesToQC = c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction"),
  
  ### Parameters for dimmension reduction plots
  MinNumberOfCellsToReducePerplexity = 150,
  ReducedPerplexity = 7,
  BaseSizeSinglePlotPdf  = 7,
  BaseSizeMultiplePlotPdfWidth  = 3.7,
  BaseSizeMultiplePlotPdfHeight = 3,
  MaxNumbLabelsPerRowInLegend   = 4,

  ### Parameters for datasets comparison
  AssaysForAverageGETables = c("RNA", "SCT", "integrated"),
  AssaysForPseudoBulk = c("RNA", "SCT")
)

### Assay types for plot and table outfiles
listAssaySuffixForOutfiles <- list(RNA="RNA", SCT="SCT", integrated="INT")

### Dimension reduction methods
DimensionReductionMethods<-list()
DimensionReductionMethods$umap$name <-"UMAP"
DimensionReductionMethods$tsne$name <-"TSNE"
DimensionReductionMethods$umap$run  <-as.function(RunUMAP)
DimensionReductionMethods$tsne$run  <-as.function(RunTSNE)

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("infiles_list", "infile_r_object", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
  }
}

####################################
### Check that optional parameters are compatible with each other
####################################
writeLines("\n*** Check that optional parameters are compatible with each other ***\n")

#### Metadata options

if (length(grep('^NA$', InfileMetadata, perl = T))) {
  writeLines(paste0("\n", "*** No metadata will be used *** \n"))
}else{
  CellPropertiesFromMetadata <- data.frame(read.table(InfileMetadata, header = T, row.names = 1, check.names = F))
  
  # This is because we are using Read10X(..., strip.suffix = T), which removes the last '-digit' from barcode ID's
  rownames(CellPropertiesFromMetadata)<-gsub(x =rownames(CellPropertiesFromMetadata), pattern = "-[0-9]+$", perl = T, replacement = "")
}

#### Gene mapping into dimension reduction plot options

if ((1 %in% RequestedApplySelectedGenes == T) | 
    (2 %in% RequestedApplySelectedGenes == T) | 
    (3 %in% RequestedApplySelectedGenes == T)
) {
  if (length(grep('^NA$', InfileSelectedGenes, perl = T))) {
    stop("ERROR!!! options `-m [1|2|3] and -g NA` are incompatible with each other")
  }else{
    writeLines(paste0("\n*** Get list of genes to map in dimension reduction plots ***\n"))
    
    ListOfGenesForDimRedPlots<-unique(readLines(con = InfileSelectedGenes, skipNul = T))
  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO LOAD DATASETS
################################################################################################################################################
################################################################################################################################################

####################################
### Load --inputs_list
####################################
writeLines("\n*** Load --inputs_list ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  if (regexpr("^NA$", MinioPath , ignore.case = T)[1] == 1) {
    
    InputsTable<-read.table(InputsList, header = F, row.names = 1, stringsAsFactors = F)
    colnames(InputsTable)<-c("PathToDataset","DatasetType","DatasetFormat","MinMitoFrac","MaxMitoFrac","MinRiboFrac","MaxRiboFrac","MinNGenes","MaxNGenes","MinNReads","MaxNReads")
    
  } else {
    MinioPaths <- as.list(strsplit(MinioPath, ",")[[1]])
    MinioDataPaths = data.frame(dataset_ID=rep(0, length(MinioPaths)), dataset_path=rep(0, length(MinioPaths)))
    
    for (i in seq_along(MinioPaths)) {
      MinioDataPaths[i, ] = c(basename(MinioPaths[[i]]), MinioPaths[[i]])
    }
    
    InputsTable0 <- read.table(InputsList, header = T, sep = ",", stringsAsFactors = F)
    
    MergedInputsTable <- merge(MinioDataPaths, InputsTable0, by="dataset_ID")
    MergeFilter <- c("name", "dataset_path", "dataset_type", "dataset_format", "mito_min", "mito_max", "ribo_min", "ribo_max", "ngenes_min", "ngenes_max", "nreads_min", "nreads_max")
    MergedInputsTableFiltered <- MergedInputsTable[MergeFilter]
    MergedInputsTableFilteredFinal <- MergedInputsTableFiltered[,-1]
    rownames(MergedInputsTableFilteredFinal) <- MergedInputsTableFiltered[,1]
    colnames(MergedInputsTableFilteredFinal) <-c("PathToDataset","DatasetType","DatasetFormat","MinMitoFrac","MaxMitoFrac","MinRiboFrac","MaxRiboFrac","MinNGenes","MaxNGenes","MinNReads","MaxNReads")
    
    InputsTable <- MergedInputsTableFilteredFinal
  }
} else {
  InputsList<-gsub("^~/",paste0(UserHomeDirectory,"/"), InputsList)
  InputsTable<-read.table(InputsList, header = F, row.names = 1, stringsAsFactors = F)
  colnames(InputsTable)<-c("PathToDataset","DatasetType","DatasetFormat","MinMitoFrac","MaxMitoFrac","MinRiboFrac","MaxRiboFrac","MinNGenes","MaxNGenes","MinNReads","MaxNReads")
}

##### Replace low dashes by dots in rownames(InputsTable) or DatasetType
rownames(InputsTable) <- gsub(x=rownames(InputsTable), pattern = "-", replacement = "_")
InputsTable[,"DatasetType"] <- gsub(x=InputsTable[,"DatasetType"], pattern = "-", replacement = "_")
NumberOfDatasets <- nrow(InputsTable)

list_DatasetToType <- list()
for (dataset in rownames(InputsTable)) {
  list_DatasetToType[[dataset]] <- InputsTable[dataset,"DatasetType"]
}

####################################
### Get's number of rows for the legend of dimension reduction plots
####################################
if (NumberOfDatasets <= DefaultParameters$MaxNumbLabelsPerRowInLegend) {
  NumbRowsForLegendPerDataset <- 1
}else{
  NumbRowsForLegendPerDataset <- round((NumberOfDatasets / DefaultParameters$MaxNumbLabelsPerRowInLegend), digits = 0)
}

NumberOfDatasetsTypes <- length(unique(InputsTable[,"DatasetType"]))
if (NumberOfDatasetsTypes <= DefaultParameters$MaxNumbLabelsPerRowInLegend) {
  NumbRowsForLegendPerDatasetType <- 1
}else{
  NumbRowsForLegendPerDatasetType <- round((NumberOfDatasetsTypes / DefaultParameters$MaxNumbLabelsPerRowInLegend), digits = 0)
}

####################################
### Load integrated datasets R object
####################################
writeLines("\n*** Load integrated datasets R object ***\n")

StopWatchStart$LoadRDSIntegratedDatasets  <- Sys.time()

seurat.object.integrated <- readRDS(InfileRobject)

StopWatchEnd$LoadRDSIntegratedDatasets  <- Sys.time()

####################################
### Outfiles for web app: save normalized count matrix as loom
####################################
if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Outfiles for web app: save normalized count matrix as loom ***\n")
  
  StopWatchStart$OutNormalizedTablesLoomFrontEnd  <- Sys.time()
  
  normalized_count_matrix <- as.matrix(seurat.object.integrated@assays[["RNA"]]@data)
  
  # all genes/features in matrix
  features_tsv <- data.frame(features = rownames(normalized_count_matrix))
  features_tsv_ordered <- as.data.frame(features_tsv[mixedorder(features_tsv$features),])
  write.table(features_tsv_ordered, file=paste0(Tempdir, "/CRESCENT_CLOUD/frontend_raw/","features.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # generating loom file of normalized count matrix
  loom_file <- paste0(Tempdir, "/CRESCENT_CLOUD/frontend_normalized/","normalized_counts.loom")
  create(loom_file, normalized_count_matrix, overwrite = T)
  
  StopWatchEnd$OutNormalizedTablesLoomFrontEnd <- Sys.time()
  
}

####################################
### Obtaining principal components
####################################
writeLines("\n*** Obtaining principal components ***\n")

StopWatchStart$RunPCA  <- Sys.time()

seurat.object.integrated <- RunPCA(seurat.object.integrated, verbose = F)

StopWatchEnd$RunPCA  <- Sys.time()

StopWatchStart$PCAPlots  <- Sys.time()

ForElbowPlot<-ElbowPlot(object = seurat.object.integrated, ndims = 50, reduction = "pca")
MaxYAxis<-as.integer(max(ForElbowPlot$data$stdev)+1)

OutfilePCElbowPlot <- paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", "PCElbowPlot", ".pdf")
pdf(file=OutfilePCElbowPlot, width = DefaultParameters$BaseSizeSinglePlotPdf, height = DefaultParameters$BaseSizeSinglePlotPdf)
print(ForElbowPlot
      + scale_x_continuous(breaks =  seq(from = 0, to = 50, by=5))
      + geom_vline(xintercept = seq(from = 0, to = 50, by=5), linetype='dotted', col="red")
      + scale_y_continuous(breaks =  seq(from = 0, to = MaxYAxis, by=0.5))
      + geom_hline(yintercept = seq(from = 0, to = MaxYAxis, by=0.5), linetype='dotted', col="red")
)
dev.off()

StopWatchEnd$PCAPlots  <- Sys.time()

####################################
### Run dimension reductions using integrated data
####################################
writeLines("\n*** Run dimension reductions using integrated data ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {
  ####################################
  ### Run non-linear dimension reductions using integrated data
  ####################################
  writeLines(paste0("\n*** Run ", DimensionReductionMethods[[dim_red_method]][["name"]], " ***\n"))
  
  StopWatchStart$DimensionReduction[[dim_red_method]]  <- Sys.time()
  
  ### NOTES:
  ### In RunTSNE: if the datasets is small user may get error:
  ### `Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large.`
  ### User can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
  ###
  ### Also using RunTSNE(..., check_duplicates = F) to skip cases where cells happen to have the same values after PCA reduction
  
  if (("tsne" %in% dim_red_method == T) & (length(colnames(seurat.object.integrated)) < DefaultParameters$MinNumberOfCellsToReducePerplexity)) {
    writeLines(paste0("\n*** Using reduced perplexity = ", DefaultParameters$ReducedPerplexity, " because found ",  length(colnames(seurat.object.integrated)), " cells", " ***\n"))
    seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, perplexity = DefaultParameters$ReducedPerplexity, check_duplicates = F)
  }else if ("tsne" %in% dim_red_method == T) {
    seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, check_duplicates = F)
  }else if ("umap" %in% dim_red_method == T) {
    seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, umap.method = "uwot")
  }
  
  StopWatchEnd$DimensionReduction[[dim_red_method]]  <- Sys.time()
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS ALL INTEGRATED DATASETS AS A WHOLE
################################################################################################################################################
################################################################################################################################################

####################################
### Globally cluster cells using integrated data
####################################
writeLines("\n*** Globally cluster cells using integrated data ***\n")

StopWatchStart$ClusterAllCells  <- Sys.time()

options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.integrated <- FindNeighbors(object = seurat.object.integrated, dims = PcaDimsUse)
seurat.object.integrated <- FindClusters(object = seurat.object.integrated, resolution = Resolution)

StopWatchEnd$ClusterAllCells  <- Sys.time()

StopWatchStart$AllCellClusterTables  <- Sys.time()

CellNames<-rownames(seurat.object.integrated@meta.data)
ClusterIdent <-seurat.object.integrated@meta.data$seurat_clusters
NumberOfClusters<-length(unique(ClusterIdent))

Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "CellClusters", ".tsv")
write.table(Headers, file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data), file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
OutfileNumbClusters<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellClusters", ".tsv")
write(x=NumberOfClusters, file = OutfileNumbClusters)
#
OutfileNumbCellsPerCluster<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerCluster", ".tsv")
Headers<-paste("Cluster", "Number_of_cells", sep="\t", collapse = "")
write.table(Headers, file = OutfileNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F)
write.table(table(ClusterIdent), file = OutfileNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
StopWatchEnd$AllCellClusterTables  <- Sys.time()

## Creates dataframe to merge for groups.tsv
if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Creates dataframe to merge for groups.tsv ***\n")
  
  StopWatchStart$SaveGlobalCellClusterIdentitiesFrontEnd  <- Sys.time()
  
  OutfileClustersDataframe  <- data.frame(NAME = rownames(seurat.object.integrated@meta.data),
                                          Seurat_Global_Clusters = seurat.object.integrated@meta.data$seurat_clusters,
                                          monochrome = 1,
                                          dataset_id = seurat.object.integrated@meta.data$dataset,
                                          dataset_type = as.character(list_DatasetToType[c(seurat.object.integrated@meta.data$dataset)])
  )
  
  StopWatchEnd$SaveGlobalCellClusterIdentitiesFrontEnd  <- Sys.time()
  
}

####################################
### Saving global cell cluster identities
####################################
writeLines("\n*** Saving global cell cluster identities ***\n")

StopWatchStart$SaveGlobalCellClusterIdentities  <- Sys.time()

seurat.object.integrated$GlobalCellClusterIdentities <- Idents(object = seurat.object.integrated)

StopWatchEnd$SaveGlobalCellClusterIdentities  <- Sys.time()

####################################
### Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset cluster (columns)
####################################
writeLines("\n*** Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset cluster (columns) ***\n")

StopWatchStart$SaveSctransformClustersInColumns  <- Sys.time()

seurat.object.integrated.list <- SplitObject(seurat.object.integrated, split.by = "dataset")

for (assay_expression in DefaultParameters$AssaysForPseudoBulk) {
  mat_for_correl_each_cluster.df <- data.frame(row.names = rownames(seurat.object.integrated.list[[1]]@assays[[assay_expression]]))
  for (DatasetId in rownames(InputsTable)) {
    if (exists(x = "seurat.object.each_dataset") == T) {
      rm(seurat.object.each_dataset)
    }
    seurat.object.each_dataset <- subset(x = seurat.object.integrated, subset = dataset == DatasetId)
    
    for (cluster_number in sort(unique(seurat.object.integrated$GlobalCellClusterIdentities))) {
      dataset_cluster <- (paste(DatasetId, cluster_number, sep = "_c", collapse = ""))
      res <- try(subset(x = seurat.object.each_dataset, subset = seurat_clusters == cluster_number), silent = TRUE)
      if (class(res) == "try-error") {
        # mat_for_correl_each_cluster.df[[dataset_cluster]] <- NA
      }else{
        if (exists(x = "seurat.object.each_dataset.each_cluster") == T) {
          rm(seurat.object.each_dataset.each_cluster)
        }
        seurat.object.each_dataset.each_cluster <- subset(x = seurat.object.each_dataset, subset = seurat_clusters == cluster_number)
        mat_for_correl_each_cluster.df[[dataset_cluster]] <- rowSums(as.matrix(seurat.object.each_dataset.each_cluster@assays[[assay_expression]][,]))
      }
    }
  }
  
  OutfilePseudoBulk <- paste0(Tempdir, "/PSEUDO_BULK/", PrefixOutfiles, ".", ProgramOutdir, "_PseudoBulk_", "EachDatasetCluster_", listAssaySuffixForOutfiles[[assay_expression]], ".tsv")
  Headers<-paste("PseudoBulk_EachDatasetCluster", paste(colnames(mat_for_correl_each_cluster.df), sep = "\t", collapse = "\t"), sep = "\t", collapse = "")
  write.table(Headers, file = OutfilePseudoBulk, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(mat_for_correl_each_cluster.df,  file = OutfilePseudoBulk, row.names = T, col.names = F, sep="\t", quote = F, append = T)
}

StopWatchEnd$SaveSctransformClustersInColumns  <- Sys.time()

####################################
### FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data to colour by: a) cell clusters, b) dataset, c) dataset type, d) selected genes and e) metadata
####################################
writeLines("\n*** FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data to colour by: a) cell clusters, b) dataset, c) dataset type, d) selected genes and e) metadata  ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {
  
  ####################################
  ### Colour dimension reduction plot by global cell clusters using integrated data
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by global cell clusters ***\n"))
  
  StopWatchStart$DimRedOPlotColourByCellCluster[[dim_red_method]]  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 10)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
  IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "GlobalClustering", "_ColourByCellClusters", ".pdf")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$DimRedOPlotColourByCellCluster[[dim_red_method]]  <- Sys.time()
  
  ####################################
  ### Write out dimension reduction plot coordinates using integrated data
  ####################################
  writeLines(paste0("\n*** Write out ", DimensionReductionMethods[[dim_red_method]][["name"]], " coordinates ***\n"))
  
  StopWatchStart$DimensionReductionWriteCoords[[dim_red_method]]  <- Sys.time()
  
  OutfileCoordinates<-paste0(Tempdir, "/DIMENSION_REDUCTION_COORDINATE_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_Coordinates", ".tsv")
  Headers<-paste("Barcode",paste(colnames(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings), sep="", collapse="\t"), sep="\t", collapse = "\t")
  write.table(Headers, file = OutfileCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$DimensionReductionWriteCoords[[dim_red_method]]  <- Sys.time()
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    ####################################
    ### Outfiles for web app: dimension reduction coordinates
    ####################################
    writeLines("\n*** Outfiles for web app: dimension reduction coordinates ***\n")
    
    StopWatchStart$DimensionReductionWriteCoordsFrontEnd[[dim_red_method]]  <- Sys.time()
    
    OutfileCoordinatesCWL<-paste0(Tempdir, "/CRESCENT_CLOUD/frontend_coordinates/",DimensionReductionMethods[[dim_red_method]][["name"]], "Coordinates.tsv")
    write.table(Headers, file = OutfileCoordinatesCWL, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinatesCWL,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
    
    StopWatchEnd$DimensionReductionWriteCoordsFrontEnd[[dim_red_method]]  <- Sys.time()
    
  }
  
  ####################################
  ### Colour dimension reduction plots by dataset using integrated data
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by dataset ***\n"))
  
  StopWatchStart$DimRedPlotsByDataset[[dim_red_method]]  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("dataset"), combine = F, reduction = dim_red_method)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = NumbRowsForLegendPerDataset, byrow = T, override.aes = list(size = 2))))
  IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "GlobalClustering", "_ColourByDataset", ".pdf")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  ### This is the same plot as above, but with only two dataset labels per row in the legend to make sure that all text is printed out
  plots <- DimPlot(seurat.object.integrated, group.by = c("dataset"), combine = F, reduction = dim_red_method)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(ncol = 2, byrow = T, override.aes = list(size = 1))))
  IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "GlobalClustering", "_ColourByDataset_FullLegend", ".pdf")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = (7 + (0.1 * NumberOfDatasets)))
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$DimRedPlotsByDataset[[dim_red_method]]  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by dataset type using integrated data
  ####################################
  
  if (NumberOfDatasetsTypes >= 1) {
    writeLines("\n*** Get dataset_type assignments replacing dataset ids using list_DatasetToType() ***\n")
    
    StopWatchStart$PrepareEachDatasetTypeForDimRedPlotByDatasetType <- Sys.time()
    
    seurat.object.integrated$dataset_type<-seurat.object.integrated$dataset
    for (dataset in rownames(InputsTable)) {
      seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
    }
    DatasetTypes <- unique(seurat.object.integrated$dataset_type)
    
    StopWatchEnd$PrepareEachDatasetTypeForDimRedPlotByDatasetType <- Sys.time()
    
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by dataset type ***\n"))
    
    Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset_type
    
    NumberOfDatasetTypes <- 0
    for (dataset_type in DatasetTypes) {
      NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
      print(NumberOfDatasetTypes)
      
      if (exists(x = "seurat.object.each_dataset_type") == T) {
        rm(seurat.object.each_dataset_type)
      }
      seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
      
      writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " coloured by dataset_type ***\n"))
      
      StopWatchStart$DimRedOPlotColourByDatasetType[[dim_red_method]]  <- Sys.time()
      
      plots <- DimPlot(seurat.object.integrated, group.by = c("dataset_type"), combine = F, reduction = dim_red_method, label = F, label.size = 10)
      plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = NumbRowsForLegendPerDatasetType, byrow = T, override.aes = list(size = 2))))
      IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "GlobalClustering", "_ColourByDatasetType", ".pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
      print(CombinePlots(plots))
      dev.off()
      
      StopWatchEnd$DimRedOPlotColourByDatasetType[[dim_red_method]]  <- Sys.time()
    }
  }
  
  ####################################
  ### Colour dimension reduction plots for all cells by selected genes (option -g)
  ####################################
  
  if (regexpr("^NA$", InfileSelectedGenes, ignore.case = T)[1] == 1) {
    print("No selected genes for dimension reduction plots")
    
  }else if (1 %in% RequestedApplySelectedGenes == T) {
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " for all cells by selected genes (option -g) ***\n"))
    
    StopWatchStart$DimRedPlotsColuredByGenes[[dim_red_method]]  <- Sys.time()
    
    ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
    ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
    seurat.object.integrated.sa <- CreateSeuratObject(seurat.object.integrated@assays$RNA@data)
    seurat.object.integrated.sa@reductions$umap <- seurat.object.integrated@reductions$umap
    seurat.object.integrated.sa@reductions$tsne <- seurat.object.integrated@reductions$tsne
    
    OutDirThisOption <- paste0(Tempdir, "/SELECTED_GENE_DIMENSION_REDUCTION_PLOTS/ALL_CELLS/", DimensionReductionMethods[[dim_red_method]][["name"]])
    dir.create(path = OutDirThisOption, recursive = T)
    sapply(ListOfGenesForDimRedPlots, FUN=function(eachGene) {
      if (eachGene %in% rownames(seurat.object.integrated.sa) == T)  {
        OutfilePathAndName <- paste0(OutDirThisOption, "/", eachGene, "_all_cells_", DimensionReductionMethods[[dim_red_method]][["name"]], ".pdf")
        pdf(file=OutfilePathAndName, width=DefaultParameters$BaseSizeMultiplePlotPdfWidth, height=DefaultParameters$BaseSizeMultiplePlotPdfHeight)
        print(FeaturePlot(object = seurat.object.integrated.sa, features = eachGene, cols = c("lightgrey", "blue"),
                          reduction = dim_red_method, order = T, slot = "data", pt.size = 0.3, min.cutoff = "q0.1", max.cutoff = "q90"))
        dev.off()
      }else{
        print(paste0("Missing ", eachGene))
      }
    })
    
    StopWatchEnd$DimRedPlotsColuredByGenes[[dim_red_method]]  <- Sys.time()
    
  }
  
  ####################################
  ### Colour dimension reduction plots by -infile_metadata using integrated data
  ####################################
  
  if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
    print("No metadata will be used for dimension reduction plots")
  }else{
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by -infile_metadata ***\n"))
    
    StopWatchStart$DimRedPlotsColuredByMetadata[[dim_red_method]]  <- Sys.time()
    
    # Note: needs data.frame(CellPropertiesFromMetadata) preloaded
    
    seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = CellPropertiesFromMetadata)
    
    # Generating outfile
    # Note DimPlot() takes the entire current device (pdf)
    # even if using layout(matrix(...)) or  par(mfrow=())
    # Thus each property plot is written to a separate page of a single *pdf outfile
    for (property in colnames(CellPropertiesFromMetadata)) {
      IntegratedDimRedPlotPdf <- paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "Metadata", "_", property, ".pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(CellPropertiesFromMetadata[,property])))))
      
      if ( (sum(CellPropertiesFromMetadata[,property] %in% 0:1 == T)) == (nrow(CellPropertiesFromMetadata)) ) { ## is binary
        CellsToHighlight <- rownames(CellPropertiesFromMetadata)[CellPropertiesFromMetadata[,property]==1]
        plots <- DimPlot(seurat.object.integrated, group.by = property, combine = F, reduction = dim_red_method, label = F, cells.highlight = CellsToHighlight)
        plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
        print(CombinePlots(plots))
        
      }else{
        plots <- DimPlot(seurat.object.integrated, group.by = property, combine = F, reduction = dim_red_method, label = T)
        plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top", legend.text.align = 0) + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
        print(CombinePlots(plots))
        
      }
      dev.off()
    }
    
    StopWatchEnd$DimRedPlotsColuredByMetadata[[dim_red_method]]  <- Sys.time()
  }
}

####################################
### Get average gene expression for each global cluster
####################################
writeLines("\n*** Get average gene expression for each global cluster ***\n")

StopWatchStart$AverageGeneExpressionGlobalClusters  <- Sys.time()

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$seurat_clusters

for (assay_expression in DefaultParameters$AssaysForAverageGETables) {
  cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F, assays = assay_expression)
  OutfileClusterAverages<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_AverageGeneExpression_", "GlobalClustering", "_AllDatasets_", listAssaySuffixForOutfiles[[assay_expression]], ".tsv")
  Headers<-paste("AVERAGE_GENE_EXPRESSION", paste("r", Resolution, "_C",  names(cluster.averages[[assay_expression]]), sep="", collapse="\t"), sep="\t", collapse = "\t")
  write.table(Headers, file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages[[assay_expression]]), file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
}

StopWatchEnd$AverageGeneExpressionGlobalClusters  <- Sys.time()

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET USING GLOBAL CELL CLUSTERS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET uses integrated/normalized counts and generates dimension reduction plots for each dataset coloured by: a) QC, b) global clusters, c) selected genes and d) metadata
####################################
writeLines("\n*** FOR EACH DATASET uses integrated/normalized counts and generates dimension reduction plots for each dataset coloured by: a) QC, b) global clusters, c) selected genes and d) metadata ***\n")

####################################
### Colour dimension reduction plots for each dataset by QC (nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction)
####################################
writeLines("\n*** Colour dimension reduction plots for each dataset by QC (nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction) ***\n")

StopWatchStart$DimRedPlotsQC  <- Sys.time()

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset

for (dataset in rownames(InputsTable)) {
  
  if (exists(x = "seurat.object.each_dataset") == T) {
    rm(seurat.object.each_dataset)
  }
  seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)

  for (dim_red_method in names(DimensionReductionMethods)) {
    OutfileDimRedPlotQc <- paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plots", "_", dataset, ".pdf")
    pdf(file=OutfileDimRedPlotQc, width = DefaultParameters$BaseSizeSinglePlotPdf * 1.5, height = DefaultParameters$BaseSizeSinglePlotPdf * 1.5)
    print(FeaturePlot(object = seurat.object.each_dataset, label = F, order = T, features = DefaultParameters$CellPropertiesToQC, cols = c("lightgrey", "blue"), reduction = dim_red_method, ncol = 2, pt.size = 1.5))
    dev.off()
  }
}

StopWatchEnd$DimRedPlotsQC  <- Sys.time()

####################################
### Merge dataset_id and global cluster identities to make dataset-specific identities
####################################
writeLines("\n*** Merge dataset_id and global cluster identities to make dataset-specific identities  ***\n")

StopWatchStart$MergeDatasetAndClusterIds  <- Sys.time()

EachDatasetGlobalCellClusters <- unlist(x = strsplit(x = paste(seurat.object.integrated$dataset, seurat.object.integrated$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetGlobalCellClusters, col.name = "EachDatasetGlobalCellClusters")
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetGlobalCellClusters

StopWatchEnd$MergeDatasetAndClusterIds  <- Sys.time()

####################################
### Write out number and fraction of cells per cluster, per dataset
####################################
writeLines("\n*** Write out number and fraction of cells per cluster, per dataset  ***\n")

StopWatchStart$SaveFractionsOfCellsPerClusterPerDataset  <- Sys.time()

OutfileNumbCellsPerClusterPerDataset<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerClusterPerDataset", ".tsv")
OutfileFracCellsPerClusterPerDataset<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "FracCellsPerClusterPerDataset", ".tsv")

seurat.object.integrated_EachDatasetGlobalCellClusters.df <- as.data.frame(table(seurat.object.integrated$EachDatasetGlobalCellClusters))
rownames(seurat.object.integrated_EachDatasetGlobalCellClusters.df) <- seurat.object.integrated_EachDatasetGlobalCellClusters.df[,1]

freq_mat <- sapply(rownames(InputsTable), function(dataset) {
  freq_clust <- sapply(sort(unique(seurat.object.integrated$seurat_clusters)), function(cluster) {
    dataset_cluster <- paste0(dataset, "_c", cluster)
    if ((dataset_cluster %in% row.names(seurat.object.integrated_EachDatasetGlobalCellClusters.df)) == TRUE)  {
      freq <- seurat.object.integrated_EachDatasetGlobalCellClusters.df[dataset_cluster,"Freq"]
    }else{
      freq <- 0
    }
    freq ## returns freq for the second sapply loop
  })
  freq_clust ## returns freq_clust for the first sapply loop
})
rownames(freq_mat) <- sort(unique(seurat.object.integrated$seurat_clusters))
frac_mat <- t(round(t(freq_mat) / colSums(freq_mat), 4)) ### note it's necessary to t(freq_ma) otherwise the colSums() are used incorrectly

Headers<-paste("cluster", paste(rownames(InputsTable) , sep="", collapse = "\t"), sep = "\t", collapse = "\t")

write.table(Headers, file = OutfileNumbCellsPerClusterPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(freq_mat, file = OutfileNumbCellsPerClusterPerDataset, row.names = T, col.names = F, sep="\t", quote = F, append = T)

write.table(Headers, file = OutfileFracCellsPerClusterPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(frac_mat , file = OutfileFracCellsPerClusterPerDataset, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$SaveFractionsOfCellsPerClusterPerDataset  <- Sys.time()

####################################
### Colour dimension reduction plots for each dataset based on global clusters, selected genes and metadata
####################################

writeLines("\n*** Colour dimension reduction plots for each dataset based on global clusters, selected genes and metadata ***\n")

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_dataset") == T) {
    rm(seurat.object.each_dataset)
  }
  seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_dataset)
  
  for (dim_red_method in names(DimensionReductionMethods)) {
    
    ####################################
    ### Colour dimension reduction plots for each dataset by global cell clusters
    ####################################
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for each dataset by global cell clusters ***\n"))
    
    StopWatchStart$DimRedPlotsByDatasetIntegratedCellClusters[[dataset]][[dim_red_method]]  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_dataset, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "GlobalClustering", "_", dataset, ".pdf")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
    
    StopWatchEnd$DimRedPlotsByDatasetIntegratedCellClusters[[dataset]][[dim_red_method]]  <- Sys.time()
    
    ####################################
    ### Colour dimension reduction plots for each dataset by selected genes (option -g)
    ####################################
    
    if (regexpr("^NA$", InfileSelectedGenes, ignore.case = T)[1] == 1) {
      print("No selected genes for each dataset dimension reduction plots")
      
    }else if (2 %in% RequestedApplySelectedGenes == T) {
      writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for each dataset by selected genes (option -g) ***\n"))
      
      StopWatchStart$DimRedPlotsColuredByGenesEachDataset[[dim_red_method]]  <- Sys.time()
      
      ### Making a new Seurat object `seurat.object.each_dataset.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
      ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
      seurat.object.each_dataset.sa <- CreateSeuratObject(seurat.object.each_dataset@assays$RNA@data)
      seurat.object.each_dataset.sa@reductions$umap <- seurat.object.each_dataset@reductions$umap
      seurat.object.each_dataset.sa@reductions$tsne <- seurat.object.each_dataset@reductions$tsne
      
      OutDirThisOption <- paste0(Tempdir, "/SELECTED_GENE_DIMENSION_REDUCTION_PLOTS/DATASETS/", DimensionReductionMethods[[dim_red_method]][["name"]])
      dir.create(path = OutDirThisOption, recursive = T)
      sapply(ListOfGenesForDimRedPlots, FUN=function(eachGene) {
        if (eachGene %in% rownames(seurat.object.each_dataset.sa) == T)  {
          OutfilePathAndName <- paste0(OutDirThisOption, "/", eachGene, "_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], ".pdf")
          pdf(file=OutfilePathAndName, width=DefaultParameters$BaseSizeMultiplePlotPdfWidth, height=DefaultParameters$BaseSizeMultiplePlotPdfHeight)
          print(FeaturePlot(object = seurat.object.each_dataset.sa, features = eachGene, cols = c("lightgrey", "blue"),
                            reduction = dim_red_method, order = T, slot = "data", pt.size = 0.3, min.cutoff = "q0.1", max.cutoff = "q90"))
          dev.off()
        }else{
          print(paste0("Missing ", eachGene))
        }
      })
      
      StopWatchEnd$DimRedPlotsColuredByGenesEachDataset[[dim_red_method]]  <- Sys.time()
      
      rm(seurat.object.each_dataset.sa)
      
    }
    
    ####################################
    ### Colour dimension reduction plots for each dataset by -infile_metadata using integrated data
    ####################################
    
    if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
      print("No metadata will be used for dimension reduction plots")
    }else{
      writeLines("\n*** Colour dimension reduction plots for each dataset by -infile_metadata using integrated data ***\n")
      
      # Note: needs data.frame(CellPropertiesFromMetadata) preloaded
      
      seurat.object.each_dataset <- AddMetaData(object = seurat.object.each_dataset, metadata = CellPropertiesFromMetadata)
      
      StopWatchStart$DimRedPlotsByMetadata[[dataset]][[dim_red_method]] <- Sys.time()
      
      # Generating outfile
      # Note DimPlot() takes the entire current device (pdf)
      # even if using layout(matrix(...)) or  par(mfrow=())
      # Thus each property plot is written to a separate page of a single *pdf outfile
      for (property in colnames(CellPropertiesFromMetadata)) {
        
        if (sum(is.na(seurat.object.each_dataset[[property]]==T)) == ncol(seurat.object.each_dataset)) {
          
          print(paste0("No extra barcode-attributes '", property, "' found for dataset ", dataset))
          
        }else{
          
          IntegratedDimRedPlotPdf <- paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "Metadata", "_", property, "_", dataset, ".pdf")
          pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(CellPropertiesFromMetadata[,property])))))
          
          if ( (sum(CellPropertiesFromMetadata[,property] %in% 0:1 == T)) == (nrow(CellPropertiesFromMetadata)) ) { ## is binary
            CellsToHighlight <- rownames(CellPropertiesFromMetadata)[CellPropertiesFromMetadata[,property]==1]
            plots <- DimPlot(seurat.object.each_dataset, group.by = property, combine = F, reduction = dim_red_method, label = F, cells.highlight = CellsToHighlight)
            plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
            print(CombinePlots(plots))
            
          }else{
            plots <- DimPlot(seurat.object.each_dataset, group.by = property, combine = F, reduction = dim_red_method, label = T)
            plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top", legend.text.align = 0) + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
            print(CombinePlots(plots))
            
          }
          
          dev.off()
        }
      }
      StopWatchEnd$DimRedPlotsByMetadata[[dataset]][[dim_red_method]]  <- Sys.time()
    }
    
  }
}

####################################
### Get average gene expression for each dataset based on global clusters
####################################
writeLines("\n*** Get average gene expression for each dataset based on global clusters ***\n")

StopWatchStart$AverageGeneExpressionEachDatasetGlobalClustersPerDataset  <- Sys.time()

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$EachDatasetGlobalCellClusters

for (assay_expression in DefaultParameters$AssaysForAverageGETables) {
  cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F, assays = assay_expression)
  OutfileClusterAverages<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_AverageGeneExpression_", "GlobalClustering", "_EachDataset_", listAssaySuffixForOutfiles[[assay_expression]], ".tsv")
  Headers<-paste("AVERAGE_GENE_EXPRESSION", paste("r", Resolution, "_C",  names(cluster.averages[[assay_expression]]), sep="", collapse="\t"), sep="\t", collapse = "\t")
  write.table(Headers, file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages[[assay_expression]]), file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
}

StopWatchEnd$AverageGeneExpressionEachDatasetGlobalClustersPerDataset  <- Sys.time()

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET USING RE-CLUSTERED CELLS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET uses integrated counts, RE-CLUSTER cells and colour dimension reduction plots
####################################

if (2 %in% RequestedClusteringInputs == T) {
  
  writeLines("\n*** FOR EACH DATASET uses integrated counts, RE-CLUSTER cells and colour dimension reduction plots ***\n")
  
  ####################################
  ### Prepares outfile headers
  ####################################
  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset
  
  # Headers for dataset-specific cell clusters
  Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
  OutfileEachDatasetClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "CellClusters", ".tsv")
  write.table(Headers, file = OutfileEachDatasetClusters, row.names = F, col.names = F, sep="\t", quote = F)
  #
  Headers<-paste("Dataset", "Number_of_clusters", sep="\t", collapse = "")
  OutfileEachDatasetNumbClusters<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "NumbCellClusters", ".tsv")
  write(x=NumberOfClusters, file = OutfileNumbClusters)
  #
  Headers<-paste("Cluster", "Number_of_cells", sep="\t", collapse = "")
  OutfileEachDatasetNumbCellsPerCluster<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "NumbCellsPerCluster", ".tsv")
  write.table(Headers, file = OutfileEachDatasetNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F)
  
  ####################################
  ### Loops each dataset
  ####################################
  NumberOfDatasets <- 0
  for (dataset in rownames(InputsTable)) {
    NumberOfDatasets <- NumberOfDatasets + 1
    print(NumberOfDatasets)
    
    if (exists(x = "seurat.object.each_dataset") == T) {
      rm(seurat.object.each_dataset)
    }
    seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
    
    ####################################
    ### Re-cluster each dataset cells using integrated data
    ####################################
    StopWatchStart$ClusterEachDatasetCellsFromInteg[[dataset]]  <- Sys.time()
    
    options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
    seurat.object.each_dataset <- FindNeighbors(object = seurat.object.each_dataset, dims = PcaDimsUse)
    seurat.object.each_dataset <- FindClusters(object = seurat.object.each_dataset, resolution = Resolution)
    ClustersThisDataset <- unlist(x = strsplit(x = paste(dataset, seurat.object.each_dataset@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
    
    StopWatchEnd$ClusterEachDatasetCellsFromInteg[[dataset]]  <- Sys.time()
    
    ####################################
    ### Write out each dataset cell re-clusters
    ####################################
    writeLines("\n*** Write out each dataset cell re-clusters ***\n")
    
    StopWatchStart$WriteClustersEachDatasetCellClusterTables[[dataset]] <- Sys.time()
    
    ### Cell cluster identities
    write.table(paste(colnames(seurat.object.each_dataset), ClustersThisDataset, sep = "\t", collapse = "\n"),
                file = OutfileEachDatasetClusters, row.names = F, col.names = F, quote = F, append = T)
    
    ### Number of clusters per dataset
    CellNames<-rownames(seurat.object.each_dataset@meta.data)
    ClusterIdent <-seurat.object.each_dataset@meta.data$seurat_clusters
    NumberOfClusters<-length(unique(ClusterIdent))
    write(x=paste(dataset, NumberOfClusters, sep = "\t", collapse = ""), file = OutfileEachDatasetNumbClusters, append = T)
    
    ### Number of cells per cluster, per dataset
    NumberOfCellPerClusterPerDataset<-table(ClusterIdent)
    rownames(NumberOfCellPerClusterPerDataset) <- paste(dataset, rownames(NumberOfCellPerClusterPerDataset), sep = "_c")
    write.table(NumberOfCellPerClusterPerDataset, file = OutfileEachDatasetNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    
    StopWatchEnd$WriteClustersEachDatasetCellClusterTables[[dataset]] <- Sys.time()
    
    ####################################
    ### Colour dimension reduction plots for each dataset by re-clustered cells
    ####################################
    writeLines("\n*** Colour dimension reduction plots for each dataset by re-clustered cells ***\n")
    
    for (dim_red_method in names(DimensionReductionMethods)) {
      
      StopWatchStart$DimRedPlotsByEachDatasetReclusteredCellCluster[[dataset]][[dim_red_method]]  <- Sys.time()
      
      plots <- DimPlot(seurat.object.each_dataset, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
      plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
      IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "EachDatasetReclustered_", dataset, ".pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
      print(CombinePlots(plots))
      dev.off()
      
      StopWatchEnd$DimRedPlotsByEachDatasetReclusteredCellCluster[[dataset]][[dim_red_method]] <- Sys.time()
      
    }
  }
  
  ####################################
  ### Load each dataset re-cluster assignments
  ####################################
  writeLines("\n*** Load each dataset re-cluster assignments ***\n")
  
  EachDatasetCellReClusters <- data.frame(read.table(OutfileEachDatasetClusters, header = T, row.names = 1, sep = "\t"))
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetCellReClusters, col.name = "EachDatasetCellReClusters")
  Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetCellReClusters
  
  ####################################
  ### Get average gene expression for each dataset re-clustered clusters
  ####################################
  writeLines("\n*** Get average gene expression for each dataset re-clustered clusters ***\n")
  
  StopWatchStart$AverageGeneExpressionEachDatasetReclustered[[dataset]]  <- Sys.time()
  
  for (assay_expression in DefaultParameters$AssaysForAverageGETables) {
    cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F, assays = assay_expression)
    
    OutfileClusterAverages<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_AverageGeneExpression_", "EachDatasetReclustered", "_", listAssaySuffixForOutfiles[[assay_expression]], ".tsv")
    Headers<-paste("AVERAGE_GENE_EXPRESSION", paste("r", Resolution, "_C",  names(cluster.averages[[assay_expression]]), sep="", collapse="\t"), sep="\t", collapse = "\t")
    write.table(Headers, file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(cluster.averages[[assay_expression]]), file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  }
  
  StopWatchEnd$AverageGeneExpressionEachDatasetReclustered[[dataset]]  <- Sys.time()
  
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET TYPE USING GLOBAL CELL CLUSTERS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET TYPE uses integrated/normalized counts and generates dimension reduction plots coloured by: a) dataset type, b) global clusters, c) selected genes and d) metadata
####################################
if (NumberOfDatasetsTypes >= 1) {
  writeLines("\n*** FOR EACH DATASET TYPE uses integrated/normalized counts and generates dimension reduction plots coloured by: a) dataset type, b) global clusters, c) selected genes and d) metadata ***\n")
  
  ####################################
  ### Prepares each dataset type object
  ####################################
  
  StopWatchStart$PrepareEachDatasetTypeForDimRedPlotsByVariousAttributes <- Sys.time()
  
  seurat.object.integrated$dataset_type<-seurat.object.integrated$dataset
  for (dataset in rownames(InputsTable)) {
    seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
  }
  DatasetTypes <- unique(seurat.object.integrated$dataset_type)
  
  StopWatchEnd$PrepareEachDatasetTypeForDimRedPlotsByVariousAttributes <- Sys.time()
  
  ####################################
  ### Write out each dataset type global cluster assignments
  ####################################
  writeLines("\n*** Write out each dataset type global cluster assignments ***\n")
  
  StopWatchStart$WriteOutDatasetTypeClusterAssignments  <- Sys.time()
  
  # Headers for dataset-type-specific cell clusters
  Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
  OutfileEachDatasetTypeGlobalClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "EachDatasetType_CellClusters", ".tsv")
  GlobalClustersByDatasetType<- unlist(x = strsplit(x = paste(seurat.object.integrated@meta.data$dataset_type, seurat.object.integrated@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
  write.table(Headers, file = OutfileEachDatasetTypeGlobalClusters, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(paste(colnames(seurat.object.integrated), GlobalClustersByDatasetType, sep = "\t", collapse = "\n"),
              file = OutfileEachDatasetTypeGlobalClusters, row.names = F, col.names = F, quote = F, append = T)
  
  OutfileNumbCellsPerClusterPerDatasetType<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerClusterPerDatasetType", ".tsv")
  Headers<-paste("Cluster", "Number_of_cells", sep="\t", collapse = "")
  write.table(Headers, file = OutfileNumbCellsPerClusterPerDatasetType, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(table(GlobalClustersByDatasetType), file = OutfileNumbCellsPerClusterPerDatasetType, row.names = F, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$WriteOutDatasetTypeClusterAssignments  <- Sys.time()
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    ####################################
    ### Outfiles for web app: each dataset type global cluster assignments
    ####################################
    writeLines("\n*** Outfiles for web app: each dataset type global cluster assignments ***\n")
    
    StopWatchStart$OutDatasetTypeClusterAssignmentsFrontEnd  <- Sys.time()
    
    OutfileEachDatasetTypeGlobalClustersDataframe  <- data.frame(NAME = rownames(seurat.object.integrated@meta.data), Seurat_Datasets_Global_Clusters = GlobalClustersByDatasetType)
    
    groupsMergedDataframe <- merge(OutfileClustersDataframe, OutfileEachDatasetTypeGlobalClustersDataframe, by="NAME")
    groupsMergedDataframeString <- sapply(groupsMergedDataframe, as.character)
    groupsMergedDataframeStringTYPE <- rbind(data.frame(NAME = "TYPE", Seurat_Global_Clusters = "group", monochrome = "group", dataset_id = "group", dataset_type = "group", Seurat_Datasets_Global_Clusters = "group"), groupsMergedDataframeString)
    colnames(groupsMergedDataframeStringTYPE) <- c("NAME", paste0("Seurat_Global_Clusters_Resolution_", Resolution), 
                                                   "monochrome",
                                                   "dataset_id",
                                                   "dataset_type",
                                                   paste0("Seurat_Datasets_Global_Clusters_Resolution_", Resolution))
    
    OutfileClustersMerged<-paste0(Tempdir, "/CRESCENT_CLOUD/frontend_groups/","groups.tsv")
    write.table(data.frame(groupsMergedDataframeStringTYPE), file = OutfileClustersMerged, row.names = F, col.names = T, sep="\t", quote = F, append = T)
    
    StopWatchEnd$OutDatasetTypeClusterAssignmentsFrontEnd  <- Sys.time()
    
  }
  
  ####################################
  ### Loops each dataset type
  ####################################
  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset_type
  
  NumberOfDatasetTypes <- 0
  for (dataset_type in DatasetTypes) {
    NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
    print(NumberOfDatasetTypes)
    
    if (exists(x = "seurat.object.each_dataset_type") == T) {
      rm(seurat.object.each_dataset_type)
    }
    seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
    
    for (dim_red_method in names(DimensionReductionMethods)) {
      
      ####################################
      ### Colour dimension reduction plots for each dataset type by global clusters
      ####################################
      
      writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for ", dataset_type, " using global cell clusters ***\n"))
      
      StopWatchStart$DimRedOPlotDatasetTypeColourByCellCluster[[dataset_type]][[dim_red_method]]  <- Sys.time()
      
      plots <- DimPlot(seurat.object.each_dataset_type, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 10)
      plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
      IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "GlobalClustering", "_", dataset_type, ".pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
      print(CombinePlots(plots))
      dev.off()
      
      StopWatchEnd$DimRedOPlotDatasetTypeColourByCellCluster[[dataset_type]][[dim_red_method]]  <- Sys.time()
      
      ####################################
      ### Colour dimension reduction plots for each dataset type by selected genes (option -g)
      ####################################
      
      if (regexpr("^NA$", InfileSelectedGenes, ignore.case = T)[1] == 1) {
        print("No selected genes for each dataset type dimension reduction plots")
        
      }else if (3 %in% RequestedApplySelectedGenes == T) {
        writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for each dataset type by selected genes (option -g) ***\n"))
        
        StopWatchStart$DimRedPlotsColuredByGenesEachDatasetType[[dim_red_method]]  <- Sys.time()
        
        ### Making a new Seurat object `seurat.object.each_dataset_type.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
        ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
        seurat.object.each_dataset_type.sa <- CreateSeuratObject(seurat.object.each_dataset_type@assays$RNA@data)
        seurat.object.each_dataset_type.sa@reductions$umap <- seurat.object.each_dataset_type@reductions$umap
        seurat.object.each_dataset_type.sa@reductions$tsne <- seurat.object.each_dataset_type@reductions$tsne
        
        OutDirThisOption <- paste0(Tempdir, "/SELECTED_GENE_DIMENSION_REDUCTION_PLOTS/DATASET_TYPES/", DimensionReductionMethods[[dim_red_method]][["name"]])
        dir.create(path = OutDirThisOption, recursive = T)
        sapply(ListOfGenesForDimRedPlots, FUN=function(eachGene) {
          if (eachGene %in% rownames(seurat.object.each_dataset_type.sa) == T)  {
            OutfilePathAndName <- paste0(OutDirThisOption, "/", eachGene, "_", dataset_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], ".pdf")
            pdf(file=OutfilePathAndName, width=DefaultParameters$BaseSizeMultiplePlotPdfWidth, height=DefaultParameters$BaseSizeMultiplePlotPdfHeight)
            print(FeaturePlot(object = seurat.object.each_dataset_type.sa, features = eachGene, cols = c("lightgrey", "blue"),
                              reduction = dim_red_method, order = T, slot = "data", pt.size = 0.3, min.cutoff = "q0.1", max.cutoff = "q90"))
            dev.off()
          }else{
            print(paste0("Missing ", eachGene))
          }
        })
        
        StopWatchEnd$DimRedPlotsColuredByGenesEachDatasetType[[dim_red_method]]  <- Sys.time()
        
        rm(seurat.object.each_dataset_type.sa)
      }
      
      ####################################
      ### Colour dimension reduction plots for each dataset type by -infile_metadata using integrated data
      ####################################
      
      if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
        print("No metadata will be used for dimension reduction plots")
      }else{
        writeLines("\n*** Colour dimension reduction plots for each dataset type by -infile_metadata using integrated data ***\n")
        
        # Note: needs data.frame(CellPropertiesFromMetadata) preloaded
        
        seurat.object.each_dataset_type <- AddMetaData(object = seurat.object.each_dataset_type, metadata = CellPropertiesFromMetadata)
        
        StopWatchStart$DimRedPlotsByMetadata[[dataset_type]][[dim_red_method]]  <- Sys.time()
        
        # Generating outfile
        # Note DimPlot() takes the entire current device (pdf)
        # even if using layout(matrix(...)) or  par(mfrow=())
        # Thus each property plot is written to a separate page of a single *pdf outfile
        for (property in colnames(CellPropertiesFromMetadata)) {
          
          if (sum(is.na(seurat.object.each_dataset_type[[property]]==T)) == ncol(seurat.object.each_dataset_type)) {
            
            print(paste0("No extra barcode-attributes '", property, "' found for dataset_type ", dataset_type))
            
          }else{
            
            IntegratedDimRedPlotPdf <- paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "Metadata", "_", property, "_", dataset_type, ".pdf")
            pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(CellPropertiesFromMetadata[,property])))))
            
            if ( (sum(CellPropertiesFromMetadata[,property] %in% 0:1 == T)) == (nrow(CellPropertiesFromMetadata)) ) { ## is binary
              CellsToHighlight <- rownames(CellPropertiesFromMetadata)[CellPropertiesFromMetadata[,property]==1]
              plots <- DimPlot(seurat.object.each_dataset_type, group.by = property, combine = F, reduction = dim_red_method, label = F, cells.highlight = CellsToHighlight)
              plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
              print(CombinePlots(plots))
            }else{
              plots <- DimPlot(seurat.object.each_dataset_type, group.by = property, combine = F, reduction = dim_red_method, label = T)
              plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top", legend.text.align = 0) + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
              print(CombinePlots(plots))
            }
            
            dev.off()
          }
        }
        StopWatchEnd$DimRedPlotsByMetadata[[dataset_type]][[dim_red_method]]  <- Sys.time()
      }
    }
  }
  
  ####################################
  ### Load each dataset type global cluster assignments
  ####################################
  writeLines("\n*** Load each dataset type global cluster assignments ***\n")
  
  StopWatchStart$LoadOutDatasetTypeClusterAssignments  <- Sys.time()
  
  EachDatasetTypeGlobalCellClusters <- data.frame(read.table(OutfileEachDatasetTypeGlobalClusters, header = T, row.names = 1, sep = "\t"))
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetTypeGlobalCellClusters, col.name = "EachDatasetTypeGlobalCellClusters")
  
  StopWatchEnd$LoadOutDatasetTypeClusterAssignments  <- Sys.time()
  
  ####################################
  ### Get average gene expression for each dataset type using global clusters
  ####################################
  writeLines("\n*** Get average gene expression for each dataset type using global clusters ***\n")
  
  StopWatchStart$AverageGeneExpressionEachDatasetTypeGlobalClusters  <- Sys.time()
  
  for (assay_expression in DefaultParameters$AssaysForAverageGETables) {
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetTypeGlobalCellClusters
    
    cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F, assays = assay_expression)
    
    OutfileClusterAverages<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_AverageGeneExpression_", "GlobalClustering", "_EachDatasetType_", listAssaySuffixForOutfiles[[assay_expression]], ".tsv")
    Headers<-paste("AVERAGE_GENE_EXPRESSION", paste("r", Resolution, "_C",  names(cluster.averages[[assay_expression]]), sep="", collapse="\t"), sep="\t", collapse = "\t")
    write.table(Headers, file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(cluster.averages[[assay_expression]]), file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
    
  }
  
  StopWatchEnd$AverageGeneExpressionEachDatasetTypeGlobalClusters  <- Sys.time()
  
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET TYPE USING RE-CLUSTERED CELLS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET TYPE uses integrated counts, RE-CLUSTER cells and colour dimension reduction plots
####################################

if (NumberOfDatasetsTypes >= 1) {
  writeLines("\n*** FOR EACH DATASET TYPE uses integrated counts, RE-CLUSTER cells and colour dimension reduction plots ***\n")
  
  ####################################
  ### Prepares each dataset type object
  ####################################
  
  StopWatchStart$PrepareEachDatasetTypeForRecusterAndDimRedPlots <- Sys.time()
  
  seurat.object.integrated$dataset_type<-seurat.object.integrated$dataset
  for (dataset in rownames(InputsTable)) {
    seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
  }
  DatasetTypes <- unique(seurat.object.integrated$dataset_type)
  
  StopWatchEnd$PrepareEachDatasetTypeForRecusterAndDimRedPlots <- Sys.time()
  
  ####################################
  ### Reclusters each dataset type cells
  ####################################
  if (3 %in% RequestedClusteringInputs == T) {
    
    ####################################
    ### Prepares outfile headers
    ####################################
    OutfileEachDatasetTypeCellReClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", "CellClusters", ".tsv")
    Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
    write.table(Headers, file = OutfileEachDatasetTypeCellReClusters, row.names = F, col.names = F, sep="\t", quote = F)
    
    ####################################
    ### Loops each dataset type
    ####################################
    NumberOfDatasetTypes <- 0
    for (dataset_type in DatasetTypes) {
      NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
      print(NumberOfDatasetTypes)
      
      ####################################
      ### Subsets seurat object per dataset type
      ####################################
      if (exists(x = "seurat.object.each_dataset_type") == T) {
        rm(seurat.object.each_dataset_type)
      }
      
      Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset_type
      
      seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
      print(seurat.object.each_dataset_type)
      
      ####################################
      ### Re-clusters each dataset type
      ####################################
      
      StopWatchStart$ReclusterEachDatasetTypeCellsFromInteg[[dataset_type]]  <- Sys.time()
      
      options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
      seurat.object.each_dataset_type <- FindNeighbors(object = seurat.object.each_dataset_type, dims = PcaDimsUse)
      seurat.object.each_dataset_type <- FindClusters(object = seurat.object.each_dataset_type, resolution = Resolution)
      ClustersThisDatasetType <- unlist(x = strsplit(x = paste(dataset_type, seurat.object.each_dataset_type@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
      
      StopWatchEnd$ReclusterEachDatasetTypeCellsFromInteg[[dataset_type]]  <- Sys.time()
      
      ####################################
      ### Write out each dataset type re-clustered cell clusters
      ####################################
      writeLines("\n*** Write out each dataset type re-clustered cell clusters ***\n")
      
      StopWatchStart$WriteReClustersEachDatasetTypeTables[[dataset_type]] <- Sys.time()
      
      write.table(paste(colnames(seurat.object.each_dataset_type), ClustersThisDatasetType, sep = "\t", collapse = "\n"),
                  file = OutfileEachDatasetTypeCellReClusters, row.names = F, col.names = F, quote = F, append = T)
      
      CellNames<-rownames(seurat.object.each_dataset_type@meta.data)
      ClusterIdent <-seurat.object.each_dataset_type@meta.data$seurat_clusters
      NumberOfClusters<-length(unique(ClusterIdent))
      OutfileNumbClusters<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", dataset_type, "_NumbCellClusters", ".tsv")
      write(x=NumberOfClusters, file = OutfileNumbClusters)
      
      StopWatchEnd$WriteReClustersEachDatasetTypeTables[[dataset_type]] <- Sys.time()
      
      ####################################
      ### Colour dimension reduction plots by each dataset type re-clustered cells
      ####################################
      
      writeLines("\n*** Colour dimension reduction plots by each dataset type re-clustered cells ***\n")
      
      for (dim_red_method in names(DimensionReductionMethods)) {
        
        StopWatchStart$DimRedPlotsReclusteredByDatasetTypeColuredByRecluster[[dim_red_method]][[dataset_type]]  <- Sys.time()
        
        plots <- DimPlot(seurat.object.each_dataset_type, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
        plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
        IntegratedUMAPPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_", "EachDatasetTypeReclustered_", dataset_type, ".pdf")
        pdf(file=IntegratedUMAPPlotPdf, width = 7, height = 8)
        print(CombinePlots(plots))
        dev.off()
        
        StopWatchEnd$DimRedPlotsReclusteredByDatasetTypeColuredByRecluster[[dim_red_method]][[dataset_type]] <- Sys.time()
      }
    }
    
    ####################################
    ### Load each dataset type re-cluster assignments
    ####################################
    writeLines("\n*** Load each dataset type re-cluster assignments ***\n")
    
    EachDatasetTypeCellReClusters <- data.frame(read.table(OutfileEachDatasetTypeCellReClusters, header = T, row.names = 1, sep = "\t"))
    seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetTypeCellReClusters, col.name = "EachDatasetTypeCellReClusters")
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetTypeCellReClusters
    
    ####################################
    ### Get average expression for each dataset type re-clustered clusters
    ####################################
    writeLines("\n*** Get average expression for each dataset type re-clustered clusters ***\n")
    
    StopWatchStart$AverageGeneExpressionEachDatasetTypeReCluster <- Sys.time()
    
    for (assay_expression in DefaultParameters$AssaysForAverageGETables) {
      cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F, assays = assay_expression)
      
      OutfileClusterAverages<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_AverageGeneExpression_", "EachDatasetTypeReclustered_", listAssaySuffixForOutfiles[[assay_expression]], ".tsv")
      Headers<-paste("AVERAGE_GENE_EXPRESSION", paste("r", Resolution, "_C",  names(cluster.averages[[assay_expression]]), sep="", collapse="\t"), sep="\t", collapse = "\t")
      write.table(Headers, file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
      write.table(data.frame(cluster.averages[[assay_expression]]), file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
    }
    
    StopWatchEnd$AverageGeneExpressionEachDatasetTypeReCluster <- Sys.time()
    
  }
}else{
  writeLines("\n*** No dataset type analyzes will be conducted because only '1' dataset type was found in -inputs_list ***\n")
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SAVE THE WHOLE RUN R_OBJECT AND LOG FILES
################################################################################################################################################
################################################################################################################################################

####################################
### Saving the full R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving the full R object ***\n")
  
  StopWatchStart$SaveRDSFull  <- Sys.time()
  
  OutfileRDS<-paste0(Tempdir, "/R_OBJECTS/", PrefixOutfiles, ".", ProgramOutdir, "_PcaClusteringDimReduction", ".rds")
  saveRDS(seurat.object.integrated, file = OutfileRDS)
  
  StopWatchEnd$SaveRDSFull  <- Sys.time()
  
}else{
  
  writeLines("\n*** Not saving the full R object ***\n")
  
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_LogFiles_", "CPUtimes", ".txt")
write(file = OutfileCPUtimes, x = paste("#Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
write(file = OutfileCPUtimes, x = paste("#MaxGlobalVariables", MaxGlobalVariables, sep = "\t", collapse = ""), append = T)

Headers<-paste("Step", "Time(minutes)", sep="\t")
write.table(Headers, file = OutfileCPUtimes, row.names = F, col.names = F, sep="\t", quote = F, append = T)

lapply(names(StopWatchStart), function(stepToClock) {
  if (regexpr("POSIXct", class(StopWatchStart[[stepToClock]]), ignore.case = T)[1] == 1) {
    TimeStart <- StopWatchStart[[stepToClock]]
    TimeEnd   <- StopWatchEnd[[stepToClock]]
    TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
    ReportTime<-c(paste(stepToClock, TimeDiff, sep = "\t", collapse = ""))
    write(file = OutfileCPUtimes, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
  }else{
    lapply(names(StopWatchStart[[stepToClock]]), function(substep) {
      if (regexpr("POSIXct", class(StopWatchStart[[stepToClock]][[substep]]), ignore.case = T)[1] == 1) {
        TimeStart <- StopWatchStart[[stepToClock]][[substep]]
        TimeEnd   <- StopWatchEnd[[stepToClock]][[substep]]
        TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
        ReportTime<-c(paste(paste0(stepToClock, "-", substep), TimeDiff, sep = "\t", collapse = ""))
        write(file = OutfileCPUtimes, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
      }
    })
  }
})

####################################
### Moving outfiles into outdir or keeping them at tempdir (if using CWL)
####################################

writeLines("\n*** Moving outfiles into outdir or keeping them at tempdir (if using CWL) ***\n")

### using two steps to copy files (`file.copy` and `file.remove`) instead of just `file.rename` to avoid issues with path to Tempdir in cluster systems
if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(Tempdir)
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    if (DirName == "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS") {
      sapply(list.dirs(TempdirWithData, full.names = F, recursive = F), FUN=function(SubDirName) {
        sapply(list.dirs(paste0(TempdirWithData, "/", SubDirName), full.names = F, recursive = F), FUN=function(SubSubDirName) {
          OutdirFinal <- gsub(pattern = Tempdir, replacement =  paste0(Outdir, "/", ProgramOutdir), x = paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName))
          dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
          sapply(list.files(paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName), pattern = ".pdf", full.names = F), FUN=function(EachFileName) {
            file.copy(from=paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
            file.remove(paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName, "/", EachFileName))
          })
        })
      })
    }else{
      OutdirFinal <- paste0(Outdir, "/", ProgramOutdir, "/", DirName)
      print(OutdirFinal)
      dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
      sapply(list.files(TempdirWithData, pattern = paste0("^", PrefixOutfiles, ".", ProgramOutdir), full.names = F), FUN=function(EachFileName) {
        file.copy(from=paste0(TempdirWithData, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
        file.remove(paste0(TempdirWithData, "/", EachFileName))
      })
    }
  })
}

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
OutfileCPUtimes <- gsub(x = OutfileCPUtimes, pattern = Tempdir, replacement = paste0(Outdir, "/", ProgramOutdir))
writeLines(paste0("END - All done!!! See:\n", OutfileCPUtimes, "\nfor computing times report"))

quit()
