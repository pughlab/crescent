####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/v3.2/sctransform_vignette.html (SCtransform normalization)
### https://satijalab.org/seurat/v3.2/integration.html (general integration)
### https://satijalab.org/seurat/v3.2/immune_alignment.html (control vs. treatment)
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads scRNA-seq data, e.g. from Cell Ranger (MTX files)
### 2) Generates QC plots for each dataset
### 3) Normalizes each dataset
### 4) Saves an R object for each normalized dataset
### 5) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_Seurat_v3_MultiDatasets_QC_Normalization.R -h'
### for help
####################################

####################################
### Tested with R v4.0.2
####################################

####################################
### Required libraries
####################################
writeLines("\n**** LOAD REQUIRED LIBRARIES ****\n")

suppressPackageStartupMessages(library(DropletUtils)) # (Bioconductor) to handle MTX/H5 format files. Note it has about the same speed than library(earlycross) which can't handle H5
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

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SETUP RUN
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** SETUP RUN ****\n")

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_Seurat_v3_MultiDatasets_QC_Normalization.R"
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
  make_option(c("-j", "--inputs_remove_barcodes"), default="NA",
              help="Path/name to a <tab> delimited list of barcodes to be removed from analysis, like:
                d1  AAACCTGAGCTCCCAG
                d2  AAACCTGTCACCATAG
                d3  AAACCTGTCAGCTTAG
                
                Note: the first column must be dataset ID and the second column the barcode
                Or type 'NA' to include all barcodes
                
                Default = 'NA'"),
  #
  make_option(c("-k", "--save_filtered_data"), default="N",
              help="Indicates if filtered raw and normalized data should be saved as MTX files. Type 'y/Y' or 'n/N'.
              
                Default = 'N'"),
  #
  make_option(c("-l", "--save_unfiltered_data"), default="N",
              help="Indicates if unfiltered raw data (exactly as inputted by --inputs_list) should be saved as MTX files. Type 'y/Y' or 'n/N'.
                This may be useful to share with collaborators along results in a single zipped file.
                
                Default = 'N'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
              
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
              
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system
              
                Default = 'MAX'"),
  #
  make_option(c("-s", "--save_r_object"), default="Y",
              help="Indicates if R objects with the data and analyzes from each dataset should be saved
                Note that this may be time consuming. Type 'y/Y' or 'n/N'
                
                Default = 'Y'"),
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
InfileRemoveBarcodes    <- opt$inputs_remove_barcodes
SaveFilteredData        <- opt$save_filtered_data
SaveUnFilteredData      <- opt$save_unfiltered_data
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
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
  dir.create(file.path("R_OBJECTS_CWL"), showWarnings = F) 
  
  FILE_TYPE_OUT_DIRECTORIES = c(
    "CRESCENT_CLOUD",
    "CRESCENT_CLOUD/frontend_qc",
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS",
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
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS",
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

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_QC_Normalization_UsedOptions", ".txt")

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

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_QC_Normalization_RSessionInfo", ".txt")
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

DefaultParameters <- list(
  ### Parameters for QC plots
  CellPropertiesToQC = c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction"),
  
  ### Parameters for Seurat filters
  MinCells = 3

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

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("infiles_list", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO LOAD AND QC DATASETS
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** LOAD AND QC DATASETS ****\n")

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
    MergeFilter <- c("name","dataset_ID","dataset_path", "dataset_type", "dataset_format", "mito_min", "mito_max", "ribo_min", "ribo_max", "ngenes_min", "ngenes_max", "nreads_min", "nreads_max")
    MergedInputsTableFiltered <- MergedInputsTable[MergeFilter]
    MergedInputsTableFilteredFinal <- MergedInputsTableFiltered[,-1]
    rownames(MergedInputsTableFilteredFinal) <- MergedInputsTableFiltered[,1]
    colnames(MergedInputsTableFilteredFinal) <-c("DatasetMinioID","PathToDataset","DatasetType","DatasetFormat","MinMitoFrac","MaxMitoFrac","MinRiboFrac","MaxRiboFrac","MinNGenes","MaxNGenes","MinNReads","MaxNReads")
    
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

####################################
### Load scRNA-seq data and generate QC plots and tables
####################################
writeLines("\n*** Load scRNA-seq data and generate QC plots and tables ***\n")

SeuratObjectsFiltered   <-list()
SeuratObjectsUnfiltered <-list()
DatasetIds              <-list()
list_DatasetToType      <-list()
list_TypeToDatasets     <-list()
list_DatasetToFormat    <-list()
list_MinMitoFrac        <-list()
list_MaxMitoFrac        <-list()
list_MinRiboFrac        <-list()
list_MaxRiboFrac        <-list()
list_MinNGenes          <-list()
list_MaxNGenes          <-list()
list_MinNReads          <-list()
list_MaxNReads          <-list()

if ((regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) & (!(regexpr("^NA$", MinioPath , ignore.case = T)[1] == 1))) {
  list_DatasetMinioIDs         <-list()
}

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  Dataset.SO <-paste0(dataset, ".so")
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    PathToDataset <- InputsTable[dataset,"PathToDataset"]
  } else {
    PathToDataset <- gsub("^~/",paste0(UserHomeDirectory,"/"), InputsTable[dataset,"PathToDataset"])
  }
  
  DatasetType   <- InputsTable[dataset,"DatasetType"]
  
  list_DatasetToType[[dataset]]      <- DatasetType
  list_DatasetToFormat[[dataset]]    <- InputsTable[dataset,"DatasetFormat"]
  list_TypeToDatasets[[DatasetType]] <- append(list_TypeToDatasets[[DatasetType]], dataset)
  list_MinMitoFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MinMitoFrac"])
  list_MaxMitoFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MaxMitoFrac"])
  list_MinRiboFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MinRiboFrac"])
  list_MaxRiboFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MaxRiboFrac"])
  list_MinNGenes[[dataset]]          <- as.numeric(InputsTable[dataset,"MinNGenes"])
  list_MaxNGenes[[dataset]]          <- as.numeric(InputsTable[dataset,"MaxNGenes"])
  list_MinNReads[[dataset]]          <- as.numeric(InputsTable[dataset,"MinNReads"])
  list_MaxNReads[[dataset]]          <- as.numeric(InputsTable[dataset,"MaxNReads"])
  
  if ((regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) & (!(regexpr("^NA$", MinioPath , ignore.case = T)[1] == 1))) {
    DatasetMinioID   <- InputsTable[dataset,"DatasetMinioID"]
    list_DatasetMinioIDs[[dataset]]      <- DatasetMinioID
  }
  
  if (regexpr("^MTX$|^TSV$|^HDF5$", list_DatasetToFormat[[dataset]], ignore.case = T, perl = T)[1] == 1) {
    
    ####################################
    ### Load scRNA-seq data
    ####################################
    
    StopWatchStart$LoadScRNAseqData[[dataset]] <- Sys.time()
    
    if (regexpr("^MTX$", list_DatasetToFormat[[dataset]], ignore.case = T)[1] == 1) {
      writeLines(paste0("\n*** Loading MTX infiles for ", dataset, " from: ", PathToDataset, " ***\n"))
      expression_matrix <- Read10X(data.dir = PathToDataset, strip.suffix = T) ### Note `strip.suffix = T` applies to Seurat v3.2 or higher
    }else if (regexpr("^TSV$", list_DatasetToFormat[[dataset]], ignore.case = T)[1] == 1) {
      writeLines(paste0("\n*** Loading matrix of genes (rows) vs. barcodes (columns) for dataset: ", dataset, " from: ", PathToDataset, " ***\n"))
      ## Note `check.names = F` is needed for both `fread` and `data.frame`
      expression_matrix <- as.matrix(data.frame(fread(PathToDataset, check.names = F), row.names=1, check.names = F))
    }else if (regexpr("^HDF5$", list_DatasetToFormat[[dataset]], ignore.case = T)[1] == 1) {
      writeLines(paste0("\n*** Loading HDF5 infile for ", dataset, " from: ", PathToDataset, " ***\n"))
      expression_matrix <- Read10X_h5(filename = PathToDataset, use.names = T, unique.features = T)
    }else{
      stop(paste0("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript ", ThisScriptName, " -h\n\n"))
    }
    dim(expression_matrix)
    
    StopWatchEnd$LoadScRNAseqData[[dataset]]  <- Sys.time()
    
    ####################################
    ### Save unfiltered data
    ####################################
    
    if (regexpr("^Y$", SaveUnFilteredData, ignore.case = T)[1] == 1) {
      writeLines("\n*** Save unfiltered data ***\n")
      
      StopWatchStart$SaveUnFilteredData[[dataset]]  <- Sys.time()
      
      ### This is to write out the input data with no filters at all
      seurat.object.tmp  <- CreateSeuratObject(counts = expression_matrix, min.cells = 0, min.features = 0, project = paste0(PrefixOutfiles, "_", dataset))
      
      OutDirUnfilteredRaw <-paste0(Tempdir, "/UNFILTERED_DATA_MATRICES/", dataset, "/RAW")
      dir.create(file.path(OutDirUnfilteredRaw), showWarnings = F, recursive = T)
      write10xCounts(path = OutDirUnfilteredRaw, x = seurat.object.tmp@assays[["RNA"]]@counts, gene.type="Raw_Gene_Expression", overwrite=T, type="sparse", version="3")
      
      rm(seurat.object.tmp)
      
      StopWatchEnd$SaveUnFilteredData[[dataset]]  <- Sys.time()
      
    }
    
    ####################################
    ### Create seurat object
    ####################################
    
    StopWatchStart$CreateSeuratObject[[dataset]]  <- Sys.time()
    
    writeLines(paste0("\n*** Create seurat object for ", dataset, " ***\n"))
    seurat.object.u <- CreateSeuratObject(counts = expression_matrix, min.cells = DefaultParameters$MinCells, min.features = list_MinNGenes[[dataset]], project = paste0(PrefixOutfiles, "_", dataset))
    seurat.object.u <- AddMetaData(object = seurat.object.u, metadata = colnames(seurat.object.u), col.name = 'barcode_sequence')
    
    StopWatchEnd$CreateSeuratObject[[dataset]]  <- Sys.time()
    
    ####################################
    ### Add dataset label
    ####################################
    writeLines(paste0("\n*** Add dataset label for ", dataset, " ***\n"))
    
    StopWatchStart$AddDatasetLabel[[dataset]]  <- Sys.time()
    
    seurat.object.u[['dataset.label']] <- dataset
    
    StopWatchEnd$AddDatasetLabel[[dataset]]  <- Sys.time()
    
    ####################################
    ### Add dataset type label
    ####################################
    writeLines(paste0("\n*** Add dataset type label for ", dataset, " ***\n"))
    
    StopWatchStart$AddDatasetTypeLabel[[dataset]]  <- Sys.time()
    
    seurat.object.u[['dataset_type.label']] <- list_DatasetToType[[dataset]]
    
    StopWatchEnd$AddDatasetTypeLabel[[dataset]]  <- Sys.time()
    
    ####################################
    ### Get mitochondrial genes
    ####################################
    writeLines(paste0("\n*** Get  mitochondrial genes for ", dataset, " ***\n"))
    
    StopWatchStart$GetMitoGenes[[dataset]]  <- Sys.time()
    
    mitoRegExpressions<- paste(c("^MT-"), collapse = "|")
    mito.features <- grep(pattern = mitoRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)
    
    if (length(mito.features)[[1]] > 0) {
      mito.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['mito.fraction']] <- mito.fraction
    }else{
      mito.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['mito.fraction']] <- mito.fraction
    }
    
    StopWatchEnd$GetMitoGenes[[dataset]]  <- Sys.time()
    
    ####################################
    ### Get ribosomal protein genes
    ####################################
    writeLines(paste0("\n*** Get ribosomal protein genes for ", dataset, " ***\n"))
    
    StopWatchStart$GetRiboGenes[[dataset]]  <- Sys.time()
    
    riboRegExpressions<- paste(c("^MRPL", "^MRPS", "^RPL", "^RPS"), collapse = "|")
    ribo.features <- grep(pattern = riboRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)
    
    if (length(ribo.features)[[1]] > 0) {
      ribo.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['ribo.fraction']] <- ribo.fraction
    }else{
      ribo.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['ribo.fraction']] <- ribo.fraction
    }
    
    StopWatchEnd$GetRiboGenes[[dataset]]  <- Sys.time()
    
    ####################################
    ### Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation
    ####################################
    writeLines(paste0("\n*** Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation for ", dataset, " ***\n"))
    
    StopWatchStart$FilterCells[[dataset]]  <- Sys.time()
    
    if (length(mito.features)[[1]] > 0) {
      seurat.object.f<-subset(x = seurat.object.u, subset = 
                                nFeature_RNA >= list_MinNGenes[[dataset]]
                              & nFeature_RNA <= list_MaxNGenes[[dataset]] 
                              & nCount_RNA   >= list_MinNReads[[dataset]]
                              & nCount_RNA   <= list_MaxNReads[[dataset]]
                              & mito.fraction >= list_MinMitoFrac[[dataset]]
                              & mito.fraction <= list_MaxMitoFrac[[dataset]]
                              & ribo.fraction >= list_MinRiboFrac[[dataset]]
                              & ribo.fraction <= list_MaxRiboFrac[[dataset]])
      
    }else{
      seurat.object.f<-subset(x = seurat.object.u, subset = 
                                nFeature_RNA >= list_MinNGenes[[dataset]]
                              & nFeature_RNA <= list_MaxNGenes[[dataset]] 
                              & nCount_RNA   >= list_MinNReads[[dataset]]
                              & nCount_RNA   <= list_MaxNReads[[dataset]]
                              & ribo.fraction >= list_MinRiboFrac[[dataset]]
                              & ribo.fraction <= list_MaxRiboFrac[[dataset]])
    }
    ### Get list of barcodes excluded by nFeature_RNA, nCount_RNA, mito.fraction or ribo.fraction
    NumberOfBarcodesExcludedByNFeature        <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA  >= list_MinNGenes[[dataset]] & nFeature_RNA <= list_MaxNGenes[[dataset]]))))
    NumberOfBarcodesExcludedByNReads          <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nCount_RNA    >= list_MinNReads[[dataset]] & nCount_RNA   <= list_MaxNReads[[dataset]]))))
    NumberOfBarcodesExcludedByMito            <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = mito.fraction >= list_MinMitoFrac[[dataset]] & mito.fraction <= list_MaxMitoFrac[[dataset]]))))
    NumberOfBarcodesExcludedByRibo            <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = ribo.fraction >= list_MinRiboFrac[[dataset]] & ribo.fraction <= list_MaxRiboFrac[[dataset]]))))
    
    print(paste0("NumberOfBarcodesExcludedByNFeature=", NumberOfBarcodesExcludedByNFeature))
    print(paste0("NumberOfBarcodesExcludedByNReads=", NumberOfBarcodesExcludedByNReads))
    print(paste0("NumberOfBarcodesExcludedByMito=", NumberOfBarcodesExcludedByMito))
    print(paste0("NumberOfBarcodesExcludedByRibo=", NumberOfBarcodesExcludedByRibo))
    
    StopWatchEnd$FilterCells[[dataset]]  <- Sys.time()
    
    ### Just reporting the summary of the UNfiltered and filtered objects
    print(paste0(dataset, "  unfiltered"))
    print(seurat.object.u)
    print(paste0(dataset, "  filtered"))
    print(seurat.object.f)
    
    ####################################
    ### QC EDA violin plots
    ####################################
    writeLines(paste0("\n*** QC EDA violin plots for ", dataset, " ***\n"))
    
    StopWatchStart$QCviolinplots[[dataset]]  <- Sys.time()
    
    ### Get unfiltered data QC statistics
    nFeature_RNA.u.df   <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
    nCount_RNA.u.df     <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
    mito.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$mito.fraction, mito.fraction = 1)
    ribo.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$ribo.fraction, ribo.fraction = 1)
    #
    mean_nFeature_RNAStats.u   <- round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
    median_nFeature_RNAStats.u <- round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
    mean_nCount_RNAStats.u     <- round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0)
    median_nCount_RNAStats.u   <- round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)
    mean_mito.fraction.u       <- round(mean(seurat.object.u@meta.data[,"mito.fraction"]),3)
    median_mito.fraction.u     <- round(median(seurat.object.u@meta.data[,"mito.fraction"]),3)
    mean_ribo.fraction.u       <- round(mean(seurat.object.u@meta.data[,"ribo.fraction"]),3)
    median_ribo.fraction.u     <- round(median(seurat.object.u@meta.data[,"ribo.fraction"]),3)
    #
    nFeature_RNAStats.u <-paste0(c("mean = ", mean_nFeature_RNAStats.u, "\n", "median = ", median_nFeature_RNAStats.u))
    nCount_RNAStats.u   <-paste0(c("mean = ", mean_nCount_RNAStats.u,   "\n", "median = ", median_nCount_RNAStats.u))
    mito.fraction.u     <-paste0(c("mean = ", mean_mito.fraction.u,     "\n", "median = ", median_mito.fraction.u))
    ribo.fraction.u     <-paste0(c("mean = ", mean_ribo.fraction.u,     "\n", "median = ", median_ribo.fraction.u))
    
    ### Get filtered data QC statistics
    nFeature_RNA.f.df   <-data.frame(Expression_level = seurat.object.f@meta.data$nFeature_RNA, nGenes = 2)
    nCount_RNA.f.df     <-data.frame(Expression_level = seurat.object.f@meta.data$nCount_RNA,   nCount_RNA = 2)
    mito.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$mito.fraction, mito.fraction = 2)
    ribo.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$ribo.fraction, ribo.fraction = 2)
    #
    mean_nFeature_RNAStats.f   <- round(mean(seurat.object.f@meta.data[,"nFeature_RNA"]),0)
    median_nFeature_RNAStats.f <- round(median(seurat.object.f@meta.data[,"nFeature_RNA"]),0)
    mean_nCount_RNAStats.f     <- round(mean(seurat.object.f@meta.data[,"nCount_RNA"]),0)
    median_nCount_RNAStats.f   <- round(median(seurat.object.f@meta.data[,"nCount_RNA"]),0)
    mean_mito.fraction.f       <- round(mean(seurat.object.f@meta.data[,"mito.fraction"]),3)
    median_mito.fraction.f     <- round(median(seurat.object.f@meta.data[,"mito.fraction"]),3)
    mean_ribo.fraction.f       <- round(mean(seurat.object.f@meta.data[,"ribo.fraction"]),3)
    median_ribo.fraction.f     <- round(median(seurat.object.f@meta.data[,"ribo.fraction"]),3)
    #
    nFeature_RNAStats.f <-paste0(c("mean = ", mean_nFeature_RNAStats.f, "\n", "median = ", median_nFeature_RNAStats.f))
    nCount_RNAStats.f   <-paste0(c("mean = ", mean_nCount_RNAStats.f,   "\n", "median = ", median_nCount_RNAStats.f))
    mito.fraction.f     <-paste0(c("mean = ", mean_mito.fraction.f,     "\n", "median = ", median_mito.fraction.f))
    ribo.fraction.f     <-paste0(c("mean = ", mean_ribo.fraction.f,     "\n", "median = ", median_ribo.fraction.f))
    
    ### Put QC statistics together
    nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
    nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
    mito.fraction.m.df  <-data.frame(rbind(mito.fraction.u.df,mito.fraction.f.df))
    ribo.fraction.m.df  <-data.frame(rbind(ribo.fraction.u.df,ribo.fraction.f.df))
    NumberOfCells<- list()
    NumberOfCells[["unfiltered"]] <- nrow(seurat.object.u@meta.data)
    NumberOfCells[["filtered"]]   <- nrow(seurat.object.f@meta.data)
    LabelUnfiltered    <-paste0("Before filters: No. of cells = ", NumberOfCells[["unfiltered"]])
    LabelFiltered      <-paste0("After filters:  No. of cells = ", NumberOfCells[["filtered"]])
    
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
      labs(x=paste0("No. of genes", "\n", "Filter: ", "min=", list_MinNGenes[[dataset]], " max=", list_MaxNGenes[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNFeature)) +
      annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])
    
    nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
      geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x=paste0("No. of reads", "\n", "Filter: ", "min=", list_MinNReads[[dataset]], " max=", list_MaxNReads[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNReads)) +
      annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])
    
    mito.fraction.plot<-ggplot(data=mito.fraction.m.df, aes(x = factor(mito.fraction), y = Expression_level)) +
      geom_violin(aes(fill = factor(mito.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x="Mitochondrial genes (fraction)") +
      labs(x=paste0("Mitochondrial genes", "\n", "Filter: ", "min=", list_MinMitoFrac[[dataset]], " max=", list_MaxMitoFrac[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByMito)) +
      annotate("text", x = 1 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.f, col = ColoursQCViolinPlots[[2]])
    
    ribo.fraction.plot<-ggplot(data=ribo.fraction.m.df, aes(x = factor(ribo.fraction), y = Expression_level)) +
      geom_violin(aes(fill = factor(ribo.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x=paste0("Ribosomal protein genes (fraction)", "\n", "Filter: ", "min=", list_MinRiboFrac[[dataset]], " max=", list_MaxRiboFrac[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByRibo)) +
      annotate("text", x = 1 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.f, col = ColoursQCViolinPlots[[2]])
    
    bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, mito.fraction.plot, ribo.fraction.plot, ncol = 4)
    
    ### Create a *pdf file with the violin ggplot's
    
    VlnPlotPdf<-paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", "ViolinPlots", "_", dataset, ".pdf")
    pdf(file=VlnPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 1.7, height = DefaultParameters$BaseSizeSinglePlotPdf)
    print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
    dev.off()
    
    StopWatchEnd$QCviolinplots[[dataset]]  <- Sys.time()
    
    if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
      ####################################
      ### Outfiles for web app: interactive qc plots
      ####################################
      writeLines("\n*** Outfiles for web app: interactive qc plots ***\n")
      
      StopWatchStart$QCviolinplotsFrontEnd[[dataset]]  <- Sys.time()
      
      # unfiltered
      interactive_qc_plot_u  <-data.frame(Barcodes = row.names(seurat.object.u@meta.data), Number_of_Genes = seurat.object.u@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.u@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.u@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.u@meta.data$ribo.fraction)
      interactive_qc_plot_u$Mitochondrial_Genes_Percentage <- interactive_qc_plot_u$Mitochondrial_Genes_Percentage * 100
      interactive_qc_plot_u$Ribosomal_Protein_Genes_Percentage <- interactive_qc_plot_u$Ribosomal_Protein_Genes_Percentage * 100
      colnames(interactive_qc_plot_u) <- c("Barcodes","Number of Genes","Number of Reads","Mitochondrial Genes Percentage","Ribosomal Protein Genes Percentage")
      write.table(interactive_qc_plot_u, paste(Tempdir, "/CRESCENT_CLOUD/frontend_qc/",dataset,"_BeforeFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
      
      # filtered
      interactive_qc_plot_f  <-data.frame(Barcodes = row.names(seurat.object.f@meta.data), Number_of_Genes = seurat.object.f@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.f@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.f@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.f@meta.data$ribo.fraction)
      interactive_qc_plot_f$Mitochondrial_Genes_Percentage <- interactive_qc_plot_f$Mitochondrial_Genes_Percentage * 100
      interactive_qc_plot_f$Ribosomal_Protein_Genes_Percentage <- interactive_qc_plot_f$Ribosomal_Protein_Genes_Percentage * 100
      colnames(interactive_qc_plot_f) <- c("Barcodes","Number of Genes","Number of Reads","Mitochondrial Genes Percentage","Ribosomal Protein Genes Percentage")
      write.table(interactive_qc_plot_f, paste(Tempdir, "/CRESCENT_CLOUD/frontend_qc/",dataset,"_AfterFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
      
      qc_tsv <- data.frame(NAME = row.names(seurat.object.f@meta.data), Number_of_Genes = seurat.object.f@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.f@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.f@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.f@meta.data$ribo.fraction)
      qc_tsv$NAME <- paste0(dataset, "_", qc_tsv$NAME)
      qc_tsv$Mitochondrial_Genes_Percentage <- qc_tsv$Mitochondrial_Genes_Percentage * 100
      qc_tsv$Ribosomal_Protein_Genes_Percentage <- qc_tsv$Ribosomal_Protein_Genes_Percentage * 100
      qc_tsv_string <- sapply(qc_tsv, as.character)
      qc_tsv_string_TYPE <- rbind(data.frame(NAME = "TYPE", Number_of_Genes = "numeric", Number_of_Reads = "numeric", Mitochondrial_Genes_Percentage = "numeric", Ribosomal_Protein_Genes_Percentage = "numeric"), qc_tsv_string)
      
      qc_outfile <-paste0(Tempdir, "/CRESCENT_CLOUD/frontend_qc/",dataset,"_qc_data.tsv")
      write.table(data.frame(qc_tsv_string_TYPE), file = qc_outfile, row.names = F, col.names = T, sep="\t", quote = F, append = T)
      
      StopWatchEnd$QCviolinplotsFrontEnd[[dataset]]  <- Sys.time()
      
    } 
    
    ####################################
    ### Feature-vs-feature scatter plot
    ####################################
    writeLines(paste0("\n*** Feature-vs-feature scatter plot for ", dataset, " ***\n"))
    
    StopWatchStart$FeatureVsFeatureplot[[dataset]]  <- Sys.time()
    
    UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                                  nGene = seurat.object.u@meta.data$nFeature_RNA,
                                  mito.fraction = seurat.object.u@meta.data$mito.fraction,
                                  filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.f))
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = ColoursQCViolinPlots[[1]])
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = F, replacement = ColoursQCViolinPlots[[2]])
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = 4)
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = F, replacement = 16)
    
    FeatureVsFeaturePlotPdf<-paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", "NumbReadsVsNumbGenesAndMito_ScatterPlot", "_", dataset, ".pdf")
    pdf(file=FeatureVsFeaturePlotPdf, width = 10, height = 5)
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
    
    StopWatchEnd$FeatureVsFeatureplot[[dataset]]  <- Sys.time()
    
    ####################################
    ### Write out filter details and number of filtered cells
    ####################################
    writeLines(paste0("\n*** Write out filter details and number of filtered cells for ", dataset, " ***\n"))
    
    StopWatchStart$OutTablesFilterDetailsAndFilteredCells[[dataset]]  <- Sys.time()
    
    Outfile.con <- bzfile(paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", "FilterDetails", "_", dataset, ".tsv.bz2"), "w")
    Headers<-paste("Step", "Filter_min", "Filter_max", "Mean_before_filter", "Median_before_filter", "Mean_after_filter", "Median_after_filter", "Excluded_cells", sep = "\t", collapse = "")
    write.table(Headers, file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F)
    
    FilterDetails.nFeature_RNA  <- paste("nFeature_RNA", list_MinNGenes[[dataset]], list_MaxNGenes[[dataset]],
                                         mean_nFeature_RNAStats.u, median_nFeature_RNAStats.u, mean_nFeature_RNAStats.f, median_nFeature_RNAStats.f, 
                                         NumberOfBarcodesExcludedByNFeature, sep = "\t", collapse = "")
    FilterDetails.nCount_RNA    <- paste("nCount_RNA", list_MinNReads[[dataset]], list_MaxNReads[[dataset]], 
                                         mean_nCount_RNAStats.u, median_nCount_RNAStats.u, mean_nCount_RNAStats.f, median_nCount_RNAStats.f, 
                                         NumberOfBarcodesExcludedByNReads, sep = "\t", collapse = "")
    FilterDetails.mito.fraction <- paste("mito.fraction", list_MinMitoFrac[[dataset]], list_MaxMitoFrac[[dataset]],
                                         mean_mito.fraction.u, median_mito.fraction.u, mean_mito.fraction.f, median_mito.fraction.f, 
                                         NumberOfBarcodesExcludedByMito, sep = "\t", collapse = "")
    FilterDetails.ribo.fraction <- paste("ribo.fraction", list_MinRiboFrac[[dataset]], list_MaxRiboFrac[[dataset]],
                                         mean_ribo.fraction.u, median_ribo.fraction.u, mean_ribo.fraction.f, median_ribo.fraction.f, 
                                         NumberOfBarcodesExcludedByRibo, sep = "\t", collapse = "")
    write.table(FilterDetails.nFeature_RNA,  file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    write.table(FilterDetails.nCount_RNA,    file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    write.table(FilterDetails.mito.fraction, file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    write.table(FilterDetails.ribo.fraction, file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    close(Outfile.con)
    
    Outfile.con <- bzfile(paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", "NumberOfFilteredCells", "_", dataset, ".tsv.bz2"), "w")
    write.table(paste("Number_of_cells_before_filters", NumberOfCells[["unfiltered"]], sep = "\t", collapse = "\n"), file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(paste("Number_of_cells_after_filters", NumberOfCells[["filtered"]],    sep = "\t", collapse = "\n"), file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    close(Outfile.con)
    
    StopWatchEnd$OutTablesFilterDetailsAndFilteredCells[[dataset]]  <- Sys.time()
    
    if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
      ####################################
      ### Outfiles for web app: write out QC metrics for each dataset
      ####################################
      writeLines("\n*** Outfiles for web app: write out QC metrics for each dataset ***\n")
      
      StopWatchStart$OutTablesFilterDetailsAndFilteredCellsFrontEnd[[dataset]]  <- Sys.time()
      
      qc_metrics_genes <- data.frame(Step = "Number of Genes", Filter_min = list_MinNGenes[[dataset]], Filter_max = list_MaxNGenes[[dataset]], Excluded_cells = NumberOfBarcodesExcludedByNFeature)
      qc_metrics_reads <- data.frame(Step = "Number of Reads", Filter_min = list_MinNReads[[dataset]], Filter_max = list_MaxNReads[[dataset]], Excluded_cells = NumberOfBarcodesExcludedByNReads)
      qc_metrics_mito <- data.frame(Step = "Percentage of Mitochondrial Genes", Filter_min = as.numeric(list_MinMitoFrac[[dataset]])*100, Filter_max = as.numeric(list_MaxMitoFrac[[dataset]])*100, Excluded_cells = NumberOfBarcodesExcludedByMito)
      qc_metrics_ribo <- data.frame(Step = "Percentage of Ribsomal Protein Genes", Filter_min = as.numeric(list_MinRiboFrac[[dataset]])*100, Filter_max = as.numeric(list_MaxRiboFrac[[dataset]])*100, Excluded_cells = NumberOfBarcodesExcludedByRibo)
      
      qc_metrics <-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",dataset,"_qc_metrics.tsv")
      write.table(qc_metrics_genes, file = qc_metrics, row.names = F, col.names = T, sep="\t", quote = F, append = T)
      write.table(qc_metrics_reads, file = qc_metrics, row.names = F, col.names = F, sep="\t", quote = F, append = T)
      write.table(qc_metrics_mito, file = qc_metrics, row.names = F, col.names = F, sep="\t", quote = F, append = T)
      write.table(qc_metrics_ribo, file = qc_metrics, row.names = F, col.names = F, sep="\t", quote = F, append = T)
      
      StopWatchEnd$OutTablesFilterDetailsAndFilteredCellsFrontEnd[[dataset]]  <- Sys.time()
      
    } 
    
    ####################################
    ### Assign data to Datasets lists
    ####################################
    writeLines(paste0("\n*** Assign data to Datasets lists: ", dataset, " ***\n"))
    
    StopWatchStart$AssignDataToDatasets  <- Sys.time()
    
    SeuratObjectsUnfiltered[[as.character(NumberOfDatasets)]]  <- seurat.object.u
    SeuratObjectsFiltered[[as.character(NumberOfDatasets)]]    <- seurat.object.f
    DatasetIds[[as.character(NumberOfDatasets)]]               <- dataset
    
    StopWatchEnd$AssignDataToDatasets  <- Sys.time()
    
    ####################################
    ### Write out QC data for each dataset
    ####################################
    writeLines("\n*** Write out QC data for each dataset ***\n")
    
    StopWatchStart$WriteOutQCData  <- Sys.time()
    
    Headers<-paste("Cell_barcode", paste(DefaultParameters$CellPropertiesToQC, sep = "", collapse = "\t") ,sep="\t")
    
    BarcodeIdsWithDatasetBeforeFilters <- colnames(SeuratObjectsUnfiltered[[NumberOfDatasets]])
    Outfile.con <- bzfile(paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", "Before_filters_QC_metadata", "_", dataset, ".tsv.bz2"), "w")
    write.table(Headers, file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(BarcodeIdsWithDatasetBeforeFilters, SeuratObjectsUnfiltered[[NumberOfDatasets]]@meta.data[,DefaultParameters$CellPropertiesToQC]),
                file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    close(Outfile.con)
    
    BarcodeIdsWithDatasetAfterFilters <- colnames(SeuratObjectsFiltered[[NumberOfDatasets]])
    Outfile.con <- bzfile(paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_QC_", "After_filters_QC_metadata", "_", dataset, ".tsv.bz2"), "w")
    write.table(Headers, file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(BarcodeIdsWithDatasetAfterFilters, SeuratObjectsFiltered[[NumberOfDatasets]]@meta.data[,DefaultParameters$CellPropertiesToQC]),
                file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    close(Outfile.con)
    
    StopWatchEnd$WriteOutQCData  <- Sys.time()
    
    ####################################
    ### Remove the Unfiltered seurat object
    ####################################
    writeLines(paste0("\n*** Remove the Unfiltered seurat object for ", dataset, " ***\n"))
    
    rm(seurat.object.u)
    rm(UnfilteredData.df)
    
  }else{
    stop(paste0("Unexpected type of input: ", DatasetType, "\n\nFor help type:\n\nRscript integrates_datasets_with_seurat.R -h\n\n"))
  }
}

####################################
### Remove barcodes by parameter -j (if applicable)
####################################

if (regexpr("^NA$", InfileRemoveBarcodes , ignore.case = T)[1] == 1) {
  
  writeLines("\n*** No barcodes are removed by option -j ***\n")
  
}else{
  
  StopWatchStart$RemoveBarcodes <- Sys.time()
  
  writeLines("\n*** Removing barcodes by parameter -j ***\n")
  
  InfileAllBarcodesToRemove<-gsub("^~/",paste0(UserHomeDirectory,"/"), InfileRemoveBarcodes)
  AllBarcodesToRemove.tab<-read.table(InfileAllBarcodesToRemove, header = F, row.names = NULL, stringsAsFactors = FALSE)
  if (ncol(AllBarcodesToRemove.tab) == 2) {
    colnames(AllBarcodesToRemove.tab) <- c("Dataset","Barcode")
  }else{
    stop(paste0("Unexpected format in file:\n", InfileAllBarcodesToRemove, "\nfor parameter -j, 2 columns were expected but found ", ncol(AllBarcodesToRemove.tab)))
  }
  
  for (SeuratObjNumb in c(1:NumberOfDatasets)) {
    print(paste0("Removing cells from:", DatasetIds[[SeuratObjNumb]]))
    seurat.object.full <- SeuratObjectsFiltered[[SeuratObjNumb]]
    seurat.object.full
    DatasetId                      <- DatasetIds[[SeuratObjNumb]]
    ThisDatasetBarcodesToRemove    <- subset(x=AllBarcodesToRemove.tab, subset = Dataset == DatasetId)[,"Barcode"]
    ThisDatasetBarcodesToKeep.log  <- !colnames(seurat.object.full) %in% ThisDatasetBarcodesToRemove
    seurat.object.subset           <- subset(seurat.object.full, cells = colnames(seurat.object.full[,ThisDatasetBarcodesToKeep.log]))
    seurat.object.subset
    SeuratObjectsFiltered[[SeuratObjNumb]] <- seurat.object.subset
    print(paste(DatasetId, paste0("Before:", ncol(seurat.object.full)), paste0("After:", ncol(seurat.object.subset)), sep = "  ", collapse = "\n"))
  }
  
  StopWatchEnd$RemoveBarcodes <- Sys.time()
  
}

####################################
### Merge Seurat objects RNA assay
####################################
writeLines("\n*** Merge Seurat objects RNA assay ***\n")

StopWatchStart$MergeSeuratObjectsFilteredRNA  <- Sys.time()
FirstSeuratObject   <- SeuratObjectsFiltered[[1]]
RestOfSeuratObjectsFiltered <- SeuratObjectsFiltered[c(2:NumberOfDatasets)]
RestOfDatasetsIds    <- unlist(DatasetIds[c(2:NumberOfDatasets)])

seurat.object.merged <- merge(FirstSeuratObject, y = RestOfSeuratObjectsFiltered, 
                              add.cell.ids = DatasetIds,
                              project = PrefixOutfiles)

seurat.object.merged <- AddMetaData(object = seurat.object.merged, metadata = seurat.object.merged@meta.data$dataset.label, col.name = "dataset")
seurat.object.list <- SplitObject(seurat.object.merged, split.by = "dataset")

StopWatchEnd$MergeSeuratObjectsFilteredRNA  <- Sys.time()

####################################
### Running SCTransform
####################################
writeLines("\n**** NORMALIZE DATASETS ****\n")

writeLines("\n*** Running SCTransform ***\n")

StopWatchStart$SCTransform  <- Sys.time()

for (i in 1:length(seurat.object.list)) {
  
  dataset <- rownames(InputsTable)[[i]]
  print(dataset)
  seurat.object.list[[i]] <- SCTransform(seurat.object.list[[i]], verbose = T)
  
  ####################################
  ### Save filtered data as MTX files
  ####################################
  writeLines("\n*** Save filtered data as MTX files ***\n")
  
  StopWatchStart$SaveFilteredData[[dataset]]  <- Sys.time()
  
  if (regexpr("^Y$", SaveFilteredData, ignore.case = T)[1] == 1) {
    
    OutDirFilteredRaw <-paste0(Tempdir, "/FILTERED_DATA_MATRICES/", dataset, "/RAW")
    dir.create(file.path(OutDirFilteredRaw), showWarnings = F, recursive = T)
    write10xCounts(path = OutDirFilteredRaw, x = seurat.object.list[[i]]@assays[["RNA"]]@counts, gene.type="Raw_Gene_Expression", overwrite=T, type="sparse", version="3")
    
    OutDirFilteredNorm<-paste0(Tempdir, "/FILTERED_DATA_MATRICES/", dataset, "/NORMALIZED_SCT")
    dir.create(file.path(OutDirFilteredNorm), showWarnings = F, recursive = T)
    write10xCounts(path = OutDirFilteredNorm, x = round(seurat.object.list[[i]]@assays[["SCT"]]@data, digits = 4), gene.type="SCTransform_Gene_Expression", overwrite=T, type="sparse", version="3")
  }
  
  StopWatchEnd$SaveFilteredData[[dataset]]  <- Sys.time()
  
}

StopWatchEnd$SCTransform  <- Sys.time()

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SAVE THE R_OBJECTS AND LOG FILES
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** SAVE THE R_OBJECTS AND LOG FILES ****\n")

####################################
### Saving each dataset R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving each dataset R object ***\n")
  
  StopWatchStart$SaveRDSEachDataset  <- Sys.time()
  
  for (i in 1:length(seurat.object.list)) {
    
    dataset <- rownames(InputsTable)[[i]]
    print(dataset)
    
    if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
      OutfileRDS<-paste0("R_OBJECTS_CWL/", list_DatasetMinioIDs[[dataset]], ".", PrefixOutfiles, ".", ProgramOutdir, "_", dataset , "_QC_Normalization.rds")
    } else {
      OutfileRDS<-paste0(Tempdir, "/R_OBJECTS/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset , "_QC_Normalization.rds")
    }
    print(OutfileRDS)
    saveRDS(seurat.object.list[[i]], file = OutfileRDS)
  }
  
  StopWatchEnd$SaveRDSEachDataset  <- Sys.time()
  
}else{
  
  writeLines("\n*** Not saving the full R object ***\n")
  
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_QC_Normalization_CPUtimes.txt")
write(file = OutfileCPUusage, x = paste("Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
write(file = OutfileCPUusage, x = paste("MaxGlobalVariables", MaxGlobalVariables, sep = "\t", collapse = ""))

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
    for (substep in rownames(InputsTable)) {
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

writeLines("\n*** Moving outfiles into outdir or keeping them at tempdir (if using CWL) ***\n")

### using two steps to copy files (`file.copy` and `file.remove`) instead of just `file.rename` to avoid issues with path to Tempdir in cluster systems
if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(Tempdir)
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    if (DirName == "FILTERED_DATA_MATRICES" | DirName == "UNFILTERED_DATA_MATRICES") {
      sapply(list.dirs(TempdirWithData, full.names = F, recursive = F), FUN=function(SubDirName) {
        sapply(list.dirs(paste0(TempdirWithData, "/", SubDirName), full.names = F, recursive = F), FUN=function(SubSubDirName) {
          OutdirFinal <- gsub(pattern = Tempdir, replacement =  paste0(Outdir, "/", ProgramOutdir), x = paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName))
          dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
          sapply(list.files(paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName), pattern = ".gz", full.names = F), FUN=function(EachFileName) {
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
writeLines(paste0("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report"))

quit()
