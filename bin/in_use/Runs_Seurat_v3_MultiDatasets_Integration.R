####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/v3.2/sctransform_vignette.html (SCtransform normalization)
### https://satijalab.org/seurat/v3.2/integration.html (general integration)
### https://satijalab.org/seurat/v3.2/immune_alignment.html (control vs. treatment)
### https://carmonalab.github.io/STACAS/tutorial.html (alternative anchor finder)
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads each normalized dataset R object, presumably from script `Runs_Seurat_v3_MultiDatasets_QC_Normalization.R'
### 2) Merges and integrates datasets correcting batch effects
### 3) Saves integrated datasets R object
### 4) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_Seurat_v3_MultiDatasets_Integration.R -h'
### for help
####################################

####################################
### THINGS TO DO:
### 1) Using pseudo-bulk correlation between datasets get dataset clusters and get representatives to be used as references
###    How does it look compared with STACAS distances?
####################################

####################################
### COMMENTS ON USING REFERENCE DATASETS VS. ALL PAIRWISE COMPARISONS TO FIND ANCHORS
### Using 4 datasets with 3K to 5K compared runs using either `-k NA` or `-k 1,2`
### Computing time was reduced in step IntegrateData() when using `-k 1,2`. Whereas FindIntegrationAnchors() was similar:
### Step/time(minutes)      Using_-k_NA   Using_-k_1,2
### FindIntegrationAnchors	3.336035	    3.118804
### IntegrateData           2.968116	    1.666707
### In both cases used `-u MAX -b 10000` in a 3.1-GHz Intel Core i5 CPU with 2 cores and 16 GB RAM
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
suppressPackageStartupMessages(library(STACAS))       # (GitHub carmonalab/STACAS) tested with v1.0.1 (compatible with Seurat v3.2.1). Needed for STACAS-based dataset integration
suppressPackageStartupMessages(library(tidyr))        # (CRAN) to handle tibbles and data.frames
suppressPackageStartupMessages(library(cluster))      # (CRAN) to cluster/sort the STACAS distances
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_Seurat_v3_MultiDatasets_Integration.R"
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
  make_option(c("-j", "--infile_r_objects"), default="NA",
              help="Path/name to a <tab> delimited file with dataset ID's and path/name to the QC and Normalization R/Seurat objects, like:
                d1  /path/to/d1_Normalization.rds
                d2  /path/to/d2_Normalization.rds

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-y", "--anchors_function"), default="Seurat",
              help="Indicates function to find integration anchors:
                'STACAS' = use function FindAnchors.STACAS() from library(STACAS)
                'Seurat' = use function FindIntegrationAnchors() from library(Seurat)
                
                Default = 'Seurat'"),
  #
  make_option(c("-z", "--reference_datasets"), default="NA",
              help="Indicates either of two options to use reference datasets for integration anchors:
                '<d1,d2,...etc>' = a <comma> delimited list of dataset ID(s), corresponding to --inputs_list column 1, to be used as references
                                   to obtain anchors (e.g. datasets expected to cover most expected cell types, with better QC, etc.)
                'NA'             = to compare all datatasets against each other to obtain anchors (i.e. no reference/s are used)
                                   Note this increases computing time and memory compared to using references.

                Note: If option '<d1,d2,...etc>' is used, anchors are first found between each non-reference and each reference dataset.
                      The references are then integrated through pairwise integration. Each non-reference is then mapped to the integrated reference.
                
                Default = 'NA'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
              
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
              
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file:
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                
                Default = '10'"),
  #
  make_option(c("-n", "--k_filter"), default="200",
              help="Integrating highly heterogenous datasets can lead to a too small number of anchors,
                resulting in an error 'Cannot find more nearest neighbours than there are points'
                This can be avoided by using a -k_filter value smaller than the default (e.g. '150')
                
                Default = '200'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system
              
                Default = 'MAX'"),
  #
  make_option(c("-s", "--save_r_object"), default="Y",
              help="Indicates if a R object with the data and analyzes from the run should be saved
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
InputRObjects           <- opt$infile_r_objects
AnchorsFunction         <- opt$anchors_function
ReferenceDatasets       <- opt$reference_datasets
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
KFilter                 <- as.numeric(opt$k_filter)
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
    "LOG_FILES",
    "R_OBJECTS",
    "PSEUDO_BULK",
    "STACAS"
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
    "LOG_FILES",
    "R_OBJECTS",
    "PSEUDO_BULK",
    "STACAS"
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

DefaultParameters <- list(

  ### Parameters for STACAS
  StacasVarGenesIntegratedN = 500,
  DigitsForRoundMedianDist = 3,
  
  ### Parameters for plots
  BaseSizeSinglePlotPdf  = 7,

  ### Parameters for datasets integration
  IntegrationNFeatures = 3000,
  
  ### Parameters for datasets comparison
  AssaysForPseudoBulk = c("RNA", "SCT")
)

### Assay types for plot and table outfiles
listAssaySuffixForOutfiles <- list(RNA="RNA", SCT="SCT", integrated="INT")

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("infiles_list", "infile_r_objects", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
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

DatasetIds <-list()
for (i in 1:nrow(InputsTable)) {
  DatasetIds[[i]] <- rownames(InputsTable)[[i]]
}

####################################
### Load each dataset R object
####################################
writeLines("\n*** Load each dataset R object ***\n")

StopWatchStart$LoadRDSEachDataset  <- Sys.time()

InputRObjects.tab <- read.table(InputRObjects, header = F, row.names = 1, stringsAsFactors = F)
colnames(InputRObjects.tab)<-c("PathToRObject")

seurat.object.list <- list()
for (dataset in rownames(InputRObjects.tab)) {
  print(dataset)
  DatasetIndexInInputsTable <- which(x = rownames(InputsTable) == dataset)
  InputRobject <- InputRObjects.tab[dataset,"PathToRObject"]
  seurat.object.list[[DatasetIndexInInputsTable]] <- readRDS(InputRobject)
}

StopWatchEnd$LoadRDSEachDataset  <- Sys.time()

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO INTEGRATE DATASETS
################################################################################################################################################
################################################################################################################################################

####################################
### Get correlation between datasets using pseudo-bulk
####################################
writeLines("\n*** Get correlation between datasets using pseudo-bulk ***\n")

StopWatchStart$GetCorrelBetweenDatasetsUsingPseudoBulk  <- Sys.time()

SeuratObjectsFilteredAndNormalized <- seurat.object.list
FirstSeuratObject   <- SeuratObjectsFilteredAndNormalized[[1]]
RestOfSeuratObjectsFiltered <- SeuratObjectsFilteredAndNormalized[c(2:NumberOfDatasets)]
RestOfDatasetsIds    <- unlist(DatasetIds[c(2:NumberOfDatasets)])

seurat.object.merged.normalized <- merge(FirstSeuratObject, y = RestOfSeuratObjectsFiltered, 
                                         add.cell.ids = DatasetIds,
                                         project = PrefixOutfiles)

seurat.object.list.normalized <- SplitObject(seurat.object.merged.normalized, split.by = "dataset")

for (assay_expression in DefaultParameters$AssaysForPseudoBulk) {
  
  ### Get pseudo-bulk matrices
  mat_for_correl_all_cells.df <- data.frame(row.names = rownames(seurat.object.list.normalized[[1]]@assays[[assay_expression]]))
  for (dataset in rownames(InputsTable)) {
    mat_for_correl_all_cells.df[[dataset]] <- rowSums(as.matrix(seurat.object.list.normalized[[dataset]]@assays[[assay_expression]][,]))
  }
  
  OutfilePseudoBulk <- paste0(Tempdir, "/PSEUDO_BULK/", PrefixOutfiles, ".", ProgramOutdir, "_PseudoBulk_EachDataset_", listAssaySuffixForOutfiles[[assay_expression]], ".tsv")
  Headers<-paste(paste0("PseudoBulk_EachDataset_", listAssaySuffixForOutfiles[[assay_expression]]), paste(colnames(mat_for_correl_all_cells.df), sep = "\t", collapse = "\t"), sep = "\t", collapse = "")
  write.table(Headers, file = OutfilePseudoBulk, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(mat_for_correl_all_cells.df,  file = OutfilePseudoBulk, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  ### Get correlation
  mat_for_correl_all_cells.cor <- round(cor(mat_for_correl_all_cells.df), digits = 3)
  OutfilePseudoBulkCor <- paste0(Tempdir, "/PSEUDO_BULK/", PrefixOutfiles, ".", ProgramOutdir, "_PseudoBulk_EachDataset_", listAssaySuffixForOutfiles[[assay_expression]], "_cor", ".tsv")
  Headers<-paste(paste0("PseudoBulk_EachDataset_", listAssaySuffixForOutfiles[[assay_expression]], "_cor"), paste(colnames(mat_for_correl_all_cells.cor), sep = "\t", collapse = "\t"), sep = "\t", collapse = "")
  write.table(Headers, file = OutfilePseudoBulkCor, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(mat_for_correl_all_cells.cor,  file = OutfilePseudoBulkCor, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
}

StopWatchEnd$GetCorrelBetweenDatasetsUsingPseudoBulk  <- Sys.time()

####################################
### Get reference datasets
####################################
writeLines("\n*** Get reference datasets ***\n")

StopWatchStart$GetReferenceDatasets <- Sys.time()

if (regexpr("^NA$", ReferenceDatasets , ignore.case = T)[1] == 1) {
  writeLines("\n*** Will compare all-vs-all datasets to get anchors ***\n")
  ReferenceDatasets.indices <- c(1:nrow(InputsTable))
}else{
  writeLines("\n*** Determine reference dataset indices ***\n")
  ReferenceDatasets.list <- unlist(strsplit(ReferenceDatasets, ","))
  NumberOfFoundReferenceDatasetIDs <- sum(ReferenceDatasets.list %in% rownames(InputsTable) == T)
  if (NumberOfFoundReferenceDatasetIDs == length(ReferenceDatasets.list)) {
    ReferenceDatasets.indices <- match(ReferenceDatasets.list, rownames(InputsTable))
  }else{
    stop(paste0("Requested ", length(ReferenceDatasets.list), " datasets as references by parameter `-z`, but found ", NumberOfFoundReferenceDatasetIDs, " in --inputs_list row headers"))
  }
  print(paste0("Will use datasets: ", paste(as.character(ReferenceDatasets.indices), sep = "", collapse = ",")))
}

StopWatchEnd$GetReferenceDatasets <- Sys.time()

####################################
### Get integration anchors
####################################
writeLines("\n*** Get integration anchors ***\n")

if (regexpr("^STACAS$", AnchorsFunction , ignore.case = T)[1] == 1) {
  
  ####################################
  ### Get anchors with STACAS
  ####################################
  writeLines("\n*** Get anchors with STACAS ***\n")
  
  StopWatchStart$FindIntegrationAnchorsStacas  <- Sys.time()
  
  seurat.object.anchors.unfiltered <- FindAnchors.STACAS(object.list = seurat.object.list, dims=PcaDimsUse, anchor.features=DefaultParameters$StacasVarGenesIntegratedN, reference = ReferenceDatasets.indices, verbose = T)
  
  StopWatchEnd$FindIntegrationAnchorsStacas  <- Sys.time()
  
  ####################################
  ### Get STACAS distances between dataset anchors
  ####################################
  writeLines("\n*** Get STACAS distances between dataset anchors ***\n")
  
  StopWatchStart$GetDistancesStacas <- Sys.time()
  
  DensityPlots <- PlotAnchors.STACAS(seurat.object.anchors.unfiltered, obj.names=names(seurat.object.list))
  
  OutfilesTsv <- paste0(Tempdir, "/STACAS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_distances.tsv")
  Headers<-paste("dataset1", "dataset2", "cell1", "cell2", "score", "dist1.2", "dist2.1", "dist.mean", sep = "\t", collapse = "")
  write.table(Headers, file = OutfilesTsv, row.names = F, col.names = F, sep="\t", quote = F)
  lapply(1:length(seurat.object.list), function(DatasetNumber) {
    write.table(file = OutfilesTsv, x = DensityPlots[[DatasetNumber]]$data[,c("dataset1", "dataset2", "cell1", "cell2", "score", "dist1.2", "dist2.1", "dist.mean")],
                quote = F, sep = "\t", row.names = F, col.names = F, append = T)
  })
  
  StopWatchEnd$GetDistancesStacas <- Sys.time()
  
  ####################################
  ### Get median of STACAS dist.mean between dataset anchors
  ####################################
  writeLines("\n*** Get median of STACAS dist.mean between dataset anchors ***\n")
  
  StopWatchStart$GetMedianDistMeanStacas <- Sys.time()
  
  Distances.df <- read.table(file = OutfilesTsv, header = T, sep = "\t", row.names = NULL)
  Distances.df$datasets <- paste0(Distances.df[,"dataset1"],"---",Distances.df[,"dataset2"])
  MedianDistances.df <- aggregate(Distances.df[,"dist.mean"], list(Distances.df$datasets), median)
  SplitDatasets <- data.frame(do.call('rbind', strsplit(as.character(MedianDistances.df[,1]),'---', fixed=TRUE)))
  MedianDistances.df <- cbind(MedianDistances.df, SplitDatasets)
  colnames(MedianDistances.df) <- c("datasets","median","dataset1","dataset2")
  MedianDistances.df<-MedianDistances.df[,c("dataset1","dataset2","median")]
  MedianDistances.df<-data.frame(pivot_wider(MedianDistances.df, id_cols = dataset1, names_from = dataset2, values_from = median))
  rownames(MedianDistances.df) <- MedianDistances.df[,1]
  MedianDistances.mat <- as.matrix(MedianDistances.df[,-1])
  MedianDistances.mat[is.na(MedianDistances.mat)] <- 1
  MedianDistances.clust <- agnes(x = 1-MedianDistances.mat, metric = "manhattan")
  MedianDistances.order <- rownames(MedianDistances.mat)[MedianDistances.clust$order]
  MedianDistances.mat   <- MedianDistances.mat[MedianDistances.order,MedianDistances.order]
  
  OutfilesTsv <- paste0(Tempdir, "/STACAS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_inv_dist.mean_median_ordered.tsv")
  write.table(data.frame("InvDistMean_median"=MedianDistances.order, round(MedianDistances.mat,DefaultParameters$DigitsForRoundMedianDist)), OutfilesTsv, row.names = F,sep="\t",quote = F)
  OutfilesTsv <- paste0(Tempdir, "/STACAS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_dist.mean_median_ordered.tsv")
  write.table(data.frame("DistMean_median"=MedianDistances.order, round(1-MedianDistances.mat,DefaultParameters$DigitsForRoundMedianDist)), OutfilesTsv, row.names = F,sep="\t",quote = F)
  
  StopWatchEnd$GetMedianDistMeanStacas <- Sys.time()
  
  ####################################
  ### Get list of STACAS anchors
  ####################################
  writeLines("\n*** Get list of STACAS anchors ***\n")
  
  StopWatchStart$GetListAnchorsStacas <- Sys.time()
  
  TableAnchorNumbers <- function(so.anchors,sos.list,OutfilesTsv,OutfilePdf,HeaderPlot) {
    anchor.stats<-table(so.anchors@anchors[,c("dataset1","dataset2")])
    rownames(anchor.stats) <- names(sos.list)
    colnames(anchor.stats) <- names(sos.list)
    
    write(x=paste0("Dataset", "\t", paste(names(sos.list), sep = "\t", collapse = "\t")), file=OutfilesTsv)
    write.table(x=anchor.stats, file=OutfilesTsv, row.names = T, col.names = F, sep="\t", quote = F, append = T)
    
    pdf(file=OutfilePdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = DefaultParameters$BaseSizeSinglePlotPdf)
    anchor.stats.norm <- apply(anchor.stats, 1, function(x) {x/sum(x)})
    toplot <- melt(anchor.stats.norm, varnames=c("Dataset1","Dataset2"), value.name = "Fraction")
    print(ggplot(toplot, aes(fill=Dataset1, y=Fraction, x=Dataset2)) +
            geom_bar(position="stack", stat="identity") +
            theme(axis.text.x = element_text(angle = 45)) +
            ggtitle(paste0("Fraction of anchors between data sets - ", HeaderPlot, " anchor filtering"))
    )
    dev.off()
  }
  
  StopWatchEnd$GetListAnchorsStacas <- Sys.time()
  
  ####################################
  ### Get all anchors
  ####################################
  writeLines("\n*** Get all anchors ***\n")
  
  StopWatchStart$GetAllAnchorsStacas <- Sys.time()
  
  OutfilesTsv <- paste0(Tempdir, "/STACAS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_all_anchor_numbers.tsv")
  OutfilePdf  <- paste0(Tempdir, "/STACAS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_all_anchor_numbers.pdf")
  so.anchors  <- seurat.object.anchors.unfiltered
  sos.list    <- seurat.object.list
  TableAnchorNumbers(so.anchors,sos.list,OutfilesTsv,OutfilePdf,"BEFORE")
  
  StopWatchEnd$GetAllAnchorsStacas <- Sys.time()
  
  ####################################
  ### Get filtered anchors
  ####################################
  writeLines("\n*** Get filtered anchors ***\n")
  
  StopWatchStart$GetFilteredAnchorsStacas <- Sys.time()
  
  seurat.object.anchors <- FilterAnchors.STACAS(seurat.object.anchors.unfiltered)
  print(seurat.object.anchors)
  
  OutfilesTsv <- paste0(Tempdir, "/STACAS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_filtered_anchor_numbers.tsv")
  OutfilePdf  <- paste0(Tempdir, "/STACAS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_filtered_anchor_numbers.pdf")
  so.anchors  <- seurat.object.anchors
  sos.list    <- seurat.object.list
  TableAnchorNumbers(so.anchors,sos.list,OutfilesTsv,OutfilePdf,"AFTER")
  
  StopWatchEnd$GetFilteredAnchorsStacas <- Sys.time()
  
  ####################################
  ### Get optimal integration tree
  ####################################
  writeLines("\n*** Get optimal integration tree ***\n")
  
  StopWatchStart$GetOptimalIntegrationTreeStacas <- Sys.time()
  
  all.genes <- row.names(seurat.object.list[[1]])
  lapply(2:length(seurat.object.list), function(DatasetNumber) {
    all.genes <- intersect(all.genes, row.names(seurat.object.list[[DatasetNumber]]))
  })
  
  SampleTree <- SampleTree.STACAS(seurat.object.anchors)
  
  StopWatchEnd$GetOptimalIntegrationTreeStacas <- Sys.time()
  
} else if (regexpr("^Seurat$", AnchorsFunction , ignore.case = T)[1] == 1) {
  
  ####################################
  ### Get anchors with Seurat
  ####################################
  writeLines("\n*** Get anchors with Seurat ***\n")
  
  StopWatchStart$SelectIntegrationFeatures  <- Sys.time()
  
  writeLines("\n*** Run SelectIntegrationFeatures() ***\n")
  seurat.object.integratedfeatures <- SelectIntegrationFeatures(object.list = seurat.object.list, nfeatures = DefaultParameters$IntegrationNFeatures)
  
  StopWatchEnd$SelectIntegrationFeatures  <- Sys.time()
  
  StopWatchStart$PrepSCTIntegration  <- Sys.time()
  
  writeLines("\n*** Run PrepSCTIntegration() ***\n")
  seurat.object.list <- PrepSCTIntegration(object.list = seurat.object.list, anchor.features = seurat.object.integratedfeatures, verbose = T)
  
  StopWatchEnd$PrepSCTIntegration  <- Sys.time()
  
  StopWatchStart$FindIntegrationAnchors  <- Sys.time()
  
  writeLines("\n*** Run FindIntegrationAnchors() ***\n")
  print(paste0("Using k.filter = ", KFilter))
  seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, k.filter = KFilter, normalization.method = "SCT", anchor.features = seurat.object.integratedfeatures, reference = ReferenceDatasets.indices, verbose = T)
  SampleTree <- NULL
  
  StopWatchEnd$FindIntegrationAnchors  <- Sys.time()
  
  print(seurat.object.anchors)
  
}

####################################
### Integrating datasets
####################################
writeLines("\n*** Integrating datasets ***\n")

StopWatchStart$IntegrateData  <- Sys.time()

writeLines("\n*** Run IntegrateData() ***\n")
seurat.object.integrated <- IntegrateData(anchorset = seurat.object.anchors, normalization.method = "SCT", sample.tree = SampleTree, preserve.order = T, verbose = T)

StopWatchEnd$IntegrateData  <- Sys.time()

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
  
  OutfileRDS<-paste0(Tempdir, "/R_OBJECTS/", PrefixOutfiles, ".", ProgramOutdir, "_Integration", ".rds")
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
    OutdirFinal <- paste0(Outdir, "/", ProgramOutdir, "/", DirName)
    print(OutdirFinal)
    dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
    sapply(list.files(TempdirWithData, pattern = paste0("^", PrefixOutfiles, ".", ProgramOutdir), full.names = F), FUN=function(EachFileName) {
      file.copy(from=paste0(TempdirWithData, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
      file.remove(paste0(TempdirWithData, "/", EachFileName))
    })
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
