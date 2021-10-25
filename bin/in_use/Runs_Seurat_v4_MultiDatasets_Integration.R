####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/articles/integration_introduction.html
### https://carmonalab.github.io/STACAS/tutorial.html (alternative anchor finder)
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads each normalized dataset R object produced by script `Runs_Seurat_v4_MultiDatasets_QC_Normalization.R'
### 2) Merges and integrates datasets correcting batch effects
### 3) Saves integrated datasets R object
### 4) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT
### For help using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_Seurat_v4_MultiDatasets_Integration.R -h'
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
writeLines("\n**** LOAD REQUIRED LIBRARIES ****\n")

suppressPackageStartupMessages(library(Seurat))       # (CRAN) tested with v4.0.2. To run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # (CRAN) to calculate the percentages of extra properties to be plotted
suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(stringr))      # (CRAN) to regex and extract matching string. Only needed if using `-w Y` and `-x`.
suppressPackageStartupMessages(library(STACAS))       # (GitHub carmonalab/STACAS) tested with v1.1.0 (compatible with Seurat v4.0.2). Needed for STACAS-based dataset integration
suppressPackageStartupMessages(library(tidyr))        # (CRAN) to handle tibbles and data.frames
suppressPackageStartupMessages(library(cluster))      # (CRAN) to cluster/sort the STACAS distances
####################################

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SETUP RUN
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** SETUP RUN ****\n")

####################################
### Turning warnings off for the sake of a cleaner output
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_Seurat_v4_MultiDatasets_Integration.R"
ProgramOutdir  <- "SEURAT"

####################################
### Get inputs from command line argumets
####################################
#####

option_list <- list(
  make_option(c("-i", "--inputs_list"), default="NA",
              help="Path/name to a <tab> delimited file with one dataset per row and the following 3 columns specifying details for each dataset:
                1) unique dataset ID (e.g. 'd1')
                2) /path_to/dataset_R_object
                3) dataset type (e.g. 'control' or 'treatment') or use 'type' to skip using dataset types

                Notes:
                (a) The order of the list of datasets in --inputs_list may influence the results, including number of clusters,
                t-SNE/UMAP and differentially expressed genes. List datasets better measured first.
                
                (b) middle dashes '-' in columns 1 or 3 will be replaced by low dashes '_'

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-j", "--infile_r_objects"), default="NA",
              help="Used with the --run_cwl to mount the --inputs_list R object into CWL containers

                Default = 'NA'"),
  #
  make_option(c("-y", "--anchors_function"), default="Seurat_CCA",
              help="Indicates function to find integration anchors:
                'STACAS'      = use function FindAnchors.STACAS() from library(STACAS)
                'Seurat_CCA'  = use function FindIntegrationAnchors(..., reduction = 'cca') from library(Seurat)
                'Seurat_RPCA' = use function FindIntegrationAnchors(..., reduction = 'rpca') from library(Seurat)
                
                Default = 'Seurat_CCA'"),
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
  make_option(c("-d", "--pca_dimensions_anchors"), default="30",
              help="Max value of PCA dimensions to use for anchor detection
                FindIntegrationAnchors(..., dims = 1:-d)

                Default = '30'"),
  #
  make_option(c("-n", "--k_filter"), default="200",
              help="Integrating highly heterogenous datasets can lead to a too small number of anchors,
                resulting in an error 'Cannot find more nearest neighbours than there are points'
                This can be avoided by using a -k_filter value smaller than the default (e.g. '150')
                
                Default = '200'"),
  #
  make_option(c("-q", "--dist_thr"), default="0.8",
              help="Only needed if using '-y STACAS'
                Distance threshold for anchor filtering. Distances are calculated by RPCA between pairs of datasets.
                If not specified, dist_thr defaults to the dist.pct percentile of the distance between the two closest datasets.
                
                Default = '0.8'"),
  #
  make_option(c("-r", "--k_weight"), default="100",
              help="Number of neighbors to consider when weighting anchors
                
                Default = '100'"),
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
  make_option(c("-w", "--run_cwl"), default="0",
              help="Indicates if this script should produce 'frontend' files for crescent.cloud and if CWL is used
                0 = no frontend files should be produced and CWL is not used
                1 = frontend files should be produced and CWL is used
                2 = frontend files should be produced but CWL is not used (e.g. run locally)

                Default = '0'"),
  #
  make_option(c("-x", "--minio_path"), default="NA",
              help="Only needed if using '-w 1' to mount input data files in 'inputs_list' to CWL containers
              
                Default = 'NA'"),
  #
  make_option(c("-a", "--max_global_variables"), default="10000",
              help="Indicates maximum allowed total size (in bytes) of global variables identified. Used by library(future) to prevent too large exports
              
                Default = '10000' for 10000 MiB")
)

#####

opt <- parse_args(OptionParser(option_list=option_list))

InputsList              <- opt$inputs_list
InfileRObjects          <- opt$infile_r_objects
AnchorsFunction         <- opt$anchors_function
AnchorsInFile           <- opt$anchors_infile
SampleTreeInFile        <- opt$sample_tree_infile
ReferenceDatasets       <- opt$reference_datasets
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions_anchors))
KFilter                 <- as.numeric(opt$k_filter)
DistThr                 <- as.numeric(opt$dist_thr)
KWeight                 <- as.numeric(opt$k_weight)
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- as.numeric(opt$run_cwl)
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
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("inputs_list", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
  }
}

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

FILE_TYPE_OUT_DIRECTORIES = c(
  "LOG_FILES",
  "ANCHORS"
)

if (RunsCwl == 0 || RunsCwl == 2) {
  FILE_TYPE_OUT_DIRECTORIES = c(
    FILE_TYPE_OUT_DIRECTORIES,
    "R_OBJECTS"
  )
}

if (RunsCwl == 1) {
  ### Using `-w 1` will make Tempdir, which takes the value of ProgramOutdir, and it will be the final out-directory
  ### for most outfiles, except R objects, which will be written into R_OBJECTS_CWL
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F)
  dir.create(file.path("R_OBJECTS_CWL"), showWarnings = F)
}else if (RunsCwl == 0 || RunsCwl == 2) {
  ## Using `-w 0` or `-w 2` will create a Tempdir/DIRECTORY for temporary storage of outfiles because sometimes
  ## long paths of outdirectories casuse R to leave outfiles unfinished
  ## 'DIRECTORY' is one of the directories specified at FILE_TYPE_OUT_DIRECTORIES
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  ## The difference between `-w 0` and `-w 2` is that the first one doesn't produce frontend outfiles for crescent.cloud
  Tempdir           <- "~/temp"
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
}else{
  stop(paste0("ERROR unexpected option '-w ", RunsCwl))
}

sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(eachdir) {
  dir.create(file.path(paste0(Tempdir, "/", eachdir)), showWarnings = F, recursive = T)
})

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

StopWatchStart$ReportUsedOptions  <- Sys.time()

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_Integration_UsedOptions", ".txt")

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

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_Integration_RSessionInfo", ".txt")
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

### Allocation of global variables allows to handle memory errors
### 4GiB (4000*1024^2) are assigned as: `options(future.globals.maxSize = 4000 * 1024^2)`
options(future.globals.maxSize = MaxGlobalVariables * 1024^2)

####################################
### Define default parameters
####################################
### Some of these default parameters are provided by Seurat developers, others are tailored empirically

DefaultParameters <- list(

  ### Parameters for STACAS
  StacasVarGenesIntegratedN = 800,
  StacasNGenes = 500,
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
### Check anchor parameters are compatible with each other
####################################
writeLines("\n*** Check anchor parameters are compatible with each other ***\n")

if (regexpr("^STACAS$|^Seurat_CCA$|^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {

  if (regexpr("^Seurat_CCA$", AnchorsFunction , ignore.case = T)[1] == 1) {
    ReductionForFindIntegrationAnchors <- "cca"
  } else if (regexpr("^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {
    ReductionForFindIntegrationAnchors <- "rpca"
  }

} else {
  stop(paste("ERROR!!! parameter -y must be 'STACAS', Seurat_CCA' or 'Seurat_RPCA'"))
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO LOAD DATASETS
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** LOAD DATASETS ****\n")

####################################
### Load --inputs_list
####################################
if (regexpr("^STACAS$|^Seurat_CCA$|^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {

  writeLines("\n*** Load --inputs_list ***\n")
  
  if (RunsCwl == 0 || RunsCwl == 2) {
    InputsList<-gsub("^~/",paste0(UserHomeDirectory,"/"), InputsList)
    InputsTable<-read.table(InputsList, header = F, row.names = 1, stringsAsFactors = F)
    colnames(InputsTable)<-c("PathToRObject","DatasetType")
    
  }else{
    MinioPaths <- as.list(strsplit(MinioPath, ",")[[1]])
    MinioDataPaths = data.frame(dataset_ID=rep(0, length(MinioPaths)), dataset_path=rep(0, length(MinioPaths)))
    
    for (i in seq_along(MinioPaths)) {
      MinioDataPaths[i, ] = c(basename(MinioPaths[[i]]), MinioPaths[[i]])
    }
    
    InputsTable0 <- read.table(InputsList, header = T, sep = ",", stringsAsFactors = F)
    
    MergedInputsTable <- merge(MinioDataPaths, InputsTable0, by="dataset_ID")
    MergeFilter <- c("name", "dataset_ID","dataset_path", "dataset_type")
    MergedInputsTableFiltered <- MergedInputsTable[MergeFilter]
    MergedInputsTableFilteredFinal <- MergedInputsTableFiltered[,-1]
    rownames(MergedInputsTableFilteredFinal) <- MergedInputsTableFiltered[,1]
    colnames(MergedInputsTableFilteredFinal) <-c("DatasetMinioID","PathToRObject","DatasetType")
    
    InputsTable <- MergedInputsTableFilteredFinal
  }
  
  ##### Replace low dashes by dots in rownames(InputsTable) or DatasetType
  rownames(InputsTable) <- gsub(x=rownames(InputsTable), pattern = "-", replacement = "_")
  InputsTable[,"DatasetType"] <- gsub(x=InputsTable[,"DatasetType"], pattern = "-", replacement = "_")
  NumberOfDatasets <- nrow(InputsTable)
  
  DatasetIds <-list()
  for (i in 1:nrow(InputsTable)) {
    DatasetIds[[i]] <- rownames(InputsTable)[[i]]
  }
}

####################################
### Load each dataset R object
####################################
if (regexpr("^STACAS$|^Seurat_CCA$|^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {
  writeLines("\n*** Load each dataset R object ***\n")
  
  StopWatchStart$LoadRDSEachDataset  <- Sys.time()

  seurat.object.list <- list()
  
  if (RunsCwl == 1) {
    RObjects <- list.files(InfileRObjects, pattern="*_QC_Normalization.rds", full.names=T)
    for (object in RObjects) {
      MinioID <- str_extract(basename(object), regex("[^.]*"))
      for (dataset in rownames(InputsTable)) {
        DatasetMinioID <- InputsTable[dataset,"DatasetMinioID"]
        if (MinioID == DatasetMinioID) {
          DatasetIndexInInputsTable <- which(x = rownames(InputsTable) == dataset)
          seurat.object.list[[DatasetIndexInInputsTable]] <- readRDS(object)
        }
      }
    }
  }else{
    for (dataset in rownames(InputsTable)) {
      print(dataset)
      DatasetIndexInInputsTable <- which(x = rownames(InputsTable) == dataset)
      InputRobject <- InputsTable[dataset,"PathToRObject"]
      seurat.object.list[[DatasetIndexInInputsTable]] <- readRDS(InputRobject)
      print(seurat.object.list[[DatasetIndexInInputsTable]])
    }
  }
  
  StopWatchEnd$LoadRDSEachDataset  <- Sys.time()
  
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO INTEGRATE DATASETS
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** INTEGRATE DATASETS ****\n")

####################################
### Merge datasets
####################################
if (regexpr("^STACAS$|^Seurat_CCA$|^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {
  writeLines("\n*** Merge datasets ***\n")
  
  StopWatchStart$MergeDatasets  <- Sys.time()
  
  FirstSeuratObject   <- seurat.object.list[[1]]
  RestOfSeuratObjectsFiltered <- seurat.object.list[c(2:NumberOfDatasets)]
  RestOfDatasetsIds    <- unlist(DatasetIds[c(2:NumberOfDatasets)])
  
  seurat.object.merged.normalized <- merge(FirstSeuratObject, y = RestOfSeuratObjectsFiltered, 
                                           add.cell.ids = DatasetIds,
                                           project = PrefixOutfiles)
  
  seurat.object.list.normalized <- SplitObject(seurat.object.merged.normalized, split.by = "dataset")

  StopWatchEnd$MergeDatasets  <- Sys.time()
}

for (i in 1:length(seurat.object.list.normalized)) {
  dataset <- rownames(InputsTable)[[i]]
  print(dataset)
  seurat.object.list.normalized[[i]] <- SCTransform(seurat.object.list.normalized[[i]], verbose = T)
}

####################################
### Get reference datasets
####################################
if (regexpr("^STACAS$|^Seurat_CCA$|^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {
  writeLines("\n*** Get reference datasets ***\n")

  StopWatchStart$GetReferenceDatasets <- Sys.time()
  
  if (regexpr("^NA$", ReferenceDatasets , ignore.case = T)[1] == 1) {
    writeLines("\n*** Will compare all-vs-all datasets to get anchors ***\n")
    ReferenceDatasets.indices <- c(1:nrow(InputsTable))
  }else{

    if (RunsCwl == 1) {
      InputsTableReferenceID <- InputsTable$DatasetMinioID
    } else {
      InputsTableReferenceID <- rownames(InputsTable)
    }
    
    writeLines("\n*** Determine reference dataset indices ***\n")
    ReferenceDatasets.list <- unlist(strsplit(ReferenceDatasets, ","))
    NumberOfFoundReferenceDatasetIDs <- sum(ReferenceDatasets.list %in% InputsTableReferenceID == T)
    if (NumberOfFoundReferenceDatasetIDs == length(ReferenceDatasets.list)) {
      ReferenceDatasets.indices <- match(ReferenceDatasets.list, InputsTableReferenceID)
    }else{
      stop(paste0("Requested ", length(ReferenceDatasets.list), " datasets as references by parameter `-z`, but found ", NumberOfFoundReferenceDatasetIDs, " in --inputs_list row headers"))
    }
    print(paste0("Will use datasets: ", paste(as.character(ReferenceDatasets.indices), sep = "", collapse = ",")))
  }
  
  StopWatchEnd$GetReferenceDatasets <- Sys.time()
}

####################################
### Get integration anchors
####################################
writeLines("\n*** Get integration anchors ***\n")

if (regexpr("^STACAS$", AnchorsFunction , ignore.case = T)[1] == 1) {
  
  ####################################
  ### Get anchors with STACAS
  ####################################
  writeLines("\n*** Get anchors with STACAS ***\n")
  
  StopWatchStart$SelectIntegrationFeatures  <- Sys.time()
  
  writeLines("\n*** Run SelectIntegrationFeatures() ***\n")
  seurat.object.integratedfeatures <- SelectIntegrationFeatures(object.list = seurat.object.list.normalized,
                                                                nfeatures = DefaultParameters$IntegrationNFeatures)
  
  StopWatchEnd$SelectIntegrationFeatures  <- Sys.time()
  
  StopWatchStart$PrepSCTIntegration  <- Sys.time()
  
  writeLines("\n*** Run PrepSCTIntegration() ***\n")
  seurat.object.list.normalized <- PrepSCTIntegration(object.list = seurat.object.list.normalized,
                                           anchor.features = seurat.object.integratedfeatures,
                                           verbose = T)
  
  StopWatchEnd$PrepSCTIntegration  <- Sys.time()
  
  StopWatchStart$FindIntegrationAnchorsStacas  <- Sys.time()
  
  seurat.object.anchors.unfiltered <- FindAnchors.STACAS(object.list = seurat.object.list.normalized,
                                                         dims=PcaDimsUse,
                                                         anchor.features=DefaultParameters$StacasVarGenesIntegratedN,
                                                         normalization.method = "SCT",
                                                         reference = ReferenceDatasets.indices,
                                                         verbose = T)
  
    
  StopWatchEnd$FindIntegrationAnchorsStacas  <- Sys.time()
  
  ####################################
  ### Get STACAS distances between dataset anchors
  ####################################
  writeLines("\n*** Get STACAS distances between dataset anchors ***\n")
  
  StopWatchStart$GetDistancesStacas <- Sys.time()
  
  DensityPlots <- PlotAnchors.STACAS(seurat.object.anchors.unfiltered, obj.names=names(seurat.object.list.normalized))
  
  OutfileStacasDistances <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_distances.tsv.bz2")
  Outfile.con <- bzfile(OutfileStacasDistances, "w")
  Headers<-paste("dataset1", "dataset2", "cell1", "cell2", "score", "dist1.2", "dist2.1", "dist.mean", sep = "\t", collapse = "")
  write.table(Headers, file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F)
  lapply(1:length(seurat.object.list.normalized), function(DatasetNumber) {
    write.table(file = Outfile.con, x = DensityPlots[[DatasetNumber]]$data[,c("dataset1", "dataset2", "cell1", "cell2", "score", "dist1.2", "dist2.1", "dist.mean")],
                quote = F, sep = "\t", row.names = F, col.names = F, append = T)
  })
  close(Outfile.con)
  
  StopWatchEnd$GetDistancesStacas <- Sys.time()
  
  ####################################
  ### Get median of STACAS dist.mean between dataset anchors
  ####################################
  writeLines("\n*** Get median of STACAS dist.mean between dataset anchors ***\n")
  
  StopWatchStart$GetMedianDistMeanStacas <- Sys.time()
  
  Distances.df <- read.table(file = OutfileStacasDistances, header = T, sep = "\t", row.names = NULL)
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
  MedianDistances.mat   <- MedianDistances.mat[as.factor(MedianDistances.order),as.factor(MedianDistances.order)]
  OutfilePathName <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_inv_dist.mean_median_ordered.tsv.bz2")
  Outfile.con <- bzfile(OutfilePathName, "w")
  write.table(data.frame("InvDistMean_median"=MedianDistances.order, round(MedianDistances.mat,DefaultParameters$DigitsForRoundMedianDist)),
              file = OutfilePathName, row.names = F,sep="\t",quote = F)
  close(Outfile.con)
  OutfilePathName <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_dist.mean_median_ordered.tsv.bz2")
  Outfile.con <- bzfile(OutfilePathName, "w")
  write.table(data.frame("DistMean_median"=MedianDistances.order, round(1-MedianDistances.mat,DefaultParameters$DigitsForRoundMedianDist)),
              file = OutfilePathName, row.names = F,sep="\t",quote = F)
  close(Outfile.con)
  
  StopWatchEnd$GetMedianDistMeanStacas <- Sys.time()
  
  ####################################
  ### Get list of STACAS anchors
  ####################################
  writeLines("\n*** Get list of STACAS anchors ***\n")
  
  StopWatchStart$GetListAnchorsStacas <- Sys.time()
  
  TableAnchorNumbers <- function(seurat.object.anchors.final,seurat.object.list.normalized,OutfilesTsv,OutfilePdf,HeaderPlot) {
    anchor.stats<-table(seurat.object.anchors.final@anchors[,c("dataset1","dataset2")])
    rownames(anchor.stats) <- names(seurat.object.list.normalized)
    colnames(anchor.stats) <- names(seurat.object.list.normalized)
    
    OutfilePathName<-OutfilesTsv
    Outfile.con <- bzfile(OutfilePathName, "w")
    write(x=paste0("Dataset", "\t", paste(names(seurat.object.list.normalized), sep = "\t", collapse = "\t")), file=Outfile.con)
    write.table(x=anchor.stats, file=Outfile.con, row.names = T, col.names = F, sep="\t", quote = F, append = T)
    close(Outfile.con)
    
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
  
  OutfilesTsv <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_all_anchor_numbers.tsv.bz2")
  OutfilePdf  <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_all_anchor_numbers.pdf")
  TableAnchorNumbers(seurat.object.anchors.unfiltered,seurat.object.list.normalized,OutfilesTsv,OutfilePdf,"BEFORE")
  
  StopWatchEnd$GetAllAnchorsStacas <- Sys.time()
  
  ####################################
  ### Get filtered anchors
  ####################################
  writeLines("\n*** Get filtered anchors ***\n")
  
  StopWatchStart$GetFilteredAnchorsStacas <- Sys.time()
  
  seurat.object.anchors.final <- FilterAnchors.STACAS(seurat.object.anchors.unfiltered, dist.thr = DistThr)

  OutfilesTsv <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_filtered_anchor_numbers.tsv.bz2")
  OutfilePdf  <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_filtered_anchor_numbers.pdf")
  TableAnchorNumbers(seurat.object.anchors.final,seurat.object.list.normalized,OutfilesTsv,OutfilePdf,"AFTER")
  
  StopWatchEnd$GetFilteredAnchorsStacas <- Sys.time()
  
  ####################################
  ### Get optimal integration tree
  ####################################
  writeLines("\n*** Get optimal integration tree ***\n")
  
  StopWatchStart$GetOptimalIntegrationTreeStacas <- Sys.time()
  
  all.genes <- row.names(seurat.object.list.normalized[[1]])
  lapply(2:length(seurat.object.list.normalized), function(DatasetNumber) {
    all.genes <- intersect(all.genes, row.names(seurat.object.list.normalized[[DatasetNumber]]))
  })
  
  SampleTree <- SampleTree.STACAS(seurat.object.anchors.final)
  
  OutfileSampleTree <- paste0(Tempdir, "/ANCHORS/", PrefixOutfiles, ".", ProgramOutdir, "_STACAS_SampleTree", ".tsv.bz2")
  Outfile.con <- bzfile(OutfileSampleTree, "w")
  write.table(SampleTree, file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F)
  close(Outfile.con)
  
  StopWatchEnd$GetOptimalIntegrationTreeStacas <- Sys.time()
  
} else if (regexpr("^Seurat_CCA$|^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {
    
  ####################################
  ### Get anchors with Seurat
  ####################################
  writeLines("\n*** Get anchors with Seurat ***\n")
  
  StopWatchStart$SelectIntegrationFeatures  <- Sys.time()
  
  writeLines("\n*** Run SelectIntegrationFeatures() ***\n")
  seurat.object.integratedfeatures <- SelectIntegrationFeatures(object.list = seurat.object.list.normalized,
                                                                nfeatures = DefaultParameters$IntegrationNFeatures)
  
  StopWatchEnd$SelectIntegrationFeatures  <- Sys.time()
  
  StopWatchStart$PrepSCTIntegration  <- Sys.time()
  
  writeLines("\n*** Run PrepSCTIntegration() ***\n")
  seurat.object.list.normalized <- PrepSCTIntegration(object.list = seurat.object.list.normalized,
                                                      anchor.features = seurat.object.integratedfeatures,
                                                      verbose = T)
  
  StopWatchEnd$PrepSCTIntegration  <- Sys.time()
  
  if (regexpr("^Seurat_RPCA$", AnchorsFunction , ignore.case = T)[1] == 1) {
    ####################################
    ### Run RunPCA() on each dataset
    ####################################
    writeLines("\n*** Run RunPCA() on each dataset ***\n")
  
    StopWatchStart$RunPCAEachDataset  <- Sys.time()
  
    ### This is not needed in Seurat's CCA method
    seurat.object.list.normalized <- lapply(X = seurat.object.list.normalized, FUN = RunPCA, features = seurat.object.integratedfeatures)
  
    StopWatchEnd$RunPCAEachDataset  <- Sys.time()
  }

  StopWatchStart$FindIntegrationAnchors  <- Sys.time()
  
  writeLines("\n*** Run FindIntegrationAnchors() ***\n")
  print(paste0("Using k.filter = ", KFilter))
  seurat.object.anchors.final <- FindIntegrationAnchors(object.list = seurat.object.list.normalized,
                                                        k.filter = KFilter,
                                                        normalization.method = "SCT",
                                                        dims=PcaDimsUse,
                                                        anchor.features = seurat.object.integratedfeatures,
                                                        reduction = ReductionForFindIntegrationAnchors,
                                                        k.anchor = 20,
                                                        reference = ReferenceDatasets.indices,
                                                        verbose = T)

  SampleTree <- NULL
  
  StopWatchEnd$FindIntegrationAnchors  <- Sys.time()
  
} else {
  
  stop(paste0("ERROR!!! unexpected parameter -y ", AnchorsFunction))
  
}

####################################
### Integrating datasets
####################################
writeLines("\n*** Integrating datasets ***\n")

StopWatchStart$IntegrateData  <- Sys.time()

seurat.object.integrated <- IntegrateData(anchorset = seurat.object.anchors.final,
                                          dims = PcaDimsUse,
                                          k.weight = KWeight,
                                          normalization.method = "SCT",
                                          sample.tree = SampleTree,
                                          preserve.order = T,
                                          verbose = T
                                          )

StopWatchEnd$IntegrateData  <- Sys.time()

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SAVE THE R_OBJECT AND LOG FILES
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** SAVE THE R_OBJECTS AND LOG FILES ****\n")

####################################
### Saving the Integration R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving the Integration R object ***\n")
  
  StopWatchStart$SaveRDSIntegration  <- Sys.time()
  
  if (RunsCwl == 1) {
    OutfileRDS<-paste0("R_OBJECTS_CWL/", PrefixOutfiles, ".", ProgramOutdir, "_Integration", ".rds")
  }else{
    OutfileRDS<-paste0(Tempdir, "/R_OBJECTS/", PrefixOutfiles, ".", ProgramOutdir, "_Integration", ".rds")
  }
  saveRDS(seurat.object.integrated, file = OutfileRDS)
  
  StopWatchEnd$SaveRDSIntegration  <- Sys.time()
}else{
  
  writeLines("\n*** Not saving the Integration R objects ***\n")
  
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_Integration_CPUtimes", ".txt")
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
if (RunsCwl == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(Tempdir)
}else{
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
