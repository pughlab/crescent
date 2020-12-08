####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/v3.2/sctransform_vignette.html (SCtransform normalization)
### https://satijalab.org/seurat/v3.2/integration.html (general integration)
### https://satijalab.org/seurat/v3.2/immune_alignment.html (control vs. treatment)
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads integrated datasets R object, presumably from script `Runs_Seurat_v3_MultiDatasets_PCA_Clustering_DimReduction.R`
### 2) Computes differential gene expression (DGE) for options specified by --diff_gene_expr_comparisons
###    Note, all clustering and metadata groups must be provided by --infile_r_object and --infile_metadata
### 3) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_Seurat_v3_MultiDatasets_DGE.R -h'
### for help
####################################

####################################
### THINGS TO DO:
### 1) See Run_Seurat_v3_SingleDataset.R for a list of functions to be implemented
### 2) In:
###    `for (dataset in rownames(InputsTable)) {
###       seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
###     }`
###   Need to avoid that the gsub replaces partial strings, for example 'SMTR03' will be replaced by list_DatasetToType[[SMTR03]] in both
###   'SMTR03' iself and in 'SMTR03t1_NonRad'. The second case is undesired.
###   Generating files like:
###   SMTR_res1.SEURAT_GlobalClustering_Rad_unknownt1_NonRad_TSNEPlot_ColourByCellClusters.pdf
###   When it should be:
###   SMTR_res1.SEURAT_GlobalClustering_Rad_unknown_TSNEPlot_ColourByCellClusters.pdf
####################################

####################################
### COMMENTS ON FINDING DIFFERENTIALLY EXPRESSED GENES
### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object.each_dataset, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object.each_dataset, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
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

ThisScriptName <- "Runs_Seurat_v3_MultiDatasets_DGE.R"
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
              help="Path/name to the Integration R/Seurat object with PCA, clustering and dimension reduction data

                Default = 'No default. It's mandatory to specify this parameter'"),
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
                d1_AAACCTGAGCGGCTTC-1  1            yes
                d2_AAACCTGAGTCGAGTG-1  1            no
                d3_AAACCTGCAAAGGAAG-1  2            yes
                Note: the barcode id must include the dataset ID

                Default = 'NA' (i.e. no metadata to be used for plots or DGE detection)"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each differentially expressed genes test returns only hits that have an UNcorrected p-value < return_thresh

                Default = '0.01'"),
  #
  make_option(c("-f", "--diff_gene_expr_comparisons"), default="1",
              help="One or more of the following options. If more than one, pass the choice <comma> delimited, e.g. '1,2,3'
                '0' = no differentially expressed genes are computed
                '1' = using global cell clusers, compares each cell cluster vs. the rest of cells. Needs `-v 1`
                '2' = using global cell clusers, for each dataset, compares each cell cluster vs. the rest of cells. Needs `-v 1`
                '3' = using global cell clusers, for each dataset, compares each cell cluster vs. the same cluster from other datasets. Needs `-v 1`
                '4' = using global cell clusers, for each dataset type, compares each cell cluster vs. the rest of cells. Needs `-v 1`
                '5' = using global cell clusers, for each dataset type, compares each cell cluster vs. the same cluster from other dataset types. Needs `-v 1`
                '6' = using re-clustered cells, for each dataset, compares each cell cluster vs. the rest of cells. Needs `-v 2`
                '7' = using re-clustered cells, for each dataset type, compares each cell cluster vs. the rest of cells. Needs `-v 3`
                '8' = using metadata annotations, compares each cell class specified by `-b` and `-c` vs. the rest of cells
                '9' = using metadata annotations, for each dataset, compares each cell class specified by `-b` and `-c` vs. the rest of cells
                '10' = using metadata annotations, for each dataset, compares each cell class specified by `-b` and `-c` vs. the same class from other datasets
                '11' = using metadata annotations, for each dataset type, compares each cell class specified by `-b` and `-c` vs. the rest of cells
                '12' = using metadata annotations, for each dataset type, compares each cell class specified by `-b` and `-c` vs. the same class from other dataset types
                '13' = using metadata annotations, for each cell class specified by `-b` and `-c`, compares each subclass vs. other subclasses

                Default = '1'"),
  #
  make_option(c("-b", "--metadata_column_names_for_dge"), default="NA",
              help="Only needed if using -f 8 to 13. It indicates <comma> delimited column names of --infile_metadata to be used for differential gene expression
                Type 'NA' if not using -f 8 to 13

                Default = 'NA'"),
  #
  make_option(c("-d", "--infile_list_metadata_subclasses"), default="NA",
              help="Only needed if using -f 13. It indicates <comma> delimited subclass pairs to compare, one pair per row.
                Type 'NA' to use all subclass pairs indicated by `-b`, `-c` and `-f 13`

                Default = 'NA'"),
  #
  make_option(c("-t", "--assay_to_use_for_dge"), default="SCT",
              help="Only needed if using -f 1 to 12. Either 'SCT' or 'RNA'

                Default = 'SCT'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system

                Default = 'MAX'"),
  #
  make_option(c("-w", "--run_cwl"), default="N",
              help="Indicates if this script is running inside a virtual machine container, such that outfiles are written directly into the 'HOME'. Type 'y/Y' or 'n/N'.
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
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
InfileMetadata          <- opt$infile_metadata
ThreshReturn            <- as.numeric(opt$return_threshold)
DiffGeneExprComparisons <- opt$diff_gene_expr_comparisons
MetadataColNamesForDge  <- opt$metadata_column_names_for_dge
InfileListSubclasses    <- opt$infile_list_metadata_subclasses
AssayForDge             <- opt$assay_to_use_for_dge
NumbCores               <- opt$number_cores
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
    "CRESCENT_CLOUD/frontend_markers",
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES",
    "LOG_FILES"
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
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES",
    "LOG_FILES"
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

RequestedDiffGeneExprComparisons = unlist(strsplit(DiffGeneExprComparisons, ","))

DefaultParameters <- list(
  ### Parameters for Cluster Biomarkers
  FindAllMarkers.MinPct     =  0.25,
  FindAllMarkers.ThreshUse  =  0.25,
  TopDGEForFrontEnd = 6
)

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
### Load integrated datasets R object
####################################
writeLines("\n*** Load integrated datasets R object ***\n")

StopWatchStart$LoadRDSIntegratedDatasets  <- Sys.time()

seurat.object.integrated <- readRDS(InfileRobject)

StopWatchEnd$LoadRDSIntegratedDatasets  <- Sys.time()

####################################
### Loading metadata from --infile_metadata
####################################

if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
  writeLines("\n*** No extra barcode-attributes will be used for dimension reduction plots ***\n")
}else{
  writeLines("\n*** Loading metadata from --infile_metadata for dimension reduction plots ***\n")
  
  CellPropertiesFromMetadata <- data.frame(read.table(InfileMetadata, header = T, row.names = 1, check.names = F))
  
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same dataset)
  # so that barcodes from --infile_metadata and --input can match each other
  rownames(CellPropertiesFromMetadata)<-gsub(x =rownames(CellPropertiesFromMetadata), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = CellPropertiesFromMetadata)
}

DatasetTypes <- unique(seurat.object.integrated$dataset_type)
NumberOfDatasetsTypes <- length(unique(InputsTable[,"DatasetType"]))

####################################
### Check that requested groups for DGE are available
####################################
writeLines("\n*** Check that requested groups for DGE are available ***\n")

listAttributesToSearchInSeuratObject <- list()
if (1 %in% RequestedDiffGeneExprComparisons == T) {
  listAttributesToSearchInSeuratObject[["seurat_clusters"]] <- 1
}
if ((2 %in% RequestedDiffGeneExprComparisons == T) | (3 %in% RequestedDiffGeneExprComparisons == T)) {
  listAttributesToSearchInSeuratObject[["EachDatasetGlobalCellClusters"]] <- 1
}
if ((4 %in% RequestedDiffGeneExprComparisons == T) | (5 %in% RequestedDiffGeneExprComparisons == T)) {
  listAttributesToSearchInSeuratObject[["EachDatasetTypeGlobalCellClusters"]] <- 1
}
if (6 %in% RequestedDiffGeneExprComparisons == T) {
  listAttributesToSearchInSeuratObject[["EachDatasetCellReClusters"]] <- 1
}
if (7 %in% RequestedDiffGeneExprComparisons == T) {
  listAttributesToSearchInSeuratObject[["EachDatasetTypeCellReClusters"]] <- 1
}

if ((8 %in% RequestedDiffGeneExprComparisons == T) |
    (9 %in% RequestedDiffGeneExprComparisons == T) |
    (10 %in% RequestedDiffGeneExprComparisons == T) |
    (11 %in% RequestedDiffGeneExprComparisons == T) |
    (12 %in% RequestedDiffGeneExprComparisons == T) |
    (13 %in% RequestedDiffGeneExprComparisons == T)
) {
  MetadataColNamesForDge.list  = unlist(strsplit(MetadataColNamesForDge, ","))
  for (property in MetadataColNamesForDge.list) {
    listAttributesToSearchInSeuratObject[[property]] <- 1
  }
}

if (all(names(listAttributesToSearchInSeuratObject) %in% names(seurat.object.integrated@meta.data)) == T) {
    writeLines("\nOk\n")
}else{
  stop(paste0(paste0("\n\nERROR!!! couldn't find all '", length(names(listAttributesToSearchInSeuratObject)), "' requested properties in Seurat object metadata:", "\n\n"),
              paste0("Requested: '"    , paste(names(listAttributesToSearchInSeuratObject), collapse = "' '"), "'\n\n"),
              paste0("Seurat_object: '", paste(names(seurat.object.integrated@meta.data),   collapse = "' '"), "'\n\n")))
}

if (regexpr("^NA$", InfileListSubclasses, ignore.case = T)[1] == 1) {
  writeLines("\nOk\n")
}else{
  if (13 %in% RequestedDiffGeneExprComparisons == T) {
    writeLines("\nOk\n")
  }else{
    stop(paste0("\n\nERROR!!! option '-d INFILE' requires '-f 13'"))
  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO COMPUTE DGE USING GLOBAL CLUSTERS
################################################################################################################################################
################################################################################################################################################

####################################
### Finding differentially expressed genes (1): using global cell clusers, compares each cell cluster vs. the rest of cells
####################################

ClusterIdent <-seurat.object.integrated@meta.data$seurat_clusters
NumberOfClusters<-length(unique(ClusterIdent))

if (1 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (1): using global cell clusers, compares each cell cluster vs. the rest of cells ***\n")

  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$seurat_clusters
  
  StopWatchStart$FindDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()

  print(paste0("Number of groups = ", length(unique(seurat.object.integrated@meta.data$seurat_clusters))))

  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  seurat.object.integrated.markers <- FindAllMarkers(object = seurat.object.integrated, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
  SimplifiedDiffExprGenes.df <- seurat.object.integrated.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
  OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "GlobalClustering", ".tsv.bz2")
  Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
  write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
  close(Outfile.con)

  StopWatchEnd$FindDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()

  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    ####################################
    ### Outfiles for web app: top differentially expressed genes
    ####################################
    writeLines("\n*** Outfiles for web app: top differentially expressed genes ***\n")

    StopWatchStart$OutTopDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()

    top_genes_by_cluster<-(seurat.object.integrated.markers %>% group_by(cluster) %>% top_n(DefaultParameters$TopDGEForFrontEnd, avg_logFC))
    globalMarkersFile <- top_genes_by_cluster[,c("gene","cluster","p_val","avg_logFC")]
    write.table(globalMarkersFile, paste0(Tempdir, "/CRESCENT_CLOUD/frontend_markers/","TopMarkersPerCluster.tsv"), row.names = F, sep="\t", quote = F)

    StopWatchEnd$OutTopDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()

  }
}

####################################
### Finding differentially expressed genes (2): using global cell clusers, for each dataset, compares each cell cluster vs. the rest of cells
####################################

if (2 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (2): using global cell clusers, for each dataset, compares each cell cluster vs. the rest of cells ***\n")

  ####################################
  ### Loops each dataset
  ####################################
  NumberOfDatasets <- 0
  for (dataset in rownames(InputsTable)) {
    
    StopWatchStart$FindDiffMarkersEachDatasetGlobalClustersVsRestOfCells[[dataset]]  <- Sys.time()

    NumberOfDatasets <- NumberOfDatasets + 1
    print(NumberOfDatasets)
    
    if (exists(x = "seurat.object.each_dataset") == T) {
      rm(seurat.object.each_dataset)
    }
    Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset
    seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
    Idents(object = seurat.object.each_dataset) <- seurat.object.each_dataset@meta.data$EachDatasetGlobalCellClusters
    
    print(paste0("Number of groups = ", length(unique(seurat.object.each_dataset@meta.data$EachDatasetGlobalCellClusters))))

    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset@meta.data))
    seurat.object.each_dataset.markers <- FindAllMarkers(object = seurat.object.each_dataset, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_dataset.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
    OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "GlobalClustering", "_dataset_", dataset, ".tsv.bz2")
    Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
    close(Outfile.con)
    
    StopWatchEnd$FindDiffMarkersEachDatasetGlobalClustersVsRestOfCells[[dataset]]  <- Sys.time()

    if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
      ####################################
      ### Outfiles for web app: top differentially expressed genes oer dataset
      ####################################
      writeLines("\n*** Outfiles for web app: top differentially expressed genes per dataset ***\n")

      StopWatchStart$OutDiffMarkersEachDatasetGlobalClustersVsRestOfCellsFrontEnd[[dataset]]  <- Sys.time()

      # top markers per dataset
      top_genes_by_cluster<-(seurat.object.each_dataset.markers %>% group_by(cluster) %>% top_n(DefaultParameters$TopDGEForFrontEnd, avg_logFC))
      datasetMarkersFile <- top_genes_by_cluster[,c("gene","cluster","p_val","avg_logFC")]
      write.table(datasetMarkersFile, paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_markers/",dataset,"_TopMarkersPerCluster.tsv"), row.names = F, sep="\t", quote = F)

      # groups.tsv per dataset
      OutfileClustersPerDatasetDataframe  <- data.frame(NAME = rownames(seurat.object.each_dataset@meta.data), Seurat_Dataset_Clusters = seurat.object.each_dataset@meta.data$seurat_clusters)
      groupsPerDatasetDataframeString <- sapply(OutfileClustersPerDatasetDataframe, as.character)
      groupsPerDatasetDataframeStringTYPE <- rbind(data.frame(NAME = "TYPE", Seurat_Dataset_Clusters = "group"), groupsPerDatasetDataframeString)
      colnames(groupsPerDatasetDataframeStringTYPE) <- c("NAME", paste("Seurat_",dataset,"_Clusters_Resolution_", Resolution, sep = "", collapse = ""))

      OutfileClustersPerDatasetFile<-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_groups/",dataset,"_groups.tsv")
      write.table(data.frame(groupsPerDatasetDataframeStringTYPE), file = OutfileClustersPerDatasetFile, row.names = F, col.names = T, sep="\t", quote = F, append = T)

      StopWatchEnd$OutDiffMarkersEachDatasetGlobalClustersVsRestOfCellsFrontEnd[[dataset]]  <- Sys.time()

    }
  }
}

####################################
### Finding differentially expressed genes (3): using global cell clusers, for each dataset, compares each cell cluster vs. the same cluster from other datasets
####################################

if (3 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (3): using global cell clusers, for each dataset, compares each cell cluster vs. the same cluster from other datasets ***\n")

  StopWatchStart$FindDiffMarkersEachDatasetGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()

  print(paste0("Number of groups = ", length(unique(seurat.object.integrated@meta.data$seurat_clusters))))
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "GlobalClustering", "_", "DatasetEquivalentClusters", ".tsv.bz2")
  Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
  HeadersOrder <- paste("cluster1", "cluster2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
  write.table(HeadersOrder, file = Outfile.con, row.names = F, col.names = F, sep="", quote = F)

  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$EachDatasetGlobalCellClusters

  for (cluster in unique(seurat.object.integrated@meta.data$seurat_clusters)) {
    for (dataset1 in rownames(InputsTable)) {
      for (dataset2 in rownames(InputsTable)) {
        Cluster1 <- paste(dataset1, cluster, sep = "_c")
        Cluster2 <- paste(dataset2, cluster, sep = "_c")
        if (Cluster1 == Cluster2) {
          ### Skip
        }else if ((sum(seurat.object.integrated@meta.data$EachDatasetGlobalCellClusters == Cluster1) >= 3) & ((sum(seurat.object.integrated@meta.data$EachDatasetGlobalCellClusters == Cluster2) >= 3))) {
          print (paste0(Cluster1, " vs. ", Cluster2))
          seurat.object.integrated.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.integrated, assay = AssayForDge, only.pos = F, ident.1 = Cluster1, ident.2 = Cluster2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
          seurat.object.integrated.each_equivalent_cluster.markers$cluster1 <- Cluster1
          seurat.object.integrated.each_equivalent_cluster.markers$cluster2 <- Cluster2
          seurat.object.integrated.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.integrated.each_equivalent_cluster.markers)
          write.table(seurat.object.integrated.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
        }else{
          print(paste0("Skip cluster ", Cluster1, " vs. ", Cluster2, " because there were not >= 3 cells in at least one of them"))
        }
      }
    }
  }
  
  close(Outfile.con)
  StopWatchEnd$FindDiffMarkersEachDatasetGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()
}

####################################
### Finding differentially expressed genes (4): using global cell clusers, for each dataset type, compares each cell cluster vs. the rest of cells
####################################
if (4 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes (4): using global cell clusers, for each dataset type, compares each cell cluster vs. the rest of cells ***\n")
  
  ####################################
  ### Loops each dataset type
  ####################################
  
  NumberOfDatasetTypes <- 0
  for (dataset_type in DatasetTypes) {
    
    StopWatchStart$FindDiffMarkersEachDatasetTypeGlobalClustersVsRestOfCells[[dataset_type]]  <- Sys.time()

    NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
    print(NumberOfDatasetTypes)
    
    if (exists(x = "seurat.object.each_dataset_type") == T) {
      rm(seurat.object.each_dataset_type)
    }
    Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset_type
    seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
    Idents(object = seurat.object.each_dataset_type) <- seurat.object.each_dataset_type@meta.data$EachDatasetTypeGlobalCellClusters

    print(paste0("Number of groups = ", length(unique(seurat.object.each_dataset_type@meta.data$EachDatasetTypeGlobalCellClusters))))

    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset_type@meta.data))
    seurat.object.each_dataset_type.markers <- FindAllMarkers(object = seurat.object.each_dataset_type, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_dataset_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
    OutfileDiffGeneExpression <- paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "GlobalClustering", "_dataset_type_", dataset_type, ".tsv.bz2")
    Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
    close(Outfile.con)
    
    StopWatchEnd$FindDiffMarkersEachDatasetTypeGlobalClustersVsRestOfCells[[dataset_type]]  <- Sys.time()
  }
}

####################################
### Finding differentially expressed genes (5): using global cell clusers, for each dataset type, compares each cell cluster vs. the same cluster from other dataset types
####################################

if (5 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes (5): using global cell clusers, compares each cell cluster from each dataset type vs. the same cluster from other dataset types ***\n")

  StopWatchStart$FindDiffMarkersEachDatasetTypeGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()
  
  print(paste0("Number of groups = ", length(unique(seurat.object.integrated@meta.data$seurat_clusters))))
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  
  OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "GlobalClustering", "_", "DatasetTypeEquivalentClusters", ".tsv.bz2")
  Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
  HeadersOrder <- paste("cluster1", "cluster2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
  write.table(HeadersOrder, file = Outfile.con, row.names = F, col.names = F, sep="", quote = F)
  
  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$EachDatasetTypeGlobalCellClusters
  
  for (cluster in unique(seurat.object.integrated@meta.data$seurat_clusters)) {
    for (dataset_type1 in unique(InputsTable[,"DatasetType"])) {
      for (dataset_type2 in unique(InputsTable[,"DatasetType"])) {
        Cluster1 <- paste(dataset_type1, cluster, sep = "_c")
        Cluster2 <- paste(dataset_type2, cluster, sep = "_c")
        if (Cluster1 == Cluster2) {
          ### Skip
        }else if ((sum(seurat.object.integrated@meta.data$EachDatasetTypeGlobalCellClusters == Cluster1) >= 3) & ((sum(seurat.object.integrated@meta.data$EachDatasetTypeGlobalCellClusters == Cluster2) >= 3))) {
          print (paste0("Diff gene expression: ", Cluster1, " vs. ", Cluster2))
          seurat.object.integrated.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.integrated, assay = AssayForDge, only.pos = F, ident.1 = Cluster1, ident.2 = Cluster2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
          seurat.object.integrated.each_equivalent_cluster.markers$cluster1 <- Cluster1
          seurat.object.integrated.each_equivalent_cluster.markers$cluster2 <- Cluster2
          seurat.object.integrated.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.integrated.each_equivalent_cluster.markers)
          write.table(seurat.object.integrated.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
          rm(seurat.object.integrated.each_equivalent_cluster.markers)
          
        }else{
          print(paste0("Skip cluster ", Cluster1, " vs. ", Cluster2, " because there were not >= 3 cells in at least one of them"))
        }
      }
    }
  }
  
  close(Outfile.con)
  StopWatchEnd$FindDiffMarkersEachDatasetTypeGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO COMPUTE DGE USING RE-CLUSTERS
################################################################################################################################################
################################################################################################################################################

####################################
### Finding differentially expressed genes (6): using re-clustered cells, for each dataset, compares each cell cluster vs. the rest of cells
####################################
if (6 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (6): using re-clustered cells, for each dataset, compares each cell cluster vs. the rest of cells ***\n")
  
  ####################################
  ### Loops each dataset
  ####################################
  NumberOfDatasets <- 0
  for (dataset in rownames(InputsTable)) {
    
    StopWatchStart$FindDiffMarkersReclusteredVsRestOfCells[[dataset]]  <- Sys.time()
    
    NumberOfDatasets <- NumberOfDatasets + 1
    print(NumberOfDatasets)

    if (exists(x = "seurat.object.each_dataset") == T) {
      rm(seurat.object.each_dataset)
    }
    Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset
    seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
    Idents(object = seurat.object.each_dataset) <- seurat.object.each_dataset@meta.data$EachDatasetCellReClusters

    print(paste0("Number of clusters = ", length(unique(seurat.object.each_dataset@meta.data$EachDatasetCellReClusters))))

    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset@meta.data))
    seurat.object.each_dataset.markers <- FindAllMarkers(object = seurat.object.each_dataset, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_dataset.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
    OutfileDiffGeneExpression <- paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "EachDatasetReclustered_", dataset, ".tsv.bz2")
    Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
    close(Outfile.con)
    
    StopWatchEnd$FindDiffMarkersReclusteredVsRestOfCells[[dataset]]  <- Sys.time()
  }
}

####################################
### Finding differentially expressed genes (7): using re-clustered cells, for each dataset type, compares each cell cluster vs. the rest of cells
####################################

if (7 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes (7): using re-clustered cells, for each dataset type, compares each cell cluster vs. the rest of cells ***\n")
  
  ####################################
  ### Loops each dataset type
  ####################################
  for (dataset_type in DatasetTypes) {
    
    StopWatchStart$FindDiffMarkersEachDatasetTypeReclusteredVsRestOfCells[[dataset_type]]  <- Sys.time()
    
    if (exists(x = "seurat.object.each_dataset_type") == T) {
      rm(seurat.object.each_dataset_type)
    }
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset_type
    seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
    Idents(object = seurat.object.each_dataset_type) <- seurat.object.each_dataset_type@meta.data$EachDatasetTypeCellReClusters

    print(paste0("Number of groups = ", length(unique(seurat.object.each_dataset_type@meta.data$EachDatasetTypeCellReClusters))))

    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset_type@meta.data))
    seurat.object.each_dataset_type.markers <- FindAllMarkers(object = seurat.object.each_dataset_type, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_dataset_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
    OutfileDiffGeneExpression <- paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "EachDatasetTypeReclustered_", dataset_type, ".tsv.bz2")
    Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
    close(Outfile.con)

    StopWatchEnd$FindDiffMarkersEachDatasetTypeReclusteredVsRestOfCells[[dataset_type]]  <- Sys.time()

  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO OBTAIN DIFFERENTIAL GENE EXPRESSION USING METADATA
################################################################################################################################################
################################################################################################################################################

####################################
### Finding differentially expressed genes (8): using metadata, compares each cell class vs. the rest of cells
####################################
if (8 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (8): using metadata, compares each cell class vs. the rest of cells ***\n")

  for (property in MetadataColNamesForDge.list) {

    NumberOfClassesInThisProperty <- length(unique(seurat.object.integrated@meta.data[[property]]))
    print(paste0("Number of groups in '", property, "' = ", NumberOfClassesInThisProperty))
    if (NumberOfClassesInThisProperty > 1) {

      StopWatchStart$FindDiffMarkersEachMetadataClassVsRestOfCells[[property]]  <- Sys.time()
      
      if (exists(x = "seurat.object.each_property") == T) {
        rm(seurat.object.each_property)
      }
      Idents(object = seurat.object.integrated) <- seurat.object.integrated[[property]]
      seurat.object.each_property <- seurat.object.integrated

      ####################################
      ### Finding differentially expressed genes
      ####################################
      FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_property@meta.data))
      seurat.object.each_property.markers <- FindAllMarkers(object = seurat.object.each_property, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
      SimplifiedDiffExprGenes.df <- seurat.object.each_property.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
      colnames(SimplifiedDiffExprGenes.df) <- c("class","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")
      OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "Metadata", "_", property, "_", "AllDatasets", ".tsv.bz2")
      Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
      write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
      close(Outfile.con)
      
      StopWatchEnd$FindDiffMarkersEachMetadataClassVsRestOfCells[[property]]  <- Sys.time()
    }
  }
}

####################################
### Finding differentially expressed genes (9): using metadata annotations, for each dataset, compares each cell class vs. the rest of cells
####################################
if (9 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (9): using metadata annotations, for each dataset, compares each cell class vs. the rest of cells ***\n")

  ####################################
  ### Loops each dataset
  ####################################
  for (dataset in rownames(InputsTable)) {

    if (exists(x = "seurat.object.each_dataset") == T) {
      rm(seurat.object.each_dataset)
    }
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset
    seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)

    for (property in MetadataColNamesForDge.list) {
      
      StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetVsRestOfCells[[dataset]][[property]] <- Sys.time()

      NumberOfClassesInThisProperty <- length(unique(seurat.object.each_dataset@meta.data[[property]]))
      print(paste0("Number of groups in '", property, " in dataset ", dataset, "' = ", NumberOfClassesInThisProperty))
      if (NumberOfClassesInThisProperty > 1) {
        
        Idents(object = seurat.object.each_dataset) <- seurat.object.each_dataset@meta.data[[property]]

        FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset@meta.data))
        seurat.object.each_dataset.markers <- FindAllMarkers(object = seurat.object.each_dataset, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
        SimplifiedDiffExprGenes.df <- seurat.object.each_dataset.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
        colnames(SimplifiedDiffExprGenes.df) <- c("class","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")
        OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "Metadata", "_", property, "_",  dataset, ".tsv.bz2")
        Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
        write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
        close(Outfile.con)
        
      }
      
      StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetVsRestOfCells[[dataset]][[property]] <- Sys.time()
    }
  }
}

####################################
### Finding differentially expressed genes (10): using metadata annotations, for each dataset, compares each cell class specified vs. the same class from other datasets
####################################
if (10 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (10): using metadata annotations, for each dataset, compares each cell class specified vs. the same class from other datasets ***\n")

  for (property in MetadataColNamesForDge.list) {

    NumberOfClassesInThisProperty <- length(unique(seurat.object.integrated@meta.data[[property]]))
    print(paste0("Number of groups in '", property, "' = ", NumberOfClassesInThisProperty))
    if (NumberOfClassesInThisProperty > 1) {

      StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetVsSameClassInOtherDatasets[[property]]  <- Sys.time()

      ####################################
      ### Subsets seurat object per property
      ####################################
      if (exists(x = "seurat.object.each_property") == T) {
        rm(seurat.object.each_property)
      }
      seurat.object.each_property <- seurat.object.integrated
      DatasetAndProperty <- unlist(x = strsplit(x = paste(seurat.object.each_property@meta.data$dataset, seurat.object.each_property@meta.data[[property]], sep = "_", collapse = "\n"), split = "\n"))
      seurat.object.each_property <- AddMetaData(object = seurat.object.each_property, metadata = DatasetAndProperty, col.name = "DatasetAndProperty")

      FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_property@meta.data))
      OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "Metadata", "_", property, "_",  "DatasetEquivalentClasses", ".tsv.bz2")
      Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
      HeadersOrder <- paste("class1", "class2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
      write.table(HeadersOrder, file = Outfile.con, row.names = F, col.names = F, sep="", quote = F)
      
      Idents(object = seurat.object.each_property) <- DatasetAndProperty
      
      for (class in unique(seurat.object.each_property@meta.data[[property]])) {
        for (dataset_1 in rownames(InputsTable)) {
          for (dataset_2 in rownames(InputsTable)) {
            Class_Dataset1 <- paste0(dataset_1, "_", class)
            Class_Dataset2 <- paste0(dataset_2, "_", class)
            N_Class_Dataset1 <- sum(seurat.object.each_property$DatasetAndProperty == Class_Dataset1)
            N_Class_Dataset2 <- sum(seurat.object.each_property$DatasetAndProperty == Class_Dataset2)

            if (Class_Dataset1 == Class_Dataset2) {
              ### Skip
            }else if (N_Class_Dataset1 >= 3 & N_Class_Dataset2 >= 3) {
              print (paste0(Class_Dataset1, " vs. ", Class_Dataset2))
              seurat.object.each_property.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.each_property, assay = AssayForDge, only.pos = F, ident.1 = Class_Dataset1, ident.2 = Class_Dataset2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
              seurat.object.each_property.each_equivalent_cluster.markers$class1 <- Class_Dataset1
              seurat.object.each_property.each_equivalent_cluster.markers$class2 <- Class_Dataset2
              seurat.object.each_property.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.each_property.each_equivalent_cluster.markers)
              write.table(seurat.object.each_property.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
              rm(seurat.object.each_property.each_equivalent_cluster.markers)
            }else{
              print(paste0("Skip class ", Class_Dataset1, " vs. ", Class_Dataset2, " because there were not >= 3 cells in at least one of them"))
            }
          }
        }
      }
      
      close(Outfile.con)
      StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetVsSameClassInOtherDatasets[[property]]  <- Sys.time()

     }
  }
}

####################################
### Finding differentially expressed genes (11): using metadata annotations, for each dataset_type, compares each cell class vs. the rest of cells
####################################
if (11 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (11): using metadata annotations, for each dataset_type, compares each cell class vs. the rest of cells ***\n")

  ####################################
  ### Loops each dataset_type
  ####################################
  for (dataset_type in DatasetTypes) {

    if (exists(x = "seurat.object.each_dataset_type") == T) {
      rm(seurat.object.each_dataset_type)
    }
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset_type
    seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)

    for (property in MetadataColNamesForDge.list) {
      
      StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetTypeVsRestOfCells[[dataset_type]][[property]] <- Sys.time()
      
      NumberOfClassesInThisProperty <- length(unique(seurat.object.each_dataset_type@meta.data[[property]]))
      print(paste0("Number of classes in '", property, " in dataset_type ", dataset_type, " = ", NumberOfClassesInThisProperty))
      if (NumberOfClassesInThisProperty > 1) {

        Idents(object = seurat.object.each_dataset_type) <- seurat.object.each_dataset_type@meta.data[[property]]

        FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset_type@meta.data))
        seurat.object.each_dataset_type.markers <- FindAllMarkers(object = seurat.object.each_dataset_type, assay = AssayForDge, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
        SimplifiedDiffExprGenes.df <- seurat.object.each_dataset_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
        colnames(SimplifiedDiffExprGenes.df) <- c("class","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")
        OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "Metadata", "_", property, "_", dataset_type, ".tsv.bz2")
        Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
        write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
        close(Outfile.con)
      }
      
      StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetTypeVsRestOfCells[[dataset_type]][[property]] <- Sys.time()
    }
  }
}

####################################
### Finding differentially expressed genes (12): using metadata annotations, for each dataset type, compares each cell class vs. the same class from other dataset types
####################################
if (12 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (12): using metadata annotations, for each dataset type, compares each cell class vs. the same class from other dataset types ***\n")

  if (NumberOfDatasetsTypes > 1) {

    for (property in MetadataColNamesForDge.list) {

      NumberOfClassesInThisProperty <- length(unique(seurat.object.integrated@meta.data[[property]]))
      print(paste0("Number of groups in '", property, "' = ", NumberOfClassesInThisProperty))
      if (NumberOfClassesInThisProperty > 1) {
        
        StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetTypeVsSameClassInOtherDatasetTypes  <- Sys.time()

        ####################################
        ### Subsets seurat object per property
        ####################################
        if (exists(x = "seurat.object.each_property") == T) {
          rm(seurat.object.each_property)
        }
        seurat.object.each_property <- seurat.object.integrated
        DatasetTypeAndProperty <- unlist(x = strsplit(x = paste(seurat.object.each_property@meta.data$dataset_type, seurat.object.each_property@meta.data[[property]], sep = "_", collapse = "\n"), split = "\n"))
        seurat.object.each_property <- AddMetaData(object = seurat.object.each_property, metadata = DatasetTypeAndProperty, col.name = "DatasetTypeAndProperty")

        FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_property@meta.data))
        OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "Metadata", "_", property, "_",  "DatasetTypeEquivalentClasses", ".tsv.bz2")
        Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
        HeadersOrder <- paste("class1", "class2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
        write.table(HeadersOrder, file = Outfile.con, row.names = F, col.names = F, sep="", quote = F)

        Idents(object = seurat.object.each_property) <- DatasetTypeAndProperty
        
        for (class in unique(seurat.object.each_property@meta.data[[property]])) {
          for (dataset_type_1 in unique(InputsTable[,"DatasetType"])) {
            for (dataset_type_2 in unique(InputsTable[,"DatasetType"])) {
              Class_DatasetType1 <- paste0(dataset_type_1, "_", class)
              Class_DatasetType2 <- paste0(dataset_type_2, "_", class)
              N_Class_DatasetType1 <- sum(seurat.object.each_property$DatasetTypeAndProperty == Class_DatasetType1)
              N_Class_DatasetType2 <- sum(seurat.object.each_property$DatasetTypeAndProperty == Class_DatasetType2)

              if (Class_DatasetType1 == Class_DatasetType2) {
                ### Skip
              }else if (N_Class_DatasetType1 >= 3 & N_Class_DatasetType2 >= 3) {
                print (paste0(Class_DatasetType1, " vs. ", Class_DatasetType2))
                seurat.object.each_property.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.each_property, assay = AssayForDge, only.pos = F, ident.1 = Class_DatasetType1, ident.2 = Class_DatasetType2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
                seurat.object.each_property.each_equivalent_cluster.markers$class1 <- Class_DatasetType1
                seurat.object.each_property.each_equivalent_cluster.markers$class2 <- Class_DatasetType2
                seurat.object.each_property.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.each_property.each_equivalent_cluster.markers)
                write.table(seurat.object.each_property.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
                rm(seurat.object.each_property.each_equivalent_cluster.markers)
              }else{
                print(paste0("Skip class ", Class_DatasetType1, " vs. ", Class_DatasetType2, " because there were not >= 3 cells in at least one of them"))
              }
            }
          }
        }
        
        close(Outfile.con)
        StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetTypeVsSameClassInOtherDatasetTypes  <- Sys.time()
      }
    }
  }else{
    print(paste0("Skip because could find only '", NumberOfDatasetsTypes, "' dataset types"))
  }
}

####################################
### Finding differentially expressed genes (13): using metadata annotations, for each cell class specified by `-b` and `-c`, compares each subclass vs. other subclasses
####################################
if (13 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes (13): using metadata annotations, for each cell class specified by `-b` and `-c`, compares each subclass vs. other subclasses ***\n")

  if (regexpr("^NA$", InfileListSubclasses, ignore.case = T)[1] == 1) {
    writeLines("\n*** All subclass pairs will be compared ***\n")
  }else{
    writeLines("\n*** Load subclass pairs to compare ***\n")
    SubclassesToUse.df<-read.table(InfileListSubclasses, header = F, row.names = NULL, stringsAsFactors = F)
    colnames(SubclassesToUse.df)<-c("SubclassPair")
    print(paste0("Will restrict comparisons to ", length(SubclassesToUse.df[,"SubclassPair"]), " subclass pairs"))
  }

  for (property in MetadataColNamesForDge.list) {
    NumberOfClassesInThisProperty <- length(unique(seurat.object.integrated@meta.data[[property]]))
    print(paste0("Number of groups in '", property, "' = ", NumberOfClassesInThisProperty))
    if (NumberOfClassesInThisProperty > 1) {
      
      StopWatchStart$FindDiffMarkersEachMetadataEachSubclassVsOtherSubclasses[[property]]  <- Sys.time()
      
      ####################################
      ### Subsets seurat object per property
      ####################################
      if (exists(x = "seurat.object.each_property") == T) {
        rm(seurat.object.each_property)
      }
      seurat.object.each_property <- seurat.object.integrated
      
      Idents(object = seurat.object.each_property) <- property

      FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_property@meta.data))
      OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_DiffExprMarkers_", "Metadata", "_", property, "_",  "SubClassesAgainstEachOther", ".tsv.bz2")
      Outfile.con <- bzfile(OutfileDiffGeneExpression, "w")
      HeadersOrder <- paste("class1", "class2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
      write.table(HeadersOrder, file = Outfile.con, row.names = F, col.names = F, sep="", quote = F)

      for (subclass1 in unique(seurat.object.each_property@meta.data[[property]])) {
        for (subclass2 in unique(seurat.object.each_property@meta.data[[property]])) {
          N_Class_1 <- sum(seurat.object.each_property[[property]] == subclass1)
          N_Class_2 <- sum(seurat.object.each_property[[property]] == subclass2)
          
          if (subclass1 == subclass2) {
            ### Skip
          }else if (N_Class_1 >= 3 & N_Class_2 >= 3) {
            SubclassPairToCompare <- paste0(subclass1, ",", subclass2)
            ToRunThisPair <- 0
            if (regexpr("^NA$", InfileListSubclasses, ignore.case = T)[1] == 1) {
              ToRunThisPair <- 1
            }else if (SubclassPairToCompare %in% SubclassesToUse.df[,"SubclassPair"] == T) {
              ToRunThisPair <- 1
            }
            if (ToRunThisPair == 1) {
              print (paste0(subclass1, " vs. ", subclass2))
              seurat.object.each_property.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.each_property, assay = AssayForDge, only.pos = F, ident.1 = subclass1, ident.2 = subclass2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
              seurat.object.each_property.each_equivalent_cluster.markers$class1 <- subclass1
              seurat.object.each_property.each_equivalent_cluster.markers$class2 <- subclass2
              seurat.object.each_property.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.each_property.each_equivalent_cluster.markers)
              write.table(seurat.object.each_property.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = Outfile.con, row.names = F, col.names = F, sep="\t", quote = F, append = T)
              rm(seurat.object.each_property.each_equivalent_cluster.markers)
            }
          }else{
            print(paste0("Skip class ", subclass1, " vs. ", subclasss2, " because there were not >= 3 cells in at least one of them"))
          }
        }
      }
      close(Outfile.con)
      StopWatchEnd$FindDiffMarkersEachMetadataEachSubclassVsOtherSubclasses[[property]]  <- Sys.time()
    }
  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SAVE THE WHOLE RUN R_OBJECT AND LOG FILES
################################################################################################################################################
################################################################################################################################################

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
