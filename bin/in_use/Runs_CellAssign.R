####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Erik Christensen - echris3@uwo.ca
### Script made based on https://irrationone.github.io/cellassign/index.html
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads scRNA-seq data 
### 2) Computes cell size factors
### 3) Loads marker gene matrix
### 4) Subsets the scRNA-seq data based on the marker gene matrix and transposes it
### 5) Runs cellassign
### 6) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_cellassign.R -h'
### for help
####################################

####################################
### THINGS TO DO:
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(DropletUtils)) # (Bioconductor) to handle MTX/H5 format files. Note it has about the same speed than library(earlycross) which can't handle H5
suppressPackageStartupMessages(library(Seurat))       # (CRAN) to run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(gtools))       # (CRAN) to do alphanumeric sorting. Only needed if using `-w Y`.
suppressPackageStartupMessages(library(scran))        # (Bioconductor) to compute cell size factors
suppressPackageStartupMessages(library(cellassign))   # (GitHub Irrationone/cellassign) to assign scRNA-seq data to known cell types
####################################

####################################
### Required external packages
####################################
### Tensorflow 2.1.0 with tensorflow-probability (required for cellassign)
###     install.packages("tensorflow")
###     tensorflow::install_tensorflow(extra_packages='tensorflow-probability', version = "2.1.0")
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_cellassign.R"
ProgramOutdir  <- "CELLASSIGN"

####################################
### Get inputs from command line argumets
####################################
#
option_list <- list(
  make_option(c("-i", "--input"), default="NA",
              help="Path/name to either a read counts matrix in either 'MTX' or 'TSV' format (see parameter -t)
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--input_type"), default="NA", 
              help="Either 'MTX', 'TSV' or 'HDF5'
                'MTX'  is the path/name to an MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
                'TSV'  is the path/name of a <tab> delimited *file* with genes in rows vs. barcodes in columns
                'HDF5' is the path/name of a *file* in hdf5 format (e.g. from Cell Ranger)
                Note 'MTX' files can be the outputs from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-m", "--input_markers"), default="NA",
              help="Path/name to a binary marker by cell type matrix in 'CSV' format.
              Default = 'No default. It's mandatory to specify this parameter'"), 
  #
  make_option(c("-x", "--design_matrix"), default="NA",
              help="Path to an N by P design matrix of covariates for any patient/batch specific effects. 
              Default = No default, optional."),
  #
  make_option(c("-l", "--learning_rate"), default="0.01",
              help="The learning rate for ADAM optimization
              Default = '0.01'"), 
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),
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
InputMarkers            <- opt$input_markers
LearningRate            <- as.numeric(opt$learning_rate)
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl
MaxGlobalVariables      <- opt$max_global_variables

####################################
### Define default parameters
####################################

DefaultParameters <- list(
  
)


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
ListMandatory<-list("input", "input_type", "input_matrix","outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
  }
}

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
    "CRESCENT_CLOUD/frontend_metadata",
    "CRESCENT_CLOUD/frontend_markers",
    "CRESCENT_CLOUD/frontend_qc",
    "CRESCENT_CLOUD/frontend_groups",
    "LOG_FILES",
    "ANNOTATIONS"
  )
  
}else{
  #PrefixOutfiles <- c(paste0(PrefixOutfiles,"_res",Resolution))
  ## Using `Tempdir/DIRECTORY` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## 'DIRECTORY' is one of the directories specified at FILE_TYPE_OUT_DIRECTORIES
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  Tempdir        <- "~/temp" 
  #
  CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
  UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
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
    "ANNOTATIONS"
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

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_UsedOptions.txt")

TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=paste0("Run started: ", TimeOfRun, "\n"))

write(file = OutfileOptionsUsed, x=c("Commands used:"), append = T)
for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

cat(file = OutfileOptionsUsed, x=paste0("\n", "One-line-commands used:", "\n", "`Rscript /path_to/", ThisScriptName), append = T)
for (optionInput in option_list) {
  cat(file = OutfileOptionsUsed, x=(paste0(" ", optionInput@short_flag, " ", opt[optionInput@dest])),append = T)
}
cat(file = OutfileOptionsUsed, x="`\n", append = T)

StopWatchEnd$ReportUsedOptions  <- Sys.time()

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo.txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)



####################################
### Load scRNA-seq data 
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

StopWatchStart$LoadScRNAseqData  <- Sys.time()

if (regexpr("^MTX$", InputType, ignore.case = T)[1] == 1) {
  writeLines(paste0("\n*** Loading MTX infiles ***\n"))
  input.matrix <- Read10X(data.dir = Input) 
}else if (regexpr("^TSV$", InputType, ignore.case = T)[1] == 1) {
  writeLines(paste0("\n*** Loading matrix of genes (rows) vs. barcodes (columns) ***\n"))
  ## Note `check.names = F` is needed for both `fread` and `data.frame`
  input.matrix <- as.matrix(data.frame(fread(Input, check.names = F), row.names=1, check.names = F))
}else if (regexpr("^HDF5$", InputType, ignore.case = T)[1] == 1) {
  writeLines(paste0("\n*** Loading HDF5 infile ***\n"))
  input.matrix <- Read10X_h5(filename = Input, use.names = T, unique.features = T) 
}else{
  stop(paste0("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n"))
}
dim(input.matrix)

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Calculate Size Factors 
####################################

writeLines("\n*** Computing Size Factors ***\n")
StopWatchStart$ComputeSizeFactors <- Sys.time()
size.factors <- calculateSumFactors(input.matrix, min.mean = 1)
StopWatchEnd$ComputeSizeFactors <- Sys.time()

####################################
### load marker gene matrix
####################################

writeLines("\n*** Loading Marker Gene Matrix ***\n")
StopWatchStart$LoadMarkerGeneMatrix <- Sys.time()
marker.matrix <- as.matrix(data.frame(fread(InputMarkers, check.names = F), row.names = 1, check.names = F))
StopWatchEnd$LoadMarkerGeneMatrix <- Sys.time()

####################################
### Subset the expression matrix with marker genes and transpose
####################################

writeLines("\n*** Removing Non Marker Genes ***\n")
StopWatchStart$RemoveNonMarkerGenes <- Sys.time()
input.matrix.final <- as.matrix(t(input.matrix[as.character(row.names(marker.matrix)),]))
StopWatchEnd$RemoveNonMarkerGenes <- Sys.time()

####################################
### Run cellassign 
####################################

writeLines("\n*** Running CellAssign ***\n")
StopWatchStart$RunCellAssign <- Sys.time()
cell.fit <- cellassign(exprs_obj = input.matrix.final, 
                        marker_gene_info = marker.matrix,
                        s = size.factors,
                        learning_rate = LearningRate,
                        verbose = FALSE)
cell.results <- data.frame(row.names(input.matrix.final), celltypes(cell.fit))
colnames(cell.results) <- c('cell', 'annotation')
StopWatchEnd$RunCellAssign <- Sys.time()

####################################
### Write cell annotations 
####################################

StopWatchStart$WriteAnnotations <- Sys.time()
OutDirAnnos <- paste0(Tempdir, '/ANNOTATIONS/', PrefixOutfiles, '.', ProgramOutdir,'_Annotations.tsv')
write.table(cell.results, OutDirAnnos, sep='\t', quote = FALSE, row.names = FALSE)
StopWatchEnd$WriteAnnotationss <- Sys.time()


####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste0(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt")
#write(file = OutfileCPUusage, x = paste("Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
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

### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(Tempdir)
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    if (DirName == "FILTERED_DATA_MATRICES" | DirName == "UNFILTERED_DATA_MATRICES") {
      sapply(list.dirs(TempdirWithData, full.names = F, recursive = F), FUN=function(SubDirName) {
        OutdirFinal <- gsub(pattern = Tempdir, replacement =  paste0(Outdir, "/", ProgramOutdir), x = paste0(TempdirWithData, "/", SubDirName))
        dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
        sapply(list.files(paste0(TempdirWithData, "/", SubDirName), pattern = ".gz", full.names = F), FUN=function(EachFileName) {
          file.copy(from=paste0(TempdirWithData, "/", SubDirName, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
          file.remove(paste0(TempdirWithData, "/", SubDirName, "/", EachFileName))
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
OutfileCPUusage <- gsub(x = OutfileCPUusage, pattern = Tempdir, replacement = Outdir)
writeLines(paste0("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report"))

quit()
