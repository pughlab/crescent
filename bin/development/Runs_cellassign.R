####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Erik Christensen - echris3@uwo.ca
### Script made based on #TODO fill in some sources
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) #TODO write an overview
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
suppressPackageStartupMessages(library(cellassign))   # TODO get other includes and describe them
#TODO see if i acutally need these
#suppressPackageStartupMessages(library(Seurat))       # (CRAN) to run QC, differential gene expression and clustering analyses
#suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
#suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
#suppressPackageStartupMessages(library(fmsb))         # (CRAN) to calculate the percentages of extra properties to be t-SNE plotted
#suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
#suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
#suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
#suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
#suppressPackageStartupMessages(library(staplr))       # (CRAN) to merge pdf files. Note it needs pdftk available. If not available use `SummaryPlots <- "N"`
#suppressPackageStartupMessages(library(gtools))       # (CRAN) to do alphanumeric sorting. Only needed if using `-w Y`.
#suppressPackageStartupMessages(library(loomR))        # (GitHub mojaveazure/loomR) needed for fron-end display of data. Only needed if using `-w Y`.
####################################

####################################
### Required external packages
####################################
#TODO Shouldn't be any?
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
  make_option(c("-t", "--input_type"), default="NA", #TODO make sure i can support this option
              help="Either 'MTX', 'TSV' or 'HDF5'
                'MTX'  is the path/name to an MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
                'TSV'  is the path/name of a <tab> delimited *file* with genes in rows vs. barcodes in columns
                'HDF5' is the path/name of a *file* in hdf5 format (e.g. from Cell Ranger)
                Note 'MTX' files can be the outputs from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-m", "--input_markers"), default="NA",
              help=""), #TODO write help for marker input
  #
  make_option(c("-s", "--size_factors"), default="NA",
              help="") #TODO write help for size factors
  #
  make_option(c("-x", "--design_matrix"), default="NA",
              help="") #TODO write help for design matrix (is it even required)
  #
  make_option(c(), default="NA",
              help="") #TODO fill other options
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-u", "--number_cores"), default="AUTO",
              help="Indicates one of three options:
                a) Type the number of cores to use for parellelization (e.g. '4')
                b) Type 'AUTO' to detect the number of cells in the sample and assign a pre-established number of cores based on that
                c) Type 'MAX' to determine and use all available cores in the system
                Default = 'AUTO'"),
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
  make_option(c("-a", "--max_global_variables"), default="AUTO",
              help="Indicates one of two options:
                a) Type the maximum allowed total size (in bytes) of global variables identified for library(future) (e.g. '4000')
                b) Type 'AUTO' to detect the number of cells in the sample and assign a pre-established size of global variables
                Default = 'AUTO'")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input                   <- opt$input
InputType               <- opt$input_type
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl
MaxGlobalVariables      <- opt$max_global_variables

####################################
### Define default parameters
####################################
#TODO define defaults for cellassign (see runs_seurat_v3_single.. for examples)

DefaultParameters <- list(
  
)

#TODO do i need these definitions
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
### Start stopwatches
####################################

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")
#TODO update list of mandatory parameters for cellassign
ListMandatory<-list("input", "input_type", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
  }
}

####################################
### Define outdirs and CWL parameters
####################################
#TODO Modify the output files
writeLines("\n*** Create outdirs ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Using `-w Y` will make Tempdir, which takes the value of ProgramOutdir, and it will be the final out-directory
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) 
  #TODO
  FILE_TYPE_OUT_DIRECTORIES = c(
    "CRESCENT_CLOUD",
    "CRESCENT_CLOUD/frontend_normalized",
    "CRESCENT_CLOUD/frontend_coordinates",
    "CRESCENT_CLOUD/frontend_raw",
    "CRESCENT_CLOUD/frontend_metadata",
    "CRESCENT_CLOUD/frontend_markers",
    "CRESCENT_CLOUD/frontend_qc",
    "CRESCENT_CLOUD/frontend_groups",
    "AVERAGE_GENE_EXPRESSION_TABLES", 
    "CELL_CLUSTER_IDENTITIES", 
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES",
    "DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS",
    "DIMENSION_REDUCTION_COORDINATE_TABLES",
    "DIMENSION_REDUCTION_PLOTS",
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS", 
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS", 
    "SUMMARY_PLOTS",
    "UNFILTERED_DATA_MATRICES"
  )
  
}else{
  PrefixOutfiles <- c(paste0(PrefixOutfiles,"_res",Resolution))
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
  #TODO
  FILE_TYPE_OUT_DIRECTORIES = c(
    "AVERAGE_GENE_EXPRESSION_TABLES", 
    "CELL_CLUSTER_IDENTITIES", 
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES",
    "DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS",
    "DIMENSION_REDUCTION_COORDINATE_TABLES",
    "DIMENSION_REDUCTION_PLOTS",
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS", 
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS", 
    "SUMMARY_PLOTS",
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

#TODO confirm that all preceding code is correct
##
## EVERYTHING ABOVE THIS IS REQUIRED FOR CRESCENT AS PER JAVIER'S EMAIL
##





####################################
### Load scRNA-seq data #TODO make sure this will work for cellassign
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

#TODO See if cellassign supports parallelisation. I'd like to keep this.
####################################
### Define number of cores for parallelization
### Note: this must run after loading scRNA-seq data to get the number of barcodes (if using `-u AUTO`)
####################################
writeLines("\n*** Define number of cores for parallelization ***\n")

NumbCoresAvailable <- as.numeric(availableCores()[[1]]) ## Number of cores available in the system
NumberOfBarcodes <- ncol(input.matrix)

### Get number of cores requested
if (regexpr("^AUTO$", NumbCores, ignore.case = T)[1] == 1) {
  if (NumberOfBarcodes <= DefaultParameters$MaxNumbCellsSmallForNumbCores) {
    NumbCoresRequested <-DefaultParameters$NumbCoresSmall
  }else if (NumbCoresAvailable < DefaultParameters$NumbCoresMedOrLarge) {
    NumbCoresRequested <- NumbCoresAvailable
  }else{
    NumbCoresRequested <-DefaultParameters$NumbCoresMedOrLarge
  }
}else if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresRequested <- NumbCoresAvailable
}else if (regexpr("^[0-9]+$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresRequested <- as.numeric(NumbCores)
}else{
  stop(paste0("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript", ThisScriptName, " -h\n\n"))
}

### Check if number of cores requested is not larger than number of cores available
if (NumbCoresAvailable < NumbCoresRequested) {
  print(paste0("WARNING: Parameter `-u` requested ", NumbCoresRequested, " cores to process ", NumberOfBarcodes, " cells. But only ", NumbCoresAvailable, " cores are available and they will be used instead"))
  NumbCoresToUse <- NumbCoresAvailable
}else{
  NumbCoresToUse <- NumbCoresRequested
}

writeLines(paste0("\nUsing ", NumbCoresToUse, " cores\n"))

#TODO see above, might have to delete this
####################################
### Define library(future) strategy for parallelization
### Note: this must run after:
###       a) loading scRNA-seq data to get the number of barcodes (if using `-a AUTO`)
###       b) defining `NumbCoresToUse`
####################################
writeLines("\n*** Define library(future) strategy for parallelization ***\n")

NumberOfBarcodes <- ncol(input.matrix)

if (regexpr("^AUTO$", MaxGlobalVariables, ignore.case = T)[1] == 1) {
  if (NumberOfBarcodes <= DefaultParameters$MaxNumbCellsSmallForGlobVars) {
    MaxGlobalVariablesToUse <- DefaultParameters$NumbGlobVarsSmall
  }else{
    MaxGlobalVariablesToUse <- DefaultParameters$NumbGlobVarsMedOrLarge
  }
}else if (regexpr("^[0-9]+$", MaxGlobalVariables, ignore.case = T)[1] == 1) {
  MaxGlobalVariablesToUse <- as.numeric(MaxGlobalVariables)
}else{
  stop(paste0("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript ", ThisScriptName, " -h\n\n"))
}

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4Gb of RAM (4000*1024^2), use:

writeLines(paste0("\nUsing ", MaxGlobalVariablesToUse, " global variables\n"))

options(future.globals.maxSize = MaxGlobalVariablesToUse * 1024^2)

plan(strategy = "multicore", workers = NumbCoresToUse)







##
## EVERYTHING AFTER THIS IS REQUIRED
##
#TODO see if any of this stuff needs changing

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste0(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt")
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
