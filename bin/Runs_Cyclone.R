####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://rdrr.io/bioc/scran/man/cyclone.html
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))       # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))     # to read tables quicker than read.table
suppressPackageStartupMessages(library(Seurat))         # to read MTX files
### Seurat v3 can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
suppressPackageStartupMessages(library(scran))          # (Bioconductor) to run Cyclone
### BiocManager::install("scran")
### used v1.14.5

####### These commands can be used to save Cyclone's pretrained human and mouse pairs
### #Getting pre-trained marker sets
### suppressPackageStartupMessages(library(scran))
### mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
### hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
### 
### Headers<-paste("Phase","first","second", sep = "\t")
### 
### OutfilePreTrained <- "~/PROGRAMS/CYCLONE/pretrained_mouse_ENSMUSG_pairs.tsv"
### write.table(Headers,file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F)
### write.table(data.frame("G1",  mm.pairs$G1),  file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F, append = T)
### write.table(data.frame("S",   mm.pairs$S),   file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F, append = T)
### write.table(data.frame("G2M", mm.pairs$G2M), file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F, append = T)
### 
### OutfilePreTrained <- "~/PROGRAMS/CYCLONE/pretrained_human_ENSG_pairs.tsv"
### write.table(Headers,file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F)
### write.table(data.frame("G1",  hs.pairs$G1),  file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F, append = T)
### write.table(data.frame("S",   hs.pairs$S),   file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F, append = T)
### write.table(data.frame("G2M", hs.pairs$G2M), file = OutfilePreTrained, row.names = F, col.names = F, sep="\t", quote = F, append = T)

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
  make_option(c("-i", "--input_counts"), default="NA",
              help="Either the path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files;
                or path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns
                Notes:
                The 'MTX' files can be for example the output from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--input_counts_type"), default="NA",
              help="Indicates if --input_counts is either a 'MTX' directory or a 'DGE' file
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-j", "--input_training"), default="NA",
              help="Indicates the path to a list of training pairs in format like:
                Phase  first   second
                G1     PSMC6   GNAI3
                S      BAX     SLC2A8
                G2M    PARP3   SETDB2
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
  make_option(c("-w", "--run_cwl"), default="N",
              help="Indicates if this script is running inside a virtual machine container, such that outfiles are written directly into the 'HOME' . Type 'y/Y' or 'n/N'.
                Note, if using 'y/Y' this supersedes option -o
                Default = 'N'")
)

opt <- parse_args(OptionParser(option_list=option_list))

InputCounts             <- opt$input_counts
InputCountsType         <- opt$input_counts_type
InputTraining           <- opt$input_training
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
RunsCwl                 <- opt$run_cwl

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

ProgramOutdir <- "CYCLONE"

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Outfiles will be stored into `ProgramOutdir` directory
  #PrefixOutfiles <- "cwl_run" 
  PrefixOutfiles  <- opt$prefix_outfiles
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) ## Note Tempdir will be the final out-directory as well
}else{
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

ListMandatory<-list("input_counts", "input_counts_type",  "input_training", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Loading training pairs
####################################
writeLines("\n*** Loading training pairs ***\n")

StopWatchStart$LoadTrainingPairs  <- Sys.time()

## Note `check.names = F` is needed for both `fread` and `data.frame`
training_pairs.df <- read.table(file = InputTraining, header = T, row.names = NULL)

training_pairs.ls <- list()
for (CellStage in unique(training_pairs.df[,"Phase"])) {
  CellStage.df<-subset(x = training_pairs.df, Phase == CellStage)
  training_pairs.ls[[CellStage]] <- CellStage.df[,c("first","second")]
  print(paste("Loaded ", CellStage, " (", nrow(CellStage.df), " pairs)", collapse = "", sep = ""))
}

StopWatchEnd$LoadTrainingPairs  <- Sys.time()

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

StopWatchStart$LoadScRNAseqData  <- Sys.time()

if (regexpr("^MTX$", InputCountsType, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  exprMatrix <- as.matrix(as.data.frame(Read10X(data.dir = InputCounts)))
}else if (regexpr("^DGE$", InputCountsType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  ## Note `check.names = F` is needed for both `fread` and `data.frame`
  exprMatrix <- as.matrix(data.frame(fread(InputCounts, check.names = F), row.names=1, check.names = F))
}else{
  stop(paste("Unexpected type of --input_counts: ", InputCountsType, "\n\nFor help type:\n\nRscript Runs_Cyclone.R -h\n\n", sep=""))
}

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Classify test dataset cells
####################################
writeLines("\n*** Classify test dataset cells ***\n")

StopWatchStart$ClassifyTestCells  <- Sys.time()

# Classify test dataset cells
assignments.cyclone <- cyclone(x = exprMatrix, training_pairs.ls)
rownames(assignments.cyclone$scores)<-colnames(exprMatrix)

StopWatchEnd$ClassifyTestCells  <- Sys.time()

####################################
### Save outfiles
####################################
writeLines("\n*** Save outfiles ***\n")

StopWatchStart$SaveOutfiles  <- Sys.time()

OutfileCycloneScores      <- paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_scores.tsv", sep="")
OutfileCyclonePhases      <- paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_phases.tsv", sep="")
OutfileCyclonePhasesTable <- paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_phases_table.tsv", sep="")
#
Headers<-paste("Cell_barcode", paste(names(assignments.cyclone$scores), sep = "", collapse = "\t") ,sep="\t")
write.table(Headers,file = OutfileCycloneScores, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(assignments.cyclone$scores),file = OutfileCycloneScores, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
Headers<-paste("Cell_barcode", "Phase", sep = "\t", collapse = "")
write.table(Headers,file = OutfileCyclonePhases, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste(colnames(exprMatrix), assignments.cyclone$phases, sep = "\t"),file = OutfileCyclonePhases, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
Headers<-paste("Phase", "Number_of_cells", sep = "\t", collapse = "")
write.table(Headers,file = OutfileCyclonePhasesTable, row.names = F, col.names = F, sep="\t", quote = F)
write.table(as.matrix(table(assignments.cyclone$phases)),file = OutfileCyclonePhasesTable, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$SaveOutfiles  <- Sys.time()

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Outdir, "/", ProgramOutdir ,"/",PrefixOutfiles,".", ProgramOutdir, "_UsedOptions.txt", sep="")

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

OutfileCPUusage<-paste(Outdir, "/", ProgramOutdir ,"/",PrefixOutfiles,".", ProgramOutdir, "_CPUusage.txt", sep="")
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