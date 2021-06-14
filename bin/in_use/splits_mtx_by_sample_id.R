####################################
### Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))       # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(Seurat))         # (CRAN) to subset matrix
suppressPackageStartupMessages(library(DropletUtils))   # (Bioconductor) to write mtx files

ThisScriptName <- "splits_mtx_by_sample_id"

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
  make_option(c("-i", "--indir_mtx"), default="NA",
              help="Path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
                The barcodes.tsv.gz file must contain barcodes with sample ID's, like:
                SampleId1_ATCGATCGATCGATCG
                SampleId1_ATTTCTCGATCGATTT
                SampleId2_AGGGGTCGATCGACCC
                SampleId2_AAAAATCGATCGAGGG
                Strings before the last '_' will be used a sample IDs

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
              
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
              
                Default = 'No default. It's mandatory to specify this parameter'")
)

opt <- parse_args(OptionParser(option_list=option_list))

InputCounts             <- opt$indir_mtx
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles

#####

opt <- parse_args(OptionParser(option_list=option_list))

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
### Define outdirs
####################################
writeLines("\n*** Create outdirs ***\n")

FILE_TYPE_OUT_DIRECTORIES = c(
  "PER_SAMPLE",
  "LOG_FILES"
)

CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Outdir<-gsub("/$", "", Outdir)
#
dir.create(file.path(Outdir), recursive = T)

sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(eachdir) {
  dir.create(file.path(paste0(Outdir, "/", eachdir)), showWarnings = F, recursive = T)
})

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

StopWatchStart$ReportUsedOptions  <- Sys.time()

OutfileOptionsUsed<-paste0(Outdir, "/LOG_FILES/", PrefixOutfiles, "_SplitMtxPerSample_UsedOptions", ".txt")

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

OutfileRSessionInfo<-paste0(Outdir, "/LOG_FILES/", PrefixOutfiles, "_SplitMtxPerSample_RSessionInfo", ".txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("indir_mtx", "outdir", "prefix_outfiles")
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

print("Loading MTX infiles")
exprRead10X <- Read10X(data.dir = InputCounts)
sampleIDs_ <- unique(regmatches(colnames(exprRead10X), regexpr(".+_",colnames(exprRead10X))))

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Write matrices for each sample ID
####################################
writeLines("\n*** Write matrices for each sample ID ***\n")

StopWatchStart$WriteSubsampledMatrix  <- Sys.time()

sapply(sampleIDs_, FUN=function(sampleID_) {
  sampleID <- gsub(pattern = "_$", replacement = "",  x = sampleID_)
  OutdirSampleID <- paste0(Outdir, "/PER_SAMPLE/", sampleID)
  write10xCounts(path = OutdirSampleID,
                 x = exprRead10X[,grepl(pattern = sampleID_, x = colnames(exprRead10X))],
                 gene.type="Gene Expression", overwrite=T, type="sparse", version="3")
})

StopWatchEnd$WriteSubsampledMatrix  <- Sys.time()

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(Outdir, "/LOG_FILES/", PrefixOutfiles, "_SplitMtxPerSample_CPUusage.txt")

Headers<-paste("Step", "Time(minutes)", sep="\t")
write.table(Headers,file = OutfileCPUtimes, row.names = F, col.names = F, sep="\t", quote = F, append = T)

for (stepToClock in names(StopWatchStart)) {
  TimeStart <- StopWatchStart[[stepToClock]]
  TimeEnd   <- StopWatchEnd[[stepToClock]]
  TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
  ReportTime<-c(paste(stepToClock, TimeDiff, sep = "\t", collapse = ""))
  write(file = OutfileCPUtimes, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
}

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
writeLines(paste0("END - All done!!! See:\n", OutfileCPUtimes, "\nfor computing times report"))

quit()
