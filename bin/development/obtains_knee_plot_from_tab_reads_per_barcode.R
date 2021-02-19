####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Entering a table with reads per barcode obtains a knee-plot and its inflection point
### indicating barcodes to be considered cells
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/obtains_knee_plot_from_tab_reads_per_barcode.R -h'
### for help
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))   # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(dropbead))   # (GitHub rajewsky-lab/dropbead) to generate the knee-plot and get inflection
suppressPackageStartupMessages(library(data.table)) # (CRAN) to read tables quicker than read.table
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "obtains_knee_plot_from_tab_reads_per_barcode.R"
ProgramOutdir  <- "KNEEPLOT"

####################################
### Get inputs from command line argumets
####################################

option_list <- list(
  make_option(c("-i", "--infile_tab"), default="NA",
              help="A path/name to a <tab> delimited *file* with number of reads per barcode, like:
                Barcode       Reads
                CTTATGGCTTTA  203177  
                TTAGCGCTTATA  196369  
                CGGCTGCTAATC  151561  
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-m", "--max_number_of_barcodes_to_plot"), default="60000",
              help="Typically 5 to 10 times the number of expected cells in the sample
                Default = '60000'"),
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

InfileTab       <- opt$infile_tab
MaxBarcodes     <- as.numeric(opt$max_number_of_barcodes_to_plot)
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles
RunsCwl         <- opt$run_cwl

####################################
### Start stopwatches
####################################

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_tab", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
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
    "KNEEPLOT",
    "LOG_FILES"
  )
  
}else{
  ## Using `Tempdir/DIRECTORY` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## 'DIRECTORY' is one of the directories specified at FILE_TYPE_OUT_DIRECTORIES
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  Tempdir        <- "~/temp" 
  #
  CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
  UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
  #
  Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
  Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
  Outdir<-gsub("/+", "/", Outdir, perl = T)
  Tempdir<-gsub("/+", "/", Tempdir, perl = T)
  Outdir<-gsub("/$", "", Outdir)
  Tempdir<-gsub("/$", "", Tempdir)
  #
  dir.create(file.path(Outdir, ProgramOutdir), recursive = T)
  dir.create(file.path(Tempdir), showWarnings = F, recursive = T)
  
  FILE_TYPE_OUT_DIRECTORIES = c(
    "KNEEPLOT",
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

OutfileOptionsUsed<-paste(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_UsedOptions.txt", sep="")

TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=paste0("Run started: ", TimeOfRun, "\n"))

write(file = OutfileOptionsUsed, x=c("Commands used:"), append = T)
for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

cat(file = OutfileOptionsUsed, x=paste("\n", "One-line-commands used:", "\n", "`Rscript /path_to/", ThisScriptName, sep = ""), append = T)
for (optionInput in option_list) {
  cat(file = OutfileOptionsUsed, x=(paste(" ", optionInput@short_flag, " ", opt[optionInput@dest], sep="", collapse = "")),append = T)
}
cat(file = OutfileOptionsUsed, x="`\n", append = T)

StopWatchEnd$ReportUsedOptions  <- Sys.time()

####################################
### Load read-counts per cell data
####################################
writeLines("\n*** Load read-counts per cell data ***\n")

StopWatchStart$LoadReadCounts  <- Sys.time()

## Note `check.names = F` is needed for both `fread` and `data.frame`
mat <- data.frame(fread(InfileTab, check.names = F), row.names=1, check.names = F)

StopWatchEnd$LoadReadCounts  <- Sys.time()

####################################
### Get knee-plot
### Note: earlier versions of this script used library(dropbead) plotCumulativeFractionOfReads() and estimateCellNumber()
###       but they failed to produce a kneeplot for dataset: ~/SINGLE_CELL/10X/SMARTER_VACCINE/DATA_FROM_PUGHLAB/SMTR04t1_NonRad/5p_RNA
###       Thus decided to program this part over
####################################
writeLines("\n*** Get knee-plot ***\n")

StopWatchStart$ReadReadCounts  <- Sys.time()

HeaderNReads <- colnames(mat)[[1]]

barcode_number <- 0
for (i in rownames(mat)) {
  barcode_number <- barcode_number+1
  ivalue <- mat[barcode_number,HeaderNReads]
  if (barcode_number == 1) {
    jvalue <- ivalue
  }else{
    jvalue <- jvalue + ivalue
  }
  mat[barcode_number,"cumulative"] <- jvalue
}

StopWatchEnd$ReadReadCounts  <- Sys.time()

####################################
### Generate outfiles
####################################
writeLines("\n*** Generate outfiles ***\n")

StopWatchStart$GenerateOutfiles  <- Sys.time()

### Kneeplot
OutPdf<-paste0(Tempdir, "/KNEEPLOT/", PrefixOutfiles, ".", ProgramOutdir, ".kneeplot.pdf")

XaxisLenght <- min(MaxBarcodes,nrow(mat))
pdf(OutPdf)
InflectionPoint<-estimateCellNumber(mat[, HeaderNReads], max.cells = min(XaxisLenght))
plot(x = c(1:XaxisLenght), y = mat[c(1:XaxisLenght),"cumulative"], type="l",
     xlab = "Sorted barcodes", ylab = paste0("Cumulative ", HeaderNReads), lwd = 2, main = paste(PrefixOutfiles, "\n", "Inflection = ", InflectionPoint, " cells", sep = "", collapse = ""))
abline(v=InflectionPoint, col = "red", lwd = 2, lty =2)
dev.off()

### Kneeplot inflection
OutInf<-paste0(Tempdir, "/KNEEPLOT/", PrefixOutfiles, ".", ProgramOutdir, ".kneeplot_inflection.txt")

write(file = OutInf, x = InflectionPoint)

### Top barcodes based on kneeplot inflection
OutTop<-paste0(Tempdir, "/KNEEPLOT/", PrefixOutfiles, ".", ProgramOutdir, ".TopCells.csv")

write.table(file = OutTop, x = "Barcode", row.names = F, col.names = F, quote = F)
write.table(file = OutTop, x = rownames(mat)[1:InflectionPoint], row.names = F, col.names = F, quote = F, append = T)

StopWatchEnd$GenerateOutfiles  <- Sys.time()

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo.txt", sep="")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt", sep="")
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
  
  writeLines(paste(Tempdir, sep="", collapse = ""))
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    OutdirFinal <- paste0(Outdir, "/", ProgramOutdir, "/", DirName)
    print(OutdirFinal)
    dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
    sapply(list.files(TempdirWithData, pattern = paste0("^", PrefixOutfiles, ".", ProgramOutdir), full.names = F), FUN=function(eachFileName) {
      file.copy(from=paste0(TempdirWithData, "/", eachFileName), to=paste0(OutdirFinal, "/", eachFileName), overwrite=T)
      file.remove(paste0(TempdirWithData, "/", eachFileName))
    })
  })
}

####################################
### Finish
####################################

print("END - All done!!!")
quit()
