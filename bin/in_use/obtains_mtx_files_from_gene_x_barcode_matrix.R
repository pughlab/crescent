####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script that takes a genes (rows) vs. barcodes (columns) sparse matrix
### and tranforms it into a mtx set of files (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz)
### or (prefix_outfiles.mtx, prefix_outfiles.mtx_rows and prefix_outfiles.mtx_cols)
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(DropletUtils)) # (Bioconductor) to handle MTX/H5 format files. Note it has about the same speed than library(earlycross) which can't handle H5
### Requires DropletUtils (tested on v1.7.1), which can be installed like:
### install.packages("remotes")
### remotes::install_github("MarioniLab/DropletUtils")
suppressPackageStartupMessages(library(Seurat))       # to create a Seurat object for earlycross
### Requires Seurat v3 (tested on v3.0.3.9023), which can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "obtains_mtx_files_from_gene_x_barcode_matrix.R"
ProgramOutdir  <- "GetMTX"

####################################
### Get inputs from command line argumets
####################################
#
option_list <- list(
  make_option(c("-i", "--input"), default="NA",
              help="Path/name to a genes (rows) vs. barcodes (columns) matrix
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-o", "--outdir"), default="selected_gene_bc_matrices",
              help="A path/name for the directory where the new mtx files will be saved
                Default = 'selected_gene_bc_matrices'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Note: if using option `-l y`, outfile names will be:
                prefix_outfiles.mtx, prefix_outfiles.mtx_rows and prefix_outfiles.mtx_cols
                Or if using `-l n`, outfile names will be:
                'barcodes.tsv.gz'  'features.tsv.gz' and 'matrix.mtx.gz'
                Default = 'NA'"),
  #
  make_option(c("-l", "--add_barcode_and_gene_numbers"), default="N",
              help="Indicates if 'barcodes.tsv' and 'genes.tsv' outfiles should have numbers in the first column (type [y/Y] or [n/N]), like:
                1	SRR3541565
                2	SRR3541564
                3	SRR3541563
                and
                1	ENSG00000000003
                2	ENSG00000000005
                3	ENSG00000000419
                Default = 'N'"),
  #
  make_option(c("-w", "--run_cwl"), default="0",
              help="Indicates if this script should produce 'frontend' files for crescent.cloud
                0 = no frontend files should be produced
                1 = frontend files should be produced and '--minio_path path' is provided
                2 = frontend files should be produced but '--minio_path path' is not provided (i.e local run)

                Default = '0'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system

                Default = 'MAX'"),
  #
  make_option(c("-a", "--max_global_variables"), default="10000",
              help="Indicates maximum allowed total size (in bytes) of global variables identified. Used by library(future) to prevent too large exports

                Default = '10000' for 10000 MiB")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input              <- opt$input
Outdir             <- opt$outdir
PrefixOutfiles     <- opt$prefix_outfiles
AddNumbers         <- opt$add_barcode_and_gene_numbers
NumbCores          <- opt$number_cores
MaxGlobalVariables <- as.numeric(opt$max_global_variables)
RunsCwl            <- as.numeric(opt$run_cwl)

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
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

FILE_TYPE_OUT_DIRECTORIES = c(
  "MTX",
  "LOG_FILES"
)

if (RunsCwl == 1) {
  ### Using `-w 1` will make Tempdir, which takes the value of ProgramOutdir, and it will be the final out-directory
  ### for most outfiles, except R objects, which will be written into R_OBJECTS_CWL
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F)
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

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_UsedOptions", ".txt")

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

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo", ".txt")
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
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Load data
####################################

StopWatchStart$LoadData <- Sys.time()

print("Loading genes (rows) vs. barcodes (columns) matrix")

## Note `check.names = F` is needed for both `fread` and `data.frame`
input.matrix <- as.matrix(data.frame(fread(Input, check.names = F), row.names=1, check.names = F))

StopWatchEnd$LoadData <- Sys.time()

StopWatchStart$CreateSeuratObject <- Sys.time()

### Create Seurat object
print("Creating Seurat object")
seurat.object  <- CreateSeuratObject(counts = input.matrix, project = PrefixOutfiles)

StopWatchEnd$CreateSeuratObject <- Sys.time()

OutdirMTX<-paste0(Tempdir, "/MTX/", PrefixOutfiles)
dir.create(file.path(OutdirMTX), showWarnings = F, recursive = T)

### Remove preexisting files
MtxFilesList <- list("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
for (file in MtxFilesList) {
  if (file.exists(paste(OutdirMTX, "/", file, sep = "", collapse = ""))) {
  system(command = paste("rm ", OutdirMTX, "/", file, sep = "", collapse = ""))
  }
}

####################################
### Write outfiles
####################################

StopWatchStart$WriteOutfiles <- Sys.time()

print("Writing MTX files")
write10xCounts(path = OutdirMTX, x = seurat.object@assays[["RNA"]]@data, gene.type="Gene Expression", overwrite=T, type="sparse", version="3")

StopWatchEnd$WriteOutfiles <- Sys.time()

####################################
### Reformat (if needed)
####################################
if (grepl(pattern = "y", ignore.case = T, x = AddNumbers) == T) {

  StopWatchStart$ReformatOutfiles <- Sys.time()
  barcodes  <- read.table(file = paste(OutdirMTX, "/barcodes.tsv.gz", sep = "", collapse = ""), row.names = 1)
  features  <- read.table(file = paste(OutdirMTX, "/features.tsv.gz", sep = "", collapse = ""), row.names = 2)
  featuresOutfile <- paste0(OutdirMTX, "/", PrefixOutfiles, ".expression_tpm.mtx_rows.gz")
  barcodesOutfile <- paste0(OutdirMTX, "/", PrefixOutfiles, ".expression_tpm.mtx_cols.gz")
  matrixOutfile   <- paste0(OutdirMTX, "/", PrefixOutfiles, ".expression_tpm.mtx.gz")
  #
  write(file = barcodesOutfile, x = paste(1:length(row.names(barcodes)), row.names(barcodes), sep = "\t", collapse = "\n"))
  write(file = featuresOutfile, x = paste(1:length(row.names(features)), row.names(features), sep = "\t", collapse = "\n"))
  #
  file.copy(from=paste(OutdirMTX, "/matrix.mtx.gz", sep = "", collapse = ""), to=matrixOutfile, overwrite=T)
  #
  file.remove(paste0(OutdirMTX, "/features.tsv.gz"))
  file.remove(paste0(OutdirMTX, "/barcodes.tsv.gz"))
  file.remove(paste0(OutdirMTX, "/matrix.mtx.gz"))
  FilesMTX <- c(paste0(PrefixOutfiles, ".expression_tpm.mtx_rows.gz"),
                paste0(PrefixOutfiles, ".expression_tpm.mtx_cols.gz"),
                paste0(PrefixOutfiles, ".expression_tpm.mtx.gz"))
  StopWatchEnd$ReformatOutfiles <- Sys.time()
}else{
  featuresOutfile <- paste0(OutdirMTX, "/features.tsv.gz")
  barcodesOutfile <- paste0(OutdirMTX, "/barcodes.tsv.gz")
  matrixOutfile   <- paste0(OutdirMTX, "/matrix.mtx.gz")
  FilesMTX <- c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_CPUtimes", ".txt")
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
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    if (DirName == "MTX") {
      OutdirFinal <- gsub(pattern = Tempdir, replacement = Outdir, x = OutdirMTX)
      dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
      sapply(FilesMTX, FUN=function(EachFileName) {
        file.copy(from=paste0(OutdirMTX, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
        file.remove(paste0(TempdirWithData, "/", EachFileName))
        print(paste0(OutdirMTX, "/", EachFileName))
        print(paste0(OutdirFinal, "/", EachFileName))
      })
    }else{
      OutdirFinal <- paste0(Outdir, "/", ProgramOutdir, "/", DirName)
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
