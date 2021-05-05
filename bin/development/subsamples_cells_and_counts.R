####################################
### Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))       # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))     # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(Seurat))         # (CRAN) to subset matrix
### Seurat v3 can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
suppressPackageStartupMessages(library(DropletUtils))   # (Bioconductor) to handle reading and writing mtx files and downsample matrix

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
  make_option(c("-j", "--input_barcodes"), default="NA",
              help="Path/name to a list of barcodes, one per row, to subset from --input_counts
                Or type 'NA' to either downsample barcodes randomly using --number_of_barcodes or to include all barcodes in --input_counts
                Default = 'NA'"),
  #
  make_option(c("-b", "--random_number_of_barcodes"), default="NA",
              help="Indicates the number of barcodes to downsample randomly from --input_counts
                Note: to query a specific list of barcodes use --input_barcodes, or type 'NA' to include all barcodes in --input_counts
                Default = 'NA'"),
  #
  make_option(c("-c", "--fraction_of_reads"), default="1",
              help="Indicates the fraction of counts to subsample from --input_counts (0 to 1)
                Default = '1'"),
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
InputBarcodes           <- opt$input_barcodes
RandNumberOfBarcodes    <- opt$random_number_of_barcodes
FractionReads           <- as.numeric(opt$fraction_of_reads)
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
RunsCwl                 <- opt$run_cwl

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

ProgramOutdir <- "SUBSAMPLE"

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Outfiles will be stored into `ProgramOutdir` directory
  #PrefixOutfiles <- "cwl_run" 
  PrefixOutfiles  <- opt$prefix_outfiles
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) ## Note Tempdir will be the final out-directory as well
}else{
  ## Using `Tempdir` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## Then at the end of the script they'll be moved into `OutdirFinal`
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
  OutdirFinal<-paste(Outdir, "/", ProgramOutdir, "/selected_gene_bc_matrices",  sep = "", collapse = "")
  #
  dir.create(file.path(OutdirFinal), recursive = T)
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

ListMandatory<-list("input_counts", "input_counts_type",  "fraction_of_reads", "outdir", "prefix_outfiles")
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

if (regexpr("^MTX$", InputCountsType, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  exprMatrix <- as.matrix(as.data.frame(Read10X(data.dir = InputCounts)))
}else if (regexpr("^DGE$", InputCountsType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  ## Note `check.names = F` is needed for both `fread` and `data.frame`
  exprMatrix <- as.matrix(data.frame(fread(InputCounts, check.names = F), row.names=1, check.names = F))
}else{
  stop(paste("Unexpected type of --input_counts: ", InputCountsType, "\n\nFor help type:\n\nsubsamples_cells_and_counts.R -h\n\n", sep=""))
}

seurat.object.u  <- CreateSeuratObject(counts = exprMatrix, project = PrefixOutfiles)

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Downsample matrix
####################################

StopWatchStart$DownsampleMatrix  <- Sys.time()

if (FractionReads == 1) {
  print("Will keep identical counts than --input_counts")
  seurat.object.downsampled <- seurat.object.u
}else{
  print(paste("Downsample matrix by fraction = ", FractionReads, sep = "", collapse = ""))
  exprMatrix.Downsampled <- downsampleMatrix(x = seurat.object.u@assays$RNA@counts, prop = FractionReads, bycol = F)
  seurat.object.downsampled <- CreateSeuratObject(counts = exprMatrix.Downsampled, project = PrefixOutfiles)
}

StopWatchEnd$DownsampleMatrix  <- Sys.time()

####################################
### Subset barcodes
####################################

StopWatchStart$SubsetBarcodes  <- Sys.time()

if (regexpr("^NA$", InputBarcodes, ignore.case = T)[1] == 1) {
  if (regexpr("^NA$", RandNumberOfBarcodes, ignore.case = T)[1] == 1) {
    print("Will keep all barcodes from --input_counts")
    seurat.object.downsampled.selected <- seurat.object.downsampled
  }else{
    RandNumberOfBarcodes <- as.numeric(RandNumberOfBarcodes)
    AllBarcodes<-colnames(seurat.object.downsampled)
    SubsampleBarcodes<-sample(x = AllBarcodes, size = RandNumberOfBarcodes, replace=FALSE)
    seurat.object.downsampled.selected <- SubsetData(object = seurat.object.downsampled, cells = as.vector(SubsampleBarcodes))
  }
}else{
  print("Downsample matrix using a list of barcodes")
  barcodesToInclude <- data.frame(read.table(InputBarcodes, header = F, row.names = NULL))
  seurat.object.downsampled.selected <- SubsetData(object = seurat.object.downsampled, cells = as.vector(barcodesToInclude[,1]))
}

StopWatchEnd$SubsetBarcodes  <- Sys.time()

####################################
### Write subsampled matrix
####################################
writeLines("\n*** Write downsubsampled matrix in MTX format ***\n")

seurat.object.downsampled.selected

StopWatchStart$WriteSubsampledMatrix  <- Sys.time()

write10xCounts(path = OutdirFinal, x = seurat.object.downsampled.selected@assays[["RNA"]]@data, gene.type="Gene Expression", overwrite=T, type="sparse", version="3")

StopWatchEnd$WriteSubsampledMatrix  <- Sys.time()

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt", sep="")

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
