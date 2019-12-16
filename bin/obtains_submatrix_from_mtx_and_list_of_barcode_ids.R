####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Uses Seurat to load the matrix and library(DropletUtils) to write the subsampled matrix
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(DropletUtils)) # (Bioconductor) to handle reading and writing mtx files
suppressPackageStartupMessages(library(Seurat))       # (CRAN) to run QC, differential gene expression and clustering analyses
### Seurat v3 can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
####################################

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
  make_option(c("-i", "--input"), default="NA",
              help="The path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
              Notes:
              The 'MTX' files can be for example the output from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
              Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
              Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-b", "--select_barcodes"), default="NA",
              help="One of three options:
              A path/name for the *file* with the list of barcodes to select, one-per-row
              Or type a *numeric value* to print barcodes with a *sum* across all genes >= *numeric value*
              Or type 'ALL' to print all barcodes from --input into outfile
              Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-g", "--select_genes"), default="NA",
              help="One of three options:
              A path/name for the *file* with the list of genes to select, one-per-row
              Or type a *numeric value* to print genes with a *sum* across all barcodes >= *numeric value*
              Or type 'ALL' to print all genes from --input into outfile
              Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the directory where the new mtx files will be saved
              Default = 'No default. It's mandatory to specify this parameter'")
  )

opt <- parse_args(OptionParser(option_list=option_list))

Input          <- opt$input
SelectBarcodes <- opt$select_barcodes
SelectGenes    <- opt$select_genes
Outdir         <- opt$outdir
PrefixOutfiles <- "selected_gene_bc_matrices"

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

####################################
### Start stopwatches
####################################

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input", "select_barcodes", "select_genes", "outdir")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs
####################################

CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
OutdirFinal<-paste(Outdir, "/selected_gene_bc_matrices",  sep = "", collapse = "")
dir.create(file.path(OutdirFinal), recursive = T)

####################################
### Load MTX data
####################################
writeLines("\n*** Load MTX files ***\n")

StopWatchStart$LoadScRNAseqData  <- Sys.time()

print("Loading MTX infiles")
input.matrix <- Read10X(data.dir = Input)
dim(input.matrix)

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Create a Seurat object
####################################
writeLines("\n*** Create a Seurat object ***\n")

StopWatchStart$CreateSeuratObject  <- Sys.time()

seurat.object.u  <- CreateSeuratObject(counts = input.matrix, project = PrefixOutfiles)
nCellsInOriginalMatrix<-length(seurat.object.u@meta.data$orig.ident)
seurat.object.u

StopWatchEnd$CreateSeuratObject  <- Sys.time()

####################################
### Determine barcodes to subsample
####################################
writeLines("\n*** Determine barcodes to subsample ***\n")

StopWatchStart$DetermineBarcodesToSubsample  <- Sys.time()

if ((grepl(pattern = "^[0-9]+$", x = SelectBarcodes)) == TRUE) {
  stop("Subsetting barcodes by number needs to be implemented")
}else if ((grepl(pattern = "^ALL$", x = SelectBarcodes)) == TRUE) {
  stop("Including ALL barcodes needs to be implemented")
}else{
  sampled.barcodes <- data.frame(read.table(SelectBarcodes, header = F, row.names = NULL))
}

StopWatchEnd$DetermineBarcodesToSubsample  <- Sys.time()

####################################
### Determine genes to subsample
####################################
writeLines("\n*** Determine genes to subsample ***\n")

StopWatchStart$DetermineGenesToSubsample  <- Sys.time()

if ((grepl(pattern = "^[0-9]+$", x = SelectGenes)) == TRUE) {
  stop("Subsetting genes needs to be implemented")
}else if ((grepl(pattern = "^ALL$", x = SelectGenes)) == TRUE) {
  print("Using all genes")
}else{
  stop("Subsetting genes needs to be implemented")
}

StopWatchEnd$DetermineGenesToSubsample  <- Sys.time()

####################################
### Get subset of barcodes and genes from Seurat object
####################################
writeLines("\n*** Get subset of barcodes and genes from Seurat object ***\n")

StopWatchStart$SubsampleMatrix  <- Sys.time()

seurat.object.subsampled <- SubsetData(object = seurat.object.u, cells = as.vector(sampled.barcodes[,1]))

StopWatchEnd$SubsampleMatrix  <- Sys.time()

####################################
### Write downsubsampled matrix
####################################
writeLines("\n*** Write downsubsampled matrix ***\n")

StopWatchStart$WriteSubsampledMatrix  <- Sys.time()

write10xCounts(path = OutdirFinal, x = seurat.object.subsampled@assays[["RNA"]]@data, gene.type="Gene Expression", overwrite=T, type="sparse", version="3")

StopWatchEnd$WriteSubsampledMatrix  <- Sys.time()

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".Subsample_UsedOptions.txt", sep="")
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

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".Subsample_CPUusage.txt", sep="")

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
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"SUBSAMPLED/",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".Subsample_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/SUBSAMPLED/",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

print(paste("Outfiles at: ", Outdir, "SUBSAMPLED/", sep = "", collapse = ""))

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
writeLines(paste("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report", sep = "", collapse = ""))

quit()
