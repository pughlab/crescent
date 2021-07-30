####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Uses Seurat to load the matrix and library(DropletUtils) to write the subsampled matrix
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(DropletUtils)) # (Bioconductor) to handle MTX/H5 format files. Note it has about the same speed than library(earlycross) which can't handle H5
suppressPackageStartupMessages(library(Seurat))       # (CRAN) to subset data. Tested with v4.0.2
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "obtains_submatrix_from_mtx_and_lists_of_barcode_ids_and_gene_ids.R"

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
              (1)
              A path/name for the *file* with the list of barcodes to select, one-per-row.
              Barcode
              ATTATCCAGACTAGAT-1
              TGCTGCTAGCTCCTTC-1
              GGACAAGTCAGTGTTG-1
              Note: the 'Barcode' header may be absent
              
              (2)
              Type a *numeric value* to print barcodes with a *sum* across all genes >= *numeric value*
              
              (3)
              Type 'ALL' to print all barcodes from --input into outfile
              
              Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-g", "--select_genes"), default="NA",
              help="One of three options:
              (1)
              A path/name for the *file* with the list of genes to select, one-per-row.
              Gene
              SAMD11
              NOC2L
              KLHL17
              Note: the 'Gene' header may be absent
              
              (2)
              Type a *numeric value* to print genes with a *sum* across all barcodes >= *numeric value*
              
              (3)
              Type 'ALL' to print all genes from --input into outfile
              
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
PrefixOutfiles <- "selected_feature_bc_matrix"

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

ListMandatory<-list("input", "select_barcodes", "select_genes", "outdir")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
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
  "LOG_FILES"
)

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
dir.create(file.path(Outdir), showWarnings = F, recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(eachdir) {
  dir.create(file.path(paste0(Tempdir, "/", eachdir)), showWarnings = F, recursive = T)
})

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

StopWatchStart$ReportUsedOptions  <- Sys.time()

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, "_SubmatrixFromMtxAndListOfBarcodes_UsedOptions", ".txt")

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

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, "_SubmatrixFromMtxAndListOfBarcodes_RSessionInfo", ".txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

####################################
### Load MTX data
####################################
writeLines("\n*** Load MTX files ***\n")

StopWatchStart$LoadScRNAseqData  <- Sys.time()

print("Loading MTX infiles")
input.matrix <- Read10X(data.dir = Input, strip.suffix = T) ### Note `strip.suffix = T` applies to Seurat v3.2 or higher
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
### Determine barcodes to subsample and subsamples Seurat object
####################################
writeLines("\n*** Determine barcodes to subsample from Seurat object ***\n")

StopWatchStart$DetermineBarcodesToSubsampleAndSubsample  <- Sys.time()

if ((grepl(pattern = "^[0-9]+$", x = SelectBarcodes)) == TRUE) {
  stop("Subsetting barcodes by number of reads needs to be implemented")
}else if ((grepl(pattern = "^ALL$", x = SelectBarcodes)) == TRUE) {
  seurat.object.subsampled.barcodes <- seurat.object.u
}else{
  sampled.barcodes <- data.frame(read.table(SelectBarcodes, header = F, row.names = NULL, check.names = FALSE))
  
  # This is because we are using Read10X(..., strip.suffix = T), which removes the last '-digit' from barcode ID's
  if ((grepl(pattern = "^barcode", ignore.case = T, x = sampled.barcodes[1,"V1"])) == TRUE) {
    sampled.barcodes <- sampled.barcodes[-1,]
    sampled.barcodes<-gsub(x =sampled.barcodes, pattern = "-[0-9]+$", perl = T, replacement = "")
  }else{
    sampled.barcodes<-gsub(x = as.factor(sampled.barcodes[,1]), pattern = "-[0-9]+$", perl = T, replacement = "")
  }
  
  seurat.object.subsampled.barcodes <- subset(x = seurat.object.u, cells = as.vector(sampled.barcodes))
}

StopWatchEnd$DetermineBarcodesToSubsampleAndSubsample  <- Sys.time()

####################################
### Determine genes to subsample and subsample Seurat object
####################################
writeLines("\n*** Determine genes to subsample from Seurat object ***\n")

StopWatchStart$DetermineGenesToSubsampleAndSubsample  <- Sys.time()

if ((grepl(pattern = "^[0-9]+$", x = SelectGenes)) == TRUE) {
  stop("Subsetting genes by number of reads needs to be implemented")
}else if ((grepl(pattern = "^ALL$", x = SelectGenes)) == TRUE) {
  seurat.object.subsampled.barcodes.genes <- seurat.object.subsampled.barcodes
}else{
  sampled.genes <- data.frame(read.table(SelectGenes, header = F, row.names = NULL, check.names = FALSE, sep = "\t"))
  sampled.genes <- as.character(sampled.genes[,1])
  sampled.genes <- gsub(x=sampled.genes, pattern = "_", replacement = "-") ## Because CreateSeuratObject() will do the same
  print(paste0("Requested ", length(sampled.genes), " genes"))
  seurat.object.subsampled.barcodes.genes <- subset(x = seurat.object.subsampled.barcodes, features = sampled.genes)
}

StopWatchEnd$DetermineGenesToSubsampleAndSubsample  <- Sys.time()

####################################
### Write subsampled matrix
####################################
writeLines("\n*** Write subsampled matrix ***\n")

StopWatchStart$WriteSubsampledMatrix  <- Sys.time()

OutdirFinal <- paste0(Outdir, "/", PrefixOutfiles)
write10xCounts(path = OutdirFinal, x = seurat.object.subsampled.barcodes.genes@assays[["RNA"]]@data, gene.type="Gene Expression",
               overwrite=T, type="sparse", version="3")

StopWatchEnd$WriteSubsampledMatrix  <- Sys.time()

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, "_SubmatrixFromMtxAndListOfBarcodes_CPUtimes", ".txt")

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
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"/", PrefixOutfiles, sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir, pattern = paste(PrefixOutfiles, "_SubmatrixFromMtxAndListOfBarcodes_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  if (Outdir == Tempdir) {
    ### Do nothing
  }else{
    print(paste(Tempdir,"/",eachFile, sep=""))
    print(paste(Outdir,"/", eachFile, sep=""))
    file.copy(from=paste(Tempdir,"/",eachFile, sep=""),to=paste(Outdir,"/", eachFile, sep=""),overwrite=T)
    file.remove(paste(Tempdir,"/",eachFile, sep=""))
  }
})

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
writeLines(paste("\nEND - All done!!! See:\n", Outdir, "/selected_feature_bc_matrix", "\nfor outfiles", sep = "", collapse = "\n\n"))

quit()
