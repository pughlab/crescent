####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))       # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))     # to read tables quicker than read.table
suppressPackageStartupMessages(library(future))         # to run parallel processes
suppressPackageStartupMessages(library(scales))         # to use opacity in scatter plots
suppressPackageStartupMessages(library(Seurat))         # to load MTX format scRNA-seq data
### Requires Seurat v3 (tested on v3.0.3.9023), which can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")

### The following packages can be installed like:
### `BiocManager::install(c("doMC", "doRNG" ... etc))`
suppressPackageStartupMessages(library(AUCell))         # to run AUCell
suppressPackageStartupMessages(library(GSEABase))       # to load gene sets from GMT file
suppressPackageStartupMessages(library(doMC))           # to support paralell execution
suppressPackageStartupMessages(library(doRNG))          # to support paralell execution


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
  make_option(c("-c", "--infile_gmt"), default="NA",
              help="A path/name to a <tab> delimited *file* of gene sets in *gmt format, like:
                GeneSet1_ID  GeneSet1_Name  Gene1 Gene2 Gene3
                GeneSet2_ID  GeneSet2_Name  Gene4 Gene5
                ... etc
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-d", "--infile_list_dim_red_files"), default="NA",
              help="A path/name to a <tab> delimited *file* with a *list* of dimention reduction coordinate files, like:
                UMAP   path_to/infile_umap.tsv
                TSNE   path_to/infile_tsne.tsv
                ... etc
                Or type 'NA' to skip this option (no plots will be produced)
                Default = 'NA'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-a", "--opacity"), default="0.2",
              help="Provides a value for the minimal opacity of gene set prediction. Use a value between 0 and 1
                Default = '0.2'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system
                Default = 'MAX'"),
  #
  make_option(c("-w", "--run_cwl"), default="N",
              help="Indicates if this script is running inside a virtual machine container, such that outfiles are written directly into the 'HOME' . Type 'y/Y' or 'n/N'.
                Note, if using 'y/Y' this supersedes option -o
                Default = 'N'")
)

opt <- parse_args(OptionParser(option_list=option_list))

InputCounts             <- opt$input_counts
InputCountsType         <- opt$input_counts_type
InfileGmt               <- opt$infile_gmt
InfileListDimRed        <- opt$infile_list_dim_red_files
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
Opacity                 <- as.numeric(opt$opacity)
NumbCores               <- opt$number_cores
RunsCwl                 <- opt$run_cwl

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

ProgramOutdir <- "AUCELL"

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
### Define number of cores and RAM for parallelization
####################################

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (is.numeric(NumbCores) == T) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

cat("Using ", NumbCoresToUse, "cores")

plan(strategy = "multicore", workers = NumbCoresToUse)

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4Gb of RAM (4000*1024^2), use:
options(future.globals.maxSize = 4000 * 1024^2)

####################################
### Define default parameters
####################################

DefaultParameters <- list(

  ### Parameters for dimmension reduction plots
  BaseSizeSinglePlotPdf  = 7,
  BaseSizeSinglePlotPng  = 480,
  BaseSizeMultiplePlotPdfWidth  = 3.7,
  BaseSizeMultiplePlotPdfHeight = 3,
  nColDimRedPlot = 4
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

ListMandatory<-list("input_counts", "input_counts_type", "input_annotations", "input_gene_order", "outdir", "prefix_outfiles")
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
  stop(paste("Unexpected type of --input_counts: ", InputCountsType, "\n\nFor help type:\n\nRscript Runs_InferCNV.R -h\n\n", sep=""))
}

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Load gene sets
####################################
writeLines("\n*** Load gene sets ***\n")

StopWatchStart$LoadGeneSets  <- Sys.time()

geneSets <- getGmt(InfileGmt)

StopWatchEnd$LoadGeneSets  <- Sys.time()

####################################
### Obtain AUCell rankings
####################################
writeLines("\n*** Obtain AUCell rankings ***\n")

StopWatchStart$ObtainAucellRankings  <- Sys.time()

cells_rankings <- AUCell_buildRankings(exprMatrix)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_NumberOfGenesDetectedPerCell.pdf", sep=""), width=DefaultParameters$BaseSizeSinglePlotPdf, height=DefaultParameters$BaseSizeSinglePlotPdf)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=NumbCoresToUse, plotStats=TRUE)
dev.off()

StopWatchEnd$ObtainAucellRankings  <- Sys.time()

####################################
### Calculate AUCell AUC
####################################
writeLines("\n*** Calculate AUCell AUC ***\n")

StopWatchStart$ObtainAucellRankings  <- Sys.time()

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

### Define number of row and columns for plot layout
if (length(geneSets) %% DefaultParameters$nColDimRedPlot == 0 ) {
  nRowDimRedPlot <- as.integer(length(geneSets)/DefaultParameters$nColDimRedPlot)
}else{
  nRowDimRedPlot <- as.integer(length(geneSets)/DefaultParameters$nColDimRedPlot)+1
}

pdfWidth  <- DefaultParameters$nColDimRedPlot * DefaultParameters$BaseSizeMultiplePlotPdfWidth
pdfHeight <- nRowDimRedPlot * DefaultParameters$BaseSizeMultiplePlotPdfHeight

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_AUCs.pdf", sep=""), width=pdfWidth, height=pdfHeight)
par(mfrow=c(nRowDimRedPlot,DefaultParameters$nColDimRedPlot))
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
dev.off()

StopWatchEnd$ObtainAucellRankings  <- Sys.time()

####################################
### Assign gene sets to cells
####################################
writeLines("\n*** Assign gene sets to cells ***\n")

StopWatchStart$AssignGeneSetsToCells  <- Sys.time()

### Assign gene sets to cells
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "GeneSet_assignment"
assignmentMat <- table(assignmentTable[,"GeneSet_assignment"], assignmentTable[,"cell"])

### Get assignment thersholds
selectedThresholds <- getThresholdSelected(cells_assignment)

### Write out gene set assignments
OutfileAssignments <- paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_GeneSetAssignments.tsv",sep="")
Headers<-paste("Cell_barcode","GeneSet_assignment", sep = "\t", collapse = "")
write.table(Headers, OutfileAssignments, row.names = F, sep="", quote = F, col.names = F)
write.table(data.frame(assignmentTable),file = OutfileAssignments, row.names = F, col.names = F, sep="\t", quote = F, append = T)

### Write out gene set assignment thersholds
OutfileThersholds <- paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_GeneSetAssignmentThresholds.tsv",sep="")
write.table(paste(t(names(selectedThresholds)),data.frame(t(selectedThresholds)), sep = "\t", collapse = "\n"), OutfileThersholds, row.names = F, sep="\t", quote = F, col.names = F)

StopWatchEnd$AssignGeneSetsToCells  <- Sys.time()

####################################
### Generate dimention reduction coloured plots
####################################

if (regexpr("^NA$", InfileListDimRed, ignore.case = T)[1] == 1) {
  writeLines("\n*** Skipping generation of dimention reduction plots ***\n")
}else{
  writeLines("\n*** Generate dimention reduction coloured plots ***\n")

  StopWatchStart$DimReductionPlots  <- Sys.time()

  ### Define colours
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("grey70","grey70"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("gold", "red"))(nBreaks)
  
  ### Load list of infiles
  ListDimRedInfiles<-read.table(InfileListDimRed, header = F, row.names = 1, stringsAsFactors = FALSE)
  colnames(ListDimRedInfiles)<-c("PathToDataset")
  
  for (dataset in rownames(ListDimRedInfiles)) {
    PathToDataset<-ListDimRedInfiles[dataset,"PathToDataset"]
    
    ## Note `check.names = F` is needed for both `fread` and `data.frame`
    cellsTsne <- as.matrix(data.frame(fread(PathToDataset, check.names = F), row.names=1, check.names = F))

    ### Generate plots
    pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".",  ProgramOutdir, "_", dataset, "Plots.pdf", sep=""), width=pdfWidth, height=pdfHeight)
    par(mfrow=c(nRowDimRedPlot,DefaultParameters$nColDimRedPlot))

    for(geneSetName in names(selectedThresholds)) {
    
      # Split cells according to their AUC value for the gene set
      passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
      if(sum(passThreshold) > 0 ) {
        aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
        # Assign cell color
        cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])),
                       setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    
        # Plot
        plot(cellsTsne, main=geneSetName,
             sub="Yellow/red cells pass the threshold",
             col=alpha(cellColor[rownames(cellsTsne)], Opacity), pch=16)
      }
    }
    dev.off()
  }
  StopWatchEnd$DimReductionPlots  <- Sys.time()
}
  
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
write(file = OutfileCPUusage, x = paste("Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
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