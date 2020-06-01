####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script 'Runs_GSVA.R' obtains Gene Set Variation Analysis enrichment scores
### for each column of --infile_mat vs. each gene set from --infile_gmt
### and selects a final gene set label for each column header from --infile_mat
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads a --infile_mat with arrays to be labeled in columns (e.g. cell clusters) and observations in rows (e.g. average gene expression)
###    The matrix can be a plain text file or an MTX format (see --infile_mat_type parameter documentation below)
### 2) Loads a file with classes (e.g. gene set signatures) in *gmt format
### 3) Runs GSVA
### 4) Cluster the GSVA enrichment scores matrix by similariry of prediction vectors and generates outputs
### 5) Obtains p-value and FDR and generates outputs, including *filtered.tsv
### 6) Determines a final label (class) for each array of --infile_mat
###    NOTE: currently the final label is simply the top-ranked class from --infile_gmt for each column of --infile_mat
###          i.e. p-value and FDR cutoffs are NOT used
### 7) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_GSVA.R -h'
### for help
####################################

####################################
### Dependencies:
####################################
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(parallel))     # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(data.table))   # (CRAN) to speed up reading matrices with 'fread' instead of 'read.table'
suppressPackageStartupMessages(library(GSA))          # (CRAN) to handle *gmt infile
suppressPackageStartupMessages(library(GSVA))         # (bioconductor) to run the gsva function
suppressPackageStartupMessages(library(qvalue))       # (bioconductor) to get FDR/q-values from GSVA's p-values
suppressPackageStartupMessages(library(cluster))      # (CRAN) to cluster/sort the *GSVA_enrichment_scores.tsv rows and columns
suppressPackageStartupMessages(library(Seurat))       # (CRAN) only needed if using `-t MTX` option. Otherwise you can comment this row adding # at the beginning
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_GSVA.R"

####################################
### Get inputs from command line argumets
####################################

option_list <- list(
  make_option(c("-i", "--infile_mat"), default="NA",
              help="A path/name to a <tab> delimited *file* with genes in rows and arrays (e.g. clusters or conditions) in columns, like:
                genes  clust1  clust2  clust3 ... etc
                RP113  0       0.0045  0.0008
                FAM14  0.0077  0.0175  0.0082
                NOC2L  0.0800  0.1532  0.0745
                ...etc
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--infile_mat_type"), default="DGE",
              help="Indicates either 'DGE' or 'MTX' if --infile_mat is either:
                a) a <tab> delimited digital gene expression 'DGE' *file* with genes in rows vs. cell barcodes or cell clusters in columns, or
                b) a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
                Default = 'DGE'"),
  #
  make_option(c("-c", "--infile_gmt"), default="NA",
              help="A path/name to a <tab> delimited *file* of gene sets in *gmt format, like:
                GeneSet1_ID  GeneSet1_Name  Gene1 Gene2 Gene3
                GeneSet2_ID  GeneSet2_Name  Gene4 Gene5
                ... etc
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
  make_option(c("-e", "--pvalue_cutoff"), default="0.05",
              help="This script produce a *filtered.tsv matrix with pairs passing -e and -f filters and unfiltered outfiles
                Default = 0.05"),
  #
  make_option(c("-f", "--fdr_cutoff"), default="0.1",
              help="Same as -e option, but for FDR scores
                Default = 0.1"),
  #
  make_option(c("-w", "--run_cwl"), default="N",
              help="Indicates if this script is running inside a virtual machine container, such that outfiles are written directly into the 'HOME' . Type 'y/Y' or 'n/N'.
                Note, if using 'y/Y' this supersedes option -o
                Default = 'N'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicates one of three options:
                a) Type the number of cores to use for parellelization (e.g. '4')
                b) Type 'MAX' to determine and use all available cores in the system
                Default = 'MAX'")

)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat       <- opt$infile_mat
InfileMatType   <- opt$infile_mat_type
InfileGmt       <- opt$infile_gmt
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles
PvalueCutoff    <- as.numeric(opt$pvalue_cutoff)
FdrCutoff       <- as.numeric(opt$fdr_cutoff)
RunsCwl         <- opt$run_cwl
NumbCores       <- opt$number_cores

####################################
### Start stopwatches
####################################

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

####################################
### Define default parameters
####################################

DefaultParameters <- list(
  DigitsForRound = 5 ### Number of digits to round up enrichment score values in outfiles
)

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("infile_mat", "infile_gmt", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

ProgramOutdir <- "GSVA"

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Using `-w Y` will make Tempdir, which takes the value of ProgramOutdir, and it will be the final out-directory
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) 
  
}else{
  ## Using `Tempdir/DIRECTORY` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## 'DIRECTORY' is one of the directories specified at FILE_TYPE_OUT_DIRECTORIES
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`

  ### Checking if the system follows a mirror /scratch vs. /home structure
  ### If that's the case, this script will use Tempdir at /scratch, not at /home
  ### This is because systems like SciNet (Compute Canada) using 'Slurm Workload Manager'
  ### only allow to write in /scratch when using slave nodes, not in /home
  ### Comment the following '4' lines if you have a mirror structure and still want to use /home (called /Users in Mac)
  ### If you don't have a mirror structure just ignore this
  ## Using `Tempdir/DIRECTORY` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## 'DIRECTORY' is one of the directories specified at FILE_TYPE_OUT_DIRECTORIES
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  
  Tempdir         <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  
  ScratchTempdir<-gsub("home","scratch", Tempdir)
  if (dir.exists(ScratchTempdir)[1] == 1) {
    Tempdir<-ScratchTempdir
  }
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
}

OutfileEnrichmentScores<-paste(Tempdir, "/GSVA_RESULTS/", PrefixOutfiles, ".GSVA_enrichment_scores.tsv", sep="")
OutfileEnrichScorsClust<-paste(Tempdir, "/GSVA_RESULTS/", PrefixOutfiles, ".GSVA_enrichment_scores_sorted.tsv", sep="")
OutfilePvalues         <-paste(Tempdir, "/GSVA_RESULTS/", PrefixOutfiles, ".GSVA_pvalues.tsv", sep="")
OutfileFdrvalues       <-paste(Tempdir, "/GSVA_RESULTS/", PrefixOutfiles, ".GSVA_fdr_values.tsv", sep="")
OutfileAllScores       <-paste(Tempdir, "/GSVA_RESULTS/", PrefixOutfiles, ".GSVA_all_scores_table.tsv", sep="")
OutfileFilteredES      <-paste(Tempdir, "/GSVA_RESULTS/", PrefixOutfiles, ".GSVA_filtered.tsv", sep="")
OutfileFinalLabel      <-paste(Tempdir, "/GSVA_RESULTS/", PrefixOutfiles, ".GSVA_final_label.tsv", sep="")

FILE_TYPE_OUT_DIRECTORIES = c(
  "GSVA_RESULTS", 
  "LOG_FILES"
)

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
### Define number of cores for parallelization
### Note: this must run after loading scRNA-seq data to get the number of barcodes (if using `-u AUTO`)
####################################
writeLines("\n*** Define number of cores for parallelization ***\n")

NumbCoresAvailable <- as.numeric(availableCores()[[1]]) ## Number of cores available in the system

### Get number of cores requested
if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresRequested <- NumbCoresAvailable
}else if (regexpr("^[0-9]+$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresRequested <- as.numeric(NumbCores)
}else{
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript", ThisScriptName, " -h\n\n", sep=""))
}

### Check if number of cores requested is not larger than number of cores available
if (NumbCoresAvailable < NumbCoresRequested) {
  print(paste("WARNING: Parameter `-u` requested ", NumbCoresRequested, " cores. But only ", NumbCoresAvailable, " cores are available and they will be used instead"))
  NumbCoresToUse <- NumbCoresAvailable
}else{
  NumbCoresToUse <- NumbCoresRequested
}

writeLines(paste("\nUsing ", NumbCoresToUse, " cores\n", sep = "", collapse = ""))

plan(strategy = "multicore", workers = NumbCoresToUse)

####################################
### Load data
####################################
writeLines("\n*** Load data ***\n")

StopWatchStart$LoadData  <- Sys.time()

### Load gene expression matrix
if (regexpr("^MTX$", InfileMatType, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  fullmat <- as.matrix(Read10X(data.dir = InfileMat))
}else if (regexpr("^DGE$", InfileMatType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  ## Note `check.names = F` is needed for both `fread` and `data.frame`
  fullmat <- as.matrix(data.frame(fread(InfileMat, check.names = F), row.names=1, check.names = F))
}else{
  stop(paste("Unexpected type of input: ", InfileMatType, "\n\nFor help type:\n\nRscript obtains_GSVA_for_MatrixColumns.R -h\n\n", sep=""))
}
rownames(fullmat) <- toupper(rownames(fullmat))

### Creates object gmt2 with the gene set memberships
gmt1<-GSA.read.gmt(InfileGmt)
gmt2<-list()
for (i in 1:length(gmt1[[1]])){
  tmp<-unlist(gmt1[[1]][i])
  tmp <- toupper(tmp)
  gmt2[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
}

StopWatchEnd$LoadData  <- Sys.time()

####################################
### Run GSVA
####################################
writeLines("\n*** Run GSVA ***\n")

StopWatchStart$RunGSVA  <- Sys.time()

### Get and print out GSVA enrichment scores
StartTimeGsva<-Sys.time()
EnrichmentScores<-gsva(expr=fullmat, gset.idx.list=gmt2, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=T, parallel.sz=0)
SortedRowNames<-rownames(EnrichmentScores)
SortedColNames<-colnames(EnrichmentScores)
#
EnrichmentScores<-EnrichmentScores[SortedRowNames,SortedColNames]
write.table(data.frame("ENRICHMENT"=colnames(EnrichmentScores), t(round(x=EnrichmentScores, digits = DefaultParameters$DigitsForRound))), OutfileEnrichmentScores, row.names = F,sep="\t",quote = F)

### maps gene set members from gmt2
Classes.list<-NULL
for (i in 1:nrow(EnrichmentScores)){
  tmp1<-unlist(strsplit(gmt1[[2]][match(rownames(EnrichmentScores)[i],gmt1[[3]])],"%"))[3]
  tmp2<-paste(unlist(gmt2[[match(rownames(EnrichmentScores)[i],names(gmt2))]]),collapse = ",")
  Classes.list<-rbind(Classes.list,c(tmp1,tmp2))
}
EndTimeGsva<-Sys.time()

StopWatchEnd$RunGSVA  <- Sys.time()

####################################
### Get and print Enrichment, P-values and Q-values (FDR)
####################################
writeLines("\n*** Get and print Enrichment, P-values and Q-values (FDR) ***\n")

StopWatchStart$GetPvalAndFdr  <- Sys.time()

### Originally used something like:
### qvalues<-qvalue(pvalues,lambda=seq(0.05,0.45,0.01))$lfdr
### But this was producing errors for some clusters given their p-values distributions, like:
### "Error in smooth.spline(lambda, pi0, df = smooth.df) : 
### missing or infinite values in inputs are not allowed"
### See https://github.com/StoreyLab/qvalue/issues/11 and 
### http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
### Javier Diaz replaced this by:
### qvalues<-qvalue(pvalues,pi0=1)$lfdr

HeaderCutoff<-paste(c("PassCutoff", "_p", PvalueCutoff, "_fdr", FdrCutoff) , sep="",collapse = "")
HeaderCutoff
HeadersForPandQvalues<-paste("CLASS","ColumnHeader","EnrichmentScore","p.Val","FDR", HeaderCutoff ,sep="\t",collapse = "")
write(x=HeadersForPandQvalues,file=OutfileAllScores)

for (columnNumber in 1:ncol(EnrichmentScores)){
  pvalues<-pnorm(-abs(scale(EnrichmentScores[,columnNumber])[,1]))
  qvalues<-qvalue(pvalues,pi0=1)$lfdr
  PassCutoffs<-ifelse((pvalues<=PvalueCutoff & qvalues<=FdrCutoff)==TRUE,1,0)
  concatenatedResults<-cbind(colnames(fullmat)[columnNumber], round(x=EnrichmentScores[,columnNumber], digits = DefaultParameters$DigitsForRound) , pvalues,qvalues,PassCutoffs)
  # Write out Table with CLASS ColumnHeader EnrichmentScore p.Val FDR PassCutoff
  write.table(x=data.frame("CLASS"=rownames(concatenatedResults),concatenatedResults),file=OutfileAllScores, row.names = F,sep="\t",quote = F,col.names = F,append = T)
}

### Index and print out GSVA p.Val and FDR
dataEPQ <- read.table(file = OutfileAllScores,row.names = NULL ,header = T)
#
ForPvaluesMat<- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"p.Val"])
PvaluesMat<-xtabs(z~x+y, data=ForPvaluesMat)
PvaluesMat<-PvaluesMat[SortedRowNames,SortedColNames]
HeadersForPvalues<-paste(c("PVALUES", colnames(PvaluesMat)), sep="\t", collapse = "\t")
write(x=HeadersForPvalues,file=OutfilePvalues)
write.table(x=PvaluesMat, file=OutfilePvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
#
ForFDRvaluesMat <- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"FDR"])
FdrvaluesMat<-xtabs(z~x+y, data=ForFDRvaluesMat)
FdrvaluesMat<-FdrvaluesMat[SortedRowNames,SortedColNames]
HeadersForFdrvalues<-paste(c("FDR_VALUES", colnames(FdrvaluesMat)), sep="\t", collapse = "\t")
write(x=HeadersForFdrvalues,file=OutfileFdrvalues)
write.table(x=FdrvaluesMat, file=OutfileFdrvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
#
FilteredESMatLogical<-(PvaluesMat<=PvalueCutoff & FdrvaluesMat<=FdrCutoff)
FilteredESMatLogical<-FilteredESMatLogical[SortedRowNames,SortedColNames]
FilteredESMatValues<-ifelse(FilteredESMatLogical==TRUE,EnrichmentScores,NA)
write.table(data.frame("ENRICHMENT_FILTERED"=rownames(EnrichmentScores), round(x=FilteredESMatValues, digits = DefaultParameters$DigitsForRound)),OutfileFilteredES, row.names = F,sep="\t",quote = F)

StopWatchEnd$GetPvalAndFdr  <- Sys.time()

####################################
### Sort *.GSVA_enrichment_scores.tsv by similarity of prediction profiles
####################################
writeLines("\n*** Sort *.GSVA_enrichment_scores.tsv by similarity of prediction profiles ***\n")

StopWatchStart$SortGSVAERmatrix  <- Sys.time()

predictions.mat<-as.matrix(data.frame(fread(OutfileEnrichmentScores, sep="\t", na.strings=c("NA")), row.names=1))
predictions.clusters.clust <- agnes(x = predictions.mat, metric = "manhattan")
predictions.clusters.order <- rownames(predictions.mat)[predictions.clusters.clust$order]

predictions.mat.t<-t(predictions.mat)
predictions.classes.clust <- agnes(x = predictions.mat.t, metric = "manhattan")
predictions.classes.order <- rownames(predictions.mat.t)[predictions.classes.clust$order]
predictions.mat.ordered   <- predictions.mat[predictions.clusters.order,predictions.classes.order]

write.table(data.frame("ENRICHMENT"=predictions.clusters.order, predictions.mat.ordered), OutfileEnrichScorsClust, row.names = F,sep="\t",quote = F)

StopWatchEnd$SortGSVAERmatrix  <- Sys.time()

####################################
### Outfile *.GSVA_final_label.tsv - currently simply taking the maximum GSVA enrichment score for each array (column of --infile_mat)
####################################

StopWatchStart$GetFinalLabel  <- Sys.time()

writeLines("\n*** Write cluster labels based on maximum GSVA enrichment scores ***\n")

row.names(predictions.mat.ordered) <- sub("C", "", row.names(predictions.mat.ordered))

write(paste(row.names(predictions.mat.ordered), colnames(predictions.mat.ordered)[max.col(predictions.mat.ordered, ties.method="first")], sep = "\t", collapse = "\n"),
      OutfileFinalLabel)

StopWatchEnd$GetFinalLabel  <- Sys.time()

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo.txt", sep="")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt", sep="")
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
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
writeLines(paste("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report", sep = "", collapse = ""))

quit()
