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
ProgramOutdir  <- "GSVA"

####################################
### Get inputs from command line argumets
####################################
#####
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
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicates one of three options:
                a) Type the number of cores to use for parellelization (e.g. '4')
                b) Type 'MAX' to determine and use all available cores in the system
                
                Default = 'MAX'"),
  #
  make_option(c("-w", "--run_cwl"), default="0",
              help="Indicates if this script should produce 'frontend' files for crescent.cloud
                0 = no frontend files should be produced
                1 = frontend files should be produced and '--minio_path path' is provided
                2 = frontend files should be produced but '--minio_path path' is not provided (i.e local run)

                Default = '0'")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat          <- opt$infile_mat
InfileMatType      <- opt$infile_mat_type
InfileGmt          <- opt$infile_gmt
Outdir             <- opt$outdir
PrefixOutfiles     <- opt$prefix_outfiles
PvalueCutoff       <- as.numeric(opt$pvalue_cutoff)
FdrCutoff          <- as.numeric(opt$fdr_cutoff)
RunsCwl            <- as.numeric(opt$run_cwl)
MaxGlobalVariables <- as.numeric(opt$max_global_variables)
NumbCores          <- opt$number_cores

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
### Define default parameters
####################################

DefaultParameters <- list(
  ### General
  DigitsForRound = 5, ### Number of digits to round up enrichment score values in outfiles
  
  ### Scoring
  MinEScoreToAssignLabel = 0,
  LabelForUndefinedScored = "UNDEFINED"
  
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

FILE_TYPE_OUT_DIRECTORIES = c(
  "GSVA",
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

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_UsedOptions.txt")

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

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_LogFiles_", "RSessionInfo", ".txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

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
EnrichmentScores<-gsva(expr=fullmat, gset.idx.list=gmt2, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=T, parallel.sz=NumbCoresToUse)
SortedRowNames<-rownames(EnrichmentScores)
SortedColNames<-colnames(EnrichmentScores)
#
EnrichmentScores<-EnrichmentScores[SortedRowNames,SortedColNames]
OutfileEnrichmentScores<-paste0(Tempdir, "/GSVA/", PrefixOutfiles, ".GSVA_enrichment_scores.tsv")
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
OutfileAllScores<-paste0(Tempdir, "/GSVA/", PrefixOutfiles, ".GSVA_all_scores_table.tsv")
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
dataEPQ <- read.table(file = OutfileAllScores,row.names = NULL ,header = T, sep = "\t")
#
ForPvaluesMat<- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"p.Val"])
PvaluesMat<-xtabs(z~x+y, data=ForPvaluesMat)
PvaluesMat<-PvaluesMat[SortedRowNames,SortedColNames]
HeadersForPvalues<-paste(c("PVALUES", colnames(PvaluesMat)), sep="\t", collapse = "\t")
OutfilePvalues<-paste0(Tempdir, "/GSVA/", PrefixOutfiles, ".GSVA_pvalues.tsv")
write(x=HeadersForPvalues,file=OutfilePvalues)
write.table(x=PvaluesMat, file=OutfilePvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
#
ForFDRvaluesMat <- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"FDR"])
FdrvaluesMat<-xtabs(z~x+y, data=ForFDRvaluesMat)
FdrvaluesMat<-FdrvaluesMat[SortedRowNames,SortedColNames]
HeadersForFdrvalues<-paste(c("FDR_VALUES", colnames(FdrvaluesMat)), sep="\t", collapse = "\t")
OutfileFdrvalues <-paste0(Tempdir, "/GSVA/", PrefixOutfiles, ".GSVA_fdr_values.tsv")
write(x=HeadersForFdrvalues,file=OutfileFdrvalues)
write.table(x=FdrvaluesMat, file=OutfileFdrvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
#
FilteredESMatLogical<-(PvaluesMat<=PvalueCutoff & FdrvaluesMat<=FdrCutoff)
FilteredESMatLogical<-FilteredESMatLogical[SortedRowNames,SortedColNames]
FilteredESMatValues<-ifelse(FilteredESMatLogical==TRUE,EnrichmentScores,NA)
OutfileFilteredES <-paste0(Tempdir, "/GSVA/", PrefixOutfiles, ".GSVA_filtered.tsv")
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

OutfileEnrichScorsClust <-paste0(Tempdir, "/GSVA/", PrefixOutfiles, ".GSVA_enrichment_scores_sorted.tsv")
write.table(data.frame("ENRICHMENT"=predictions.clusters.order, predictions.mat.ordered), OutfileEnrichScorsClust, row.names = F,sep="\t",quote = F)

StopWatchEnd$SortGSVAERmatrix  <- Sys.time()

####################################
### Outfile *.GSVA_final_label.tsv - currently simply taking the maximum GSVA enrichment score for each array (column of --infile_mat)
### and making enrichment scores <= 0 to show 'UNDEFINED' predictions
####################################

StopWatchStart$GetFinalLabel  <- Sys.time()

writeLines("\n*** Write cluster labels based on maximum GSVA enrichment scores ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  row.names(predictions.mat.ordered) <- sub("C", "", row.names(predictions.mat.ordered))
}

HighestScoreLabels <- colnames(predictions.mat.ordered)[max.col(predictions.mat.ordered, ties.method="first")]
HighestScoreLabels.df <- as.data.frame(cbind(row.names(predictions.mat.ordered), HighestScoreLabels))
colnames(HighestScoreLabels.df)[1] <- "cluster"

# just to clean-up OutfileFinalLabel if there is a pre-existing run
OutfileClusterFinalLabel<-paste0(Tempdir, "/GSVA/", PrefixOutfiles, ".GSVA_cluster_final_label.tsv")
close(file(OutfileClusterFinalLabel, open="w"))

sapply(row.names(HighestScoreLabels.df), FUN=function(eachClusterN) {
  ClusterId     <- as.character(HighestScoreLabels.df[eachClusterN,"cluster"])
  LabelAssigned <- as.character(HighestScoreLabels.df[eachClusterN,"HighestScoreLabels"])
  EnrichmentScore <- predictions.mat.ordered[ClusterId,LabelAssigned]

  if (EnrichmentScore > DefaultParameters$MinEScoreToAssignLabel) {
    write(file = OutfileClusterFinalLabel, x = paste(ClusterId, LabelAssigned, paste0(ClusterId, ".", LabelAssigned), sep = "\t"), append = T)
  }else{
    write(file = OutfileClusterFinalLabel, x = paste(ClusterId, DefaultParameters$LabelForUndefinedScored, paste0(ClusterId, ".", DefaultParameters$LabelForUndefinedScored), sep = "\t"), append = T)
  }
})

StopWatchEnd$GetFinalLabel  <- Sys.time()

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

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

writeLines("\n*** Moving outfiles into outdir or keeping them at tempdir (if using CWL) ***\n")

### using two steps to copy files (`file.copy` and `file.remove`) instead of just `file.rename` to avoid issues with path to Tempdir in cluster systems
if (RunsCwl == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(Tempdir)
}else{
  writeLines("\n*** Moving outfiles into outdir ***\n")
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    OutdirFinal <- paste0(Outdir, "/", ProgramOutdir, "/", DirName)
    print(OutdirFinal)
    dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
    sapply(list.files(TempdirWithData, pattern = paste0("^", PrefixOutfiles, ".", ProgramOutdir), full.names = F), FUN=function(EachFileName) {
      file.copy(from=paste0(TempdirWithData, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
      file.remove(paste0(TempdirWithData, "/", EachFileName))
    })
    sapply(list.files(TempdirWithData, pattern = paste0("^", "infercnv"), full.names = F), FUN=function(EachFileName) {
      file.copy(from=paste0(TempdirWithData, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
      file.remove(paste0(TempdirWithData, "/", EachFileName))
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
