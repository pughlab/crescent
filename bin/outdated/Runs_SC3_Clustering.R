####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html
### Things missing from this tutorial:
### 1) Automated selection of granularity of clusters
### 2) Assigning cell type identity to clusters (needs supervised annotations, maybe based on Gene Set Enrichment analysis)
### 3) Iterative shiny sc3_interactive() application, because it needs a web navigator
### 4) Use knitr() to produce better html plot layout (https://yihui.name/knitr/demo/stitch/)
### 5) See note on SC3 developers define spike-in controls at the end of this script
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))               # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))             # to speed up reading matrices with 'fread' instead of 'read.table'
suppressPackageStartupMessages(library(SingleCellExperiment))   # (bioconductor) to create SingleCellExperiment objects
suppressPackageStartupMessages(library(SC3))                    # (bioconductor) to do the SC3 clustering
suppressPackageStartupMessages(library(scater))                 # (bioconductor) to generate scatter plots and calculate single-cell QC metrics
####################################

####################################
### Get inputs from command line argumets
####################################

option_list <- list(
  make_option(c("-i", "--infile_mat"), default="NA",
              help="A path/name to a <tab> delimited *file* with cell barcodes in columns and genes in rows"),

  make_option(c("-t", "--infile_cell_types"), default="NA",
              help="A path/name to a <tab> delimited *file* with two columns, like:
                barcode             cell_type
                AAACCTGCACCGAAAG.1  cell_type_label_1
                AAACCTGTCCTCATTA.1  cell_type_label_1
                GGACGGGGTACTTAGC.1  cell_type_label_2
                AAAGATGGTTGCTCCT.1  cell_type_label_3

                Note: if not applicable just type 'NA'"),

  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID"),
  
  make_option(c("-l", "--min_ks"), default="2",
              help="Minimum number of clusters k used for SC3 clustering, e.g. '2'"),
  
  make_option(c("-u", "--max_ks"), default="10",
              help="Maximum number of clusters k used for SC3 clustering, e.g. '10'"),
  
  make_option(c("-c", "--max_cores"), default="2",
              help="Maximum number of computer cores to use for SC3 clustering, e.g. '2'")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat       <- opt$infile_mat
InfileCellTypes <- opt$infile_cell_types
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles
MinKs           <- opt$min_ks
MaxKs           <- opt$max_ks
MaxCores        <- opt$max_cores
Tempdir         <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

####################################
### Define tailored parameters
####################################
### Some of these parameters are the defaults provided by SC3 developers
###
### Parameters for cell spike-ins
SpikeInLabels  <- "ERCC-" ### It's important to use '-' to distingish spike-ins from Excision Repair Cross-Complementing proteins

### Parameters for data normalization
AddPseudocount  <- 1e-99 ### Default is 1, but using this value to match what was used in 'Run_Seurat_Clustering.R' script. Also see https://goo.gl/3VzQ3L

StartTimeOverall <-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_mat", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs
####################################
writeLines("\n*** Create outdirs ***\n")
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
Outdir<-gsub("/$", "", Outdir)
Tempdir<-gsub("/$", "", Tempdir)
#
dir.create(file.path(Outdir, "SC3"), showWarnings = F, recursive = T)
dir.create(file.path(Tempdir),        showWarnings = F, recursive = T)

####################################
### Load data
####################################

fullmat<-data.frame(fread(InfileMat, sep="\t", na.strings=c("NA")), row.names=1)
if (InfileCellTypes == "NA") {
## Do nothing
  }else{
  celltypes<-data.frame(fread(InfileCellTypes, sep="\t", na.strings=c("NA")), row.names=1)
  CellTypeHeader<-colnames(celltypes)[1]
}

####################################
### Create a SingleCellExperiment object
####################################

if (InfileCellTypes == "NA") {
  singlecellexp.object <- SingleCellExperiment(assays = list(counts = as.matrix(fullmat), logcounts = log2(as.matrix(fullmat) + AddPseudocount)))
}else{
  singlecellexp.object <- SingleCellExperiment(assays = list(counts = as.matrix(fullmat), logcounts = log2(as.matrix(fullmat) + AddPseudocount)), colData = celltypes)
}
singlecellexp.object

####################################
### Define feature names (gene_IDs or gene names) in feature_symbol column
####################################

rowData(singlecellexp.object)$feature_symbol <- rownames(singlecellexp.object)
### Remove features with duplicated names
### This shoulndn't affect the singlecellexp.object if we are using unique Gene_IDs
singlecellexp.object <- singlecellexp.object[!duplicated(rowData(singlecellexp.object)$feature_symbol), ]

### Define spike-ins
isSpike(singlecellexp.object, SpikeInLabels) <- grepl(SpikeInLabels, rowData(singlecellexp.object)$feature_symbol)

####################################
### PCA plot
####################################

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SC3_plotPCA.pdf", sep=""),width=7, height=7)
if (InfileCellTypes == "NA") {
  plotPCA(singlecellexp.object)
}else{
  plotPCA(singlecellexp.object, colour_by = CellTypeHeader)
}
dev.off()

####################################
### Clustering
####################################
### See note on CPU/time usage note if including 'biology' step at the end of this script
####################################

StartTimeClustering<-Sys.time()
singlecellexp.object <- sc3(singlecellexp.object, ks = MinKs:MaxKs, biology = FALSE, n_cores=MaxCores)
EndTimeClustering<-Sys.time()

StartTimeBiol<-Sys.time()
singlecellexp.object <- sc3_calc_biology(object = singlecellexp.object, ks=MinKs:MaxKs)
EndTimeBiol<-Sys.time()

### Clustering results are added to the colData slot and can be found by "sc3_" prefix
### This is the same table exported by command 'sc3_export_results_xls' into an Excel file
clusters_data<-(colData(singlecellexp.object)[ , grep("sc3_", colnames(colData(singlecellexp.object)))])

OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles,".SC3_clusters.tsv", sep="")
write.table(data.frame("CLUSTERS"=rownames(clusters_data),clusters_data),file = OutfileClusters, row.names = F, sep="\t", quote = F)

####################################
### Finding differentially expressed genes (cluster biomarkers)
####################################

MarkersData <- rowData(singlecellexp.object)
OutfileMarkersPerCluster<-paste(Tempdir,"/",PrefixOutfiles,".SC3_MarkersPerCluster.tsv", sep="")
write.table(MarkersData, file = OutfileMarkersPerCluster, row.names = F, col.names = T, sep="\t", quote = F)

####################################
###  PCA plot coloured by clustering results
####################################

### Having SC3 results colour and size of dots can be used for a scatter plot
ColourBy<-paste("sc3_",MaxKs,"_clusters",sep = "")
SizeBy<-paste("sc3_",MaxKs,"_log2_outlier_score",sep = "")

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SC3_plotPCA_ColourAndSizeByCluster.pdf", sep=""),width=7, height=7)
plotPCA(
  singlecellexp.object,
  colour_by = ColourBy,
  size_by = SizeBy
)
dev.off()

###################################
## Report used options
###################################
OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".SC3_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report time used
####################################
EndTimeOverall<-Sys.time()

TookTimeClustering <-format(difftime(EndTimeClustering,StartTimeClustering, units = "min"))
TookTimeBiology    <-format(difftime(EndTimeBiol,      StartTimeBiol,       units = "min"))
TookTimeOverall    <-format(difftime(EndTimeOverall,   StartTimeOverall,    units = "min"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".SC3_CPUusage.tsv", sep="")
ReportTime<-c(
  paste("clustering",TookTimeClustering,collapse = "\t"),
  paste("biology",TookTimeBiology,collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".SC3_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/SC3/",eachFile,sep=""),overwrite=T)
  file.remove(from=paste(Tempdir,"/",eachFile,sep=""))
})


####################################
### Finish
####################################

print("END - All done!!! Took time:")
print(ReportTime)

quit()

#######################################################
#######################################################
### CPU/time usage note if including 'biology' step ###
#######################################################
### Using a matrix of 33694 genes x 690 cells with a MacBookPro14,2 (Intel Core i5, 3.1 GHz, 2 cores, 16 GB RAM)
###
### Clustering alone by:
### singlecellexp.object <- sc3(singlecellexp.object, ks = MinKs:MaxKs, biology = FALSE, n_cores=MaxCores)
### Takes ~3.91 min
### Then running the sc3_calc_biology() alone by:
### sc3_calc_biology(object = singlecellexp.object, ks=MinKs:MaxKs)
### Takes ~1.59 min. Thus the two steps are ~5.5 min
###
### In contrast, running clustering and biology together by:
### singlecellexp.object <- sc3(singlecellexp.object, ks = MinKs:MaxKs, biology = TRUE, n_cores=MaxCores)
### Takes ~4.99 min
#######################################################
#######################################################

#######################################################
#######################################################
### Spike-in definition by SC3 developers #############
#######################################################
### SC3 tutorial:
### https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html
### uses data 'ann' and the following command to get spike-ins
### `isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)`
### 
### ERCC stands for External RNA Controls Consortium, which are the spike-ins and their labels are like:
### ERCC-00074, ERCC-00004, ERCC-00113, etc.
### 
### But ERCC also stands for Excision repair cross-complementing, which are a set of proteins involved in DNA repair.
### In humans, ERCC proteins are transcribed from the following genes: ERCC1, ERCC2, ERCC3, ERCC4, ERCC5, ERCC6, and ERCC8
### The labels produced with the command above using the 'ann' data from the SC3 tutorial are the last one, not the spike-ins
### Thus the assumption that they are spike-ins may be incorrect. Better to use:
### `isSpike(sce, "ERCC-") <- grepl("ERCC-", rowData(sce)$feature_symbol)`
#######################################################
#######################################################
