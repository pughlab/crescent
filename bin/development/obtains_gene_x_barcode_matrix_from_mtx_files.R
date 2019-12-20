####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script that takes a mtx set of files (barcodes.tsv, genes.tsv and matrix.mtx)
### and transforms it into a genes (rows) vs. barcodes (columns) sparse matrix
###
### NEW IMPLEMENTATIONS SINCE Seurat v2:
### 1) Rewritten with Seurat v3 commands (including ability to read output from Cell Ranger v3)
###    Main differences vs. Seurat v2 include:
###    a) new function names
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
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
              help="The path/name to a 10X *directory* with barcodes.tsv, genes.tsv and matrix.mtx files"),
#
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
#
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID")
)
opt <- parse_args(OptionParser(option_list=option_list))

Input          <- opt$input
Outdir         <- opt$outdir
PrefixOutfiles <- opt$prefix_outfiles

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

StartTimeOverall<-Sys.time()

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
### Create outdirs
####################################
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
Outdir<-gsub("/$", "", Outdir)
Tempdir<-gsub("/$", "", Tempdir)
#
dir.create(file.path(Outdir, "GENE_VS_BARCODE_MATRIX"), recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")
input.matrix <- Read10X(data.dir = Input)
dim(input.matrix)

####################################
### Create a Seurat object
####################################
writeLines("\n*** Create a Seurat object ***\n")
seurat.object  <- CreateSeuratObject(counts = input.matrix, project = PrefixOutfiles)

####################################
### Write gene vs. barcode outfile
####################################
OutfileDGE<-paste(Tempdir,"/",PrefixOutfiles,".gene_vs_barcode.tsv", sep="")
Headers<-paste("GENES",paste(colnames(GetAssayData(object = seurat.object, slot = 'counts')),sep="",collapse = "\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileDGE, row.names = F, col.names = F, sep="\t", quote = F)
write.table(x=as.matrix(GetAssayData(object = seurat.object, slot = 'counts')), file = OutfileDGE, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Report used options
####################################
OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".gene_vs_barcode_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report time used
####################################
EndTimeOverall<-Sys.time()

TookTimeOverall <-format(difftime(EndTimeOverall, StartTimeOverall, units = "min"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".gene_vs_barcode_CPUusage.tsv", sep="")
ReportTime<-c(
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".gene_vs_barcode", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/GENE_VS_BARCODE_MATRIX/",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################

print("END - All done!!! Took time:")
print(ReportTime)

quit()

