####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script that takes a mtx set of files (barcodes.tsv, genes.tsv and matrix.mtx)
### and transforms it into a genes (rows) vs. barcodes (columns) sparse matrix
####################################

####################################
### Required libraries
####################################
### 'optparse'   to handle one-line-commands
### 'Seurat'     to run QC, differential gene expression and clustering analyses
### 'dplyr'      needed by Seurat for data manupulation
####################################

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

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
Tempdir<-gsub("/$", "", Outdir)
#
dir.create(file.path(Outdir, "DGE"), recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Load data
####################################
print("Loading 10X infiles")
input.matrix <- Read10X(data.dir = Input)
dim(input.matrix)

####################################
### Create a Seurat object
####################################

seurat.object  <- CreateSeuratObject(raw.data = input.matrix)
seurat.object

####################################
### Write DGE out
####################################
OutfileDGE<-paste(Tempdir,"/",PrefixOutfiles,".DGE.tsv", sep="")

Headers<-paste("GENES",paste(colnames(seurat.object@data),sep="",collapse = "\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileDGE, row.names = F, col.names = F, sep="\t", quote = F)
write.table(x=seurat.object@data, file = OutfileDGE, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Report used options
####################################
OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".DGE_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report time used
####################################
EndTimeOverall<-Sys.time()

TookTimeOverall        <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "min"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".DGE_CPUusage.tsv", sep="")
ReportTime<-c(
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".DGE_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){ 
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/DGE/",eachFile,sep=""),overwrite=T)
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

