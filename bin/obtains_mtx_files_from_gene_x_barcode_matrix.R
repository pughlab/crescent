####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script that takes a genes (rows) vs. barcodes (columns) sparse matrix
### and tranforms it into a mtx set of files (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz)
### or (prefix_outfiles.mtx, prefix_outfiles.mtx_rows and prefix_outfiles.mtx_cols)
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(earlycross))   # to handle reading and writing mtx files
### Which can be installed like:
### install.packages('devtools')
### devtools::install_github("daskelly/earlycross")
suppressPackageStartupMessages(library(Seurat))       # to create a Seurat object for earlycross
### Requires Seurat v3 (tested on v3.0.3.9023), which can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table
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
              help="Path/name to a genes (rows) vs. barcodes (columns) matrix"),
  #
  make_option(c("-o", "--outdir"), default="selected_gene_bc_matrices",
              help="A path/name for the directory where the new mtx files will be saved"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
              Note: if using option `-l y`, outfile names will be:
              prefix_outfiles.mtx, prefix_outfiles.mtx_rows and prefix_outfiles.mtx_cols
              Or if using `-l n`, outfile names will be:
              'barcodes.tsv', 'genes.tsv' and 'matrix.mtx'"),
  #
  make_option(c("-l", "--add_barcode_and_gene_numbers"), default="N",
              help="Indicates if 'barcodes.tsv' and 'genes.tsv' outfiles should have numbers in the first column (type [y/Y] or [n/N]), like:
              1	SRR3541565
              2	SRR3541564
              3	SRR3541563
              and
              1	ENSG00000000003
              2	ENSG00000000005
              3	ENSG00000000419")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input          <- opt$input
Outdir         <- opt$outdir
PrefixOutfiles <- opt$prefix_outfiles
AddNumbers     <- opt$add_barcode_and_gene_numbers

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Load data
####################################

print("Loading genes (rows) vs. barcodes (columns) matrix")
input.matrix <- data.frame(fread(Input),row.names=1)

print("Creating Seurat object")
seurat.object  <- CreateSeuratObject(counts = input.matrix, project = PrefixOutfiles)

OutdirFinal<-paste(Outdir, "/MTX",  sep = "", collapse = "")
dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)

### Remove preexisting files
MtxFilesList <- list("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
for (file in MtxFilesList) {
  if (file.exists(paste(OutdirFinal, "/", file, sep = "", collapse = ""))) {
  system(command = paste("rm ", OutdirFinal, "/", file, sep = "", collapse = ""))
  }
}

print("Writing MTX files")
Write10X(obj = seurat.object, dir = OutdirFinal)

####################################
### Reformat (if needed)
####################################
if (grepl(pattern = "y", ignore.case = T, x = AddNumbers) == T) {
  barcodes <- read.table(file = paste(OutdirFinal, "/barcodes.tsv.gz", sep = "", collapse = ""), row.names = 1)
  genes    <- read.table(file = paste(OutdirFinal, "/features.tsv.gz", sep = "", collapse = ""), row.names = 2)
  genesOutfile    <- paste(OutdirFinal, "/", PrefixOutfiles, ".expression_tpm.mtx_rows.gz", sep = "", collapse = "")
  barcodesOutfile <- paste(OutdirFinal, "/", PrefixOutfiles, ".expression_tpm.mtx_cols.gz", sep = "", collapse = "")
  #
  write(file = barcodesOutfile, x = paste(1:length(row.names(barcodes)), row.names(barcodes), sep = "\t", collapse = "\n"))
  write(file = genesOutfile,    x = paste(1:length(row.names(genes)), row.names(genes), sep = "\t", collapse = "\n"))
  #
  file.copy(from=paste(OutdirFinal, "/matrix.mtx.gz", sep = "", collapse = ""), to=paste(OutdirFinal, "/", PrefixOutfiles, ".expression_tpm.mtx.gz", sep = "", collapse = "") , overwrite=T)
  #
  file.remove(paste(OutdirFinal, "/barcodes.tsv.gz", sep = "", collapse = ""))
  file.remove(paste(OutdirFinal, "/features.tsv.gz", sep = "", collapse = ""))
  file.remove(paste(OutdirFinal, "/matrix.mtx.gz",   sep = "", collapse = ""))
}

####################################
### Report used options
####################################
OutfileOptionsUsed<-paste(Outdir,"/",PrefixOutfiles,".mtx_from_gene_vs_barcode_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report time used
####################################
EndTimeOverall<-Sys.time()

TookTimeOverall <-format(difftime(EndTimeOverall, StartTimeOverall, units = "secs"))
OutfileCPUusage<-paste(Outdir,"/",PrefixOutfiles,".mtx_from_gene_vs_barcode_CPUusage.txt", sep="")
ReportTime<-c(
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################

print("END - All done!!!")

quit()
