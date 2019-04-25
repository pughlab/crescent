####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script that takes a genes (rows) vs. barcodes (columns) sparse matrix
### and tranforms it into a mtx set of files (barcodes.tsv, genes.tsv and matrix.mtx) or (prefix_outfiles.mtx, prefix_outfiles.mtx_rows and prefix_outfiles.mtx_cols)
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(DropletUtils)) # to habdle reading and writing mtx files
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

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

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
OutdirFinal<-paste(Outdir, "/MTX",  sep = "", collapse = "")
dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Load data
####################################

print("Loading genes (rows) vs. barcodes (columns) matrix")
input.matrix        <- as.matrix(data.frame(fread(Input),row.names=1))
input.matrix.sparse <- as(input.matrix, "dgCMatrix")
dim(input.matrix.sparse)
barcode.ids  <- colnames(input.matrix.sparse)
gene.ids     <- rownames(input.matrix.sparse)
gene.symbols <- rownames(input.matrix.sparse)

####################################
### Write output
####################################
## Be carefull with this command write10xCounts(..., overwrite=T)
## as it removes any preesixting contents of 'path' out directory
write10xCounts(path = OutdirFinal, x = input.matrix.sparse , barcodes=barcode.ids, gene.id=gene.ids,
               gene.symbol=gene.symbols, overwrite=T)

####################################
### Reformat (if needed)
####################################
if (grepl(pattern = "y", ignore.case = T, x = AddNumbers) == T) {
  barcodes <- read.table(file = paste(OutdirFinal, "/barcodes.tsv", sep = "", collapse = ""), row.names = 1)
  genes    <- read.table(file = paste(OutdirFinal, "/genes.tsv", sep = "", collapse = ""), row.names = 1)
  genesOutfile    <- paste(OutdirFinal, "/", PrefixOutfiles, ".expression_tpm.mtx_rows", sep = "", collapse = "")
  barcodesOutfile <- paste(OutdirFinal, "/", PrefixOutfiles, ".expression_tpm.mtx_cols", sep = "", collapse = "")
  #
  write(file = barcodesOutfile, x = paste(1:length(row.names(barcodes)), row.names(barcodes), sep = "\t", collapse = "\n"))
  write(file = genesOutfile,    x = paste(1:length(row.names(genes)), row.names(genes), sep = "\t", collapse = "\n"))
  #
  file.copy(from=paste(OutdirFinal, "/matrix.mtx", sep = "", collapse = ""), to=paste(OutdirFinal, "/", PrefixOutfiles, ".expression_tpm.mtx", sep = "", collapse = "") , overwrite=T)
  #
  file.remove(paste(OutdirFinal, "/barcodes.tsv", sep = "", collapse = ""))
  file.remove(paste(OutdirFinal, "/genes.tsv", sep = "", collapse = ""))
  file.remove(paste(OutdirFinal, "/matrix.mtx", sep = "", collapse = ""))
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
