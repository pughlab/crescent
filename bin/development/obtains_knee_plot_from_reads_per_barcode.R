####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script to obtain either:
### (a) read-counts from a possorted_genome_bam.bam file and then (b) a knee-plot and it's inflection, or
### entering (a), will obtain (b)
###
### Note (a) is pressumably the output from Drop-seq_tools script 'BAMTagHistogram' (http://mccarrolllab.com/dropseq/)
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse)) # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(dropbead)) # to generate the knee-plot and get inflection. Install with: library(devtools), then install_github("rajewsky-lab/dropbead")
####################################

####################################
### Required external dependencies
####################################
### Dropseq-tools 'BAMTagHistogram' if using '-t BAM'
####################################

### Change this path as needed, if using '-t BAM'
BAMTagHistogram <- "~/PROGRAMS/DROPSEQ/Drop-seq_tools-1.13/BAMTagHistogram"

####################################
### Get inputs from command line argumets
####################################

option_list <- list(
  make_option(c("-i", "--infile"), default="NA",
              help="A path/name to either a <tab> delimited *file* with number of reads per barcode, like:
                #INPUT=/path_to_file/file.bam  TAG=XC FILTER_PCR...
                203177  CTTATGGCTTTA
                196369  TTAGCGCTTATA
                151561  CGGCTGCTAATC
                Presumably produced by Dropseq-tools 'BAMTagHistogram'
                Or
                A possorted_genome_bam.bam file resumably produced by 'cellranger count'
                In the last case 'BAMTagHistogram' will be run first
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--input_type"), default="NA",
              help="Indicates if input is either a 'TAB' delimited file or a 'BAM' file
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-m", "--max_number_of_barcodes_to_plot"), default="60000",
              help="A path/name for the results directory
                Default = '60000'"),
  #   
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
 #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-c", "--cluster_name"), default="NA",
              help="Name of either of the following clusters where process are running
                'h4h'
                'local'")
)

opt <- parse_args(OptionParser(option_list=option_list))

Infile          <- opt$infile
InputType       <- opt$input_type
MaxBarcodes     <- as.numeric(opt$max_number_of_barcodes_to_plot)
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles
ClusterName     <- opt$cluster_name
Tempdir         <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile", "input_type", "outdir", "prefix_outfiles", "cluster_name")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs and get path of files
### Basically replacing any '~/' by the User HOME, because otherwise it conflicts with system(...)
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
Infile<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Infile)
BAMTagHistogram<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), BAMTagHistogram)
#
dir.create(file.path(Outdir),   showWarnings = F, recursive = T)
dir.create(file.path(Tempdir),  showWarnings = F, recursive = T)

OutPdf<-paste(Tempdir,"/", PrefixOutfiles, ".kneeplot.pdf",collapse="",sep="")
OutInf<-paste(Tempdir,"/", PrefixOutfiles, ".kneeplot_inflection.txt",collapse="",sep="")

####################################
### Obtain read-counts per cell
####################################

if(regexpr("^TAB$", InputType, ignore.case = T)[1] == 1) {
  print("Loading *tab infile")
  InfileTab <- Infile
}else if (regexpr("^BAM$", InputType, ignore.case = T)[1] == 1) {
  print("Processing possorted_genome_bam.bam infile")
  print(Sys.time())

  InfileTab <- paste(Tempdir, "/" , PrefixOutfiles, ".kneeplot_reads_per_barcode.tsv",collapse="",sep="")

  CommandsForBAMTagHistogram <- c("module load java/8",
                                  paste(BAMTagHistogram, " \\" , sep = "", collapse = ""),
                                  paste("I=", Infile, " \\", sep = "", collapse = ""),
                                  paste("O=", InfileTab, " \\", sep = "", collapse = ""),
                                  "TAG=CB")

  OutfileLocalCommandsForBAMTagHistogram<-paste(Tempdir,"/",PrefixOutfiles,".kneeplot_CommandsForBAMTagHistogram.sh", sep="")
  
  write(file = OutfileLocalCommandsForBAMTagHistogram, x = CommandsForBAMTagHistogram)
        
  Chmod<-paste("chmod +x ", OutfileLocalCommandsForBAMTagHistogram, sep = "", collapse = "")
  system(command = Chmod[[1]],  input = NULL, intern = TRUE)
  
  ### Runnig BAMTagHistogram
  ### Note in clusters we have to wait for the 'qsub' (PBS) process to finish
  ### We'll monitor this by waiting until InfileTab is written
    if (file.exists(InfileTab)) {
    file.remove(InfileTab) ## will remove preexisting InfileTab (if any)
    }

    if(regexpr("^local$", ClusterName, ignore.case = T)[1] == 1) {
      CommandToGetKneePlot<-OutfileLocalCommandsForBAMTagHistogram
      
    }else if (regexpr("^h4h$", ClusterName, ignore.case = T)[1] == 1) {
      HeadersSlurm <- "#!/bin/bash
#SBATCH -t 72:00:00  #	<-- walltime for your job in format of days-hours:minutes:seconds
#SBATCH --mem=60G    #	<-- amount of memory your job needs to run, in megabytes by default
#SBATCH -J BAMTagHist  #	<-- name of your job
#SBATCH -p himem     #  <-- partition you job wants to run
#SBATCH -c 1 	     #  <-- number of CPUs your job needs to run
#SBATCH -N 1 	     #  <-- number of nodes your job needs to run. Please use 1 node unless you are running mpi jobs.
#SBATCH -o %x-%j.out #  <-- redirect job output (both stdout and stderr) to a file called “<job name>-<job id>.out”. The Slurm default output file is called “Slurm-<job id>.out”
"
      OutfileH4HCommandsForBAMTagHistogram<-paste(Tempdir,"/",PrefixOutfiles,".kneeplot_CommandsForBAMTagHistogram.sh", sep="")
      write(file = OutfileH4HCommandsForBAMTagHistogram, x = c(HeadersSlurm, CommandsForBAMTagHistogram))
      
      #### Here need to find how to send sbatch
      
      
      SendSlurm<-paste("sbatch ", OutfileH4HCommandsForBAMTagHistogram, sep = "", collapse = "")
      # stop(SendSlurm)
      print(SendSlurm[[1]])
      system(command = SendSlurm[[1]],  input = NULL, intern = TRUE)

    }else{
      stop(paste("Unexpected name of cluster: ", ClusterName, "\n\nFor help type:\n\nRscript obtains_knee_plot_from_reads_per_barcode.R -h\n\n", sep=""))
    }
 
 print(CommandToGetKneePlot)
 system(command = CommandToGetKneePlot, input = NULL, intern = TRUE)

    ### Wait for InfileTab to be written
    while (!file.exists(InfileTab)) {
      Sys.sleep(1)
    }
  print("Done - processing possorted_genome_bam.bam infile")
  print(paste("Created file:", InfileTab), sep="", collapse="")
  print(Sys.time())

}else{
 stop(paste("Unexpected type of infile: ", InputType, "\n\nFor help type:\n\nRscript obtains_knee_plot_from_reads_per_barcode.R -h\n\n", sep=""))
}

####################################
### Load read-counts per cell data and get knee-plot
### Note: earlier versions of this script used library(dropbead) plotCumulativeFractionOfReads() and estimateCellNumber()
###       but they failed to produce a kneeplot for dataset: ~/SINGLE_CELL/10X/SMARTER_VACCINE/DATA_FROM_PUGHLAB/SMTR04t1_NonRad/5p_RNA
###       Thus decided to program this part over
####################################
print(paste("Start reading table of read-counts from:",InfileTab))
print(Sys.time())

mat<-read.table(InfileTab, row.names=NULL, header = FALSE)
colnames(mat)<-c("reads","barcode")

### Get cumulative

barcode_number <- 0
for (i in rownames(mat)) {
  barcode_number <- barcode_number+1
  ivalue <- mat[barcode_number,"reads"]
  if (barcode_number == 1) {
    jvalue <- ivalue
  }else{
    jvalue <- jvalue + ivalue
  }
  mat[barcode_number,"cumulative"] <- jvalue
}

####################################
### Generate outfiles
####################################

### Kneeplot
XaxisLenght <- min(MaxBarcodes,nrow(mat))
pdf(OutPdf)
InflectionPoint<-estimateCellNumber(mat[, "reads"], max.cells = min(XaxisLenght))
plot(x = c(1:XaxisLenght), y = mat[c(1:XaxisLenght),"cumulative"], type="l",
     xlab = "Sorted barcodes", ylab = "Cumulative reads", lwd = 2, main = paste(PrefixOutfiles, "\n", "Inflection = ", InflectionPoint, " cells", sep = "", collapse = ""))
abline(v=InflectionPoint, col = "red", lwd = 2, lty =2)
dev.off()

### Kneeplot inflection
write(file = OutInf, x = InflectionPoint)

### Top barcodes based on kneeplot inflection
OutTop<-paste(Tempdir, "/", PrefixOutfiles, "_Top", InflectionPoint, "cells.csv", sep="", collapse = "")
write.table(file = OutTop, x = "Barcode", row.names = F, col.names = F, quote = F)
write.table(file = OutTop, x = mat[c(1:InflectionPoint),"barcode"], row.names = F, col.names = F, quote = F, append = T)

####################################
### Report used options
####################################
OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".kneeplot_UsedOptions.txt", sep="")
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

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".kneeplot_CPUusage.txt", sep="")
ReportTime<-c(
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir or keeping them at tempdir (if using CWL)
####################################

writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir, sep="", collapse = ""))
outfiles_to_move <- list.files(Tempdir, pattern = PrefixOutfiles, full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile) {
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""), to=paste(Outdir, "/", eachFile, sep=""), overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

####################################
### Finish
####################################

print("END - All done!!!")
quit()
