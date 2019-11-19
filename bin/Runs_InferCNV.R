####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://github.com/broadinstitute/inferCNV/wiki#quick_start
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))       # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))     # to read tables quicker than read.table
suppressPackageStartupMessages(library(future))         # To run parallel processes
suppressPackageStartupMessages(library(infercnv))       # (Bioconductor) to create 'infercnv_obj'
### if (!requireNamespace("BiocManager", quietly = TRUE))
###  install.packages("BiocManager")
### BiocManager::install("infercnv")
suppressPackageStartupMessages(library(Seurat))         # to load MTX format scRNA-seq data
### Requires Seurat v3 (tested on v3.0.3.9023), which can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")

####################################
### Required external packages
####################################
### 'jags'    JAGS (Just Another Gibbs Sampler)
###           https://sourceforge.net/projects/mcmc-jags/
### 'ShaidyMapGen.jar' For rendering NGCHMs. Required if using NGCHM R package to generate standalone NGCHMs
###           https://www.ngchm.net/Downloads/index.html
PathToShaidyMapGenJar <- "~/PROGRAMS/SHAIDYMAPGEN/ShaidyMapGen.jar"
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
PathToShaidyMapGenJar<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), PathToShaidyMapGenJar)
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
  make_option(c("-a", "--input_annotations"), default="NA",
              help="Path/name of a table with annotations file which indicates which cells are tumor vs. normal
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-g", "--input_gene_coordinates"), default="NA",
              help="Path/name of a table with gene/chromosome positions
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-m", "--gene_mean_count_cutoff"), default="0.1",
              help="Indicates which genes will be used for the infercnv analysis.
                For smart-seq a value of '1' works well. For 10x, where the count matrix tends to be more sparse, a value of '0.1' work well
                Default = '0.1'"),
  #
  make_option(c("-n", "--noise_filter"), default="0.1",
              help="Values +- from the reference cell mean will be set to zero (whitening effect) default
                'NA', instead will use sd_amplifier
                Default = '0.1'"),
  #
  make_option(c("-s", "--sd_amplifier"), default="0.15",
              help="Used only if --noise_filter is set to 'NA'
                Noise is defined as mean(reference_cells) +- sdev(reference_cells) * sd_amplifier
                Default = '0.15'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),
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
InputAnnotations        <- opt$input_annotations
InputCoordinates        <- opt$input_gene_coordinates
GeneMeanCountCutoff     <- as.numeric(opt$gene_mean_count_cutoff)
NoiseFilter             <- as.numeric(opt$noise_filter)
SdAmplifier             <- as.numeric(opt$sd_amplifier)
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
NumbCores               <- opt$number_cores
RunsCwl                 <- opt$run_cwl

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

ProgramOutdir <- "INFERCNV"

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
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript Runs_InferCNV.R -h\n\n", sep=""))
}

cat("Using ", NumbCoresToUse, "cores")

plan(strategy = "multicore", workers = NumbCoresToUse)

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4Gb of RAM (4000*1024^2), use:
options(future.globals.maxSize = 4000 * 1024^2)

####################################
### Define default parameters
####################################

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
  counts.df <- as.data.frame(Read10X(data.dir = InputCounts))
}else if (regexpr("^DGE$", InputCountsType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  counts.df <- data.frame(fread(InputCounts),row.names=1, check.names = FALSE)
}else{
  stop(paste("Unexpected type of --input_counts: ", InputCountsType, "\n\nFor help type:\n\nRscript Runs_InferCNV.R -h\n\n", sep=""))
}

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Load annotations and gene coordinates data
####################################
writeLines("\n*** Load annotations and gene coordinates data ***\n")

StopWatchStart$CreateInfercnvObject  <- Sys.time()

annots.df<-as.data.frame(read.table(file = InputAnnotations, header = F, row.names = 1, sep = "\t", check.names = FALSE))
coords.df<-as.data.frame(read.table(file = InputCoordinates, header = F, row.names = 1, sep = "\t", check.names = FALSE))

### needed because Seurat's Read10X() replaces '-' by '.' in barcode IDs
if (regexpr("^MTX$", InputCountsType, ignore.case = T)[1] == 1) {
  rownames(annots.df)<-gsub(pattern = "-", replacement = ".", x = rownames(annots.df))
}

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts.df,
                                    annotations_file=annots.df,
                                    gene_order_file=coords.df,
                                    ref_group_names=c("Reference"
                                    ))


StopWatchEnd$CreateInfercnvObject  <- Sys.time()

####################################
### Run infercnv
####################################
writeLines("\n*** Run infercnv ***\n")

StopWatchStart$InfercnvRun  <- Sys.time()

### Note using paths such as "~/temp" or "/Users/UserName/temp" for out_dir produces:
### `(Error in nls(y ~ .logistic_midpt_slope(x, midpt = x0, slope = k), data = df, : 
### number of iterations exceeded maximum of 50), couldn't fit logistic, but no worries, going to use a spline
### (Error in nls(y ~ .logistic_midpt_slope(x, midpt = x0, slope = k), data = df, : 
### number of iterations exceeded maximum of 50), couldn't fit logistic, but no worries, going to use a spline
### Error in if (runif(1) <= padj) { : missing value where TRUE/FALSE needed`
### Thus, we stick to `tempfile()``

if (length(grep('^NA$',opt[["noise_filter"]], perl = T))) {
  HeatmapHeader<-paste(" gene_cutoff=", GeneMeanCountCutoff, " noise_filter=NA sd_amplifier=", as.numeric(opt$sd_amplifier), " ", sep = "", collapse = "")
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=GeneMeanCountCutoff,
                               out_dir=tempfile(), 
                               cluster_by_groups=TRUE,
                               plot_steps = TRUE,
                               scale_data = TRUE,
                               denoise=TRUE,
                               sd_amplifier=as.numeric(opt$sd_amplifier),
                               HMM = TRUE,
                               HMM_type = "i6"
  )

}else{
  HeatmapHeader<-paste("gene_cutoff=", GeneMeanCountCutoff, " noise_filter=", as.numeric(opt$noise_filter), " sd_amplifier=NA", sep = "", collapse = "")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=GeneMeanCountCutoff,
                               out_dir=tempfile(), 
                               cluster_by_groups=TRUE,
                               plot_steps = TRUE,
                               scale_data = TRUE,
                               denoise=TRUE,
                               noise_filter=as.numeric(opt$noise_filter),
                               HMM = TRUE,
                               HMM_type = "i6"
  )
}

StopWatchEnd$InfercnvRun  <- Sys.time()

####################################
### Generate matrix and heatmap
####################################
writeLines("\n*** Generate matrix and heatmap ***\n")

StopWatchStart$PlotCnv  <- Sys.time()

plot_cnv(infercnv_obj = infercnv_obj, out_dir = Tempdir, title = HeatmapHeader, output_format = "pdf" )

StopWatchEnd$PlotCnv  <- Sys.time()

####################################
### Transform *txt files into *tsv
####################################
writeLines("\n*** Transform *txt files into *tsv ***\n")

StopWatchStart$TransformTxtToTsv  <- Sys.time()

ListOfTxtFilesToTsv<- c(paste(Tempdir, "/infercnv.references.txt", sep = "", collapse = ""),
                        paste(Tempdir, "/infercnv.observations.txt", sep = "", collapse = ""),
                        paste(Tempdir, "/infercnv.observation_groupings.txt", sep = "", collapse = "")
)

sapply(ListOfTxtFilesToTsv,FUN=function(eachFile) {
  intxt  <- eachFile
  outtsv <- gsub(".txt$",".tsv", x= intxt, ignore.case = T, perl = T)
  in.mat <- data.frame(fread(intxt),row.names=1)
  Headers<-paste("DATA", paste(names(in.mat), sep = "", collapse = "\t"), sep="\t", collapse = "\t")
  write.table(Headers,file = outtsv, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(in.mat, file = outtsv, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  file.remove(intxt)
})

StopWatchEnd$TransformTxtToTsv  <- Sys.time()

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
  
  outfiles_to_move <- list.files(Tempdir, pattern = "infercnv", full.names = F)
  sapply(outfiles_to_move,FUN=function(eachFile) {
    ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
    inFile<-eachFile
    outFile<-gsub("infercnv", PrefixOutfiles, inFile)
    inFile<-paste(Tempdir, "/", inFile, sep = "", collapse = "")
    outFile<-paste(Outdir,"/", ProgramOutdir, "/", outFile, sep = "", collapse = "")
    file.copy(from=inFile, to=outFile, overwrite=T)
    file.remove(inFile)
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