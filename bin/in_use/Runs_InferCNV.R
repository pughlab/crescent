####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(optparse))       # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))     # to read tables quicker than read.table
suppressPackageStartupMessages(library(future))         # To run parallel processes
suppressPackageStartupMessages(library(infercnv))       # (Bioconductor) to create 'infercnv_obj'
suppressPackageStartupMessages(library(Seurat))         # to load MTX format scRNA-seq data

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_InferCNV.R"
ProgramOutdir  <- "INFERCNV"

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
  make_option(c("-j", "--input_annotations"), default="NA",
              help="Path/name of a <tab> delimited file with cell type annotations, like:
                Barcode1  Brain_Cortex
                Barcode2  Brain_Cerebellum
                Barcode3  MGH31_tumor
                Barcode4  MGH32_tumor

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-k", "--normal_cell_types"), default="NA",
              help="A <comma> delimited list of cell types from --input_annotations to use as normal cells, like:
                Brain_Cortex,Brain_Cerebellum
              
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-g", "--input_gene_coordinates"), default="NA",
              help="Path/name of a <tab> delimited file with gene, chromosome and coordinates, like:
                DDX11L1  chr1  11869  14412
                WASH7P   chr1  14363  29806
                OR4G4P   chr1  52473  54936

                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-m", "--cutoff"), default="0.1",
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
  make_option(c("-w", "--run_cwl"), default="0",
              help="Indicates if this script should produce 'frontend' files for crescent.cloud
                0 = no frontend files should be produced
                1 = frontend files should be produced and '--minio_path path' is provided
                2 = frontend files should be produced but '--minio_path path' is not provided (i.e local run)

                Default = '0'"),
  #
  make_option(c("-a", "--max_global_variables"), default="10000",
              help="Indicates maximum allowed total size (in bytes) of global variables identified. Used by library(future) to prevent too large exports
              
                Default = '10000' for 10000 MiB")
)

opt <- parse_args(OptionParser(option_list=option_list))

InputCounts             <- opt$input_counts
InputCountsType         <- opt$input_counts_type
InputAnnotations        <- opt$input_annotations
NormalCellTypes         <- opt$normal_cell_types
InputCoordinates        <- opt$input_gene_coordinates
GeneMeanCountCutoff     <- as.numeric(opt$cutoff)
NoiseFilter             <- as.numeric(opt$noise_filter)
SdAmplifier             <- as.numeric(opt$sd_amplifier)
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
NumbCores               <- opt$number_cores
RunsCwl                 <- as.numeric(opt$run_cwl)
MaxGlobalVariables      <- as.numeric(opt$max_global_variables)

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
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

FILE_TYPE_OUT_DIRECTORIES = c(
  "INFERCNV",
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

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles,".", ProgramOutdir, "_INFERCNV_UsedOptions", ".txt")

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

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_INFERCNV_RSessionInfo", ".txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

####################################
### Define number of cores and RAM for parallelization
####################################

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (regexpr("^[0-9]+$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste0("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript ", ThisScriptName, " -h\n\n"))
}

if (NumbCoresToUse == 1) {
  plan(strategy = "sequential", workers = NumbCoresToUse)
  writeLines(paste0("\n", "*** Running: ", ThisScriptName, " in 'sequential' mode with ", NumbCoresToUse, " core ***", "\n"))
}else if (NumbCoresToUse > 1) {
  plan(strategy = "multicore", workers = NumbCoresToUse)
  writeLines(paste0("\n", "*** Running: ", ThisScriptName, " in 'multicore' mode with ", NumbCoresToUse, " cores ***", "\n"))
}else{
  stop(paste0("Unexpected --number_cores = ", NumbCoresToUse))
}

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4GiB of global variables identified (4000*1024^2), use: `options(future.globals.maxSize = 4000 * 1024^2)`
options(future.globals.maxSize = MaxGlobalVariables * 1024^2)

####################################
### Define default parameters
####################################

### No default parameters are defined

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
  counts.df <- as.data.frame(Read10X(data.dir = InputCounts, strip.suffix = T)) ### Note `strip.suffix = T` applies to Seurat v3.2 or higher
}else if (regexpr("^DGE$", InputCountsType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  counts.df <- data.frame(fread(InputCounts),row.names=1, check.names = FALSE)
}else{
  stop(paste("Unexpected type of --input_counts: ", InputCountsType, "\n\nFor help type:\n\nRscript Runs_InferCNV.R -h\n\n", sep=""))
}

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Creating InferCNV object
####################################
writeLines("\n*** Creating InferCNV object ***\n")

StopWatchStart$CreateInfercnvObject  <- Sys.time()

NormalCellTypes.ls = unlist(strsplit(NormalCellTypes, ","))

print("The following cell types will be used as normal:")
print(NormalCellTypes.ls)

print("Creatinf The following cell types will be used as normal:")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts.df,
                                    annotations_file=InputAnnotations,
                                    delim = "\t",
                                    gene_order_file=InputCoordinates,
                                    ref_group_names=c(NormalCellTypes.ls)
                                    )


StopWatchEnd$CreateInfercnvObject  <- Sys.time()

####################################
### Run infercnv
####################################
writeLines("\n*** Run infercnv ***\n")

StopWatchStart$InfercnvRun  <- Sys.time()

### Note using paths such as "~/temp" or "/Users/UserName/temp" for out_dir produces:
### (Error in nls(y ~ .logistic_midpt_slope(x, midpt = x0, slope = k), data = df, : 
### number of iterations exceeded maximum of 50), couldn't fit logistic, but no worries, going to use a spline
### Error in if (runif(1) <= padj) { : missing value where TRUE/FALSE needed`
### Thus, we'll use tempfile()

if (length(grep('^NA$',opt[["noise_filter"]], perl = T))) {
  HeatmapHeader<-paste(" cutoff=", GeneMeanCountCutoff, " noise_filter=NA sd_amplifier=", as.numeric(opt$sd_amplifier), " ", sep = "", collapse = "")
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=GeneMeanCountCutoff,
                               out_dir=tempfile(), 
                               cluster_by_groups=T,
                               plot_steps = T,
                               scale_data = T,
                               denoise=T,
                               sd_amplifier=as.numeric(opt$sd_amplifier),
                               HMM = T,
                               HMM_type = "i6"
                               )

}else{
  HeatmapHeader<-paste("cutoff=", GeneMeanCountCutoff, " noise_filter=", as.numeric(opt$noise_filter), " sd_amplifier=NA", sep = "", collapse = "")

  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=GeneMeanCountCutoff,
                               out_dir=tempfile(), 
                               cluster_by_groups=T,
                               plot_steps = T,
                               scale_data = T,
                               denoise=T,
                               noise_filter=as.numeric(opt$noise_filter),
                               HMM = T,
                               HMM_type = "i6",
                               num_threads = NumbCoresToUse
                               )
}

StopWatchEnd$InfercnvRun  <- Sys.time()

####################################
### Generate matrix and heatmap
####################################
writeLines("\n*** Generate matrix and heatmap ***\n")

StopWatchStart$PlotCnv  <- Sys.time()

plot_cnv(infercnv_obj = infercnv_obj, out_dir = paste0(Tempdir, "/INFERCNV"), title = HeatmapHeader, output_format = "pdf" )

StopWatchEnd$PlotCnv  <- Sys.time()

####################################
### Transform *txt files into *tsv
####################################
writeLines("\n*** Transform *txt files into *tsv ***\n")

StopWatchStart$TransformTxtToTsv  <- Sys.time()

ListOfTxtFilesToTsv<- c(paste0(Tempdir, "/INFERCNV/", "infercnv.references.txt"),
                        paste0(Tempdir, "/INFERCNV/", "infercnv.observations.txt"),
                        paste0(Tempdir, "/INFERCNV/", "infercnv.observation_groupings.txt")
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
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_INFERCNV_CPUtimes.txt")
write(file = OutfileCPUusage, x = paste("Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
write(file = OutfileCPUusage, x = paste("MaxGlobalVariables", MaxGlobalVariables, sep = "\t", collapse = ""))

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
    for (substep in rownames(InputsTable)) {
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
writeLines(paste0("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report"))

quit()
