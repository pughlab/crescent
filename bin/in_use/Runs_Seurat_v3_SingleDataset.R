####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
### and https://satijalab.org/seurat/v3.0/sctransform_vignette.html
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads scRNA-seq data and generate QC plots
### 2) Normalizes and scales measurements
### 3) Runs dimension reduction
### 4) Determines cell clusters
### 5) Runs average gene expression
### 6) Draws UMAP/tSNE plots using cell clusters, requested genes and metadata
### 7) Runs differential gene expression (DGE)
### 8) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_Seurat_v3_SingleDataset.R -h'
### for help
####################################

####################################
### THINGS TO DO:
### 1) Parallelize FindAllMarkers() to send each cluster to a separate core
###    Although it seems that FindMarkers() is already parallelized by Seurat:
###    https://satijalab.org/seurat/v3.0/future_vignette.html
### 2) Determine automatically the inflection point from the ElbowPlot() and use it as `dims` for FindNeighbors(), RunUMAP() and RunTSNE()
### 3) Determine automatically the optimal resolution from FindClusters() following Iness and Bader (F1000Research, 2019)
###    By picking the number of clusters based on differentially expressed genes. See https://f1000research.com/articles/7-1522/v1
### 4) Add a lists of ENSEMBL Ids for mitochondrial genes instead of just MT- and mt- (at gene names)
###    Need to do it for both Human and Mouse
### 5) To implement a flag to see if inputted dataset matches expected mito.fraction
###    In particular, if inputting a whole cell sample and using `-m 0,0.05` likely will filterout most cells
###    and the script may crash at step 'Perform linear dimensional reduction by PCA'
###    because the filtered matrix will be too sparse and small to get the gene PC's
### 6) Make *MarkersPerCluster.tsv and *TopTwoMarkersPerCluster.tsv a single file
###
### THINGS NICE TO HAVE:
### 1) Assigning cell type identity to clusters (needs supervised annotations, maybe based on GSVA)
### 2) Use knitr() to produce html plot layouts (https://yihui.name/knitr/demo/stitch/)
### 3) In "Load data" we use Seurat's Read10X() function. In this command, when all barcodes come from the same sample (i.e. finish with the same digit), like:
###    CTCTACGCAAGAGGCT-1
###    CTGAAACCAAGAGGCT-1
###    CTGAAACCAAGAGGCT-1
###    ... etc
###    Read10X will remove the '-digit'
###
###    Hence we need to implement code to remove the '-digit' from --input_clusters barcode ID's as well WHEN all barcodes come from the same sample
###    For now, this script is removing the digit always. And user must provide the inputs like:
###    1-CTCTACGCAAGAGGCT
###    2-CTCGAAAAGCTAACAA
###    3-CTGCCTAGTGCAGGTA
###    Instead of:
###    CTCTACGCAAGAGGCT-1
###    CTCGAAAAGCTAACAA-2
###    CTGCCTAGTGCAGGTA-3
####################################

####################################
### COMMENTS ON FINDING DIFFERENTIALLY EXPRESSED GENES
### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object.each_sample, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object.each_sample, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster5.markers, n = 5))
###
### NOTES:
### 1) FindAllMarkers() uses return.thresh = 0.01 as defaults, but FindMarkers() displays all genes passing previous filters.
###    Thus to make the outputs between these two commands identical to each other use return.thresh = 1
###
### 2) Default pseudocount.use=1, which sounds high for a Log level correction
###    An earlier version of this script was using 1e-99, but it was probably too small
###    Now using the inverse of the number of cells in the data.
###    This is sufficiently small as to not compress logGER magnitudes,
###    while keeping comparisons with zero reasonably close to the range of potential logGER values (Innes and Bader, 2018, F1000 Research)
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(DropletUtils)) # (Bioconductor) to handle MTX/H5 format files. Note it has about the same speed than library(earlycross) which can't handle H5
suppressPackageStartupMessages(library(Seurat))       # (CRAN) to run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # (CRAN) to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(staplr))       # (CRAN) to merge pdf files. Note it needs pdftk available. If not available use `SummaryPlots <- "N"`
suppressPackageStartupMessages(library(gtools))       # (CRAN) to do alphanumeric sorting. Only needed if using `-w Y`.
suppressPackageStartupMessages(library(loomR))        # (GitHub mojaveazure/loomR) needed for fron-end display of data. Only needed if using `-w Y`.
####################################

####################################
### Required external packages
####################################
### 'pdftk'   to merge selected *pdf files into *summary_plots.pdf using libary(staplr)
###           and to add statistics to violin plots overlapping pdf files
###           in Mac you can install it with something like:
###           'sudo port install pdftk'
###           https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/
### 'convert' to transofrm *pdf files into *png
###           it's part of ImageMagic (https://imagemagick.org/index.php)
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_Seurat_v3_SingleDataset.R"
ProgramOutdir  <- "SEURAT"

####################################
### Get inputs from command line argumets
####################################
#
option_list <- list(
  make_option(c("-i", "--input"), default="NA",
              help="Path/name to either a read counts matrix in either 'MTX' or 'TSV' format (see parameter -t)
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--input_type"), default="NA",
              help="Either 'MTX', 'TSV' or 'HDF5'
                'MTX'  is the path/name to an MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
                'TSV'  is the path/name of a <tab> delimited *file* with genes in rows vs. barcodes in columns
                'HDF5' is the path/name of a *file* in hdf5 format (e.g. from Cell Ranger)
                Note 'MTX' files can be the outputs from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-j", "--inputs_remove_barcodes"), default="NA",
              help="Path/name to a <tab> delimited list of barcodes to be removed from analysis, like:
                dataset1_id  AAACCTGAGCTCCCAG
                dataset2_id  AAACCTGTCACCATAG
                dataset3_id  AAACCTGTCAGCTTAG
                Or type 'NA' to include all barcodes
                Default = 'NA'"),
  #
  make_option(c("-b", "--normalize_and_scale_sample"), default="2",
              help="Indicates one of three options:
                a) Type '1' if --input contains raw counts (e.g. from cellranger) and this script should use Seurat's NormalizeData() function
                b) Type '2' if --input contains raw counts (e.g. from cellranger) and this script should use Seurat's SCTransform() function
                c) Type '3' if --input contains pre-normalized counts (e.g. transcripts per million counts [TPMs]) and the normalization step should be skipped
                Default = '2'"),
  #
  make_option(c("-k", "--save_filtered_data"), default="N",
              help="Indicates if filtered raw and normalized data should be saved as MTX files. Type 'y/Y' or 'n/N'
                Default = 'N'"),
  #
  make_option(c("-l", "--save_unfiltered_data"), default="N",
              help="Indicates if unfiltered raw data (exactly as inputted by --input) should be saved as MTX files. Type 'y/Y' or 'n/N'
                This may be useful to share with collaborators along results in a single zipped file
                Default = 'N'"),
  #
  make_option(c("-r", "--resolution"), default="1",
              help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of cell clusters
                Default = '1'"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-c", "--infile_metadata"), default="NA",
              help="A <tab> delimited table of barcodes and discrete properties to colour the dimension reduction plots, like:
                Barcode                CellClass    InOtherDatasets
                d1_AAACCTGAGCGGCTTC-1   1            yes
                d2_AAACCTGAGTCGAGTG-1   1            no
                d3_AAACCTGCAAAGGAAG-1   2            yes
                Note: the barcode id must include the dataset ID
                Default = 'NA' (i.e. no --infile_metadata is provided)"),
  #
  make_option(c("-g", "--list_genes"), default="NA",
              help="A <comma> delimited list of gene identifiers whose expression will be mapped into the dimension reduction plots
                Default = 'NA' (no --list_genes provided)"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Number of PCA dimensions to use for clustering and dimension reduction functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file *PCElbowPlot.pdf,
                use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                Default = '10'"),
  #
  make_option(c("-f", "--apply_cell_filters"), default="Y",
              help="Indicates if parameters `-m, -q, -n, -v` should be applied. Type 'y/Y' or 'n/N')
                Default = 'Y'"),
  #
  make_option(c("-m", "--percent_mito"), default="0,0.05",
              help="<comma> delimited min,max fraction of gene counts of mitochondrial origin a cell to be included in normalization and clustering analyses
                For example, for whole cell scRNA-seq use '0,0.2', or for Nuc-seq use '0,0.05'
                For negative values (e.g. if using TPM in log scale refer negative values with a double backslash '\\', like this '\\-1,0.5')
                Default = '0,0.05'"),
  #
  make_option(c("-q", "--percent_ribo"), default="0,0.75",
              help="<comma> delimited min,max fraction of gene counts of ribosomal proteins to be included in normalization and clustering analyses
                For negative values (e.g. if using TPM in log scale refer negative values with a double backslash '\\', like this '\\-1,0.5')
                Default = '0,0.75'"),
  #
  make_option(c("-n", "--n_genes"), default="50,8000",
              help="<comma> delimited min,max number of unique genes measured in a cell to be included in normalization and clustering analyses
                Default = '50,8000'"),
  #
  make_option(c("-v", "--n_reads"), default="1,80000",
              help="<comma> delimited min,max number of reads measured in a cell to be included in normalization and clustering analyses
                Default = '1,80000'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each differentially expressed genes test returns only hits that have an UNcorrected p-value < return_thresh
                Default = '0.01'"),
  #
  make_option(c("-u", "--number_cores"), default="AUTO",
              help="Indicates one of three options:
                a) Type the number of cores to use for parellelization (e.g. '4')
                b) Type 'AUTO' to detect the number of cells in the sample and assign a pre-established number of cores based on that
                c) Type 'MAX' to determine and use all available cores in the system
                Default = 'AUTO'"),
  #
  make_option(c("-s", "--save_r_object"), default="N",
              help="Indicates if a R object with the data and analyzes from the run should be saved
                Note that this may be time consuming. Type 'y/Y' or 'n/N'
                Default = 'N'"),
  #
  make_option(c("-w", "--run_cwl"), default="N",
              help="Indicates if this script is running inside a virtual machine container, such that outfiles are written directly into the 'HOME' . Type 'y/Y' or 'n/N'.
                Note, if using 'y/Y' this supersedes option -o
                Default = 'N'"),
  #
  make_option(c("-a", "--max_global_variables"), default="AUTO",
              help="Indicates one of two options:
                a) Type the maximum allowed total size (in bytes) of global variables identified for library(future) (e.g. '4000')
                b) Type 'AUTO' to detect the number of cells in the sample and assign a pre-established size of global variables
                Default = 'AUTO'")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input                   <- opt$input
InputType               <- opt$input_type
InfileRemoveBarcodes    <- opt$inputs_remove_barcodes
NormalizeAndScale       <- as.numeric(opt$normalize_and_scale_sample)
SaveFilteredData        <- opt$save_filtered_data 
SaveUnFilteredData      <- opt$save_unfiltered_data 
Resolution              <- as.numeric(opt$resolution)
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
InfileMetadata          <- opt$infile_metadata
ListGenes               <- opt$list_genes
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
ApplyCellFilters        <- opt$apply_cell_filters
ListPMito               <- opt$percent_mito
ListNGenes              <- opt$n_genes
ListPRibo               <- opt$percent_ribo
ListNReads              <- opt$n_reads
ThreshReturn            <- as.numeric(opt$return_threshold)
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl
MaxGlobalVariables      <- opt$max_global_variables

####################################
### Define default parameters
####################################
### Some of these default parameters are provided by Seurat developers,
### others are tailored according to clusters/t-SNE granularity

### To be able to take negative values for from cell filtering parameters
### entered using optparse (e.g. -m \\-1,1 means mito.percentage from -1 to 1)
ListNReads <- gsub("\\\\","", ListNReads, ignore.case = T)
ListNGenes <- gsub("\\\\","", ListNGenes, ignore.case = T)
ListPMito  <- gsub("\\\\","", ListPMito, ignore.case = T)
ListPRibo  <- gsub("\\\\","", ListPRibo, ignore.case = T)

ListNGenes = unlist(strsplit(ListNGenes, ","))
MinGenes   = as.numeric(ListNGenes[1])
MaxGenes   = as.numeric(ListNGenes[2])
#
ListPMito  = unlist(strsplit(ListPMito,  ","))
MinPMito   = as.numeric(ListPMito[1])
MaxPMito   = as.numeric(ListPMito[2])
#
ListNReads = unlist(strsplit(ListNReads, ","))
MinReads   = as.numeric(ListNReads[1])
MaxReads   = as.numeric(ListNReads[2])
#
ListPRibo  = unlist(strsplit(ListPRibo,  ","))
MinPRibo   = as.numeric(ListPRibo[1])
MaxPRibo   = as.numeric(ListPRibo[2])

DefaultParameters <- list(
  
  ### Parameters for QC plots
  CellPropertiesToQC = c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction"),
  
  ### Parameters for Seurat filters
  MinCells = 3,
  MinReads = MinReads,
  MaxReads = MaxReads,
  MinGenes = MinGenes,
  MaxGenes = MaxGenes,
  MinPMito = MinPMito,
  MaxPMito = MaxPMito,
  MinPRibo = MinPRibo,
  MaxPRibo = MaxPRibo,
  
  ### Parameters to assign number of cores
  NumbCoresSmall = 4,
  NumbCoresMedOrLarge = 10,
  MaxNumbCellsSmallForNumbCores = 30000,
  
  ### Parameters to assign size of global variables
  NumbGlobVarsSmall = 4000,
  NumbGlobVarsMedOrLarge = 16000,
  MaxNumbCellsSmallForGlobVars = 30000,
  
  ### Parameters for Seurat normalization
  ScaleFactor = 1000000, ### Using 1000000 to set scale.factor as counts per million (CPM)
  NormalizationMethod = "LogNormalize",
  
  ### Parameters for Seurat variable gene detection
  VarGeneDetectSelectMethod = "mean.var.plot",
  XLowCutoff = 0.0125,
  XHighCutoff = 3,
  YCutoff = 0.5,
  
  ### Parameters for PCA
  PrintPCA.PcsPrint = 1:5,
  PrintPCA.GenesPrint = 5,
  VizPCA.PcsUse = 1:6,
  VizPCA.nGenesToPlot = 20,
  PCHeatmapCellsUse = 300,
  PCHeatmapComponentsToPlot = 18,
  
  ### Parameters for Cluster Biomarkers
  FindAllMarkers.MinPct     =  0.25,
  FindAllMarkers.ThreshUse  =  0.25,
  NumberOfGenesPerClusterToPlotTsne  =  2,
  NumberOfGenesPerClusterToPlotHeatmap  =  10,
  
  ### Parameters for dimmension reduction plots
  MinNumberOfCellsToReducePerplexity = 150,
  ReducedPerplexity = 7,
  BaseSizeSinglePlotPdf  = 7,
  BaseSizeSinglePlotPng  = 480,
  BaseSizeMultiplePlotPdfWidth  = 3.7,
  BaseSizeMultiplePlotPdfHeight = 3,
  MaxNumbLabelsPerRowInLegend   = 4
  
)

### Colour definitions
ColourDefinitions<-list("orange"        = "#E69F00",
                        "bluishgreen"   = "#009E73",
                        "reddishpurple" = "#CC79A7",
                        "dodgerblue"    = "#1E90FF",
                        "vermillion"    = "#D55E00",
                        "snow4"         = "#8B8989",
                        "yellow"        = "#FFD700",
                        "seagreen3"     = "#43CD80",
                        "pink"          = "#FFC0CB",
                        "skyblue"       = "#56B4E9",
                        "orchid4"       = "#8B4789",
                        "blue"          = "#0072B2",
                        "black"         = "#000000"
)
ColoursQCViolinPlots <- c(ColourDefinitions[["skyblue"]][[1]], ColourDefinitions[["orange"]][[1]])

### Dimension reduction methods
DimensionReductionMethods<-list()
DimensionReductionMethods$umap$name <-"UMAP"
DimensionReductionMethods$tsne$name <-"TSNE"
DimensionReductionMethods$umap$run  <-as.function(RunUMAP)
DimensionReductionMethods$tsne$run  <-as.function(RunTSNE)

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

ListMandatory<-list("input", "input_type", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
  }
}

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Using `-w Y` will make Tempdir, which takes the value of ProgramOutdir, and it will be the final out-directory
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) 
  
  FILE_TYPE_OUT_DIRECTORIES = c(
    "CRESCENT_CLOUD",
    "CRESCENT_CLOUD/frontend_normalized",
    "CRESCENT_CLOUD/frontend_coordinates",
    "CRESCENT_CLOUD/frontend_raw",
    "CRESCENT_CLOUD/frontend_metadata",
    "CRESCENT_CLOUD/frontend_markers",
    "CRESCENT_CLOUD/frontend_qc",
    "CRESCENT_CLOUD/frontend_groups",
    "AVERAGE_GENE_EXPRESSION_TABLES", 
    "CELL_CLUSTER_IDENTITIES", 
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES",
    "DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS",
    "DIMENSION_REDUCTION_COORDINATE_TABLES",
    "DIMENSION_REDUCTION_PLOTS",
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS", 
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS", 
    "SUMMARY_PLOTS",
    "UNFILTERED_DATA_MATRICES"
  )
  
}else{
  PrefixOutfiles <- c(paste0(PrefixOutfiles,"_res",Resolution))
  ## Using `Tempdir/DIRECTORY` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## 'DIRECTORY' is one of the directories specified at FILE_TYPE_OUT_DIRECTORIES
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  Tempdir        <- "~/temp" 
  #
  CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
  UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
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
  
  FILE_TYPE_OUT_DIRECTORIES = c(
    "AVERAGE_GENE_EXPRESSION_TABLES", 
    "CELL_CLUSTER_IDENTITIES", 
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES",
    "DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS",
    "DIMENSION_REDUCTION_COORDINATE_TABLES",
    "DIMENSION_REDUCTION_PLOTS",
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS", 
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS", 
    "SUMMARY_PLOTS",
    "UNFILTERED_DATA_MATRICES"
  )
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

cat(file = OutfileOptionsUsed, x=paste0("\n", "One-line-commands used:", "\n", "`Rscript /path_to/", ThisScriptName), append = T)
for (optionInput in option_list) {
  cat(file = OutfileOptionsUsed, x=(paste0(" ", optionInput@short_flag, " ", opt[optionInput@dest])),append = T)
}
cat(file = OutfileOptionsUsed, x="`\n", append = T)

StopWatchEnd$ReportUsedOptions  <- Sys.time()

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/" , PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo.txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

StopWatchStart$LoadScRNAseqData  <- Sys.time()

if (regexpr("^MTX$", InputType, ignore.case = T)[1] == 1) {
  writeLines(paste0("\n*** Loading MTX infiles ***\n"))
  input.matrix <- Read10X(data.dir = Input)
}else if (regexpr("^TSV$", InputType, ignore.case = T)[1] == 1) {
  writeLines(paste0("\n*** Loading matrix of genes (rows) vs. barcodes (columns) ***\n"))
  ## Note `check.names = F` is needed for both `fread` and `data.frame`
  input.matrix <- as.matrix(data.frame(fread(Input, check.names = F), row.names=1, check.names = F))
}else if (regexpr("^HDF5$", InputType, ignore.case = T)[1] == 1) {
  writeLines(paste0("\n*** Loading HDF5 infile ***\n"))
  input.matrix <- Read10X_h5(filename = Input, use.names = T, unique.features = T)
}else{
  stop(paste0("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n"))
}
dim(input.matrix)

StopWatchEnd$LoadScRNAseqData  <- Sys.time()

####################################
### Save unfiltered data
####################################

if (regexpr("^Y$", SaveUnFilteredData, ignore.case = T)[1] == 1) {
  writeLines("\n*** Save unfiltered data ***\n")
  
  StopWatchStart$SaveUnFilteredData  <- Sys.time()

  ### This is to write out the input data with no filters at all
  seurat.object.tmp  <- CreateSeuratObject(counts = input.matrix, min.cells = 0, min.features = 0, project = PrefixOutfiles)
  
  OutDirUnfilteredRaw <-paste0(Tempdir, "/UNFILTERED_DATA_MATRICES/", "RAW")
  dir.create(file.path(OutDirUnfilteredRaw), showWarnings = F, recursive = T)
  write10xCounts(path = OutDirUnfilteredRaw, x = seurat.object.tmp@assays[["RNA"]]@counts, gene.type="Raw_Gene_Expression", overwrite=T, type="sparse", version="3")

  StopWatchEnd$SaveUnFilteredData  <- Sys.time()

}

####################################
### Define number of cores for parallelization
### Note: this must run after loading scRNA-seq data to get the number of barcodes (if using `-u AUTO`)
####################################
writeLines("\n*** Define number of cores for parallelization ***\n")

NumbCoresAvailable <- as.numeric(availableCores()[[1]]) ## Number of cores available in the system
NumberOfBarcodes <- ncol(input.matrix)

### Get number of cores requested
if (regexpr("^AUTO$", NumbCores, ignore.case = T)[1] == 1) {
  if (NumberOfBarcodes <= DefaultParameters$MaxNumbCellsSmallForNumbCores) {
    NumbCoresRequested <-DefaultParameters$NumbCoresSmall
  }else if (NumbCoresAvailable < DefaultParameters$NumbCoresMedOrLarge) {
    NumbCoresRequested <- NumbCoresAvailable
  }else{
    NumbCoresRequested <-DefaultParameters$NumbCoresMedOrLarge
  }
}else if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresRequested <- NumbCoresAvailable
}else if (regexpr("^[0-9]+$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresRequested <- as.numeric(NumbCores)
}else{
  stop(paste0("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript", ThisScriptName, " -h\n\n"))
}

### Check if number of cores requested is not larger than number of cores available
if (NumbCoresAvailable < NumbCoresRequested) {
  print(paste0("WARNING: Parameter `-u` requested ", NumbCoresRequested, " cores to process ", NumberOfBarcodes, " cells. But only ", NumbCoresAvailable, " cores are available and they will be used instead"))
  NumbCoresToUse <- NumbCoresAvailable
}else{
  NumbCoresToUse <- NumbCoresRequested
}

writeLines(paste0("\nUsing ", NumbCoresToUse, " cores\n"))

####################################
### Define library(future) strategy for parallelization
### Note: this must run after:
###       a) loading scRNA-seq data to get the number of barcodes (if using `-a AUTO`)
###       b) defining `NumbCoresToUse`
####################################
writeLines("\n*** Define library(future) strategy for parallelization ***\n")

NumberOfBarcodes <- ncol(input.matrix)

if (regexpr("^AUTO$", MaxGlobalVariables, ignore.case = T)[1] == 1) {
  if (NumberOfBarcodes <= DefaultParameters$MaxNumbCellsSmallForGlobVars) {
    MaxGlobalVariablesToUse <- DefaultParameters$NumbGlobVarsSmall
  }else{
    MaxGlobalVariablesToUse <- DefaultParameters$NumbGlobVarsMedOrLarge
  }
}else if (regexpr("^[0-9]+$", MaxGlobalVariables, ignore.case = T)[1] == 1) {
  MaxGlobalVariablesToUse <- as.numeric(MaxGlobalVariables)
}else{
  stop(paste0("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript ", ThisScriptName, " -h\n\n"))
}

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4Gb of RAM (4000*1024^2), use:

writeLines(paste0("\nUsing ", MaxGlobalVariablesToUse, " global variables\n"))

options(future.globals.maxSize = MaxGlobalVariablesToUse * 1024^2)

plan(strategy = "multicore", workers = NumbCoresToUse)

####################################
### Create a Seurat object
####################################
writeLines("\n*** Create a Seurat object ***\n")

StopWatchStart$CreateSeuratObject  <- Sys.time()

seurat.object.u  <- CreateSeuratObject(counts = input.matrix, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = PrefixOutfiles)
nCellsInOriginalMatrix<-length(seurat.object.u@meta.data$orig.ident)

StopWatchEnd$CreateSeuratObject  <- Sys.time()

####################################
### Get mitochondrial genes
####################################
writeLines("\n*** Get  mitochondrial genes ***\n")

StopWatchStart$GetMitoGenes  <- Sys.time()

mitoRegExpressions<- paste(c("^MT-"),collapse = "|")
mito.features <- grep(pattern = mitoRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)

if (length(mito.features)[[1]] > 0) {
  mito.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['mito.fraction']] <- mito.fraction
}else{
  mito.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['mito.fraction']] <- mito.fraction
}

StopWatchEnd$GetMitoGenes  <- Sys.time()

####################################
### Get ribosomal protein genes
####################################
writeLines("\n*** Get ribosomal protein genes ***\n")

StopWatchStart$GetRiboGenes  <- Sys.time()

riboRegExpressions<- paste(c("^MRPL", "^MRPS", "^RPL", "^RPS"),collapse = "|")
ribo.features <- grep(pattern = riboRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)

if (length(ribo.features)[[1]] > 0) {
  ribo.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['ribo.fraction']] <- ribo.fraction
}else{
  ribo.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
  seurat.object.u[['ribo.fraction']] <- ribo.fraction
}

StopWatchEnd$GetRiboGenes  <- Sys.time()

####################################
### Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation (if applicable)
####################################

if (regexpr("^Y$", ApplyCellFilters, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation ***\n")
  
  StopWatchStart$FilterCells  <- Sys.time()
  
  if (length(mito.features)[[1]] > 0) {
    seurat.object.f<-subset(x = seurat.object.u, subset = 
                              nFeature_RNA >= DefaultParameters$MinGenes
                            & nFeature_RNA <= DefaultParameters$MaxGenes 
                            & nCount_RNA   >= DefaultParameters$MinReads
                            & nCount_RNA   <= DefaultParameters$MaxReads 
                            & mito.fraction >= DefaultParameters$MinPMito
                            & mito.fraction <= DefaultParameters$MaxPMito
                            & ribo.fraction >= DefaultParameters$MinPRibo
                            & ribo.fraction <= DefaultParameters$MaxPRibo)
    
  }else{
    seurat.object.f<-subset(x = seurat.object.u, subset = 
                              nFeature_RNA >= DefaultParameters$MinGenes
                            & nFeature_RNA <= DefaultParameters$MaxGenes 
                            & nCount_RNA   >= DefaultParameters$MinReads
                            & nCount_RNA   <= DefaultParameters$MaxReads 
                            & ribo.fraction >= DefaultParameters$MinPRibo
                            & ribo.fraction <= DefaultParameters$MaxPRibo)
  }
    
  ### Get list of barcodes excluded by nFeature_RNA, nCount_RNA or ribo.fraction
  NumberOfBarcodesExcludedByNFeature <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA >= DefaultParameters$MinGenes & nFeature_RNA <= DefaultParameters$MaxGenes))))
  NumberOfBarcodesExcludedByNReads   <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nCount_RNA   >= DefaultParameters$MinReads & nCount_RNA   <= DefaultParameters$MaxReads))))
  NumberOfBarcodesExcludedByMito     <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = mito.fraction >= DefaultParameters$MinPMito & mito.fraction <= DefaultParameters$MaxPMito))))
  NumberOfBarcodesExcludedByRibo     <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = ribo.fraction >= DefaultParameters$MinPRibo & ribo.fraction <= DefaultParameters$MaxPRibo))))

  print(paste0("NumberOfBarcodesExcludedByNFeature=", NumberOfBarcodesExcludedByNFeature))
  print(paste0("NumberOfBarcodesExcludedByNReads=", NumberOfBarcodesExcludedByNReads))
  print(paste0("NumberOfBarcodesExcludedByMito=", NumberOfBarcodesExcludedByMito))
  print(paste0("NumberOfBarcodesExcludedByRibo=", NumberOfBarcodesExcludedByRibo))

}else{
  writeLines("\n*** Not applying cell filters ***\n")
  
  seurat.object.f <- seurat.object.u
  
  NumberOfBarcodesExcludedByNFeature <- 0
  NumberOfBarcodesExcludedByNReads   <- 0
  NumberOfBarcodesExcludedByMito     <- 0
  NumberOfBarcodesExcludedByRibo     <- 0
  
  MinGenes <- "NA"
  MaxGenes <- "NA"
  MinReads <- "NA"
  MaxReads <- "NA"
  MinPMito <- "NA"
  MaxPMito <- "NA"
  MinPRibo <- "NA"
  MaxPRibo <- "NA"
  ListNGenes <- c(MinGenes, MaxGenes)
  ListNReads <- c(MinReads, MaxReads)
  ListPMito  <- c(MinPMito, MaxPMito)
  ListPRibo  <- c(MinPRibo, MaxPRibo)
  DefaultParameters$MinGenes <- MinGenes
  DefaultParameters$MaxGenes <- MaxGenes
  DefaultParameters$MinReads <- MinReads
  DefaultParameters$MaxReads <- MaxReads
  DefaultParameters$MinPMito <- MinPMito
  DefaultParameters$MaxPMito <- MaxPMito
  DefaultParameters$MinPRibo <- MinPRibo
  DefaultParameters$MaxPRibo <- MaxPRibo
  
}

StopWatchEnd$FilterCells  <- Sys.time()

### Just reporting the summary of the UNfiltered and filtered objects
print("Unfiltered object")
seurat.object.u
print("Filtered object")
seurat.object.f

####################################
### QC EDA violin plots
####################################
writeLines("\n*** QC EDA violin plots ***\n")

StopWatchStart$QCviolinplots  <- Sys.time()

QCStats <- list()
QCStats$unfiltered$mean$nFeature_RNA    <-round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$unfiltered$median$nFeature_RNA  <-round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
QCStats$unfiltered$mean$nCount_RNA      <-round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$unfiltered$median$nCount_RNA    <-round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)
QCStats$unfiltered$mean$mito.fraction   <-round(mean(seurat.object.u@meta.data[,"mito.fraction"]),3)
QCStats$unfiltered$median$mito.fraction <-round(median(seurat.object.u@meta.data[,"mito.fraction"]),3)
QCStats$unfiltered$mean$ribo.fraction   <-round(mean(seurat.object.u@meta.data[,"ribo.fraction"]),3)
QCStats$unfiltered$median$ribo.fraction <-round(median(seurat.object.u@meta.data[,"ribo.fraction"]),3)
#
QCStats$filtered$mean$nFeature_RNA      <-round(mean(seurat.object.f@meta.data[,"nFeature_RNA"]),0)
QCStats$filtered$median$nFeature_RNA    <-round(median(seurat.object.f@meta.data[,"nFeature_RNA"]),0)
QCStats$filtered$mean$nCount_RNA        <-round(mean(seurat.object.f@meta.data[,"nCount_RNA"]),0)
QCStats$filtered$median$nCount_RNA      <-round(median(seurat.object.f@meta.data[,"nCount_RNA"]),0)
QCStats$filtered$mean$mito.fraction     <-round(mean(seurat.object.f@meta.data[,"mito.fraction"]),3)
QCStats$filtered$median$mito.fraction   <-round(median(seurat.object.f@meta.data[,"mito.fraction"]),3)
QCStats$filtered$mean$ribo.fraction     <-round(mean(seurat.object.f@meta.data[,"ribo.fraction"]),3)
QCStats$filtered$median$ribo.fraction   <-round(median(seurat.object.f@meta.data[,"ribo.fraction"]),3)

### Get unfiltered data QC statistics
nFeature_RNA.u.df   <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
nCount_RNA.u.df     <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
mito.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$mito.fraction, mito.fraction = 1)
ribo.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$ribo.fraction, ribo.fraction = 1)
#
nFeature_RNAStats.u <-paste0(c(" mean = ", QCStats$unfiltered$mean$nFeature_RNA,"\n", "median = ",  QCStats$unfiltered$median$nFeature_RNA))
nCount_RNAStats.u   <-paste0(c(" mean = ", QCStats$unfiltered$mean$nCount_RNA,  "\n", "median = ",  QCStats$unfiltered$median$nCount_RNA))
mito.fraction.u     <-paste0(c(" mean = ", QCStats$unfiltered$mean$mito.fraction,"\n", "median = ", QCStats$unfiltered$median$mito.fraction))
ribo.fraction.u     <-paste0(c(" mean = ", QCStats$unfiltered$mean$ribo.fraction,"\n", "median = ", QCStats$unfiltered$median$ribo.fraction))

### Get filtered data QC statistics
nFeature_RNA.f.df   <-data.frame(Expression_level = seurat.object.f@meta.data$nFeature_RNA, nGenes = 2)
nCount_RNA.f.df     <-data.frame(Expression_level = seurat.object.f@meta.data$nCount_RNA,   nCount_RNA = 2)
mito.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$mito.fraction, mito.fraction = 2)
ribo.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$ribo.fraction, ribo.fraction = 2)
#
nFeature_RNAStats.f <-paste0(c(" mean = ", QCStats$filtered$mean$nFeature_RNA,"\n", "median = ", QCStats$filtered$median$nFeature_RNA))
nCount_RNAStats.f   <-paste0(c(" mean = ", QCStats$filtered$mean$nCount_RNA,  "\n", "median = ", QCStats$filtered$median$nCount_RNA))
mito.fraction.f     <-paste0(c(" mean = ", QCStats$filtered$mean$mito.fraction,"\n", "median = ", QCStats$filtered$median$mito.fraction))
ribo.fraction.f     <-paste0(c(" mean = ", QCStats$filtered$mean$ribo.fraction,"\n", "median = ", QCStats$filtered$median$ribo.fraction))

### Put QC statistics together
nFeature_RNA.m.df   <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
nCount_RNA.m.df     <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
mito.fraction.m.df  <-data.frame(rbind(mito.fraction.u.df,mito.fraction.f.df))
ribo.fraction.m.df  <-data.frame(rbind(ribo.fraction.u.df,ribo.fraction.f.df))
LabelUnfiltered     <-paste0("Before filters: No. of cells = ", nrow(seurat.object.u@meta.data))
LabelFiltered       <-paste0("After filters:  No. of cells = ", nrow(seurat.object.f@meta.data))
NumberOfCells<- list()
NumberOfCells[["unfiltered"]] <- nrow(seurat.object.u@meta.data)
NumberOfCells[["filtered"]]   <- nrow(seurat.object.f@meta.data)

### Commands for violin ggplot's
DataForHeader.df<-data.frame(forx = c(0.4,0.4), fory = c(0.09,0.03), label = c(LabelFiltered,LabelUnfiltered))
Headers.plot<-ggplot(data=DataForHeader.df, aes(x = forx, y = fory)) + theme_void() + 
  geom_point(colour = ColoursQCViolinPlots, size = 7) + xlim(0,1) + ylim(0,0.12) +
  geom_text(hjust = -0.08, label = c(LabelUnfiltered,LabelFiltered), size = 4)

nFeature_RNA.plot<-ggplot(data=nFeature_RNA.m.df, aes(x = factor(nGenes), y = Expression_level)) +
  geom_violin(aes(fill = factor(nGenes))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x=paste0("No. of genes", "\n", "Filter: ", "min=", ListNGenes[[1]], " max=", ListNGenes[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNFeature)) +
  annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])

nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
  geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x=paste0("No. of reads", "\n", "Filter: ", "min=", ListNReads[[1]], " max=", ListNReads[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNReads)) +
  annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])

mito.fraction.plot<-ggplot(data=mito.fraction.m.df, aes(x = factor(mito.fraction), y = Expression_level)) +
  geom_violin(aes(fill = factor(mito.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x="Mitochondrial genes (fraction)") +
  labs(x=paste0("Mitochondrial genes", "\n", "Filter: ", "min=", ListPMito[[1]], " max=", ListPMito[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByMito)) +
  annotate("text", x = 1 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.f, col = ColoursQCViolinPlots[[2]])

ribo.fraction.plot<-ggplot(data=ribo.fraction.m.df, aes(x = factor(ribo.fraction), y = Expression_level)) +
  geom_violin(aes(fill = factor(ribo.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  scale_fill_manual(values = ColoursQCViolinPlots) +
  labs(x=paste0("Ribosomal protein genes (fraction)", "\n", "Filter: ", "min=", ListPRibo[[1]], " max=", ListPRibo[[2]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByRibo)) +
  annotate("text", x = 1 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.f, col = ColoursQCViolinPlots[[2]])

bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, mito.fraction.plot, ribo.fraction.plot, ncol = 4)

### Create a *pdf file with the violin ggplot's
VlnPlotPdf<-paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_QC_VlnPlot.pdf")
pdf(file=VlnPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 1.7, height = DefaultParameters$BaseSizeSinglePlotPdf)
print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
dev.off()

StopWatchEnd$QCviolinplots  <- Sys.time()

####################################
### CWL interactive qc plots
####################################
writeLines("\n*** CWL interactive qc plots ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  # unfiltered
  interactive_qc_plot_u  <-data.frame(Barcodes = row.names(seurat.object.u@meta.data), Number_of_Genes = seurat.object.u@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.u@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.u@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.u@meta.data$ribo.fraction)
  interactive_qc_plot_u$Mitochondrial_Genes_Percentage <- interactive_qc_plot_u$Mitochondrial_Genes_Percentage * 100
  interactive_qc_plot_u$Ribosomal_Protein_Genes_Percentage <- interactive_qc_plot_u$Ribosomal_Protein_Genes_Percentage * 100
  colnames(interactive_qc_plot_u) <- c("Barcodes","Number of Genes","Number of Reads","Mitochondrial Genes Percentage","Ribosomal Protein Genes Percentage")
  write.table(interactive_qc_plot_u, paste(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",PrefixOutfiles,"_BeforeFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
  
  # filtered
  interactive_qc_plot_f  <-data.frame(Barcodes = row.names(seurat.object.f@meta.data), Number_of_Genes = seurat.object.f@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.f@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.f@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.f@meta.data$ribo.fraction)
  interactive_qc_plot_f$Mitochondrial_Genes_Percentage <- interactive_qc_plot_f$Mitochondrial_Genes_Percentage * 100
  interactive_qc_plot_f$Ribosomal_Protein_Genes_Percentage <- interactive_qc_plot_f$Ribosomal_Protein_Genes_Percentage * 100
  colnames(interactive_qc_plot_f) <- c("Barcodes","Number of Genes","Number of Reads","Mitochondrial Genes Percentage","Ribosomal Protein Genes Percentage")
  write.table(interactive_qc_plot_f, paste(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",PrefixOutfiles,"_AfterFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
  
  qc_tsv <- data.frame(NAME = row.names(seurat.object.f@meta.data), Number_of_Genes = seurat.object.f@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.f@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.f@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.f@meta.data$ribo.fraction)
  qc_tsv$Mitochondrial_Genes_Percentage <- qc_tsv$Mitochondrial_Genes_Percentage * 100
  qc_tsv$Ribosomal_Protein_Genes_Percentage <- qc_tsv$Ribosomal_Protein_Genes_Percentage * 100
  qc_tsv_string <- sapply(qc_tsv, as.character)
  qc_tsv_string_TYPE <- rbind(data.frame(NAME = "TYPE", Number_of_Genes = "numeric", Number_of_Reads = "numeric", Mitochondrial_Genes_Percentage = "numeric", Ribosomal_Protein_Genes_Percentage = "numeric"), qc_tsv_string)
  
  qc_outfile <-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",PrefixOutfiles,"_qc_data.tsv")
  write.table(data.frame(qc_tsv_string_TYPE),file = qc_outfile, row.names = F, col.names = T, sep="\t", quote = F, append = T)
} 

####################################
### Feature-vs-feature scatter plot
####################################
writeLines("\n*** Feature-vs-feature scatter plot ***\n")

StopWatchStart$FeatureVsFeatureplot  <- Sys.time()

UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                              nGene = seurat.object.u@meta.data$nFeature_RNA,
                              mito.fraction = seurat.object.u@meta.data$mito.fraction,
                              filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.f))
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = ColoursQCViolinPlots[[1]])
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = F, replacement = ColoursQCViolinPlots[[2]])
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = 4)
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = F, replacement = 16)

FeatureVsFeaturePlotPdf<-paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_NumbReadsVsNumbGenesAndMito_ScatterPlot.pdf")
pdf(file=FeatureVsFeaturePlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 2, height = DefaultParameters$BaseSizeSinglePlotPdf)
par(mfrow=c(1,2))
## No. of reads vs. Mitochond. %
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$mito.fraction, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "Mitochond. %")
legend("topright", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))

## No. of reads vs. No. of genes
plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$nGene, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
     xlab ="No. of reads", ylab = "No. of genes")
legend("topleft", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))
dev.off()

StopWatchEnd$FeatureVsFeatureplot  <- Sys.time()

####################################
### Write out filter details and number of filtered cells
####################################
writeLines(paste0("\n*** Write out filter details and number of filtered cells ***\n"))

StopWatchStart$OutTablesFilterDetailsAndFilteredCells  <- Sys.time()

OutTableFilterDetails<-paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_FilterDetails.tsv")
Headers<-paste("Step", "Filter_min", "Filter_max", "Mean_before_filter", "Median_before_filter", "Mean_after_filter", "Median_after_filter", "Excluded_cells", sep = "\t", collapse = "")
write.table(Headers, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F)

FilterDetails.nFeature_RNA  <- paste("nFeature_RNA", ListNGenes[[1]], ListNGenes[[2]],
                                     QCStats$unfiltered$mean$nFeature_RNA,  QCStats$unfiltered$median$nFeature_RNA,
                                     QCStats$filtered$mean$nFeature_RNA,    QCStats$filtered$median$nFeature_RNA, 
                                     NumberOfBarcodesExcludedByNFeature, sep = "\t", collapse = "")
FilterDetails.nCount_RNA    <- paste("nCount_RNA", ListNReads[[1]], ListNReads[[2]], 
                                     QCStats$unfiltered$mean$nCount_RNA,    QCStats$unfiltered$median$nCount_RNA,
                                     QCStats$filtered$mean$nCount_RNA,      QCStats$filtered$median$nCount_RNA, 
                                     NumberOfBarcodesExcludedByNReads, sep = "\t", collapse = "")
FilterDetails.mito.fraction <- paste("mito.fraction", ListPMito[[1]], ListPMito[[2]],
                                     QCStats$unfiltered$mean$mito.fraction, QCStats$unfiltered$median$mito.fraction,
                                     QCStats$filtered$mean$mito.fraction,   QCStats$filtered$median$mito.fraction, 
                                     NumberOfBarcodesExcludedByMito, sep = "\t", collapse = "")
FilterDetails.ribo.fraction <- paste("ribo.fraction", ListPRibo[[1]], ListPRibo[[2]],
                                     QCStats$unfiltered$mean$ribo.fraction, QCStats$unfiltered$median$ribo.fraction,
                                     QCStats$filtered$mean$ribo.fraction,   QCStats$filtered$median$ribo.fraction, 
                                     NumberOfBarcodesExcludedByRibo, sep = "\t", collapse = "")

write.table(FilterDetails.nFeature_RNA,  file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(FilterDetails.nCount_RNA,    file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(FilterDetails.mito.fraction, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
write.table(FilterDetails.ribo.fraction, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)

OutTableFilteredCells<-paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_NumberOfFilteredCells.tsv")
write.table(paste("Number_of_cells_before_filters", NumberOfCells[["unfiltered"]], sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste("Number_of_cells_after_filters", NumberOfCells[["filtered"]],    sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$OutTablesFilterDetailsAndFilteredCells  <- Sys.time()

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  qc_metrics_genes <- data.frame(Step = "Number of Genes", Filter_min = ListNGenes[[1]], Filter_max = ListNGenes[[2]], Excluded_cells = NumberOfBarcodesExcludedByNFeature)
  qc_metrics_reads <- data.frame(Step = "Number of Reads", Filter_min = ListNReads[[1]], Filter_max = ListNReads[[2]], Excluded_cells = NumberOfBarcodesExcludedByNReads)
  qc_metrics_mito <- data.frame(Step = "Percentage of Mitochondrial Genes", Filter_min = as.numeric(ListPMito[[1]])*100, Filter_max = as.numeric(ListPMito[[2]])*100, Excluded_cells = NumberOfBarcodesExcludedByMito)
  qc_metrics_ribo <- data.frame(Step = "Percentage of Ribsomal Protein Genes", Filter_min = as.numeric(ListPRibo[[1]])*100, Filter_max = as.numeric(ListPRibo[[2]])*100, Excluded_cells = NumberOfBarcodesExcludedByRibo)
  
  qc_metrics <-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",PrefixOutfiles,"_qc_metrics.tsv")
  write.table(qc_metrics_genes, file = qc_metrics, row.names = F, col.names = T, sep="\t", quote = F, append = T)
  write.table(qc_metrics_reads, file = qc_metrics, row.names = F, col.names = F, sep="\t", quote = F, append = T)
  write.table(qc_metrics_mito, file = qc_metrics, row.names = F, col.names = F, sep="\t", quote = F, append = T)
  write.table(qc_metrics_ribo, file = qc_metrics, row.names = F, col.names = F, sep="\t", quote = F, append = T)
}

####################################
### Write out QC data
####################################
writeLines("\n*** Write out QC data ***\n")

StopWatchStart$WriteOutQCData  <- Sys.time()

Headers<-paste("Cell_barcode", paste(DefaultParameters$CellPropertiesToQC, sep = "", collapse = "\t") ,sep="\t")

BarcodeIdsBeforeFilters <- unlist(x = colnames(seurat.object.u))
OutfileQCMetadataBeforeFilters<-paste0(Tempdir,"/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", "Before_filters_QC_metadata.tsv")
write.table(Headers, file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(BarcodeIdsBeforeFilters, seurat.object.u@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)

BarcodeIdsAfterFilters <- unlist(x = colnames(seurat.object.f))
OutfileQCMetadataAfterFilters<-paste0(Tempdir,"/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", "After_filters_QC_metadata.tsv")
write.table(Headers, file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(BarcodeIdsAfterFilters, seurat.object.f@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$WriteOutQCData  <- Sys.time()

####################################
### Remove barcodes by parameter -j (if applicable)
####################################

if (regexpr("^NA$", InfileRemoveBarcodes , ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Ignoring option -j ***\n")
  
}else{
  
  StopWatchStart$RemoveBarcodes <- Sys.time()
  
  writeLines("\n*** Removing barcodes by parameter -j ***\n")
  
  InfileAllBarcodesToRemove<-gsub("^~/",paste0(c(UserHomeDirectory,"/")), InfileRemoveBarcodes)
  AllBarcodesToRemove.tab<-read.table(InfileAllBarcodesToRemove, header = F, row.names = NULL, stringsAsFactors = FALSE)
  colnames(AllBarcodesToRemove.tab) <- c("Barcode")
  
  seurat.object.full    <- seurat.object.f
  seurat.object.subset  <- subset(seurat.object.full, cells = setdiff(colnames(seurat.object.full),BarcodesToSubset.tab[,"Barcode"]))
  print(paste(paste0("Before:", ncol(seurat.object.full)), paste0("After:", ncol(seurat.object.subset)), sep = "  ", collapse = "\n"))
  seurat.object.f       <- seurat.object.subset
  
  StopWatchEnd$RemoveBarcodes <- Sys.time()
  
}

####################################
### Remove temporary Seurat objects
####################################
writeLines("\n*** Remove temporary Seurat objects ***\n")

rm(seurat.object.u)
rm(UnfilteredData.df)
rm(seurat.object.full)
rm(seurat.object.subset)


####################################
### Normalize data (if applicable)
####################################
if (NormalizeAndScale == 1) {
  
  writeLines("\n*** Normalize data using NormalizeData() ***\n")
  
  StopWatchStart$NormalizeData  <- Sys.time()
  
  seurat.object.f <- NormalizeData(object = seurat.object.f, normalization.method = DefaultParameters$NormalizationMethod, scale.factor = DefaultParameters$ScaleFactor)
  
  StopWatchEnd$NormalizeData  <- Sys.time()
  
}else if (NormalizeAndScale == 2) {
  
  writeLines("\n*** Normalize data using SCTransform() ***\n")
  
  StopWatchStart$SCTransform  <- Sys.time()
  
  seurat.object.f <- SCTransform(object = seurat.object.f, verbose = T)
  
  StopWatchEnd$SCTransform  <- Sys.time()
  
}else if (NormalizeAndScale == 3) {
  
  writeLines("\n*** Data assumed to be already normalized ***\n")
  
} else {
  stop(paste0("Unexpected option -b", NormalizeAndScale, ". Only options '1', '2', or '3' are allowed. \n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n"))
}

####################################
### Save filtered data
####################################

if (regexpr("^Y$", SaveFilteredData, ignore.case = T)[1] == 1) {
  writeLines("\n*** Save filtered data ***\n")
  
  OutDirFilteredRaw <-paste0(Tempdir, "/FILTERED_DATA_MATRICES/", "RAW")
  dir.create(file.path(OutDirFilteredRaw), showWarnings = F, recursive = T)
  write10xCounts(path = OutDirFilteredRaw, x = seurat.object.f@assays[["RNA"]]@counts, gene.type="Raw_Gene_Expression", overwrite=T, type="sparse", version="3")

  if (NormalizeAndScale == 1) {
    OutDirFilteredNorm<-paste0(Tempdir, "/FILTERED_DATA_MATRICES/", "NORMALIZED_LOG")
    dir.create(file.path(OutDirFilteredNorm), showWarnings = F, recursive = T)
    write10xCounts(path = OutDirFilteredNorm, x = round(seurat.object.f@assays[["RNA"]]@data, digits = 4), gene.type="LogNormalize_Gene_Expression", overwrite=T, type="sparse", version="3")
  }else if (NormalizeAndScale == 2) {
    OutDirFilteredNorm<-paste0(Tempdir, "/FILTERED_DATA_MATRICES/", "NORMALIZED_SCT")
    dir.create(file.path(OutDirFilteredNorm), showWarnings = F, recursive = T)
    write10xCounts(path = OutDirFilteredNorm, x = round(seurat.object.f@assays[["SCT"]]@data, digits = 4), gene.type="SCTransform_Gene_Expression", overwrite=T, type="sparse", version="3")
  }
}

####################################
### Save normalized count matrix as loom
####################################

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Save normalized count matrix as loom ***\n")
  
  normalized_count_matrix <- as.matrix(seurat.object.f@assays[["RNA"]]@data)
  
  # all genes/features in matrix
  features_tsv <- data.frame(features = rownames(normalized_count_matrix))
  features_tsv_ordered <- as.data.frame(features_tsv[mixedorder(features_tsv$features),])
  write.table(features_tsv_ordered, file=paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_raw/","features.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  # generating loom file of normalized count matrix
  loom_file <- paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_normalized/","normalized_counts.loom")
  create(loom_file, normalized_count_matrix)
  
}

####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes ***\n")

StopWatchStart$FindVariableGenes  <- Sys.time()

### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object.f <- FindVariableFeatures(object = seurat.object.f, selection.method = DefaultParameters$VarGeneDetectSelectMethod, mean.cutoff = c(DefaultParameters$XLowCutoff, DefaultParameters$XHighCutoff), dispersion.cutoff = c(DefaultParameters$YCutoff, Inf))
VariableGenes<-VariableFeatures(object = seurat.object.f)

StopWatchEnd$FindVariableGenes  <- Sys.time()

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation ***\n")

StopWatchStart$ScaleData  <- Sys.time()

seurat.object.f <- ScaleData(object = seurat.object.f, vars.to.regress = c("nCount_RNA", "mito.fraction"), features = rownames(x = seurat.object.f), verbose = T)

StopWatchEnd$ScaleData  <- Sys.time()

####################################
### Perform linear dimensional reduction by PCA
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA ***\n")

StopWatchStart$RunPCA  <- Sys.time()

seurat.object.f <- RunPCA(object = seurat.object.f, features = VariableGenes, verbose = T, do.print = T, ndims.print = DefaultParameters$PrintPCA.PcsPrint, nfeatures.print = DefaultParameters$PrintPCA.GenesPrint)

StopWatchEnd$RunPCA  <- Sys.time()

StopWatchStart$PCAPlots  <- Sys.time()

seurat.object.f <- ProjectDim(object = seurat.object.f, overwrite = T, verbose = T, nfeatures.print = 10)

StopWatchEnd$PCAPlots  <- Sys.time()

####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components ***\n")

StopWatchStart$GetSignificantPCs  <- Sys.time()

ForElbowPlot<-ElbowPlot(object = seurat.object.f, ndims = 50, reduction = "pca")
MaxYAxis<-as.integer(max(ForElbowPlot$data$stdev)+1)

pdf(file=paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_PCElbowPlot.pdf"))
print(ForElbowPlot
      + scale_x_continuous(breaks =  seq(from = 0, to = 50, by=5))
      + geom_vline(xintercept = seq(from = 0, to = 50, by=5), linetype='dotted', col="red")
      + scale_y_continuous(breaks =  seq(from = 0, to = MaxYAxis, by=0.5))
      + geom_hline(yintercept = seq(from = 0, to = MaxYAxis, by=0.5), linetype='dotted', col="red")
)
dev.off()

StopWatchEnd$GetSignificantPCs  <- Sys.time()

####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells ***\n")

StopWatchStart$ClusterCells  <- Sys.time()

options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.f <- FindNeighbors(object = seurat.object.f, dims = PcaDimsUse) ## This step was part of FindClusters() in Seurat v2
seurat.object.f <- FindClusters(object = seurat.object.f, resolution = Resolution)

StopWatchEnd$ClusterCells  <- Sys.time()

StopWatchStart$CellClusterTables  <- Sys.time()

ClusterIdent<-seurat.object.f@meta.data$seurat_clusters
NumberOfClusters<-length(unique(ClusterIdent))

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  metadata_tsv  <-data.frame(NAME = row.names(seurat.object.f@meta.data), seurat_clusters = seurat.object.f@meta.data$seurat_clusters)
  metadata_tsv_string <- sapply(metadata_tsv, as.character)
  metadata_tsv_string_TYPE <- rbind(data.frame(NAME = "TYPE", seurat_clusters = "group"), metadata_tsv_string)
  ### Note paste0() didn't work here. Use paste(...,  sep = "", collapse = "") instead
  colnames(metadata_tsv_string_TYPE) <- c("NAME", paste("Seurat_Clusters_Resolution", Resolution, sep = "", collapse = ""))
  OutfileClusters<-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_groups/","groups.tsv")
  write.table(data.frame(metadata_tsv_string_TYPE),file = OutfileClusters, row.names = F, col.names = T, sep="\t", quote = F, append = T)
  
}

CellNames<-rownames(seurat.object.f@meta.data)
Headers<-paste("Cell_barcode", paste0("seurat_cluster_resolution", Resolution) ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
OutfileClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles,".", ProgramOutdir, "_CellClusters.tsv")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
OutfileNumbClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles,".", ProgramOutdir, "_NumbCellClusters", ".tsv")
write(x=NumberOfClusters,file = OutfileNumbClusters)

StopWatchEnd$CellClusterTables  <- Sys.time()

####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene ***\n")

StopWatchStart$AverageGeneExpression  <- Sys.time()

cluster.averages<-AverageExpression(object = seurat.object.f, use.scale = F, use.counts = F)
OutfileClusterAverages<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles,".", ProgramOutdir, "_AverageGeneExpressionPerCluster.tsv")
Headers<-paste("AVERAGE_GENE_EXPRESSION", paste("C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$AverageGeneExpression  <- Sys.time()

####################################
### Run and plot dimension reductions
####################################
writeLines("\n*** Run and plot dimension reductions ***\n")

names(DimensionReductionMethods) <- c("umap", "tsne")

for (dim_red_method in names(DimensionReductionMethods)) {
  
  ####################################
  ### Run non-linear dimensional reductions
  ####################################
  writeLines(paste0("\n*** Run ", DimensionReductionMethods[[dim_red_method]][["name"]], " ***\n"))
  
  StopWatchStart$DimensionReduction$dim_red_method  <- Sys.time()
  
  ### NOTES:
  ### In RunTSNE: if the datasets is small user may get error:
  ### `Error in .check_tsne_params(nrow(X), dims = dims, perplexity = perplexity,  : 
  ### perplexity is too large for the number of samples`
  ### User can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
  ###
  ### Also using RunTSNE(..., check_duplicates = F) to skip cases where cells happen to have the same values after PCA reduction
  
  if (("tsne" %in% dim_red_method) & (length(colnames(seurat.object.f)) < DefaultParameters$MinNumberOfCellsToReducePerplexity)) {
    writeLines(paste0("\n*** Using reduced perplexity = ", DefaultParameters$ReducedPerplexity, " because found ",  length(colnames(seurat.object.f)), " cells", " ***\n"))
    seurat.object.f <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.f, dims = PcaDimsUse, perplexity = DefaultParameters$ReducedPerplexity, check_duplicates = F)
  }else if ("tsne" %in% dim_red_method) {
    seurat.object.f <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.f, dims = PcaDimsUse, check_duplicates = F)
  }else if ("umap" %in% dim_red_method) {
    seurat.object.f <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.f, dims = PcaDimsUse, umap.method = "uwot")
  }
  
  StopWatchEnd$DimensionReduction$dim_red_method  <- Sys.time()
  
  StopWatchStart$DimensionReductionPlot$dim_red_method  <- Sys.time()
  
  pdf(file=paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf"))
  print(DimPlot(object = seurat.object.f, reduction = dim_red_method, group.by = 'ident', label = T, label.size=10))
  dev.off()
  
  StopWatchEnd$DimensionReductionPlot$dim_red_method  <- Sys.time()
  
  ####################################
  ### Write out coordinates
  ####################################
  writeLines(paste0("\n*** Write out ", DimensionReductionMethods[[dim_red_method]][["name"]], " coordinates ***\n"))
  
  StopWatchStart$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  Headers<-paste("Barcode",paste(colnames(seurat.object.f@reductions[[dim_red_method]]@cell.embeddings),sep="",collapse="\t"), sep="\t", collapse = "\t")
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    OutfileCoordinatesCWL<-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_coordinates/",DimensionReductionMethods[[dim_red_method]][["name"]], "Coordinates.tsv")
    write.table(Headers,file = OutfileCoordinatesCWL, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(seurat.object.f@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinatesCWL,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
  }
  
  OutfileCoordinates<-paste0(Tempdir, "/DIMENSION_REDUCTION_COORDINATE_TABLES/", PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Coordinates.tsv")
  write.table(Headers,file = OutfileCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(seurat.object.f@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction ***\n"))
  
  StopWatchStart$QCDimRedPlots$dim_red_method  <- Sys.time()
  
  CellPropertiesToColour<-c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction")
  pdf(file=paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByQC.pdf"), width = DefaultParameters$BaseSizeSinglePlotPdf * 1.5, height = DefaultParameters$BaseSizeSinglePlotPdf * 1.5)
  print(FeaturePlot(object = seurat.object.f, label = T, order = T, features = CellPropertiesToColour, cols = c("lightgrey", "blue"), reduction = dim_red_method, ncol = 2, pt.size = 1.5))
  dev.off()
  
  StopWatchEnd$QCDimRedPlots$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by -infile_metadata
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by -infile_metadata ***\n"))
  
  if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
    print("No metadata will be used for dimension reduction plots")
  }else{
    
    StopWatchStart$DimRedPlotsColuredByOptC$dim_red_method  <- Sys.time()
    
    seurat.object.meta.data<-seurat.object.f@meta.data
    ExtraCellProperties <- data.frame(read.table(InfileMetadata, header = T, row.names = 1))
    print(head(ExtraCellProperties))
    
    # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
    # so that barcodes from --infile_colour_dim_red_plots and --input can match each other
    rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
    seurat.object.f <- AddMetaData(object = seurat.object.f, metadata = ExtraCellProperties)
    
    # Generating outfile
    # Note DimPlot() takes the entire current device (pdf)
    # even if using layout(matrix(...)) or  par(mfrow=())
    # Thus each property plot is written to a separate page of a single *pdf outfile
    pdf(file=paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByMetadata.pdf"), width = DefaultParameters$BaseSizeSinglePlotPdf, height = DefaultParameters$BaseSizeSinglePlotPdf)
    for (property in colnames(ExtraCellProperties)) {
      print(DimPlot(object = seurat.object.f, reduction = dim_red_method, group.by = property, combine = T, legend = "none") + ggtitle(property))
    }
    dev.off()
    
    StopWatchEnd$DimRedPlotsColuredByOptC$dim_red_method  <- Sys.time()
    
  }
  
  ####################################
  ### Colour dimension reduction plots showing each requested gene
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot showing each requested gene ***\n"))
  
  ### To program layout() for more than 3 genes in multiple rows
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("No selected genes for dimension reduction plots")
  }else{
    
    StopWatchStart$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
    
    ListOfGenesForDimRedPlots<-unlist(strsplit(ListGenes, ","))
    if (length(ListOfGenesForDimRedPlots) <= 4) {
      pdfWidth  <- (length(ListOfGenesForDimRedPlots) * DefaultParameters$BaseSizeMultiplePlotPdfWidth)
      pdfHeight <- DefaultParameters$BaseSizeMultiplePlotPdfHeight
      nColFeaturePlot <- length(ListOfGenesForDimRedPlots)
    }else{
      pdfWidth  <- 4 * DefaultParameters$BaseSizeMultiplePlotPdfWidth
      pdfHeight <- (as.integer(length(ListOfGenesForDimRedPlots) / 4) + 1) * DefaultParameters$BaseSizeMultiplePlotPdfHeight
      nColFeaturePlot <- 4
    }
    
    pdf(file=paste0(Tempdir, "/SELECTED_GENE_DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySelectedGenes.pdf"), width=pdfWidth, height=pdfHeight)
    print(FeaturePlot(object = seurat.object.f, ncol = nColFeaturePlot, features = c(ListOfGenesForDimRedPlots), cols = c("lightgrey", "blue"), reduction = dim_red_method, order = T))
    dev.off()
    
    StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
    
  }
}

####################################
### Finding differentially expressed genes for each cell cluster
####################################
writeLines("\n*** Finding differentially expressed genes for each cell cluster ***\n")

print(paste0("NumberOfClusters=", NumberOfClusters))

StopWatchStart$FindDiffMarkers  <- Sys.time()

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.f@meta.data))

seurat.object.markers <- FindAllMarkers(object = seurat.object.f, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)

SimplifiedDiffExprGenes.df <- seurat.object.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC")]
write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles,".", ProgramOutdir, "_MarkersPerCluster.tsv"), row.names = F, sep="\t", quote = F)

### Get top-2 genes sorted by cluster, then by avg_logFC
top_genes_by_cluster_for_tsne<-(seurat.object.markers %>% group_by(cluster) %>% top_n(DefaultParameters$NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers[["cluster"]]))

StopWatchEnd$FindDiffMarkers  <- Sys.time()

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  top_6_genes_by_cluster<-(seurat.object.markers %>% group_by(cluster) %>% top_n(6, avg_logFC))
  markers_file <- top_6_genes_by_cluster[,c("gene","cluster","p_val","avg_logFC")]
  write.table(markers_file, paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_markers/","TopTwoMarkersPerCluster.tsv"),row.names = F,sep="\t",quote = F)
} 

####################################
### Saving the R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving the R object ***\n")
  
  StopWatchStart$SaveRDS  <- Sys.time()
  
  OutfileRDS<-paste0(Tempdir, "/R_OBJECTS/", PrefixOutfiles,".", ProgramOutdir, "_object.rds")
  saveRDS(seurat.object.f, file = OutfileRDS)
  
  StopWatchEnd$SaveRDS  <- Sys.time()
  
}else{
  
  writeLines("\n*** Not saving R object ***\n")
  
}

####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes ***\n")

StopWatchStart$DiffMarkerViolinPlots  <- Sys.time()

NumberOfPanesForFeaturesPlot<-(NumberOfClusters*DefaultParameters$NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5

### Here need to program a better way to control the y-axis labels
pdf(file=paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_VlnPlot_CountsLog10_AfterClusters.pdf"), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), slot = 'counts', log = T, adjust = 1, pt.size = 0.5))
dev.off()

pdf(file=paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_VlnPlot_Norm_AfterClusters.pdf"), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), adjust = 1, pt.size = 0.5))
dev.off()

StopWatchEnd$DiffMarkerViolinPlots  <- Sys.time()

####################################
### Plot dimension reductions showing each cluster top genes
####################################
writeLines("\n*** Plot dimension reductions showing each cluster top genes ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {
  
  ####################################
  ### Dimension reduction plots showing each cluster top genes
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by each cluster top genes ***\n"))
  
  pdfWidth  <- 4 * DefaultParameters$BaseSizeMultiplePlotPdfWidth
  pdfHeight <- NumberOfClusters * DefaultParameters$BaseSizeMultiplePlotPdfHeight / 2
  pdf(file=paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS/", PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_EachClusterTop2DEGs.pdf"), width=pdfWidth, height=pdfHeight)
  print(FeaturePlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), cols = c("lightgrey", "blue"), reduction = dim_red_method, order = T))
  dev.off()
  
}

####################################
### Create summary plots outfile
####################################
writeLines("\n*** Create summary plots outfile ***\n")

### Make this option "N" if library(staplr) or pdftk are not available
SummaryPlots <- "Y"

StopWatchStart$SummaryPlots  <- Sys.time()

if (regexpr("^y$", SummaryPlots, ignore.case = T)[1] == 1) {
  SummaryPlotsPdf<-paste0(Tempdir, "/SUMMARY_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_summary_plots.pdf")
  ListOfPdfFilesToMerge<-c(paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_QC_VlnPlot.pdf"),
                           paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_UMAPPlot_ColourByCellClusters.pdf"),
                           paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_TSNEPlot_ColourByCellClusters.pdf"),
                           paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_VlnPlot_CountsLog10_AfterClusters.pdf"),
                           paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_UMAPPlot_EachClusterTop2DEGs.pdf"),
                           paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_TSNEPlot_EachClusterTop2DEGs.pdf")
  )
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("Create summary file")
  }else{
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste0(Tempdir, "/", PrefixOutfiles,".", ProgramOutdir, "_UMAPPlot_ColourBySelectedGenes.pdf"))
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste0(Tempdir, "/", PrefixOutfiles,".", ProgramOutdir, "_TSNEPlot_ColourBySelectedGenes.pdf"))
    print("Create summary file including dimension reduction plots for selected genes")
  }
  
  ### Stappling *pdf's
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge, output_filepath = SummaryPlotsPdf)
}

StopWatchEnd$SummaryPlots  <- Sys.time()

####################################
### Transform select *pdf files into *png
####################################
writeLines("\n*** Transform select *pdf files into *png ***\n")

StopWatchStart$TransformPdfToPng  <- Sys.time()

ListOfPdfFilesToPng<-c(paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_PCElbowPlot.pdf")
)

sapply(ListOfPdfFilesToPng,FUN=function(eachFile) {
  inpdf  <- eachFile
  outpng <- gsub(".pdf$",".png", x= inpdf, ignore.case = T, perl = T)
  CommandConvert <- paste0("convert -density 150 ", inpdf, " -quality 300 ", outpng)
  system(command = CommandConvert, input = NULL, wait = T)
})

StopWatchEnd$TransformPdfToPng  <- Sys.time()

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

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

### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(Tempdir)
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    if (DirName == "FILTERED_DATA_MATRICES" | DirName == "UNFILTERED_DATA_MATRICES") {
      sapply(list.dirs(TempdirWithData, full.names = F, recursive = F), FUN=function(SubDirName) {
        OutdirFinal <- gsub(pattern = Tempdir, replacement =  paste0(Outdir, "/", ProgramOutdir), x = paste0(TempdirWithData, "/", SubDirName))
        dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
        sapply(list.files(paste0(TempdirWithData, "/", SubDirName), pattern = ".gz", full.names = F), FUN=function(EachFileName) {
          file.copy(from=paste0(TempdirWithData, "/", SubDirName, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
          file.remove(paste0(TempdirWithData, "/", SubDirName, "/", EachFileName))
        })
      })
    }else{
      OutdirFinal <- paste0(Outdir, "/", ProgramOutdir, "/", DirName)
      print(OutdirFinal)
      dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
      sapply(list.files(TempdirWithData, pattern = paste0("^", PrefixOutfiles, ".", ProgramOutdir), full.names = F), FUN=function(EachFileName) {
        file.copy(from=paste0(TempdirWithData, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
        file.remove(paste0(TempdirWithData, "/", EachFileName))
      })
    }
  })
}

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
OutfileCPUusage <- gsub(x = OutfileCPUusage, pattern = Tempdir, replacement = Outdir)
writeLines(paste0("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report"))

quit()
