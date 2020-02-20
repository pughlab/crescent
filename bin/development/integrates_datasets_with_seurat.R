####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/v3.1/integration.html
###
### THINGS TO DO:
### 1) See Run_Seurat_v3.R for a list of functions to be implemented
### 2) Clean up outfiles, keep only those needed for front end and publication style plots
###    For example, saving the `seurat.object.integrated.sa` may not be necessary because the general `seurat.object.integrated.sa` includes the same info
### 3) In:
###    `for (dataset in rownames(InputsTable)) {
###       seurat.object.integrated$sample_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$sample_type)
###     }`
###   Need to avoid that the gsub replaces partial strings, for example 'SMTR03' will be replaced by list_DatasetToType[[SMTR03]] in both
###   'SMTR03' iself and in 'SMTR03t1_NonRad'. The second case is undesired.
###   Generating files like:
###   SMTR_res1.SEURAT_GlobalClustering_Rad_unknownt1_NonRad_TSNEPlot_ColourByCellClusters.pdf
###   When it should be:
###   SMTR_res1.SEURAT_GlobalClustering_Rad_unknown_TSNEPlot_ColourByCellClusters.pdf
####################################

####################################
### COMMENTS ON DIFFERENT SUB-VERSIONS OF SEURAT v3
### Parameter -v allows to use either sub-version of Seurat v3.1.1 or v3.0.3.923
### In principle users should use Seurat v3.1.1.
### Javier Diaz uses v3.0.3.923 for reproducibility of a study since the two versions produce different solutions with the exact same inputs and parameters
### In future only v3.1.1 or higher will be used
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
### COMMENTS ON USING REFERENCE DATASETS VS. ALL PAIRWISE COMPARISONS TO FIND ANCHORS
### Using 4 datasets with 3K to 5K compared runs using either `-k NA` or `-k 1,2`
### The main cimputing time was in step IntegrateData(). Whereas FindIntegrationAnchors() was similar:
### Step/time(minutes)      Using_-k_NA   Using_-k_1,2
### FindIntegrationAnchors	3.336035	    3.118804
### IntegrateData           2.968116	    1.666707
### In both cases used `-u MAX -b 10000` in a 3.1-GHz Intel Core i5 CPU with 2 cores and 16 GB RAM
####################################

####################################
### Required libraries
####################################
### Package 'Seurat' version 3 is needed to run QC, differential gene expression and clustering analyses
### Seurat's latest stable version can be installed like:
### install.packages('Seurat')
###
### Or the development version can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
###
### Or a specific older version (e.g. v3.0.3.9023) can be installed from GitHub like:
### `PathForV3_0_3_9023Libs<-paste(.libPaths(), "/Seurat_V3_0_3_9023", sep = "", collapse = "")`
### `dir.create(file.path(PathForV3_0_3_9023Libs), showWarnings = F, recursive = T)``
### `devtools::install_github("satijalab/seurat", ref = "556e598", lib = PathForV3_0_3_9023Libs)`
### ### To check that v3.0.3.9023 is installed
### `library("Seurat",lib.loc = PathForV3_0_3_9023Libs)``
### `packageVersion("Seurat")`
### ### should return:
### [1] ‘3.0.3.9023’
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # (CRAN) to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(cluster))      # (CRAN) to cluster/sort the sample correlations
####################################

####################################
### Required external packages
####################################
### 'UMAP' can be installed using `pip install umap-learn`
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
  make_option(c("-i", "--inputs_list"), default="NA",
              help="Path/name to a <tab> delimited file with the list of dataset IDs, infile path/names, and dataset type (e.g. batch) like:
                dataset1_id_d1  /path_to/d1  d1_type  d1_format  d1_mito_fraction  d1_ribo_fraction  d1_number_of_genes d1_number_of_reads
                dataset2_id_d2  /path_to/d2  d2_type  d2_format  d2_mito_fraction  d2_ribo_fraction  d2_number_of_genes d2_number_of_reads
                dataset3_id_d3  /path_to/d3  d3_type  d3_format  d3_mito_fraction  d3_ribo_fraction  d3_number_of_genes d3_number_of_reads
                ...etc
                Default = 'No default. It's mandatory to specify this parameter'

                Notes:
                (1)
                The order of the list of datasets in --inputs_list influences the results,
                including number of clusters, t-SNE/UMAP and differentially expressed genes
                
                (2)
                Column 4 indicates the dataset format. It must be in either 'MTX' or 'DGE':
                a) the path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files; or
                b) the path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns
                Notes: The 'MTX' files can be the outputs from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                       Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                
                (3)    
                Column 5 indicates the min,max cutoff values for the fraction of mitochondrial gene counts in each cell. For example:
                a) for whole cell scRNA-seq use '0,0.2'
                b) for Nuc-seq use '0,0.05'
                c) for negative value counts (e.g. if using TPM in log scale) refer negative values with an 'n', like 'n1,0.5'

                (4)    
                Column 6 indicates the min,max cutoff values for the fraction of ribosomal protein gene counts in each cell.
                For example: '0,0.75'
                c) for negative value counts (e.g. if using TPM in log scale) refer negative values with an 'n', like 'n1,0.5'
                
                (5)
                Column 7 indicates the min,max cutoff values for the number of unique genes measured in each cell. For example:
                '50,8000'
                
                (6)
                Column 8 indicates the min,max cutoff values for the number of reads measured in each cell. For example:
                '1,80000'
                
                (7)
                Datasets will be normalized using SCTransform and three levels of integration will be used for clustering:
                1) cluster cells from each dataset
                2) integrate cells by 'dataset_type' (column 3) and cluster them
                3) integrate cells from all datasets and cluster them"),
  
  #
  make_option(c("-k", "--reference_datasets"), default="NA",
              help="<comma> delimited number of row(s) in --inputs_list of datasets to be used as reference(s) for integration
                Or type 'NA' to run all-vs-all dataset pairwise comparisons (this is more time and memory consuming than using references).
                If references are used, then references will be integrated and anchors identified between non-reference and reference datasets,
                but not anchors will be identified between non-reference datasets (this saves time and memory)
                Default = 'NA'"),
  #
  make_option(c("-j", "--inputs_remove_barcodes"), default="NA",
              help="Path/name to a <tab> delimited list of barcodes to be removed from analysis, like:
                dataset1_id_d1  AAACCTGAGCTCCCAG
                dataset2_id_d2  AAACCTGTCACCATAG
                dataset3_id_d3  AAACCTGTCAGCTTAG
                Note: the barcode id must include the dataset ID
                Or type 'NA' to include all barcodes
                Default = 'NA'"),
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
  make_option(c("-c", "--infile_colour_dim_red_plots"), default="NA",
              help="A <tab> delimited table of barcodes and discrete properties to colour the t-SNE, like:
                Barcode                             CellClass    InOtherDatasets
                dataset1_id_d1_AAACCTGAGCGGCTTC-1   1            yes
                dataset2_id_d2_AAACCTGAGTCGAGTG-1   1            no
                dataset3_id_d3_AAACCTGCAAAGGAAG-1   2            yes
                Note: the barcode id must include the dataset ID
                Default = 'NA' (i.e. no --infile_colour_dim_red_plots is provided)"),
  #
  make_option(c("-g", "--list_genes"), default="NA",
              help="A <comma> delimited list of gene identifiers whose expression will be mapped into the t-SNE plots
                Default = 'NA' (no --list_genes are provided)"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file:
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                Default = '10'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh
                Default = '0.01'"),
  #
  make_option(c("-f", "--diff_gene_expr_comparisons"), default="1",
              help="One or more of the following options. If more than one, pass the choice <comma> delimited, e.g. '1,2,3'
                '0' = no differentially expressed genes are computed
                '1' = using global cell clusers, compares each cell cluster vs. the rest of cells
                '2' = using global cell clusers, compares each cell cluster from one sample vs. the same cluster from other samples
                '3' = using global cell clusers, for each sample type, compares each cell cluster vs. the rest of cells
                '4' = using global cell clusers, for each sample type, compares each cell cluster vs. the same cluster from other sample types
                '5' = using each sample re-clustered, compares each cell cluster vs. the rest of cells
                '6' = using each sample type re-clustered, compares each cell cluster vs. the rest of cells
                Default = '1'"),
  #
  make_option(c("-u", "--number_cores"), default="MAX",
              help="Indicate the number of cores to use for parellelization (e.g. '4') or type 'MAX' to determine and use all available cores in the system
                Default = 'MAX'"),
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
  make_option(c("-v", "--seurat_version"), default="1",
              help="Indicates the sub-version of Seurat to use:
                '1' to use the default version (likely v3.1.1)
                '2' to use v3.0.3.923
                Default = '1'"),
  #
  make_option(c("-b", "--max_global_variables"), default="4000",
              help="Indicates maximum allowed total size (in bytes) of global variables identified.
                 Used by library(future) to prevent too large exports
                 Default = '4000' for 4000 MiB")
)

opt <- parse_args(OptionParser(option_list=option_list))

InputsList              <- opt$inputs_list
ReferenceDatasets       <- opt$reference_datasets
InfileRemoveBarcodes    <- opt$inputs_remove_barcodes
Resolution              <- as.numeric(opt$resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
InfileColourDimRedPlots <- opt$infile_colour_dim_red_plots
ListGenes               <- opt$list_genes
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
ThreshReturn            <- as.numeric(opt$return_threshold)
DiffGeneExprComparisons <- opt$diff_gene_expr_comparisons
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl
SeuratVersion           <- as.numeric(opt$seurat_version)
MaxGlobalVariables      <- as.numeric(opt$max_global_variables)

####################################
### Load Seurat
####################################
writeLines("\n*** Create outdirs ***\n")

if (SeuratVersion == 1) {
  suppressPackageStartupMessages(library("Seurat"))
} else if (SeuratVersion == 2) {
  PathForV3_0_3_9023Libs<-paste(.libPaths(), "/Seurat_V3_0_3_9023", sep = "", collapse = "") ### PathForV3_0_3_9023Libs must be the path of the sub-version installation
  suppressPackageStartupMessages(library("Seurat", lib.loc = PathForV3_0_3_9023Libs)) 
}else{
  stop("Couldn't determine Seurat version to use")
}

####################################
### Define outdirs and CWL parameters
####################################
writeLines("\n*** Create outdirs ***\n")

ProgramOutdir <- "SEURAT"

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  ### Outfiles will be stored into `ProgramOutdir` directory
  #PrefixOutfiles <- "cwl_run" 
  PrefixOutfiles  <- opt$prefix_outfiles
  Tempdir         <- ProgramOutdir
  dir.create(file.path(Tempdir), showWarnings = F) ## Note Tempdir will be the final out-directory as well
}else{
  PrefixOutfiles <- c(paste(PrefixOutfiles,"_res",Resolution,sep=""))
  ## Using `Tempdir` for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished
  ## Then at the end of the script they'll be moved into `Outdir/ProgramOutdir`
  Tempdir        <- "~/temp" 
  #
  CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
  UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
  #
  Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
  Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
  Outdir<-gsub("/+", "/", Outdir, perl = T)
  Tempdir<-gsub("/+", "/", Tempdir, perl = T)
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
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript integrates_datasets_with_seurat.R -h\n\n", sep=""))
}

cat("Using ", NumbCoresToUse, "cores")

plan(strategy = "multicore", workers = NumbCoresToUse)

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4Gb of RAM (4000*1024^2), use: `options(future.globals.maxSize = 4000 * 1024^2)`
options(future.globals.maxSize = MaxGlobalVariables * 1024^2)

####################################
### Define default parameters
####################################
### Some of these default parameters are provided by Seurat developers,
### others are tailored according to clusters/t-SNE granularity

RequestedDiffGeneExprComparisons = unlist(strsplit(DiffGeneExprComparisons, ","))

DefaultParameters <- list(
  ### Parameters for QC plots
  CellPropertiesToQC = c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction"),

  ### Parameters for Seurat filters
  MinCells = 3,

  ### Parameters for Cluster Biomarkers
  FindAllMarkers.MinPct     =  0.25,
  FindAllMarkers.ThreshUse  =  0.25,

  ### Parameters for dimmension reduction plots
  BaseSizeSinglePlotPdf  = 7,
  BaseSizeSinglePlotPng  = 480,
  BaseSizeMultiplePlotPdfWidth  = 3.7,
  BaseSizeMultiplePlotPdfHeight = 3,
  MaxNumbLabelsPerRowInLegend   = 4,
  
  ### Parameters for datasets integration
  IntegrationNFeatures = 3000

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

ListMandatory<-list("infiles_list", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Load --inputs_list
####################################
writeLines("\n*** Load --inputs_list ***\n")

InputsList<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), InputsList)
InputsTable<-read.table(InputsList, header = F, row.names = 1, stringsAsFactors = F)
colnames(InputsTable)<-c("PathToDataset","DatasetType","DatasetFormat","MitoFrac","RiboFrac","NGenes","NReads")

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

SeuratObjectsFiltered   <-list()
SeuratObjectsUnfiltered <-list()
DatasetIds              <-list()
list_DatasetToType      <-list()
list_TypeToDatasets     <-list()
list_DatasetToFormat    <-list()
list_MinMitoFrac        <-list()
list_MaxMitoFrac        <-list()
list_MinRiboFrac        <-list()
list_MaxRiboFrac        <-list()
list_MinNGenes          <-list()
list_MaxNGenes          <-list()
list_MinNReads          <-list()
list_MaxNReads          <-list()

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  Dataset.SO <-paste(dataset, ".so",    sep = "", collapse = "")
  
  PathToDataset <- InputsTable[dataset,"PathToDataset"]
  DatasetType   <- InputsTable[dataset,"DatasetType"]
  DatasetFormat <- InputsTable[dataset,"DatasetFormat"]
  ListMitoFrac  <- InputsTable[dataset,"MitoFrac"]
  ListRiboFrac  <- InputsTable[dataset,"RiboFrac"]
  ListNGenes    <- InputsTable[dataset,"NGenes"]
  ListNReads    <- InputsTable[dataset,"NReads"]
  
  PathToDataset <- gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), PathToDataset)
  #
  ListMitoFrac = unlist(strsplit(ListMitoFrac,  ","))
  MinMitoFrac  = as.numeric(ListMitoFrac[1])
  MaxMitoFrac  = as.numeric(ListMitoFrac[2])
  #
  ListRiboFrac = unlist(strsplit(ListRiboFrac,  ","))
  MinRiboFrac  = as.numeric(ListRiboFrac[1])
  MaxRiboFrac  = as.numeric(ListRiboFrac[2])
  #
  ListNGenes = unlist(strsplit(ListNGenes, ","))
  MinNGenes  = as.numeric(ListNGenes[1])
  MaxNGenes  = as.numeric(ListNGenes[2])
  #
  ListNReads = unlist(strsplit(ListNReads, ","))
  MinNReads  = as.numeric(ListNReads[1])
  MaxNReads  = as.numeric(ListNReads[2])
  
  list_DatasetToType[[dataset]]      <- DatasetType
  list_DatasetToFormat[[dataset]]    <- DatasetFormat
  list_TypeToDatasets[[DatasetType]] <- append(list_TypeToDatasets[[DatasetType]], dataset)
  list_MinMitoFrac[[dataset]]        <- MinMitoFrac
  list_MaxMitoFrac[[dataset]]        <- MaxMitoFrac
  list_MinRiboFrac[[dataset]]        <- MinRiboFrac
  list_MaxRiboFrac[[dataset]]        <- MaxRiboFrac
  list_MinNGenes[[dataset]]          <- MinNGenes
  list_MaxNGenes[[dataset]]          <- MaxNGenes
  list_MinNReads[[dataset]]          <- MinNReads
  list_MaxNReads[[dataset]]          <- MaxNReads
  
  if (regexpr("^MTX$|^DGE$", DatasetFormat, ignore.case = T, perl = T)[1] == 1) {

    ####################################
    ### Loading MTX or DGE infiles
    ####################################
  
    StopWatchStart$LoadScRNAseqData$dataset <- Sys.time()
  
    if (regexpr("^MTX$", DatasetFormat, ignore.case = T)[1] == 1) {
      writeLines(paste("\n*** Loading MTX infiles for ", dataset, " from: ", PathToDataset, " ***\n", sep = "", collapse = ""))
      expression_matrix <- Read10X(data.dir = PathToDataset)
    }else if (regexpr("^DGE$", DatasetFormat, ignore.case = T)[1] == 1) {
      print("Loading Digital Gene Expression matrix")
      writeLines(paste("\n*** Loading Digital Gene Expression matrix for ", dataset, " from: ", PathToDataset, " ***\n", sep = "", collapse = ""))
      ## Note `check.names = F` is needed for both `fread` and `data.frame`
      expression_matrix <- as.matrix(data.frame(fread(PathToDataset, check.names = F), row.names=1, check.names = F))
    }
    dim(expression_matrix)
  
    StopWatchEnd$LoadScRNAseqData$dataset  <- Sys.time()

    ####################################
    ### Create seurat object
    ####################################
    
    StopWatchStart$CreateSeuratObject$dataset  <- Sys.time()
    
    writeLines(paste("\n*** Create seurat object for ", dataset, " ***\n", sep = "", collapse = ""))
    seurat.object.u  <- CreateSeuratObject(counts = expression_matrix, min.cells = DefaultParameters$MinCells, min.features = list_MinNGenes[[dataset]], project = paste(PrefixOutfiles, "_", dataset, sep = "", collapse = ""))
    
    StopWatchEnd$CreateSeuratObject$dataset  <- Sys.time()
    
    ####################################
    ### Add dataset label
    ####################################
    writeLines(paste("\n*** Add dataset label for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$AddDatasetLabel$dataset  <- Sys.time()
    
    seurat.object.u[['dataset.label']] <- dataset
    
    StopWatchEnd$AddDatasetLabel$dataset  <- Sys.time()

    ####################################
    ### Add dataset type label
    ####################################
    writeLines(paste("\n*** Add dataset type label for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$AddDatasetLabel$dataset_type  <- Sys.time()
    
    seurat.object.u[['dataset_type.label']] <- list_DatasetToType[[dataset]]
    
    StopWatchEnd$AddDatasetLabel$dataset_type  <- Sys.time()

    ####################################
    ### Get mitochondrial genes
    ####################################
    writeLines(paste("\n*** Get  mitochondrial genes for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$GetMitoGenes$dataset  <- Sys.time()
    
    mitoRegExpressions<- paste(c("^MT-"), collapse = "|")
    mito.features <- grep(pattern = mitoRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)
    
    if (length(mito.features)[[1]] > 0) {
      mito.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['mito.fraction']] <- mito.fraction
    }else{
      mito.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['mito.fraction']] <- mito.fraction
    }
    
    StopWatchEnd$GetMitoGenes$dataset  <- Sys.time()
    
    ####################################
    ### Get ribosomal protein genes
    ####################################
    writeLines(paste("\n*** Get ribosomal protein genes for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$GetRiboGenes$dataset  <- Sys.time()
    
    riboRegExpressions<- paste(c("^MRPL", "^MRPS", "^RPL", "^RPS"),collapse = "|")
    ribo.features <- grep(pattern = riboRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)
    
    if (length(ribo.features)[[1]] > 0) {
      ribo.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['ribo.fraction']] <- ribo.fraction
    }else{
      ribo.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['ribo.fraction']] <- ribo.fraction
    }
    
    StopWatchEnd$GetRiboGenes$dataset  <- Sys.time()
    
    ####################################
    ### Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation
    ####################################
    writeLines(paste("\n*** Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$FilterCells$dataset  <- Sys.time()
    
    if (length(mito.features)[[1]] > 0) {
      seurat.object.f<-subset(x = seurat.object.u, subset = 
                                nFeature_RNA >= list_MinNGenes[[dataset]]
                              & nFeature_RNA <= list_MaxNGenes[[dataset]] 
                              & nCount_RNA   >= list_MinNReads[[dataset]]
                              & nCount_RNA   <= list_MaxNReads[[dataset]]
                              & mito.fraction >= list_MinMitoFrac[[dataset]]
                              & mito.fraction <= list_MaxMitoFrac[[dataset]]
                              & ribo.fraction >= list_MinRiboFrac[[dataset]]
                              & ribo.fraction <= list_MaxRiboFrac[[dataset]])

      ### Get list of barcodes excluded by nFeature_RNA, nCount_RNA, mito.fraction or ribo.fraction
      BarcodesExcludedByNFeature        <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA  >= list_MinNGenes[[dataset]] & nFeature_RNA <= list_MaxNGenes[[dataset]])))
      BarcodesExcludedByNReads          <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nCount_RNA    >= list_MinNReads[[dataset]] & nCount_RNA   <= list_MaxNReads[[dataset]])))
      BarcodesExcludedByMito            <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = mito.fraction >= list_MinMitoFrac[[dataset]] & mito.fraction <= list_MaxMitoFrac[[dataset]])))
      BarcodesExcludedByRibo            <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = ribo.fraction >= list_MinRiboFrac[[dataset]] & ribo.fraction <= list_MaxRiboFrac[[dataset]])))

    }else{
      seurat.object.f<-subset(x = seurat.object.u, subset = 
                                nFeature_RNA >= list_MinNGenes[[dataset]]
                              & nFeature_RNA <= list_MaxNGenes[[dataset]] 
                              & nCount_RNA   >= list_MinNReads[[dataset]]
                              & nCount_RNA   <= list_MaxNReads[[dataset]]
                              & ribo.fraction >= list_MinRiboFrac[[dataset]]
                              & ribo.fraction <= list_MaxRiboFrac[[dataset]])
      
      ### Get list of barcodes excluded by nFeature_RNA, nCount_RNA, mito.fraction or ribo.fraction
      BarcodesExcludedByNFeature        <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA  >= list_MinNGenes[[dataset]] & nFeature_RNA <= list_MaxNGenes[[dataset]])))
      BarcodesExcludedByNReads          <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nCount_RNA    >= list_MinNReads[[dataset]] & nCount_RNA   <= list_MaxNReads[[dataset]])))
      BarcodesExcludedByMito            <- 0
      BarcodesExcludedByRibo            <- setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = ribo.fraction >= list_MinRiboFrac[[dataset]] & ribo.fraction <= list_MaxRiboFrac[[dataset]])))
    }
    NumberOfBarcodesExcludedByNFeature <- length(BarcodesExcludedByNFeature)
    NumberOfBarcodesExcludedByNReads   <- length(BarcodesExcludedByNReads)
    NumberOfBarcodesExcludedByMito     <- length(BarcodesExcludedByMito)
    NumberOfBarcodesExcludedByRibo     <- length(BarcodesExcludedByRibo)

    StopWatchEnd$FilterCells$dataset  <- Sys.time()
       
    ### Just reporting the summary of the UNfiltered and filtered objects
    print(paste(dataset, "  unfiltered", sep = "", collapse = ""))
    print(seurat.object.u)
    print(paste(dataset, "  filtered", sep = "", collapse = ""))
    print(seurat.object.f)
    
    ####################################
    ### QC EDA violin plots
    ####################################
    writeLines(paste("\n*** QC EDA violin plots for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$QCviolinplots$dataset  <- Sys.time()
    
    ### Get unfiltered data QC statistics
    nFeature_RNA.u.df   <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
    nCount_RNA.u.df     <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
    mito.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$mito.fraction, mito.fraction = 1)
    ribo.fraction.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$ribo.fraction, ribo.fraction = 1)
    #
    mean_nFeature_RNAStats.u   <- round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
    median_nFeature_RNAStats.u <- round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)
    mean_nCount_RNAStats.u     <- round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0)
    median_nCount_RNAStats.u   <- round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)
    mean_mito.fraction.u       <- round(mean(seurat.object.u@meta.data[,"mito.fraction"]),3)
    median_mito.fraction.u     <- round(median(seurat.object.u@meta.data[,"mito.fraction"]),3)
    mean_ribo.fraction.u       <- round(mean(seurat.object.u@meta.data[,"ribo.fraction"]),3)
    median_ribo.fraction.u     <- round(median(seurat.object.u@meta.data[,"ribo.fraction"]),3)
    #
    nFeature_RNAStats.u <-paste(c("mean = ", mean_nFeature_RNAStats.u, "\n", "median = ", median_nFeature_RNAStats.u), sep = "", collapse="")
    nCount_RNAStats.u   <-paste(c("mean = ", mean_nCount_RNAStats.u,   "\n", "median = ", median_nCount_RNAStats.u),   sep = "", collapse="")
    mito.fraction.u     <-paste(c("mean = ", mean_mito.fraction.u,     "\n", "median = ", median_mito.fraction.u),     sep = "", collapse="")
    ribo.fraction.u     <-paste(c("mean = ", mean_ribo.fraction.u,     "\n", "median = ", median_ribo.fraction.u),     sep = "", collapse="")
    
    ### Get filtered data QC statistics
    nFeature_RNA.f.df   <-data.frame(Expression_level = seurat.object.f@meta.data$nFeature_RNA, nGenes = 2)
    nCount_RNA.f.df     <-data.frame(Expression_level = seurat.object.f@meta.data$nCount_RNA,   nCount_RNA = 2)
    mito.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$mito.fraction, mito.fraction = 2)
    ribo.fraction.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$ribo.fraction, ribo.fraction = 2)
    #
    mean_nFeature_RNAStats.f   <- round(mean(seurat.object.f@meta.data[,"nFeature_RNA"]),0)
    median_nFeature_RNAStats.f <- round(median(seurat.object.f@meta.data[,"nFeature_RNA"]),0)
    mean_nCount_RNAStats.f     <- round(mean(seurat.object.f@meta.data[,"nCount_RNA"]),0)
    median_nCount_RNAStats.f   <- round(median(seurat.object.f@meta.data[,"nCount_RNA"]),0)
    mean_mito.fraction.f       <- round(mean(seurat.object.f@meta.data[,"mito.fraction"]),3)
    median_mito.fraction.f     <- round(median(seurat.object.f@meta.data[,"mito.fraction"]),3)
    mean_ribo.fraction.f       <- round(mean(seurat.object.f@meta.data[,"ribo.fraction"]),3)
    median_ribo.fraction.f     <- round(median(seurat.object.f@meta.data[,"ribo.fraction"]),3)
    #
    nFeature_RNAStats.f <-paste(c("mean = ", mean_nFeature_RNAStats.f, "\n", "median = ", median_nFeature_RNAStats.f), sep = "", collapse="")
    nCount_RNAStats.f   <-paste(c("mean = ", mean_nCount_RNAStats.f,   "\n", "median = ", median_nCount_RNAStats.f),   sep = "", collapse="")
    mito.fraction.f     <-paste(c("mean = ", mean_mito.fraction.f,     "\n", "median = ", median_mito.fraction.f),     sep = "", collapse="")
    ribo.fraction.f     <-paste(c("mean = ", mean_ribo.fraction.f,     "\n", "median = ", median_ribo.fraction.f),     sep = "", collapse="")
    
    ### Put QC statistics together
    nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
    nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
    mito.fraction.m.df  <-data.frame(rbind(mito.fraction.u.df,mito.fraction.f.df))
    ribo.fraction.m.df  <-data.frame(rbind(ribo.fraction.u.df,ribo.fraction.f.df))
    NumberOfCells<- list()
    NumberOfCells[["unfiltered"]] <- nrow(seurat.object.u@meta.data)
    NumberOfCells[["filtered"]]   <- nrow(seurat.object.f@meta.data)
    LabelUnfiltered    <-paste("Before filters: No. of cells = ", NumberOfCells[["unfiltered"]], sep ="", collapse = "")
    LabelFiltered      <-paste("After filters:  No. of cells = ", NumberOfCells[["filtered"]],   sep ="", collapse = "")

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
      labs(x=paste("No. of genes", "\n", "Filter: ", "min=", list_MinNGenes[[dataset]], " max=", list_MaxNGenes[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNFeature, sep = "", collapse = "")) +
      annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])
    
    nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
      geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x=paste("No. of reads", "\n", "Filter: ", "min=", list_MinNReads[[dataset]], " max=", list_MaxNReads[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNReads, sep = "", collapse = "")) +
      annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])

    mito.fraction.plot<-ggplot(data=mito.fraction.m.df, aes(x = factor(mito.fraction), y = Expression_level)) +
      geom_violin(aes(fill = factor(mito.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x="Mitochondrial genes (fraction)") +
      labs(x=paste("Mitochondrial genes", "\n", "Filter: ", "min=", list_MinMitoFrac[[dataset]], " max=", list_MaxMitoFrac[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByMito, sep = "", collapse = "")) +
      annotate("text", x = 1 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.f, col = ColoursQCViolinPlots[[2]])
    
    ribo.fraction.plot<-ggplot(data=ribo.fraction.m.df, aes(x = factor(ribo.fraction), y = Expression_level)) +
      geom_violin(aes(fill = factor(ribo.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x=paste("Ribosomal protein genes (fraction)", "\n", "Filter: ", "min=", list_MinRiboFrac[[dataset]], " max=", list_MaxRiboFrac[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByRibo, sep = "", collapse = "")) +
      annotate("text", x = 1 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.f, col = ColoursQCViolinPlots[[2]])

    bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, mito.fraction.plot, ribo.fraction.plot, ncol = 4)
    
    ### Create a *pdf file with the violin ggplot's
    
    VlnPlotPdf<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_QC_VlnPlot.pdf", sep="", collapse = )
    pdf(file=VlnPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 1.7, height = DefaultParameters$BaseSizeSinglePlotPdf)
    print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
    dev.off()
    
    StopWatchEnd$QCviolinplots$dataset  <- Sys.time()
    
    ####################################
    ### Feature-vs-feature scatter plot
    ####################################
    writeLines(paste("\n*** Feature-vs-feature scatter plot for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$FeatureVsFeatureplot$dataset  <- Sys.time()
    
    UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                                  nGene = seurat.object.u@meta.data$nFeature_RNA,
                                  mito.fraction = seurat.object.u@meta.data$mito.fraction,
                                  filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.f))
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = ColoursQCViolinPlots[[1]])
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = F, replacement = ColoursQCViolinPlots[[2]])
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = 4)
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = F, replacement = 16)
    
    FeatureVsFeaturePlotPdf<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_NumbReadsVsNumbGenesAndMito_ScatterPlot.pdf", sep="", collapse = "")
    pdf(file=FeatureVsFeaturePlotPdf, width = 10, height = 5)
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
    
    StopWatchEnd$FeatureVsFeatureplot$dataset  <- Sys.time()
    
    ####################################
    ### Write out filter details and number of filtered cells
    ####################################
    writeLines(paste("\n*** Write out filter details and number of filtered cells for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$OutTablesFilterDetailsAndFilteredCells$dataset  <- Sys.time()
    
    OutTableFilterDetails<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_FilterDetails.tsv", sep="", collapse = "")
    Headers<-paste("Step", "Filter_min", "Filter_max", "Mean_before_filter", "Median_before_filter", "Mean_after_filter", "Median_after_filter", "Excluded_cells", sep = "\t", collapse = "")
    write.table(Headers, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F)
    
    FilterDetails.nFeature_RNA  <- paste("nFeature_RNA", list_MinNGenes[[dataset]], list_MaxNGenes[[dataset]],
                                      mean_nFeature_RNAStats.u, median_nFeature_RNAStats.u, mean_nFeature_RNAStats.f, median_nFeature_RNAStats.f, 
                                      NumberOfBarcodesExcludedByNFeature, sep = "\t", collapse = "")
    FilterDetails.nCount_RNA    <- paste("nCount_RNA", list_MinNReads[[dataset]], list_MaxNReads[[dataset]], 
                                      mean_nCount_RNAStats.u, median_nCount_RNAStats.u, mean_nCount_RNAStats.f, median_nCount_RNAStats.f, 
                                      NumberOfBarcodesExcludedByNReads, sep = "\t", collapse = "")
    FilterDetails.mito.fraction <- paste("mito.fraction", list_MinMitoFrac[[dataset]], list_MaxMitoFrac[[dataset]],
                                      mean_mito.fraction.u, median_mito.fraction.u, mean_mito.fraction.f, median_mito.fraction.f, 
                                      NumberOfBarcodesExcludedByMito, sep = "\t", collapse = "")
    FilterDetails.ribo.fraction <- paste("ribo.fraction", list_MinRiboFrac[[dataset]], list_MaxRiboFrac[[dataset]],
                                      mean_ribo.fraction.u, median_ribo.fraction.u, mean_ribo.fraction.f, median_ribo.fraction.f, 
                                      NumberOfBarcodesExcludedByRibo, sep = "\t", collapse = "")
    
    write.table(FilterDetails.nFeature_RNA,  file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    write.table(FilterDetails.nCount_RNA,    file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    write.table(FilterDetails.mito.fraction, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    write.table(FilterDetails.ribo.fraction, file = OutTableFilterDetails, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    
    OutTableFilteredCells<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_NumberOfFilteredCells.tsv", sep="", collapse = "")
    write.table(paste("Number_of_cells_before_filters", NumberOfCells[["unfiltered"]], sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(paste("Number_of_cells_after_filters", NumberOfCells[["filtered"]],    sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F, append = T)

    StopWatchEnd$OutTablesFilterDetailsAndFilteredCells$dataset  <- Sys.time()

    ####################################
    ### Assign data to Datasets lists
    ####################################
    writeLines(paste("\n*** Assign data to Datasets lists: ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$AssignDataToDatasets  <- Sys.time()
    
    SeuratObjectsUnfiltered[[as.character(NumberOfDatasets)]]  <- seurat.object.u
    SeuratObjectsFiltered[[as.character(NumberOfDatasets)]]    <- seurat.object.f
    DatasetIds[[as.character(NumberOfDatasets)]]               <- dataset

    StopWatchEnd$AssignDataToDatasets  <- Sys.time()

    ####################################
    ### Write out QC data for each sample
    ####################################
    writeLines("\n*** Write out QC data for each sample ***\n")
    
    StopWatchStart$WriteOutQCData  <- Sys.time()
    
    Headers<-paste("Cell_barcode", paste(DefaultParameters$CellPropertiesToQC, sep = "", collapse = "\t") ,sep="\t")
    
    BarcodeIdsWithDatasetBeforeFilters <- unlist(x = strsplit(x = paste(dataset, colnames(SeuratObjectsUnfiltered[[NumberOfDatasets]]), sep = "_", collapse = "\n"), split = "\n"))
    OutfileQCMetadataBeforeFilters<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", "Before_filters_QC_metadata.tsv", sep = "", collapse = "")
    write.table(Headers, file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(BarcodeIdsWithDatasetBeforeFilters, SeuratObjectsUnfiltered[[NumberOfDatasets]]@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)

    BarcodeIdsWithDatasetAfterFilters <- unlist(x = strsplit(x = paste(dataset, colnames(SeuratObjectsFiltered[[NumberOfDatasets]]), sep = "_", collapse = "\n"), split = "\n"))
    OutfileQCMetadataAfterFilters<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", "After_filters_QC_metadata.tsv", sep = "", collapse = "")
    write.table(Headers, file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(BarcodeIdsWithDatasetAfterFilters, SeuratObjectsFiltered[[NumberOfDatasets]]@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    
    StopWatchEnd$WriteOutQCData  <- Sys.time()
    
    ####################################
    ### Remove the Unfiltered seurat object
    ####################################
    writeLines(paste("\n*** Remove the Unfiltered seurat object for ", dataset, " ***\n", sep = "", collapse = ""))
    
    rm(seurat.object.u)
    rm(UnfilteredData.df)
  
  }else{
    stop(paste("Unexpected type of input: ", DatasetType, "\n\nFor help type:\n\nRscript integrates_datasets_with_seurat.R -h\n\n", sep=""))
  }
}

####################################
### Get's number of rows for the legend of dimension reduction plots
####################################
if (NumberOfDatasets <= DefaultParameters$MaxNumbLabelsPerRowInLegend) {
  NumbRowsForLegendPerSample <- 1
}else{
  NumbRowsForLegendPerSample <- round((NumberOfDatasets / DefaultParameters$MaxNumbLabelsPerRowInLegend), digits = 0)
}

NumberOfDatasetsTypes <- length(unique(names(list_TypeToDatasets)))
if (NumberOfDatasetsTypes <= DefaultParameters$MaxNumbLabelsPerRowInLegend) {
  NumbRowsForLegendPerSampleType <- 1
}else{
  NumbRowsForLegendPerSampleType <- round((NumberOfDatasetsTypes / DefaultParameters$MaxNumbLabelsPerRowInLegend), digits = 0)
}

####################################
### Remove barcodes by parameter -j (if applicable)
####################################

if (regexpr("^NA$", InfileRemoveBarcodes , ignore.case = T)[1] == 1) {

    writeLines("\n*** No barcodes are removed by option -j ***\n")

}else{

  StopWatchStart$RemoveBarcodes <- Sys.time()
  
  writeLines("\n*** Removing barcodes by parameter -j ***\n")
  
  InfileAllBarcodesToRemove<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), InfileRemoveBarcodes)
  AllBarcodesToRemove.tab<-read.table(InfileAllBarcodesToRemove, header = F, row.names = NULL, stringsAsFactors = FALSE)
  if (ncol(AllBarcodesToRemove.tab) == 2) {
    colnames(AllBarcodesToRemove.tab) <- c("Dataset","Barcode")
  }else{
    stop(paste("Unexpected format in file:\n", InfileAllBarcodesToRemove, "\nfor parameter -j, 2 columns were expected but found ", ncol(AllBarcodesToRemove.tab), sep = "", collapse = ""))
  }
  
  for (SeuratObjNumb in c(1:NumberOfDatasets)) {
    print(paste("Removing cells from:", DatasetIds[[SeuratObjNumb]] , sep = "", collapse = ""))
    seurat.object.full <- SeuratObjectsFiltered[[SeuratObjNumb]]
    seurat.object.full
    DatasetId          <- DatasetIds[[SeuratObjNumb]]
    ThisDatasetBarcodesToRemove    <- subset(x=AllBarcodesToRemove.tab, subset = Dataset == DatasetId)[,"Barcode"]
    ThisDatasetBarcodesToKeep.log  <- !colnames(seurat.object.full) %in% ThisDatasetBarcodesToRemove
    seurat.object.subset           <- subset(seurat.object.full, cells = colnames(seurat.object.full[,ThisDatasetBarcodesToKeep.log]))
    seurat.object.subset
    SeuratObjectsFiltered[[SeuratObjNumb]] <- seurat.object.subset
    print(paste(DatasetId, paste("Before:", ncol(seurat.object.full), sep = "", collapse =""), paste("After:", ncol(seurat.object.subset), sep = "", collapse =""), sep = "  ", collapse = "\n"))
  }
  
  StopWatchEnd$RemoveBarcodes <- Sys.time()
  
}

####################################
### Merge Seurat objects
####################################
writeLines("\n*** Merge Seurat objects ***\n")

StopWatchStart$MergeSeuratObjectsFiltered  <- Sys.time()

FirstSeuratObject   <- SeuratObjectsFiltered[[1]]
FirstSampleId       <- DatasetIds[[1]]
RestOfSeuratObjectsFiltered <- SeuratObjectsFiltered[c(2:NumberOfDatasets)]
RestOfSamplesIds    <- unlist(DatasetIds[c(2:NumberOfDatasets)])

seurat.object.merged <- merge(FirstSeuratObject, y = RestOfSeuratObjectsFiltered, 
                              add.cell.ids = DatasetIds,
                              project = PrefixOutfiles)

dataset.label.metadata<-cbind.data.frame(sample=seurat.object.merged@meta.data$dataset.label)
rownames(dataset.label.metadata)<-colnames(seurat.object.merged)

seurat.object.merged.withmetadata <- CreateSeuratObject(seurat.object.merged@assays$RNA@counts, meta.data = dataset.label.metadata)
seurat.object.list <- SplitObject(seurat.object.merged.withmetadata, split.by = "sample")

StopWatchEnd$MergeSeuratObjectsFiltered  <- Sys.time()

####################################
### Running SCTransform
####################################
writeLines("\n*** Running SCTransform ***\n")

StopWatchStart$SCTransform  <- Sys.time()

for (i in 1:length(seurat.object.list)) {
  seurat.object.list[[i]] <- SCTransform(seurat.object.list[[i]], verbose = T)
}

StopWatchEnd$SCTransform  <- Sys.time()

####################################
### Integrating datasets
####################################
writeLines("\n*** Integrating datasets ***\n")

StopWatchStart$SelectIntegrationFeatures  <- Sys.time()

writeLines("\n*** Run SelectIntegrationFeatures() ***\n")
seurat.object.integratedfeatures <- SelectIntegrationFeatures(object.list = seurat.object.list, nfeatures = DefaultParameters$IntegrationNFeatures)

StopWatchEnd$SelectIntegrationFeatures  <- Sys.time()

StopWatchStart$PrepSCTIntegration  <- Sys.time()

writeLines("\n*** Run PrepSCTIntegration() ***\n")
seurat.object.list <- PrepSCTIntegration(object.list = seurat.object.list, anchor.features = seurat.object.integratedfeatures, verbose = T)

StopWatchEnd$PrepSCTIntegration  <- Sys.time()

StopWatchStart$FindIntegrationAnchors  <- Sys.time()

writeLines("\n*** Run FindIntegrationAnchors() ***\n")
if (regexpr("^NA$", ReferenceDatasets , ignore.case = T)[1] == 1) {
  writeLines("\n*** No reference datasets will be used. Finding anchors in all-vs-all dataset pairwise comparisons ***\n")
  seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, normalization.method = "SCT", anchor.features = seurat.object.integratedfeatures, verbose = T)

}else if (regexpr("[0-9]", ReferenceDatasets , ignore.case = T, perl = T)[1] == 1) {
  writeLines(paste("\n*** Dataset number(s): ", ReferenceDatasets, " will be used as reference(s) ***\n", sep = "", collapse = ""))
  seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, normalization.method = "SCT", anchor.features = seurat.object.integratedfeatures, reference = c(as.numeric(unlist(strsplit(ReferenceDatasets, ",")))), verbose = T)

}else{
  stop(paste("Unexpected format in --reference_datasets ", ReferenceDatasets, sep = "", collapse = ""))
}

StopWatchEnd$FindIntegrationAnchors  <- Sys.time()

StopWatchStart$IntegrateData  <- Sys.time()

writeLines("\n*** Run IntegrateData() ***\n")
seurat.object.integrated <- IntegrateData(anchorset = seurat.object.anchors, normalization.method = "SCT", verbose = T)

StopWatchEnd$IntegrateData  <- Sys.time()

####################################
### Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset (columns)
####################################
writeLines("\n*** Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset (columns) ***\n")

seurat.object.integrated.list <- SplitObject(seurat.object.integrated, split.by = "sample")

mat_for_correl_all_cells.df <- data.frame(row.names = rownames(seurat.object.integrated.list[[1]]@assays$SCT))
for (dataset in rownames(InputsTable)) {
  mat_for_correl_all_cells.df[[dataset]] <- rowSums(as.matrix(seurat.object.integrated.list[[dataset]]@assays$SCT[,]))
}

OutfileSCTransform <- paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_SCTransform_SamplesInColumns.tsv", sep = "", collapse = "")
Headers<-paste("SCTransform", paste(colnames(mat_for_correl_all_cells.df), sep = "\t", collapse = "\t"), sep = "\t", collapse = "")
write.table(Headers, file = OutfileSCTransform, row.names = F, col.names = F, sep="\t", quote = F)
write.table(mat_for_correl_all_cells.df,  file = OutfileSCTransform, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Saving the R object up to dataset integration
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving the R object up to dataset integration ***\n")
  
  StopWatchStart$SaveRDS  <- Sys.time()
  
  OutfileRDS<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_integrated_object_incl_integrated_datasets.rds", sep="")
  saveRDS(seurat.object.integrated, file = OutfileRDS)
  
  StopWatchEnd$SaveRDS  <- Sys.time()
  
}else{
  
  writeLines("\n*** Saving the R object up to dataset integration ***\n")
  
}

####################################
### Obtaining principal components
####################################
writeLines("\n*** Obtaining principal components ***\n")

StopWatchStart$RunPCA  <- Sys.time()

seurat.object.integrated <- RunPCA(seurat.object.integrated, verbose = F)

StopWatchEnd$RunPCA  <- Sys.time()

StopWatchStart$PCAPlots  <- Sys.time()

ForElbowPlot<-ElbowPlot(object = seurat.object.integrated, ndims = 50, reduction = "pca")
MaxYAxis<-as.integer(max(ForElbowPlot$data$stdev)+1)

pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_PCElbowPlot_Integrated.pdf", sep=""))
print(ForElbowPlot
      + scale_x_continuous(breaks =  seq(from = 0, to = 50, by=5))
      + geom_vline(xintercept = seq(from = 0, to = 50, by=5), linetype='dotted', col="red")
      + scale_y_continuous(breaks =  seq(from = 0, to = MaxYAxis, by=0.5))
      + geom_hline(yintercept = seq(from = 0, to = MaxYAxis, by=0.5), linetype='dotted', col="red")
)
dev.off()

StopWatchEnd$PCAPlots  <- Sys.time()

####################################
### Run dimension reductions using integrated data
####################################
writeLines("\n*** Run dimension reductions using integrated data ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {
  ####################################
  ### Run non-linear dimension reductions using integrated data
  ####################################
  writeLines(paste("\n*** Run ", DimensionReductionMethods[[dim_red_method]][["name"]], " ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimensionReduction$dim_red_method  <- Sys.time()
  
  ### NOTES:
  ### In RunTSNE: if the datasets is small user may get error:
  ### `Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large.`
  ### User can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
  
  seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, do.fast = T)
  
  StopWatchEnd$DimensionReduction$dim_red_method  <- Sys.time()
}

####################################
### Colour dimension reduction plots for each sample by QC (nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction)
####################################
writeLines("\n*** Colour dimension reduction plots for each sample by QC (nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction) ***\n")

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  
  ### Note: need to AddMetaData() mito.fraction and ribo.fraction from OutfileQCMetadata to generate FeaturePlot()
  ### because they are not inherited in seurat.object.integrated by IntegrateData()
  InfileQCMetadataAfterFilters<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", "After_filters_QC_metadata.tsv", sep = "", collapse = "")
  
  QCMetadata <- data.frame(read.table(InfileQCMetadataAfterFilters, header = T, row.names = 1))
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = QCMetadata)
  
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_sample") == T) {
    rm(seurat.object.each_sample)
  }
  seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_sample)
  
  for (dim_red_method in names(DimensionReductionMethods)) {
    pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByQC.pdf", sep=""), width = DefaultParameters$BaseSizeSinglePlotPdf * 1.5, height = DefaultParameters$BaseSizeSinglePlotPdf * 1.5)
    print(FeaturePlot(object = seurat.object.each_sample, label = F, order = T, features = DefaultParameters$CellPropertiesToQC, cols = c("lightgrey", "blue"), reduction = dim_red_method, ncol = 2, pt.size = 1.5))
    dev.off()
  }
}

####################################
### Globally cluster cells using integrated data
####################################
writeLines("\n*** Globally cluster cells using integrated data ***\n")

StopWatchStart$ClusterAllCells  <- Sys.time()

## Uncomment this once fixed
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.integrated <- FindNeighbors(object = seurat.object.integrated, dims = PcaDimsUse)
seurat.object.integrated <- FindClusters(object = seurat.object.integrated, resolution = Resolution)

StopWatchEnd$ClusterAllCells  <- Sys.time()

StopWatchStart$AllCellClusterTables  <- Sys.time()

CellNames<-rownames(seurat.object.integrated@meta.data)
ClusterIdent <-seurat.object.integrated@meta.data$seurat_clusters
NumberOfClusters<-length(unique(ClusterIdent))

Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)
#
OutfileNumbCellsPerCluster<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerCluster.tsv", sep="")
Headers<-paste("Cluster", "Number_of_cells" ,sep="\t", collapse = "")
write.table(Headers,file = OutfileNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F)
write.table(table(ClusterIdent),file = OutfileNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
StopWatchEnd$AllCellClusterTables  <- Sys.time()

####################################
### Saving global cell cluster identities
####################################
writeLines("\n*** Saving global cell cluster identities ***\n")

seurat.object.integrated$GlobalCellClusterIdentities <- Idents(object = seurat.object.integrated)

####################################
### Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset cluster (columns)
####################################
writeLines("\n*** Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset cluster (columns) ***\n")

seurat.object.integrated.list <- SplitObject(seurat.object.integrated, split.by = "sample")

mat_for_correl_each_cluster.df <- data.frame(row.names = rownames(seurat.object.integrated.list[[1]]@assays$SCT))
for (dataset in rownames(InputsTable)) {
  if (exists(x = "seurat.object.each_sample") == T) {
    rm(seurat.object.each_sample)
  }
  seurat.object.each_sample <- subset(x = seurat.object.integrated, subset = sample == dataset)
  
  for (cluster_number in sort(unique(seurat.object.integrated$GlobalCellClusterIdentities))) {
    dataset_cluster <- (paste(dataset, cluster_number, sep = "_c", collapse = ""))
    res <- try(subset(x = seurat.object.each_sample, subset = seurat_clusters == cluster_number), silent = TRUE)
    if (class(res) == "try-error") {
      # mat_for_correl_each_cluster.df[[dataset_cluster]] <- NA
      }else{
      if (exists(x = "seurat.object.each_sample.each_cluster") == T) {
        rm(seurat.object.each_sample.each_cluster)
      }
      seurat.object.each_sample.each_cluster <- subset(x = seurat.object.each_sample, subset = seurat_clusters == cluster_number)
      mat_for_correl_each_cluster.df[[dataset_cluster]] <- rowSums(as.matrix(seurat.object.each_sample.each_cluster@assays$SCT[,]))
    }
  }
}

OutfileSCTransform <- paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_SCTransform_SampleClustersInColumns.tsv", sep = "", collapse = "")
Headers<-paste("SCTransform", paste(colnames(mat_for_correl_each_cluster.df), sep = "\t", collapse = "\t"), sep = "\t", collapse = "")
write.table(Headers, file = OutfileSCTransform, row.names = F, col.names = F, sep="\t", quote = F)
write.table(mat_for_correl_each_cluster.df,  file = OutfileSCTransform, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Get average gene expression for each global cluster
####################################
writeLines("\n*** Get average gene expression for each global cluster ***\n")

StopWatchStart$AverageGeneExpression  <- Sys.time()

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)
#
OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllSamples_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllSamples_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllSamples_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("r", Resolution, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$AverageGeneExpression  <- Sys.time()

####################################
### Finding differentially expressed genes: using global cell clusers, compares each cell cluster vs. the rest of cells
####################################

if (1 %in% RequestedDiffGeneExprComparisons == T) {

  writeLines("\n*** Finding differentially expressed genes: using global cell clusers, compares each cell cluster vs. the rest of cells ***\n")
  
  print(paste("NumberOfClusters=", NumberOfClusters, sep = "", collapse = ""))
  
  StopWatchStart$FindDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  seurat.object.integrated.markers <- FindAllMarkers(object = seurat.object.integrated, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
  SimplifiedDiffExprGenes.df <- seurat.object.integrated.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC")]
  write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "MarkersPerCluster.tsv", sep=""), row.names = F, sep="\t", quote = F)
  
  StopWatchEnd$FindDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()

}

####################################
### FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data
####################################
writeLines("\n*** FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {

  ####################################
  ### Colour dimension reduction plot by global cell clusters using integrated data
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by global cell clusters ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimRedOPlotColourByCellCluster$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 10)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
  IntegratedDimRedPlotPdf<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir,  "_GlobalClustering_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$DimRedOPlotColourByCellCluster$dim_red_method  <- Sys.time()

  ####################################
  ### Write out dimension reduction plot coordinates using integrated data
  ####################################
  writeLines(paste("\n*** Write out ", DimensionReductionMethods[[dim_red_method]][["name"]], " coordinates ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  OutfileCoordinates<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_Coordinates.tsv", sep="", collapse = "")
  Headers<-paste("Barcode",paste(colnames(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
  write.table(Headers,file = OutfileCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by dataset using integrated data
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by dataset ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimRedPlotsByDataset$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("sample"), combine = F, reduction = dim_red_method)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = NumbRowsForLegendPerSample, byrow = T, override.aes = list(size = 2))))
  IntegratedDimRedPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByDataset.pdf", sep="")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  ### This is the same plot as above, but with only two dataset labels per row in the legend to make sure that all text is printed out
  plots <- DimPlot(seurat.object.integrated, group.by = c("sample"), combine = F, reduction = dim_red_method)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(ncol = 2, byrow = T, override.aes = list(size = 1))))
  IntegratedDimRedPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByDataset_FullLegend.pdf", sep="")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = (7 + (0.1 * NumberOfDatasets)))
  print(CombinePlots(plots))
  dev.off()

  StopWatchEnd$DimRedPlotsByDataset$dim_red_method  <- Sys.time()

  ####################################
  ### Colour dimension reduction plots by -infile_colour_dim_red_plots using integrated data
  ####################################
  
  if (regexpr("^NA$", InfileColourDimRedPlots, ignore.case = T)[1] == 1) {
    print("No extra barcode-attributes will be used for dimension reduction plots")
  }else{
    writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by -infile_colour_dim_red_plots ***\n", sep = "", collapse = ""))
    
    StopWatchStart$DimRedPlotsColuredByMetadata$dim_red_method  <- Sys.time()
    
    seurat.object.meta.data<-seurat.object.integrated@meta.data
    ExtraCellProperties <- data.frame(read.table(InfileColourDimRedPlots, header = T, row.names = 1))

    # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
    # so that barcodes from --infile_colour_dim_red_plots and --input can match each other
    rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
    seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = ExtraCellProperties)
    
    # Generating outfile
    # Note DimPlot() takes the entire current device (pdf)
    # even if using layout(matrix(...)) or  par(mfrow=())
    # Thus each property plot is written to a separate page of a single *pdf outfile
    for (property in colnames(ExtraCellProperties)) {
      IntegratedDimRedPlotPdf <- paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBy_", property, ".pdf", sep="")
      pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(ExtraCellProperties[,property])))))
      
      if ( (sum(ExtraCellProperties[,property] %in% 0:1 == T)) == (nrow(ExtraCellProperties)) ) { ## is binary
        CellsToHighlight <- rownames(ExtraCellProperties)[ExtraCellProperties[,property]==1]
        plots <- DimPlot(seurat.object.integrated, group.by = property, combine = F, reduction = dim_red_method, label = F, cells.highlight = CellsToHighlight)
        plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
        print(CombinePlots(plots))
        
      }else{
        plots <- DimPlot(seurat.object.integrated, group.by = property, combine = F, reduction = dim_red_method, label = T)
        plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top", legend.text.align = 0) + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
        print(CombinePlots(plots))
        
      }
      dev.off()
    }
    
    StopWatchEnd$DimRedPlotsColuredByMetadata$dim_red_method  <- Sys.time()
  }

  ####################################
  ### Colour dimension reduction plots showing each requested gene using integrated data
  ####################################
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot showing each requested gene ***\n", sep = "", collapse = ""))
  
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
    
    ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
    ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
    seurat.object.integrated.sa <- CreateSeuratObject(seurat.object.integrated@assays$RNA@data)
    seurat.object.integrated.sa@reductions$umap <- seurat.object.integrated@reductions$umap
    seurat.object.integrated.sa@reductions$tsne <- seurat.object.integrated@reductions$tsne
    
    pdf(file=paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
    print(FeaturePlot(object = seurat.object.integrated.sa, ncol = nColFeaturePlot, features = ListOfGenesForDimRedPlots, cols = c("lightgrey", "blue"),
                      reduction = dim_red_method, order = T, slot = "data", pt.size = 1.5, min.cutoff = "q0.1", max.cutoff = "q90"))
    dev.off()
    
    #### This object is mainly for the front end
    if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
      
      writeLines("\n*** Saving the R object with only RNA data ***\n")
      
      StopWatchStart$SaveRDSOnlyRNA$DimensionReductionMethods[[dim_red_method]][["name"]]  <- Sys.time()
      
      OutfileRDS<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "_integrated_object_sa_assays_RNA_data.rds", sep="")
      saveRDS(seurat.object.integrated.sa, file = OutfileRDS)
      
      StopWatchEnd$SaveRDSOnlyRNA$DimensionReductionMethods[[dim_red_method]][["name"]]  <- Sys.time()
      
    }

    rm(seurat.object.integrated.sa)
    
    StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
    
  }
}

####################################
### FOR EACH SAMPLE merge sample_id and GLOBAL CLUSTER identities to make sample-specific identities
### Then get DGE and colour dimension reduction plots
####################################
writeLines("\n*** FOR EACH SAMPLE merge sample_id and GLOBAL CLUSTER to make sample-specific identities. Then get DGE and colour dimension reduction plots  ***\n")

EachSampleGlobalCellClusters <- unlist(x = strsplit(x = paste(seurat.object.integrated$sample, seurat.object.integrated$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleGlobalCellClusters, col.name = "EachSampleGlobalCellClusters")

# switch the identity class of all cells to reflect sample-specific identities
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleGlobalCellClusters

####################################
### Write out number and fraction of cells per cluster, per sample
####################################
writeLines("\n*** Write out number and fraction of cells per cluster, per sample  ***\n")

OutfileNumbCellsPerClusterPerSample<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerClusterPerSample.tsv", sep="")
OutfileFracCellsPerClusterPerSample<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "FracCellsPerClusterPerSample.tsv", sep="")

seurat.object.integrated_EachSampleGlobalCellClusters.df <- as.data.frame(table(seurat.object.integrated$EachSampleGlobalCellClusters))
rownames(seurat.object.integrated_EachSampleGlobalCellClusters.df) <- seurat.object.integrated_EachSampleGlobalCellClusters.df[,1]

freq_mat <- sapply(rownames(InputsTable), function(dataset) {
  freq_clust <- sapply(sort(unique(seurat.object.integrated$seurat_clusters)), function(cluster) {
    dataset_cluster <- paste(dataset, "_c", cluster, sep = "", collapse = "")
    if ((dataset_cluster %in% row.names(seurat.object.integrated_EachSampleGlobalCellClusters.df)) == TRUE)  {
      freq <- seurat.object.integrated_EachSampleGlobalCellClusters.df[dataset_cluster,"Freq"]
    }else{
      freq <- 0
    }
    freq ## returns freq for the second sapply loop
  })
  freq_clust ## returns freq_clust for the first sapply loop
})

Headers<-paste("cluster", paste(rownames(InputsTable) , sep="", collapse = "\t"), sep = "\t", collapse = "\t")

write.table(Headers, file = OutfileNumbCellsPerClusterPerSample, row.names = F, col.names = F, sep="\t", quote = F)
write.table(freq_mat, file = OutfileNumbCellsPerClusterPerSample, row.names = T, col.names = F, sep="\t", quote = F, append = T)

write.table(Headers, file = OutfileFracCellsPerClusterPerSample, row.names = F, col.names = F, sep="\t", quote = F)
write.table(round((freq_mat / colSums(freq_mat)), 4) , file = OutfileFracCellsPerClusterPerSample, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Get average gene expression for each sample based on global clusters
####################################
writeLines("\n*** Get average gene expression for each sample based on global clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSample_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSample_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSample_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Finding differentially expressed genes: using global cell clusers, compares each cell cluster from one sample vs. the same cluster from other samples
####################################

if (2 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes: using global cell clusers, compares each cell cluster from one sample vs. the same cluster from other samples ***\n")
  
  print(paste("NumberOfClusters=", NumberOfClusters, sep = "", collapse = ""))
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  
  OutfileDiffGeneExpression<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "MarkersPerSampleEquivalentClusters.tsv", sep="")
  HeadersOrder <- paste("cluster1", "cluster2", "gene", "p_val","p_val_adj","avg_logFC", sep = "\t")
  write.table(HeadersOrder,file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="", quote = F)

  for (cluster in unique(seurat.object.integrated$seurat_clusters)) {
    for (dataset1 in rownames(InputsTable)) {
      for (dataset2 in rownames(InputsTable)) {
        Cluster1 <- paste(dataset1, cluster, sep = "_c")
        Cluster2 <- paste(dataset2, cluster, sep = "_c")
        if (Cluster1 == Cluster2) {
        ### Skip
        }else if ((sum(seurat.object.integrated$EachSampleGlobalCellClusters == Cluster1) >= 3) & ((sum(seurat.object.integrated$EachSampleGlobalCellClusters == Cluster2) >= 3))) {
          print (paste(Cluster1, " vs. ", Cluster2, sep = "", collapse = ""))
          seurat.object.integrated.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.integrated, only.pos = F, ident.1 = Cluster1, ident.2 = Cluster2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
          seurat.object.integrated.each_equivalent_cluster.markers$cluster1 <- Cluster1
          seurat.object.integrated.each_equivalent_cluster.markers$cluster2 <- Cluster2
          seurat.object.integrated.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.integrated.each_equivalent_cluster.markers)
          write.table(seurat.object.integrated.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="\t", quote = F, append = T)
        }else{
          print(paste("Skip cluster ", Cluster1, " vs.", Cluster2, "because there were not >= 3 cells in at least one of them", sep = "", collapse = ""))
        }
      }
    }
  }
}

####################################
### Colour dimension reduction plots for each sample based on global clusters and selected genes
####################################

writeLines("\n*** Colour dimension reduction plots for each sample based on global clusters and selected genes ***\n")

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_sample") == T) {
    rm(seurat.object.each_sample)
  }
  seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_sample)

  for (dim_red_method in names(DimensionReductionMethods)) {
    
    ####################################
    ### Colour dimension reduction plots by global cell clusters
    ####################################
    writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by global cell clusters ***\n", sep = "", collapse = ""))

    StopWatchStart$DimRedPlotsByDatasetIntegratedCellClusters$dataset$dim_red_method  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_sample, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
    
    StopWatchEnd$DimRedPlotsByDatasetIntegratedCellClusters$dataset$dim_red_method  <- Sys.time()

    ####################################
    ### Colour dimension reduction plots by selected genes (option -g)
    ####################################
    
    ### To program layout() for more than 3 genes in multiple rows
    if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
      print("No selected genes for dimension reduction plots")
    }else{
      writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by selected genes (option -g) ***\n", sep = "", collapse = ""))
      
      StopWatchStart$DimRedPlotsByDatasetColuredByGenes$dim_red_method  <- Sys.time()
      
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
      
      ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
      ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
      seurat.object.each_sample.sa <- CreateSeuratObject(seurat.object.each_sample@assays$RNA@data)
      seurat.object.each_sample.sa@reductions$umap <- seurat.object.each_sample@reductions$umap
      seurat.object.each_sample.sa@reductions$tsne <- seurat.object.each_sample@reductions$tsne
      
      pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
      print(FeaturePlot(object = seurat.object.each_sample.sa, ncol = nColFeaturePlot, features = ListOfGenesForDimRedPlots, cols = c("lightgrey", "blue"),
                        reduction = dim_red_method, order = T, slot = "data", pt.size = 1.5, min.cutoff = "q0.1", max.cutoff = "q90"))
      dev.off()
      
      StopWatchEnd$DimRedPlotsByDatasetColuredByGenes$dim_red_method  <- Sys.time()
  
      rm(seurat.object.each_sample.sa)
      
    }
  }
}

####################################
### Colour dimension reduction plots for each sample by -infile_colour_dim_red_plots using integrated data
####################################

if (regexpr("^NA$", InfileColourDimRedPlots, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for dimension reduction plots")
}else{
  writeLines("\n*** Colour dimension reduction plots for each sample by -infile_colour_dim_red_plots using integrated data ***\n")
    
  seurat.object.meta.data<-seurat.object.integrated@meta.data
  ExtraCellProperties <- data.frame(read.table(InfileColourDimRedPlots, header = T, row.names = 1))

  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_dim_red_plots and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = ExtraCellProperties)

  NumberOfDatasets <- 0
  for (dataset in rownames(InputsTable)) {
    NumberOfDatasets <- NumberOfDatasets + 1
    print(NumberOfDatasets)
    
    if (exists(x = "seurat.object.each_sample") == T) {
      rm(seurat.object.each_sample)
    }
    seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
    print(seurat.object.each_sample)
    
    for (dim_red_method in names(DimensionReductionMethods)) {
      
      StopWatchStart$DimRedPlotsByMetadata$dataset$dim_red_method  <- Sys.time()
  
      # Generating outfile
      # Note DimPlot() takes the entire current device (pdf)
      # even if using layout(matrix(...)) or  par(mfrow=())
      # Thus each property plot is written to a separate page of a single *pdf outfile
      for (property in colnames(ExtraCellProperties)) {

        if (sum(is.na(seurat.object.each_sample[[property]]==T)) == ncol(seurat.object.each_sample)) {
          
          print(paste("No extra barcode-attributes '", property, "' found for sample ", dataset))
          
        }else{
          
          IntegratedDimRedPlotPdf <- paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBy_", property, ".pdf", sep="")
          pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(ExtraCellProperties[,property])))))
          
          if ( (sum(ExtraCellProperties[,property] %in% 0:1 == T)) == (nrow(ExtraCellProperties)) ) { ## is binary
            CellsToHighlight <- rownames(ExtraCellProperties)[ExtraCellProperties[,property]==1]
            plots <- DimPlot(seurat.object.integrated, group.by = property, combine = F, reduction = dim_red_method, label = F, cells.highlight = CellsToHighlight)
            plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
            print(CombinePlots(plots))
            
          }else{
            plots <- DimPlot(seurat.object.integrated, group.by = property, combine = F, reduction = dim_red_method, label = T)
            plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top", legend.text.align = 0) + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
            print(CombinePlots(plots))
            
          }
          
          dev.off()
        }
      }
      
      StopWatchEnd$DimRedPlotsByMetadata$dataset$dim_red_method  <- Sys.time()
      
    }
  }
}

####################################
### FOR EACH SAMPLE uses integrated counts and RE-CLUSTER cells
### Then get DGE and colour dimension reduction plots 
####################################
writeLines("\n*** FOR EACH SAMPLE uses integrated counts and RE-CLUSTER cells. Then get DGE and colour dimension reduction plots ***\n")

# switch the identity class of all cells to reflect sample
Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample

# Headers for sample-specific cell clusters
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = ""), sep="\t")
OutfileEachSampleClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "CellClusters.tsv", sep="")
write.table(Headers,file = OutfileEachSampleClusters, row.names = F, col.names = F, sep="\t", quote = F)
#
Headers<-paste("Sample", "Number_of_clusters" ,sep="\t", collapse = "")
OutfileEachSampleNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)
#
Headers<-paste("Cluster", "Number_of_cells" ,sep="\t", collapse = "")
OutfileEachSampleNumbCellsPerCluster<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "NumbCellsPerCluster.tsv", sep="")
write.table(Headers,file = OutfileEachSampleNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F)

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_sample") == T) {
    rm(seurat.object.each_sample)
  }
  seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_sample)

  StopWatchStart$ClusterEachSampleCellsFromInteg$dataset  <- Sys.time()
  
  ####################################
  ### Re-cluster each sample cells using integrated data
  ####################################
  
  options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
  seurat.object.each_sample <- FindNeighbors(object = seurat.object.each_sample, dims = PcaDimsUse)
  seurat.object.each_sample <- FindClusters(object = seurat.object.each_sample, resolution = Resolution)
  ClustersThisSample <- unlist(x = strsplit(x = paste(dataset, seurat.object.each_sample@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))

  StopWatchEnd$ClusterEachSampleCellsFromInteg$dataset  <- Sys.time()
  
  ####################################
  ### Write out each sample cell clusters
  ####################################
  writeLines("\n*** Write out each sample cell clusters ***\n")
  
  StopWatchStart$WriteClustersEachSampleCellClusterTables$dataset <- Sys.time()
  
  ### Cell cluster identities
  write.table(paste(colnames(seurat.object.each_sample),ClustersThisSample, sep = "\t", collapse = "\n"),
              file = OutfileEachSampleClusters, row.names = F, col.names = F, quote = F, append = T)

  ### Number of clusters per sample
  CellNames<-rownames(seurat.object.each_sample@meta.data)
  ClusterIdent <-seurat.object.each_sample@meta.data$seurat_clusters
  NumberOfClusters<-length(unique(ClusterIdent))
  write(x=paste(dataset, NumberOfClusters, sep = "\t", collapse = ""), file = OutfileEachSampleNumbClusters, append = T)
  
  ### Number of cells per cluster, per sample
  NumberOfCellPerClusterPerSample<-table(ClusterIdent)
  rownames(NumberOfCellPerClusterPerSample) <- paste(dataset, rownames(NumberOfCellPerClusterPerSample), sep = "_c")
  write.table(NumberOfCellPerClusterPerSample, file = OutfileEachSampleNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F, append = T)

  StopWatchEnd$WriteClustersEachSampleCellClusterTables$dataset <- Sys.time()
  
  if (5 %in% RequestedDiffGeneExprComparisons == T) {
    
    ####################################
    ### Finding differentially expressed genes: using each sample re-clustered, compares each cell cluster vs. the rest of cells
    ####################################
    writeLines("\n*** Finding differentially expressed genes: using each sample re-clustered, compares each cell cluster vs. the rest of cells ***\n")
    
    print(paste("NumberOfClusters=", NumberOfClusters, sep = "", collapse = ""))
    
    StopWatchStart$FindDiffMarkersReclusteredVsRestOfCells$dataset  <- Sys.time()
    
    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_sample@meta.data))
    seurat.object.each_sample.markers <- FindAllMarkers(object = seurat.object.each_sample, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_sample.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC")]
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", dataset, "_MarkersPerCluster.tsv", sep=""), row.names = F, sep="\t", quote = F)
  
    StopWatchEnd$FindDiffMarkersReclusteredVsRestOfCells$dataset  <- Sys.time()
  
  }
  
  ####################################
  ### Colour dimension reduction plots for each sample re-clustered cells
  ####################################
  writeLines("\n*** Colour dimension reduction plots for each sample re-clustered cells ***\n")
  
  for (dim_red_method in names(DimensionReductionMethods)) {
    
    StopWatchStart$DimRedPlotsByEachSampleReclusteredCellCluster$dataset$dim_red_method  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_sample, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
    
    StopWatchEnd$DimRedPlotsByEachSampleReclusteredCellCluster$dataset$dim_red_method <- Sys.time()
    
  }
}

####################################
### Load each sample re-cluster assignments
####################################
writeLines("\n*** Load each sample re-cluster assignments ***\n")

EachSampleCellReClusters <- data.frame(read.table(OutfileEachSampleClusters, header = T, row.names = 1, sep = "\t"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleCellReClusters, col.name = "EachSampleCellReClusters")

# switch the identity class of all cells to reflect each sample clusters
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleCellReClusters

####################################
### Get average gene expression for each sample re-clustered clusters
####################################
writeLines("\n*** Get average gene expression for each sample re-clustered clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "PerSample_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "PerSample_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir, "/", PrefixOutfiles, ".", ProgramOutdir, "_EachSampleReclustered_", "PerSample_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### FOR EACH SAMPLE TYPE uses integrated/normalized counts and generates dimension reduction plots for each sample type using:
### a) sample type
### b) global clusters
### c) re-clustered cells of each sample - also gets DGE
####################################
writeLines("\n*** FOR EACH SAMPLE TYPE uses integrated/normalized counts and generates dimension reduction plots for each sample type using:
a) sample type
b) global clusters
c) re-clustered cells of each sample - also gets DGE  ***\n")

# Get sample_type assignments replacing sample ids using list_DatasetToType()
seurat.object.integrated$sample_type<-seurat.object.integrated$sample
for (dataset in rownames(InputsTable)) {
seurat.object.integrated$sample_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$sample_type)
}
SampleTypes <- unique(seurat.object.integrated$sample_type)

####################################
### Colour dimension reduction plots by sample type
####################################
writeLines("\n*** Colour dimension reduction plots by sample type ***\n")

# switch identity to sample type
Idents(object = seurat.object.integrated) <- seurat.object.integrated$sample_type

for (dim_red_method in names(DimensionReductionMethods)) {
  
  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " coloured by sample_type ***\n", sep = "", collapse = ""))
  
  StopWatchStart$DimRedOPlotColourBySampleType$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("sample_type"), combine = F, reduction = dim_red_method, label = F, label.size = 10)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = NumbRowsForLegendPerSampleType, byrow = T, override.aes = list(size = 2))))
  IntegratedDimRedPlotPdf<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySampleType.pdf", sep="")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  print(IntegratedDimRedPlotPdf)
  
  StopWatchEnd$DimRedOPlotColourBySampleType$dim_red_method  <- Sys.time()
  
}

####################################
### Write out each sample type global cluster assignments
####################################
writeLines("\n*** Write out each sample type global cluster assignments ***\n")

# Headers for sample-type-specific cell clusters
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = ""), sep="\t")
OutfileEachSampleTypeGlobalClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "EachSampleType_CellClusters.tsv", sep="")
GlobalClustersBySampleType<- unlist(x = strsplit(x = paste(seurat.object.integrated@meta.data$sample_type, seurat.object.integrated@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
write.table(Headers,file = OutfileEachSampleTypeGlobalClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(paste(colnames(seurat.object.integrated),GlobalClustersBySampleType, sep = "\t", collapse = "\n"),
            file = OutfileEachSampleTypeGlobalClusters, row.names = F, col.names = F, quote = F, append = T)

OutfileNumbCellsPerClusterPerSampleType<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "OutfileNumbCellsPerClusterPerSampleType.tsv", sep="")
Headers<-paste("Cluster", "Number_of_cells" ,sep="\t", collapse = "")
write.table(Headers, file = OutfileNumbCellsPerClusterPerSampleType, row.names = F, col.names = F, sep="\t", quote = F)
write.table(table(GlobalClustersBySampleType),file = OutfileNumbCellsPerClusterPerSampleType, row.names = F, col.names = F, sep="\t", quote = F, append = T)

####################################
### Load each sample type global cluster assignments
####################################
writeLines("\n*** Load each sample type global cluster assignments ***\n")

EachSampleTypeGlobalCellClusters <- data.frame(read.table(OutfileEachSampleTypeGlobalClusters, header = T, row.names = 1, sep = "\t"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleTypeGlobalCellClusters, col.name = "EachSampleTypeGlobalCellClusters")

# switch the identity class of all cells to reflect each sample clusters
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleTypeGlobalCellClusters

####################################
### Get average gene expression for each sample type using global clusters
####################################
writeLines("\n*** Get average gene expression for each sample type using global clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSampleType_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSampleType_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerSampleType_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Finding differentially expressed genes: using global cell clusers, for each sample type, compares each cell cluster vs. the same cluster from other sample types
####################################

if (4 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes: using global cell clusers, for each sample type, compares each cell cluster vs. the same cluster from other sample types ***\n")
  
  NumberOfClusters <- length(unique(seurat.object.integrated$seurat_clusters))
  
  print(paste("NumberOfClusters=", NumberOfClusters, sep = "", collapse = ""))
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  
  OutfileDiffGeneExpression<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "MarkersPerSampleTypeEquivalentClusters.tsv", sep="")
  HeadersOrder <- paste("cluster1", "cluster2", "gene", "p_val","p_val_adj","avg_logFC", sep = "\t")
  write.table(HeadersOrder,file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="", quote = F)
  
  for (cluster in unique(seurat.object.integrated$seurat_clusters)) {
    for (dataset_type1 in unique(InputsTable[,"DatasetType"])) {
      for (dataset_type2 in unique(InputsTable[,"DatasetType"])) {
        Cluster1 <- paste(dataset_type1, cluster, sep = "_c")
        Cluster2 <- paste(dataset_type2, cluster, sep = "_c")
        if (Cluster1 == Cluster2) {
          ### Skip
        }else if ((sum(seurat.object.integrated$EachSampleTypeGlobalCellClusters == Cluster1) >= 3) & ((sum(seurat.object.integrated$EachSampleTypeGlobalCellClusters == Cluster2) >= 3))) {
          print (paste(Cluster1, " vs. ", Cluster2, sep = "", collapse = ""))
          seurat.object.integrated.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.integrated, only.pos = F, ident.1 = Cluster1, ident.2 = Cluster2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
          seurat.object.integrated.each_equivalent_cluster.markers$cluster1 <- Cluster1
          seurat.object.integrated.each_equivalent_cluster.markers$cluster2 <- Cluster2
          seurat.object.integrated.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.integrated.each_equivalent_cluster.markers)
          write.table(seurat.object.integrated.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="\t", quote = F, append = T)
          rm(seurat.object.integrated.each_equivalent_cluster.markers)
        }else{
          print(paste("Skip cluster ", Cluster1, " vs.", Cluster2, "because there were not >= 3 cells in at least one of them", sep = "", collapse = ""))
        }
      }
    }
  }
}

####################################
### Colour dimension reduction plots for each sample type by glabal clusters and gets DGE
####################################
writeLines("\n*** Colour dimension reduction plots for each sample type by glabal clusters and gets DGE ***\n")

# switch identity to sample type global cluster
Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample_type

NumberOfSampleTypes <- 0
for (sample_type in SampleTypes) {
  NumberOfSampleTypes <- NumberOfSampleTypes + 1
  print(NumberOfSampleTypes)
  
  if (exists(x = "seurat.object.each_sample_type") == T) {
    rm(seurat.object.each_sample_type)
  }
  seurat.object.each_sample_type <- subset(x = seurat.object.integrated, idents = sample_type)

  ####################################
  ### Colour dimension reduction plots by each sample type using global cell clusters
  ####################################

  writeLines(paste("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for ", sample_type, " using global cell clusters ***\n", sep = "", collapse = ""))

  for (dim_red_method in names(DimensionReductionMethods)) {

    StopWatchStart$DimRedOPlotSampleTypeColourByCellCluster$sample_type$dim_red_method  <- Sys.time()

    plots <- DimPlot(seurat.object.each_sample_type, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 10)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste(Tempdir,"/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", sample_type,"_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()

    StopWatchEnd$DimRedOPlotSampleTypeColourByCellCluster$sample_type$dim_red_method  <- Sys.time()

  }

  if (3 %in% RequestedDiffGeneExprComparisons == T) {
    
    ####################################
    ### Finding differentially expressed genes: using global cell clusers, for each sample type, compares each cell cluster vs. the rest of cells
    ####################################
    writeLines("\n*** Finding differentially expressed genes: using global cell clusers, for each sample type, compares each cell cluster vs. the rest of cells ***\n")
    
    NumberOfClusters <- length(unique(seurat.object.each_sample_type$seurat_clusters))
    
    print(paste("NumberOfClusters=", NumberOfClusters, sep = "", collapse = ""))
  
    Idents(object = seurat.object.each_sample_type) <- seurat.object.each_sample_type@meta.data$EachSampleTypeGlobalCellClusters
  
    StopWatchStart$FindDiffMarkersGlobalClustersVsRestOfCells$sample_type  <- Sys.time()
    
    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_sample_type@meta.data))
    seurat.object.each_sample_type.markers <- FindAllMarkers(object = seurat.object.each_sample_type, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_sample_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC")]
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_" , sample_type, "_MarkersPerCluster.tsv", sep=""), row.names = F, sep="\t", quote = F)
  
    StopWatchEnd$FindDiffMarkersGlobalClustersVsRestOfCells$sample_type  <- Sys.time()
  }
}

####################################
### Reclusters each sample type cells and gets DGE
####################################

### Write out headers for sample type cell re-cluster outfile
OutfileEachSampleTypeReClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "CellClusters.tsv", sep="")
Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = ""), sep="\t")
write.table(Headers,file = OutfileEachSampleTypeReClusters, row.names = F, col.names = F, sep="\t", quote = F)

NumberOfSampleTypes <- 0
for (sample_type in SampleTypes) {
  NumberOfSampleTypes <- NumberOfSampleTypes + 1
  print(NumberOfSampleTypes)
  
  if (exists(x = "seurat.object.each_sample_type") == T) {
    rm(seurat.object.each_sample_type)
  }
  seurat.object.each_sample_type <- subset(x = seurat.object.integrated, idents = sample_type)
  print(seurat.object.each_sample_type)
  
  ####################################
  ### Re-cluster each sample type cells using integrated data
  ####################################
  writeLines("\n*** Re-cluster each sample type cells using integrated data ***\n")
  
  StopWatchStart$ClusterEachSampleTypeCellsFromInteg$sample_type  <- Sys.time()

  options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
  seurat.object.each_sample_type <- FindNeighbors(object = seurat.object.each_sample_type, dims = PcaDimsUse)
  seurat.object.each_sample_type <- FindClusters(object = seurat.object.each_sample_type, resolution = Resolution)
  ClustersThisSampleType <- unlist(x = strsplit(x = paste(sample_type, seurat.object.each_sample_type@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
  
  StopWatchEnd$ClusterEachSampleTypeCellsFromInteg$sample_type  <- Sys.time()
  
  ####################################
  ### Write out each sample type re-clustered cell clusters
  ####################################
  writeLines("\n*** Write out each sample type re-clustered cell clusters ***\n")

  StopWatchStart$WriteClustersEachSampleTypeCellClusterTables$sample_type <- Sys.time()
  
  write.table(paste(colnames(seurat.object.each_sample_type), ClustersThisSampleType, sep = "\t", collapse = "\n"),
              file = OutfileEachSampleTypeReClusters, row.names = F, col.names = F, quote = F, append = T)
  
  CellNames<-rownames(seurat.object.each_sample_type@meta.data)
  ClusterIdent <-seurat.object.each_sample_type@meta.data$seurat_clusters
  NumberOfClusters<-length(unique(ClusterIdent))
  OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", sample_type, "_NumbCellClusters", ".tsv", sep="")
  write(x=NumberOfClusters,file = OutfileNumbClusters)
  
  StopWatchEnd$WriteClustersEachSampleTypeCellClusterTables$sample_type <- Sys.time()
  
  if (6 %in% RequestedDiffGeneExprComparisons == T) {
    
    ####################################
    ### Finding differentially expressed genes: using each sample type re-clustered, compares each cell cluster vs. the rest of cells
    ####################################
    writeLines("\n*** Finding differentially expressed genes: using each sample type re-clustered, compares each cell cluster vs. the rest of cells ***\n")
    
    print(paste("NumberOfClusters=", NumberOfClusters, sep = "", collapse = ""))
  
    StopWatchStart$FindDiffMarkersReclusteredVsRestOfCells$sample_type  <- Sys.time()
    
    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_sample_type@meta.data))
    seurat.object.each_sample_type.markers <- FindAllMarkers(object = seurat.object.each_sample_type, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_sample_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC")]
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", sample_type, "_MarkersPerCluster.tsv", sep=""), row.names = F, sep="\t", quote = F)
  
    StopWatchEnd$FindDiffMarkersReclusteredVsRestOfCells$sample_type  <- Sys.time()
  
  }
  
  ####################################
  ### Colour dimension reduction plots by each sample type re-clustered cells
  ####################################

  writeLines("\n*** Colour dimension reduction plots by each sample type re-clustered cells ***\n")
  
  for (dim_red_method in names(DimensionReductionMethods)) {

    StopWatchStart$UmapAndTsneColourByEachSampleTypeCellCluster$sample_type  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_sample_type, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
    IntegratedUMAPPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", sample_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf", sep="")
    pdf(file=IntegratedUMAPPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
  
    StopWatchEnd$UmapAndTsneColourByEachSampleTypeCellCluster$sample_type <- Sys.time()
    
    ### To program layout() for more than 3 genes in multiple rows
    if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
      print("No selected genes for dimension reduction plots")
    }else{
      
      StopWatchEnd$DimRedPlotsBySampleTypeColuredByGenes$dim_red_method  <- Sys.time()
  
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
      
      ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
      ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
      seurat.object.each_sample_type.sa <- CreateSeuratObject(seurat.object.each_sample_type@assays$RNA@data)
      seurat.object.each_sample_type.sa@reductions$umap <- seurat.object.each_sample_type@reductions$umap
      seurat.object.each_sample_type.sa@reductions$tsne <- seurat.object.each_sample_type@reductions$tsne
      
      pdf(file=paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_", sample_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBySelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
      print(FeaturePlot(object = seurat.object.each_sample_type.sa, ncol = nColFeaturePlot, features = ListOfGenesForDimRedPlots, cols = c("lightgrey", "blue"),
                        reduction = dim_red_method, order = T, slot = "data", pt.size = 1.5, min.cutoff = "q0.1", max.cutoff = "q90"))
      dev.off()
      
      rm(seurat.object.each_sample_type.sa)
      
      StopWatchEnd$DimRedPlotsBySampleTypeColuredByGenes$dim_red_method  <- Sys.time()
    
    }
  }
}

####################################
### Load each sample type re-cluster assignments
####################################
writeLines("\n*** Load each sample type re-cluster assignments ***\n")

EachSampleTypeReclusteredCellClusters <- data.frame(read.table(OutfileEachSampleTypeReClusters, header = T, row.names = 1, sep = "\t"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleTypeReclusteredCellClusters, col.name = "EachSampleTypeReclusteredCellClusters")

# switch the identity class of all cells to reflect each sample_type clusters
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleTypeReclusteredCellClusters

####################################
### Get average expression for each sample type re-clustered clusters
####################################
writeLines("\n*** Get average expression for each sample type re-clustered clusters ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "PerSampleType_AverageGeneExpression_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "PerSampleType_AverageGeneExpression_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".", ProgramOutdir, "_EachSampleTypeReclustered_", "PerSampleType_AverageGeneExpression_integrated.tsv", sep="")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Saving the full R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {

  writeLines("\n*** Saving the full R object ***\n")
  
  StopWatchStart$SaveRDS  <- Sys.time()
  
  OutfileRDS<-paste(Tempdir,"/",PrefixOutfiles,".", ProgramOutdir, "_integrated_object_full.rds", sep="")
  saveRDS(seurat.object.integrated, file = OutfileRDS)
  
  StopWatchEnd$SaveRDS  <- Sys.time()

}else{
  
  writeLines("\n*** Not saving the full R object ***\n")
  
}

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo.txt", sep="")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste(Tempdir, "/" , PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt", sep="")
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

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(paste(Tempdir, sep="", collapse = ""))
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  writeLines(paste(Outdir,"/", ProgramOutdir, "/", sep="", collapse = ""))
  outfiles_to_move <- list.files(Tempdir, pattern = PrefixOutfiles, full.names = F)
  sapply(outfiles_to_move,FUN=function(eachFile) {
    ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
    file.copy(from=paste(Tempdir,"/",eachFile,sep=""), to=paste(Outdir, "/", ProgramOutdir, "/", eachFile,sep=""), overwrite=T)
    file.remove(paste(Tempdir,"/",eachFile,sep=""))
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
