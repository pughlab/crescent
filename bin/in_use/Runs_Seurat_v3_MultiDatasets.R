####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/v3.1/integration.html
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads scRNA-seq data and generate QC plots for each dataset
### 2) Normalizes and scales each dataset
### 3) Merges and integrates datasets correcting batch effects
### 4) Process integrated datasets as a whole, including:
###    - dimension reduction
###    - 'global' cell clustering
###    - UMAP/tSNE plots by global cell clusters, sample, sample type, requested genes and metadata
###    - differential gene expression (DGE)
###    - average gene expression
### 5) For each dataset on its own uses global clusters and repeats analysis from step 3
### 6) For each dataset on its own re-clusters cells and repeats analysis from step 3
### 7) For each dataset type on its own uses global clusters and repeats analysis from step 3
### 8) For each dataset type on its own re-clusters cells and repeats analysis from step 3
### 9) Obtains DGE for metadata classes
### 10) Saves log files
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/Runs_Seurat_v3_MultiDatasets.R -h'
### for help
####################################

####################################
### THINGS TO DO:
### 1) See Run_Seurat_v3_SingleDataset.R for a list of functions to be implemented
### 2) In:
###    `for (dataset in rownames(InputsTable)) {
###       seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
###     }`
###   Need to avoid that the gsub replaces partial strings, for example 'SMTR03' will be replaced by list_DatasetToType[[SMTR03]] in both
###   'SMTR03' iself and in 'SMTR03t1_NonRad'. The second case is undesired.
###   Generating files like:
###   SMTR_res1.SEURAT_GlobalClustering_Rad_unknownt1_NonRad_TSNEPlot_ColourByCellClusters.pdf
###   When it should be:
###   SMTR_res1.SEURAT_GlobalClustering_Rad_unknown_TSNEPlot_ColourByCellClusters.pdf
### 3) To program a loop to avoid an error when integrating very heretogeneous datasets--where few anchors are expected:
###    Error in nn2(data = c(-0.0268782296438318, -0.0333332648557412, -0.032492403275445,  : 
###       Cannot find more nearest neighbours than there are points
###    To avoid this error one can modify `k.filter` from 200 (default) to lower values, e.g. `k.filter = 160` worked out for
###    Pugh/dePerrot labe SMARTER datasets, where sample SMTR05t1_NonRad integrated vs. SMTR0[2-4] produced the above error
###    Documented here: https://github.com/satijalab/seurat/issues/997
####################################

####################################
### COMMENTS ON FINDING DIFFERENTIALLY EXPRESSED GENES
### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object.each_dataset, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object.each_dataset, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
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
### Computing time was reduced in step IntegrateData() when using `-k 1,2`. Whereas FindIntegrationAnchors() was similar:
### Step/time(minutes)      Using_-k_NA   Using_-k_1,2
### FindIntegrationAnchors	3.336035	    3.118804
### IntegrateData           2.968116	    1.666707
### In both cases used `-u MAX -b 10000` in a 3.1-GHz Intel Core i5 CPU with 2 cores and 16 GB RAM
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(DropletUtils)) # (Bioconductor) to handle MTX/H5 format files. Note it has about the same speed than library(earlycross) which can't handle H5
suppressPackageStartupMessages(library(Seurat))       # (CRAN) to run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # (CRAN) needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # (CRAN) to calculate the percentages of extra properties to be plotted
suppressPackageStartupMessages(library(data.table))   # (CRAN) to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # (CRAN) to run parallel processes
suppressPackageStartupMessages(library(gtools))       # (CRAN) to do alphanumeric sorting. Only needed if using `-w Y`.
suppressPackageStartupMessages(library(loomR))        # (GitHub mojaveazure/loomR) needed for fron-end display of data. Only needed if using `-w Y`.     
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

ThisScriptName <- "Runs_Seurat_v3_MultiDatasets.R"
ProgramOutdir  <- "SEURAT"

####################################
### Get inputs from command line argumets
####################################
#
option_list <- list(
  make_option(c("-i", "--inputs_list"), default="NA",
              help="Path/name to a <tab> delimited file with one dataset per row and
                the following 12 columns specifying filters/details for each dataset:
                1) unique dataset ID (e.g. 'd1')
                2) /path_to/dataset
                3) dataset type (e.g. 'control' or 'treatment') or use 'type' to skip using dataset types
                4) dataset format (either 'MTX', 'TSV' or 'HDF5')
                5) dataset minimum mitochondrial fraction (e.g. '0')
                6) dataset maximum mitochondrial fraction (e.g. '0.05' for NucSeq or '0.2' for whole cell scRNA-seq)
                7) dataset minimum ribosomal fraction (e.g. '0')
                8) dataset maximum ribosomal fraction (e.g. '0.75')
                9) dataset minimum number of genes (e.g. '50')
                10) dataset maximum number of genes (e.g. '8000')
                11) dataset minimum number of reads (e.g. '1')
                12) dataset maximum number of reads (e.g. '80000')
                Default = 'No default. It's mandatory to specify this parameter'

                Notes:
                (a) The order of the list of datasets in --inputs_list may influence the results, including number of clusters,
                t-SNE/UMAP and differentially expressed genes. List datasets better measured first.
                
                (b) Column 4 indicates the dataset format. It must be in either:
                'MTX'  and column 2 must be the path/name to an MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files.
                'TSV'  and column 2 must be the path/name to a <tab> delimited *file* with genes in rows vs. barcodes in columns.
                'HDF5' and column 2 must be the path/name of a *file* in hdf5 format (e.g. from Cell Ranger)
                Note 'MTX' files can be the outputs from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz"),
  #
  make_option(c("-j", "--inputs_remove_barcodes"), default="NA",
              help="Path/name to a <tab> delimited list of barcodes to be removed from analysis, like:
                d1  AAACCTGAGCTCCCAG
                d2  AAACCTGTCACCATAG
                d3  AAACCTGTCAGCTTAG
                Note: the barcode id must include the dataset ID
                Or type 'NA' to include all barcodes
                Default = 'NA'"),
  #
  make_option(c("-k", "--save_filtered_data"), default="N",
              help="Indicates if filtered raw and normalized data should be saved as MTX files. Type 'y/Y' or 'n/N'
                Default = 'N'"),
  #
  make_option(c("-l", "--save_unfiltered_data"), default="N",
              help="Indicates if unfiltered raw data (exactly as inputted by --inputs_list) should be saved as MTX files. Type 'y/Y' or 'n/N'
                This may be useful to share with collaborators along results in a single zipped file
                Default = 'N'"),
  #
  make_option(c("-z", "--reference_datasets"), default="NA",
              help="A list of <comma> delimited number of row(s) from --inputs_list of datasets to be used as reference(s) for integration
                Or type 'NA' to run all-vs-all dataset pairwise comparisons (this is more time and memory consuming than using references).
                If references are used, then references will be integrated and anchors identified between non-reference and reference datasets,
                but not anchors will be identified between non-reference datasets (this saves time and memory)
                Default = 'NA'"),
  #
  make_option(c("-r", "--resolution"), default="1",
              help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of cell clusters
                Default = '1'"),
  #
  make_option(c("-v", "--clustering_inputs"), default="1",
              help="One or more of the following options. If more than one, pass the choice <comma> delimited, e.g. '1,2,3'
                '1' = globaly clusters the integrated datasets (mandatory in this script)
                '2' = re-clusters each datasets split from the integrated datasets. Needed for `-f 6`.
                '3' = re-clusters each dataset type (column 3 from -inputs_list) from the integrated datasets. Needed for `-f 7`.
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
                Default = 'NA' (i.e. no metadata to be used for plots or DGE detection)"),
  #
  make_option(c("-g", "--infile_selected_genes"), default="NA",
              help="A path/name to a file with a list of gene identifiers (one-per-row) whose expression will be mapped into the dimension reduction plots
                Default = 'NA' (i.e. no genes to be plotted)"),
  #
  make_option(c("-m", "--apply_list_genes"), default="0",
              help="One or more of the following options. If more than one, pass the choice <comma> delimited, e.g. '1,2,3'
                '0' = no genes to be mapped into dimension reduction plots
                '1' = using cells from all integrated datasets to map genes requested by `-g`
                '2' = using cells from each dataset to map genes requested by `-g`
                '3' = using cells from each dataset type to map genes requested by `-g`
                Default = '0'"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check file:
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                Default = '10'"),
  #
  make_option(c("-n", "--k_filter"), default="200",
              help="Integrating highly heterogenous datasets can lead to a too small number of anchors,
                resulting in an error 'Cannot find more nearest neighbours than there are points'
                This can be avoided by using a -k_filter value smaller than the default (e.g. '150')
                Default = '200'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each differentially expressed genes test returns only hits that have an UNcorrected p-value < return_thresh
                Default = '0.01'"),
  #
  make_option(c("-f", "--diff_gene_expr_comparisons"), default="1",
              help="One or more of the following options. If more than one, pass the choice <comma> delimited, e.g. '1,2,3'
                '0' = no differentially expressed genes are computed
                '1' = using global cell clusers, compares each cell cluster vs. the rest of cells. Needs `-v 1`.
                '2' = using global cell clusers, for each dataset, compares each cell cluster vs. the rest of cells. Needs `-v 1`.
                '3' = using global cell clusers, for each dataset, compares each cell cluster vs. the same cluster from other datasets. Needs `-v 1`.
                '4' = using global cell clusers, for each dataset type, compares each cell cluster vs. the rest of cells. Needs `-v 1`.
                '5' = using global cell clusers, for each dataset type, compares each cell cluster vs. the same cluster from other dataset types. Needs `-v 1`.
                '6' = using re-clustered cells, for each dataset, compares each cell cluster vs. the rest of cells. Needs `-v 2`.
                '7' = using re-clustered cells, for each dataset type, compares each cell cluster vs. the rest of cells. Needs `-v 3`.
                '8' = using metadata annotations, compares each cell class specified by `-b` and `-c` vs. the rest of cells
                '9' = using metadata annotations, for each dataset, compares each cell class specified by `-b` and `-c` vs. the rest of cells
                '10' = using metadata annotations, for each dataset, compares each cell class specified by `-b` and `-c` vs. the same class from other datasets
                '11' = using metadata annotations, for each dataset type, compares each cell class specified by `-b` and `-c` vs. the rest of cells
                '12' = using metadata annotations, for each dataset type, compares each cell class specified by `-b` and `-c` vs. the same class from other dataset types
                Default = '1'"),
  #
  make_option(c("-b", "--metadata_column_names_for_dge"), default="NA",
              help="Only needed if using -f 8 to 12. It indicates <comma> delimited column names of --infile_metadata to be used for differential gene expression
                Type 'NA' if not using -f 8 to 12
                Default = 'NA'"),
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
  make_option(c("-x", "--minio_path"), default="NA",
              help="Used with the 'run_cwl' option above for mounting input data files in 'inputs_list' to CWL containers. 
                Default = 'NA'"),
  #
  make_option(c("-a", "--max_global_variables"), default="4000",
              help="Indicates maximum allowed total size (in bytes) of global variables identified.
                Used by library(future) to prevent too large exports
                Default = '4000' for 4000 MiB")
)

opt <- parse_args(OptionParser(option_list=option_list))

InputsList              <- opt$inputs_list
InfileRemoveBarcodes    <- opt$inputs_remove_barcodes
ReferenceDatasets       <- opt$reference_datasets
SaveFilteredData        <- opt$save_filtered_data
SaveUnFilteredData      <- opt$save_unfiltered_data 
Resolution              <- as.numeric(opt$resolution)
ClusteringInputs        <- opt$clustering_inputs
Outdir                  <- opt$outdir
PrefixOutfiles          <- opt$prefix_outfiles
InfileMetadata          <- opt$infile_metadata
InfileSelectedGenes     <- opt$infile_selected_genes
ApplySelectedGenes      <- opt$apply_list_genes
PcaDimsUse              <- c(1:as.numeric(opt$pca_dimensions))
KFilter                 <- as.numeric(opt$k_filter)
ThreshReturn            <- as.numeric(opt$return_threshold)
DiffGeneExprComparisons <- opt$diff_gene_expr_comparisons
MetadataColNamesForDge  <- opt$metadata_column_names_for_dge
NumbCores               <- opt$number_cores
SaveRObject             <- opt$save_r_object
RunsCwl                 <- opt$run_cwl
MaxGlobalVariables      <- as.numeric(opt$max_global_variables)
MinioPath               <- opt$minio_path

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
    "CELL_FRACTIONS", 
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES", 
    "DIMENSION_REDUCTION_COORDINATE_TABLES", 
    "DIMENSION_REDUCTION_PLOTS",
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS",
    "SCTRANSFORM_PSEUDO_BULK",
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS",
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
    "CELL_FRACTIONS", 
    "DIFFERENTIAL_GENE_EXPRESSION_TABLES", 
    "DIMENSION_REDUCTION_COORDINATE_TABLES", 
    "DIMENSION_REDUCTION_PLOTS", 
    "FILTERED_DATA_MATRICES",
    "LOG_FILES",
    "QC_PLOTS", 
    "QC_TABLES", 
    "R_OBJECTS",
    "SCTRANSFORM_PSEUDO_BULK",
    "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS",
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

OutfileOptionsUsed<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_UsedOptions.txt")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_RSessionInfo.txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)

####################################
### Define number of cores and RAM for parallelization
####################################

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (regexpr("^[0-9]+$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste0("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nintegrates_datasets_with_seurat_dev.R -h\n\n"))
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
RequestedClusteringInputs        = unlist(strsplit(ClusteringInputs, ","))
RequestedApplySelectedGenes      = unlist(strsplit(ApplySelectedGenes, ","))

DefaultParameters <- list(
  ### Parameters for QC plots
  CellPropertiesToQC = c("nFeature_RNA", "nCount_RNA", "mito.fraction", "ribo.fraction"),
  
  ### Parameters for Seurat filters
  MinCells = 3,
  
  ### Parameters for Cluster Biomarkers
  FindAllMarkers.MinPct     =  0.25,
  FindAllMarkers.ThreshUse  =  0.25,
  
  ### Parameters for dimmension reduction plots
  MinNumberOfCellsToReducePerplexity = 150,
  ReducedPerplexity = 7,
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
    stop(paste0("Parameter -", param, " can't be 'NA' (default). Use option -h for help."))
  }
}

####################################
### Check that optional parameters are compatible with each other
####################################
writeLines("\n*** Check that optional parameters are compatible with each other ***\n")

#### Differential gene expression options

if ((1 %in% RequestedDiffGeneExprComparisons == T) | (2 %in% RequestedDiffGeneExprComparisons == T) | (3 %in% RequestedDiffGeneExprComparisons == T) | (4 %in% RequestedDiffGeneExprComparisons == T) | (5 %in% RequestedDiffGeneExprComparisons == T))  {
  if (1 %in% RequestedClusteringInputs == T) {
    writeLines("\nOk\n")
  }else{
    stop("ERROR!!! option `-f [1|2|3|4|5]` needs option `-v 1`")
  }
}

if (6 %in% RequestedDiffGeneExprComparisons == T) {
  if (2 %in% RequestedClusteringInputs == T) {
    writeLines("\nOk\n")
  }else{
    stop("ERROR!!! options `-f 6` needs option `-v 2`")
  }
}

if (7 %in% RequestedDiffGeneExprComparisons == T) {
  if (3 %in% RequestedClusteringInputs == T) {
    writeLines("\nOk\n")
  }else{
    stop("ERROR!!! option `-f 7` needs option `-v 3`")
  }
}

#### Metadata options

if (length(grep('^NA$', InfileMetadata, perl = T))) {
  writeLines("\nNo metadata will be used\n")
}else{
  CellPropertiesFromMetadata <- data.frame(read.table(InfileMetadata, header = T, row.names = 1, check.names = F))
}

if ((8 %in% RequestedDiffGeneExprComparisons == T) | 
    (9 %in% RequestedDiffGeneExprComparisons == T) | 
    (10 %in% RequestedDiffGeneExprComparisons == T) |
    (11 %in% RequestedDiffGeneExprComparisons == T) | 
    (12 %in% RequestedDiffGeneExprComparisons == T)
)  {
  if (length(grep('^NA$', InfileMetadata, perl = T))) {
    stop("ERROR!!! options `-f [8|9|10|11|12] and -c NA` are incompatible with each other")
  }else if (length(grep('^NA$', MetadataColNamesForDge, perl = T))) {
    stop("ERROR!!! options `-f [8|9|10|11|12] and -b NA` are incompatible with each other")
  }else if (file.exists(InfileMetadata)) {
    writeLines("\nOk\n")
  }else{
    stop(paste("ERROR!!! couldnt find file:", InfileMetadata, " for option `-f [8|9|10|11|12]`", sep = "", collapse = "\n"))
  }
  
  #### Check that requested metadata columns exist in InfileMetadata
  
  writeLines(paste0("\n*** Check that requested --metadata_column_names_for_dge exist in --infile_metadata ***\n"))
  
  MetadataColNamesForDge.list  = unlist(strsplit(MetadataColNamesForDge, ","))
  if (sum(MetadataColNamesForDge.list %in% colnames(CellPropertiesFromMetadata)) == length(MetadataColNamesForDge.list)) {
    writeLines("\nOk\n")
  }else{
    stop(paste("ERROR!!! couldnt find all '", length(MetadataColNamesForDge.list), "' requested properties in --infile_metadata", sep = "", collapse = "\n"))
  }
}

#### Gene mapping into dimension reduction plot options

if ((1 %in% RequestedApplySelectedGenes == T) | 
    (2 %in% RequestedApplySelectedGenes == T) | 
    (3 %in% RequestedApplySelectedGenes == T)
    ) {
  if (length(grep('^NA$', InfileSelectedGenes, perl = T))) {
    stop("ERROR!!! options `-l [1|2|3] and -g NA` are incompatible with each other")
  }else{
    writeLines(paste0("\n*** Get list of genes to map in dimension reduction plots ***\n"))
    
    ListOfGenesForDimRedPlots<-unique(readLines(con = InfileSelectedGenes, skipNul = T))
  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO LOAD, QC AND INTEGRATE DATASETS
################################################################################################################################################
################################################################################################################################################

####################################
### Load --inputs_list
####################################
writeLines("\n*** Load --inputs_list ***\n")

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  MinioPaths <- as.list(strsplit(MinioPath, ",")[[1]])
  MinioDataPaths = data.frame(dataset_ID=rep(0, length(MinioPaths)), dataset_path=rep(0, length(MinioPaths)))
  
  for (i in seq_along(MinioPaths)) {
    MinioDataPaths[i, ] = c(basename(MinioPaths[[i]]), MinioPaths[[i]])
  }
  
  InputsTable0 <- read.table(InputsList, header = T, sep = ",", stringsAsFactors = F)
  
  MergedInputsTable <- merge(MinioDataPaths, InputsTable0, by="dataset_ID")
  MergeFilter <- c("name", "dataset_path", "dataset_type", "dataset_format", "mito_min", "mito_max", "ribo_min", "ribo_max", "ngenes_min", "ngenes_max", "nreads_min", "nreads_max")
  MergedInputsTableFiltered <- MergedInputsTable[MergeFilter]
  MergedInputsTableFilteredFinal <- MergedInputsTableFiltered[,-1]
  rownames(MergedInputsTableFilteredFinal) <- MergedInputsTableFiltered[,1]
  colnames(MergedInputsTableFilteredFinal) <-c("PathToDataset","DatasetType","DatasetFormat","MinMitoFrac","MaxMitoFrac","MinRiboFrac","MaxRiboFrac","MinNGenes","MaxNGenes","MinNReads","MaxNReads")
  
  InputsTable <- MergedInputsTableFilteredFinal
  
} else {
  InputsList<-gsub("^~/",paste0(UserHomeDirectory,"/"), InputsList)
  InputsTable<-read.table(InputsList, header = F, row.names = 1, stringsAsFactors = F)
  colnames(InputsTable)<-c("PathToDataset","DatasetType","DatasetFormat","MinMitoFrac","MaxMitoFrac","MinRiboFrac","MaxRiboFrac","MinNGenes","MaxNGenes","MinNReads","MaxNReads")
}
####################################
### Load scRNA-seq data and generate QC plots and tables
####################################
writeLines("\n*** Load scRNA-seq data and generate QC plots and tables ***\n")

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
  Dataset.SO <-paste0(dataset, ".so")
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    PathToDataset <- InputsTable[dataset,"PathToDataset"]
  } else {
    PathToDataset <- gsub("^~/",paste0(UserHomeDirectory,"/"), InputsTable[dataset,"PathToDataset"])
  }
  
  DatasetType   <- InputsTable[dataset,"DatasetType"]
  
  list_DatasetToType[[dataset]]      <- DatasetType
  list_DatasetToFormat[[dataset]]    <- InputsTable[dataset,"DatasetFormat"]
  list_TypeToDatasets[[DatasetType]] <- append(list_TypeToDatasets[[DatasetType]], dataset)
  list_MinMitoFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MinMitoFrac"])
  list_MaxMitoFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MaxMitoFrac"])
  list_MinRiboFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MinRiboFrac"])
  list_MaxRiboFrac[[dataset]]        <- as.numeric(InputsTable[dataset,"MaxRiboFrac"])
  list_MinNGenes[[dataset]]          <- as.numeric(InputsTable[dataset,"MinNGenes"])
  list_MaxNGenes[[dataset]]          <- as.numeric(InputsTable[dataset,"MaxNGenes"])
  list_MinNReads[[dataset]]          <- as.numeric(InputsTable[dataset,"MinNReads"])
  list_MaxNReads[[dataset]]          <- as.numeric(InputsTable[dataset,"MaxNReads"])
  
  if (regexpr("^MTX$|^TSV$|^HDF5$", list_DatasetToFormat[[dataset]], ignore.case = T, perl = T)[1] == 1) {
    
    ####################################
    ### Loading MTX or TSV infiles
    ####################################
    
    StopWatchStart$LoadScRNAseqData$dataset <- Sys.time()
    
    if (regexpr("^MTX$", list_DatasetToFormat[[dataset]], ignore.case = T)[1] == 1) {
      writeLines(paste0("\n*** Loading MTX infiles for ", dataset, " from: ", PathToDataset, " ***\n"))
      expression_matrix <- Read10X(data.dir = PathToDataset)
    }else if (regexpr("^TSV$", list_DatasetToFormat[[dataset]], ignore.case = T)[1] == 1) {
      writeLines(paste0("\n*** Loading matrix of genes (rows) vs. barcodes (columns) for dataset: ", dataset, " from: ", PathToDataset, " ***\n"))
      ## Note `check.names = F` is needed for both `fread` and `data.frame`
      expression_matrix <- as.matrix(data.frame(fread(PathToDataset, check.names = F), row.names=1, check.names = F))
    }else if (regexpr("^HDF5$", list_DatasetToFormat[[dataset]], ignore.case = T)[1] == 1) {
      writeLines(paste0("\n*** Loading HDF5 infile for ", dataset, " from: ", PathToDataset, " ***\n"))
      expression_matrix <- Read10X_h5(filename = PathToDataset, use.names = T, unique.features = T)
    }else{
      stop(paste0("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n"))
    }
    dim(expression_matrix)
    
    StopWatchEnd$LoadScRNAseqData$dataset  <- Sys.time()
    
    ####################################
    ### Save unfiltered data
    ####################################
    
    if (regexpr("^Y$", SaveUnFilteredData, ignore.case = T)[1] == 1) {
      writeLines("\n*** Save unfiltered data ***\n")
      
      StopWatchStart$SaveUnFilteredData$dataset  <- Sys.time()
      
      ### This is to write out the input data with no filters at all
      seurat.object.tmp  <- CreateSeuratObject(counts = expression_matrix, min.cells = 0, min.features = 0, project = paste0(PrefixOutfiles, "_", dataset))
      
      OutDirUnfilteredRaw <-paste0(Tempdir, "/UNFILTERED_DATA_MATRICES/", dataset, "/RAW")
      dir.create(file.path(OutDirUnfilteredRaw), showWarnings = F, recursive = T)
      write10xCounts(path = OutDirUnfilteredRaw, x = seurat.object.tmp@assays[["RNA"]]@counts, gene.type="Raw_Gene_Expression", overwrite=T, type="sparse", version="3")
      
      rm(seurat.object.tmp)
      
      StopWatchEnd$SaveUnFilteredData$dataset  <- Sys.time()
      
    }
    
    ####################################
    ### Create seurat object
    ####################################
    
    StopWatchStart$CreateSeuratObject$dataset  <- Sys.time()
    
    writeLines(paste0("\n*** Create seurat object for ", dataset, " ***\n"))
    seurat.object.u  <- CreateSeuratObject(counts = expression_matrix, min.cells = DefaultParameters$MinCells, min.features = list_MinNGenes[[dataset]], project = paste0(PrefixOutfiles, "_", dataset))
    
    StopWatchEnd$CreateSeuratObject$dataset  <- Sys.time()
    
    ####################################
    ### Add dataset label
    ####################################
    writeLines(paste0("\n*** Add dataset label for ", dataset, " ***\n"))
    
    StopWatchStart$AddDatasetLabel$dataset  <- Sys.time()
    
    seurat.object.u[['dataset.label']] <- dataset
    
    StopWatchEnd$AddDatasetLabel$dataset  <- Sys.time()
    
    ####################################
    ### Add dataset type label
    ####################################
    writeLines(paste0("\n*** Add dataset type label for ", dataset, " ***\n"))
    
    StopWatchStart$AddDatasetLabel$dataset_type  <- Sys.time()
    
    seurat.object.u[['dataset_type.label']] <- list_DatasetToType[[dataset]]
    
    StopWatchEnd$AddDatasetLabel$dataset_type  <- Sys.time()
    
    ####################################
    ### Get mitochondrial genes
    ####################################
    writeLines(paste0("\n*** Get  mitochondrial genes for ", dataset, " ***\n"))
    
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
    writeLines(paste0("\n*** Get ribosomal protein genes for ", dataset, " ***\n"))
    
    StopWatchStart$GetRiboGenes$dataset  <- Sys.time()
    
    riboRegExpressions<- paste(c("^MRPL", "^MRPS", "^RPL", "^RPS"), collapse = "|")
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
    writeLines(paste0("\n*** Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation for ", dataset, " ***\n"))
    
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
      
    }else{
      seurat.object.f<-subset(x = seurat.object.u, subset = 
                                nFeature_RNA >= list_MinNGenes[[dataset]]
                              & nFeature_RNA <= list_MaxNGenes[[dataset]] 
                              & nCount_RNA   >= list_MinNReads[[dataset]]
                              & nCount_RNA   <= list_MaxNReads[[dataset]]
                              & ribo.fraction >= list_MinRiboFrac[[dataset]]
                              & ribo.fraction <= list_MaxRiboFrac[[dataset]])
    }
    ### Get list of barcodes excluded by nFeature_RNA, nCount_RNA, mito.fraction or ribo.fraction
    NumberOfBarcodesExcludedByNFeature        <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nFeature_RNA  >= list_MinNGenes[[dataset]] & nFeature_RNA <= list_MaxNGenes[[dataset]]))))
    NumberOfBarcodesExcludedByNReads          <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = nCount_RNA    >= list_MinNReads[[dataset]] & nCount_RNA   <= list_MaxNReads[[dataset]]))))
    NumberOfBarcodesExcludedByMito            <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = mito.fraction >= list_MinMitoFrac[[dataset]] & mito.fraction <= list_MaxMitoFrac[[dataset]]))))
    NumberOfBarcodesExcludedByRibo            <- length(setdiff(colnames(seurat.object.u), colnames(subset(x = seurat.object.u, subset = ribo.fraction >= list_MinRiboFrac[[dataset]] & ribo.fraction <= list_MaxRiboFrac[[dataset]]))))
    
    print(paste0("NumberOfBarcodesExcludedByNFeature=", NumberOfBarcodesExcludedByNFeature))
    print(paste0("NumberOfBarcodesExcludedByNReads=", NumberOfBarcodesExcludedByNReads))
    print(paste0("NumberOfBarcodesExcludedByMito=", NumberOfBarcodesExcludedByMito))
    print(paste0("NumberOfBarcodesExcludedByRibo=", NumberOfBarcodesExcludedByRibo))
    
    StopWatchEnd$FilterCells$dataset  <- Sys.time()
    
    ### Just reporting the summary of the UNfiltered and filtered objects
    print(paste0(dataset, "  unfiltered"))
    print(seurat.object.u)
    print(paste0(dataset, "  filtered"))
    print(seurat.object.f)
    
    ####################################
    ### QC EDA violin plots
    ####################################
    writeLines(paste0("\n*** QC EDA violin plots for ", dataset, " ***\n"))
    
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
    nFeature_RNAStats.u <-paste0(c("mean = ", mean_nFeature_RNAStats.u, "\n", "median = ", median_nFeature_RNAStats.u))
    nCount_RNAStats.u   <-paste0(c("mean = ", mean_nCount_RNAStats.u,   "\n", "median = ", median_nCount_RNAStats.u))
    mito.fraction.u     <-paste0(c("mean = ", mean_mito.fraction.u,     "\n", "median = ", median_mito.fraction.u))
    ribo.fraction.u     <-paste0(c("mean = ", mean_ribo.fraction.u,     "\n", "median = ", median_ribo.fraction.u))
    
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
    nFeature_RNAStats.f <-paste0(c("mean = ", mean_nFeature_RNAStats.f, "\n", "median = ", median_nFeature_RNAStats.f))
    nCount_RNAStats.f   <-paste0(c("mean = ", mean_nCount_RNAStats.f,   "\n", "median = ", median_nCount_RNAStats.f))
    mito.fraction.f     <-paste0(c("mean = ", mean_mito.fraction.f,     "\n", "median = ", median_mito.fraction.f))
    ribo.fraction.f     <-paste0(c("mean = ", mean_ribo.fraction.f,     "\n", "median = ", median_ribo.fraction.f))
    
    ### Put QC statistics together
    nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
    nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
    mito.fraction.m.df  <-data.frame(rbind(mito.fraction.u.df,mito.fraction.f.df))
    ribo.fraction.m.df  <-data.frame(rbind(ribo.fraction.u.df,ribo.fraction.f.df))
    NumberOfCells<- list()
    NumberOfCells[["unfiltered"]] <- nrow(seurat.object.u@meta.data)
    NumberOfCells[["filtered"]]   <- nrow(seurat.object.f@meta.data)
    LabelUnfiltered    <-paste0("Before filters: No. of cells = ", NumberOfCells[["unfiltered"]])
    LabelFiltered      <-paste0("After filters:  No. of cells = ", NumberOfCells[["filtered"]])
    
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
      labs(x=paste0("No. of genes", "\n", "Filter: ", "min=", list_MinNGenes[[dataset]], " max=", list_MaxNGenes[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNFeature)) +
      annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])
    
    nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
      geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x=paste0("No. of reads", "\n", "Filter: ", "min=", list_MinNReads[[dataset]], " max=", list_MaxNReads[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByNReads)) +
      annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])
    
    mito.fraction.plot<-ggplot(data=mito.fraction.m.df, aes(x = factor(mito.fraction), y = Expression_level)) +
      geom_violin(aes(fill = factor(mito.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x="Mitochondrial genes (fraction)") +
      labs(x=paste0("Mitochondrial genes", "\n", "Filter: ", "min=", list_MinMitoFrac[[dataset]], " max=", list_MaxMitoFrac[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByMito)) +
      annotate("text", x = 1 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(mito.fraction.m.df$Expression_level)*1.1, label = mito.fraction.f, col = ColoursQCViolinPlots[[2]])
    
    ribo.fraction.plot<-ggplot(data=ribo.fraction.m.df, aes(x = factor(ribo.fraction), y = Expression_level)) +
      geom_violin(aes(fill = factor(ribo.fraction))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x=paste0("Ribosomal protein genes (fraction)", "\n", "Filter: ", "min=", list_MinRiboFrac[[dataset]], " max=", list_MaxRiboFrac[[dataset]], "\n", "Excluded cells: ", NumberOfBarcodesExcludedByRibo)) +
      annotate("text", x = 1 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(ribo.fraction.m.df$Expression_level)*1.1, label = ribo.fraction.f, col = ColoursQCViolinPlots[[2]])
    
    bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, mito.fraction.plot, ribo.fraction.plot, ncol = 4)
    
    ### Create a *pdf file with the violin ggplot's
    
    VlnPlotPdf<-paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_QC_VlnPlot.pdf")
    pdf(file=VlnPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf * 1.7, height = DefaultParameters$BaseSizeSinglePlotPdf)
    print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
    dev.off()
    
    StopWatchEnd$QCviolinplots$dataset  <- Sys.time()


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
      write.table(interactive_qc_plot_u, paste(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",dataset,"_BeforeFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
      
      # filtered
      interactive_qc_plot_f  <-data.frame(Barcodes = row.names(seurat.object.f@meta.data), Number_of_Genes = seurat.object.f@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.f@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.f@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.f@meta.data$ribo.fraction)
      interactive_qc_plot_f$Mitochondrial_Genes_Percentage <- interactive_qc_plot_f$Mitochondrial_Genes_Percentage * 100
      interactive_qc_plot_f$Ribosomal_Protein_Genes_Percentage <- interactive_qc_plot_f$Ribosomal_Protein_Genes_Percentage * 100
      colnames(interactive_qc_plot_f) <- c("Barcodes","Number of Genes","Number of Reads","Mitochondrial Genes Percentage","Ribosomal Protein Genes Percentage")
      write.table(interactive_qc_plot_f, paste(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",dataset,"_AfterFiltering.tsv",sep=""),row.names = F,sep="\t",quote = F)
      
      qc_tsv <- data.frame(NAME = row.names(seurat.object.f@meta.data), Number_of_Genes = seurat.object.f@meta.data$nFeature_RNA, Number_of_Reads = seurat.object.f@meta.data$nCount_RNA, Mitochondrial_Genes_Percentage = seurat.object.f@meta.data$mito.fraction, Ribosomal_Protein_Genes_Percentage = seurat.object.f@meta.data$ribo.fraction)
      qc_tsv$Mitochondrial_Genes_Percentage <- qc_tsv$Mitochondrial_Genes_Percentage * 100
      qc_tsv$Ribosomal_Protein_Genes_Percentage <- qc_tsv$Ribosomal_Protein_Genes_Percentage * 100
      qc_tsv_string <- sapply(qc_tsv, as.character)
      qc_tsv_string_TYPE <- rbind(data.frame(NAME = "TYPE", Number_of_Genes = "numeric", Number_of_Reads = "numeric", Mitochondrial_Genes_Percentage = "numeric", Ribosomal_Protein_Genes_Percentage = "numeric"), qc_tsv_string)
      
      qc_outfile <-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_qc/",dataset,"_qc_data.tsv")
      write.table(data.frame(qc_tsv_string_TYPE),file = qc_outfile, row.names = F, col.names = T, sep="\t", quote = F, append = T)
    } 
    
    ####################################
    ### Feature-vs-feature scatter plot
    ####################################
    writeLines(paste0("\n*** Feature-vs-feature scatter plot for ", dataset, " ***\n"))
    
    StopWatchStart$FeatureVsFeatureplot$dataset  <- Sys.time()
    
    UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                                  nGene = seurat.object.u@meta.data$nFeature_RNA,
                                  mito.fraction = seurat.object.u@meta.data$mito.fraction,
                                  filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.f))
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = ColoursQCViolinPlots[[1]])
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = F, replacement = ColoursQCViolinPlots[[2]])
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = T, replacement = 4)
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = F, replacement = 16)
    
    FeatureVsFeaturePlotPdf<-paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_NumbReadsVsNumbGenesAndMito_ScatterPlot.pdf")
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
    writeLines(paste0("\n*** Write out filter details and number of filtered cells for ", dataset, " ***\n"))
    
    StopWatchStart$OutTablesFilterDetailsAndFilteredCells$dataset  <- Sys.time()
    
    OutTableFilterDetails<-paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_FilterDetails.tsv")
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
    
    OutTableFilteredCells<-paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_NumberOfFilteredCells.tsv")
    write.table(paste("Number_of_cells_before_filters", NumberOfCells[["unfiltered"]], sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(paste("Number_of_cells_after_filters", NumberOfCells[["filtered"]],    sep = "\t", collapse = "\n"), file = OutTableFilteredCells, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    
    StopWatchEnd$OutTablesFilterDetailsAndFilteredCells$dataset  <- Sys.time()
    
    ####################################
    ### Assign data to Datasets lists
    ####################################
    writeLines(paste0("\n*** Assign data to Datasets lists: ", dataset, " ***\n"))
    
    StopWatchStart$AssignDataToDatasets  <- Sys.time()
    
    SeuratObjectsUnfiltered[[as.character(NumberOfDatasets)]]  <- seurat.object.u
    SeuratObjectsFiltered[[as.character(NumberOfDatasets)]]    <- seurat.object.f
    DatasetIds[[as.character(NumberOfDatasets)]]               <- dataset
    
    StopWatchEnd$AssignDataToDatasets  <- Sys.time()
    
    ####################################
    ### Write out QC data for each dataset
    ####################################
    writeLines("\n*** Write out QC data for each dataset ***\n")
    
    StopWatchStart$WriteOutQCData  <- Sys.time()
    
    Headers<-paste("Cell_barcode", paste(DefaultParameters$CellPropertiesToQC, sep = "", collapse = "\t") ,sep="\t")
    
    BarcodeIdsWithDatasetBeforeFilters <- unlist(x = strsplit(x = paste(dataset, colnames(SeuratObjectsUnfiltered[[NumberOfDatasets]]), sep = "_", collapse = "\n"), split = "\n"))
    OutfileQCMetadataBeforeFilters<-paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", "Before_filters_QC_metadata.tsv")
    write.table(Headers, file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(BarcodeIdsWithDatasetBeforeFilters, SeuratObjectsUnfiltered[[NumberOfDatasets]]@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataBeforeFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    
    BarcodeIdsWithDatasetAfterFilters <- unlist(x = strsplit(x = paste(dataset, colnames(SeuratObjectsFiltered[[NumberOfDatasets]]), sep = "_", collapse = "\n"), split = "\n"))
    OutfileQCMetadataAfterFilters<-paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", "After_filters_QC_metadata.tsv")
    write.table(Headers, file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(BarcodeIdsWithDatasetAfterFilters, SeuratObjectsFiltered[[NumberOfDatasets]]@meta.data[,DefaultParameters$CellPropertiesToQC]), file = OutfileQCMetadataAfterFilters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    
    StopWatchEnd$WriteOutQCData  <- Sys.time()
    
    ####################################
    ### Remove the Unfiltered seurat object
    ####################################
    writeLines(paste0("\n*** Remove the Unfiltered seurat object for ", dataset, " ***\n"))
    
    rm(seurat.object.u)
    rm(UnfilteredData.df)
    
  }else{
    stop(paste0("Unexpected type of input: ", DatasetType, "\n\nFor help type:\n\nRscript integrates_datasets_with_seurat.R -h\n\n"))
  }
}

####################################
### Get's number of rows for the legend of dimension reduction plots
####################################
if (NumberOfDatasets <= DefaultParameters$MaxNumbLabelsPerRowInLegend) {
  NumbRowsForLegendPerDataset <- 1
}else{
  NumbRowsForLegendPerDataset <- round((NumberOfDatasets / DefaultParameters$MaxNumbLabelsPerRowInLegend), digits = 0)
}

NumberOfDatasetsTypes <- length(unique(names(list_TypeToDatasets)))
if (NumberOfDatasetsTypes <= DefaultParameters$MaxNumbLabelsPerRowInLegend) {
  NumbRowsForLegendPerDatasetType <- 1
}else{
  NumbRowsForLegendPerDatasetType <- round((NumberOfDatasetsTypes / DefaultParameters$MaxNumbLabelsPerRowInLegend), digits = 0)
}

####################################
### Remove barcodes by parameter -j (if applicable)
####################################

if (regexpr("^NA$", InfileRemoveBarcodes , ignore.case = T)[1] == 1) {
  
  writeLines("\n*** No barcodes are removed by option -j ***\n")
  
}else{
  
  StopWatchStart$RemoveBarcodes <- Sys.time()
  
  writeLines("\n*** Removing barcodes by parameter -j ***\n")
  
  InfileAllBarcodesToRemove<-gsub("^~/",paste0(UserHomeDirectory,"/"), InfileRemoveBarcodes)
  AllBarcodesToRemove.tab<-read.table(InfileAllBarcodesToRemove, header = F, row.names = NULL, stringsAsFactors = FALSE)
  if (ncol(AllBarcodesToRemove.tab) == 2) {
    colnames(AllBarcodesToRemove.tab) <- c("Dataset","Barcode")
  }else{
    stop(paste0("Unexpected format in file:\n", InfileAllBarcodesToRemove, "\nfor parameter -j, 2 columns were expected but found ", ncol(AllBarcodesToRemove.tab)))
  }
  
  for (SeuratObjNumb in c(1:NumberOfDatasets)) {
    print(paste0("Removing cells from:", DatasetIds[[SeuratObjNumb]]))
    seurat.object.full <- SeuratObjectsFiltered[[SeuratObjNumb]]
    seurat.object.full
    DatasetId                      <- DatasetIds[[SeuratObjNumb]]
    ThisDatasetBarcodesToRemove    <- subset(x=AllBarcodesToRemove.tab, subset = Dataset == DatasetId)[,"Barcode"]
    ThisDatasetBarcodesToKeep.log  <- !colnames(seurat.object.full) %in% ThisDatasetBarcodesToRemove
    seurat.object.subset           <- subset(seurat.object.full, cells = colnames(seurat.object.full[,ThisDatasetBarcodesToKeep.log]))
    seurat.object.subset
    SeuratObjectsFiltered[[SeuratObjNumb]] <- seurat.object.subset
    print(paste(DatasetId, paste0("Before:", ncol(seurat.object.full)), paste0("After:", ncol(seurat.object.subset)), sep = "  ", collapse = "\n"))
  }
  
  StopWatchEnd$RemoveBarcodes <- Sys.time()
  
}

####################################
### Merge Seurat objects
####################################
writeLines("\n*** Merge Seurat objects ***\n")

StopWatchStart$MergeSeuratObjectsFiltered  <- Sys.time()

FirstSeuratObject   <- SeuratObjectsFiltered[[1]]
FirstDatasetId       <- DatasetIds[[1]]
RestOfSeuratObjectsFiltered <- SeuratObjectsFiltered[c(2:NumberOfDatasets)]
RestOfDatasetsIds    <- unlist(DatasetIds[c(2:NumberOfDatasets)])

seurat.object.merged <- merge(FirstSeuratObject, y = RestOfSeuratObjectsFiltered, 
                              add.cell.ids = DatasetIds,
                              project = PrefixOutfiles)

dataset.label.metadata<-cbind.data.frame(dataset=seurat.object.merged@meta.data$dataset.label)
rownames(dataset.label.metadata)<-colnames(seurat.object.merged)

seurat.object.merged.withmetadata <- CreateSeuratObject(seurat.object.merged@assays$RNA@counts, meta.data = dataset.label.metadata)
seurat.object.list <- SplitObject(seurat.object.merged.withmetadata, split.by = "dataset")

StopWatchEnd$MergeSeuratObjectsFiltered  <- Sys.time()

####################################
### Running SCTransform
####################################
writeLines("\n*** Running SCTransform ***\n")

StopWatchStart$SCTransform  <- Sys.time()

for (i in 1:length(seurat.object.list)) {
  
  dataset <- rownames(InputsTable)[[i]]
  print(dataset)
  seurat.object.list[[i]] <- SCTransform(seurat.object.list[[i]], verbose = T)
  
  ####################################
  ### Save filtered data
  ####################################
  writeLines("\n*** Save filtered data ***\n")
  
  StopWatchStart$SaveFilteredData$dataset  <- Sys.time()
  
  if (regexpr("^Y$", SaveFilteredData, ignore.case = T)[1] == 1) {
    
    OutDirFilteredRaw <-paste0(Tempdir, "/FILTERED_DATA_MATRICES/", dataset, "/RAW")
    dir.create(file.path(OutDirFilteredRaw), showWarnings = F, recursive = T)
    write10xCounts(path = OutDirFilteredRaw, x = seurat.object.list[[i]]@assays[["RNA"]]@counts, gene.type="Raw_Gene_Expression", overwrite=T, type="sparse", version="3")
    
    OutDirFilteredNorm<-paste0(Tempdir, "/FILTERED_DATA_MATRICES/", dataset, "/NORMALIZED_SCT")
    dir.create(file.path(OutDirFilteredNorm), showWarnings = F, recursive = T)
    write10xCounts(path = OutDirFilteredNorm, x = round(seurat.object.list[[i]]@assays[["SCT"]]@data, digits = 4), gene.type="SCTransform_Gene_Expression", overwrite=T, type="sparse", version="3")
  }
  
  StopWatchEnd$SaveFilteredData$dataset  <- Sys.time()
  
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
  seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, k.filter = KFilter, normalization.method = "SCT", anchor.features = seurat.object.integratedfeatures, verbose = T)
  
}else if (regexpr("[0-9]", ReferenceDatasets , ignore.case = T, perl = T)[1] == 1) {
  writeLines(paste0("\n*** Dataset number(s): ", ReferenceDatasets, " will be used as reference(s) ***\n"))
  seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, k.filter = KFilter, normalization.method = "SCT", anchor.features = seurat.object.integratedfeatures, reference = c(as.numeric(unlist(strsplit(ReferenceDatasets, ",")))), verbose = T)
  
}else{
  stop(paste0("Unexpected format in --reference_datasets ", ReferenceDatasets))
}

StopWatchEnd$FindIntegrationAnchors  <- Sys.time()

StopWatchStart$IntegrateData  <- Sys.time()

writeLines("\n*** Run IntegrateData() ***\n")
seurat.object.integrated <- IntegrateData(anchorset = seurat.object.anchors, normalization.method = "SCT", verbose = T)

StopWatchEnd$IntegrateData  <- Sys.time()

####################################
### Saving the R object up to dataset integration
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving the R object up to dataset integration ***\n")
  
  StopWatchStart$SaveRDSUpToIntegration  <- Sys.time()
  
  OutfileRDS<-paste0(Tempdir, "/R_OBJECTS/", PrefixOutfiles, ".", ProgramOutdir, "_integrated_object_incl_integrated_datasets.rds")
  saveRDS(seurat.object.integrated, file = OutfileRDS)
  
  StopWatchEnd$SaveRDSUpToIntegration  <- Sys.time()
  
}else{
  
  writeLines("\n*** Saving the R object up to dataset integration ***\n")
  
}

####################################
### Save normalized count matrix as loom
####################################

if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Save normalized count matrix as loom ***\n")
  
  normalized_count_matrix <- as.matrix(seurat.object.integrated@assays[["RNA"]]@data)
  
  # all genes/features in matrix
  features_tsv <- data.frame(features = rownames(normalized_count_matrix))
  features_tsv_ordered <- as.data.frame(features_tsv[mixedorder(features_tsv$features),])
  write.table(features_tsv_ordered, file=paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_raw/","features.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # generating loom file of normalized count matrix
  loom_file <- paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_normalized/","normalized_counts.loom")
  create(loom_file, normalized_count_matrix)
  
}

####################################
### Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset (columns)
####################################
writeLines("\n*** Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset (columns) ***\n")

StopWatchStart$SaveSctransformDatasetsInColumns  <- Sys.time()

seurat.object.integrated.list <- SplitObject(seurat.object.integrated, split.by = "dataset")

mat_for_correl_all_cells.df <- data.frame(row.names = rownames(seurat.object.integrated.list[[1]]@assays$SCT))
for (dataset in rownames(InputsTable)) {
  mat_for_correl_all_cells.df[[dataset]] <- rowSums(as.matrix(seurat.object.integrated.list[[dataset]]@assays$SCT[,]))
}

OutfileSCTransform <- paste(Tempdir, "/SCTRANSFORM_PSEUDO_BULK/", PrefixOutfiles, ".", ProgramOutdir, "_SCTransform_DatasetsInColumns.tsv", sep = "", collapse = "")
Headers<-paste("SCTransform", paste(colnames(mat_for_correl_all_cells.df), sep = "\t", collapse = "\t"), sep = "\t", collapse = "")
write.table(Headers, file = OutfileSCTransform, row.names = F, col.names = F, sep="\t", quote = F)
write.table(mat_for_correl_all_cells.df,  file = OutfileSCTransform, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$SaveSctransformDatasetsInColumns  <- Sys.time()

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

pdf(file=paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_PCElbowPlot_Integrated.pdf"))
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
  writeLines(paste0("\n*** Run ", DimensionReductionMethods[[dim_red_method]][["name"]], " ***\n"))
  
  StopWatchStart$DimensionReduction$dim_red_method  <- Sys.time()
  
  ### NOTES:
  ### In RunTSNE: if the datasets is small user may get error:
  ### `Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large.`
  ### User can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10
  ###
  ### Also using RunTSNE(..., check_duplicates = F) to skip cases where cells happen to have the same values after PCA reduction
  
  if (("tsne" %in% dim_red_method) & (length(colnames(seurat.object.integrated)) < DefaultParameters$MinNumberOfCellsToReducePerplexity)) {
    writeLines(paste0("\n*** Using reduced perplexity = ", DefaultParameters$ReducedPerplexity, " because found ",  length(colnames(seurat.object.integrated)), " cells", " ***\n"))
    seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, perplexity = DefaultParameters$ReducedPerplexity, check_duplicates = F)
  }else if ("tsne" %in% dim_red_method) {
    seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, check_duplicates = F)
  }else if ("umap" %in% dim_red_method) {
    seurat.object.integrated <- DimensionReductionMethods[[dim_red_method]][["run"]](object = seurat.object.integrated, dims = PcaDimsUse, umap.method = "uwot")
  }
  
  StopWatchEnd$DimensionReduction$dim_red_method  <- Sys.time()
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS ALL INTEGRATED DATASETS AS A WHOLE
################################################################################################################################################
################################################################################################################################################

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

Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "CellClusters.tsv")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)

# dataframe to merge for groups.tsv 
if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  OutfileClustersDataframe  <- data.frame(NAME = rownames(seurat.object.integrated@meta.data), Seurat_Global_Clusters = seurat.object.integrated@meta.data$seurat_clusters)
}

#
OutfileNumbClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellClusters", ".tsv")
write(x=NumberOfClusters,file = OutfileNumbClusters)
#
OutfileNumbCellsPerCluster<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerCluster.tsv")
Headers<-paste("Cluster", "Number_of_cells", sep="\t", collapse = "")
write.table(Headers,file = OutfileNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F)
write.table(table(ClusterIdent),file = OutfileNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
StopWatchEnd$AllCellClusterTables  <- Sys.time()

####################################
### Saving global cell cluster identities
####################################
writeLines("\n*** Saving global cell cluster identities ***\n")

StopWatchStart$SaveGlobalCellClusterIdentities  <- Sys.time()

seurat.object.integrated$GlobalCellClusterIdentities <- Idents(object = seurat.object.integrated)

StopWatchEnd$SaveGlobalCellClusterIdentities  <- Sys.time()

####################################
### Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset cluster (columns)
####################################
writeLines("\n*** Generate a matrix with each gene (rows) marginals of SCTransform values from all cells for each dataset cluster (columns) ***\n")

StopWatchStart$SaveSctransformClustersInColumns  <- Sys.time()

seurat.object.integrated.list <- SplitObject(seurat.object.integrated, split.by = "dataset")

mat_for_correl_each_cluster.df <- data.frame(row.names = rownames(seurat.object.integrated.list[[1]]@assays$SCT))
for (DatasetId in rownames(InputsTable)) {
  if (exists(x = "seurat.object.each_dataset") == T) {
    rm(seurat.object.each_dataset)
  }
  seurat.object.each_dataset <- subset(x = seurat.object.integrated, subset = dataset == DatasetId)
  
  for (cluster_number in sort(unique(seurat.object.integrated$GlobalCellClusterIdentities))) {
    dataset_cluster <- (paste(DatasetId, cluster_number, sep = "_c", collapse = ""))
    res <- try(subset(x = seurat.object.each_dataset, subset = seurat_clusters == cluster_number), silent = TRUE)
    if (class(res) == "try-error") {
      # mat_for_correl_each_cluster.df[[dataset_cluster]] <- NA
    }else{
      if (exists(x = "seurat.object.each_dataset.each_cluster") == T) {
        rm(seurat.object.each_dataset.each_cluster)
      }
      seurat.object.each_dataset.each_cluster <- subset(x = seurat.object.each_dataset, subset = seurat_clusters == cluster_number)
      mat_for_correl_each_cluster.df[[dataset_cluster]] <- rowSums(as.matrix(seurat.object.each_dataset.each_cluster@assays$SCT[,]))
    }
  }
}

OutfileSCTransform <- paste0(Tempdir, "/SCTRANSFORM_PSEUDO_BULK/", PrefixOutfiles, ".", ProgramOutdir, "_SCTransform_DatasetClustersInColumns.tsv")
Headers<-paste("SCTransform", paste(colnames(mat_for_correl_each_cluster.df), sep = "\t", collapse = "\t"), sep = "\t", collapse = "")
write.table(Headers, file = OutfileSCTransform, row.names = F, col.names = F, sep="\t", quote = F)
write.table(mat_for_correl_each_cluster.df,  file = OutfileSCTransform, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$SaveSctransformClustersInColumns  <- Sys.time()

####################################
### FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data to colour by: a) cell clusters, b) dataset, c) dataset type, d) selected genes and e) metadata
####################################
writeLines("\n*** FOR EACH DIMENSION REDUCTION TYPE colour plots using integrated data to colour by: a) cell clusters, b) dataset, c) dataset type, d) selected genes and e) metadata  ***\n")

for (dim_red_method in names(DimensionReductionMethods)) {
  
  ####################################
  ### Colour dimension reduction plot by global cell clusters using integrated data
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by global cell clusters ***\n"))
  
  StopWatchStart$DimRedOPlotColourByCellCluster$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 10)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
  IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir,  "_GlobalClustering_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$DimRedOPlotColourByCellCluster$dim_red_method  <- Sys.time()
  
  ####################################
  ### Write out dimension reduction plot coordinates using integrated data
  ####################################
  writeLines(paste0("\n*** Write out ", DimensionReductionMethods[[dim_red_method]][["name"]], " coordinates ***\n"))
  
  StopWatchStart$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  OutfileCoordinates<-paste0(Tempdir, "/DIMENSION_REDUCTION_COORDINATE_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_Coordinates.tsv")
  Headers<-paste("Barcode",paste(colnames(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings), sep="", collapse="\t"), sep="\t", collapse = "\t")
  write.table(Headers,file = OutfileCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    OutfileCoordinatesCWL<-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_coordinates/",DimensionReductionMethods[[dim_red_method]][["name"]], "Coordinates.tsv")
    write.table(Headers,file = OutfileCoordinatesCWL, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(seurat.object.integrated@reductions[[dim_red_method]]@cell.embeddings, file = OutfileCoordinatesCWL,  row.names = T, col.names = F, sep="\t", quote = F, append = T)
  }
  
  StopWatchEnd$DimensionReductionWriteCoords$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by dataset using integrated data
  ####################################
  writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by dataset ***\n"))
  
  StopWatchStart$DimRedPlotsByDataset$dim_red_method  <- Sys.time()
  
  plots <- DimPlot(seurat.object.integrated, group.by = c("dataset"), combine = F, reduction = dim_red_method)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = NumbRowsForLegendPerDataset, byrow = T, override.aes = list(size = 2))))
  IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByDataset.pdf")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  ### This is the same plot as above, but with only two dataset labels per row in the legend to make sure that all text is printed out
  plots <- DimPlot(seurat.object.integrated, group.by = c("dataset"), combine = F, reduction = dim_red_method)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(ncol = 2, byrow = T, override.aes = list(size = 1))))
  IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByDataset_FullLegend.pdf")
  pdf(file=IntegratedDimRedPlotPdf, width = 7, height = (7 + (0.1 * NumberOfDatasets)))
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$DimRedPlotsByDataset$dim_red_method  <- Sys.time()
  
  ####################################
  ### Colour dimension reduction plots by dataset type using integrated data
  ####################################
  
  if (NumberOfDatasetsTypes >= 1) {
    writeLines("\n*** Get dataset_type assignments replacing dataset ids using list_DatasetToType() ***\n")
    
    StopWatchStart$PrepareEachDatasetTypeForDimRedPlotByDatasetType <- Sys.time()
    
    seurat.object.integrated$dataset_type<-seurat.object.integrated$dataset
    for (dataset in rownames(InputsTable)) {
      seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
    }
    DatasetTypes <- unique(seurat.object.integrated$dataset_type)
    
    StopWatchEnd$PrepareEachDatasetTypeForDimRedPlotByDatasetType <- Sys.time()
    
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by dataset type ***\n"))
    
    # switch identity to dataset type global cluster
    Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset_type
    
    head(seurat.object.integrated@meta.data)
    
    NumberOfDatasetTypes <- 0
    for (dataset_type in DatasetTypes) {
      NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
      print(NumberOfDatasetTypes)
      
      if (exists(x = "seurat.object.each_dataset_type") == T) {
        rm(seurat.object.each_dataset_type)
      }
      seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
      
      writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " coloured by dataset_type ***\n"))
      
      StopWatchStart$DimRedOPlotColourByDatasetType$dim_red_method  <- Sys.time()
      
      plots <- DimPlot(seurat.object.integrated, group.by = c("dataset_type"), combine = F, reduction = dim_red_method, label = F, label.size = 10)
      plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = NumbRowsForLegendPerDatasetType, byrow = T, override.aes = list(size = 2))))
      IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByDatasetType.pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
      print(CombinePlots(plots))
      dev.off()
      
      StopWatchEnd$DimRedOPlotColourByDatasetType$dim_red_method  <- Sys.time()
    }
  }
  
  ####################################
  ### Colour dimension reduction plots for all cells by selected genes (option -g)
  ####################################
  
  if (regexpr("^NA$", InfileSelectedGenes, ignore.case = T)[1] == 1) {
    print("No selected genes for dimension reduction plots")
    
  }else if (1 %in% RequestedApplySelectedGenes == T) {
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " for all cells by selected genes (option -g) ***\n"))
    
    StopWatchStart$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
    
    ### Making a new Seurat object `seurat.object.integrated.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
    ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
    seurat.object.integrated.sa <- CreateSeuratObject(seurat.object.integrated@assays$RNA@data)
    seurat.object.integrated.sa@reductions$umap <- seurat.object.integrated@reductions$umap
    seurat.object.integrated.sa@reductions$tsne <- seurat.object.integrated@reductions$tsne
    
    OutDirThisOption <- paste0(Tempdir, "/SELECTED_GENE_DIMENSION_REDUCTION_PLOTS/ALL_CELLS/", DimensionReductionMethods[[dim_red_method]][["name"]]) 
    dir.create(path = OutDirThisOption, recursive = T)
    sapply(ListOfGenesForDimRedPlots, FUN=function(eachGene) {
      OutfilePathAndName <- paste0(OutDirThisOption, "/", eachGene, "_all_cells_", DimensionReductionMethods[[dim_red_method]][["name"]], ".pdf")
      pdf(file=OutfilePathAndName, width=DefaultParameters$BaseSizeMultiplePlotPdfWidth, height=DefaultParameters$BaseSizeMultiplePlotPdfHeight)
      print(FeaturePlot(object = seurat.object.integrated.sa, features = eachGene, cols = c("lightgrey", "blue"),
                        reduction = dim_red_method, order = T, slot = "data", pt.size = 0.3, min.cutoff = "q0.1", max.cutoff = "q90"))
      dev.off()
    })
    
    StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
    
  }
  
  ####################################
  ### Colour dimension reduction plots by -infile_metadata using integrated data
  ####################################
  
  if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
    print("No metadata will be used for dimension reduction plots")
  }else{
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot by -infile_metadata ***\n"))
    
    StopWatchStart$DimRedPlotsColuredByMetadata$dim_red_method  <- Sys.time()
    
    # Note: needs data.frame(CellPropertiesFromMetadata) preloaded
    
    # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same dataset)
    # so that barcodes from --infile_metadata and --input can match each other
    rownames(CellPropertiesFromMetadata)<-gsub(x =rownames(CellPropertiesFromMetadata), pattern = "-[0-9]+$", perl = T, replacement = "")
    seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = CellPropertiesFromMetadata)
    
    # Generating outfile
    # Note DimPlot() takes the entire current device (pdf)
    # even if using layout(matrix(...)) or  par(mfrow=())
    # Thus each property plot is written to a separate page of a single *pdf outfile
    for (property in colnames(CellPropertiesFromMetadata)) {
      IntegratedDimRedPlotPdf <- paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBy_", property, ".pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(CellPropertiesFromMetadata[,property])))))
      
      if ( (sum(CellPropertiesFromMetadata[,property] %in% 0:1 == T)) == (nrow(CellPropertiesFromMetadata)) ) { ## is binary
        CellsToHighlight <- rownames(CellPropertiesFromMetadata)[CellPropertiesFromMetadata[,property]==1]
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
}

####################################
### Get average gene expression for each global cluster
####################################
writeLines("\n*** Get average gene expression for each global cluster ***\n")

StopWatchStart$AverageGeneExpressionGlobalClusters  <- Sys.time()

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$seurat_clusters

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)
#
OutfileClusterAveragesRNA<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllDatasets_AverageGeneExpression_RNA.tsv")
OutfileClusterAveragesSCT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllDatasets_AverageGeneExpression_SCT.tsv")
OutfileClusterAveragesINT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "AllDatasets_AverageGeneExpression_INT.tsv")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION", paste("r", Resolution, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$AverageGeneExpressionGlobalClusters  <- Sys.time()

####################################
### Finding differentially expressed genes (1): using global cell clusers, compares each cell cluster vs. the rest of cells
####################################

if (1 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes (1): using global cell clusers, compares each cell cluster vs. the rest of cells ***\n")
  
  StopWatchStart$FindDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()
  
  print(paste0("Number of clusters = ", NumberOfClusters))
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  seurat.object.integrated.markers <- FindAllMarkers(object = seurat.object.integrated, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
  SimplifiedDiffExprGenes.df <- seurat.object.integrated.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
  OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "MarkersPerCluster.tsv")
  write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
  
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    top_6_genes_by_cluster<-(seurat.object.integrated.markers %>% group_by(cluster) %>% top_n(6, avg_logFC))
    globalMarkersFile <- top_6_genes_by_cluster[,c("gene","cluster","p_val","avg_logFC")]
    write.table(globalMarkersFile, paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_markers/","TopTwoMarkersPerCluster.tsv"), row.names = F, sep="\t", quote = F)
  } 
  
  StopWatchEnd$FindDiffMarkersGlobalClustersVsRestOfCells  <- Sys.time()
  
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET USING GLOBAL CELL CLUSTERS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET uses integrated/normalized counts and generates dimension reduction plots for each dataset coloured by: a) QC, b) global clusters, c) selected genes and d) metadata
####################################
writeLines("\n*** FOR EACH DATASET uses integrated/normalized counts and generates dimension reduction plots for each dataset coloured by: a) QC, b) global clusters, c) selected genes and d) metadata ***\n")

####################################
### Colour dimension reduction plots for each dataset by QC (nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction)
####################################
writeLines("\n*** Colour dimension reduction plots for each dataset by QC (nFeature_RNA, nCount_RNA, mito.fraction and ribo.fraction) ***\n")

StopWatchStart$DimRedPlotsQC  <- Sys.time()

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  
  ### Note: need to AddMetaData() mito.fraction and ribo.fraction from OutfileQCMetadata to generate FeaturePlot()
  ### because they are not inherited in seurat.object.integrated by IntegrateData()
  InfileQCMetadataAfterFilters<-paste0(Tempdir, "/QC_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", "After_filters_QC_metadata.tsv")
  
  QCMetadata <- data.frame(read.table(InfileQCMetadataAfterFilters, header = T, row.names = 1))
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = QCMetadata)
  
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_dataset") == T) {
    rm(seurat.object.each_dataset)
  }
  seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_dataset)
  
  for (dim_red_method in names(DimensionReductionMethods)) {
    pdf(file=paste0(Tempdir, "/QC_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset,"_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByQC.pdf"), width = DefaultParameters$BaseSizeSinglePlotPdf * 1.5, height = DefaultParameters$BaseSizeSinglePlotPdf * 1.5)
    print(FeaturePlot(object = seurat.object.each_dataset, label = F, order = T, features = DefaultParameters$CellPropertiesToQC, cols = c("lightgrey", "blue"), reduction = dim_red_method, ncol = 2, pt.size = 1.5))
    dev.off()
  }
}

StopWatchEnd$DimRedPlotsQC  <- Sys.time()

####################################
### Merge dataset_id and global cluster identities to make dataset-specific identities
####################################
writeLines("\n*** Merge dataset_id and global cluster identities to make dataset-specific identities  ***\n")

StopWatchStart$MergeDatasetAndClusterIds  <- Sys.time()

EachDatasetGlobalCellClusters <- unlist(x = strsplit(x = paste(seurat.object.integrated$dataset, seurat.object.integrated$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetGlobalCellClusters, col.name = "EachDatasetGlobalCellClusters")

# switch the identity class of all cells to reflect dataset-specific identities
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetGlobalCellClusters

StopWatchEnd$MergeDatasetAndClusterIds  <- Sys.time()

####################################
### Write out number and fraction of cells per cluster, per dataset
####################################
writeLines("\n*** Write out number and fraction of cells per cluster, per dataset  ***\n")

StopWatchStart$SaveFractionsOfCellsPerClusterPerDataset  <- Sys.time()

OutfileNumbCellsPerClusterPerDataset<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerClusterPerDataset.tsv")
OutfileFracCellsPerClusterPerDataset<-paste0(Tempdir, "/CELL_FRACTIONS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "FracCellsPerClusterPerDataset.tsv")

seurat.object.integrated_EachDatasetGlobalCellClusters.df <- as.data.frame(table(seurat.object.integrated$EachDatasetGlobalCellClusters))
rownames(seurat.object.integrated_EachDatasetGlobalCellClusters.df) <- seurat.object.integrated_EachDatasetGlobalCellClusters.df[,1]

freq_mat <- sapply(rownames(InputsTable), function(dataset) {
  freq_clust <- sapply(sort(unique(seurat.object.integrated$seurat_clusters)), function(cluster) {
    dataset_cluster <- paste0(dataset, "_c", cluster)
    if ((dataset_cluster %in% row.names(seurat.object.integrated_EachDatasetGlobalCellClusters.df)) == TRUE)  {
      freq <- seurat.object.integrated_EachDatasetGlobalCellClusters.df[dataset_cluster,"Freq"]
    }else{
      freq <- 0
    }
    freq ## returns freq for the second sapply loop
  })
  freq_clust ## returns freq_clust for the first sapply loop
})
rownames(freq_mat) <- sort(unique(seurat.object.integrated$seurat_clusters))
frac_mat <- t(round(t(freq_mat) / colSums(freq_mat), 4)) ### note it's necessary to t(freq_ma) otherwise the colSums() are used incorrectly

Headers<-paste("cluster", paste(rownames(InputsTable) , sep="", collapse = "\t"), sep = "\t", collapse = "\t")

write.table(Headers, file = OutfileNumbCellsPerClusterPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(freq_mat, file = OutfileNumbCellsPerClusterPerDataset, row.names = T, col.names = F, sep="\t", quote = F, append = T)

write.table(Headers, file = OutfileFracCellsPerClusterPerDataset, row.names = F, col.names = F, sep="\t", quote = F)
write.table(frac_mat , file = OutfileFracCellsPerClusterPerDataset, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$SaveFractionsOfCellsPerClusterPerDataset  <- Sys.time()

####################################
### Colour dimension reduction plots for each dataset based on global clusters, selected genes and metadata
####################################

writeLines("\n*** Colour dimension reduction plots for each dataset based on global clusters, selected genes and metadata ***\n")

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset

NumberOfDatasets <- 0
for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_dataset") == T) {
    rm(seurat.object.each_dataset)
  }
  seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_dataset)
  
  for (dim_red_method in names(DimensionReductionMethods)) {
    
    ####################################
    ### Colour dimension reduction plots for each dataset by global cell clusters
    ####################################
    writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for each dataset by global cell clusters ***\n"))
    
    StopWatchStart$DimRedPlotsByDatasetIntegratedCellClusters$dataset$dim_red_method  <- Sys.time()
    
    plots <- DimPlot(seurat.object.each_dataset, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
    plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
    IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf")
    pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
    print(CombinePlots(plots))
    dev.off()
    
    StopWatchEnd$DimRedPlotsByDatasetIntegratedCellClusters$dataset$dim_red_method  <- Sys.time()
    
    ####################################
    ### Colour dimension reduction plots for each dataset by selected genes (option -g)
    ####################################
    
    if (regexpr("^NA$", InfileSelectedGenes, ignore.case = T)[1] == 1) {
      print("No selected genes for each dataset dimension reduction plots")
      
    }else if (2 %in% RequestedApplySelectedGenes == T) {
      writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for each dataset by selected genes (option -g) ***\n"))
      
      StopWatchStart$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
      
      ### Making a new Seurat object `seurat.object.each_dataset.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
      ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
      seurat.object.each_dataset.sa <- CreateSeuratObject(seurat.object.each_dataset@assays$RNA@data)
      seurat.object.each_dataset.sa@reductions$umap <- seurat.object.each_dataset@reductions$umap
      seurat.object.each_dataset.sa@reductions$tsne <- seurat.object.each_dataset@reductions$tsne
      
      OutDirThisOption <- paste0(Tempdir, "/SELECTED_GENE_DIMENSION_REDUCTION_PLOTS/DATASETS/", DimensionReductionMethods[[dim_red_method]][["name"]]) 
      dir.create(path = OutDirThisOption, recursive = T)
      sapply(ListOfGenesForDimRedPlots, FUN=function(eachGene) {
        OutfilePathAndName <- paste0(OutDirThisOption, "/", eachGene, "_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], ".pdf")
        pdf(file=OutfilePathAndName, width=DefaultParameters$BaseSizeMultiplePlotPdfWidth, height=DefaultParameters$BaseSizeMultiplePlotPdfHeight)
        print(FeaturePlot(object = seurat.object.each_dataset.sa, features = eachGene, cols = c("lightgrey", "blue"),
                          reduction = dim_red_method, order = T, slot = "data", pt.size = 0.3, min.cutoff = "q0.1", max.cutoff = "q90"))
        dev.off()
      })
      
      StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
      
      rm(seurat.object.each_dataset.sa)
      
    }
    
    ####################################
    ### Colour dimension reduction plots for each dataset by -infile_metadata using integrated data
    ####################################
    
    if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
      print("No metadata will be used for dimension reduction plots")
    }else{
      writeLines("\n*** Colour dimension reduction plots for each dataset by -infile_metadata using integrated data ***\n")
      
      # Note: needs data.frame(CellPropertiesFromMetadata) preloaded
      
      # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same dataset)
      # so that barcodes from --infile_metadata and --input can match each other
      rownames(CellPropertiesFromMetadata)<-gsub(x =rownames(CellPropertiesFromMetadata), pattern = "-[0-9]+$", perl = T, replacement = "")
      seurat.object.each_dataset <- AddMetaData(object = seurat.object.each_dataset, metadata = CellPropertiesFromMetadata)
      
      StopWatchStart$DimRedPlotsByMetadata$dataset$dim_red_method  <- Sys.time()
      
      # Generating outfile
      # Note DimPlot() takes the entire current device (pdf)
      # even if using layout(matrix(...)) or  par(mfrow=())
      # Thus each property plot is written to a separate page of a single *pdf outfile
      for (property in colnames(CellPropertiesFromMetadata)) {
        
        if (sum(is.na(seurat.object.each_dataset[[property]]==T)) == ncol(seurat.object.each_dataset)) {
          
          print(paste0("No extra barcode-attributes '", property, "' found for dataset ", dataset))
          
        }else{
          
          IntegratedDimRedPlotPdf <- paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBy_", property, ".pdf")
          pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(CellPropertiesFromMetadata[,property])))))
          
          if ( (sum(CellPropertiesFromMetadata[,property] %in% 0:1 == T)) == (nrow(CellPropertiesFromMetadata)) ) { ## is binary
            CellsToHighlight <- rownames(CellPropertiesFromMetadata)[CellPropertiesFromMetadata[,property]==1]
            plots <- DimPlot(seurat.object.each_dataset, group.by = property, combine = F, reduction = dim_red_method, label = F, cells.highlight = CellsToHighlight)
            plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
            print(CombinePlots(plots))
            
          }else{
            plots <- DimPlot(seurat.object.each_dataset, group.by = property, combine = F, reduction = dim_red_method, label = T)
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
### Get average gene expression for each dataset based on global clusters
####################################
writeLines("\n*** Get average gene expression for each dataset based on global clusters ***\n")

StopWatchStart$AverageGeneExpressionEachDatasetGlobalClustersPerDataset  <- Sys.time()

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerDataset_AverageGeneExpression_RNA.tsv")
OutfileClusterAveragesSCT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerDataset_AverageGeneExpression_SCT.tsv")
OutfileClusterAveragesINT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerDataset_AverageGeneExpression_INT.tsv")
#
Headers<-paste("AVERAGE_GENE_EXPRESSION", paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
#
write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
#
write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$AverageGeneExpressionEachDatasetGlobalClustersPerDataset  <- Sys.time()

####################################
### Finding differentially expressed genes (2): using global cell clusers, for each dataset, compares each cell cluster vs. the rest of cells
####################################

if (2 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes (2): using global cell clusers, for each dataset, compares each cell cluster vs. the rest of cells ***\n")
  
  print(paste0("Number of clusters = ", NumberOfClusters))
  
  ####################################
  ### Loops each dataset
  ####################################
  NumberOfDatasets <- 0
  for (dataset in rownames(InputsTable)) {
    
    StopWatchStart$FindDiffMarkersEachDatasetGlobalClustersVsRestOfCells$dataset  <- Sys.time()
    
    NumberOfDatasets <- NumberOfDatasets + 1
    print(NumberOfDatasets)
    
    if (exists(x = "seurat.object.each_dataset") == T) {
      rm(seurat.object.each_dataset)
    }
    
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset
    seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$seurat_clusters
    Idents(object = seurat.object.each_dataset) <- seurat.object.each_dataset$seurat_clusters
    print(seurat.object.each_dataset)
    
    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset@meta.data))
    seurat.object.each_dataset.markers <- FindAllMarkers(object = seurat.object.each_dataset, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    SimplifiedDiffExprGenes.df <- seurat.object.each_dataset.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
    OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", dataset, "_MarkersPerCluster.tsv")
    write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
    
    StopWatchEnd$FindDiffMarkersEachDatasetGlobalClustersVsRestOfCells$dataset  <- Sys.time()
    
  }
}

####################################
### Finding differentially expressed genes (3): using global cell clusers, for each dataset, compares each cell cluster vs. the same cluster from other datasets
####################################

if (3 %in% RequestedDiffGeneExprComparisons == T) {
  
  writeLines("\n*** Finding differentially expressed genes (3): using global cell clusers, for each dataset, compares each cell cluster vs. the same cluster from other datasets ***\n")
  
  StopWatchStart$FindDiffMarkersEachDatasetGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()
  
  print(paste0("Number of clusters = ", NumberOfClusters))
  
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
  
  OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "MarkersPerDatasetEquivalentClusters.tsv")
  HeadersOrder <- paste("cluster1", "cluster2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
  write.table(HeadersOrder,file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="", quote = F)
  
  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$EachDatasetGlobalCellClusters
  
  for (cluster in unique(seurat.object.integrated$seurat_clusters)) {
    for (dataset1 in rownames(InputsTable)) {
      for (dataset2 in rownames(InputsTable)) {
        Cluster1 <- paste(dataset1, cluster, sep = "_c")
        Cluster2 <- paste(dataset2, cluster, sep = "_c")
        if (Cluster1 == Cluster2) {
          ### Skip
        }else if ((sum(seurat.object.integrated$EachDatasetGlobalCellClusters == Cluster1) >= 3) & ((sum(seurat.object.integrated$EachDatasetGlobalCellClusters == Cluster2) >= 3))) {
          print (paste0(Cluster1, " vs. ", Cluster2))
          seurat.object.integrated.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.integrated, only.pos = F, ident.1 = Cluster1, ident.2 = Cluster2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
          seurat.object.integrated.each_equivalent_cluster.markers$cluster1 <- Cluster1
          seurat.object.integrated.each_equivalent_cluster.markers$cluster2 <- Cluster2
          seurat.object.integrated.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.integrated.each_equivalent_cluster.markers)
          write.table(seurat.object.integrated.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="\t", quote = F, append = T)
        }else{
          print(paste0("Skip cluster ", Cluster1, " vs. ", Cluster2, " because there were not >= 3 cells in at least one of them"))
        }
      }
    }
  }
  
  StopWatchEnd$FindDiffMarkersEachDatasetGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()
  
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET USING RE-CLUSTERED CELLS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET uses integrated counts and RE-CLUSTER cells. Then colour dimension reduction plots and get DGE (if requested)
####################################

if (2 %in% RequestedClusteringInputs == T) {
  
  writeLines("\n*** FOR EACH DATASET uses integrated counts and RE-CLUSTER cells. Then colour dimension reduction plots and get DGE (if requested) ***\n")
  
  ####################################
  ### Prepares outfile headers
  ####################################
  
  # switch the identity class of all cells to reflect dataset
  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset
  
  # Headers for dataset-specific cell clusters
  Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
  OutfileEachDatasetClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "CellClusters.tsv")
  write.table(Headers,file = OutfileEachDatasetClusters, row.names = F, col.names = F, sep="\t", quote = F)
  #
  Headers<-paste("Dataset", "Number_of_clusters", sep="\t", collapse = "")
  OutfileEachDatasetNumbClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "NumbCellClusters", ".tsv")
  write(x=NumberOfClusters,file = OutfileNumbClusters)
  #
  Headers<-paste("Cluster", "Number_of_cells", sep="\t", collapse = "")
  OutfileEachDatasetNumbCellsPerCluster<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "NumbCellsPerCluster.tsv")
  write.table(Headers,file = OutfileEachDatasetNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F)
  
  ####################################
  ### Loops each dataset
  ####################################
  NumberOfDatasets <- 0
  for (dataset in rownames(InputsTable)) {
    NumberOfDatasets <- NumberOfDatasets + 1
    print(NumberOfDatasets)
    
    if (exists(x = "seurat.object.each_dataset") == T) {
      rm(seurat.object.each_dataset)
    }
    seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
    print(seurat.object.each_dataset)
    
    StopWatchStart$ClusterEachDatasetCellsFromInteg$dataset  <- Sys.time()
    
    ####################################
    ### Re-cluster each dataset cells using integrated data
    ####################################
    
    options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
    seurat.object.each_dataset <- FindNeighbors(object = seurat.object.each_dataset, dims = PcaDimsUse)
    seurat.object.each_dataset <- FindClusters(object = seurat.object.each_dataset, resolution = Resolution)
    ClustersThisDataset <- unlist(x = strsplit(x = paste(dataset, seurat.object.each_dataset@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
    
    StopWatchEnd$ClusterEachDatasetCellsFromInteg$dataset  <- Sys.time()
    
    ####################################
    ### Write out each dataset cell re-clusters
    ####################################
    writeLines("\n*** Write out each dataset cell re-clusters ***\n")
    
    StopWatchStart$WriteClustersEachDatasetCellClusterTables$dataset <- Sys.time()
    
    ### Cell cluster identities
    write.table(paste(colnames(seurat.object.each_dataset), ClustersThisDataset, sep = "\t", collapse = "\n"),
                file = OutfileEachDatasetClusters, row.names = F, col.names = F, quote = F, append = T)
    
    ### Number of clusters per dataset
    CellNames<-rownames(seurat.object.each_dataset@meta.data)
    ClusterIdent <-seurat.object.each_dataset@meta.data$seurat_clusters
    NumberOfClusters<-length(unique(ClusterIdent))
    write(x=paste(dataset, NumberOfClusters, sep = "\t", collapse = ""), file = OutfileEachDatasetNumbClusters, append = T)
    
    ### Number of cells per cluster, per dataset
    NumberOfCellPerClusterPerDataset<-table(ClusterIdent)
    rownames(NumberOfCellPerClusterPerDataset) <- paste(dataset, rownames(NumberOfCellPerClusterPerDataset), sep = "_c")
    write.table(NumberOfCellPerClusterPerDataset, file = OutfileEachDatasetNumbCellsPerCluster, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    
    StopWatchEnd$WriteClustersEachDatasetCellClusterTables$dataset <- Sys.time()
    
    ####################################
    ### Colour dimension reduction plots for each dataset by re-clustered cells
    ####################################
    writeLines("\n*** Colour dimension reduction plots for each dataset by re-clustered cells ***\n")
    
    for (dim_red_method in names(DimensionReductionMethods)) {
      
      StopWatchStart$DimRedPlotsByEachDatasetReclusteredCellCluster$dataset$dim_red_method  <- Sys.time()
      
      plots <- DimPlot(seurat.object.each_dataset, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
      plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
      IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", dataset, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
      print(CombinePlots(plots))
      dev.off()
      
      StopWatchEnd$DimRedPlotsByEachDatasetReclusteredCellCluster$dataset$dim_red_method <- Sys.time()
      
    }
    
    ####################################
    ### Finding differentially expressed genes (6): using re-clustered cells, for each dataset, compares each cell cluster vs. the rest of cells
    ####################################
    if (6 %in% RequestedDiffGeneExprComparisons == T) {
      
      writeLines("\n*** Finding differentially expressed genes (6): using re-clustered cells, for each dataset, compares each cell cluster vs. the rest of cells ***\n")
      
      print(paste0("Number of clusters = ", NumberOfClusters))
      
      StopWatchStart$FindDiffMarkersReclusteredVsRestOfCells$dataset  <- Sys.time()
      
      FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset@meta.data))
      seurat.object.each_dataset.markers <- FindAllMarkers(object = seurat.object.each_dataset, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
      SimplifiedDiffExprGenes.df <- seurat.object.each_dataset.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
      OutfileDiffGeneExpression <- paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", dataset, "_MarkersPerCluster.tsv")
      write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
      
      StopWatchEnd$FindDiffMarkersReclusteredVsRestOfCells$dataset  <- Sys.time()
      
    }
  }
  
  ####################################
  ### Load each dataset re-cluster assignments
  ####################################
  writeLines("\n*** Load each dataset re-cluster assignments ***\n")
  
  EachDatasetCellReClusters <- data.frame(read.table(OutfileEachDatasetClusters, header = T, row.names = 1, sep = "\t"))
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetCellReClusters, col.name = "EachDatasetCellReClusters")
  
  # switch the identity class of all cells to reflect each dataset clusters
  Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetCellReClusters
  
  ####################################
  ### Get average gene expression for each dataset re-clustered clusters
  ####################################
  writeLines("\n*** Get average gene expression for each dataset re-clustered clusters ***\n")
  
  StopWatchStart$AverageGeneExpressionEachDatasetReclustered$dataset  <- Sys.time()
  
  cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)
  
  OutfileClusterAveragesRNA<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "PerDataset_AverageGeneExpression_RNA.tsv")
  OutfileClusterAveragesSCT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "PerDataset_AverageGeneExpression_SCT.tsv")
  OutfileClusterAveragesINT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetReclustered_", "PerDataset_AverageGeneExpression_INT.tsv")
  #
  Headers<-paste("AVERAGE_GENE_EXPRESSION", paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
  #
  write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  #
  write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  #
  write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$AverageGeneExpressionEachDatasetReclustered$dataset  <- Sys.time()
  
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET TYPE USING GLOBAL CELL CLUSTERS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET TYPE uses integrated/normalized counts and generates dimension reduction plots coloured by: a) dataset type, b) global clusters, c) selected genes and d) metadata
####################################
if (NumberOfDatasetsTypes >= 1) {
  writeLines("\n*** FOR EACH DATASET TYPE uses integrated/normalized counts and generates dimension reduction plots coloured by: a) dataset type, b) global clusters, c) selected genes and d) metadata ***\n")
  
  ####################################
  ### Prepares each dataset type object
  ####################################
  
  StopWatchStart$PrepareEachDatasetTypeForDimRedPlotsByVariousAttributes <- Sys.time()
  
  seurat.object.integrated$dataset_type<-seurat.object.integrated$dataset
  for (dataset in rownames(InputsTable)) {
    seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
  }
  DatasetTypes <- unique(seurat.object.integrated$dataset_type)
  
  StopWatchEnd$PrepareEachDatasetTypeForDimRedPlotsByVariousAttributes <- Sys.time()
  
  ####################################
  ### Write out each dataset type global cluster assignments
  ####################################
  writeLines("\n*** Write out each dataset type global cluster assignments ***\n")
  
  StopWatchStart$WriteOutDatasetTypeClusterAssignments  <- Sys.time()
  
  # Headers for dataset-type-specific cell clusters
  Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
  OutfileEachDatasetTypeGlobalClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "EachDatasetType_CellClusters.tsv")
  GlobalClustersByDatasetType<- unlist(x = strsplit(x = paste(seurat.object.integrated@meta.data$dataset_type, seurat.object.integrated@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
  write.table(Headers,file = OutfileEachDatasetTypeGlobalClusters, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(paste(colnames(seurat.object.integrated), GlobalClustersByDatasetType, sep = "\t", collapse = "\n"),
              file = OutfileEachDatasetTypeGlobalClusters, row.names = F, col.names = F, quote = F, append = T)
  
  # dataframe to merge for groups.tsv 
  if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
    OutfileEachDatasetTypeGlobalClustersDataframe  <- data.frame(NAME = colnames(seurat.object.integrated), Seurat_Datasets_Global_Clusters = GlobalClustersByDatasetType)
    
    groupsMergedDataframe <- merge(OutfileClustersDataframe, OutfileEachDatasetTypeGlobalClustersDataframe, by="NAME")
    groupsMergedDataframeString <- sapply(groupsMergedDataframe, as.character)
    groupsMergedDataframeStringTYPE <- rbind(data.frame(NAME = "TYPE", Seurat_Global_Clusters = "group", Seurat_Datasets_Global_Clusters = "group"), groupsMergedDataframeString)
    colnames(groupsMergedDataframeStringTYPE) <- c("NAME", paste("Seurat_Global_Clusters_Resolution", Resolution, sep = "", collapse = ""), paste("Seurat_Datasets_Global_Clusters_Resolution", Resolution, sep = "", collapse = ""))
    
    OutfileClustersMerged<-paste0(Tempdir,"/","CRESCENT_CLOUD/frontend_groups/","groups.tsv")
    write.table(data.frame(groupsMergedDataframeStringTYPE), file = OutfileClustersMerged, row.names = F, col.names = T, sep="\t", quote = F, append = T)
  }
  
  OutfileNumbCellsPerClusterPerDatasetType<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "NumbCellsPerClusterPerDatasetType.tsv")
  Headers<-paste("Cluster", "Number_of_cells", sep="\t", collapse = "")
  write.table(Headers, file = OutfileNumbCellsPerClusterPerDatasetType, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(table(GlobalClustersByDatasetType),file = OutfileNumbCellsPerClusterPerDatasetType, row.names = F, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$WriteOutDatasetTypeClusterAssignments  <- Sys.time()
  
  ####################################
  ### Loops each dataset type
  ####################################
  # switch identity to dataset type global cluster
  Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset_type
  
  NumberOfDatasetTypes <- 0
  for (dataset_type in DatasetTypes) {
    NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
    print(NumberOfDatasetTypes)
    
    if (exists(x = "seurat.object.each_dataset_type") == T) {
      rm(seurat.object.each_dataset_type)
    }
    seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
    
    for (dim_red_method in names(DimensionReductionMethods)) {
      
      ####################################
      ### Colour dimension reduction plots for each dataset type by global clusters
      ####################################
      
      writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for ", dataset_type, " using global cell clusters ***\n"))
      
      StopWatchStart$DimRedOPlotDatasetTypeColourByCellCluster$dataset_type$dim_red_method  <- Sys.time()
      
      plots <- DimPlot(seurat.object.each_dataset_type, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 10)
      plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 3))))
      IntegratedDimRedPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", dataset_type,"_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf")
      pdf(file=IntegratedDimRedPlotPdf, width = 7, height = 8)
      print(CombinePlots(plots))
      dev.off()
      
      StopWatchEnd$DimRedOPlotDatasetTypeColourByCellCluster$dataset_type$dim_red_method  <- Sys.time()
      
      ####################################
      ### Colour dimension reduction plots for each dataset type by selected genes (option -g)
      ####################################
      
      if (regexpr("^NA$", InfileSelectedGenes, ignore.case = T)[1] == 1) {
        print("No selected genes for each dataset type dimension reduction plots")
        
      }else if (3 %in% RequestedApplySelectedGenes == T) {
        writeLines(paste0("\n*** Colour ", DimensionReductionMethods[[dim_red_method]][["name"]], " plot for each dataset type by selected genes (option -g) ***\n"))
        
        StopWatchStart$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
        
        ### Making a new Seurat object `seurat.object.each_dataset_type.sa` with one single slot `@assays$RNA@data` to avoid `FeaturePlot` calling:
        ### Warning: Found the following features in more than one assay, excluding the default. We will not include these in the final dataframe:...
        seurat.object.each_dataset_type.sa <- CreateSeuratObject(seurat.object.each_dataset_type@assays$RNA@data)
        seurat.object.each_dataset_type.sa@reductions$umap <- seurat.object.each_dataset_type@reductions$umap
        seurat.object.each_dataset_type.sa@reductions$tsne <- seurat.object.each_dataset_type@reductions$tsne
        
        OutDirThisOption <- paste0(Tempdir, "/SELECTED_GENE_DIMENSION_REDUCTION_PLOTS/DATASET_TYPES/", DimensionReductionMethods[[dim_red_method]][["name"]])
        dir.create(path = OutDirThisOption, recursive = T)
        sapply(ListOfGenesForDimRedPlots, FUN=function(eachGene) {
          OutfilePathAndName <- paste0(OutDirThisOption, "/", eachGene, "_", dataset_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], ".pdf")
          pdf(file=OutfilePathAndName, width=DefaultParameters$BaseSizeMultiplePlotPdfWidth, height=DefaultParameters$BaseSizeMultiplePlotPdfHeight)
          print(FeaturePlot(object = seurat.object.each_dataset_type.sa, features = eachGene, cols = c("lightgrey", "blue"),
                            reduction = dim_red_method, order = T, slot = "data", pt.size = 0.3, min.cutoff = "q0.1", max.cutoff = "q90"))
          dev.off()
        })
        
        StopWatchEnd$DimRedPlotsColuredByGenes$dim_red_method  <- Sys.time()
        
        rm(seurat.object.each_dataset_type.sa)
      }
      
      ####################################
      ### Colour dimension reduction plots for each dataset type by -infile_metadata using integrated data
      ####################################
      
      if (regexpr("^NA$", InfileMetadata, ignore.case = T)[1] == 1) {
        print("No metadata will be used for dimension reduction plots")
      }else{
        writeLines("\n*** Colour dimension reduction plots for each dataset type by -infile_metadata using integrated data ***\n")
        
        # Note: needs data.frame(CellPropertiesFromMetadata) preloaded
        
        # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same dataset_type)
        # so that barcodes from --infile_metadata and --input can match each other
        rownames(CellPropertiesFromMetadata)<-gsub(x =rownames(CellPropertiesFromMetadata), pattern = "-[0-9]+$", perl = T, replacement = "")
        seurat.object.each_dataset_type <- AddMetaData(object = seurat.object.each_dataset_type, metadata = CellPropertiesFromMetadata)
        
        StopWatchStart$DimRedPlotsByMetadata$dataset_type$dim_red_method  <- Sys.time()
        
        # Generating outfile
        # Note DimPlot() takes the entire current device (pdf)
        # even if using layout(matrix(...)) or  par(mfrow=())
        # Thus each property plot is written to a separate page of a single *pdf outfile
        for (property in colnames(CellPropertiesFromMetadata)) {
          
          if (sum(is.na(seurat.object.each_dataset_type[[property]]==T)) == ncol(seurat.object.each_dataset_type)) {
            
            print(paste0("No extra barcode-attributes '", property, "' found for dataset_type ", dataset_type))
            
          }else{
            
            IntegratedDimRedPlotPdf <- paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourBy_", property, ".pdf")
            pdf(file=IntegratedDimRedPlotPdf, width = DefaultParameters$BaseSizeSinglePlotPdf, height = (DefaultParameters$BaseSizeSinglePlotPdf + (0.12 * length(unique(CellPropertiesFromMetadata[,property])))))
            
            if ( (sum(CellPropertiesFromMetadata[,property] %in% 0:1 == T)) == (nrow(CellPropertiesFromMetadata)) ) { ## is binary
              CellsToHighlight <- rownames(CellPropertiesFromMetadata)[CellPropertiesFromMetadata[,property]==1]
              plots <- DimPlot(seurat.object.each_dataset_type, group.by = property, combine = F, reduction = dim_red_method, label = F, cells.highlight = CellsToHighlight)
              plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
              print(CombinePlots(plots))
              
            }else{
              plots <- DimPlot(seurat.object.each_dataset_type, group.by = property, combine = F, reduction = dim_red_method, label = T)
              plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top", legend.text.align = 0) + labs(title = property) + guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))))
              print(CombinePlots(plots))
              
            }
            
            dev.off()
          }
        }
        StopWatchEnd$DimRedPlotsByMetadata$dataset_type$dim_red_method  <- Sys.time()
      }
    }
  }
  
  ####################################
  ### Load each dataset type global cluster assignments
  ####################################
  writeLines("\n*** Load each dataset type global cluster assignments ***\n")
  
  StopWatchStart$LoadOutDatasetTypeClusterAssignments  <- Sys.time()
  
  EachDatasetTypeGlobalCellClusters <- data.frame(read.table(OutfileEachDatasetTypeGlobalClusters, header = T, row.names = 1, sep = "\t"))
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetTypeGlobalCellClusters, col.name = "EachDatasetTypeGlobalCellClusters")
  
  # switch the identity class of all cells to reflect each dataset clusters
  Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetTypeGlobalCellClusters
  
  StopWatchEnd$LoadOutDatasetTypeClusterAssignments  <- Sys.time()
  
  ####################################
  ### Get average gene expression for each dataset type using global clusters
  ####################################
  writeLines("\n*** Get average gene expression for each dataset type using global clusters ***\n")
  
  StopWatchStart$AverageGeneExpressionEachDatasetTypeGlobalClusters  <- Sys.time()
  
  cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)
  
  OutfileClusterAveragesRNA<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerDatasetType_AverageGeneExpression_RNA.tsv")
  OutfileClusterAveragesSCT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerDatasetType_AverageGeneExpression_SCT.tsv")
  OutfileClusterAveragesINT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "PerDatasetType_AverageGeneExpression_INT.tsv")
  #
  Headers<-paste("AVERAGE_GENE_EXPRESSION", paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
  #
  write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  #
  write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  #
  write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
  StopWatchEnd$AverageGeneExpressionEachDatasetTypeGlobalClusters  <- Sys.time()
  
  
  ####################################
  ### Finding differentially expressed genes (4): using global cell clusers, for each dataset type, compares each cell cluster vs. the rest of cells
  ####################################
  if (4 %in% RequestedDiffGeneExprComparisons == T) {
    
    writeLines("\n*** Finding differentially expressed genes (4): using global cell clusers, for each dataset type, compares each cell cluster vs. the rest of cells ***\n")
    
    # switch identity to dataset type global cluster
    Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$dataset_type
    
    NumberOfDatasetTypes <- 0
    for (dataset_type in DatasetTypes) {
      NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
      print(NumberOfDatasetTypes)
      
      if (exists(x = "seurat.object.each_dataset_type") == T) {
        rm(seurat.object.each_dataset_type)
      }
      seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
      
      StopWatchStart$FindDiffMarkersEachDatasetTypeGlobalClustersVsRestOfCells$dataset_type  <- Sys.time()
      
      NumberOfClusters <- length(unique(seurat.object.each_dataset_type$seurat_clusters))
      
      print(paste0("Number of clusters = ", NumberOfClusters))
      
      Idents(object = seurat.object.each_dataset_type) <- seurat.object.each_dataset_type@meta.data$EachDatasetTypeGlobalCellClusters
      
      FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset_type@meta.data))
      seurat.object.each_dataset_type.markers <- FindAllMarkers(object = seurat.object.each_dataset_type, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
      SimplifiedDiffExprGenes.df <- seurat.object.each_dataset_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
      OutfileDiffGeneExpression <- paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_" , dataset_type, "_MarkersPerCluster.tsv")
      write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
      
      StopWatchEnd$FindDiffMarkersEachDatasetTypeGlobalClustersVsRestOfCells$dataset_type  <- Sys.time()
    }
  }
  
  ####################################
  ### Finding differentially expressed genes (5): using global cell clusers, for each dataset type, compares each cell cluster vs. the same cluster from other dataset types
  ####################################
  
  if (5 %in% RequestedDiffGeneExprComparisons == T) {
    
    writeLines("\n*** Finding differentially expressed genes (5): using global cell clusers, compares each cell cluster from each dataset type vs. the same cluster from other dataset types ***\n")
    
    StopWatchStart$FindDiffMarkersEachDatasetTypeGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()
    
    NumberOfClusters <- length(unique(seurat.object.integrated$seurat_clusters))
    
    print(paste0("Number of clusters = ", NumberOfClusters))
    
    # switch the identity class of all cells to reflect each dataset clusters
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetTypeGlobalCellClusters
    
    FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
    
    OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_GlobalClustering_", "MarkersPerDatasetTypeEquivalentClusters.tsv")
    HeadersOrder <- paste("cluster1", "cluster2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
    write.table(HeadersOrder,file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="", quote = F)
    
    for (cluster in unique(seurat.object.integrated$seurat_clusters)) {
      for (dataset_type1 in unique(InputsTable[,"DatasetType"])) {
        for (dataset_type2 in unique(InputsTable[,"DatasetType"])) {
          Cluster1 <- paste(dataset_type1, cluster, sep = "_c")
          Cluster2 <- paste(dataset_type2, cluster, sep = "_c")
          if (Cluster1 == Cluster2) {
            ### Skip
          }else if ((sum(seurat.object.integrated$EachDatasetTypeGlobalCellClusters == Cluster1) >= 3) & ((sum(seurat.object.integrated$EachDatasetTypeGlobalCellClusters == Cluster2) >= 3))) {
            print (paste0(Cluster1, " vs. ", Cluster2))
            seurat.object.integrated.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.integrated, only.pos = F, ident.1 = Cluster1, ident.2 = Cluster2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
            seurat.object.integrated.each_equivalent_cluster.markers$cluster1 <- Cluster1
            seurat.object.integrated.each_equivalent_cluster.markers$cluster2 <- Cluster2
            seurat.object.integrated.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.integrated.each_equivalent_cluster.markers)
            write.table(seurat.object.integrated.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="\t", quote = F, append = T)
            rm(seurat.object.integrated.each_equivalent_cluster.markers)
          }else{
            print(paste0("Skip cluster ", Cluster1, " vs. ", Cluster2, " because there were not >= 3 cells in at least one of them"))
          }
        }
      }
    }
    
    StopWatchEnd$FindDiffMarkersEachDatasetTypeGlobalClustersVsSameClusterInOtherDatasets  <- Sys.time()
  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO PROCESS EACH DATASET TYPE USING RE-CLUSTERED CELLS
################################################################################################################################################
################################################################################################################################################

####################################
### FOR EACH DATASET TYPE uses integrated counts and RE-CLUSTER cells. Then colour dimension reduction plots and get DGE (if requested)
####################################

if (NumberOfDatasetsTypes >= 1) {
  writeLines("\n*** FOR EACH DATASET TYPE uses integrated counts and RE-CLUSTER cells. Then colour dimension reduction plots and get DGE (if requested) ***\n")
  
  ####################################
  ### Prepares each dataset type object
  ####################################
  
  StopWatchStart$PrepareEachDatasetTypeForRecusterAndDimRedPlots <- Sys.time()
  
  seurat.object.integrated$dataset_type<-seurat.object.integrated$dataset
  for (dataset in rownames(InputsTable)) {
    seurat.object.integrated$dataset_type<-mapply(gsub, pattern = dataset, replacement = list_DatasetToType[[dataset]], seurat.object.integrated$dataset_type)
  }
  DatasetTypes <- unique(seurat.object.integrated$dataset_type)
  
  StopWatchEnd$PrepareEachDatasetTypeForRecusterAndDimRedPlots <- Sys.time()
  
  ####################################
  ### Reclusters each dataset type cells and gets DGE
  ####################################
  if (3 %in% RequestedClusteringInputs == T) {
    
    ####################################
    ### Prepares outfile headers
    ####################################
    OutfileEachDatasetTypeReClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", "CellClusters.tsv")
    Headers<-paste("Cell_barcode", paste0("seurat_cluster_r", Resolution), sep="\t")
    write.table(Headers,file = OutfileEachDatasetTypeReClusters, row.names = F, col.names = F, sep="\t", quote = F)
    
    ####################################
    ### Loops each dataset type
    ####################################
    NumberOfDatasetTypes <- 0
    for (dataset_type in DatasetTypes) {
      NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
      print(NumberOfDatasetTypes)
      
      ####################################
      ### Subsets seurat object per dataset type
      ####################################
      if (exists(x = "seurat.object.each_dataset_type") == T) {
        rm(seurat.object.each_dataset_type)
      }
      
      Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset_type
      
      seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
      print(seurat.object.each_dataset_type)
      
      ####################################
      ### Re-clusters each dataset type
      ####################################
      
      StopWatchStart$ReclusterEachDatasetTypeCellsFromInteg$dataset_type  <- Sys.time()
      
      options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
      seurat.object.each_dataset_type <- FindNeighbors(object = seurat.object.each_dataset_type, dims = PcaDimsUse)
      seurat.object.each_dataset_type <- FindClusters(object = seurat.object.each_dataset_type, resolution = Resolution)
      ClustersThisDatasetType <- unlist(x = strsplit(x = paste(dataset_type, seurat.object.each_dataset_type@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
      
      StopWatchEnd$ReclusterEachDatasetTypeCellsFromInteg$dataset_type  <- Sys.time()
      
      ####################################
      ### Write out each dataset type re-clustered cell clusters
      ####################################
      writeLines("\n*** Write out each dataset type re-clustered cell clusters ***\n")
      
      StopWatchStart$WriteReClustersEachDatasetTypeTables$dataset_type <- Sys.time()
      
      write.table(paste(colnames(seurat.object.each_dataset_type), ClustersThisDatasetType, sep = "\t", collapse = "\n"),
                  file = OutfileEachDatasetTypeReClusters, row.names = F, col.names = F, quote = F, append = T)
      
      CellNames<-rownames(seurat.object.each_dataset_type@meta.data)
      ClusterIdent <-seurat.object.each_dataset_type@meta.data$seurat_clusters
      NumberOfClusters<-length(unique(ClusterIdent))
      OutfileNumbClusters<-paste0(Tempdir, "/CELL_CLUSTER_IDENTITIES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", dataset_type, "_NumbCellClusters", ".tsv")
      write(x=NumberOfClusters,file = OutfileNumbClusters)
      
      StopWatchEnd$WriteReClustersEachDatasetTypeTables$dataset_type <- Sys.time()
      
      ####################################
      ### Colour dimension reduction plots by each dataset type re-clustered cells
      ####################################
      
      writeLines("\n*** Colour dimension reduction plots by each dataset type re-clustered cells ***\n")
      
      for (dim_red_method in names(DimensionReductionMethods)) {
        
        StopWatchStart$DimRedPlotsReclusteredByDatasetTypeColuredByRecluster$dim_red_method$dataset_type  <- Sys.time()
        
        plots <- DimPlot(seurat.object.each_dataset_type, group.by = c("seurat_clusters"), combine = F, reduction = dim_red_method, label = T, label.size = 5)
        plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3))))
        IntegratedUMAPPlotPdf<-paste0(Tempdir, "/DIMENSION_REDUCTION_PLOTS/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", dataset_type, "_", DimensionReductionMethods[[dim_red_method]][["name"]], "Plot_ColourByCellClusters.pdf")
        pdf(file=IntegratedUMAPPlotPdf, width = 7, height = 8)
        print(CombinePlots(plots))
        dev.off()
        
        StopWatchEnd$DimRedPlotsReclusteredByDatasetTypeColuredByRecluster$dim_red_method$dataset_type <- Sys.time()
      }
      ####################################
      ### Finding differentially expressed genes (7): using re-clustered cells, for each dataset type, compares each cell cluster vs. the rest of cells
      ####################################
      
      if (7 %in% RequestedDiffGeneExprComparisons == T) {
        
        writeLines("\n*** Finding differentially expressed genes (7): using re-clustered cells, for each dataset type, compares each cell cluster vs. the rest of cells ***\n")
        
        StopWatchStart$FindDiffMarkersEachDatasetTypeReclusteredVsRestOfCells$dataset_type  <- Sys.time()
        
        print(paste0("Number of clusters = ", NumberOfClusters))
        
        FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset_type@meta.data))
        seurat.object.each_dataset_type.markers <- FindAllMarkers(object = seurat.object.each_dataset_type, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
        SimplifiedDiffExprGenes.df <- seurat.object.each_dataset_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
        OutfileDiffGeneExpression <- paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", dataset_type, "_MarkersPerCluster.tsv")
        write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
        
        StopWatchEnd$FindDiffMarkersEachDatasetTypeReclusteredVsRestOfCells$dataset_type  <- Sys.time()
        
      }
    }
    
    ####################################
    ### Load each dataset type re-cluster assignments
    ####################################
    writeLines("\n*** Load each dataset type re-cluster assignments ***\n")
    
    EachDatasetTypeReclusteredCellClusters <- data.frame(read.table(OutfileEachDatasetTypeReClusters, header = T, row.names = 1, sep = "\t"))
    seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachDatasetTypeReclusteredCellClusters, col.name = "EachDatasetTypeReclusteredCellClusters")
    
    # switch the identity class of all cells to reflect each dataset_type clusters
    Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachDatasetTypeReclusteredCellClusters
    
    ####################################
    ### Get average expression for each dataset type re-clustered clusters
    ####################################
    writeLines("\n*** Get average expression for each dataset type re-clustered clusters ***\n")
    
    StopWatchStart$AverageGeneExpressionEachDatasetTypeReCluster <- Sys.time()
    
    cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)
    
    OutfileClusterAveragesRNA<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", "PerDatasetType_AverageGeneExpression_RNA.tsv")
    OutfileClusterAveragesSCT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", "PerDatasetType_AverageGeneExpression_SCT.tsv")
    OutfileClusterAveragesINT<-paste0(Tempdir, "/AVERAGE_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachDatasetTypeReclustered_", "PerDatasetType_AverageGeneExpression_INT.tsv")
    #
    Headers<-paste("AVERAGE_GENE_EXPRESSION", paste(names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
    #
    write.table(Headers,file = OutfileClusterAveragesRNA, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAveragesRNA, row.names = T, col.names = F, sep="\t", quote = F, append = T)
    #
    write.table(Headers,file = OutfileClusterAveragesSCT, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(cluster.averages$SCT),file = OutfileClusterAveragesSCT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
    #
    write.table(Headers,file = OutfileClusterAveragesINT, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(data.frame(cluster.averages$integrated),file = OutfileClusterAveragesINT, row.names = T, col.names = F, sep="\t", quote = F, append = T)
    
    StopWatchEnd$AverageGeneExpressionEachDatasetTypeReCluster <- Sys.time()
    
  }
}else{
  writeLines("\n*** No dataset type analyzes will be conducted because only '1' dataset type was found in -inputs_list ***\n")
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO OBTAIN DIFFERENTIAL GENE EXPRESSION USING METADATA
################################################################################################################################################
################################################################################################################################################

####################################
### Finding differentially expressed genes: using meatadata classes
####################################
if ((8 %in% RequestedDiffGeneExprComparisons == T) | 
    (9 %in% RequestedDiffGeneExprComparisons == T) | 
    (10 %in% RequestedDiffGeneExprComparisons == T) |
    (11 %in% RequestedDiffGeneExprComparisons == T) | 
    (12 %in% RequestedDiffGeneExprComparisons == T)
    )  {
  
  ####################################
  ### Define metadata for DGE
  ####################################
  
  writeLines(paste0("\n*** Define metadata for DGE ***\n"))
  
  StopWatchStart$LoadMetadataForDGE  <- Sys.time()
  
  MetadataColNamesForDge.list  = unlist(strsplit(MetadataColNamesForDge, ","))
  
  # Note: needs data.frame(CellPropertiesFromMetadata) preloaded
  
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same dataset)
  # so that barcodes from --infile_metadata and --input can match each other
  rownames(CellPropertiesFromMetadata)<-gsub(x =rownames(CellPropertiesFromMetadata), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = CellPropertiesFromMetadata)
  
  StopWatchEnd$LoadMetadataForDGE  <- Sys.time()
  
  ####################################
  ### Finding differentially expressed genes (8): using metadata, compares each cell class vs. the rest of cells
  ####################################
  if (8 %in% RequestedDiffGeneExprComparisons == T) {
    
    writeLines("\n*** Finding differentially expressed genes (8): using metadata, compares each cell class vs. the rest of cells ***\n")
    
    StopWatchStart$FindDiffMarkersEachMetadataClassVsRestOfCells  <- Sys.time()
    
    for (property in MetadataColNamesForDge.list) {
      
      NumberOfClassesInThisProperty <- length(unique(seurat.object.integrated@meta.data[[property]]))
      
      print(paste0("Number of classes in '", property, "' = ", NumberOfClassesInThisProperty))
      
      if (NumberOfClassesInThisProperty > 1) {
        
        ####################################
        ### Subsets seurat object per property
        ####################################
        if (exists(x = "seurat.object.each_property") == T) {
          rm(seurat.object.each_property)
        }
        # switch the identity class of all cells to reflect property identities
        Idents(object = seurat.object.integrated) <- seurat.object.integrated[[property]]
        seurat.object.each_property <- seurat.object.integrated
        
        ####################################
        ### Finding differentially expressed genes
        ####################################
        FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_property@meta.data))
        seurat.object.each_property.markers <- FindAllMarkers(object = seurat.object.each_property, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
        SimplifiedDiffExprGenes.df <- seurat.object.each_property.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
        colnames(SimplifiedDiffExprGenes.df) <- c("class","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")
        OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachMetadataClass_", property, "_MarkersPerClass.tsv")
        write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
      }
    }
    
    StopWatchEnd$FindDiffMarkersEachMetadataClassVsRestOfCells  <- Sys.time()
    
  }
  
  ####################################
  ### Finding differentially expressed genes (9): using metadata annotations, for each dataset, compares each cell class vs. the rest of cells
  ####################################
  if (9 %in% RequestedDiffGeneExprComparisons == T) {
    
    writeLines("\n*** Finding differentially expressed genes (9): using metadata annotations, for each dataset, compares each cell class vs. the rest of cells ***\n")
    
    ####################################
    ### Loops each dataset
    ####################################
    NumberOfDatasets <- 0
    for (dataset in rownames(InputsTable)) {
      NumberOfDatasets <- NumberOfDatasets + 1
      print(NumberOfDatasets)
      
      if (exists(x = "seurat.object.each_dataset") == T) {
        rm(seurat.object.each_dataset)
      }
      
      StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetVsRestOfCells$dataset <- Sys.time()
      
      Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset
      seurat.object.each_dataset <- subset(x = seurat.object.integrated, idents = dataset)
      Idents(object = seurat.object.integrated) <- seurat.object.integrated$seurat_clusters
      Idents(object = seurat.object.each_dataset) <- seurat.object.each_dataset@meta.data[[property]]
      print(seurat.object.each_dataset)
      
      for (property in MetadataColNamesForDge.list) {
        NumberOfClassesInThisProperty <- length(unique(seurat.object.each_dataset@meta.data[[property]]))
        
        print(paste0("Number of classes in '", property, " in dataset ", dataset, " = ", NumberOfClassesInThisProperty))
        
        if (NumberOfClassesInThisProperty > 1) {
          print(paste0("**** ", property, " ****\n"))
          print(table(seurat.object.each_dataset@active.ident))
          
          FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset@meta.data))
          seurat.object.each_dataset.markers <- FindAllMarkers(object = seurat.object.each_dataset, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
          SimplifiedDiffExprGenes.df <- seurat.object.each_dataset.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
          colnames(SimplifiedDiffExprGenes.df) <- c("class","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")
          OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset, "_EachMetadataClass_", property, "_MarkersPerClass.tsv")
          write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
        }
      }
      
      StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetVsRestOfCells$dataset <- Sys.time()
      
    }
  }
  
  ####################################
  ### Finding differentially expressed genes (10): using metadata annotations, for each dataset, compares each cell class specified vs. the same class from other datasets
  ####################################
  if (10 %in% RequestedDiffGeneExprComparisons == T) {
    
    writeLines("\n*** Finding differentially expressed genes (10): using metadata annotations, for each dataset, compares each cell class specified vs. the same class from other datasets ***\n")
    
    StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetVsSameClassInOtherDatasets  <- Sys.time()
    
    for (property in MetadataColNamesForDge.list) {
      
      NumberOfClassesInThisProperty <- length(unique(seurat.object.integrated@meta.data[[property]]))
      
      print(paste0("Number of classes in '", property, "' = ", NumberOfClassesInThisProperty))
      print(paste0("Number of datasets = ", nrow(InputsTable)))
      
      if (NumberOfClassesInThisProperty > 1) {
        
        ####################################
        ### Subsets seurat object per property
        ####################################
        if (exists(x = "seurat.object.each_property") == T) {
          rm(seurat.object.each_property)
        }
        # switch the identity class of all cells to reflect dataset_property identities
        seurat.object.each_property <- seurat.object.integrated
        DatasetAndProperty <- unlist(x = strsplit(x = paste(seurat.object.each_property@meta.data$dataset, seurat.object.each_property@meta.data[[property]], sep = "_", collapse = "\n"), split = "\n"))
        seurat.object.each_property <- AddMetaData(object = seurat.object.each_property, metadata = DatasetAndProperty, col.name = "DatasetAndProperty")
        
        Idents(object = seurat.object.each_property) <- DatasetAndProperty
        
        FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_property@meta.data))
        
        OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachMetadataClass_", property, "_MarkersPerDatasetEquivalentClasses.tsv")
        HeadersOrder <- paste("class1", "class2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
        write.table(HeadersOrder,file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="", quote = F)
        
        for (class in unique(seurat.object.each_property@meta.data[[property]])) {
          for (dataset_1 in rownames(InputsTable)) {
            for (dataset_2 in rownames(InputsTable)) {
              Class_Dataset1 <- paste(dataset_1, class, sep = "_")
              Class_Dataset2 <- paste(dataset_2, class, sep = "_")
              N_Class_Dataset1 <- sum(seurat.object.each_property$DatasetAndProperty == Class_Dataset1)
              N_Class_Dataset2 <- sum(seurat.object.each_property$DatasetAndProperty == Class_Dataset2)
              
              if (Class_Dataset1 == Class_Dataset2) {
                ### Skip
              }else if (N_Class_Dataset1 >= 3 & N_Class_Dataset2 >= 3) {
                print (paste0(Class_Dataset1, " vs. ", Class_Dataset2))
                seurat.object.each_property.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.each_property, only.pos = F, ident.1 = Class_Dataset1, ident.2 = Class_Dataset2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
                seurat.object.each_property.each_equivalent_cluster.markers$class1 <- Class_Dataset1
                seurat.object.each_property.each_equivalent_cluster.markers$class2 <- Class_Dataset2
                seurat.object.each_property.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.each_property.each_equivalent_cluster.markers)
                write.table(seurat.object.each_property.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="\t", quote = F, append = T)
                rm(seurat.object.each_property.each_equivalent_cluster.markers)
              }else{
                print(paste0("Skip class ", Class_Dataset1, " vs. ", Class_Dataset2, " because there were not >= 3 cells in at least one of them"))
              }
            }
          }
        }
      }
    }
    StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetVsSameClassInOtherDatasets  <- Sys.time()
  }
  
  ####################################
  ### Finding differentially expressed genes (11): using metadata annotations, for each dataset_type, compares each cell class vs. the rest of cells
  ####################################
  if (11 %in% RequestedDiffGeneExprComparisons == T) {
    
    writeLines("\n*** Finding differentially expressed genes (10): using metadata annotations, for each dataset_type, compares each cell class vs. the rest of cells ***\n")
    
    ####################################
    ### Loops each dataset_type
    ####################################
    NumberOfDatasetTypes <- 0
    for (dataset_type in DatasetTypes) {
      NumberOfDatasetTypes <- NumberOfDatasetTypes + 1
      print(NumberOfDatasetTypes)
      
      if (exists(x = "seurat.object.each_dataset_type") == T) {
        rm(seurat.object.each_dataset_type)
      }
      
      StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetTypeVsRestOfCells$dataset_type <- Sys.time()
      
      Idents(object = seurat.object.integrated) <- seurat.object.integrated$dataset_type
      seurat.object.each_dataset_type <- subset(x = seurat.object.integrated, idents = dataset_type)
      Idents(object = seurat.object.integrated) <- seurat.object.integrated$seurat_clusters
      Idents(object = seurat.object.each_dataset_type) <- seurat.object.each_dataset_type@meta.data[[property]]
      print(seurat.object.each_dataset_type)
      
      for (property in MetadataColNamesForDge.list) {
        NumberOfClassesInThisProperty <- length(unique(seurat.object.each_dataset_type@meta.data[[property]]))
        
        print(paste0("Number of classes in '", property, " in dataset_type ", dataset_type, " = ", NumberOfClassesInThisProperty))
        
        if (NumberOfClassesInThisProperty > 1) {
          print(paste0("**** ", property, " ****\n"))
          print(table(seurat.object.each_dataset_type@active.ident))
          
          FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_dataset_type@meta.data))
          seurat.object.each_dataset_type.markers <- FindAllMarkers(object = seurat.object.each_dataset_type, only.pos = F, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
          SimplifiedDiffExprGenes.df <- seurat.object.each_dataset_type.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]
          colnames(SimplifiedDiffExprGenes.df) <- c("class","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")
          OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_", dataset_type, "_EachMetadataClass_", property, "_MarkersPerClass.tsv")
          write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = OutfileDiffGeneExpression, row.names = F, sep="\t", quote = F)
        }
      }
      
      StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetTypeVsRestOfCells$dataset_type <- Sys.time()
    }
  }
  
  ####################################
  ### Finding differentially expressed genes (12): using metadata annotations, for each dataset type, compares each cell class vs. the same class from other dataset types
  ####################################
  if (12 %in% RequestedDiffGeneExprComparisons == T) {
    
    writeLines("\n*** Finding differentially expressed genes (12): using metadata annotations, for each dataset type, compares each cell class vs. the same class from other dataset types ***\n")
    
    StopWatchStart$FindDiffMarkersEachMetadataClassEachDatasetTypeVsSameClassInOtherDatasetTypes  <- Sys.time()
    
    if (NumberOfDatasetsTypes > 1) {
      
      for (property in MetadataColNamesForDge.list) {
        
        NumberOfClassesInThisProperty <- length(unique(seurat.object.integrated@meta.data[[property]]))
        
        print(paste0("Number of classes in '", property, "' = ", NumberOfClassesInThisProperty))
        print(paste0("Number of dataset types = ", NumberOfDatasetsTypes))
        
        if (NumberOfClassesInThisProperty > 1) {
          
          ####################################
          ### Subsets seurat object per property
          ####################################
          if (exists(x = "seurat.object.each_property") == T) {
            rm(seurat.object.each_property)
          }
          # switch the identity class of all cells to reflect dataset_type_property identities
          seurat.object.each_property <- seurat.object.integrated
          DatasetTypeAndProperty <- unlist(x = strsplit(x = paste(seurat.object.each_property@meta.data$dataset_type, seurat.object.each_property@meta.data[[property]], sep = "_", collapse = "\n"), split = "\n"))
          seurat.object.each_property <- AddMetaData(object = seurat.object.each_property, metadata = DatasetTypeAndProperty, col.name = "DatasetTypeAndProperty")
          
          Idents(object = seurat.object.each_property) <- DatasetTypeAndProperty
          
          FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_property@meta.data))
          
          OutfileDiffGeneExpression<-paste0(Tempdir, "/DIFFERENTIAL_GENE_EXPRESSION_TABLES/", PrefixOutfiles, ".", ProgramOutdir, "_EachMetadataClass_", property, "_MarkersPerDatasetTypeEquivalentClasses.tsv")
          HeadersOrder <- paste("class1", "class2", "gene", "p_val","p_val_adj","avg_logFC","pct.1","pct.2", sep = "\t")
          write.table(HeadersOrder,file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="", quote = F)
          
          for (class in unique(seurat.object.each_property@meta.data[[property]])) {
            for (dataset_type_1 in unique(InputsTable[,"DatasetType"])) {
              for (dataset_type_2 in unique(InputsTable[,"DatasetType"])) {
                Class_DatasetType1 <- paste(dataset_type_1, class, sep = "_")
                Class_DatasetType2 <- paste(dataset_type_2, class, sep = "_")
                N_Class_DatasetType1 <- sum(seurat.object.each_property$DatasetTypeAndProperty == Class_DatasetType1)
                N_Class_DatasetType2 <- sum(seurat.object.each_property$DatasetTypeAndProperty == Class_DatasetType2)
                
                if (Class_DatasetType1 == Class_DatasetType2) {
                  ### Skip
                }else if (N_Class_DatasetType1 >= 3 & N_Class_DatasetType2 >= 3) {
                  print (paste0(Class_DatasetType1, " vs. ", Class_DatasetType2))
                  seurat.object.each_property.each_equivalent_cluster.markers <- data.frame(FindMarkers(object = seurat.object.each_property, only.pos = F, ident.1 = Class_DatasetType1, ident.2 = Class_DatasetType2, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount))
                  seurat.object.each_property.each_equivalent_cluster.markers$class1 <- Class_DatasetType1
                  seurat.object.each_property.each_equivalent_cluster.markers$class2 <- Class_DatasetType2
                  seurat.object.each_property.each_equivalent_cluster.markers$gene     <- rownames(seurat.object.each_property.each_equivalent_cluster.markers)
                  write.table(seurat.object.each_property.each_equivalent_cluster.markers[,unlist(strsplit(HeadersOrder, "\t"))], file = OutfileDiffGeneExpression, row.names = F, col.names = F, sep="\t", quote = F, append = T)
                  rm(seurat.object.each_property.each_equivalent_cluster.markers)
                }else{
                  print(paste0("Skip class ", Class_DatasetType1, " vs. ", Class_DatasetType2, " because there were not >= 3 cells in at least one of them"))
                }
              }
            }
          }
        }
      }
    }else{
      print(paste0("Skip because could find only '", NumberOfDatasetsTypes, "' dataset types"))
    }
    StopWatchEnd$FindDiffMarkersEachMetadataClassEachDatasetTypeVsSameClassInOtherDatasetTypes  <- Sys.time()
  }
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SAVE THE WHOLE RUN R_OBJECT AND LOG FILES
################################################################################################################################################
################################################################################################################################################

####################################
### Saving the full R object
####################################

if (regexpr("^Y$", SaveRObject, ignore.case = T)[1] == 1) {
  
  writeLines("\n*** Saving the full R object ***\n")
  
  StopWatchStart$SaveRDSFull  <- Sys.time()
  
  OutfileRDS<-paste0(Tempdir, "/R_OBJECTS/", PrefixOutfiles, ".", ProgramOutdir, "_integrated_object_full.rds")
  saveRDS(seurat.object.integrated, file = OutfileRDS)
  
  StopWatchEnd$SaveRDSFull  <- Sys.time()
  
}else{
  
  writeLines("\n*** Not saving the full R object ***\n")
  
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUusage<-paste0(Tempdir, "/LOG_FILES/", PrefixOutfiles, ".", ProgramOutdir, "_CPUusage.txt")
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

### using two steps to copy files (`file.copy` and `file.remove`) instead of just `file.rename` to avoid issues with path to Tempdir in cluster systems
if (regexpr("^Y$", RunsCwl, ignore.case = T)[1] == 1) {
  writeLines("\n*** Keeping files at: ***\n")
  writeLines(Tempdir)
} else {
  writeLines("\n*** Moving outfiles into outdir ***\n")
  sapply(FILE_TYPE_OUT_DIRECTORIES, FUN=function(DirName) {
    TempdirWithData <- paste0(Tempdir, "/", DirName)
    if (DirName == "SELECTED_GENE_DIMENSION_REDUCTION_PLOTS") {
      sapply(list.dirs(TempdirWithData, full.names = F, recursive = F), FUN=function(SubDirName) {
        sapply(list.dirs(paste0(TempdirWithData, "/", SubDirName), full.names = F, recursive = F), FUN=function(SubSubDirName) {
          OutdirFinal <- gsub(pattern = Tempdir, replacement =  paste0(Outdir, "/", ProgramOutdir), x = paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName))
          dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
          sapply(list.files(paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName), pattern = ".pdf", full.names = F), FUN=function(EachFileName) {
            file.copy(from=paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
            file.remove(paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName, "/", EachFileName))
          })
        })
      })
    }else if (DirName == "FILTERED_DATA_MATRICES" | DirName == "UNFILTERED_DATA_MATRICES") {
      sapply(list.dirs(TempdirWithData, full.names = F, recursive = F), FUN=function(SubDirName) {
        sapply(list.dirs(paste0(TempdirWithData, "/", SubDirName), full.names = F, recursive = F), FUN=function(SubSubDirName) {
          OutdirFinal <- gsub(pattern = Tempdir, replacement =  paste0(Outdir, "/", ProgramOutdir), x = paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName))
          dir.create(file.path(OutdirFinal), showWarnings = F, recursive = T)
          sapply(list.files(paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName), pattern = ".gz", full.names = F), FUN=function(EachFileName) {
            file.copy(from=paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName, "/", EachFileName), to=paste0(OutdirFinal, "/", EachFileName), overwrite=T)
            file.remove(paste0(TempdirWithData, "/", SubDirName, "/", SubSubDirName, "/", EachFileName))
          })
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
writeLines(paste0("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report"))

quit()
