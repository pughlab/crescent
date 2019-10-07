####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/v3.0/merge_vignette.html
### and https://satijalab.org/seurat/v3.0/sctransform_vignette.html
####################################

####################################
### Required libraries
####################################
suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
### Seurat v3 can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = "satijalab/seurat", ref = "develop")
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
suppressPackageStartupMessages(library(future))       # To run parallel processes
### library(staplr)     only if using option '-s y', note it needs pdftk
####################################

####################################
### Other packages required
####################################
### UMAP
### can be installed using `pip install umap-learn`
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
              help="Path/name to a <tab> delimited file with the list of dataset IDs and infile path/names, like:
                Dataset1  /path_to/dataset1_mtx
                Dataset2  /path_to/dataset2_mtx
                Dataset3  /path_to/dataset3_mtx
                ...etc
                Each infile must be either the path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files;
                or path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns
                Notes:
                The 'MTX' files can be for example the output from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--inputs_type"), default="NA",
              help="Indicates if the inputs are either a 'MTX' directories or 'DGE' files
                Default = 'No default. It's mandatory to specify this parameter'"),
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
  make_option(c("-a", "--opacity"), default="0.1",
              help="If using a --list_genes, this parameter provides a value for the minimal opacity of gene expression. Use a value between 0 and 1
                Default = 'y'"),
  #
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check these two files:
                *JackStraw*pdf, use the number of PC's where the solid curve shows a plateau along the dotted line, and
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers
                Default = '10'"),
  #
  make_option(c("-m", "--percent_mito"), default="0,0.05",
              help="<comma> delimited min,max number of percentage of mitochondrial gene counts in a cell to be included in normalization and clustering analyses
               For example, for whole cell scRNA-seq use '0,0.2', or for Nuc-seq use '0,0.05'
               Default = '0,0.05'"),
  #
  make_option(c("-n", "--n_genes"), default="50,8000",
              help="<comma> delimited min,max number of unique gene counts in a cell to be included in normalization and clustering analyses
                Default = '50,8000'"),
  #
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh
                Default = '0.01'"),
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

InputsList       <- opt$inputs_list
InputsType       <- opt$inputs_type
Resolution       <- as.numeric(opt$resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir           <- opt$outdir
PrefixOutfiles   <- opt$prefix_outfiles
Opacity          <- as.numeric(opt$opacity)
PcaDimsUse       <- c(1:as.numeric(opt$pca_dimensions))
ListPMito        <- opt$percent_mito
ListNGenes       <- opt$n_genes
ThreshReturn     <- as.numeric(opt$return_threshold)
NumbCores        <- opt$number_cores
RunsCwl          <- opt$run_cwl

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

####################################
### Define number of cores and RAM for parallelization
####################################

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (is.numeric(NumbCores) == T) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste("Unexpected format for --number_cores: ", NumbCores, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}

cat("Using ", NumbCoresToUse, "cores")

plan(strategy = "multicore", workers = NumbCoresToUse)

### To avoid a memmory error with getGlobalsAndPackages() while using ScaleData()
### allocate 4Gb of RAM (4000*1024^2), use:
options(future.globals.maxSize = 4000 * 1024^2)

####################################
### Define default parameters
####################################
### Some of these default parameters are provided by Seurat developers,
### others are tailored according to clusters/t-SNE granularity

ListNGenes = unlist(strsplit(ListNGenes, ","))
MinGenes   = as.numeric(ListNGenes[1])
MaxGenes   = as.numeric(ListNGenes[2])
#
ListPMito  = unlist(strsplit(ListPMito,  ","))
MinPMito   = as.numeric(ListPMito[1])
MaxPMito   = as.numeric(ListPMito[2])

DefaultParameters <- list(
  
  ### Parameters for Seurat filters
  MinCells = 3,
  MinGenes = MinGenes,
  MaxGenes = MaxGenes,
  MinPMito = MinPMito,
  MaxPMito = MaxPMito,
  
  ### Parameters for Seurat normalization
  ScaleFactor = 10000,
  
  ### Parameters for Seurat variable gene detection
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
  
  ### Parameters for t-SNE plots
  BaseSizeSingleTnePlot  = 7,
  BaseSizeMultipleWidth  = 3.7,
  BaseSizeMultipleHeight = 3
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

ListMandatory<-list("infiles_list", "inputs_type", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs
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
dir.create(file.path(Outdir, "SEURAT"), recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

InputsList<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), InputsList)
InputsList
InputsTable<-read.table(InputsList, header = F, row.names = 1)

SeuratObjects   <- list()
DatasetIds      <- list()
NumberOfDatasets <- 0

for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  Dataset.SO <-paste(dataset, ".so",    sep = "", collapse = "")
  
  PathToDataset<-InputsTable[dataset,1]
  PathToDataset<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), PathToDataset)
  
  if (regexpr("^MTX$", InputsType, ignore.case = T)[1] == 1) {
    
    ####################################
    ### Loading MTX infiles
    ####################################
    
    StopWatchStart$LoadScRNAseqData$dataset <- Sys.time()
    
    writeLines(paste("\n*** Loading MTX infiles for ", dataset, " from: ", PathToDataset, " ***\n", sep = "", collapse = ""))
    expression_matrix <- Read10X(data.dir = PathToDataset)
    
    StopWatchEnd$LoadScRNAseqData$dataset  <- Sys.time()
    
    ####################################
    ### Create seurat object
    ####################################
    
    StopWatchStart$CreateSeuratObject$dataset  <- Sys.time()
    
    writeLines(paste("\n*** Create seurat object for ", dataset, " ***\n", sep = "", collapse = ""))
    seurat.object.u  <- CreateSeuratObject(counts = expression_matrix, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = paste(PrefixOutfiles, "_", dataset, sep = "", collapse = ""))
    
    StopWatchEnd$CreateSeuratObject$dataset  <- Sys.time()
    
    ####################################
    ### Add dataset label
    ####################################
    writeLines(paste("\n*** Add dataset label for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$AddDatasetLabel$dataset  <- Sys.time()
    
    seurat.object.u[['dataset.label']] <- dataset
    
    StopWatchEnd$AddDatasetLabel$dataset  <- Sys.time()
    
    
    ####################################
    ### Get mitochondrial genes
    ####################################
    writeLines(paste("\n*** Get  mitochondrial genes for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$GetMitoGenes$dataset  <- Sys.time()
    
    mitoRegExpressions<- paste(c("^MT-"), collapse = "|")
    mito.features <- grep(pattern = mitoRegExpressions, ignore.case = T, x = rownames(x = seurat.object.u), value = T)
    
    if (length(mito.features)[[1]] > 0) {
      percent.mito <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.mito']] <- percent.mito
    }else{
      percent.mito <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.mito']] <- percent.mito
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
      percent.ribo <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[ribo.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.ribo']] <- percent.ribo
    }else{
      percent.ribo <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
      seurat.object.u[['percent.ribo']] <- percent.ribo
    }
    
    StopWatchEnd$GetRiboGenes$dataset  <- Sys.time()
    
    ####################################
    ### Filter cells based gene counts and mitochondrial representation
    ####################################
    writeLines(paste("\n*** Filter cells based gene counts and mitochondrial representation for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$FilterCells$dataset  <- Sys.time()
    
    if (length(mito.features)[[1]] > 0) {
      seurat.object.integrated<-subset(x = seurat.object.u, subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes & percent.mito > DefaultParameters$MinPMito & percent.mito < DefaultParameters$MaxPMito)
    }else{
      seurat.object.integrated<-subset(x = seurat.object.u, subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes)
    }
    
    StopWatchEnd$FilterCells$dataset  <- Sys.time()
    
    ### Just reporting the summary of the UNfiltered and filtered objects
    seurat.object.u
    seurat.object.integrated
    
    ####################################
    ### QC EDA violin plots
    ####################################
    writeLines(paste("\n*** QC EDA violin plots for ", dataset, " ***\n", sep = "", collapse = ""))
    
    StopWatchStart$QCviolinplots$dataset  <- Sys.time()
    
    ### Get unfiltered data QC statistics
    nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
    nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
    percent.mito.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$percent.mito, percent.mito = 1)
    percent.ribo.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$percent.ribo, percent.ribo = 1)
    #
    nFeature_RNAStats.u<-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
    nCount_RNAStats.u  <-paste(c( "mean = ",round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
    percent.mito.u     <-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object.u@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
    percent.ribo.u     <-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object.u@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")
    
    ### Get filtered data QC statistics
    nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object.integrated@meta.data$nFeature_RNA, nGenes = 2)
    nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object.integrated@meta.data$nCount_RNA,   nCount_RNA = 2)
    percent.mito.f.df  <-data.frame(Expression_level = seurat.object.integrated@meta.data$percent.mito, percent.mito = 2)
    percent.ribo.f.df  <-data.frame(Expression_level = seurat.object.integrated@meta.data$percent.ribo, percent.ribo = 2)
    #
    nFeature_RNAStats.f<-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object.integrated@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
    nCount_RNAStats.f  <-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object.integrated@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
    percent.mito.f     <-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object.integrated@meta.data[,"percent.mito"]),3)), sep = "", collapse="")
    percent.ribo.f     <-paste(c(" mean = ",round(mean(seurat.object.integrated@meta.data[,"percent.ribo"]),3),"\n", "median = ",round(median(seurat.object.integrated@meta.data[,"percent.ribo"]),3)), sep = "", collapse="")
    
    ### Put QC statistics together
    nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
    nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
    percent.mito.m.df  <-data.frame(rbind(percent.mito.u.df,percent.mito.f.df))
    percent.ribo.m.df  <-data.frame(rbind(percent.ribo.u.df,percent.ribo.f.df))
    LabelUnfiltered    <-paste("Before filters: No. of cells = ", nrow(seurat.object.u@meta.data), sep ="", collapse = "")
    LabelFiltered      <-paste("After filters:  No. of cells = ", nrow(seurat.object.integrated@meta.data), sep ="", collapse = "")
    #LabelFilters       <-paste("Filters: No. of genes per cell = ", ListNGenes, "Fraction mitochondrial protein genes per cell = ", ListPMito), sep ="", collapse = "")
    
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
      labs(x="No. of genes") +
      annotate("text", x = 1 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nFeature_RNA.m.df$Expression_level)*1.1, label = nFeature_RNAStats.f, col = ColoursQCViolinPlots[[2]])
    
    nCount_RNA.plot<-ggplot(data=nCount_RNA.m.df, aes(x = factor(nCount_RNA), y = Expression_level)) +
      geom_violin(aes(fill = factor(nCount_RNA))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x="No. of reads") +
      annotate("text", x = 1 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(nCount_RNA.m.df$Expression_level)*1.1, label = nCount_RNAStats.f, col = ColoursQCViolinPlots[[2]])
    
    percent.mito.plot<-ggplot(data=percent.mito.m.df, aes(x = factor(percent.mito), y = Expression_level)) +
      geom_violin(aes(fill = factor(percent.mito))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x="Mitochondrial genes (fraction)") +
      annotate("text", x = 1 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.f, col = ColoursQCViolinPlots[[2]])
    
    percent.ribo.plot<-ggplot(data=percent.ribo.m.df, aes(x = factor(percent.ribo), y = Expression_level)) +
      geom_violin(aes(fill = factor(percent.ribo))) + geom_jitter(height = 0, width = 0.1) + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = ColourDefinitions["medium_grey"][[1]]), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = ColoursQCViolinPlots) +
      labs(x="Ribosomal protein genes (fraction)") +
      annotate("text", x = 1 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.u, col = ColoursQCViolinPlots[[1]]) +
      annotate("text", x = 2 , y = max(percent.ribo.m.df$Expression_level)*1.1, label = percent.ribo.f, col = ColoursQCViolinPlots[[2]])
    
    bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, percent.mito.plot, percent.ribo.plot, ncol = 4)
    
    ### Create a *pdf file with the violin ggplot's
    
    VlnPlotPdf<-paste(Tempdir, "/" , PrefixOutfiles, "_", dataset, ".SEURAT_QC_VlnPlot.pdf", sep="", collapse = )
    pdf(file=VlnPlotPdf, width = 12, height = 7)
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
                                  percent.mito = seurat.object.u@meta.data$percent.mito,
                                  filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.integrated))
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = ColoursQCViolinPlots[[1]])
    UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = FALSE, replacement = ColoursQCViolinPlots[[2]])
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = 4)
    UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = FALSE, replacement = 16)
    
    FeatureVsFeaturePlotPdf<-paste(Tempdir, "/", PrefixOutfiles, "_", dataset, ".SEURAT_NumbReadsVsNumbGenesAndMito_VlnPlot.pdf", sep="", collapse = "")
    pdf(file=FeatureVsFeaturePlotPdf, width = 10, height = 5)
    par(mfrow=c(1,2))
    ## No. of reads vs. Mitochond. %
    plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$percent.mito, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
         xlab ="No. of reads", ylab = "Mitochond. %")
    legend("topright", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))
    
    ## No. of reads vs. No. of genes
    plot(x = UnfilteredData.df$nCount_RNA, UnfilteredData.df$nGene, col = UnfilteredData.df$DotColour, pch = as.integer(UnfilteredData.df$DotPch),
         xlab ="No. of reads", ylab = "No. of genes")
    legend("topleft", legend = c("No", "Yes"), title = "Filtered cells", col = ColoursQCViolinPlots, pch = c(4,16))
    dev.off()
    
    StopWatchEnd$FeatureVsFeatureplot$dataset  <- Sys.time()
    
    ####################################
    ### Remove the Unfiltered seurat object
    ####################################
    writeLines(paste("\n*** Remove the Unfiltered seurat object for ", dataset, " ***\n", sep = "", collapse = ""))
    
    rm(seurat.object.u)
    rm(UnfilteredData.df)
    
    ####################################
    ### Assign data to Datasets lists
    ####################################
    writeLines(paste("\n*** Assign data to Datasets lists: ", dataset, " ***\n", sep = "", collapse = ""))
    
    SeuratObjects[[as.character(NumberOfDatasets)]]  <- seurat.object.integrated
    DatasetIds[[as.character(NumberOfDatasets)]]     <- dataset
    
  }else{
    stop(paste("Unexpected type of input: ", InputsType, "\n\nFor help type:\n\nnormalizes_and_merges_datasets_with_seurat.R -h\n\n", sep=""))
  }
}

####################################
### Merge Seurat objects
####################################
writeLines("\n*** Merge Seurat objects ***\n")

#### Here, need to automate handling datasets

StopWatchStart$MergeSeuratObjects  <- Sys.time()

FirstSeuratObject <- SeuratObjects[[1]]
FirstSampleId     <- DatasetIds[[1]]
RestOfSeuratObjects <- SeuratObjects[c(2:NumberOfDatasets)]
RestOfSamplesIds    <- unlist(DatasetIds[c(2:NumberOfDatasets)])

seurat.object.merged <- merge(FirstSeuratObject, y = RestOfSeuratObjects, 
                              add.cell.ids = DatasetIds,
                              project = PrefixOutfiles)

StopWatchEnd$MergeSeuratObjects  <- Sys.time()

dataset.label.metadata<-cbind.data.frame(sample=seurat.object.merged@meta.data$dataset.label)
rownames(dataset.label.metadata)<-colnames(seurat.object.merged)

seurat.object.merged.withmetadata <- CreateSeuratObject(seurat.object.merged@assays$RNA@counts, meta.data = dataset.label.metadata)
seurat.object.list <- SplitObject(seurat.object.merged.withmetadata, split.by = "sample")

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

StopWatchStart$Integration  <- Sys.time()

seurat.object.integratedfeatures <- SelectIntegrationFeatures(object.list = seurat.object.list, nfeatures = 3000)
seurat.object.list <- PrepSCTIntegration(object.list = seurat.object.list, anchor.features = seurat.object.integratedfeatures, verbose = FALSE)
seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, normalization.method = "SCT", anchor.features = seurat.object.integratedfeatures, verbose = FALSE)
seurat.object.integrated <- IntegrateData(anchorset = seurat.object.anchors, normalization.method = "SCT", verbose = FALSE)

StopWatchEnd$Integration  <- Sys.time()

####################################
### Reducing dimensions
####################################
writeLines("\n*** Reducing dimensions ***\n")

StopWatchStart$RunPCA  <- Sys.time()

seurat.object.integrated <- RunPCA(seurat.object.integrated, verbose = FALSE)

StopWatchEnd$RunPCA  <- Sys.time()

StopWatchStart$PCAPlots  <- Sys.time()

ForElbowPlot<-ElbowPlot(object = seurat.object.integrated, ndims = 50, reduction = "pca")
MaxYAxis<-as.integer(max(ForElbowPlot$data$stdev)+1)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCElbowPlot_Integrated.pdf", sep=""))
print(ForElbowPlot
      + scale_x_continuous(breaks =  seq(from = 0, to = 50, by=5))
      + geom_vline(xintercept = seq(from = 0, to = 50, by=5), linetype='dotted', col="red")
      + scale_y_continuous(breaks =  seq(from = 0, to = MaxYAxis, by=0.5))
      + geom_hline(yintercept = seq(from = 0, to = MaxYAxis, by=0.5), linetype='dotted', col="red")
)
dev.off()

StopWatchEnd$PCAPlots  <- Sys.time()

NumberOfDimsBasedOnElbowPlot<-17

####################################
### Run UMAP
####################################
writeLines("\n*** Run UMAP ***\n")

StopWatchStart$RunUmap  <- Sys.time()

seurat.object.integrated <- RunUMAP(seurat.object.integrated, dims = 1:NumberOfDimsBasedOnElbowPlot)

StopWatchEnd$RunUmap  <- Sys.time()

####################################
### Write out UMAP coordinates
####################################
writeLines("\n*** Write out UMAP coordinates ***\n")

StopWatchStart$WriteUmapCoords  <- Sys.time()

OutfileUMAPCoordinates<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_UMAPCoordinates.tsv", sep="")
Headers<-paste("Barcode",paste(colnames(seurat.object.integrated@reductions$umap@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileUMAPCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
write.table(seurat.object.integrated@reductions$umap@cell.embeddings, file = OutfileUMAPCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$WriteUmapCoords  <- Sys.time()

####################################
### Run t-SNE
####################################
writeLines("\n*** Run t-SNE ***\n")

StopWatchStart$RunTsne  <- Sys.time()

seurat.object.integrated <- RunTSNE(seurat.object.integrated, dims = 1:NumberOfDimsBasedOnElbowPlot)

StopWatchEnd$RunTsne  <- Sys.time()

####################################
### Write out t-SNE coordinates
####################################
writeLines("\n*** Write out t-SNE coordinates ***\n")

StopWatchStart$WriteTsneCoords <- Sys.time()

OutfileTsneCoordinates<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNECoordinates.tsv", sep="")
Headers<-paste("Barcode", paste(colnames(seurat.object.integrated@reductions$tsne@cell.embeddings), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers, file = OutfileTsneCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
write.table(seurat.object.integrated@reductions$tsne@cell.embeddings, file = OutfileTsneCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

StopWatchEnd$WriteTsneCoords <- Sys.time()

####################################
### Colour UMAP plot by dataset
####################################
writeLines("\n*** Colour UMAP plot by dataset ***\n")

StopWatchStart$UmapColourByDataset  <- Sys.time()

plots <- DimPlot(seurat.object.integrated, group.by = c("sample"), combine = FALSE, reduction = "umap")
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
IntegratedUMAPPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_UMAP_ColourByDataset.pdf", sep="")
pdf(file=IntegratedUMAPPlotPdf, width = 7, height = 8)
print(CombinePlots(plots))
dev.off()

StopWatchEnd$UmapColourByDataset  <- Sys.time()

####################################
### Colour t-SNE plot by dataset
####################################
writeLines("\n*** Colour t-SNE plot by dataset ***\n")

StopWatchStart$TsneColourByDataset  <- Sys.time()

plots <- DimPlot(seurat.object.integrated, group.by = c("sample"), combine = FALSE, reduction = "tsne")
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
IntegratedTSNEPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_ColourByDataset.pdf", sep="")
pdf(file=IntegratedTSNEPlotPdf, width = 7, height = 8)
print(CombinePlots(plots))
dev.off()

StopWatchEnd$TsneColourByDataset  <- Sys.time()

####################################
### Cluster all cells from integrated data
####################################
writeLines("\n*** Cluster all cells from integrated data ***\n")

StopWatchStart$ClusterAllCells  <- Sys.time()

options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.integrated <- FindNeighbors(object = seurat.object.integrated, dims = 1:NumberOfDimsBasedOnElbowPlot)
seurat.object.integrated <- FindClusters(object = seurat.object.integrated, reduction.type = "pca", resolution = Resolution, print.output = 0, save.SNN = T)

StopWatchEnd$ClusterAllCells  <- Sys.time()

StopWatchStart$AllCellClusterTables  <- Sys.time()

CellNames<-rownames(seurat.object.integrated@meta.data)
ClusterIdent <-seurat.object.integrated@meta.data$seurat_clusters

Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = "") ,sep="\t")
clusters_data<-paste(CellNames, ClusterIdent, sep="\t")
#
OutfileClusters<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CellClusters.tsv", sep="")
write.table(Headers,file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(clusters_data),file = OutfileClusters, row.names = F, col.names = F, sep="\t", quote = F, append = T)
#
NumberOfClusters<-length(unique(ClusterIdent))
OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_NumbCellClusters", ".tsv", sep="")
write(x=NumberOfClusters,file = OutfileNumbClusters)

StopWatchEnd$AllCellClusterTables  <- Sys.time()

####################################
### Colour UMAP plot by all cell clusters
####################################
writeLines("\n*** Colour UMAP plot by all cell clusters ***\n")

StopWatchStart$UmapColourByCellCluster  <- Sys.time()

plots <- DimPlot(seurat.object.integrated, group.by = c("seurat_clusters"), combine = FALSE, reduction = "umap", label = T, label.size = 5)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
IntegratedUMAPPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_UMAP_ColourByCellClusters.pdf", sep="")
pdf(file=IntegratedUMAPPlotPdf, width = 7, height = 8)
print(CombinePlots(plots))
dev.off()

StopWatchEnd$UmapColourByCellCluster  <- Sys.time()

####################################
### Colour t-SNE plot by all cell clusters
####################################
writeLines("\n*** Colour t-SNE plot by all cell clusters ***\n")

StopWatchStart$TsneColourByCellCluster  <- Sys.time()

plots <- DimPlot(seurat.object.integrated, group.by = c("seurat_clusters"), combine = FALSE, reduction = "tsne", label = T, label.size = 5)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
IntegratedTSNEPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_ColourByCellClusters.pdf", sep="")
pdf(file=IntegratedTSNEPlotPdf, width = 7, height = 8)
print(CombinePlots(plots))
dev.off()

StopWatchEnd$TsneColourByCellCluster  <- Sys.time()

####################################
### Colour UMAP plots by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo
####################################
writeLines("\n*** Colour UMAP plots by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo ***\n")

StopWatchStart$UmapColourByQC  <- Sys.time()

CellPropertiesToUmap<-c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_UMAPlot_QC.pdf", sep=""), width = 21, height = 5)
p <- FeaturePlot(object = seurat.object.integrated, combine = FALSE, label = T, order = T, features = CellPropertiesToUmap, cols = c("lightgrey", "blue"), reduction = "umap", ncol = 4, pt.size = 1.5)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()
}
print(cowplot::plot_grid(plotlist = p))
dev.off()

StopWatchEnd$UmapColourByQC  <- Sys.time()

####################################
### Colour t-SNE plots by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo
####################################
writeLines("\n*** Colour t-SNE plots by nFeature_RNA, nCount_RNA, percent.mito and percent.ribo ***\n")

StopWatchStart$TsneColourByQC  <- Sys.time()

CellPropertiesToTsne<-c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 21, height = 5)
p <- FeaturePlot(object = seurat.object.integrated, combine = FALSE, label = T, order = T, features = CellPropertiesToTsne, cols = c("lightgrey", "blue"), reduction = "tsne", ncol = 4, pt.size = 1.5)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()
}
print(cowplot::plot_grid(plotlist = p))
dev.off()

StopWatchEnd$TsneColourByQC  <- Sys.time()

####################################
### Get average expression for each cluster, for each gene
####################################
writeLines("\n*** Get average expression for each cluster, for each gene ***\n")

StopWatchStart$AverageGeneExpression  <- Sys.time()

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)
#
OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster_integrated.tsv", sep="")
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
### Create metadata with each sample ID and global clustering cluster assignments
####################################
writeLines("\n*** Create metadata with each sample ID and global clustering cluster assignments ***\n")

seurat.object.integrated$AllSampleClusters <- Idents(object = seurat.object.integrated)

EachSampleGlobalClusteredCellClusters <- unlist(x = strsplit(x = paste(seurat.object.integrated$sample, seurat.object.integrated$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleGlobalClusteredCellClusters, col.name = "EachSampleGlobalClusteredCellClusters")

# switch the identity class of all cells to reflect "clusters"
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleGlobalClusteredCellClusters

####################################
### Get average expression for each sample, for each cluster, for each gene (global clustering)
####################################
writeLines("\n*** Get average expression for each sample, for each cluster, for each gene  (global clustering) ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_AverageGeneExpressionPerSampleGlobalClustering_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_AverageGeneExpressionPerSampleGlobalClustering_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_AverageGeneExpressionPerSampleGlobalClustering_integrated.tsv", sep="")
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
### Cluster each sample cells, colour UMAP and t-SNE plots by cell cluster
### and finding differentially expressed genes from integrated data
####################################
writeLines("\n*** Cluster each sample cells, colour UMAP and t-SNE plots by cell cluster and finding differentially expressed genes from integrated data ***\n")

Idents(object = seurat.object.integrated) <- seurat.object.integrated@meta.data$sample

Headers<-paste("Cell_barcode", paste("seurat_cluster_r", Resolution, sep = "", collapse = ""), sep="\t")
#
OutfileEachSampleClusters<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_EachSampleCellClusters.tsv", sep="")
write.table(Headers,file = OutfileEachSampleClusters, row.names = F, col.names = F, sep="\t", quote = F)

for (dataset in rownames(InputsTable)) {
  NumberOfDatasets <- NumberOfDatasets + 1
  print(NumberOfDatasets)
  
  if (exists(x = "seurat.object.each_sample") == T) {
    rm(seurat.object.each_sample)
  }
  seurat.object.each_sample <- subset(x = seurat.object.integrated, idents = dataset)
  print(seurat.object.each_sample)
  
  StopWatchStart$ClusterEachSampleCellsFromInteg$dataset  <- Sys.time()
  
  options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
  seurat.object.each_sample <- FindNeighbors(object = seurat.object.each_sample, dims = 1:NumberOfDimsBasedOnElbowPlot)
  seurat.object.each_sample <- FindClusters(object = seurat.object.each_sample, reduction.type = "pca", resolution = Resolution, print.output = 0, save.SNN = T)
  ClustersThisSample <- unlist(x = strsplit(x = paste(dataset, seurat.object.each_sample@meta.data$seurat_clusters, sep = "_c", collapse = "\n"), split = "\n"))
  
  StopWatchEnd$ClusterEachSampleCellsFromInteg$dataset  <- Sys.time()
  
  StopWatchStart$ClusterEachSampleCellClusterTables$dataset <- Sys.time()
  
  write.table(paste(colnames(seurat.object.each_sample),ClustersThisSample, sep = "\t", collapse = "\n"),
              file = OutfileEachSampleClusters, row.names = F, col.names = F, quote = F, append = T)

  NumberOfClusters<-length(unique(ClusterIdent))
  OutfileNumbClusters<-paste(Tempdir,"/",PrefixOutfiles, "_", dataset, ".SEURAT_NumbCellClusters", ".tsv", sep="")
  write(x=NumberOfClusters,file = OutfileNumbClusters)
  
  StopWatchEnd$ClusterEachSampleCellClusterTables$dataset <- Sys.time()
  
  ####################################
  ### Colour UMAP and t-SNE plots by each sample cell clusters
  ####################################
  writeLines("\n*** Colour UMAP and t-SNE plots by each sample cell clusters ***\n")
  
  StopWatchStart$UmapAndTsneColourByEachSampleCellCluster$dataset  <- Sys.time()
  
  plots <- DimPlot(seurat.object.each_sample, group.by = c("seurat_clusters"), combine = FALSE, reduction = "umap", label = T, label.size = 5)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
  IntegratedUMAPPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", dataset, ".SEURAT_UMAP_ColourByCellClusters.pdf", sep="")
  pdf(file=IntegratedUMAPPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  plots <- DimPlot(seurat.object.each_sample, group.by = c("seurat_clusters"), combine = FALSE, reduction = "tsne", label = T, label.size = 5)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))))
  IntegratedTSNEPlotPdf<-paste(Tempdir,"/",PrefixOutfiles, "_", dataset, ".SEURAT_TSNEPlot_ColourByCellClusters.pdf", sep="")
  pdf(file=IntegratedTSNEPlotPdf, width = 7, height = 8)
  print(CombinePlots(plots))
  dev.off()
  
  StopWatchEnd$UmapAndTsneColourByEachSampleCellCluster$dataset <- Sys.time()

  ####################################
  ### Finding differentially expressed genes for each cell cluster (each sample reclustered)
  ####################################
  writeLines("\n*** Finding differentially expressed genes for each cell cluster (each sample reclustered) ***\n")
  
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
  
  StopWatchStart$FindDiffMarkers$dataset  <- Sys.time()
    
  FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.each_sample@meta.data))
    
  seurat.object.each_sample.markers <- FindAllMarkers(object = seurat.object.each_sample, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
    
  write.table(data.frame("GENE"=rownames(seurat.object.each_sample.markers),seurat.object.each_sample.markers),paste(Tempdir,"/",PrefixOutfiles, "_", dataset, ".SEURAT_MarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)

  StopWatchEnd$FindDiffMarkers$dataset  <- Sys.time()

}

####################################
### Load each sample cluster assignments
####################################
writeLines("\n*** Load each sample cluster assignments ***\n")

EachSampleReclusteredCellClusters <- data.frame(read.table(OutfileEachSampleClusters, header = T, row.names = 1, sep = "\t"))
seurat.object.integrated <- AddMetaData(object = seurat.object.integrated, metadata = EachSampleReclusteredCellClusters, col.name = "EachSampleReclusteredCellClusters")

# switch the identity class of all cells to reflect "clusters"
Idents(object = seurat.object.integrated) <- seurat.object.integrated$EachSampleReclusteredCellClusters
  
####################################
### Get average expression for each sample, for each cluster, for each gene (reclustered)
####################################
writeLines("\n*** Get average expression for each sample, for each cluster, for each gene (reclustered) ***\n")

cluster.averages<-AverageExpression(object = seurat.object.integrated, use.scale = F, use.counts = F)

OutfileClusterAveragesRNA<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_AverageGeneExpressionPerSampleReclustered_RNA.tsv", sep="")
OutfileClusterAveragesSCT<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_AverageGeneExpressionPerSampleReclustered_SCT.tsv", sep="")
OutfileClusterAveragesINT<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_AverageGeneExpressionPerSampleReclustered_integrated.tsv", sep="")
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
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir, "/" , PrefixOutfiles, ".SEURAT_UsedOptions.txt", sep="")
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

OutfileCPUusage<-paste(Tempdir, "/" , PrefixOutfiles, ".SEURAT_CPUusage.txt", sep="")
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
    ### Because the StopWatch[Start|End]$LoadScRNAseqData$dataset
    ### are using the word `dataset` itself as key instead of the dataset 'ID'
    for (dataset in rownames(InputsTable)) {
      if (regexpr("POSIXct", class(StopWatchStart[[stepToClock]][[dataset]]), ignore.case = T)[1] == 1) {
        TimeStart <- StopWatchStart[[stepToClock]][[dataset]]
        TimeEnd   <- StopWatchEnd[[stepToClock]][[dataset]]
        TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
        ReportTime<-c(paste(stepToClock, TimeDiff, dataset, sep = "\t", collapse = ""))
        write(file = OutfileCPUusage, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
      }
    }
  }
}

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")

writeLines(paste(Outdir,"/SEURAT/", sep="", collapse = ""))

outfiles_to_move <- list.files(Tempdir, pattern = PrefixOutfiles, full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""), to=paste(Outdir,"/SEURAT/", eachFile,sep=""), overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################
writeLines(paste("END - All done!!! See:\n", OutfileCPUusage, "\nfor computing times report", sep = "", collapse = ""))

quit()
