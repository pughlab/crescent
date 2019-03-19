####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script made based on https://satijalab.org/seurat/pbmc3k_tutorial_v3.html
###
### NEW IMPLEMENTATIONS SINCE Seurat v2:
### 1) Rewritten with Seurat v3 commands (including ability to read output from Cell Ranger v3)
###    Main differences vs. Seurat v2 include:
###    a) new function names
###    b) a new function DimPlot() is used for 'uma', 'tsne' and 'pca' plots, instead of PCAPlot() and TSNEPlot()
###       Hence, now we define 'dimensionss' instead of principal components
###    c) since Cell Ranger v3 allows now to have multiple features (not only genes),
###       in general all references to 'genes' in v2 are now called 'features'
### 2) Default parameters are indexed using list(DefaultParameters) instead of variable names. This allows list(DefaultParameters)
###    to be reported in a *log file
### 3) Implemented QC violin plots directly in ggplots instead of Seurat's VlnPlot() to have more control of layout, legends, titles, etc.
### 4) Using R base plot() to create feature-vs-feature scatter plots instead of FeatureScatter()
### 5) Implemented all ggplots and Seurat plotting functions (which are based on ggplots) with print() function, like:
###    `print(FeaturePlot(...))` instead of `FeaturePlot(...)` alone
###    Otherwise using ggplots and Seurat plots inside if/else loops cause errors
###    https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-do-lattice_002ftrellis-graphics-not-work_003f
### 
### THINGS TO DO:
### 1) Pick the right number of dimension components
###    E.g. try to automatically get the inflection point from the PCElbowPlot() function output
### 2) Pick the right resolution from FindClusters()
###    E.g. implement Iness and Bader paper https://f1000research.com/articles/7-1522/v1 approach
###    By picking the number of clusters based on differentially expressed genes
### 3) Add a lists of ENSEMBL Ids for mitochondrial genes instead of just MT- and mt- (at gene names)
###    Need to do it for both Human and Mouse
### 6) In 'Colour t-SNE plots showing each requested gene' need to automate to get a layout() for more than 3 genes in multiple rows
### 7) In Seurat v2, Suluxan reported that the -c example Javier provided called a duplicated row names error
###
### THINGS NICE TO HAVE:
### 1) Assigning cell type identity to clusters (needs supervised annotations, maybe based on GSVA)
### 2) Use knitr() to produce html plot layouts (https://yihui.name/knitr/demo/stitch/)
### 3) In "Load data" we use Seurat(Read10X). In this command, when all barcodes come from the same sample (i.e. finish with the same digit), like:
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
### Required libraries
####################################
suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
### Seurat v3 can be installed like:
### install.packages('devtools')
### devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(fmsb))         # to calculate the percentages of extra properties to be t-SNE plotted
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table - only needed is using '-t DGE'
suppressPackageStartupMessages(library(ggplot2))      # (CRAN) to generate QC violin plots
suppressPackageStartupMessages(library(cowplot))      # (CRAN) to arrange QC violin plots and top legend
### library(staplr)     only if using option '-s y', note it needs pdftk
####################################

####################################
### Required external packages
####################################
### 'pdftk'   to merge selected *pdf files into *summary_plots.pdf using libary(staplr)
###           and to add statistics to violin plots overlapping pdf files
###           in Mac you can install it with something like:
###           'sudo port install pdftk'
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
              help="Either the path/name to a 10X *directory* with barcodes.tsv, genes.tsv and matrix.mtx files;
                or path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns"),
#
  make_option(c("-t", "--input_type"), default="NA",
              help="Indicates if input is either a '10X' directory or a 'DGE' file"),
#
  make_option(c("-r", "--resolution"), default="1",
              help="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain
                a larger (smaller) number of communities"),
#
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
#
  make_option(c("-p", "--prefix_outfiles"), default="NA",
            help="A prefix for outfile names, e.g. your project ID"),
#
  make_option(c("-s", "--summary_plots"), default="y",
              help="Indicates if a *summary_plots.pdf file should be generated [use 'y'] or not [use 'n']
                Note this needs 'pdftk' and R library(staplr)"),
#
  make_option(c("-c", "--infile_colour_tsne"), default="NA",
              help="A <tab> delimited table of barcodes and discrete properties to colour the t-SNE, like:
                Barcode              CellClass    InOtherDatasets
                AAACCTGAGCGGCTTC-1   1            yes
                AAACCTGAGTCGAGTG-1   1            no
                AAACCTGCAAAGGAAG-1   2            yes
                AAACCTGGTCTCATCC-1   2            no"),
#
  make_option(c("-g", "--list_genes"), default="NA",
              help="A <comma> delimited list of gene identifiers whose expression will be mapped into the t-SNE plots"),
#
  make_option(c("-a", "--opacity"), default="0.1",
              help="If using a --list_genes, this parameter provides a value for the minimal opacity of gene expression. Use a value between 0 and 1"),
#
  make_option(c("-d", "--pca_dimensions"), default="10",
              help="Max value of PCA dimensions to use for clustering and t-SNE functions
                FindClusters(..., dims.use = 1:-d) and RunTSNE(..., dims.use = 1:-d)
                Typically '10' is enough, if unsure use '10' and afterwards check these two files:
                *JackStraw*pdf, use the number of PC's where the solid curve shows a plateau along the dotted line, and
                *PCElbowPlot.pdf, use the number of PC's where the elbow shows a plateau along the y-axes low numbers"),
#
 make_option(c("-m", "--percent_mito"), default="0,0.05",
             help="<comma> delimited min,max number of percentage of mitochondrial gene counts in a cell to be included in normalization and clustering analyses
                For example, for regular scRNA-seq use '0,0.2', or for Nuc-seq use '0,0.05'"),
#
 make_option(c("-n", "--n_genes"), default="200,8000",
             help="<comma> delimited min,max number of unique gene counts in a cell to be included in normalization and clustering analyses. E.g '200,8000'"),
#
  make_option(c("-e", "--return_threshold"), default="0.01",
              help="For each cluster only return markers that have a p-value < return_thresh,  e.g. '0.01'")
)

opt <- parse_args(OptionParser(option_list=option_list))

Input            <- opt$input
InputType        <- opt$input_type
Resolution       <- as.numeric(opt$resolution) ## using as.numeric avoids FindClusters() to crash by inputting it as.character [default from parse_args()]
Outdir           <- opt$outdir
PrefixOutfiles   <- opt$prefix_outfiles
SummaryPlots     <- opt$summary_plots
InfileColourTsne <- opt$infile_colour_tsne
ListGenes        <- opt$list_genes
Opacity          <- as.numeric(opt$opacity)
PcaDimsUse       <- c(1:as.numeric(opt$pca_dimensions))
ListPMito        <- opt$percent_mito
ListNGenes        <- opt$n_genes
ThreshReturn     <- as.numeric(opt$return_threshold)

PrefixOutfiles <- c(paste(PrefixOutfiles,"_res",Resolution,sep=""))
Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

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
### Start stopwatch
####################################

StartTimeOverall <-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input", "input_type", "outdir", "prefix_outfiles")
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
if (regexpr("^10X$", InputType, ignore.case = T)[1] == 1) {
  print("Loading 10X infiles")
  input.matrix <- Read10X(data.dir = Input)
}else if (regexpr("^DGE$", InputType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  input.matrix <- data.frame(fread(Input),row.names=1)
}else{
  stop(paste("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}
dim(input.matrix)

####################################
### Create a Seurat object
####################################
writeLines("\n*** Create a Seurat object ***\n")
seurat.object.u  <- CreateSeuratObject(counts = input.matrix, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = PrefixOutfiles)
nCellsInOriginalMatrix<-length(seurat.object.u@meta.data$orig.ident)

####################################
### Get  mitochondrial genes
####################################
writeLines("\n*** Get  mitochondrial genes ***\n")
mitoRegExpressions<- paste(c("^MT", "^mt"),collapse = "|")
mito.features <- grep(pattern = mitoRegExpressions, x = rownames(x = seurat.object.u), value = T)
percent.mito <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
seurat.object.u[['percent.mito']] <- percent.mito

####################################
### Filter cells based gene counts and mitochondrial representation
####################################
writeLines("\n*** Filter cells based gene counts and mitochondrial representation ***\n")
seurat.object.f<-subset(x = seurat.object.u, subset = nFeature_RNA > DefaultParameters$MinGenes & nFeature_RNA < DefaultParameters$MaxGenes & percent.mito > DefaultParameters$MinPMito & percent.mito < DefaultParameters$MaxPMito)
seurat.object.u
seurat.object.f

####################################
### QC EDA violin plots
####################################
writeLines("\n*** QC EDA violin plots ***\n")

### Get unfiltered data QC statistics
nFeature_RNA.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$nFeature_RNA, nGenes = 1)
nCount_RNA.u.df    <-data.frame(Expression_level = seurat.object.u@meta.data$nCount_RNA,   nCount_RNA = 1)
percent.mito.u.df  <-data.frame(Expression_level = seurat.object.u@meta.data$percent.mito, percent.mito = 1)
#
nFeature_RNAStats.u<-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object.u@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.u  <-paste(c( "mean = ",round(mean(seurat.object.u@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object.u@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.u     <-paste(c(" mean = ",round(mean(seurat.object.u@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object.u@meta.data[,"percent.mito"]),3)), sep = "", collapse="")

### Get filtered data QC statistics
nFeature_RNA.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$nFeature_RNA, nGenes = 2)
nCount_RNA.f.df    <-data.frame(Expression_level = seurat.object.f@meta.data$nCount_RNA,   nCount_RNA = 2)
percent.mito.f.df  <-data.frame(Expression_level = seurat.object.f@meta.data$percent.mito, percent.mito = 2)
#
nFeature_RNAStats.f<-paste(c(" mean = ",round(mean(seurat.object.f@meta.data[,"nFeature_RNA"]),0),"\n", "median = ",round(median(seurat.object.f@meta.data[,"nFeature_RNA"]),0)), sep = "", collapse="")
nCount_RNAStats.f  <-paste(c(" mean = ",round(mean(seurat.object.f@meta.data[,"nCount_RNA"]),0),  "\n", "median = ",round(median(seurat.object.f@meta.data[,"nCount_RNA"]),0)),   sep = "", collapse="")
percent.mito.f     <-paste(c(" mean = ",round(mean(seurat.object.f@meta.data[,"percent.mito"]),3),"\n", "median = ",round(median(seurat.object.f@meta.data[,"percent.mito"]),3)), sep = "", collapse="")

### Put QC statistics together
nFeature_RNA.m.df  <-data.frame(rbind(nFeature_RNA.u.df,nFeature_RNA.f.df))
nCount_RNA.m.df    <-data.frame(rbind(nCount_RNA.u.df,nCount_RNA.f.df))
percent.mito.m.df  <-data.frame(rbind(percent.mito.u.df,percent.mito.f.df))
LabelUnfiltered    <-paste("Before filters: No. of cells = ", nrow(seurat.object.u@meta.data), sep ="", collapse = "")
LabelFiltered      <-paste("After filters:  No. of cells = ", nrow(seurat.object.f@meta.data), sep ="", collapse = "")

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
  labs(x="Mitochond. %") +
  annotate("text", x = 1 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.u, col = ColoursQCViolinPlots[[1]]) +
  annotate("text", x = 2 , y = max(percent.mito.m.df$Expression_level)*1.1, label = percent.mito.f, col = ColoursQCViolinPlots[[2]])

bottom_row<-plot_grid(nFeature_RNA.plot, nCount_RNA.plot, percent.mito.plot, ncol = 3)

### Create a *pdf file with the violin ggplot's

VlnPlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_QC_VlnPlot.pdf", sep="")
pdf(file=VlnPlotPdf, width = 10, height = 7)
print(plot_grid(Headers.plot, bottom_row, ncol = 1, rel_heights = c(0.2,1)))
dev.off()

####################################
### Feature-vs-feature scatter plots
####################################
writeLines("\n*** Feature-vs-feature scatter plots ***\n")

UnfilteredData.df<-data.frame(nCount_RNA = seurat.object.u@meta.data$nCount_RNA,
                              nGene = seurat.object.u@meta.data$nFeature_RNA,
                              percent.mito = seurat.object.u@meta.data$percent.mito,
                              filtered_out = colnames(seurat.object.u) %in% colnames(seurat.object.f))
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = ColoursQCViolinPlots[[1]])
UnfilteredData.df$DotColour<-gsub(x=UnfilteredData.df$DotColour,    pattern = FALSE, replacement = ColoursQCViolinPlots[[2]])
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$filtered_out, pattern = TRUE,  replacement = 4)
UnfilteredData.df$DotPch   <-gsub(x=UnfilteredData.df$DotPch,       pattern = FALSE, replacement = 16)

FeatureVsFeaturePlotPdf<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_NumbReadsVsNumbGenesAndMito_VlnPlot.pdf", sep="")
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

####################################
### Remove the Unfiltered seurat object
####################################
writeLines("\n*** Remove the Unfiltered seurat object ***\n")
rm(seurat.object.u)
rm(UnfilteredData.df)

####################################
### Normalize data
####################################
writeLines("\n*** Normalize data ***\n")
seurat.object.f <- NormalizeData(object = seurat.object.f, normalization.method = "LogNormalize", scale.factor = DefaultParameters$ScaleFactor)

####################################
### Detect, save list and plot variable genes
####################################
writeLines("\n*** Detect, save list and plot variable genes ***\n")
### Note: Seurat developers recommend to set default parameters to mark visual outliers on the dispersion plot
### x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
### but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy.
seurat.object.f <- FindVariableFeatures(object = seurat.object.f, selection.method = 'mean.var.plot', mean.cutoff = c(DefaultParameters$XLowCutoff, DefaultParameters$XHighCutoff), dispersion.cutoff = c(DefaultParameters$YCutoff, Inf))
VariableGenes<-VariableFeatures(object = seurat.object.f)
length(VariableGenes)

write(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VariableGenes.txt", sep=""), x=VariableGenes)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VariableGenes.pdf", sep=""))
print(VariableFeaturePlot(object = seurat.object.f, cols = c("blue", "red")))
dev.off()

####################################
### Scale data and remove unwanted sources of variation such as:
### cell cycle stage,  batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data)
####################################
writeLines("\n*** Scale data and remove unwanted sources of variation ***\n")
seurat.object.f <- ScaleData(object = seurat.object.f, vars.to.regress = c("nCount_RNA", "percent.mito"), features = rownames(x = seurat.object.f), display.progress=F)

####################################
### Perform linear dimensional reduction by PCA
### Examine and visualize PCA results a few different ways
### Note: DimPlot can now handle 'umap' and 'tsne' in addition to 'pca', but for 'umap'
### you must first install the umap-learn python package (e.g. via pip install umap-learn)
### https://github.com/satijalab/seurat/blob/master/R/dimensional_reduction.R
####################################
writeLines("\n*** Perform linear dimensional reduction by PCA ***\n")
seurat.object.f <- RunPCA(object = seurat.object.f, features = VariableGenes, verbose = T, do.print = T, ndims.print = DefaultParameters$PrintPCA.PcsPrint, nfeatures.print = DefaultParameters$PrintPCA.GenesPrint)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VizPCA.pdf", sep=""), width=7, height=10)
print(VizDimLoadings(object = seurat.object.f, reduction = "pca", dims = DefaultParameters$VizPCA.PcsUse, nfeatures = DefaultParameters$VizPCA.nGenesToPlot))
dev.off()

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCAPlot.pdf", sep=""))
print(DimPlot(object = seurat.object.f, dims = c(1,2), reduction = "pca") + theme(legend.position="none"))
dev.off()

seurat.object.f <- ProjectDim(object = seurat.object.f, overwrite = T, verbose = T, nfeatures.print = 10)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCHeatmap.C1to",DefaultParameters$PCHeatmapComponentsToPlot,".pdf", sep=""), width=7, height=12)
print(DimHeatmap(object = seurat.object.f, dims = 1:DefaultParameters$PCHeatmapComponentsToPlot, cells = DefaultParameters$PCHeatmapCellsUse, balanced = T))
dev.off()

####################################
### Determine statistically significant principal components
####################################
writeLines("\n*** Determine statistically significant principal components ***\n")
### NOTE: JackStraw() process can take a long time for big datasets
### More approximate techniques such PCElbowPlot() can be used to reduce computation time

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_PCElbowPlot.pdf", sep=""))
print(ElbowPlot(object = seurat.object.f))
dev.off()

####################################
### Cluster the cells
####################################
writeLines("\n*** Cluster the cells ***\n")
StartTimeClustering<-Sys.time()
options(scipen=10) ## Needed to avoid an 'Error in file(file, "rt") : cannot open the connection'
seurat.object.f <- FindNeighbors(object = seurat.object.f, dims = PcaDimsUse) ## This step was part of FindClusters() in Seurat v2
seurat.object.f <- FindClusters(object = seurat.object.f, reduction.type = "pca", resolution = Resolution, print.output = 0, save.SNN = T)
EndTimeClustering<-Sys.time()

CellNames<-rownames(seurat.object.f@meta.data)
ClusterIdent<-seurat.object.f@meta.data$RNA_snn_res.1
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

####################################
### Get average expression for each cluster for each gene
####################################
writeLines("\n*** Get average expression for each cluster for each gene ***\n")

cluster.averages<-AverageExpression(object = seurat.object.f, use.scale = F, use.counts = F)
OutfileClusterAverages<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster.tsv", sep="")
Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("r", Resolution, "_C", names(cluster.averages$RNA), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
write.table(data.frame(cluster.averages$RNA),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Run Non-linear dimensional reduction (tSNE)
####################################
writeLines("\n*** Run Non-linear dimensional reduction (tSNE) ***\n")
### NOTE: if the datasets is small you may get
### "Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : Perplexity is too large."
### One can try tunning down the default RunTSNE(..., perplexity=30) to say 5 or 10

seurat.object.f <- RunTSNE(object = seurat.object.f, dims.use = PcaDimsUse, do.fast = T)

pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot.pdf", sep=""), width = 7, height = 7)
print(DimPlot(object = seurat.object.f, reduction = 'tsne', group.by = 'ident', label = T, label.size=10))
dev.off()

####################################
### Saving the R object
####################################
writeLines("\n*** Saving the R object ***\n")
OutfileRDS<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_object.rds", sep="")
saveRDS(seurat.object.f, file = OutfileRDS)

####################################
### Write out t-SNE coordinates
####################################
writeLines("\n*** Write out t-SNE coordinates ***\n")
OutfileTsneCoordinates<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNECoordinates.tsv", sep="")

Headers<-paste("Barcode",paste(colnames(seurat.object.f@reductions$tsne@cell.embeddings),sep="",collapse="\t"),sep="\t",collapse = "\t")
write.table(Headers,file = OutfileTsneCoordinates, row.names = F, col.names = F, sep="\t", quote = F)
write.table(seurat.object.f@reductions$tsne@cell.embeddings, file = OutfileTsneCoordinates,  row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Colour t-SNE by nFeature_RNA, nCount_RNA, and percent.mito
####################################
writeLines("\n*** Colour t-SNE by nFeature_RNA, nCount_RNA, and percent.mito ***\n")
CellPropertiesToTsne<-c("nFeature_RNA", "nCount_RNA", "percent.mito")
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_QC.pdf", sep=""), width = 16, height = 5)
print(FeaturePlot(object = seurat.object.f, features = CellPropertiesToTsne, cols = c("lightgrey", "blue"), reduction = "tsne", ncol = 3, pt.size = 1.5))
dev.off()

####################################
### Colour t-SNE by -infile_colour_tsne
####################################
writeLines("\n*** Colour t-SNE by -infile_colour_tsne ***\n")
if (regexpr("^NA$", InfileColourTsne, ignore.case = T)[1] == 1) {
  print("No extra barcode-attributes will be used for t-SNE plots")
}else{
  seurat.object.meta.data<-seurat.object.f@meta.data
  ExtraCellProperties <- data.frame(read.table(InfileColourTsne, header = T, row.names = 1))
  head(ExtraCellProperties)

  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --infile_colour_tsne and --input can match each other
  rownames(ExtraCellProperties)<-gsub(x =rownames(ExtraCellProperties), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object.f <- AddMetaData(object = seurat.object.f, metadata = ExtraCellProperties)

  # Generating outfile
  # Note DimPlot() takes the entire current device (pdf)
  # even if using layout(matrix(...)) or  par(mfrow=())
  # Thus each property t-SNE is written to a separate page of a single *pdf outfile
  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_ExtraProperties.pdf", sep=""), height = DefaultParameters$BaseSizeSingleTnePlot, width = DefaultParameters$BaseSizeSingleTnePlot)
  for (property in colnames(ExtraCellProperties)) {
    print(DimPlot(object = seurat.object.f, reduction = 'tsne', group.by = property, combine = T, legend = "none") + ggtitle(property))
    }
  dev.off()
}

####################################
### Colour t-SNE plots showing each requested gene
####################################
writeLines("\n*** Colour t-SNE plots showing each requested gene ***\n")

### To program layout() for more than 3 genes in multiple rows
if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
  print("No selected genes for t-SNE plots")
}else{
  ListOfGenesForTsnes<-unlist(strsplit(ListGenes, ","))
  if (length(ListOfGenesForTsnes) <= 4) {
    pdfWidth  <- (length(ListOfGenesForTsnes) * DefaultParameters$BaseSizeMultipleWidth)
    pdfHeight <- DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- length(ListOfGenesForTsnes)
  }else{
    pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
    pdfHeight <- (as.integer(length(ListOfGenesForTsnes) / 4) + 1) * DefaultParameters$BaseSizeMultipleHeight
    nColFeaturePlot <- 4
  }

  pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_SelectedGenes.pdf", sep=""), width=pdfWidth, height=pdfHeight)
  print(FeaturePlot(object = seurat.object.f, ncol = nColFeaturePlot, features = c(ListOfGenesForTsnes), cols = c("lightgrey", "blue"), reduction = "tsne"))
  dev.off()
}

####################################
### Finding differentially expressed genes (cluster biomarkers)
####################################
writeLines("\n*** Finding differentially expressed genes (cluster biomarkers) ***\n")
### Finding markers for every cluster compared to all remaining cells
### only.pos allows to report only the positive ones
### Using min.pct and thresh.use (renamed logfc.threshold in latest Seurat versions) to speed comparisons up. Other options to further speed are min.diff.pct, and  max.cells.per.ident
### See http://satijalab.org/seurat/de_vignette.html
#########
### Other options
### find all markers of cluster 1
### cluster1.markers <- FindMarkers(object = seurat.object.f, ident.1 = 1, min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster1.markers, n = 5))
###
### find all markers distinguishing cluster 5 from clusters 0 and 3
### cluster5.markers <- FindMarkers(object = seurat.object.f, ident.1 = 5, ident.2 = c(0, 3), min.pct = DefaultParameters$FindAllMarkers.MinPct)
### print(x = head(x = cluster5.markers, n = 5))
########
### NOTES:
### 1) FindAllMarkers() uses return.thresh = 0.01 as defaults, but FindMarkers() displays all genes passing previous filters.
###    Thus to make the outputs between these two commands identical to each other use return.thresh = 1
###
### 2) Default pseudocount.use=1, which sounds high for a Log level correction
###    An earlier version of this script was using 1e-99, but it was probably too small
###    Now using the inverse of the number of cells in the data.
###    This is sufficiently small as to not compress logGER magnitudes,
###    while keeping comparisons with zero reasonably close to the range of potential logGER values (Innes and Bader, 2018, F1000 Research)
#########

FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.f@meta.data))

StartTimeFindAllMarkers<-Sys.time()
seurat.object.markers <- FindAllMarkers(object = seurat.object.f, only.pos = T, min.pct = DefaultParameters$FindAllMarkers.MinPct, return.thresh = ThreshReturn, logfc.threshold = DefaultParameters$FindAllMarkers.ThreshUse, pseudocount.use = FindMarkers.Pseudocount)
EndTimeFindAllMarkers<-Sys.time()

write.table(data.frame("GENE"=rownames(seurat.object.markers),seurat.object.markers),paste(Tempdir,"/",PrefixOutfiles,".SEURAT_MarkersPerCluster.tsv",sep=""),row.names = F,sep="\t",quote = F)
### Get top-2 genes sorted by cluster, then by p-value
top_genes_by_cluster_for_tsne<-(seurat.object.markers %>% group_by(cluster) %>% top_n(DefaultParameters$NumberOfGenesPerClusterToPlotTsne, avg_logFC))
NumberOfClusters<-length(unique(seurat.object.markers[["cluster"]]))

####################################
### Violin plots for top genes
####################################
writeLines("\n*** Violin plots for top genes ***\n")
NumberOfPanesForFeaturesPlot<-(NumberOfClusters*DefaultParameters$NumberOfGenesPerClusterToPlotTsne)
top_genes_by_cluster_for_tsne.list<-top_genes_by_cluster_for_tsne[c(1:NumberOfPanesForFeaturesPlot),"gene"][[1]]
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5

### Here need to program a better way to control the y-axis labels
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_VlnPlot_AfterClusters.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(VlnPlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), slot = 'counts', log = T, adjust = 1, pt.size = 0.5))
dev.off()

####################################
### t-SNE plots showing each cluster top genes
####################################
writeLines("\n*** t-SNE plots showing each cluster top genes ***\n")
pdfWidth  <- 4 * DefaultParameters$BaseSizeMultipleWidth
pdfHeight <- NumberOfClusters * DefaultParameters$BaseSizeMultipleHeight / 2
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_EachTopGene.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(FeaturePlot(object = seurat.object.f, ncol = 4, features = c(top_genes_by_cluster_for_tsne.list), cols = c("lightgrey", "blue"), reduction = "tsne"))
dev.off()

####################################
### Cell clusters heatmap
####################################
writeLines("\n*** Cell clusters heatmap ***\n")
top_genes_by_cluster_for_heatmap <- seurat.object.markers %>% group_by(cluster) %>% top_n(n = DefaultParameters$NumberOfGenesPerClusterToPlotHeatmap, wt = avg_logFC)
pdfWidth<-7
pdfHeight<-NumberOfClusters*1.5
pdf(file=paste(Tempdir,"/",PrefixOutfiles,".SEURAT_Heatmap.pdf", sep=""), width=pdfWidth, height=pdfHeight)
print(DoHeatmap(object = seurat.object.f, features = top_genes_by_cluster_for_heatmap$gene, label = T, group.bar = T, raster = F, angle = 0) + NoLegend() + ggtitle("Cell clusters"))
dev.off()

####################################
### Create summary plots outfile
####################################
writeLines("\n*** Create summary plots outfile ***\n")

if (regexpr("^y$", SummaryPlots, ignore.case = T)[1] == 1) {
  suppressPackageStartupMessages(library(staplr))
  SummaryPlotsPdf<-paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_summary_plots.pdf", sep = "", collapse = "")
  ListOfPdfFilesToMerge<-c(paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_QC_VlnPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_Heatmap.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_TSNEPlot.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_VlnPlot_AfterClusters.pdf", sep = "", collapse = ""),
                           paste(Tempdir,"/",PrefixOutfiles, ".SEURAT_TSNEPlot_EachTopGene.pdf", sep = "", collapse = "")
                           )
  if (regexpr("^NA$", ListGenes, ignore.case = T)[1] == 1) {
    print("Create summary file")
  }else{
    ListOfPdfFilesToMerge<-c(ListOfPdfFilesToMerge, paste(Tempdir,"/",PrefixOutfiles,".SEURAT_TSNEPlot_SelectedGenes.pdf", sep="", collapse = ""))
    print("Create summary file including t-SNE's for selected genes")
  }
  
  ### Stappling *pdf's
  staple_pdf(input_directory = NULL, input_files = ListOfPdfFilesToMerge, output_filepath = SummaryPlotsPdf)
}

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")
OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report time used
####################################
writeLines("\n*** Report time used ***\n")

EndTimeOverall<-Sys.time()

TookTimeClustering     <-format(difftime(EndTimeClustering,     StartTimeClustering,     units = "min"))
TookTimeFindAllMarkers <-format(difftime(EndTimeFindAllMarkers, StartTimeFindAllMarkers, units = "min"))
TookTimeOverall        <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "min"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_CPUusage.txt", sep="")
ReportTime<-c(
  paste("clustering",TookTimeClustering,collapse = "\t"),
  paste("FindAllMarkers",TookTimeFindAllMarkers,collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"/SEURAT/",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir, pattern = paste(PrefixOutfiles, ".SEURAT_", sep=""), full.names = F)
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
print("END - All done!!! Took time:")
print(ReportTime)

quit()

