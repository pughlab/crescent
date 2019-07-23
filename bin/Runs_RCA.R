
suppressPackageStartupMessages(library(RCA))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(WGCNA))

pre_clean <- function(sce) {

  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  del_gns <- c(which(rowData(sce1)$is_feature_control_rbp), which(rowData(sce1)$is_feature_control_mt))
  sce1 <- sce1[-del_gns,]
  sce1 <- scater::normalize(sce1)

  return(sce1)
}

main <- function() {

  option_list <- list(
    make_option(c("-i", "--input"), default="NA",
            help="A RData file containing a SingleCellExperiment object"),
    #
    make_option(c("-p", "--panel"), default="GlobalPanel",
            help="A reference panel for the projection step"),

    make_option(c("-n", "--normalize"), default="no_norm",
            help="A normalize method, either 'no_norm' or 'scQ' "),

    make_option(c("-m", "--minSize"), default=5,
            help="Minimum size of the resulting clusters"),

    make_option(c("-s", "--split"), default=1,
            help="Integer indicating how deep to cut dendrogram, possible values are 0,1,2,3 or 4"),

    make_option(c("-o", "--outdir"), default="NA",
                help="A path/name for the results directory")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  Input<- opt$input
  Outdir<- opt$outdir
  Panel <- opt$panel
  Normaize <- opt$normalize
  MinSize <- opt$minSize
  Split <- opt$split

  ####################################
  ### Check that mandatory parameters are not 'NA' (default)
  ####################################

  ListMandatory<-list("input", "outdir")
  for (param in ListMandatory) {
    if (length(grep('^NA$',opt[[param]], perl = T))) {
      stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h       for help.", sep = "", collapse = ""))
     }
   }
  load(Input)

  sce1 <- pre_clean(sce)

  obj <- as.matrix(calculateCPM(sce1))
  rownames(obj) <- rowData(sce1)$mat.gene_name
  data_obj <- dataConstruct(obj)
  data_obj <- geneFilt(obj_in = data_obj, method = "default")
  data_obj <- cellNormalize(data_obj, method = "no_norm")
  data_obj <- dataTransform(data_obj, method = "log10")
  data_obj <- featureConstruct(data_obj, method = Panel)
  set.seed(20742579)
  data_obj <- cellClust(data_obj, method = "hclust", deepSplit_wgcna = Split,
                        min_group_Size_wgcna = MinSize)

  res <- data.frame(row.names(colData(sce)), data_obj$group_labels_color$groupLabel)
  colnames(res) <- c("cell", "rca_clusters")
  write.csv(res, file=Outdir, quote=FALSE, row.names = FALSE, col.names=TRUE)

}

main()
