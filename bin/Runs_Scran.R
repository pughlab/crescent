
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dynamicTreeCut))

pre_clean <- function(sce) {

  ave.counts <- rowMeans(as.matrix(counts(sce)))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  del_gns <- c(which(rowData(sce1)$is_feature_control_rbp), which(rowData(sce1)$is_feature_control_mt))
  sce1 <- sce1[-del_gns,]

  return(sce1)
}

main <- function() {
  option_list <- list(
    make_option(c("-i", "--input"), default="NA",
                help="A RData file containing a SingleCellExperiment object"),
    #
    make_option(c("-o", "--outdir"), default="NA",
                help="A path/name for the results directory")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  Input<- opt$input
  Outdir<- opt$outdir

  load(Input)

  sce1 <- pre_clean(sce)
  clusters <- scran::quickCluster(counts(sce1), use.ranks=FALSE)
  sce1 <- scran::computeSumFactors(sce1, cluster=clusters)
  sce1 <- normalize(sce1)
  my.dist <- dist(t(counts(sce1)))
  set.seed(13115645)
  my.tree <- hclust(my.dist, method = "ward.D2")
  my.clusters <- unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0))

  colData(sce)$scran<-my.clusters
  res <- data.frame(row.names(colData(sce)), my.clusters)
  colnames(res) <- c("cell", "scran_clusters")

  write.csv(res, file=Outdir, quote=FALSE, row.names = FALSE, col.names=TRUE)

}

main()
