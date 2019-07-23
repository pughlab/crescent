
suppressPackageStartupMessages(library(SIMLR))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(RSpectra))
suppressPackageStartupMessages(library(RcppEigen))
suppressPackageStartupMessages(library(optparse))

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

  set.seed(345454654)
  cluster_num <- SIMLR_Estimate_Number_of_Clusters(X = as.matrix(counts(sce)), NUMC = 2:20, cores.ratio = 1)
  cluster_num <- (min(which.min(cluster_num$K1), which.min(cluster_num$K2))+1)
  output <-  SIMLR_Large_Scale(X = as.matrix(counts(sce)), c = cluster_num, k = 10, kk = 100)

  colData(sce)$SIMLR <- output$y$cluster
  res <- data.frame(row.names(colData(sce)),output$y$cluster)
  colnames(res) <- c("cell", "simlr_clusters")

  write.csv(res, file=Outdir, quote=FALSE, row.names = FALSE, col.names=TRUE)

}

main()
