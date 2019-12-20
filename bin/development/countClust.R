suppressPackageStartupMessages(library(CountClust))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(optparse))

main <- function() {
  
  option_list <- list(
    make_option(c("-i", "--input"), default="NA",
                help="A RData file containing a SingleCellExperiment object"),
    make_option(c("-o", "--outdir"), default="NA",
                help="A path/name for the results directory"),
    make_option(c("-c", "--num_clusters"), default=8,
                help="Number of clusters")
    
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  Input<- opt$input
  Outdir<- opt$outdir
  ClusterNum<-opt$num_clusters
  
  load(Input)

  input <- as.matrix(counts(sce))
  outs <- FitGoM(t(input), K = ClusterNum, tol = 0.1, path_rda = NULL)

  colData(sce)$countClust <- unlist(apply(outs$fit$omega, 1,  which.max))
  res <- res<-data.frame(row.names(colData(sce)), unlist(apply(outs$fit$omega, 1,which.max)))
  
  write.csv(res, file=Outdir, quote=FALSE, row.names = FALSE, col.names=TRUE)

}

main()
