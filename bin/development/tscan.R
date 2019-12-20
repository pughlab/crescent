suppressPackageStartupMessages(library(TSCAN))
suppressPackageStartupMessages(library(SingleCellExperiment))
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

  procdata <- preprocess(as.matrix(counts(sce)), minexpr_percent = 0.01, 
                         clusternum = NULL, takelog = TRUE, logbase = 2,
                         pseudocount = 1, minexpr_value = 1,
                         cvcutoff = 1)
  set.seed(4343646)
  lpsmclust <- exprmclust(procdata, clusternum = 2:20, modelNames = "VVV", reduce = T)
  colData(sce)$TSCAN <- lpsmclust$clusterid
  res <- data.frame(row.names(colData(sce)), lpsmclust$clusterid)

  write.csv(res, file=Outdir, quote=FALSE, row.names = FALSE, col.names=TRUE)

}

main()
