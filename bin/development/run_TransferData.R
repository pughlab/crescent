# Integrate the bone marrow data and use them as reference

cat("\014")
rm(list = ls())

suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(Matrix))

# options(future.globals.maxSize = 16 * 1024 * 1024^2)

option_list <- list(
  make_option(c("-r", "--Ref_input"), default="NA",
              help="The path/name to the reference *directory* with seurat object files"),

  make_option(c("-q", "--Query_input"), default="NA",
              help="The path/name to the query *directory* with seurat object list files"),

  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
)
opt <- parse_args(OptionParser(option_list=option_list))

Ref_input <- opt$Ref_input
Query_input <- opt$Query_input
Outdir <- opt$outdir

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input", "outdir")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs
####################################
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Outdir<-gsub("/$", "", Outdir)

dir.create(file.path(Outdir, "ALIGNED_MATRIX"), recursive = T)

ref <- readRDS(Ref_input)

query.list <- readRDS(Query_input)

ref <- SCTransform(ref, verbose = F)

for (dataset in names(query.list)) {
  query.list[[dataset]] <- SCTransform(query.list[[dataset]], verbose = FALSE)
}

for (dataset in names(query.list)) {
  anchor.t <- FindTransferAnchors(ref, query.list[[dataset]], normalization.method='SCT',
                                  verbose = FALSE)

  temp <- TransferData(anchor.t, ref@assays$SCT@data, verbose = FALSE)

  #save the transferred data
  mtx <- temp@data
  mtx.name <- paste(Outdir, "/", dataset, ".mtx", sep="")
  writeMM(mtx, mtx.name)
  gzip(mtx.name)
}

mtx <- ref@assays$SCT@data
mtx <- as(mtx, "dgCMatrix")
writeMM(mtx, paste(Outdir, "/", "Ref.mtx", sep=""))
gzip(paste(Outdir, "/", "Ref.mtx", sep=""))

Ref.label <- ref@meta.data$cell.type

write.table(Ref.label, file.path(Outdir, 'Labels.tsv'),
            quote = F, row.names = F, col.names = F)

write.table(names(query.list), file.path(Outdir, 'Query_names.tsv'),
            quote = F, row.names = F, col.names = F)
