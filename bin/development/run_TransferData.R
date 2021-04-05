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

#Read the reference and query Suerat object

ref <- readRDS(Ref_input)
query.list <- readRDS(Query_input)

#Do SCTransform
ref <- SCTransform(ref, verbose = F)

for (dataset in names(query.list)) {
  query.list[[dataset]] <- SCTransform(query.list[[dataset]], verbose = FALSE)
}

#Align each query sample to the reference dataset, and save the aligned matrix in mtx format
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

# Save the reference matrix in mtx format. Note that we don't need to save this
# file every time we want to predict the labels of any query dataset. For those public
# datasets we suggest users to use as reference, we can save a copy of this mtx file
# in the server and use it directly in the future.
mtx <- ref@assays$SCT@data
mtx <- as(mtx, "dgCMatrix")
writeMM(mtx, paste(Outdir, "/", "Ref.mtx", sep=""))
gzip(paste(Outdir, "/", "Ref.mtx", sep=""))

# Obtain the labels of reference dataset and save it.
Ref.label <- ref@meta.data$cell.type

write.table(Ref.label, file.path(Outdir, 'Labels.tsv'),
            quote = F, row.names = F, col.names = F)

# The Query_names.tsv file contains the name of the 
write.table(names(query.list), file.path(Outdir, 'Query_names.tsv'),
            quote = F, row.names = F, col.names = F)
