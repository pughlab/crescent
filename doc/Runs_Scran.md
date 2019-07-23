Script name
================
`Runs_Scran.R`

Description
================
Runs scran clustering on a SingleCellExperiment object.
The clustering procedure is based on this RCA tutorial https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html

It allows to provide one-line commands.

Outfiles
================
A table with the cell clusters.<br />

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Load data from a SingleCellExperiment object
  2. Normalize data using a deconvolution strategy for scaling normalization
  3. Compute distance matrix
  4. Hierarchical Clustering
  5. Generated cell clusters using Dynamic Tree Cutting
  6. Reports output
  7. Ends

  Example commands
  ================
  This example works with files provided in folder ~/examples/INPUTS<br />

  To display help commands type: <br />
  `Rscript Runs_Scran.R -h`

  To run the script type something like:<br />
  `Rscript ~/bin/Runs_Scran.R -i ~/path_to_/example_file.RData -o ~/example/scran_clusters.csv`

  Inputs Description
  ================

  a *RData file* with a SingleCellExperiment object

  Example infiles are provided in folder 'examples'

  Outputs Description
  ================
  | First Header |  Second Header |
  | ---------------------------------------- |  ------------------------------------------- |
  | *scran_clusters.csv                        |  Cell barcodes/ID's per cluster              |


Dependencies
================

**R and the following R packages** <br /><br />
**SingleCellExperiment** <br />
Can be installed in R console with `source("https://bioconductor.org/biocLite.R")` <br />
`biocLite("SingleCellExperiment")` <br />
It's a dependency of SC3 <br /><br />
**SCRAN** <br />
Can be installed in R console with `BiocManager::install("scran")` <br />
**scater** <br />
Can be installed in R console with `install.packages('scater')`<br /><br />
**optparse**<br />
Can be installed in R console with `install.packages('optparse')`<br />
It's used to handle one-line commands<br /><br />
