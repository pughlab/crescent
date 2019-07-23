Script name
================
`Runs_RCA.R`

Description
================
Runs RCA clustering on a SingleCellExperiment object.
The clustering procedure is based on this RCA tutorial https://cran.r-project.org/web/packages/RCA/RCA.pdf

It allows to provide one-line commands.

Outfiles
================
A table with the cell clusters.<br />

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Load data from a SingleCellExperiment object
  2. Filter out lowly expressed genes
  3. Normalize data
  4. Log transform data
  5. Project the expression data into reference component data. Can choose reference of "GlobalPanel", "ColonEpitheliumPanel" or "SelfProjection"
  6. Generate cell clusters using hierarchical clustering
  7. Reports output
  8. Ends

  Example commands
  ================
  This example works with files provided in folder ~/examples/INPUTS<br />

  To display help commands type: <br />
  `Rscript Runs_RCA.R -h`

  To run the script type something like:<br />
  `Rscript ~/bin/Runs_RCA.R -i ~/path_to_/example_file.RData -o ~/example/rca_clusters.csv`

  Inputs Description
  ================

  a *RData file* with a SingleCellExperiment object

  Example infiles are provided in folder 'examples'

  Outputs Description
  ================
  | First Header |  Second Header |
  | ---------------------------------------- |  ------------------------------------------- |
  | *rca_clusters.csv                        |  Cell barcodes/ID's per cluster              |


Dependencies
================

**R and the following R packages** <br /><br />
**SingleCellExperiment** <br />
Can be installed in R console with `source("https://bioconductor.org/biocLite.R")` <br />
`biocLite("SingleCellExperiment")` <br />
It's a dependency of SC3 <br /><br />
**RCA** <br />
Can be installed in R console with `install.packages("RCA", repo="https://cran.rstudio.com/")` <br />
**scater** <br />
Can be installed in R console with `install.packages('scater')`<br /><br />
**optparse**<br />
Can be installed in R console with `install.packages('optparse')`<br />
It's used to handle one-line commands<br /><br />
