
  
Script name
================
`Runs_SC3_Clustering.R`


Description
================
Runs SC3 clustering on a *file* with cell barcodes in columns and genes in rows.
The clustering procedure is based on this SC3 tutorial https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html

It allows to provide one-line commands.

Parameters for SC3 are defined based on the SC3 tutorial or from empirical observations.<br />
For example, a *pseudocount* for SingleCellExperiment() is set to 1e-99, instead of the default which is 1 (see https://goo.gl/3VzQ3L)

Parameters can be changed in section "Tailored parameters".

**Note: SC3 conducts heuristic evaluations, thus running this script with the same inputs and parameters may produce slightly different results between runs.**

 
Outfiles
================
A table with the cell clusters.<br />
\*pdf files from the PCNA and t-SNE clustering etc., as shown in the SC3 tutorial.<br />
Future versions will have an html report.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Reads data from a matrix of genes vs. cell-barcodes
  2. Flags spike-in genes
  3. Performs linear dimensional reduction by PCA
  4. Clusters the cells
  5. Runs Non-linear dimensional reduction (tSNE)
  6. Reports output
  7. Ends

Example commands
================
This example works with files provided in folder ~/examples/INPUTS<br />

To display help commands type: <br />
`Rscript Runs_SC3_Clustering.R -h`

To run the script type something like:<br />
`Rscript ~/bin/Runs_SC3_Clustering.R -i ~/path_to_/example_infile.tsv -t ~/path_to_/example_celltype.tsv -o ~/example/outfiles -p example_name -l 2 -u 10 -c 2`

Inputs Description
================

Either:<br />
a) a *file* with cell barcodes in columns and genes in rows

Example infiles are provided in folder 'examples'

Outputs Description
================
| First Header |  Second Header |
| ---------------------------------------- |  ------------------------------------------- |
| *SC3_clusters.tsv                        |  Cell barcodes/ID's per cluster              |
| *SC3_MarkersPerCluster.tsv               |  Gene markers per cluster                    |
| *SC3_plotPCA.pdf                         |  PCA plot before clustering                  |
| *SC3_plotPCA_ColourAndSizeByCluster.pdf  |  PCA plot coloured by SC3 clustering results |
| *SC3_CPUusage.tsv                        |  Stats on CPU usage, memory and time         |
| *SC3_UsedOptions.txt                     |  Options used in the run                     |


Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**R and the following R packages** <br /><br />
**SingleCellExperiment** <br />
Can be installed in R console with `source("https://bioconductor.org/biocLite.R")` <br />
`biocLite("SingleCellExperiment")` <br />
It's a dependency of SC3 <br /><br />
**SC3** <br />
Can be installed in R console with `source("https://bioconductor.org/biocLite.R")` <br />
`biocLite("SC3")` <br />
It's a wrapper of SC3 packages <br /><br />
**scater** <br />
Can be installed in R console with `install.packages('scater')`<br /><br />
**optparse**<br />
Can be installed in R console with `install.packages('optparse')`<br />
It's used to handle one-line commands<br /><br />
**data.table**<br />
Can be installed in R console with `install.packages('data.table')`<br />
It's used to read matrices faster then read.table()
