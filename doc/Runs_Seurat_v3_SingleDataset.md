Script name
================
`Runs_Seurat_v3_SingleDataset.R`

Description
================
Runs Seurat version 3 scRNA-seq data QC, normalization, dimension reduction, cell clustering and differential gene expression. The input is a matrix in either MTX or TSV format (descriptions below). 

The clustering procedure is based on this Seurat tutorial https://satijalab.org/seurat/pbmc3k_tutorial_v3.html

It allows the user provide parameters as one-line commands via R library(optparse).

Default parameters are set as defults based on the Seurat tutorial or from empirical observations and can be changed in list(DefaultParameters).<br />

Outfiles
================
Tables with QC data, cell clusters, differential gene expression (DGE), t-SNE and UMAP plots, DGE violin plots, etc.
See 'Outputs Description' below for details.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Loads scRNA-seq data and generate QC plots
  2. Normalizes and scales measurements
  3. Runs dimension reduction
  4. Determines cell clusters
  5. Runs average gene expression
  6. Draws UMAP/tSNE plots using cell clusters, requested genes and metadata
  7. Runs DGE
  8. Saves log files
  9. Ends

Example commands
================
This example works with files provided in folder ~/examples/INPUTS<br />

To display help commands type: <br />
`Rscript Runs_Seurat_v3.R -h`

To run the script type something like:<br />
`Rscript ~/bin/Runs_Seurat_v3.R -i ~/path_to/filtered_feature_bc_matrix/ -t MTX -j NA -b 2 -k Y -r 1 -o ~/example/outfiles -p sample_ID -c ~/path_to/example_cell_type.tsv -g ~/path_to/example_selected_genes.txt -d 10 -f Y -m 0,0.5 -q 0.075 -n 50,8000 -v 1,80000 -e 0.01 -u AUTO -s Y -w N -a AUTO`

Inputs Description
================

One of the following input types:<br />
a) a *directory* with 10X files from Cell Ranger v2 (barcodes.tsv, genes.tsv and matrix.mtx) <br />
b) a *directory* with 10X files from Cell Ranger v3 (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz) <br />
c) a *file* in TSV <tab> delimited format with cell-barcodes in columns and genes in rows (e.g. from DropSeq)

Example infiles are provided in folder 'examples'



Outputs Description
================
|               Extension                  |                        Contents                        |
| --------------------------------------   |  ----------------------------------------------------- |
| QC_PLOTS                                 | violin and t-SNE/UMAP plots showing QC metrics         | 
| QC_TABLES                                | tables underlying QC_PLOTS |
| FILTERED_DATA_MATRICES                   | MTX files with filtered raw and normalized values |
| DIMENSION_REDUCTION_PLOTS                | t-SNE/UMAP plots showing cell clusters and metadata | 
| DIMENSION_REDUCTION_COORDINATE_TABLES    | tables underlying DIMENSION_REDUCTION_PLOTS | 
| SELECTED_GENE_DIMENSION_REDUCTION_PLOTS  | t-SNE/UMAP plots showing selected genes| 
| CELL_CLUSTER_IDENTITIES                  | table with cell cluster identities | 
| AVERAGE_GENE_EXPRESSION_TABLES           | table with each gene average expression for each cell cluster | 
| DIFFERENTIAL_GENE_EXPRESSION_TABLES      | table with DGE for each cell cluster vs. rest of cells in the dataset | 
| DIFFERENTIAL_GENE_EXPRESSION_TOP_2_GENE_PLOTS  | t-SNE/UMAP plots showing each cluster top-2 DGE genes|  
| SUMMARY_PLOTS                            | selection of plots typically used for publications |  
| R_OBJECTS                                | R object files | 
| LOG_FILES                                | tables with run commands, computing times and R libraries used |

Note: if the run uses `-w Y` then other directories with prefix `frontend_` needed by CReSCENT's graphic user 
interface will be produced. 


Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies [brackets indicate tested versions]
================

**R [3.6.1] and the following R packages** <br /><br />
**Seurat [3.1.1]** <br />
Can be installed in R console with: <br />
`install.packages('devtools')`<br />
`devtools::install_github("satijalab/seurat@v3.1.1")`<br /><br />
**DropletUtils [1.6.1]** <br />
Can be installed in R console with: <br />
`if (!requireNamespace("BiocManager", quietly = TRUE))`<br />
`install.packages("BiocManager")`<br />
`BiocManager::install("DropletUtils")`<br /><br />
**dplyr [0.8.3]** <br />
Can be installed in R console with: <br />
`install.packages('dplyr')`<br /><br />
**optparse [1.6.4]**<br />
Can be installed in R console with: <br />
`install.packages('optparse')`<br />
It's used to handle one-line commands<br /><br />
**data.table [1.12.2]**<br />
Can be installed in R console with: <br />
`install.packages('data.table')`<br />
It's used to read Gene Expression matrices faster than read.table()<br /><br />
**staplr [2.9.0]**<br />
Can be installed in R console with: <br />
`install.packages('staplr')`<br />
It's used to manipulate \*pdf files. It requires pdftk, which can be obtained from<br />
https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/<br /><br />
**fmsb [0.6.3]**<br />
Can be installed in R console with: <br />
`install.packages('fmsb')`<br />
It's used to calculate the percentages of extra properties to be t-SNE plotted<br /><br />
**ggplot2 [3.2.1]**<br />
Can be installed in R console with: <br />
`install.packages('ggplot2')`<br />
It's used to generate QC violin plots<br /><br />
**cowplot [1.0.0]**<br />
Can be installed in R console with: <br />
`install.packages('cowplot')`<br />
It's used arrange QC violin plots and top legend<br /><br />
**future [1.15.1]**<br />
Can be installed in R console with: <br />
`install.packages('future')`<br />
It's used parallelize time consuming functions like: RunPCA<br /><br />
**loomR [0.2.0]**<br />
Can be installed in R console with: <br />
`devtools::install_github(repo = "mojaveazure/loomR")` <br />
It's needed for front-end display of data. Only needed if using `-w Y`
