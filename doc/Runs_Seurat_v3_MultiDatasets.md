Script name
================
`Runs_Seurat_v3_MultiDatasets.R`

Description
================
Runs Seurat version 3 scRNA-seq data normalization, integration, batch effect correction, dimension reduction, 
cell clustering and differentil gene expression. The input is a table with paths to datasets and parameters to integrate 
datasets. Datasets can be in either MTX or TSV format (descriptions below).

The pipeline is based on this Seurat tutorial https://satijalab.org/seurat/v3.1/integration.html

It allows the user provide parameters as one-line commands via R library(optparse).

Default parameters are set as defults based on the Seurat tutorial or from empirical observations and can be changed in list(DefaultParameters).<br />

Outfiles
================
Tables with QC data, cell clusters, differential gene expression (DGE), t-SNE and UMAP plots, DGE violin plots, etc.
See 'Outputs Description' below for details.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Loads definition of paths to scRNA-seq datasets and parameters to integrate them
  2. Loads scRNA-seq data and generate QC plots for each dataset
  3. Normalizes and scales each dataset
  4. Merges and integrates datasets correcting batch effects
  5. Process integrated datasets as a whole, including: dimension reduction, 'global' cell clustering, UMAP/tSNE plots, DGE, average gene expression
  6. For each dataset on its own uses global clusters and repeats analysis from step 3
  7. For each dataset on its own re-clusters cells and repeats analysis from step 3
  8. For each dataset type on its own uses global clusters and repeats analysis from step 3
  9. For each dataset type on its own re-clusters cells and repeats analysis from step 3
  10. Obtains DGE for metadata classes
  11. Saves log files
  12. Ends

Example commands
================
This example works with files provided in folder ~/examples/INPUTS<br />

To display help commands type: <br />
`Rscript Runs_Seurat_v3_MultiDatasets.R -h`

To run the script type something like:<br />
`Rscript ~/r_programs/Runs_Seurat_v3_MultiDatasets.R -i ~/path_to/list_of_datasets_and_params.tsv -j ~/path_to/list_of_barcodes_to_remove -k reference_datasets_list -r 0.4 -v 1,2,3 -o ~/path_to/outdir -p outfiles_prefix  -c ~/path_to/metadata.tsv -g GENE1,GENE2 -d 10 -e 0.01 -f 1,2,3,4,5,6,7,8,9,10,11,12 -b r0.4_cell_type -u MAX -s Y -w N -a 4000`

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
| AVERAGE_GENE_EXPRESSION_TABLES           | table with each gene average expression for each cell cluster | 
| CELL_CLUSTER_IDENTITIES                  | table with cell cluster identities | 
| CELL_FRACTIONS                           | table with cells for each cluster, for each dataset | 
| DIFFERENTIAL_GENE_EXPRESSION_TABLES      | table with DGE for each cell cluster vs. rest of cells in the dataset | 
| DIMENSION_REDUCTION_COORDINATE_TABLES    | tables underlying DIMENSION_REDUCTION_PLOTS | 
| DIMENSION_REDUCTION_PLOTS                | t-SNE/UMAP plots showing cell clusters and metadata | 
| LOG_FILES                                | tables with run commands, computing times and R libraries used |
| QC_PLOTS                                 | violin and t-SNE/UMAP plots showing QC metrics         | 
| QC_TABLES                                | tables underlying QC_PLOTS |
| R_OBJECTS                                | R object files | 
| SELECTED_GENE_DIMENSION_REDUCTION_PLOTS  | t-SNE/UMAP plots showing selected genes| 

Note: if the run uses `-w Y` then other directories with prefix `frontend_` needed by CReSCENT's graphic user interface will be produced. 

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
