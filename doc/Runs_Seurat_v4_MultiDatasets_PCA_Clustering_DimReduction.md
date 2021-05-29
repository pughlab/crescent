Script name
================
`Runs_Seurat_v4_MultiDatasets_PCA_Clustering_DimReduction.R`

Description
================
Runs Seurat version 4 scRNA-seq data PCA, dimension reduction (UMAP and t-SNE), and Clustering<br />
The input is an R object from 'Runs_Seurat_v4_MultiDatasets_Integration.R' script

The pipeline is based on these Seurat tutorials:<br />
https://satijalab.org/seurat/articles/integration_introduction.html<br />

It allows the user to provide parameters as one-line commands via R library(optparse).

Default parameters are set as defults based on the Seurat tutorial or from empirical observations and can be changed in list(DefaultParameters).<br />

Outfiles
================
Tables with integration data and an R objects for downstream analyzes<br />
See 'Outputs Description' below for details.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Loads integrated datasets R object produced by script `Runs_Seurat_v4_MultiDatasets_Integration.R'
  2. Process integrated datasets as a whole, including:
     - dimension reduction
     - 'global' cell clustering
     - UMAP/tSNE plots by global cell clusters, sample, sample type, requested genes and metadata
     - average gene expression
  3. For each dataset on its own, uses global clusters and repeats analysis from step 2
  4. For each dataset on its own, re-clusters cells and repeats analysis from step 2
  5. For each dataset type on its own, uses global clusters and repeats analysis from step 2
  6. For each dataset type on its own, re-clusters cells and repeats analysis from step 2
  7. Saves datasets R object with PCA, clustering and dimension reduction
  8. Saves log files
  9. Ends

Example commands
================
For example inputs see (SEURAT section): <br />
https://github.com/pughlab/crescent/tree/master/examples <br />

To display help commands type: <br />
`Rscript Runs_Seurat_v4_MultiDatasets_PCA_Clustering_DimReduction.R -h`

To run the script type something like:<br />
`Rscript /path_to/Runs_Seurat_v4_MultiDatasets_DGE.R -i /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/LISTS_AND_COMMANDS/list_10X_1kPBMCs_Input_ForDGE.tsv -j /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/SEURAT/R_OBJECTS/pbmcs.SEURAT_PCA_Clustering_DimReduction.rds -o /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/ -p pbmcs -c /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/METADATA/pbmcs.metadata.tsv -b SampleType_Manual -f 8 -e 0.01 -g 0.25 -d /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/METADATA/pbmcs.subclasses_to_compare.tsv -t SCT -k Y -u MAX -s N -w 2 -x NA -a 10000`

Inputs Description
================
a) An R object from script Runs_Seurat_v3_MultiDatasets_Integration.R<br />
b) A table indicating dataset types

Outputs Description
================
|               Extension                  |                        Contents                        |
| --------------------------------------   |  ----------------------------------------------------- |
| AVERAGE_GENE_EXPRESSION_TABLES	   | Tables with average gene expression for each gene, for each cell cluster |
| CELL_CLUSTER_IDENTITIES 		   | Cell cluster identity tables |
| CELL_FRACTIONS 			   | Tables with number and fraction of cells per cluster, per dataset |
| DIMENSION_REDUCTION_COORDINATE_TABLES    | Tables with UMAP and t-SNE coordinates |
| DIMENSION_REDUCTION_PLOTS		   | UMAP and t-SNE plots coloured by cell clusters and metadata |
| QC_PLOTS				   | UMAP and t-SNE plots coloured by QC metrics |
| PSEUDO_BULK				   | Tables with a dataset comparison at pseudo-bulk level for each dataset |
| SELECTED_GENE_DIMENSION_REDUCTION_PLOTS  | UMAP and t-SNE plots coloured by selected genes |
| R_OBJECTS				   | An R object with data from QC, Normalization, Integration, and this script |
| R_OBJECTS_CWL				   | An R object for the CReSCENT WebApp |
| LOOM_FILES_CWL			   | Loom files with raw and normalized data for the CReSCENT WebApp|
| CRESCENT_CLOUD			   | Other files for the CReSCENT WebApp |
| LOG_FILES				   | Log files with parameters, libraries and versions used, and computing time for each step |


Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**R [4.0.2] and the following R packages** <br /><br />
Base packages:<br />
parallel, stats4, stats, graphics, grDevices, utils, datasets, methods, base

Other attached packages:<br />
tidyr_1.1.32, loomR_0.2.1.9000, hdf5r_1.3.3, R6_2.4.1, 
future_1.21.01, cowplot_1.1.1, ggplot2_3.3.3, data.table_1.14.0, 
fmsb_0.7.1, optparse_1.6.6, dplyr_1.0.6, Seurat_4.0.2
