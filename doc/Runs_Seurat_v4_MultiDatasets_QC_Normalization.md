Script name
================
`Runs_Seurat_v4_MultiDatasets_QC_Normalization.R`

Description
================
Runs Seurat version 4 scRNA-seq data QC and normalization<br />
The input is a table with paths to datasets and parameters for QC<br />
Datasets can be in either MTX or TSV format (descriptions below).

The pipeline is based on these Seurat tutorials:<br />
https://satijalab.org/seurat/articles/integration_introduction.html (control vs. treatment)<br />

It allows the user to provide parameters as one-line commands via R library(optparse).

Default parameters are set as defults based on the Seurat tutorial or from empirical observations and can be changed in list(DefaultParameters).<br />

Outfiles
================
Tables and plots with QC data, R objects for downstream analyzes, and files for CReSCENT's WebApp QC plots<br />
See 'Outputs Description' below for details.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Loads definition of paths to scRNA-seq datasets and QC parameters
  2. Loads scRNA-seq data and generate QC plots for each dataset
  3. Saves outfiles
  4. Saves log files
  5. Ends

Example commands
================
For example inputs see (SEURAT section): <br />
https://github.com/pughlab/crescent/tree/master/examples <br />

To display help commands type: <br />
`Rscript Runs_Seurat_v4_MultiDatasets_QC_Normalization.R -h`

To run the script type something like:<br />
`Rscript /path_to/Runs_Seurat_v4_MultiDatasets_PCA_Clustering_DimReduction.R -i /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/LISTS_AND_COMMANDS/list_10X_1kPBMCs_Input_ForPcaClusteringDimReduction.tsv -j /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/SEURAT/R_OBJECTS/pbmcs.SEURAT_Integration.rds -q NA -r 0.1 -v 1 -o /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/ -p pbmcs -c NA -g NA -m NA -y RNA,SCT -d 10 -u MAX -s Y -t RNA,SCT -w 2 -x NA -a 10000`

Inputs Description
================
One of the following input types:<br />
a) For each dataset, a *directory* with 10X files from Cell Ranger v3 (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz) <br />
b) For each dataset, a *file* in TSV <tab> delimited format with cell-barcodes in columns and genes in rows

Outputs Description
================
|               Extension                  |                        Contents                        |
| --------------------------------------   |  ----------------------------------------------------- |
| UNFILTERED_DATA_MATRICES                 | Inputted datasets in MTX format before any QC filters was applied |
| FILTERED_DATA_MATRICES                   | Inputted datasets in MTX format after QC filters were applied |
| QC_PLOTS		                   | Quality control violin plots for each dataset |
| QC_TABLES		                   | Quality control tables for each dataset |
| R_OBJECTS		                   | R objects for each dataset for downstream analysis |
| R_OBJECTS_CWL		                   | R objects for the CReSCENT WebApp |
| CRESCENT_CLOUD	                   | Files for the CReSCENT WebApp |
| LOG_FILES		                   | Log files with parameters, libraries and versions used, and computing time for each step |

Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**R [4.0.2] and the following R packages** <br /><br />
Base packages:<br />
parallel, stats4, stats, graphics, grDevices, utils, datasets, methods, base

Other attached packages:<br />
tidyr_1.1.3, future_1.21.0,
cowplot_1.1.1, ggplot2_3.3.3,
data.table_1.14.0, fmsb_0.7.1,
optparse_1.6.6, dplyr_1.0.6,
SeuratObject_4.0.1, Seurat_4.0.2,
DropletUtils_1.10.3, SingleCellExperiment_1.12.0,
SummarizedExperiment_1.20.0, Biobase_2.50.0,
GenomicRanges_1.42.0, GenomeInfoDb_1.26.7,
IRanges_2.24.1, S4Vectors_0.28.1,
BiocGenerics_0.36.1, MatrixGenerics_1.2.1,
matrixStats_0.58.0
