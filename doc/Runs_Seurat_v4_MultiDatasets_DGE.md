Script name
================
`Runs_Seurat_v4_MultiDatasets_DGE.R`

Description
================
Runs Seurat version 4 scRNA-seq data QC and normalization<br />
The input is a table with paths to datasets and parameters for QC<br />
Datasets can be in either MTX or TSV format (descriptions below).

The pipeline is based on these Seurat tutorials:<br />
https://satijalab.org/seurat/articles/integration_introduction.html<br />

It allows the user to provide parameters as one-line commands via R library(optparse).

Default parameters are set as defults based on the Seurat tutorial or from empirical observations and can be changed in list(DefaultParameters).<br />

Outfiles
================
Tables with differentially expressed genes<br />
See 'Outputs Description' below for details.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Loads integrated datasets R object produced by script `Runs_Seurat_v4_MultiDatasets_PCA_Clustering_DimReduction.R`
  2. Computes differential gene expression (DGE) for options specified by --diff_gene_expr_comparisons
     Note, all clustering and metadata groups must be provided by --infile_r_object and --infile_metadata
  3. Saves log files
  5. Ends

Example commands
================
For example inputs see (SEURAT section): <br />
https://github.com/pughlab/crescent/tree/master/examples <br />

To display help commands type: <br />
`Rscript Runs_Seurat_v4_MultiDatasets_DGE.R -h`

To run the script type something like:<br />
`Rscript /path_to/Runs_Seurat_v4_MultiDatasets_DGE.R -i /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/LISTS_AND_COMMANDS/list_10X_1kPBMCs_Input_ForDGE.tsv -j /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/SEURAT/R_OBJECTS/pbmcs.SEURAT_PCA_Clustering_DimReduction.rds -o /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/ -p pbmcs -c /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/METADATA/pbmcs.metadata.tsv -b SampleType_Manual -f 8 -e 0.01 -g 0.25 -d /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/METADATA/pbmcs.subclasses_to_compare.tsv -t SCT -k Y -u MAX -s N -w 2 -x NA -a 10000`

Inputs Description
================
a) An R object from script Runs_Seurat_v3_MultiDatasets_PCA_Clustering_DimReduction.R<br />
b) A table indicating dataset types

Outputs Description
================
|               Extension                  |                        Contents                        |
| --------------------------------------   |  ----------------------------------------------------- |
| DIFFERENTIAL_GENE_EXPRESSION_TABLES      | Inputted datasets in MTX format before any QC filters was applied |
| R_OBJECTS				   | An R object with data from QC, Normalization, Integration, PCA, Dim. Reduction, Clustering and this script |
| R_OBJECTS_CWL				   | An R object for the CReSCENT WebApp |
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
tidyr_1.1.3, future_1.21.0, data.table_1.14.0, fmsb_0.7.1, 
optparse_1.6.6, dplyr_1.0.6, Seurat_4.0.2
