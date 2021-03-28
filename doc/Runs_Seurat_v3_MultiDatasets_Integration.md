Script name
================
`Runs_Seurat_v3_MultiDatasets_Integration.R`

Description
================
Runs Seurat version 3 scRNA-seq data integration using CCA<br />
The input is a table with paths to R objects from 'Runs_Seurat_v3_MultiDatasets_QC_Normalization.R' script

The pipeline is based on these Seurat tutorials:<br />
https://satijalab.org/seurat/v3.2/sctransform_vignette.html (SCtransform normalization)<br />
https://satijalab.org/seurat/v3.2/integration.html (general integration)<br />
https://satijalab.org/seurat/v3.2/immune_alignment.html (control vs. treatment)<br />
https://carmonalab.github.io/STACAS/tutorial.html (alternative anchor finder)<br />

It allows the user to provide parameters as one-line commands via R library(optparse).

Default parameters are set as defults based on the Seurat tutorial or from empirical observations and can be changed in list(DefaultParameters).<br />

Outfiles
================
Tables with integration data and an R objects for downstream analyzes<br />
See 'Outputs Description' below for details.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Loads each normalized dataset R object produced by script 'Runs_Seurat_v3_MultiDatasets_QC_Normalization.R'
  2. Merges and integrates datasets correcting batch effects
  3. Saves integrated datasets R object
  4. Saves log files
  5. Ends

Example commands
================
For example inputs see (SEURAT section): <br />
https://github.com/pughlab/crescent/tree/master/examples <br />

To display help commands type: <br />
`Rscript Runs_Seurat_v3_MultiDatasets_Integration.R -h`

To run the script type something like:<br />
`Rscript /path_to/Runs_Seurat_v3_MultiDatasets_QC_Normalization.R -i /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/LISTS_AND_COMMANDS/list_10X_1kPBMCs_Input_ForQcAndNormalization.tsv -j NA -k N -l N -o /path_to/SINGLE_CELL/CRESCENT/TESTING_SOFTWARE/ -p pbmcs -u MAX -s Y -w 2 -x NA -a 10000`

Inputs Description
================
a) For each dataset, an R object from script Runs_Seurat_v3_MultiDatasets_QC_Normalization.R<br />
b) A table indicating dataset types

Outputs Description
================
|  Extension     |                        Contents                        |
| ------------   |  ----------------------------------------------------- |
| PSEUDO_BULK    | Tables with a dataset comparison at pseudo-bulk level, averaging expression of genes across all cells |
| ANCHORS        | Anchors found by STACAS or Seurat |
| R_OBJECTS      | An R object with data from QC, Normalization, and this script |
| R_OBJECTS_CWL	 | An R object for the CReSCENT WebApp |
| LOG_FILES      | Log files with parameters, libraries and versions used, and computing time for each step |


Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**R [4.0.2] and the following R packages** <br /><br />
Base packages:<br />
parallel, stats4, stats, graphics, grDevices, utils, datasets, methods, base

Other attached packages:<br />
cluster_2.1.0, tidyr_1.1.2, STACAS_1.0.1, stringr_1.4.0, 
future_1.19.1, cowplot_1.1.0, ggplot2_3.3.2, data.table_1.13.0, 
fmsb_0.7.0, optparse_1.6.6, dplyr_1.0.2, Seurat_3.2.1
