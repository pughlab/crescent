Script name
================
`Runs_InferCNV.R`

Description
================
Runs InferCNV to predict copy number variants (CNVs) in a mixed sample of cells, using library(InferCNV) <br />
The script is based on https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV<br />

It allows the user to provide parameters as one-line commands via R library(optparse).

Outfiles
================
Tables and heatmaps with InferCNV Modified Gene Expression per chromosome and cell class (cell cluster or cell type)<br />
See 'Outputs Description' below for details.

General workflow
================
  1. Loads scRNA-seq data
  2. Loads cell classes
  3. Computes InferCNV
  4. Ends

Example commands
================
For example inputs see (INFERCNV section): <br />
https://github.com/pughlab/crescent/tree/master/examples <br />

To display help commands type: <br />
`Rscript Runs_InferCNV.R -h`

To run the script type something like:<br />
`Rscript /path_to/Runs_InferCNV.R -i INPUTS/glio.wGtexBrain.counts.matrix.gz -t DGE -j INPUTS/glio.wGtexBrain.sample_annots.txt -k Brain_Cerebellum,Brain_Caudate_basal_ganglia,Brain_Cortex,Brain_Nucleus_accumbens_basal_ganglia,Brain_Cerebellar_Hemisphere,Brain_Frontal_Cortex_BA9,Brain_Hippocampus -g INPUTS/gencode_v19_gene_pos.txt -m 0.1 -n 0.1 -s 0.15 -o OUTPUTS -p glio.wGtexBrain -u MAX -w 0 -a 10000`

Note: a run with ~64,000 cells requires at least 300GB of RAM, 10 cores and 75hrs of CPU usage.

Inputs Description
================
a) scRNA-seq raw data<br />
b) table indicating cell classes<br />
c) list of cell classes to be used as 'normal' cells<br />

Outputs Description
================
|               Extension                  |                        Contents                        |
| --------------------------------------   |  ----------------------------------------------------- |
| infercnv.pdf | Heatmap with Modified Gene Expression |
| *tsv and *txt | Text files with values and colours used for the heatmap |
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
Seurat_3.2.1, infercnv_1.4.0, future_1.19.1, data.table_1.13.0, optparse_1.6.6
