Script name
================
`Runs_GSVA.R`

Description
================
Runs GSVA to predict cell type labels to cell clusters using average gene expression and cell type gene sets. Infiles are described in section 'Inputs Description'<br />
The script is based on Diaz-Mejia F1000Res 2019 (https://f1000research.com/articles/8-296/v3)<br />

It allows the user to provide parameters as one-line commands via R library(optparse).

Outfiles
================
Tables with cell type predictions<br />
See 'Outputs Description' below for details.

General workflow
================
  1. Loads average gene expression
  2. Loads gene sets
  3. Computes GSVA
  4. Ends

Example commands
================
For example inputs see (GSVA section): <br />
https://github.com/pughlab/crescent/tree/master/examples <br />

To display help commands type: <br />
`Rscript Runs_GSVA.R -h`

To run the script type something like:<br />
`Rscript /path_to/Runs_GSVA.R -i INFILES/pbmcs.SEURAT_AverageGeneExpression_GlobalClustering_AllDatasets_SCT.tsv.bz2 -t DGE -c INFILES/LM22_signature.cutoff3000.gmt -o OUTFILES -p pbmcs -e 0.05 -f 0.1 -u MAX -w 0`

Inputs Description
================
a) Average gene expression table with genes in rows and cell clusters in columns, e.g. from `Runs_Seurat_v3_MultiDatasets_PCA_Clustering_DimReduction.R`<br />
b) Gene sets in GMT format, with one cell type per row<br />

Outputs Description
================
|               Extension                  |                        Contents                        |
| --------------------------------------   |  ----------------------------------------------------- |
| GSVA/*all_scores_table.tsv | Table with all scores, including GSVA enrichment score (ES), p-value and FDR |
| GSVA/*enrichment_scores.tsv | Matrix with ES for cell cluster (rows) and cell types (columns) |
| GSVA/*enrichment_scores_sorted.tsv | Same data as GSVA_enrichment_scores.tsv but sorted by similarity of ES profiles |
| GSVA/*cluster_final_label.tsv | Cell type label from top ES for each cell cluster | 
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
Seurat_3.2.1, cluster_2.1.0, qvalue_2.20.0, GSVA_1.36.2, GSA_1.03.1, data.table_1.13.0, future_1.19.1, optparse_1.6.6
