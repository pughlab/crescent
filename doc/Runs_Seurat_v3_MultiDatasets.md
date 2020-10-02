Script name
================
`Runs_Seurat_v3_MultiDatasets.R`

Description
================
Runs Seurat version 3 scRNA-seq data normalization, integration, batch effect correction, dimension reduction, 
cell clustering and differentil gene expression. The input is a table with paths to datasets and parameters to integrate 
datasets. Datasets can be in either MTX or TSV format (descriptions below).

The pipeline is based on these Seurat tutorials:<br />
https://satijalab.org/seurat/v3.2/sctransform_vignette.html (SCtransform normalization)<br />
https://satijalab.org/seurat/v3.2/integration.html (general integration)<br />
https://satijalab.org/seurat/v3.2/immune_alignment.html (e.g. control vs. treatment)<br />
https://carmonalab.github.io/STACAS/tutorial.html (alternative anchor finder)<br />

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
`Rscript ~/r_programs/Runs_Seurat_v3_MultiDatasets.R -i ~/path_to/list_of_datasets_and_params.tsv -j ~/path_to/list_of_barcodes_to_remove -k N -l N -y Seurat -z list_of_referece_datasets -r 0.4 -v 1,2,3 -o ~/path_to/outdir -p outfiles_prefix  -c ~/path_to/metadata.tsv -g /path_to/infile_selected_genes -m 1,2,3 -d 10 -e 0.01 -f 1,2,3,4,5,6,7,8,9,10,11,12 -b r0.4_cell_type -u MAX -s Y -w N -x NA -a 4000`

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
| FILTERED_DATA_MATRICES                   | tables with raw and normalized counts, after filtering datasets by QC parameters in -i | 
| LOG_FILES                                | tables with run commands, computing times and R libraries used |
| QC_PLOTS                                 | violin and t-SNE/UMAP plots showing QC metrics         | 
| QC_TABLES                                | tables underlying QC_PLOTS |
| R_OBJECTS                                | R object files | 
| PSEUDO_BULK                              | tables with marginal counts for each gene in each dataset |
| SELECTED_GENE_DIMENSION_REDUCTION_PLOTS  | t-SNE/UMAP plots showing selected genes|
| STACAS                                   | results from using library(STACAS) to select anchors for dataset integration |
| UNFILTERED_DATA_MATRICES                 | similar to FILTERED_DAT_MATRICES, but for unfiltered datasets |

Note: if the run uses `-w Y` then other directories with prefix `frontend_` needed by CReSCENT's graphic user interface will be produced. 

Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**R [4.0.2] and the following R packages** <br /><br />
Base packages:<br />
"parallel", "stats4", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base"

Other attached packages:<br />
"cluster_2.1.0", "tidyr_1.1.2", "STACAS_1.0.1", "loomR_0.2.0", "itertools_0.1-3", "iterators_1.0.12", "hdf5r_1.3.3", "R6_2.4.1", "gtools_3.8.2", "future_1.19.1", "cowplot_1.1.0", "ggplot2_3.3.2", "data.table_1.13.0", "fmsb_0.7.0", "optparse_1.6.6", "dplyr_1.0.2", "Seurat_3.2.1", "DropletUtils_1.8.0", "SingleCellExperiment_1.10.1", "SummarizedExperiment_1.18.2", "DelayedArray_0.14.1", "matrixStats_0.57.0", "Biobase_2.48.0", "GenomicRanges_1.40.0", "GenomeInfoDb_1.24.2", "IRanges_2.22.2", "S4Vectors_0.26.1", "BiocGenerics_0.34.0"
