Script name
================
`Runs_Seurat_v3.R`


Description
================
Runs Seurat version 3 scRNA-seq data normalization, dimension reduction and cell clustering of
Cell Ranger Matrix Market files or a Digital Expression Matrix (see *Inputs Description* below)

The clustering procedure is based on this Seurat tutorial https://satijalab.org/seurat/pbmc3k_tutorial_v3.html

It allows the user provide parameters as one-line commands.

Other parameters are set as defults based on the Seurat tutorial or from empirical observations and can be changed in list(DefaultParameters).<br />

Outfiles
================
A table with the cell-clusters, a table with differentially expressed genes on each cell-cluster, plots provided as
\*pdf files from the t-SNE, violin plots, clustering etc. See 'Outputs Description' below for details.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Reads data from 10X directory or a DGE matrix file
  2. Flags mitochondrial genes
  3. Creates Quality Control (QC) 'violin' and scatter plots before and after filtering data
  4. Creates scatter plots of No. of reads vs. mitochondrial representation and, No. of reads vs. No. of genes
  4. Normalizes data
  5. Detects variable genes
  6. Performs linear dimensional reduction by PCA
  7. Determines statistically significant principal components
  8. Clusters the cells
  9. Gets average gene expression per cluster
  10. Runs Non-linear dimensional reduction (tSNE)
  11. Finds differentially expressed genes (gene markers for each cell cluster)
  12. Creates 'violin' plots for top differentially expressed genes for each cell cluster
  13. Creates heatmaps for top differentially expressed genes for each cell cluster
  14. Creates t-SNE plots for top differentially expressed genes for each cell cluster
  15. Creates t-SNE plots for requested barcode-attributes (optional)
  16. Creates t-SNE plots for selected genes (optional)
  17. Creates a file with summary plots
  18. Reports used options
  19. Ends

Example commands
================
This example works with files provided in folder ~/examples/INPUTS<br />

To display help commands type: <br />
`Rscript Runs_Seurat_Clustering.R -h`

To run the script type something like:<br />
`Rscript ~/bin/Runs_Seurat_v3.R -i ~/path_to_/filtered_feature_bc_matrix/ -t 10X -r 1 -o ~/example/outfiles 
-l 10X_PBMC -s y -c ~/path_to_/example_cell_type.tsv -g MALAT1,GAPDH -a 0.3 -d 10 -n 50,8000 -e 0.01 -p sample_ID`

Inputs Description
================

Either:<br />
a) a *directory* with 10X files from Cell Ranger v2 (barcodes.tsv, genes.tsv and matrix.mtx) <br />
b) a *directory* with 10X files from Cell Ranger v3 (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz) <br />
c) a *file* with the Digital Gene Expression (DGE) matrix with cell-barcodes in columns and genes in rows (e.g. from DropSeq)

Example infiles are provided in folder 'examples'

Outputs Description
================
| Extension |  Contents |
| -------------------------------------- |  ----------------------------------------------------- |
| *CellClusters.tsv                      |  Cell clusters                                         |
| *AverageGeneExpressionPerCluster.tsv   |  Average gene expression per cell cluster              |
| *MarkersPerCluster.tsv                 |  Markers per cluster                                   |
| *CPUusage.tsv                          |  CPU time usage                                        |
| *UsedOptions.txt                       |  Options used in run                                   |
| *GenePlot.pdf                          |  Gene plot                                             |
| *GenePlot.seurat_filtered.pdf          |  Gene plot after filters                               |
| *Heatmap.pdf                           |  Heatmap top genes, all clusters                       |
| *JackStraw.C1toC12.pdf                 |  JackStraw plot                                        |
| *PCAPlot.pdf                           |  PCA plot                                              |
| *PCElbowPlot.pdf                       |  PCElbow plot                                          |
| *PCHeatmap.C1.pdf                      |  PCHeatmap cluster 1                                   |
| *PCHeatmap.C1toN.pdf                   |  PCHeatmap all clusters                                |
| *TSNEPlot.pdf                          |  t-SNE plot                                            |
| *TSNEPlot_EachTopGene.pdf              |  t-SNE plot mapping top-2 genes for each cluster       |
| *TSNECoordinates.tsv                   |  t-SNE plot coordinates                                |
| *VariableGenes.pdf                     |  Variable genes plot                                   |
| *VariableGenes.txt                     |  Variable genes list                                   |
| *VizPCA.pdf                            |  Vizualize genes assciated with PCA                    |
| *VlnPlot.pdf                           |  QC violin plot                                        |
| *VlnPlot.seurat_filtered.pdf           |  QC violin plot after applying FilterCells()           |
| *VlnPlot_AfterClusters.pdf             |  Violing plot of top-2 genes for each cluster          |
| *TSNEPlot_SelectedGenes.pdf            |  t-SNE plot mapping genes selected by option -g        |
| *summary_plots.pdf                     |  Summary plots                                         |


Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**R and the following R packages** <br />
**Seurat version 3** <br />
Can be installed in R console with: <br />
`install.packages('devtools')`<br />
`devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')`<br /><br />
**dplyr** <br />
Can be installed in R console with: <br />
`install.packages('dplyr')`<br /><br />
**optparse**<br />
Can be installed in R console with: <br />
`install.packages('optparse')`<br />
It's used to handle one-line commands<br /><br />
**data.table**<br />
Can be installed in R console with: <br />
`install.packages('data.table')`<br />
It's used to read Gene Expression matrices faster than read.table()<br /><br />
**staplr**<br />
Can be installed in R console with: <br />
`install.packages('staplr')`<br />
It's used to manipulate *pdf files. It requires pdftk, which can be obtained from<br />
https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/<br /><br />
**fmsb**<br />
Can be installed in R console with: <br />
`install.packages('fmsb')`<br />
It's used to calculate the percentages of extra properties to be t-SNE plotted<br /><br />
**ggplot2**<br />
Can be installed in R console with: <br />
`install.packages('ggplot2')`<br />
It's used to generate QC violin plots<br /><br />
**cowplot**<br />
Can be installed in R console with: <br />
`install.packages('cowplot')`<br />
It's used arrange QC violin plots and top legend<br /><br />