
  
Script name
================
`Runs_Seurat_Clustering.R`


Description
================
Runs Seurat's clustering on either 10X files from cellranger (barcodes.tsv, genes.tsv and matrix.mtx)
or a *file* with the Gene Expression Matrix with cell-barcodes in columns and genes in rows.
The clustering procedure is based on this Seurat tutorial http://satijalab.org/seurat/pbmc3k_tutorial.html

It allows to provide one-line commands.

Parameters for Seurat are defined based on the Seurat tutorial or from empirical observations.<br />
For example, *pseudocount.use* for FindAllMarkers() is set to 1e-99, instead of the default is 1 (see https://goo.gl/3VzQ3L)

Parameters can be changed in section "Tailored parameters".
 
Outfiles
================
A table with the cell-clusters, a table with differentially expressed genes on each cell-cluster.<br />
Plots provided as \*pdf files from the t-SNE, violin plots, clustering etc., as shown in the Seurat tutorial.<br />
See 'Outputs Description' below for details

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Read data from 10X directory or a matrix file
  2. Flag mitochondrial genes
  3. Create 'violin' and scatter plots before and after filtering data (by gene counts and mitochondrial representation)
  4. Normalize data
  5. Detect variable genes
  6. Perform linear dimensional reduction by PCA
  7. Determine statistically significant principal components
  8. Cluster the cells
  9. Run Non-linear dimensional reduction (tSNE)
  10. Find differentially expressed genes (each cluster markers)
  11. Create 'violin' plots for top differentially expressed genes for each cluster
  12. Create heatmaps of differentially expressed genes
  13. Create a file with summary plots
  13 . Report output
  14. End

Example commands
================
This example works with files provided in folder ~/examples/INPUTS<br />

To display help commands type: <br />
`Rscript Runs_Seurat_Clustering.R -h`

To run the script type something like:<br />
`Rscript ~/bin/Runs_Seurat_Clustering.R -i ~/path_to_/filtered_gene_bc_matrices -t 10X -o ~/example/outfiles -p example_10X -r 1 -e 0.01 -d 10 -s y -g MALAT1`

Inputs Description
================

Either:<br />
a) a 10X *directory* cointaining outfiles from cellranger (barcodes.tsv, genes.tsv and matrix.mtx), or<br />
b) a Dropseq *file* with cell barcodes in columns and genes in rows

Example infiles are provided in folder 'examples'

Outputs Description
================
| Extension |  Contents |
| ------------------------------ |  -----------------------------------------------------   |
| *CellClusters.tsv              |  Cell clusters                                           |
| *MarkersPerCluster.tsv         |  Markers per cluster                                     |
| *CPUusage.tsv                  |  CPU time usage                                          |
| *UsedOptions.txt               |  Options used in run                                     |
| *GenePlot.pdf                  |  Gene plot                                               |
| *GenePlot.seurat_filtered.pdf  |  Gene plot after filters                                 |
| *Heatmap.pdf                   |  Heatmap top genes, all clusters                         |
| *JackStraw.C1toC12.pdf         |  JackStraw plot                                          |
| *PCAPlot.pdf                   |  PCA plot                                                |
| *PCElbowPlot.pdf               |  PCElbow plot                                            |
| *PCHeatmap.C1.pdf              |  PCHeatmap cluster 1                                     |
| *PCHeatmap.C1toN.pdf           |  PCHeatmap all clusters                                  |
| *TSNEPlot.pdf                  |  t-SNE plot                                              |
| *TSNEPlot_EachTopGene.pdf      |  t-SNE plot mapping top-2 genes for each cluster         |
| *VariableGenes.pdf             |  Variable genes plot                                     |
| *VariableGenes.txt             |  Variable genes list                                     |
| *VizPCA.pdf                    |  Vizualize genes assciated with PCA                      |
| *VlnPlot.pdf                   |  QC violin plot                                          |
| *VlnPlot.seurat_filtered.pdf   |  QC violin plot after applying FilterCells()             |
| *VlnPlot_AfterClusters.pdf     |  Violing plot of top-2 genes for each cluster            |
| *TSNEPlot_SelectedGenes.pdf    |  t-SNE plot mapping genes selected by option -g          |
| *summary_plots.pdf             |  Summary plots                                           |


Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**R and the following R packages** <br />
**Seurat** <br />
Can be installed in R console with `install.packages('Seurat')`<br /><br />
**dplyr** <br />
Can be installed in R console with `install.packages('dplyr')`<br /><br />
**optparse**<br />
Can be installed in R console with `install.packages('optparse')`<br />
It's used to handle one-line commands<br /><br />
**data.table**<br />
Can be installed in R console with `install.packages('data.table')`<br />
It's used to read Gene Expression matrices faster than read.table()<br /><br />
**staplr**<br />
Can be installed in R console with `install.packages('staplr')`<br />
It's used to manipulate *pdf files. It requires pdftk, which can be obtained from<br />
https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/<br /><br />
