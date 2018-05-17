
  
Script name
================
`Runs_Seurat_Clustering.R`


Description
================
Runs Seurat's clustering on either 10X files from cellranger (barcodes.tsv, genes.tsv and matrix.mtx)
or a Dropseq *file* with cell barcodes in columns and genes in rows.
The clustering procedure is based on this Seurat tutorial http://satijalab.org/seurat/pbmc3k_tutorial.html

It allows to provide one-line commands.

Parameters for Seurat are defined based on the Seurat tutorial or from empirical observations.<br />
For example, *pseudocount.use* for FindAllMarkers() is set to 1e-99, instead of the default is 1 (see https://goo.gl/3VzQ3L)

Parameters can be changed in section "Tailored parameters".
 
Outfiles
================
A table with the cell clusters.<br />
Are \*pdf files from the t-SNE, violin plots, clustering etc., as shown in the Seurat tutorial.<br />
Future versions will have an html report.

General workflow
================
This code is a wrapper library written in R and the general workflow of this script is as follows:
  1. Reads data from 10X or Dropseq
  2. Flags mitochondrial genes
  3. Creates 'violin' and scatter plots before and after filtering data (by gene counts and mitochondrial representation)
  4. Normalizes data
  5. Detects variable genes
  6. Performs linear dimensional reduction by PCA
  7. Determine statistically significant principal components
  8. Clusters the cells
  9. Runs Non-linear dimensional reduction (tSNE)
  10. Findins differentially expressed genes (cluster biomarkers)
  11. Creates 'violin' plots for top differentially expressed genes for each cluster
  12. Creates heatmaps of differentially expressed genes
  13 . Reports output
  14. Ends

Example commands
================
This example works with files provided in folder ~/examples/INPUTS<br />

To display help commands type: <br />
`Rscript Runs_Seurat_Clustering.R -h`

To run the script type something like:<br />
`Rscript ~/bin/Runs_Seurat_Clustering.R -i ~/path_to_/filtered_gene_bc_matrices -t 10X -o ~/example/outfiles -p example_10X`

Inputs Description
================

Either:<br />
a) a 10X *directory* cointaining outfiles from cellranger (barcodes.tsv, genes.tsv and matrix.mtx), or<br />
b) a Dropseq *file* with cell barcodes in columns and genes in rows

Example infiles are provided in folder 'examples'

Outputs Description
================
| Extension |  Contents |
| ------------------------------ |  ---------------------  |
| *CellClusters.tsv              |  Cell clusters          |
| *MarkersPerCluster.tsv         |  Markers per cluster    |
| *CPUusage.tsv                  |  CPU time usage         |
| *UsedOptions.txt               |  Options used in run    |
| *GenePlot.pdf                  |  GenePlot               |
| *GenePlot.seurat_filtered.pdf  |  GenePlot after filters |
| *Heatmap.pdf                   |  Heatmap                |
| *JackStraw.C1toC12.pdf         |  JackStraw              |
| *PCAPlot.pdf                   |  PCAPlot                |
| *PCElbowPlot.pdf               |  PCElbowPlot            |
| *PCHeatmap.C1.pdf              |  PCHeatmap cluster 1    |
| *PCHeatmap.C1toN.pdf           |  PCHeatmap all clusters |
| *TSNEPlot.pdf                  |  TSNEPlot               |
| *TSNEPlot_EachTopGene.pdf      |  TSNEPlot_EachTopGene   |
| *VariableGenes.pdf             |  VariableGenes          |
| *VariableGenes.txt             |  VariableGenes          |
| *VizPCA.pdf                    |  VizPCA                 |
| *VlnPlot.pdf                   |  VlnPlot                |
| *VlnPlot.seurat_filtered.pdf   |  VlnPlot                |
| *VlnPlot_AfterClusters.pdf     |  VlnPlot_AfterClusters  |


Authors
================

**Javier Diaz (https://github.com/jdime)**

Dependencies
================

**3) R and the following R packages** <br />
**Seurat** <br />
Can be installed in R console with `install.packages('Seurat')`<br /><br />
**dplyr** <br />
Can be installed in R console with `install.packages('dplyr')`<br /><br />
**optparse**<br />
Can be installed in R console with `install.packages('optparse')`<br />
It's used to handle one-line commands<br /><br />
**data.table**<br />
Can be installed in R console with `install.packages('data.table')`<br />
It's used to read Dropseq matrices faster then read.table()
