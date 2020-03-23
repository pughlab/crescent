# Example files for CReSCENT

CELL RANGER
================

We are using Cell Ranger v3.X

Documentation of cellranger can be found here:
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

Inputs: <br />
Takes fastq files, each lane has three fastq files (R1, R2 and I1), which can be Gb in size. <br />
You can download example files from 10X: <br />
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

Relevant outputs:
Include three files with the sparse matrix in Matrix Market format with the mapped read counts per gene, per cell barcode: <br />
`outs/filtered_feature_bc_matrices/` <br />
barcodes.tsv.gz <br />
features.tsv.gz <br />
matrix.mtx.gz <br />


SEURAT
================

Script name: Runs_Seurat_v3.R

Infiles: <br />
Either the Matrix Market (MTX) files from Cell Ranger `outs/filtered_feature_bc_matrices/`, or a matrix with genes (rows) and cell barcodes (columns), called Digital Gene Expression (DGE) matrix.

Outfiles: <br />
See `https://github.com/jdime/crescent/doc/Runs_Seurat_v3.md`


INTER-CONVERTING MTX AND DGE FORMAT FILES
================

Download example MTX files from (10X):
http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz

_How to run scripts:_

To convert a DGE file into MTX files use: <br />
`Rscript ~/path_to/obtains_mtx_files_from_gene_x_barcode_matrix.R -i ~/path_to/gene_vs_barcode_matrix.tsv -o ~/path_to_store_outfiles/ -p prefix_for_log_files -l n` <br />
Note: the DGE file may be gzipped

To convert MTX files into a DGE file use: <br />
`Rscript ~/path_to/obtains_gene_x_barcode_matrix_from_mtx_files.R -i ~/path_to/filtered_feature_bc_matrices/ -o ~/path_to_store_outfiles/ -p prefix_for_log_files_outfiles` <br />
Note: Cell Ranger v2 produces unzipped files and the features.tsv.gz file is called genes.tsv


OLDER VERSIONS
================

10X provides example datasets: <br />
https://support.10xgenomics.com/single-cell-gene-expression/datasets/

