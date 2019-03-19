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
`outs/filtered_gene_bc_matrices/Reference_Transcriptome/` <br />
barcodes.tsv.gz <br />
features.tsv.gz <br />
matrix.mtx.gz <br />

SEURAT
================

Script name: Runs_Seurat_v3.R

Infiles: <br />
Either the Matrix Market files from Cell Ranger `outs/filtered_gene_bc_matrices/Reference_Transcriptome/` or a matrix with genes (rows) and cell barcodes (columns), called Digital Gene Expression (DGE) matrix.

Outfiles: <br />
See `https://github.com/jdime/crescent/doc/Runs_Seurat_v3.md`

OLDER VERSIONS
================

10X provides v1 and v2 chemistry example datasets:
https://support.10xgenomics.com/single-cell-gene-expression/datasets/

