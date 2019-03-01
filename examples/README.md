# Example files for CReSCENT

CELL RANGER
================

We are using Cell Ranger v3.X

Documentation of cellranger can be found here:
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

Inputs: <br />
Takes fastq files, each lane has three fastq files (R1, R2 and I1), which can be Gb in size. <br />
You can download example files from 10X: <br />
https://support.10xgenomics.com/single-cell-gene-expression/datasets

Outputs:
Include three files with the sparse matrix in Matrix Market (*mtx) format with the mapped read counts per gene, per cell barcode: <br />
`outs/filtered_gene_bc_matrices/Reference_Transcriptome/` <br />
barcodes.tsv <br />
genes.tsv (features.tsv in Cell Ranger v3) <br />
matrix.mtx <br />

SEURAT
================

Script name: Runs_Seurat_Clustering.R

Infiles: <br />
Either the `outs/filtered_gene_bc_matrices/Reference_Transcriptome/` from Cell Ranger or a matrix with genes (rows) and cell barcodes (columns)

Outfiles: <br />
See `crescent/doc/Runs_Seurat_Clustering.md`
