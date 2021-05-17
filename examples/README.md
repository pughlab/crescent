# Example files for CReSCENT

CReSCENT SEURAT WRAPPERS
================

Our four Seurat scripts are meant to run in the following order:
1. `Runs_Seurat_v3_MultiDatasets_QC_Normalization.R`
2. `Runs_Seurat_v3_MultiDatasets_Integration.R`
3. `Runs_Seurat_v3_MultiDatasets_PCA_Clustering_DimReduction.R`
4. `Runs_Seurat_v3_MultiDatasets_DGE.R`

They can run as one-line-commad tools and each script has its own documentation, which can be found here:
https://github.com/pughlab/crescent/tree/master/doc

`Runs_Seurat_v3_MultiDatasets_QC_Normalization.R` uses MTX files as inputs (e.g. from 10X Cell Ranger); whereas the other three scripts use R objects produced in previous steps. Users can obtain example MTX files and input parameters from here:<br />
https://zenodo.org/record/4642759/files/crescent_v2.0_pbmc_example_infiles_and_commands.tar.bz2?download=1

And example outfiles producing these files and parameters:<br />
https://zenodo.org/record/4642759/files/crescent_v2.0_pbmc_example_outfiles.tar.bz2?download=1


INTER-CONVERTING MTX AND GENE-VS-BARCODE FORMAT FILES
================

- MTX files are three files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz), e.g. from Cell Ranger's filtered_feature_bc_matrix directory.
- GENE-VS-BARCODE files are <tab> delimited files with genes in rows and barcodes in columns.
The main advantage of MTX over GENE-VS-BARCODE files, are that MTX files don't store the 0's and hence save hard drive space.


Download example MTX files from (10X):
http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz

_How to run scripts:_

To convert a gene-vs-barcode file into MTX files use: <br />
`Rscript ~/path_to/obtains_mtx_files_from_gene_x_barcode_matrix.R -i ~/path_to/gene_vs_barcode_matrix.tsv -o ~/path_to_store_outfiles/ -p prefix_for_log_files -l n` <br />
Note: the gene-vs-barcode file may be gzipped

To convert MTX files into a gene-vs-barcode file use: <br />
`Rscript ~/path_to/obtains_gene_x_barcode_matrix_from_mtx_files.R -i ~/path_to/filtered_feature_bc_matrices/ -o ~/path_to_store_outfiles/ -p prefix_for_log_files_outfiles` <br />


INFERCNV WRAPPER (detect copy number variants)
================

Example infiles and outfiles for our R wrapper to run InferCNV are provided.

The wrapper can run as a one-line-commad tool and its documentation can be found here:
https://github.com/pughlab/crescent/tree/master/doc

`Runs_InferCNV.R` uses an infile matrix with all cells (normal and cancer) in either format MTX or gene-vs-barcodes (DGE). Example infiles from glioblastoma can be downloaded from here:<br />
https://zenodo.org/record/4766424/files/crescent_infercnv_glio_example_infiles_and_commands.tar.bz2?download=1

And example outfiles producing these files and parameters:<br />
https://zenodo.org/record/4766424/files/crescent_infercnv_glio_example_outfiles_and_commands.tar.bz2?download=1

Example commands using the infiles above:<br />
`Rscript ~/r_programs/Runs_InferCNV.R -i INPUTS/glio.wGtexBrain.counts.matrix.gz -t DGE -j INPUTS/glio.wGtexBrain.sample_annots.txt -k Brain_Cerebellum,Brain_Caudate_basal_ganglia,Brain_Cortex,Brain_Nucleus_accumbens_basal_ganglia,Brain_Cerebellar_Hemisphere,Brain_Frontal_Cortex_BA9,Brain_Hippocampus -g INPUTS/gencode_v19_gene_pos.txt -m 0.1 -n 0.1 -s 0.15 -o OUTPUTS -p glio.wGtexBrain -u MAX -w 0 -a 10000`


GSVA WRAPPER (assign cell type labels to cell clusters)
================

Example infiles and outfiles for our R wrapper to run GSVA to predic cell-cluster labels are provided.

The wrapper can run as a one-line-commad tool and its documentation can be found here:
https://github.com/pughlab/crescent/tree/master/doc

`Runs_GSVA.R` uses an infile matrix with the average gene expression of cell clusters, e.g. \*AverageGeneExpression_GlobalClustering_AllDatasets\* files from `Runs_Seurat_v3_MultiDatasets_PCA_Clustering_DimReduction.R`, and a file with cell type gene sets in \*gmt format. Example infiles from PBMCs can be downloaded from here:<br />
https://zenodo.org/record/4766424/files/crescent_gsva_pbmcs_example_infiles_and_commands.tar.bz2?download=1

And example outfiles producing these files and parameters:<br />
https://zenodo.org/record/4766424/files/crescent_gsva_pbmcs_example_outfiles_and_commands.tar.bz2?download=1

Example commands using the infiles above:<br />
`Rscript ~/r_programs/Runs_GSVA.R -i INFILES/pbmcs.SEURAT_AverageGeneExpression_GlobalClustering_AllDatasets_SCT.tsv.bz2 -t DGE -c INFILES/LM22_signature.cutoff3000.gmt -o OUTFILES -p pbmcs -u MAX -w 0`


