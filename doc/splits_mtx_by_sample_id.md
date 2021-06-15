Script name
================
`splits_mtx_by_sample_id.R`

Description
================
Loads an MTX set of files (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz and splits the barcodes into samples<br />
barcodes.tsv.gz file must contain barcodes in format like:<br />
SampleId1_ATCGATCGATCGATCG<br />
SampleId1_ATTTCTCGATCGATTT<br />
SampleId2_AGGGGTCGATCGACCC<br />
SampleId2_AAAAATCGATCGAGGG<br />
It will use the string before the last `_` as sample ID

It allows the user to provide parameters as one-line commands via R library(optparse).

Outfiles
================
One new set of MTX files for each sample<br />

General workflow
================
  1. Loads original MTX set of files
  2. Identified sample IDs
  3. Generates new MTX sets of files, one for each sample
  4. Ends

Example commands
================
To display help commands type: <br />
`Rscript /path_to/splits_mtx_by_sample_id.R -h`

To run the script type something like:<br />
`Rscript /path_to/splits_mtx_by_sample_id.R -i /path_to/original_mtx_directory -o /path_to/new_split_mtx -p string_prefix_for_log_files`

Inputs Description
================
a) A set of MTX files (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz)<br />

Outputs Description
================
|               Extension                  |                        Contents                        |
| --------------------------------------   |  ----------------------------------------------------- |
| New barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz | One folder with MTX files for each sample |
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
DropletUtils_1.8.0 Seurat_3.2.1 optparse_1.6.6
