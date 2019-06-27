cwlVersion: v1.0

class: CommandLineTool

baseCommand: [Rscript]

inputs:
  R_script:
    type: File
    inputBinding:
      position: 0

  sc_input:
    type: Directory
    inputBinding:
      position: 1
      prefix: -i

  sc_input_type:
    type: string
    inputBinding:
      position: 2
      prefix: -t

  resolution:
    type: int?
    inputBinding:
      position: 3
      prefix: -r

  project_id:
    type: string
    inputBinding:
      position: 5
      prefix: -p

  summary_plots:
    type: string?
    inputBinding:
      position: 6
      prefix: -s

  colour_tsne_discrete:
    type: File?
    inputBinding:
      position: 7
      prefix: -c

  list_genes:
    type: string?
    inputBinding:
      position: 8
      prefix: -g

  opacity:
    type: float?
    inputBinding:
      position: 9
      prefix: -a

  pca_dimensions:
    type: int?
    inputBinding:
      position: 10
      prefix: -d

  percent_mito:
    type: string?
    inputBinding:
      position: 11
      prefix: -m

  number_genes:
    type: string?
    inputBinding:
      position: 12
      prefix: -n

  return_threshold:
    type: float?
    inputBinding:
      position: 13
      prefix: -e

  number_cores:
    type: string?
    inputBinding:
      position: 14
      prefix: -u

outputs:
  seurat_output:
    type: Directory
    outputBinding:
      glob: SEURAT
