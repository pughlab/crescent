cwlVersion: v1.0

class: CommandLineTool

baseCommand: [cellranger, count]

inputs:
  run_name:
    type: string
    inputBinding:
      position: 1
      prefix: --id=
      separate: false

  transcriptome:
    type: Directory
    inputBinding:
      position: 2
      prefix: --transcriptome=
      separate: false

  fastqs:
    type: Directory
    inputBinding:
      position: 3
      prefix: --fastqs=
      separate: false

  sample:
    type: string
    inputBinding:
      position: 4
      prefix: --sample=
      separate: false

  jobmode:
    type: string?
    inputBinding:
      position: 5
      prefix: --jobmode=
      separate: false

outputs:
  cellranger_meta_output:
    type: Directory
    outputBinding:
      glob: $(inputs.run_name)
  cellranger_output:
    type: Directory
    outputBinding:
      #glob: $(inputs.run_name)
      glob: $(inputs.run_name)/outs/filtered_gene_bc_matrices/*
