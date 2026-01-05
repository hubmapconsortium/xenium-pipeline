cwlVersion: v1.2
class: CommandLineTool
label: Convert Alevin sparse output to anndata.AnnData object, save as h5ad
requirements:
  DockerRequirement:
    dockerPull: hubmap/spatial-transcriptomics-format-conversion:latest
  MultipleInputFeatureRequirement: {}
baseCommand: /opt/format_conversion.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  data_dir:
    type: Directory
    inputBinding:
      position: 1
outputs:
  count_matrix_h5ad:
    type: File
    outputBinding:
      glob: expr.h5ad
  sdata_zarr:
    type: Directory
    outputBinding:
      glob: Xenium.zarr