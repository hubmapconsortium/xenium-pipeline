cwlVersion: v1.2
class: CommandLineTool
label: Compute QC metrics
requirements:
  DockerRequirement:
    dockerPull: hubmap/spatial-transcriptomics-analysis:0.1.3
baseCommand: /opt/compute_qc_metrics.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  primary_matrix_path:
    type: File
    inputBinding:
      position: 1
outputs:
  scanpy_qc_results:
    type: File
    outputBinding:
      glob: qc_results.hdf5
