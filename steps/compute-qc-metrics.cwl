cwlVersion: v1.2
class: CommandLineTool
label: Compute QC metrics
requirements:
  DockerRequirement:
    dockerPull: hubmap/spatial-transcriptomics-analysis:latest
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
  salmon_dir:
    type: Directory
    inputBinding:
      position: 3
outputs:
  scanpy_qc_results:
    type: File
    outputBinding:
      glob: qc_results.hdf5
