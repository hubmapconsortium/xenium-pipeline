#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.2
label: spatial transcriptomics pipeline including analysis with scanpy and squidpy
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
inputs:
  data_dir:
    label: "Directory containing FASTQ files"
    type: Directory
  assay:
    label: "scRNA-seq assay"
    type: string
outputs:
  count_matrix_h5ad:
    outputSource: convert_formats/count_matrix_h5ad
    type: File
    label: "Count matrix converted to h5ad"
#  scanpy_qc_results:
#    outputSource: compute_qc_results/scanpy_qc_results
#    type: File
#    label: "Quality control metrics from Scanpy"
  dispersion_plot:
    outputSource: scanpy_analysis/dispersion_plot
    type: File
    label: "Gene expression dispersion plot"
  umap_plot:
    outputSource: scanpy_analysis/umap_plot
    type: File
    label: "UMAP dimensionality reduction plot"
  umap_density_plot:
    outputSource: scanpy_analysis/umap_density_plot
    type: File
    label: "UMAP dimensionality reduction plot, colored by cell density"
  spatial_plot:
    outputSource: scanpy_analysis/spatial_plot
    type: File?
    label: "Slide-seq bead plot, colored by Leiden cluster"
  filtered_data_h5ad:
    outputSource: scanpy_analysis/filtered_data_h5ad
    type: File
    label: Full data set of filtered results
    doc: >-
      Full data set of filtered results: expression matrix, coordinates in
      dimensionality-reduced space (PCA and UMAP), cluster assignments via
      the Leiden algorithm, and marker genes for one cluster vs. rest
  marker_gene_plot_t_test:
    outputSource: scanpy_analysis/marker_gene_plot_t_test
    type: File
    label: "Cluster marker genes, t-test"
  marker_gene_plot_logreg:
    outputSource: scanpy_analysis/marker_gene_plot_logreg
    type: File
    label: "Cluster marker genes, logreg method"
  scvelo_annotated_h5ad:
    outputSource: scvelo_analysis/annotated_h5ad_file
    type: File?
    label: "scVelo-annotated h5ad file, including cell RNA velocity"
  scvelo_embedding_grid_plot:
    outputSource: scvelo_analysis/embedding_grid_plot
    type: File?
    label: "scVelo velocity embedding grid plot"
  squidpy_annotated_h5ad:
    outputSource: squidpy_analysis/squidpy_annotated_h5ad
    type: File?
  neighborhood_enrichment_plot:
    outputSource: squidpy_analysis/neighborhood_enrichment_plot
    type: File?
  co_occurrence_plot:
    outputSource: squidpy_analysis/co_occurrence_plot
    type: File?
  interaction_matrix_plot:
    outputSource: squidpy_analysis/interaction_matrix_plot
    type: File?
  centrality_scores_plot:
    outputSource: squidpy_analysis/centrality_scores_plot
    type: File?
  ripley_plot:
    outputSource: squidpy_analysis/ripley_plot
    type: File?
  squidpy_spatial_plot:
    outputSource: squidpy_analysis/spatial_plot
    type: File?
steps:
  convert_formats:
    in:
      data_dir:
        source: data_dir
      assay:
        source: assay
    out:
      - count_matrix_h5ad
    run: steps/convert-formats.cwl
  scanpy_analysis:
    in:
      assay:
        source: assay
      h5ad_file:
        source: convert_formats/count_matrix_h5ad
    out:
      - filtered_data_h5ad
      - umap_plot
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
      - dispersion_plot
      - umap_density_plot
      - spatial_plot
    run: salmon-rnaseq/steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
  squidpy_analysis:
    in:
      assay:
        source: assay
      h5ad_file:
        source: scanpy_analysis/filtered_data_h5ad
      img_dir:
        source: img_dir
    out:
      - squidpy_annotated_h5ad
      - neighborhood_enrichment_plot
      - co_occurrence_plot
      - interaction_matrix_plot
      - ripley_plot
      - centrality_scores_plot
      - spatial_plot
    run: salmon-rnaseq/steps/squidpy-analysis.cwl
    label: "Spatial analysis via SquidPy"
#  compute_qc_results:
#    in:
#      assay:
#        source: assay
#      primary_matrix_path:
#        source: salmon_quantification/count_matrix_h5ad
#      secondary_matrix_path:
#        source: scanpy_analysis/filtered_data_h5ad
#      salmon_dir:
#        source: salmon_quantification/salmon_output
#    out:
#      - scanpy_qc_results
#      - qc_metrics
#    run: steps/compute-qc-metrics.cwl
#    label: "Compute QC metrics"
