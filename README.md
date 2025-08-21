# -------------------------------------------------------------
# Mouse Mtb Innate Atlas: DEG and GO:BP GSEA Viewer
# -------------------------------------------------------------
# This R script serves as documentation for the Shiny application.
# It describes installation, dependencies, usage, acknowledgements,
# and references in plain R comment format.

# -------------------------------------------------------------
# Overview
# -------------------------------------------------------------
# This repository contains a Shiny application for interactive
# visualization of Differential Gene Expression (DEG) and
# Gene Set Enrichment Analysis (GSEA) results from mouse models
# of Mycobacterium tuberculosis (Mtb) infection.
#
# The app integrates:
#  - DEG results (Seurat::FindMarkers)
#  - GSEA (GO Biological Processes) results
#  - Redundancy-reduced GSEA results (rrvgo)
#
# Features:
#  - Interactive filtering across tissue, infection status,
#    timepoint, and cell type
#  - Visualization options: Volcano, Dot, Bar, Heatmap
#  - Collapsible acknowledgements/references block
#  - Styled tabs for Plot and Table (lavender background)
#  - Downloadable plots and data

# -------------------------------------------------------------
# Installation
# -------------------------------------------------------------
# Clone the repository:
#   git clone https://github.com/eshamskh/Mouse-Mtb-Innate-Atlas.git
#   cd Mouse-Mtb-Innate-Atlas.git
#
# Install required R packages:
#   install.packages(c(
#     "shiny", "ggplot2", "dplyr", "readr", "DT", "tibble",
#     "tidyr", "purrr", "scales", "plotly", "openxlsx", "stringr"
#   ))
#
# If you are running the full DEG/GSEA pipeline, also install:
#   install.packages(c("Seurat", "fgsea", "rrvgo", "clusterProfiler", "org.Mm.eg.db"))

# -------------------------------------------------------------
# Running the App
# -------------------------------------------------------------
# From within R:
#   shiny::runApp()
#
# Or deploy via Posit Connect Cloud:
# https://eshamskh-mouse-mtb-innate-atlas.share.connect.posit.cloud


# -------------------------------------------------------------
# Data
# -------------------------------------------------------------
# The app expects merged results CSVs in the following directories:
#   merged_DEG_results/
#   merged_GSEA_results/
#   merged_GSEA_reduced_results/
#
# Each folder should contain comparison-specific CSVs
# (e.g., merged_DEG_0_15.csv).

# -------------------------------------------------------------
# Acknowledgements & References
# -------------------------------------------------------------
# Made in collaboration between the Gerner Lab
# (University of Washington, Department of Immunology)
# and the Urdahl Lab (Seattle Children’s Research Institute, CGIDR).
#
# For use in submissions and publications, and for further information
# regarding methods and analysis, please reference:
#
#   "Monocytic Niches Enable Mycobacterium tuberculosis Persistence in Lymph Nodes"
#   (Manuscript in Progress)
#   Authors: Elya A. Shamskhou, Fergal R. Duffy, Lauren M. Cross,
#            Vitaly V. Ganusov, Courtney R. Plumlee, Benjamin H. Gern,
#            Alan H. Dierks, Sara B. Cohen, Kevin B. Urdahl*,
#            Michael Y. Gerner*

# -------------------------------------------------------------
# License
# -------------------------------------------------------------
# This repository is provided for research and academic use.
