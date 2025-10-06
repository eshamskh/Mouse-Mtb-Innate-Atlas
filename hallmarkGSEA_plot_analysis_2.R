# Hallmark GSEA Bubble Plot: D15/D28 vs D0 using ImmGenGroup ####

# Load required libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(fgsea)
library(msigdbr)
library(conflicted)

# Prefer base intersect to avoid ambiguity
conflicts_prefer(base::intersect)

# Load source functions and data
source("deg_functions.R")

# Define analysis column
analysisMarkerCol <- "ImmGenGroup"

# Add ImmGenGroup metadata
allCells <- addGroupings(allCells, groupingDefs = IMMGEN_GROUPINGS_3, groupLabel = analysisMarkerCol)

# Define tissues and celltypes shared at Day 0
shared_tissues <- intersect(unique(allCells$tissue[allCells$day == 0]), unique(allCells$tissue))
shared_celltypes <- intersect(
  unique(allCells@meta.data[[analysisMarkerCol]][allCells@meta.data$day == 0]),
  unique(allCells@meta.data[[analysisMarkerCol]])
)

# DEG comparisons: D15/D28 Mtb+ or Mtb- vs D0 Mtb- baseline
mtb_groups <- c("Mtb+", "Mtb-")
deg_comparisons <- expand_grid(
  test_day = c(15, 28),
  tissue = shared_tissues,
  celltype = shared_celltypes,
  Mtb = mtb_groups
)

# Run DEG analysis comparing test group vs matched D0 Mtb- cells
deg_results <- pmap_dfr(deg_comparisons, function(test_day, tissue, celltype, Mtb) {
  test_cells <- colnames(allCells)[
    allCells$day == test_day &
      allCells$Mtb == Mtb &
      allCells$tissue == tissue &
      allCells[[analysisMarkerCol]] == celltype
  ]
  ref_cells <- colnames(allCells)[
    allCells$day == 0 &
      allCells$Mtb == "Mtb-" &
      allCells$tissue == tissue &
      allCells[[analysisMarkerCol]] == celltype
  ]
  
  if (length(test_cells) < 10 || length(ref_cells) < 10) return(NULL)
  
  combined <- subset(allCells, cells = c(test_cells, ref_cells))
  combined$comparison_group <- factor(ifelse(Cells(combined) %in% test_cells, test_day, 0), levels = c(0, test_day))
  Idents(combined) <- combined$comparison_group
  
  message(sprintf("DEGs: %s D%d (%s) vs D0 (Mtb-) in %s", celltype, test_day, Mtb, tissue))
  
  FindMarkers(combined,
              ident.1 = test_day,
              ident.2 = 0,
              assay = "SCT",
              test.use = "negbinom",
              min.pct = 0.2,
              logfc.threshold = 0.5,
              recorrect_umi = FALSE) %>%
    rownames_to_column("genes") %>%
    mutate(tissue = tissue,
           celltype = celltype,
           comparison = sprintf("D%d vs D0, %s", test_day, Mtb))
})

write_csv(time_hallmark_gsea_ImmGenGroup_renamed, sprintf("analysis/time_hallmark_deg_results_%s.csv", analysisMarkerCol))

# Run GSEA on Hallmark gene sets
hallmarkGsea <- degGseaEnrichment(
  degResults = deg_results,
  groupingFactors = c("tissue", "celltype", "comparison"),
  saveFile = sprintf("analysis/time_hallmark_gsea_%s.rds", analysisMarkerCol)
)


# recode_transDC_to_tDC.R
# 1) Load your GSEA results
time_hallmark_gsea_ImmGenGroup <- readRDS("path/to/time_hallmark_gsea_ImmGenGroup.rds")

# 2) Mutate celltype: transDC → tDC
library(dplyr)
time_hallmark_deg_results_ImmGenGroup <- time_hallmark_deg_results_ImmGenGroup %>%
  mutate(celltype = if_else(celltype == "transDC", "tDC", celltype))

# 3) (Optional) Save the updated object
saveRDS(time_hallmark_gsea_ImmGenGroup, 
        file = "analysis/time_hallmark_gsea_results_ImmGenGroup.rds")

# 2) Mutate celltype: transDC → tDC
library(dplyr)
time_hallmark_gsea_ImmGenGroup <- time_hallmark_gsea_ImmGenGroup %>%
  mutate(celltype = if_else(celltype == "transDC", "tDC", celltype))

# 3) (Optional) Save the updated object
saveRDS(time_hallmark_gsea_ImmGenGroup, 
        file = "analysis/time_hallmark_gsea_ImmGenGroup_renamed.rds")

write_csv(time_hallmark_deg_results_ImmGenGroup_renamed, sprintf("analysis/time_hallmark_deg_results_%s.csv", analysisMarkerCol))


# Define cell type order and filter
ordered_celltypes <- c("Neutrophil", "pMO", "iMO1", "iMO2", "iMO3", "MC", "AM", 
                       "migDC1", "migDC2", "tDC", "resDC1", "resDC2", "Basophil")
valid_celltypes <- intersect(ordered_celltypes, unique(time_hallmark_deg_results_ImmGenGroup_renamed$celltype))

analysisMarkerCol <- "ImmGenGroup"

# Generate Hallmark bubble heatmap using updated cell type order
hallmarkPlot <- clustDotplotGseaHeatmap(
  gseaRes = filter(time_hallmark_gsea_ImmGenGroup_renamed, grepl("^HALLMARK", pathway) & !celltype %in% c("unk")),
  allCellTypes = valid_celltypes,
  mtbPosCellTypes = valid_celltypes,
  fdrThresh = 0.05
)

# Save final plot
ggsave(plot = hallmarkPlot,
       filename = sprintf("analysis/no2time_hallmark_gsea_heatmap_%s.pdf", analysisMarkerCol),
       width = 20, height = 10)

#!/usr/bin/env Rscript
# plot_hallmark_heatmap_exclude_AM_blood_mln.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# 1) cell‐type ordering + intersection
ordered_celltypes <- c(
  "Neutrophil", "pMO", "iMO1", "iMO2", "iMO3", "MC", "AM", 
  "migDC1", "migDC2", "tDC", "resDC1", "resDC2", "Basophil"
)

# assume you’ve already loaded your GSEA results and DEG results:
# time_hallmark_gsea_ImmGenGroup <- readRDS("time_hallmark_gsea_ImmGenGroup_mod.rds")
# time_hallmark_deg_results_ImmGenGroup <- readRDS("…")

valid_celltypes <- intersect(
  ordered_celltypes,
  unique(time_hallmark_gsea_ImmGenGroup_renamed$celltype)
)

analysisMarkerCol <- "ImmGenGroup"

# 2) filter out:
#    • only HALLMARK pathways
#    • drop any “unk” celltypes
#    • **add**: drop AM in blood & mln
gsea_filter <- time_hallmark_gsea_ImmGenGroup_renamed %>%
  filter(
    grepl("^HALLMARK", pathway),
    !celltype %in% "unk",
    !(celltype == "AM" & tissue %in% c("blood", "mln"))
  )

# 3) generate the heatmap with your existing function
hallmarkPlot <- clustDotplotGseaHeatmap(
  gseaRes           = gsea_filter,
  allCellTypes      = valid_celltypes,
  mtbPosCellTypes   = valid_celltypes,
  fdrThresh         = 0.05
)

# 4) save
ggsave(
  plot    = hallmarkPlot,
  filename= sprintf("analysis/no3time_hallmark_gsea_heatmap_%s.pdf", analysisMarkerCol),
  width   = 20,
  height  = 10
)

