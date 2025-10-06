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

write_csv(deg_results, sprintf("analysis/time_hallmark_deg_results_%s.csv", analysisMarkerCol))

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
saveRDS(time_hallmark_deg_results_ImmGenGroup, 
        file = "analysis/time_hallmark_deg_results_ImmGenGroup_renamed.rds")

# 2) Mutate celltype: transDC → tDC
library(dplyr)
time_hallmark_gsea_ImmGenGroup <- time_hallmark_gsea_ImmGenGroup %>%
  mutate(celltype = if_else(celltype == "transDC", "tDC", celltype))

# 3) (Optional) Save the updated object
saveRDS(time_hallmark_gsea_ImmGenGroup, 
        file = "analysis/time_hallmark_gsea_ImmGenGroup_renamed.rds")



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

###################
# Your desired Hallmark pathway order (without the "HALLMARK_" prefix)
my_pathways <- c(
  "ANGIOGENESIS",
  "MYOGENESIS",
  "EPITHELIAL_MESENCHYMAL_TRANSITION",
  "COAGULATION",
  "MTORC1_SIGNALING",
  "GLYCOLYSIS",
  "FATTY_ACID_METABOLISM",
  "CHOLESTEROL_HOMEOSTASIS",
  "OXIDATIVE_PHOSPHORYLATION",
  "XENOBIOTIC_METABOLISM",
  "HYPOXIA",
  "HEME_METABOLISM",
  "APOPTOSIS",
  "P53_PATHWAY",
  "IL2_STAT5_SIGNALING",
  "TNFA_SIGNALING_VIA_NFKB",
  "COMPLEMENT",
  "ALLOGRAFT_REJECTION",
  "INFLAMMATORY_RESPONSE",
  "IL6_JAK_STAT3_SIGNALING",
  "INTERFERON_GAMMA_RESPONSE",
  "INTERFERON_ALPHA_RESPONSE"
)

# Your desired facet/block order:
# blood Mtb- D15 D28,
# lung  Mtb- D15 D28 THEN Mtb+ D15 D28,
# mln   Mtb- D15 D28 THEN Mtb+ D15 D28
my_blocks <- c(
  "blood Mtb- D15", "blood Mtb- D28",
  "lung Mtb- D15",  "lung Mtb- D28",
  "lung Mtb+ D15",  "lung Mtb+ D28",
  "mln Mtb- D15",   "mln Mtb- D28",
  "mln Mtb+ D15",   "mln Mtb+ D28"
)

# Filter your GSEA results to Hallmark only, and drop "unk" cell types & AM in blood/mln if desired
gsea_filter <- time_hallmark_gsea_ImmGenGroup_renamed %>%
  dplyr::filter(
    grepl("^HALLMARK", pathway),
    !celltype %in% "unk"
    # , !(celltype == "AM" & tissue %in% c("blood", "mln")) # <- keep or drop; you had this version in one script
  )

hallmarkPlot <- clustDotplotGseaHeatmap(
  gseaRes           = gsea_filter,
  allCellTypes      = valid_celltypes,
  mtbPosCellTypes   = valid_celltypes,
  fdrThresh         = 0.05,
  trimPathNames     = TRUE,
  pointSizeRange    = c(2,5),
  pathway_order     = my_pathways,
  facet_block_order = my_blocks,
  bold_axes         = TRUE,
  y_expand_mult     = 0.03   # tight top/bottom padding
)

rows <- attr(hallmarkPlot, "n_pathways")
# Tweak inches_per_row to taste (smaller = tighter). 0.15–0.18 works well.
inches_per_row <- 0.22
height_in <- max(4, min(12, rows * inches_per_row))

hallmarkPlot <- hallmarkPlot +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey60", linewidth = 0.1),
    panel.grid.major.x = element_line(color = "grey60", linewidth = 0.1)
  )

ggsave(
  plot     = hallmarkPlot,
  filename = sprintf("analysis/no5_time_hallmark_gsea_heatmap_%s_customOrder.pdf", "ImmGenGroup"),
  width    = 20,
  height   = height_in
)

# If you haven't loaded it:
library(grDevices)  # cairo_pdf lives here

# Optional: pick a common font to avoid substitution
hallmarkPlot <- hallmarkPlot + theme(text = element_text(family = "Helvetica"))

# optional: check if Cairo is available
if (capabilities("cairo")) {
  ggsave(
    "analysis/hallmark_heatmap_vector.pdf",
    plot   = hallmarkPlot,
    width  = 20, height = height_in, units = "in",
    device = grDevices::cairo_pdf,
    bg     = "white",
  )
} else {
  # fallback to base PDF (still vector)
  ggsave(
    "analysis/hallmark_heatmap_vector.pdf",
    plot   = hallmarkPlot,
    width  = 20, height = height_in, units = "in",
    device = grDevices::pdf,
    bg     = "white",
    useDingbats = FALSE
  )
}



