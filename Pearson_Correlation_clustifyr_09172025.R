suppressPackageStartupMessages({
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(dplyr)
})

DefaultAssay(allCells) <- "SCT"

# --- helpers ---
get_avg <- function(obj, group.by, slot = "scale.data") {
  # use SCT z-scores/residuals to remove global level effects
  Seurat::AverageExpression(obj, assays = "SCT", slot = slot, group.by = group.by)$SCT
}

# Correlate *deviation from pooled mean* (per gene) to kill the "all-positive" floor
center_cor <- function(A, B, genes, method = "pearson") {
  X <- as.matrix(A[genes, , drop = FALSE])  # genes x labels_A
  Y <- as.matrix(B[genes, , drop = FALSE])  # genes x labels_B
  mu <- rowMeans(cbind(X, Y))               # pooled gene-wise mean across both sets
  Xd <- sweep(X, 1, mu, "-")
  Yd <- sweep(Y, 1, mu, "-")
  cor(Xd, Yd, method = method, use = "pairwise.complete.obs")
}

# --- build averaged matrices from SCT scale.data ---
avg_IG   <- get_avg(allCells, "ImmGenGroup",             slot = "scale.data")
avg_fine <- get_avg(allCells, "ImmgenClustersFine",      slot = "scale.data")
avg_ct   <- get_avg(allCells, "russell_mac_celltypist_model", slot = "scale.data")

# --- gene universe: intersect, drop housekeeping, keep top variable genes ---
g <- Reduce(intersect, list(rownames(avg_IG), rownames(avg_fine), rownames(avg_ct)))

# Drop common low-informative genes (mt-, ribosomal, hemoglobin, etc.)
drop_re <- "^(mt-|MT-|Rpl|Rps|RPL|RPS|Hb-|Hba|Hbb|Malat1|Mrpl|Mrps)"
g <- g[!grepl(drop_re, g)]

# Keep top-N most variable genes across the pooled columns
pool <- cbind(avg_IG[g, , drop = FALSE], avg_fine[g, , drop = FALSE], avg_ct[g, , drop = FALSE])
v <- matrixStats::rowVars(as.matrix(pool))
n_top <- min(3000, length(v))  # tune this (e.g., 2000–5000)
keep <- names(sort(v, decreasing = TRUE))[seq_len(n_top)]

# --- build centered correlation matrices (IG vs Fine, IG vs Celltypist) ---
cmA <- center_cor(avg_IG, avg_fine, keep, method = "pearson")
cmB <- center_cor(avg_IG, avg_ct,   keep, method = "pearson")

# harmonize row names to your canonical labels
rownames(cmA) <- dplyr::recode(rownames(cmA), "transDC" = "tDC", "unk" = "Unk")
rownames(cmB) <- dplyr::recode(rownames(cmB), "transDC" = "tDC", "unk" = "Unk")
# --- anchor columns to your ImmGenGroup order and sort by that anchor ---

# Helper: order columns by the ImmGenGroup row they correlate to best
order_cols_by_anchor <- function(cor_mat, row_order) {
  have_rows <- intersect(row_order, rownames(cor_mat))
  M <- cor_mat[have_rows, , drop = FALSE]
  # best matching row per column
  best_row <- apply(M, 2, function(x) names(which.max(x)))
  # map row name -> its position in your canonical order
  row_pos <- setNames(seq_along(row_order), row_order)
  anchor_idx <- row_pos[best_row]
  # tie-breaker within the same anchor: stronger max corr first
  best_val <- mapply(function(col, br) if (is.na(br)) -Inf else M[br, col],
                     colnames(M), best_row, USE.NAMES = FALSE)
  # any NA anchors go to the end; keep column name as final tiebreaker
  na_anchor <- is.na(anchor_idx)
  ord <- order(na_anchor, anchor_idx, -best_val, colnames(M))
  colnames(M)[ord]
}

# (Optional) get the name of the anchor row for a top annotation
anchor_group <- function(cor_mat, row_order) {
  have_rows <- intersect(row_order, rownames(cor_mat))
  M <- cor_mat[have_rows, , drop = FALSE]
  apply(M, 2, function(x) names(which.max(x)))
}

# Build complete row orders (your preferred first, then any extras)
row_orderA <- c(preferred[preferred %in% rownames(cmA)], setdiff(rownames(cmA), preferred))
row_orderB <- c(preferred[preferred %in% rownames(cmB)], setdiff(rownames(cmB), preferred))

# Compute the anchored column orders
col_orderA <- order_cols_by_anchor(cmA, preferred)
col_orderB <- order_cols_by_anchor(cmB, preferred)

# (Optional) column anchor labels for a top annotation bar
col_anchorA <- anchor_group(cmA, preferred)
col_anchorB <- anchor_group(cmB, preferred)

# 1) Your new canonical order
preferred <- c("pMO","iMO1","iMO2","iMO3","MC","AM",
               "migDC1","migDC2","tDC","resDC1","resDC2","pDC","Basophil","Neutrophil","Unk")

# 2) Recompute row + column orders anchored to that sequence
row_orderA <- c(preferred[preferred %in% rownames(cmA)], setdiff(rownames(cmA), preferred))
row_orderB <- c(preferred[preferred %in% rownames(cmB)], setdiff(rownames(cmB), preferred))

col_orderA <- order_cols_by_anchor(cmA, preferred)
col_orderB <- order_cols_by_anchor(cmB, preferred)

# (optional) anchor labels for top annotation AND column blocks
col_anchorA <- anchor_group(cmA, preferred)
col_anchorB <- anchor_group(cmB, preferred)

# Color scale tuned to your actual correlation spread (quantile-based for contrast)
qA <- quantile(cmA, c(0.05, 0.50, 0.95), na.rm = TRUE)
qB <- quantile(cmB, c(0.05, 0.50, 0.95), na.rm = TRUE)
col_funA <- circlize::colorRamp2(qA, c("#2166AC", "#FFFFFF", "#B2182B"))
col_funB <- circlize::colorRamp2(qB, c("#2166AC", "#FFFFFF", "#B2182B"))
library(ComplexHeatmap)
library(grid)

# Unicode: ρ = \u03C1, Δ = \u0394
legend_title <- "Pearson r (centered delta)"
ff <- "sans"  # safe generic
legend_title <- "Pearson r (centered delta)"

haA <- Heatmap(
  cmA,
  name = NULL,
  heatmap_legend_param = list(
    title = legend_title,
    title_gp  = gpar(fontfamily = ff),
    labels_gp = gpar(fontfamily = ff)
  ),
  col = col_funA,
  cluster_rows = FALSE,
  cluster_columns = FALSE,        # <-- we control columns now
  row_order = row_orderA,
  column_order = col_orderA,
  top_annotation = HeatmapAnnotation(
    anchor = col_anchorA[col_orderA],
    annotation_legend_param = list(title = "Anchors to IG")
  ),
  column_title = "ImmGenGroup vs ImmgenClustersFine"
)

haB <- Heatmap(
  cmB,
  name = NULL,
  heatmap_legend_param = list(
    title = legend_title,
    title_gp  = gpar(fontfamily = ff),
    labels_gp = gpar(fontfamily = ff)
  ),
  col = col_funB,
  cluster_rows = FALSE,
  cluster_columns = FALSE,        # <-- we control columns now
  row_order = row_orderB,
  column_order = col_orderB,
  top_annotation = HeatmapAnnotation(
    anchor = col_anchorB[col_orderB],
    annotation_legend_param = list(title = "Anchors to IG")
  ),
  column_title = "ImmGenGroup vs Celltypist"
)

# PDF via Cairo (good Unicode support)
if (capabilities("cairo")) {
  pdf("analysis/ImmGenGroup_vs_ImmgenClustersFine_cor_SCT_centeredpearson.pdf",
      width = 10, height = 8, family = ff, useDingbats = FALSE)
  draw(haA); dev.off()
  
  pdf("analysis/ImmGenGroup_vs_Celltypist_cor_SCT_centeredpearson.pdf",
      width = 10, height = 8, family = ff, useDingbats = FALSE)
  draw(haB); dev.off()
} else {
  # SVG fallback (also great with Unicode)
  if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
  svglite::svglite("analysis/ImmGenGroup_vs_ImmgenClustersFine_cor_SCT_centeredpearson.svg",
                   width = 10, height = 8, system_fonts = list(sans = ff))
  draw(haA); dev.off()
  
  svglite::svglite("analysis/ImmGenGroup_vs_Celltypist_cor_SCT_centeredpearson.svg",
                   width = 10, height = 8, system_fonts = list(sans = ff))
  draw(haB); dev.off()
}


# Also export matrices
write.csv(cmA, "analysis/cor_ImmGenGroup_vs_ImmgenClustersFine_SCT_centeredpearson.csv")
write.csv(cmB, "analysis/cor_ImmGenGroup_vs_Celltypist_SCT_centeredpearson.csv")
