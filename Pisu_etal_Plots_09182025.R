###RUSSELL PLOTS
plotDir <- "~/ondemand/scRNAseq/FD/scRNAseq/plots"

Russell_COLORS <- c(
  "Monocytes"          = "maroon3",
  "IM_1"               = "red",
  "IM_2"               = "burlywood4",
  "IM_3"               = "burlywood1",
  "AM_1"               = "pink4",
  "AM_2"               = "lightblue1",
  "AM_3"               = "pink2",
  "AM_4"               = "green1",
  "DC.103+11B-"        = "blue3",
  "DC"                 = "green4",
  "Neutrophils"        = "gray87",
  "Macrophages/B-cells"= "yellow1",
  "Macrophages/T-cells"= "yellow4"
)


umap_russell <- DimPlot(
  allCells,
  reduction = "umap",
  group.by = "russell_mac_celltypist_model",
  cols = Russell_COLORS
) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),       # removes box border
    axis.line = element_line(color = "black")  # keep axis lines
  ) +
  labs(title = "UMAP by Russell CellTypist Annotations")

ggsave(
  plot = umap_russell,
  filename = file.path(plotDir, "umap_russell_celltypist.pdf"),
  width = 8, height = 6
)


# Alluvial plot
alluvial_russell <- allCells@meta.data %>%
  select(ImmGenGroup, russell_mac_celltypist_model) %>%
  group_by(across(everything())) %>%
  summarize(Freq = n(), .groups = "drop") %>%
  ggplot(aes(
    axis1 = ImmGenGroup,
    axis2 = russell_mac_celltypist_model,
    y = Freq,
    fill = russell_mac_celltypist_model
  )) +
  geom_stratum(color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
  scale_x_discrete(
    limits = c("ImmGenGroup", "Russell CellTypist"),
    expand = c(.05, .05)
  ) +
  scale_fill_manual(values = Russell_COLORS) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(title = "ImmGenGroup ↔ Russell CellTypist Correspondence")

# Save as PDF
ggsave(
  plot = alluvial_russell,
  filename = file.path(plotDir, "alluvial_russell_celltypist.pdf"),
  width = 10, height = 12
)

library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)

#----------------------------------------------------------
# 1. Define your Russell label column
#----------------------------------------------------------
# Replace with the exact column name in allCells@meta.data
russellCol <- "russell_mac_celltypist_model"

# 👇 make sure RUSSELL_LABEL_COLORS exists, e.g.:
# RUSSELL_LABEL_COLORS <- c("pMO"="gold", "iMO1"="darkred", "iMO2"="red", ...)
# same mapping you used for UMAPs

#----------------------------------------------------------
# 2. Generate cluster counts per sample/condition
#----------------------------------------------------------
russellClusterCounts <- allCells@meta.data %>%
  count(sample, condition, replicate, tissue, infect, day, .data[[russellCol]]) %>%
  group_by(sample) %>%
  mutate(nPerSample = sum(n)) %>%
  ungroup() %>%
  mutate(cellsPerThousand = round((n / nPerSample) * 1000, 2))

#----------------------------------------------------------
# 3. Plot abundance with consistent colors
#----------------------------------------------------------
russellAbundancePlot <- russellClusterCounts %>%
  mutate(across(all_of(russellCol), ~ifelse(is.na(.), "unk", .))) %>%
  ggplot(aes(x = factor(replicate),
             y = cellsPerThousand,
             fill = .data[[russellCol]])) +
  geom_col() +
  geom_text(aes(label = str_trunc(.data[[russellCol]], 6, ellipsis = "..")),
            position = position_stack(vjust = 0.5), size = 2) +
  facet_grid(~ tissue + day + infect, scales = "free", space = "free") +
  scale_fill_manual(values = Russell_COLORS) +   # <- color consistency!
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = NULL, y = "Cells per thousand", title = "Russell labels")

ggsave(
  plot = russellAbundancePlot,
  filename = file.path(plotDir, "abundance_colPlot_russell.pdf"),
  width = 12, height = 6
)


# Collapse across replicates -> one pie per tissue+day+infect
russellPieByCond <- russellClusterCounts %>%
  group_by(tissue, day, infect, !!sym(russellCol)) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(tissue, day, infect) %>%
  mutate(total = sum(n),
         freq = n / total * 1000) %>%  # still per-thousand scaling if desired
  ungroup()

# Pie chart per condition
russellPiePlot <- ggplot(russellPieByCond,
       aes(x = "", y = freq, fill = .data[[russellCol]])) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  facet_grid(~ tissue + day + infect) +
  scale_fill_manual(values = Russell_COLORS) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold")
  ) +
  labs(title = "Russell labels (Pie Chart per Condition)")

ggsave(
  plot     = russellPiePlot,
  filename = file.path(plotDir, "abundance_piePlot_russell_byCondition.pdf"),
  width    = 14, height = 6
)

# Collapse across replicates -> one pie per tissue+day+infect
russellPieByCond <- russellClusterCounts %>%
  group_by(tissue, day, infect, !!sym(russellCol)) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(tissue, day, infect) %>%
  mutate(total = sum(n),
         perc = n / total * 100) %>%
  ungroup()

# Pie chart per condition, with percentages inside
russellPiePlot <- ggplot(russellPieByCond,
                         aes(x = "", y = perc, fill = .data[[russellCol]])) +
  geom_col(width = 1, color = "black") +
  geom_text(
    aes(label = paste0(round(perc, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ tissue + day + infect, ncol = 5) +
  scale_fill_manual(values = Russell_COLORS) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(size = 5, face = "bold")
  )+ guides(fill = guide_legend(ncol = 3))+
  labs(title = "Russell labels (Pie Chart per Condition)")

ggsave(
  plot     = russellPiePlot,
  filename = file.path(plotDir, "abundance_piePlot_russell_byCondition.pdf"),
  width    = 14, height = 10
)

# Ensure russellCol has fixed factor levels in the order of your palette
allCells@meta.data[[russellCol]] <- factor(
  allCells@meta.data[[russellCol]],
  levels = names(Russell_COLORS)
)

# Summarize & compute percentages + label positions
russellPieByCond <- russellClusterCounts %>%
  group_by(tissue, day, infect, !!sym(russellCol)) %>%
  summarise(Freq = sum(cellsPerThousand), .groups = "drop") %>%
  group_by(tissue, day, infect) %>%
  mutate(
    total = sum(Freq),
    perc = Freq / total * 100,
    ypos = cumsum(perc) - 0.5 * perc   # 👈 midpoint of each slice
  ) %>%
  ungroup()

# Pie chart with properly centered labels (only ≥10%)
russellPiePlot <- ggplot(russellPieByCond,
                         aes(x = "", y = perc, fill = .data[[russellCol]])) +
  geom_col(width = 1, color = "black") +
  geom_text(
    data = subset(russellPieByCond, perc >= 10),
    aes(y = ypos, 
        label = paste0(.data[[russellCol]], "\n", round(perc, 1), "%")),
    size = 3
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ tissue + day + infect, ncol = 5) +
  scale_fill_manual(values = Russell_COLORS, drop = FALSE) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    aspect.ratio = 1,
    panel.spacing = unit(2, "lines")
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(title = "Russell labels (Pie Chart per Condition, ≥10% labeled)")

ggsave(file.path(plotDir, "pie_russell_celltypist.pdf"),
       plot = russellPiePlot,
       width = 12, height = 10)

library(dplyr)
library(ggplot2)

# Summarize & compute percentages + explicit slice midpoints
russellPieByCond <- russellClusterCounts %>%
  group_by(tissue, day, infect, !!sym(russellCol)) %>%
  summarise(Freq = sum(cellsPerThousand), .groups = "drop") %>%
  group_by(tissue, day, infect) %>%
  arrange(desc(Freq)) %>%  # 👈 important: order slices
  mutate(
    total = sum(Freq),
    perc = Freq / total * 100,
    cumperc = cumsum(perc),
    ypos = cumperc - perc / 2   # 👈 midpoint of each slice
  ) %>%
  ungroup()

# Pie chart
russellPiePlot <- ggplot(russellPieByCond,
                         aes(x = "", y = perc, fill = .data[[russellCol]])) +
  geom_col(width = 1, color = "black") +
  geom_text(
    data = subset(russellPieByCond, perc >= 10),
    aes(y = ypos,
        label = paste0(.data[[russellCol]], "\n", round(perc, 1), "%")),
    size = 3
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ tissue + day + infect, ncol = 5) +
  scale_fill_manual(values = Russell_COLORS, drop = FALSE) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    aspect.ratio = 1,
    panel.spacing = unit(2, "lines")
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(title = "Russell labels (Pie Chart per Condition, ≥10% labeled)")

ggsave(file.path(plotDir, "pie_russell_celltypist.pdf"),
       plot = russellPiePlot,
       width = 12, height = 10)
