library(tidyverse)
library(ggplot2)



#D15 and D28 iMO2 and MC Mtb pos vs D0

# Genes of interest (same as before)
genes_of_interest <- c(
  "Irf1",
  "Il6",
  "Il1rn",
  "Il1b",
  "Il12b",
  "Mmp2",
  "H2-D1",
  "H2-T23",
  "H2-K1",
  "Gbp7",
  "Stat1",
  "Nos2",
  "Fbxl5",
  "Hmox1",
  "Ctsz",
  "Ctss",
  "Ctsc",
  "Slc7a11",
  "Vegfa",
  "Lyz2",
  "Lamp1",
  "Hif1a",
  "Cd68",
  "Bcl2a1b",
  "Vcam1",
  "Il18bp",
  "H2-Q4",
  "Cd274",
  "Inhba",
  "Il1a",
  "Ly6i",
  "H2-M2",
  "Ccl5",
  "Apoe",
  "Acod1",
  "Isg15",
  "Gbp5",
  "Gbp2",
  "Cxcl9",
  "Cxcl10",
  "Tnf",
  "Procr",
  "Mertk",
  "Sash1",
  "Ctsd",
  "Zbp1",
  "Socs1",
  "Irf7",
  "Oasl2",
  "Oasl1",
  "Il15ra",
  "Gm4951",
  "Cd40"
)


# Load merged time_vs0 DEGs
deg_time <- read_csv("merged_DEG_results/merged_DEG_time.csv", show_col_types = FALSE)

# Pre-filter for Mtb+ vs D0 in MC & iMO2, our genes
df_base <- deg_time %>%
  filter(
    comparison_family == "time_vs0",
    group1       %in% c("15","28"),
    Mtb          == "pos",
    celltype     %in% c("MC","iMO2"),
    gene         %in% genes_of_interest
  ) %>%
  mutate(
    day    = factor(group1, levels = c("15","28")),
    size   = -log10(p_val_adj),
    log2FC = avg_log2FC,
    gene   = factor(gene, levels = genes_of_interest)
  )

for (tiss in c("lung","mln")) {
  df <- df_base %>% filter(tissue == tiss)
  if (!nrow(df)) next
  
  # dynamic height (~0.3" per gene; min 6")
  h <- max(6, 0.3 * n_distinct(df$gene))
  
  p <- ggplot(df, aes(x = celltype, y = gene)) +
    geom_point(aes(size = size, color = log2FC), alpha = .8) +
    facet_grid(. ~ day, scales = "free_y", space = "free_y") +
    scale_color_gradient2(
      low      = "blue", mid = "white", high = "red", midpoint = 0,
      name     = "log₂FC"
    ) +
    scale_size_continuous(range = c(2,8), name = "-log₁₀(adj p)") +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill="grey20", color=NA),
      plot.background  = element_rect(fill="white", color=NA),
      panel.grid.major = element_line(color="grey30"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill="grey30", color=NA),
      strip.text       = element_text(color="white", face="bold"),
      axis.text.x      = element_text(color="black", angle=45, hjust=1),
      axis.text.y      = element_text(color="black"),        # <-- genes in black
      axis.title       = element_text(color="black"),
      plot.margin      = margin(5,5,5,5)
    ) +
    labs(
      title = sprintf("Mtb⁺ vs D0 (MC & iMO2) — %s", tiss),
      x     = "Cell Type",
      y     = "Gene"
    )
  
  
  out_file <- sprintf("plots/2bubble_Mtbpos_vsD0_MC_iMO2_%s.pdf", tiss)
  ggsave(out_file, p, width = 12, height = h, limitsize = FALSE)
  message("Saved: ", out_file)
}


###15 vs 28 analysis

# — adjust these to your own paths if needed —
merged_time_csv <- "merged_DEG_results/merged_DEG_time.csv"
out_dir        <- "plots"
dir.create(out_dir, showWarnings = FALSE)




# load your merged time‐DEG table
deg_time <- read_csv(merged_time_csv, show_col_types = FALSE)

# filter for D15 vs D28, Mtb+, those two celltypes, lung+mln, our genes, sig only
df <- deg_time %>%
  filter(
    comparison_family == "time",
    group1 == "15", group2 == "28",
    Mtb == "pos",
    tissue %in% c("lung","mln"),
    celltype %in% c("MC","iMO2"),
    gene %in% genes_of_interest,
    p_val_adj < 0.05
  ) %>%
  mutate(
    log2FC    = avg_log2FC,
    neglog10p = -log10(p_val_adj),
    # preserve your ordering:
    gene     = factor(gene, levels = genes_of_interest),
    celltype = factor(celltype, levels = c("MC","iMO2"))
  )

# plotting loop
for(tiss in unique(df$tissue)) {
  df_tiss <- df %>% filter(tissue == tiss)
  
  p <- ggplot(df_tiss, aes(x = celltype, y = gene)) +
    geom_point(aes(size = neglog10p, color = log2FC), alpha = .8) +
    facet_grid(. ~ day, scales = "free_y", space = "free_y") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          name = expression(log[2]*FC)) +
    scale_size_continuous(range = c(2, 8), name = expression(-log[10]*(adj~p))) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "grey20", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey30"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey30", color = NA),
      strip.text       = element_text(color = "white", face = "bold"),
      axis.text.x      = element_text(color = "black", angle = 45, hjust = 1),
      axis.text.y      = element_text(color = "black"),
      axis.title       = element_text(color = "black"),
      plot.title       = element_text(color = "black", hjust = 0.5, size = 16)
    ) +
    labs(
      title = sprintf("%s: D15 vs D28 (Mtb⁺)", tools::toTitleCase(tiss)),
      x     = "Cell Type",
      y     = "Gene"
    )
  
  out_pdf <- file.path(out_dir, paste0("2bubble_D15vs28_Mtbpos_", tiss, ".pdf"))
  ggsave(out_pdf, plot = p, width = 12, height = 10, limitsize = FALSE)
  message("Saved: ", out_pdf)
}
