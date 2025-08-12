### SET WD and LOAD allCells from EAS_run_DEG_GSEA_All ###

# ================================
# SCT-based DEG + GO:BP GSEA + rrvgo (mt- excluded; no split-axis plotting)
# ================================
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(fgsea)
  library(GO.db)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(rrvgo)
  library(ggplot2)
})

# ---------- helpers ----------
`%||%` <- function(a, b) if (!is.null(a)) a else b
.now <- function() Sys.time()
.log_dur <- function(start, label) {
  dur <- round(as.numeric(difftime(.now(), start, units = "secs")), 2)
  message(sprintf("%s took %ss", label, dur))
}
safe_file <- function(..., ext = "csv") {
  base <- paste(..., sep = "_")
  base <- gsub("[^A-Za-z0-9._-]+", "-", base)
  paste0(base, ".", ext)
}
is_mito_gene <- function(x) grepl("^mt-", x, ignore.case = TRUE)

# ---------- parameters ----------
padj_cutoff         <- 0.05
lfc_cutoff          <- 0.5
min_cells_grp       <- 10
top_go_keep         <- 200  # limit sent to rrvgo
force_recompute_gsea <- TRUE

# ---------- metadata hygiene ----------
if ("Mtb" %in% colnames(allCells@meta.data)) {
  allCells$Mtb <- as.character(allCells$Mtb)
  sel <- allCells$Mtb %in% c("Mtb-","Mtb+")
  allCells$Mtb[sel] <- c("neg","pos")[match(allCells$Mtb[sel], c("Mtb-","Mtb+"))]
  allCells$Mtb <- factor(allCells$Mtb, levels = c("neg","pos"))
}
allCells$day    <- factor(as.character(allCells$day))
allCells$tissue <- factor(as.character(allCells$tissue))

celltypes   <- sort(unique(allCells$ImmGenGroup))
days_all    <- levels(allCells$day)
tissues_all <- levels(allCells$tissue)
mtb_levels  <- levels(allCells$Mtb) %||% c("neg","pos")

# ---------- output folders ----------
dir.create("per_comparison/DEG",             showWarnings = FALSE, recursive=TRUE)
dir.create("per_comparison/GSEA_GOBP",       showWarnings = FALSE, recursive=TRUE)
dir.create("per_comparison/GSEA_GOBP_rrvgo", showWarnings = FALSE, recursive=TRUE)
dir.create("merged_DEG_results", showWarnings = FALSE)
dir.create("merged_GSEA_results", showWarnings = FALSE)
dir.create("merged_GSEA_reduced_results", showWarnings = FALSE)
dir.create("GSEA_plots", showWarnings = FALSE)

# ---------- GO:BP definitions (cache) ----------
message("Preparing GO:BP gene sets (using cache if available)...")
mouse_go_cache_file <- "mouse_go_cache.rds"
if (file.exists(mouse_go_cache_file)) {
  mouse_go <- readRDS(mouse_go_cache_file)
} else {
  mouse_go <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = keys(org.Mm.eg.db, keytype = "SYMBOL"),
    columns = c("GOALL","ONTOLOGY"),
    keytype = "SYMBOL"
  ) %>% distinct()
  saveRDS(mouse_go, mouse_go_cache_file)
}
bp_root        <- "GO:0008150"
bp_ids         <- unique(c(bp_root, unlist(as.list(GOBPOFFSPRING))))
bp_go_ids_mouse<- intersect(unique(mouse_go$GOALL), bp_ids)
go_terms_df    <- AnnotationDbi::select(
  GO.db,
  keys    = bp_go_ids_mouse,
  columns = "TERM",
  keytype = "GOID"
) %>% distinct(GOID,TERM)
mouse_go_bp    <- mouse_go %>% filter(ONTOLOGY=="BP", GOALL%in%bp_go_ids_mouse) %>% distinct(SYMBOL,GOALL)
gmt_go_bp      <- split(mouse_go_bp$SYMBOL, mouse_go_bp$GOALL)
gs_id_to_name  <- go_terms_df %>% filter(GOID%in%names(gmt_go_bp)) %>% distinct(GOID,TERM) %>% deframe()
gmt_gobp       <- gmt_go_bp

# ---------- rrvgo reducer ----------
reduce_gsea_terms <- function(gsea_df, similarity_cutoff = 0.7) {
  if (is.null(gsea_df)||!nrow(gsea_df)) return(gsea_df%>%mutate(parentTerm=NA_character_))
  key_col <- if ("gs_id"%in%names(gsea_df)) "gs_id" else "pathway"
  go_terms <- unique(gsea_df[[key_col]][grepl("^GO:",gsea_df[[key_col]])])
  if (length(go_terms)<2) return(gsea_df%>%mutate(parentTerm=NA_character_))
  sim <- tryCatch(rrvgo::calculateSimMatrix(go_terms, orgdb="org.Mm.eg.db", ont="BP", method="Rel"), error=function(e) NULL)
  if (is.null(sim)) return(gsea_df%>%mutate(parentTerm=NA_character_))
  scores <- setNames(abs(gsea_df$NES[match(go_terms, gsea_df[[key_col]])]), go_terms)
  red <- tryCatch(rrvgo::reduceSimMatrix(sim, scores, threshold=similarity_cutoff, orgdb=org.Mm.eg.db), error=function(e) NULL)
  if (is.null(red)||!nrow(red)) return(gsea_df%>%mutate(parentTerm=NA_character_))
  if (!"parentTerm"%in%names(red)) {
    red$parentTerm <- gs_id_to_name[red$parent]; red$parentTerm[is.na(red$parentTerm)] <- red$parent
  }
  out <- left_join(gsea_df, red[,c("go","parentTerm")], by=setNames("go",key_col))
  if (all(is.na(out$parentTerm))) message("rrvgo: no matches")
  out
}

# ---------- per-comparison writer ----------
write_per_comparison <- function(comp_name, deg=NULL, gsea=NULL, gsea_red=NULL) {
  if (!is.null(deg)     && nrow(deg))     write_csv(deg,     file.path("per_comparison/DEG",         safe_file(comp_name))) 
  if (!is.null(gsea)    && nrow(gsea))    write_csv(gsea,    file.path("per_comparison/GSEA_GOBP",   safe_file(comp_name))) 
  if (!is.null(gsea_red)&& nrow(gsea_red))write_csv(gsea_red, file.path("per_comparison/GSEA_GOBP_rrvgo", safe_file(comp_name))) 
}

# ---------- Bubble plot ----------
plot_gsea_bubble <- function(gsea_red, output_path) {
  req_cols <- c("parentTerm","NES","pval")
  if (!all(req_cols%in%colnames(gsea_red))) return(NULL)
  d <- gsea_red %>% filter(!is.na(parentTerm), pval>0)
  if (!nrow(d)) return(NULL)
  title_str <- if ("comparison"%in%names(d)) unique(d$comparison)[1] else basename(output_path)
  ord <- d %>% group_by(parentTerm) %>% summarize(m=mean(NES)) %>% arrange(m) %>% pull(parentTerm)
  d$parentTerm <- factor(d$parentTerm, levels=ord)
  h <- max(4, min(0.35*length(ord), 12))
  p <- ggplot(d, aes(NES,parentTerm)) +
    geom_point(aes(size=-log10(pval),color=NES)) +
    scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0) +
    scale_size_continuous(range=c(2,10)) +
    theme_minimal(base_size=12) +
    labs(title=title_str, x="NES", y="GO term", color="NES", size="-log10(p)") +
    theme(plot.title=element_text(hjust=0))
  ggsave(output_path, p, width=10, height=h, limitsize=FALSE)
}

# ---------- DEG+GSEA helper ----------
run_deg_gsea <- function(seu, id_col,id1,id2, meta_cols) {
  comp_name <- paste0(meta_cols$comparison_family %||% id_col, "_",
                      meta_cols$tissue  %||% "", "_",
                      meta_cols$day     %||% "", "_",
                      meta_cols$Mtb     %||% "", "_",
                      meta_cols$celltype%||% "", "_",
                      id_col,"-",id1,"vs",id2) %>%
    gsub("__+","_",.) %>% gsub("_$","",.)
  deg_p  <- file.path("per_comparison/DEG",         safe_file(comp_name))
  gsea_p <- file.path("per_comparison/GSEA_GOBP",   safe_file(comp_name))
  rrv_p  <- file.path("per_comparison/GSEA_GOBP_rrvgo", safe_file(comp_name))
  if (file.exists(deg_p)&&file.exists(gsea_p)&&file.exists(rrv_p)&&!force_recompute_gsea) {
    message("✔ skip ",comp_name); return(NULL)
  }
  message("---- running ",comp_name)
  # 1) DEG
  if (file.exists(deg_p)) {
    deg <- read_csv(deg_p,show_col_types=FALSE)
  } else {
    idv <- as.character(seu[[id_col,drop=TRUE]]); Idents(seu)<-factor(idv)
    if (any(is.na(match(c(id1,id2),levels(Idents(seu))))) ||
        min(table(Idents(seu))[c(id1,id2)])<min_cells_grp) {
      message(" ✘ not enough cells → skip "); return(NULL)
    }
    t0 <- .now(); seu2 <- PrepSCTFindMarkers(seu); .log_dur(t0,"PrepSCTFindMarkers")
    t1 <- .now()
    deg <- FindMarkers(seu2,ident.1=id1,ident.2=id2,assay="SCT",test.use="negbinom",
                       min.pct=0.2,logfc.threshold=lfc_cutoff,recorrect_umi=FALSE)%>%
      rownames_to_column("gene")
    .log_dur(t1,sprintf("FindMarkers %s vs %s",id1,id2))
    deg <- deg %>% filter(!is_mito_gene(gene))
    for (nm in names(meta_cols)) deg[[nm]] <- meta_cols[[nm]]
    write_per_comparison(comp_name, deg=deg)
  }
  # 2) GSEA
  if (file.exists(gsea_p)&&!file.exists(rrv_p)) {
    gsea <- read_csv(gsea_p,show_col_types=FALSE)
  } else if (!file.exists(gsea_p)) {
    deg_sig <- deg %>% filter(p_val_adj<padj_cutoff, abs(avg_log2FC)>=lfc_cutoff)
    if (nrow(deg_sig)<15) { message(" ✘ too few DEGs"); return(list(deg=deg,gsea=NULL,gsea_red=NULL)) }
    stats <- deframe(deg_sig %>% filter(!is_mito_gene(gene)) %>% select(gene,avg_log2FC))
    t2 <- .now()
    gsea <- fgseaMultilevel(gmt_gobp,stats=stats)
    .log_dur(t2,"fgsea")
    if (!nrow(gsea)) return(list(deg=deg,gsea=NULL,gsea_red=NULL))
    gsea <- gsea %>% mutate(leadingEdge=sapply(leadingEdge,paste,collapse=",")) %>%
      mutate(term=gs_id_to_name[pathway]) %>%
      mutate(!!!meta_cols) %>%
      slice_head(n=top_go_keep)
    write_per_comparison(comp_name, gsea=gsea)
  } else {
    message("✔ GSEA exists & rrvgo exists → skip")
    return(NULL)
  }
  # 3) rrvgo
  if (!file.exists(rrv_p)) {
    t3 <- .now()
    gsea_red <- reduce_gsea_terms(gsea)
    .log_dur(t3,"rrvgo")
    if (nrow(gsea_red)>0) write_per_comparison(comp_name, gsea_red=gsea_red)
  }
  list(deg=deg, gsea=gsea, gsea_red=if(exists("gsea_red")) gsea_red else NULL)
}


# ───────────────────────────────────────────────────────────────────────────────
# 0) Baseline vs D0 (only Mtb="pos" at test day vs Mtb="neg" day 0)
# ───────────────────────────────────────────────────────────────────────────────
message("Running baseline vs D0 (Mtb pos @15/28 vs neg @0)...")
for (ct in celltypes) {
  for (tiss in tissues_all) {
    # only test Mtb="pos"
    for (test_day in c("15","28")) {
      idx_test <- WhichCells(allCells,
                             expression=ImmGenGroup==ct & tissue==tiss & Mtb=="pos" & day==test_day)
      idx_ref  <- WhichCells(allCells,
                             expression=ImmGenGroup==ct & tissue==tiss & Mtb=="neg" & day=="0")
      if (length(idx_test)<min_cells_grp || length(idx_ref)<min_cells_grp) next
      base <- subset(allCells, cells=c(idx_test,idx_ref))
      res <- run_deg_gsea(
        seu     = base,
        id_col  = "day",
        id1     = test_day,
        id2     = "0",
        meta_cols = list(
          comparison_family = "time_vs0",
          group1            = test_day,
          group2            = "0",
          tissue            = tiss,
          Mtb               = "pos",
          celltype          = ct
        )
      )
      if (is.null(res)) next
    }
  }
}

# ───────────────────────────────────────────────────────────────────────────────
# 1) Timepoint 15 vs 28, 2) Tissue pairs, 3) Mtb pos vs neg
# (unchanged from your script)
# ───────────────────────────────────────────────────────────────────────────────
# ... [your existing loops for comparison_family="time", "tissue", "mtb"] ...

# ───────────────────────────────────────────────────────────────────────────────
# Plot all per‐comparison rrvgo results
# ───────────────────────────────────────────────────────────────────────────────
message("Plotting bubble‐plots…")
csvs <- list.files("per_comparison/GSEA_GOBP_rrvgo", pattern="\\.csv$", full.names=TRUE)
for(f in csvs){
  df <- read_csv(f, show_col_types=FALSE)
  out  <- file.path("GSEA_plots", paste0(tools::file_path_sans_ext(basename(f)),"_bubble.pdf"))
  plot_gsea_bubble(df,out)
}

# ───────────────────────────────────────────────────────────────────────────────
# Merge & write final CSVs (including time_vs0 into the “time” file)
# ───────────────────────────────────────────────────────────────────────────────
merge_configs <- list(
  DEG = list(in_dir="per_comparison/DEG",             out="merged_DEG_results/merged_DEG_time.csv",            by="comparison_family", keep=c("time","time_vs0")),
  TIS = list(in_dir="per_comparison/GSEA_GOBP",       out="merged_GSEA_results/merged_GSEA_GOBP_time.csv",    by="comparison_family", keep=c("time","time_vs0")),
  RRV = list(in_dir="per_comparison/GSEA_GOBP_rrvgo", out="merged_GSEA_reduced_results/merged_GSEA_reduced_GOBP_time.csv", by="comparison_family", keep=c("time","time_vs0"))
)
for(cfg in merge_configs){
  files <- list.files(cfg$in_dir,pattern="\\.csv$",full.names=TRUE)
  dfs   <- lapply(files,function(f){
    df<-read_csv(f,show_col_types=FALSE)
    if("group1"%in%names(df)) df$group1<-as.character(df$group1)
    if("group2"%in%names(df)) df$group2<-as.character(df$group2)
    df
  })
  big   <- bind_rows(dfs)
  sel   <- big %>% filter(.data[[cfg$by]] %in% cfg$keep)
  write_csv(sel, cfg$out)
  message("Wrote ",cfg$out)
}

message("✅ Finished all comparisons, merges, and plots.")
