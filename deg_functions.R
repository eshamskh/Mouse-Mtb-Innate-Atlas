#Subset and DEG functions
library("Seurat")
library("tidyverse")
library("fgsea")
library("tmod")
library("msigdbr")
library("ggplot2")
library("ggrepel")
library("patchwork")
library("clustree")
library("ggalluvial")
library("conflicted")
library("readxl")
library("UCell")
conflicts_prefer(dplyr::count, dplyr::filter, dplyr::select, base::setdiff)

# List of required packages
required_packages <- c(
  "Seurat", "tidyverse", "fgsea", "tmod", "msigdbr", "ggplot2",
  "ggrepel", "patchwork", "clustree", "ggalluvial",
  "conflicted", "readxl", "UCell"
)

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Function to install and load packages
install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "..."))
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }, error = function(e) {
      message(paste("Failed to install", pkg, ":", e$message))
    })
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Apply function to all packages
invisible(lapply(required_packages, install_and_load))


IMMGEN_GROUPINGS_1 <- list(Neutrophil=c("Neutrophils (GN.ARTH)", "Neutrophils (GN.Thio)", "Neutrophils (GN)"),
                           MDM = c("Macrophages (MF.103-11B+24-)", "Macrophages (MF.II+480LO)" ),
                           AM = "Macrophages (MF)",
                           Monocyte1 = c("Monocytes (MO.6C-II-)"),
                           Monocyte2 = c("Monocytes (MO.6C-IIINT)", "Monocytes (MO.6C+II-)"),
                           DC1="DC (DC.8-)",
                           DC2=c("DC (DC.8-4-11B+)", "DC (DC)"),
                           DCother=c("DC (DC.103-11B+F4-80LO.KD)", "DC (DC.PDC.8+)", "DC (DC.8-4-11B-)"),
                           Basophil="Basophils (BA)")

IMMGEN_GROUPINGS_2 <- list(Neutrophil=c("Neutrophils (GN.ARTH)", "Neutrophils (GN.Thio)", "Neutrophils (GN)"),
                           MDM = c("Macrophages (MF.103-11B+24-)", "Macrophages (MF.II+480LO)" ),
                           AM = "Macrophages (MF)",
                           Monocyte1 = c("Monocytes (MO.6C-II-)"),
                           Monocyte2a = list(ImmgenClustersFine=c("Monocytes (MO.6C-IIINT)", "Monocytes (MO.6C+II-)"), 
                                             seurat_clusters=c(3,4,9)),
                           Monocyte2b = list(ImmgenClustersFine = c("Monocytes (MO.6C-IIINT)", "Monocytes (MO.6C+II-)"),
                                             seurat_clusters = c(0,11,14)),
                           DC1="DC (DC.8-)",
                           DC2=c("DC (DC.8-4-11B+)", "DC (DC)"),
                           DCa="DC (DC.103-11B+F4-80LO.KD)",
                           DCb="DC (DC.8-4-11B-)",
                           pDC="DC (DC.PDC.8+)",
                           Basophil="Basophils (BA)")

IMMGEN_GROUPINGS_3 <- list(Neutrophil=c("Neutrophils (GN.ARTH)", "Neutrophils (GN.Thio)", "Neutrophils (GN)"),
                           MC = c("Macrophages (MF.103-11B+24-)", "Macrophages (MF.II+480LO)" ),
                           AM = "Macrophages (MF)",
                           pMO = c("Monocytes (MO.6C-II-)"),
                           iMO1 = list(ImmgenClustersFine=c("Monocytes (MO.6C-IIINT)", "Monocytes (MO.6C+II-)"), 
                                       seurat_clusters =c(3,4,9)),
                           iMO2 = list(ImmgenClustersFine = c("Monocytes (MO.6C-IIINT)", "Monocytes (MO.6C+II-)"),
                                       seurat_clusters = c(0,11)),
                           iMO3 = list(ImmgenClustersFine = c("Monocytes (MO.6C-IIINT)", "Monocytes (MO.6C+II-)"),
                                       seurat_clusters = c(14)),
                           resDC1="DC (DC.8-)",
                           migDC1="DC (DC)",
                           migDC2="DC (DC.8-4-11B+)",
                           resDC2="DC (DC.103-11B+F4-80LO.KD)",
                           tDC="DC (DC.8-4-11B-)",
                           pDC="DC (DC.PDC.8+)",
                           Basophil="Basophils (BA)",
                           Unk="Unk")

IMMGEN_3_COLORS <- c(iMO1="maroon3", iMO2="red", iMO3 = "pink", MC="burlywood2",
                     AM="pink4", pMO="gold2", 
                     resDC2="darkblue", transDC="plum3", pDC="turquoise4",
                     migDC2="lightblue", resDC1="darkgreen", migDC1="lightgreen",
                     Neutrophil="gray87", Basophil="yellow", Unk="grey1")
ImmGenGroup_COLORS <- c(iMO1="maroon3", iMO2="red", iMO3 = "pink", MC="burlywood2",
                     AM="pink4", pMO="gold2", 
                     resDC2="darkblue", tDC="plum3", pDC="turquoise4",
                     migDC2="lightblue", resDC1="darkgreen", migDC1="lightgreen",
                     Neutrophil="gray87", Basophil="yellow", Unk="grey1")

Russell_COLORS <- c(
  "Monocytes"          = "maroon3",
  "IM_1"               = "red",
  "IM_2"               = "burlywood2",
  "IM_3"               = "burlywood2",
  "AM_1"               = "pink4",
  "AM_2"               = "lightblue1",
  "AM_3"               = "pink4",
  "AM_4"               = "green3",
  "DC.103+11B-"        = "lightblue1",
  "DC"                 = "green3",
  "Neutrophils"        = "gray87",
  "Macrophages/B-cells"= "yellow",
  "Macrophages/T-cells"= "yellow"
)

DimPlot(allCells, group.by = "russell_mac_celltypist_model", cols = Russell_COLORS)

getMarkers <- function(scObj, celltypeCol, saveFile) {
  tryCatch({
    read_csv(saveFile)
  }, error=\(e){
    message("Could not read file ", saveFile, " calculating markers...")
    message(as.character(e))
    oldIdents <- Idents(scObj)
    Idents(scObj) <- celltypeCol
    markers <- FindAllMarkers(scObj,
                              test.use="roc",
                              min.pct=0.2,
                              verbose=T) %>%
      data.frame(check.names = F) %>%
      rownames_to_column("rownames")
    Idents(scObj) <- oldIdents
    write_csv(markers, saveFile)
    read_csv(saveFile)
  })
}

celltypeVsCelltypeMarkers <- function(scObj, celltypePairs, celltypeCol, saveFile) {
  tryCatch({
    read_csv(saveFile)
  }, error=\(e){
    message("Could not read file ", saveFile, " calculating markers...")
    message(as.character(e))
    map(celltypePairs, \(ct){
      stopifnot(length(ct)==2 & all(ct %in% scObj[[celltypeCol,drop=T]]))
    })
    oldIdents <- Idents(scObj)
    Idents(scObj) <- celltypeCol
    markers <- map_dfr(celltypePairs, \(celltypes){
      message(sprintf("Finding markers for %s vs %s", celltypes[1], celltypes[2]))
      FindMarkers(scObj, 
                  ident.1=celltypes[1], 
                  ident.2=celltypes[2],
                  test.use="roc",
                  min.pct=0.2,
                  verbose=T) %>%
        data.frame(check.names = F) %>%
        rownames_to_column("rownames") %>%
        mutate(celltype.1=celltypes[1],
               celltype.2=celltypes[2])
    })

    Idents(scObj) <- oldIdents
    write_csv(markers, saveFile)
    read_csv(saveFile)
  })
}

addGroupings <- function(scObj, groupingDefs, groupLabel="ImmGenGroup", seuratClustCol="integrated_snn_res.0.5") {

  groupingMap <- imap_dfr(groupingDefs, \(cellTypes, group){
    if (is.character(cellTypes)) {
      data.frame(ImmgenClustersFine=cellTypes, ImmGenGroup=group)
    } else {
      combs <- expand.grid(cellTypes$ImmgenClustersFine, cellTypes$seurat_clusters) %>%
        {paste(.[[1]], .[[2]], sep="|")}
      data.frame(ImmgenClustersFine=combs, ImmGenGroup=group)
    }
  }) %>% {set_names(.$ImmGenGroup, .$ImmgenClustersFine)}
  shortLabel <- scObj$ImmgenClustersFine
  fullLabel <- paste(scObj$ImmgenClustersFine, as.numeric(as.character(scObj[[seuratClustCol, drop=T]])), sep="|")
  newGroupLabel <- ifelse(is.na(groupingMap[shortLabel]), groupingMap[fullLabel], groupingMap[shortLabel]) %>%
    {ifelse(is.na(.), "unk", .)} %>%
    set_names(colnames(scObj))
  AddMetaData(scObj,
              newGroupLabel,
              col.name=groupLabel)
}

getTissueDegs <- function(cells,
                         compTissues,
                         times,
                         infectStatuses,
                         cellTypeCol = "lung_mouse_mtb_immune_celltypist_model",
                         saveFile=NULL) {
  tryCatch({
    read_csv(saveFile)
  }, error=\(e){
    stopifnot(length(compTissues)==2 && all(compTissues %in% cells$tissue))
    comps <- expand_grid(time=times,
                         celltype=unique(cells[[cellTypeCol, drop=T]]),
                         infection=infectStatuses)
    cells <- PrepSCTFindMarkers(cells)
    markers <- pmap_dfr(comps, \(infection, time, celltype){
      message(sprintf("Calculating DEGs for %s: %s vs %s, Mtb %s, day %d", celltype, compTissues[1], compTissues[2], infection, time))
      
      filtCells <- tryCatch({
        seleCells <- cells %>%
          {.[,.$tissue %in% compTissues & .$day==time & .[[cellTypeCol,drop=T]]==celltype]}
        if (!all(is.na(infection))) {
          seleCells <- seleCells[,seleCells$infect %in% infection]
        }
        seleCells
      }, error=\(e){
        message(as.character(e))
        NULL
      })
      #Need both Mtb+ and Mtb- cells to look for DEGs
      if (is.null(filtCells) || is.na(table(filtCells$tissue)[compTissues[1]]) || is.na(table(filtCells$tissue)[compTissues[2]]) || min(table(filtCells$tissue)) < 10) {
        return(data.frame())
      }
      Idents(filtCells) <- "tissue"
      mrk <- FindMarkers(filtCells,
                         ident.1=compTissues[1],
                         ident.2=compTissues[2],
                         assay="SCT",
                         test.use="negbinom",
                         min.pct=0.2,
                         logfc.threshold = 0.5,
                         recorrect_umi=FALSE) %>%
        data.frame(check.names = F) %>%
        rownames_to_column("genes") %>%
        mutate(tissues=paste(compTissues, collapse="_"), time=time, celltype=celltype, infect=infection)
      print(head(mrk))
      mrk
    })
    if (!is.null(saveFile)) {
      write_csv(markers, saveFile)
      return(read_csv(saveFile))
    }
    markers
  })
}

getMtbDegs <- function(cells, 
                       mtbTissues,
                       mtbTimes,
                       cellTypeCol = "lung_mouse_mtb_immune_celltypist_model",
                       saveFile=NULL) {
  tryCatch({
    read_csv(saveFile)
  }, error=\(e){
    comps <- expand_grid(tissue=mtbTissues, 
                         time=mtbTimes, 
                         celltype=unique(cells[[cellTypeCol, drop=T]]))  
    cells <- PrepSCTFindMarkers(cells)
    markers <- pmap_dfr(comps, \(tissue, time, celltype){
      message("Calculating DEGs for :", tissue, ", ", time, ", ", celltype)
      
      filtCells <- tryCatch({
        cells %>%
          {.[,.$tissue==tissue & .$day==time & .[[cellTypeCol,drop=T]]==celltype]}
      }, error=\(e){
        message(as.character(e))
        NULL
      })
      #Need both Mtb+ and Mtb- cells to look for DEGs
      if (is.null(filtCells) || is.na(table(filtCells$infect)["neg"]) || is.na(table(filtCells$infect)["pos"]) || min(table(filtCells$infect)) < 10) {
        return(data.frame())
      }
      Idents(filtCells) <- "infect"
      mrk <- FindMarkers(filtCells,
                         ident.1="pos",
                         ident.2="neg",
                         assay="SCT",
                         test.use="negbinom",
                         min.pct=0.2, 
                         logfc.threshold = 0.5,
                         recorrect_umi=FALSE) %>%
        data.frame(check.names = F) %>%
        rownames_to_column("genes") %>%
        mutate(tissue=tissue, time=time, celltype=celltype)
      print(head(mrk))
      mrk
    })
    if (!is.null(saveFile)) {
      write_csv(markers, saveFile)
      return(read_csv(saveFile))
    }
    markers
  })
}

getTimeDegs <- function(cells,
                        compTimes,
                        tissues,
                        infectStatuses,
                        cellTypeCol = "lung_mouse_mtb_immune_celltypist_model",
                        saveFile=NULL) {
  tryCatch({
    read_csv(saveFile)
  }, error=\(e){
    stopifnot(length(compTimes)==2 && all(compTimes %in% cells$day))
    comps <- expand_grid(tissue=tissues,
                         celltype=unique(cells[[cellTypeCol, drop=T]]))
    cells <- PrepSCTFindMarkers(cells)
    markers <- pmap_dfr(comps, \(tissue, celltype, infection){
      message(sprintf("Calculating DEGs for %s: %s %d vs %d, Mtb %s", celltype, tissue, compTimes["t1"], compTimes["t2"], paste(infectStatuses, collapse="/")))
      
      filtCells <- tryCatch({
        cells %>%
          {.[,.$day %in% compTimes & .$tissue==tissue & .[[cellTypeCol,drop=T]]==celltype & .$infect %in% infectStatuses]}
      }, error=\(e){
        message(as.character(e))
        NULL
      })
      #Need both Mtb+ and Mtb- cells to look for DEGs
      if (is.null(filtCells) || is.na(table(filtCells$day)[as.character(compTimes["t1"])]) || is.na(table(filtCells$day)[as.character(compTimes["t2"])]) || min(table(filtCells$day)) < 10) {
        message("Not enough cells, could not calculate marker genes...")
        return(data.frame())
      }
      Idents(filtCells) <- "day"
      mrk <- FindMarkers(filtCells,
                         ident.1=compTimes["t1"],
                         ident.2=compTimes["t2"],
                         assay="SCT",
                         test.use="negbinom",
                         min.pct=0.2,
                         logfc.threshold = 0.5,
                         recorrect_umi=FALSE) %>%
        data.frame(check.names = F) %>%
        rownames_to_column("genes") %>%
        mutate(times=paste(compTimes, collapse="_"), 
               tissue=tissue, 
               celltype=celltype, 
               infect=paste(infectStatuses, collapse="/"))
      print(head(mrk))
      mrk
    })
    if (!is.null(saveFile)) {
      write_csv(markers, saveFile)
      return(read_csv(saveFile))
    }
    markers
  })
}

downsampleSeuratClusters <- function(scObj, clustCol, nPerClust=1000) {
  keepIdxs <- data.frame(clust=scObj[[clustCol,drop=T]]) %>%
    mutate(cellIdx=1:n()) %>%
    group_by(.data[[clustCol]]) %>%
    slice_sample(n=nPerClust) %>%
    {.$cellIdx}
  scObj[,keepIdxs]
}

humanGenesToMouse <- function(humanGenes, homologyFile="mouse_homology/HOM_MouseHumanSequence.txt") {
  mapping <- read_tsv(homologyFile) 
  humanKeys <- mapping %>%
    filter(`Common Organism Name` =="human" & Symbol %in% humanGenes) %>%
    {.$`DB Class Key`}
  mapping %>%
    filter(`Common Organism Name`=="mouse, laboratory" & `DB Class Key` %in% humanKeys) %>%
    {unique(.$Symbol)}
  
}

getMouseBTMs <- function(maxMissingGeneFrac=0.2) {
  data("tmod")
  humanBtms <- tmod$gs2gv %>% 
    map(\(gIds)tmod$gv[gIds]) %>% 
    set_names(paste(tmod$gs$ID, tmod$gs$Title)) 
  imap(humanBtms, \(gs, gsName){
    mgs <- humanGenesToMouse(gs)
    matchedFrac <- (length(mgs)/length(gs))
    if (1-matchedFrac > maxMissingGeneFrac) {
      message(sprintf("Dropping %s BTM, not enough human -> mouse homologs", gsName))
      return(c())
    }
    mgs
  }) %>% compact()
}

getMousePathways <- function(saveFile="analysis/mouse_pathways.rds") {
  tryCatch({
    read_rds(saveFile)
  }, error=\(e){
    message(as.character(e))
    message("Compiling mouse pathway definitions and saving...")
    hallmark <- msigdbr(species="Mus musculus", category="H")
    gobp <- msigdbr(species="Mus musculus", category="C5", subcategory="GO:BP")
    #wp <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:WIKIPATHWAYS")
    #reactome <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:REACTOME")
    #kegg <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG")
    #mouseBtms <- getMouseBTMs()
    mouseMsigDBpathways <- bind_rows(hallmark, gobp) %>%
      select(gs_name, gene_symbol) %>%
      distinct() %>%
      split(.$gs_name) %>%
      map(\(df)df$gene_symbol)
    pathways <- c(mouseMsigDBpathways)
    write_rds(pathways, saveFile)
    read_rds(saveFile)
  })
}

degGseaEnrichment <- function(degResults, groupingFactors, saveFile) {
  tryCatch({
    read_rds(saveFile)
  }, error=\(e){
    message(as.character(e))
    message("Running GSEA...")
    pathways <- getMousePathways()
    res <- degResults %>% 
      group_by(across(one_of(groupingFactors))) %>% 
      group_modify(\(df, grp){
        stats <- df %>%
          arrange(dplyr::desc(avg_log2FC)) %>%
          {set_names(.$avg_log2FC, .$genes)}
        fgsea(pathways, stats)
      })
    write_rds(res, saveFile)
    read_rds(saveFile)
  })

}


clustDotplotGseaHeatmap <- function(gseaRes, 
                                    allCellTypes,
                                    mtbPosCellTypes,
                                    fdrThresh=0.05, 
                                    trimPathNames=T,
                                    pointSizeRange=c(2,5)) {
  
  
  plotDF <- gseaRes %>%
    {if (trimPathNames) mutate(., pathway=sub("^[A-Z]+_", "", pathway)) else .} %>%
    filter(celltype %in% allCellTypes) %>%
    mutate(padj=p.adjust(pval, "fdr")) %>%
    mutate(infect=ifelse(grepl("Mtb+", comparison, fixed=T), "Mtb+", "Mtb-")) %>%
    mutate(celltype=factor(celltype, levels=allCellTypes)) %>%
    ungroup() %>%
    complete(pathway, tissue, celltype, nesting(infect, comparison)) %>%
    filter(infect=="Mtb-" | celltype %in% mtbPosCellTypes) %>%
    filter(!is.na(comparison) & !(tissue=="blood" & infect=="Mtb+")) %>%
    mutate(time=str_split_i(comparison, ", ", 1))
  
  plotClust <- plotDF %>%
    filter(!is.na(NES)) %>%
    mutate(label=paste(celltype, tissue, time, infect)) %>%
    select(pathway, NES, label) %>%
    pivot_wider(names_from=label, values_from=NES, values_fill=0) %>%
    column_to_rownames("pathway") %>%
    dist() %>%
    hclust(method = "ward.D2")
  
  plotDF <- filter(plotDF, padj<fdrThresh)
    
  pathwayOrder <- base::intersect(plotClust$labels[plotClust$order],
                            unique(plotDF$pathway))
    
  ggplot(plotDF, aes(x=celltype, y=pathway, colour=NES, size=-log10(padj))) +
    geom_point() +
    facet_grid(~tissue+time+infect, scales="free_x", space="free") +
    scale_colour_distiller(palette="RdBu", limits=c(-3, 3), oob=scales::squish) +
    scale_size_continuous(range=pointSizeRange) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.background = element_rect(fill="grey20"),
          panel.grid = element_blank()) +
    scale_y_discrete(limits=pathwayOrder) +
    labs(x=NULL, y=NULL, caption=sprintf("Showing FDR < %.2f", fdrThresh))
}


####New GSEA Heatmap
clustDotplotGseaHeatmap <- function(
    gseaRes, 
    allCellTypes,
    mtbPosCellTypes,
    fdrThresh = 0.05, 
    trimPathNames = TRUE,
    pointSizeRange = c(2,5),
    # NEW: custom pathway order (character vec, no "HALLMARK_" prefix)
    pathway_order = NULL,
    # NEW: custom block order (character vec of "tissue Mtb± D##")
    facet_block_order = NULL,
    # NEW: make axis tick labels bold
    bold_axes = TRUE,
    y_expand_mult = 0.001  # keep rows tight without changing text size
) {
  # Prepare base DF (same as before, with a couple of extras)
  plotDF <- gseaRes %>%
    { if (trimPathNames) mutate(., pathway = sub("^[A-Z]+_", "", pathway)) else . } %>%
    filter(celltype %in% allCellTypes) %>%
    mutate(padj  = p.adjust(pval, "fdr")) %>%
    mutate(infect = ifelse(grepl("Mtb\\+", comparison, perl = TRUE), "Mtb+", "Mtb-")) %>%
    mutate(celltype = factor(celltype, levels = allCellTypes)) %>%
    ungroup() %>%
    complete(pathway, tissue, celltype, nesting(infect, comparison)) %>%
    filter(infect == "Mtb-" | celltype %in% mtbPosCellTypes) %>%
    filter(!is.na(comparison) & !(tissue == "blood" & infect == "Mtb+")) %>%
    # Extract a simple time like "D15", "D28" from the comparison field
    mutate(time_full = stringr::str_split_i(comparison, ", ", 1),
           time = stringr::str_extract(time_full, "D\\d+")) %>%
    # Build a single “block” variable for faceting, e.g. "lung Mtb- D15"
    mutate(block = paste(tissue, infect, time))
  
  # Filter by FDR
  plotDF <- filter(plotDF, !is.na(NES) & padj < fdrThresh)
  
  # Apply custom pathway order if provided
  if (!is.null(pathway_order)) {
    # keep only requested pathways, in the specified order
    plotDF <- dplyr::filter(plotDF, pathway %in% pathway_order) %>%
      mutate(pathway = factor(pathway, levels = rev(pathway_order)))  # top-to-bottom in your requested order
  } else {
    # fall back to the old clustering-based order
    pathwayOrder <- plotDF %>%
      mutate(label = paste(celltype, tissue, time, infect)) %>%
      select(pathway, NES, label) %>%
      tidyr::pivot_wider(names_from = label, values_from = NES, values_fill = 0) %>%
      tibble::column_to_rownames("pathway") %>%
      dist() %>%
      hclust(method = "ward.D2") %>%
      {\(hc) base::intersect(hc$labels[hc$order], unique(plotDF$pathway))}()
    plotDF <- mutate(plotDF, pathway = factor(pathway, levels = pathwayOrder))
  }
  
  # Apply custom facet block order if provided
  if (!is.null(facet_block_order)) {
    plotDF <- plotDF %>%
      filter(block %in% facet_block_order) %>%
      mutate(block = factor(block, levels = facet_block_order))
  }
  
  p <- ggplot(plotDF, aes(x = celltype, y = pathway, colour = NES, size = -log10(padj))) +
    geom_point() +
    # Facet **by the single ordered block** variable
    facet_grid(~ block, scales = "free_x", space = "free") +
    scale_colour_distiller(palette = "RdBu", limits = c(-3, 3), oob = scales::squish) +
    scale_size_continuous(range = pointSizeRange) +
    # ↓ tighten top/bottom padding for discrete y (does not change label size)
    scale_y_discrete(expand = expansion(mult = c(y_expand_mult, y_expand_mult))) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 face = if (bold_axes) "bold" else "plain"),
      axis.text.y = element_text(face = if (bold_axes) "bold" else "plain"),
      panel.background = element_rect(fill = "grey20"),
      panel.grid = element_blank(),
      strip.text.x = element_text(face = if (bold_axes) "bold" else "plain")
    ) +
    labs(x = NULL, y = NULL, caption = sprintf("Showing FDR < %.2f", fdrThresh))
  
  # attach how many y rows are actually plotted so we can size the figure
  n_unique_pathways <- nlevels(plotDF$pathway)
  attr(p, "n_pathways") <- n_unique_pathways
  
  p
}



clustDotGeneHeatmap <- function(degDF,
                                genesToPlot,
                                allCellTypes,
                                mtbPosCellTypes,
                                minFCdiff=2,
                                minPctDiff=0.1,
                                inflabels= c(`uninf/pos`="Mtb+", 
                                             `neg/uninf/inf`="Mtb-"),
                                timelabels=c(`15_0`="D15 vs D0", 
                                             `28_0`="D28 vs D0", 
                                             `28_15`="D28 vs D15"),
                                pointSizeRange=c(2,5)) {
  plotData <- filter(degDF, genes %in% genesToPlot & celltype %in% allCellTypes)
  
  if (!all(is.na(inflabels))) {
    plotData <- mutate(plotData, inflabel=inflabels[infect]) 
  } else {
    plotData <- mutate(plotData, inflabel=infect) 
  }
  
  if (!all(is.na(timelabels))) {
    plotData <- mutate(plotData, time=timelabels[times]) 
  }
  
  plotData <- filter(plotData, inflabel=="Mtb-" | celltype %in% mtbPosCellTypes)
  
  plotClust <- plotData %>%
    mutate(group=paste(tissue, celltype, inflabel, time)) %>%
    select(genes, avg_log2FC, group) %>%
    pivot_wider(names_from="group", values_from=avg_log2FC, values_fill = 0) %>%
    column_to_rownames("genes") %>%
    dist() %>%
    hclust(method="ward.D2")
    
  plotLabels <- plotClust$labels[plotClust$order]
  
  plotData %>%
    mutate(genes=factor(genes, levels=plotLabels)) %>%
    filter(abs(avg_log2FC)>minFCdiff & abs(pct.1-pct.2)>minPctDiff) %>%
    ggplot(aes(x=celltype, y=genes, colour=avg_log2FC, size=abs(pct.1-pct.2))) +
      geom_point() +
      facet_grid(~tissue+time+inflabel, space="free", scales="free_x") +
      scale_colour_distiller(palette="RdBu", limits=c(-5, 5), oob=scales::squish) +
      scale_size_continuous(range=pointSizeRange) +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            panel.background = element_rect(fill="grey20"),
            panel.grid = element_blank()) +
      labs(x=NULL, y=NULL,
           caption=sprintf("Showing only avg log2 FC > %.2f and diff %%expressing > %.2f", minFCdiff, minPctDiff))
}

getUCellSigScores <- function(seuratObj, signatures, sigSaveFile) {
  signatureColumns <- c(paste0(names(signatures), "_UCell"))
  
  if (all(signatureColumns %in% colnames(seuratObj@meta.data))) {
    message("Signatures found in exising Seurat object, returning...")
    return(seuratObj)
  }  
  
  scores <- data.frame()
  if (file.exists(sigSaveFile)) {
    scores <- read_csv(sigSaveFile) %>%
      data.frame() %>%
      column_to_rownames("cellID")
  } 
  
  if (all(signatureColumns %in% colnames(scores))) {
    message(sprintf("Found complete score file %s, adding scores...", sigSaveFile))
    seuratObj <- AddMetaData(seuratObj, scores)
  } else {
    message(sprintf("Calculating UCell signature scores: %s", paste(names(signatures), collapse=", ")))
    seuratObj <- AddModuleScore_UCell(seuratObj, features=signatures, storeRanks=T)
    scores <- seuratObj@meta.data %>%
      mutate(cellID=colnames(seuratObj)) %>%
      select(cellID, ends_with("_UCell"))
    message(sprintf("Saving scores to %s", sigSaveFile))
    write_csv(scores, sigSaveFile)
  }
  
  seuratObj
}
