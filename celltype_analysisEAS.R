# Functions ####
source("deg_functions.R")

trimName <- function(name, maxLen=6) {
  ifelse(str_length(name)<=6, name, paste0(substr(name, 1, maxLen-2), ".."))
}

strSignif <- function(pvals) {
  case_when(pvals<0.005~"***",
            pvals<0.01~"**",
            pvals<0.05~"*",
            TRUE~"n.s.")
}

sortClusterFactorLevels <- function(cl) {
  sort(as.numeric(as.character(unique(cl))))
}

# ImmGenGroupings <- list(Neutrophil=c("Neutrophils (GN.ARTH)", "Neutrophils (GN.Thio)", "Neutrophils (GN)"),
#                         MDM = c("Macrophages (MF.103-11B+24-)", "Macrophages (MF.II+480LO)" ),
#                         AM = "Macrophages (MF)",
#                         Monocyte1 = c("Monocytes (MO.6C-II-)"),
#                         Monocyte2 = c("Monocytes (MO.6C-IIINT)", "Monocytes (MO.6C+II-)"),
#                         DC1="DC (DC.8-)",
#                         DC2=c("DC (DC.8-4-11B+)", "DC (DC)"),
#                         DCother=c("DC (DC.103-11B+F4-80LO.KD)", "DC (DC.PDC.8+)", "DC (DC.8-4-11B-)"),
#                         Basophil="Basophils (BA)")

ImmGenGroupings <- IMMGEN_GROUPINGS_1


#Load the data ####

allCells <- read_rds("~/ondemand/scRNAseq/innate.int.rds")
allCells<-annotated_seurat
#Res 0.5 clusters were used to assign ImmGen labels
Idents(allCells) <- "integrated_snn_res.0.5"
# allCells <- AddMetaData(allCells, 
#                         set_names(ImmGenGroupingMap[allCells$ImmgenClustersFine], colnames(allCells)), 
#                         col.name="ImmGenGroup")

allCells <- allCells %>%
  addGroupings(groupingDefs=IMMGEN_GROUPINGS_1, groupLabel="ImmGenGroup") %>%
  addGroupings(groupingDefs=IMMGEN_GROUPINGS_2, groupLabel="ImmGenGroup_Refined") %>%
  addGroupings(groupingDefs=IMMGEN_GROUPINGS_3, groupLabel="ImmGenGroup_Refinedv2") %>%
  AddMetaData(set_names(IMMGEN_3_RENAMES[.[["ImmGenGroup_Refinedv2",drop=T]]], colnames(.)), col.name="ImmGenGroup_Refinedv2_rename")

ImmGenGroupMarkers <- getMarkers(allCells, "ImmGenGroup", "ImmgenGroupMarkers.csv")


markerList <- list(ImmGenGroup=ImmgenGroupMarkers)

### Select cell grouping column ####
analysisMarkerCol <- "ImmGenGroup"

markers <-  markerList[[analysisMarkerCol]]

monocyteTypes <- grep("^Mono", unique(allCells[[analysisMarkerCol,drop=T]]), value=T)

monocyteMarkers <- celltypeVsCelltypeMarkers(scObj=allCells,
                                            celltypePairs=combn(monocyteTypes, 2, simplify=F),
                                            celltypeCol=analysisMarkerCol,
                                            saveFile=sprintf("%s_MonocyteMarkers.csv", analysisMarkerCol))


topMarkers <- markers  %>%
  filter(!is.na(avg_diff) & !startsWith(gene, "mt-")) %>%
  group_by(cluster) %>%
  slice_max(order_by=power, n=5) %>%
  {unique(.$gene)}

monocyteTopMarkers <- monocyteMarkers %>%
  filter(!is.na(avg_diff) & !is.infinite(avg_diff) & !startsWith(rownames, "mt-")) %>%
  group_by(celltype.1, celltype.2) %>%
  slice_max(order_by=power, n=5) %>%
  {unique(.$rownames)}
#Also compare monocytes


#Plot markers and clusterings ####
plotDir <- "~/ondemand/scRNAseq/FD/scRNAseq/plots"

DotPlot(allCells[,!allCells[[analysisMarkerCol,drop=T]] %in% "unk"], 
        unique(c(topMarkers)), 
        group.by = "ImmGenGroup", assay="SCT",
        cols="YlGnBu", dot.min=0.01, dot.scale=4) +
  theme(panel.background = element_rect(fill="grey20"),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x=NULL, y=NULL) +
  coord_flip()
ggsave(sprintf(file.path(plotDir, "rename_clusterMarker_dotplot_%s.pdf"), "ImmGenGroup"), 
       width=6, height=12)

#Custom gene list
# Custom genes of interest
custom_genes <- c("Siglecf", "Sirpa", "Itgam", "Fcgr1", "Ccr2",  
                  "Ly6c2", "Itgax","Flt3","Xcr1","Ccr7",
                  "Irf4", "Cx3cr1", "Csf1r", "Nr4a1")

# Make the dot plot
DotPlot(allCells[, !allCells[[analysisMarkerCol, drop = TRUE]] %in% "unk"],
        features = custom_genes,
        group.by = analysisMarkerCol,
        assay = "SCT",
        cols = "YlGnBu",
        dot.min = 0.01,
        dot.scale = 4) +
  theme(panel.background = element_rect(fill = "grey20"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL) +
  coord_flip()

# Save to file
ggsave(file.path(plotDir, "customGene_dotplot_ImmGenGroup.pdf"),
       width = 6, height = 5)


celltypistLabels <- list.files("celltypist-predictions", full.names = T) %>%
  set_names(basename(gsub(".csv$", "", .))) %>%
  imap_dfc(\(fname, modelName){
    read_csv(fname) %>%
      select(cellBarcodeID, majority_voting) %>%
      column_to_rownames("cellBarcodeID") %>%
      rename_with(\(.)modelName, .cols=majority_voting)
  })

allCells <- AddMetaData(allCells, celltypistLabels)

selectedGoTerms <- read_file("enriched_geneset_ap_clusters_wg_result1707337515.txt") %>%
  str_split_1(pattern="[\\t\\n]") %>%
  {.[.!=""]} 

selectedGoTermNames <- msigdbr(species="Mus musculus", category="C5", subcategory="GO:BP") %>%
  select(gs_name, gs_exact_source) %>%
  filter(gs_exact_source %in% selectedGoTerms) %>%
  distinct() 

colMap <- list(ImmGenGroup_Refinedv2_rename=IMMGEN_3_RENAME_COLORS)

celltypistModels <- colnames(celltypistLabels)
clusteringNames <- c("integrated_snn_res.0.5", 
                     "ImmGenGroup") %>%
  set_names()

celltypeUmaps <- map(clusteringNames, \(clustering){
  colours <- if (clustering %in% names(colMap)) colMap[[clustering]] else NULL
  DimPlot(allCells, label=T, repel=T, cols=colours, group.by=clustering, raster = T)
}) 

imap(celltypeUmaps, \(plt, name){
  ggsave(plot=plt & guides(colour="none"), 
         filename=file.path(plotDir, sprintf("celltype_umap_%s.pdf", name)))
})

(celltypeUmaps$integrated_snn_res.0.5 + celltypeUmaps$ImmGenGroup_Refinedv2_rename) & guides(colour="none")

celltypeUmapPlot <- wrap_plots(celltypeUmaps) & 
  guides(colour="none")
ggsave(plot = celltypeUmapPlot, filename=file.path(plotDir, "celltype_umaps.pdf"), width=20, height=15)


for (clustering in clusteringNames){
  colours <- if (clustering %in% names(colMap)) colMap[[clustering]] else NULL
  condUmap <- DimPlot(allCells, 
                      cols=colours,
                      label=T, 
                      repel=T, 
                      label.size=2.5,
                      group.by=clustering, 
                      split.by = "condition", 
                      ncol = 4)
  ggsave(plot=condUmap, filename=file.path(plotDir, sprintf("condition_umap_by_%s.pdf", clustering)), width=16, height=16)
  condUmapNolabel <- DimPlot(allCells, 
                             cols=colours,
                              label=F, 
                              repel=T, 
                              label.size=2.5,
                              group.by=clustering, 
                              split.by = "condition", 
                              ncol = 4)
  ggsave(plot=condUmapNolabel, filename=file.path(plotDir, sprintf("condition_umap_nolabel_by_%s.pdf", clustering)), width=16, height=16)
  
}

#Add Mtb grouping
allCells$Mtb <- dplyr::recode(allCells$infect,
                              !!!setNames(rep("Mtb-", length(infectGroups$`Mtb-`)), infectGroups$`Mtb-`),
                              pos = "Mtb+")

#### Fig 3 B - UMAP per condition ####
infectGroups <- list(`Mtb-`=c("inf", "uninf", "neg"), `Mtb+`="pos")
timeGroups <- list(`Preinf (D0)`=0, `D15`=15, `D28`=28)

fig3bParts <- imap(infectGroups, \(infectGroups, infectGroupName){
  imap(timeGroups, \(time, timeName){
    tryCatch({
      cells <- allCells[,allCells$infect %in% infectGroups & 
                          allCells$day==time & 
                          !is.na(allCells$ImmGenGroup)]
      UMAPPlot(cells, 
               split.by="tissue", 
               group.by="ImmGenGroup",
               cols=IMMGEN_3_COLORS,
               ncol=1) + labs(title=paste(infectGroupName, timeName))
    }, error=\(e)c())
  })
}) %>% unlist(recursive=F) %>% 
  compact()


wrap_plots(fig3bParts) +
  plot_layout(guides = "collect",
              design = "ABCDE
                        ABCDE
                        ABCDE
                        ABC##")
ggsave(filename=file.path(plotDir, "fig3_b.pdf"), 
       width=16, height=8)


clustByCondUmap <- DimPlot(allCells, 
        label=F, 
        repel=T, 
        label.size=2.5,
        group.by="condition", 
        split.by = analysisMarkerCol, 
        ncol = 4)
ggsave(plot=clustByCondUmap, filename=file.path(plotDir, sprintf("clusterByCondition_umap_nolabel_by_%s.pdf", analysisMarkerCol)), width=16, height=16)

####Altered 3B#####
# Reusable theme: no axes, no titles, no facet strip labels
umap_theme_clean <- theme(
  axis.title = element_blank(),
  axis.text  = element_blank(),
  axis.ticks = element_blank(),
  axis.line  = element_blank(),
  panel.grid = element_blank(),
  plot.title = element_blank(),
  strip.text = element_blank(),
  strip.background = element_blank()
)
plotDir <- "~/ondemand/scRNAseq/FD/scRNAseq/plots"
#### Fig 3 B - UMAP per condition ####
infectGroups <- list(`Mtb-`=c("inf", "uninf", "neg"), `Mtb+`="pos")
timeGroups   <- list(`Preinf (D0)`=0, `D15`=15, `D28`=28)

fig3bParts <- imap(infectGroups, \(infectGroups, infectGroupName){
  imap(timeGroups, \(time, timeName){
    tryCatch({
      cells <- allCells[, allCells$infect %in% infectGroups &
                          allCells$day == time &
                          !is.na(allCells$ImmGenGroup)]
      UMAPPlot(
        cells,
        split.by = "tissue",
        group.by = "ImmGenGroup",
        cols     = IMMGEN_3_COLORS,
        ncol     = 1
      ) +
        umap_theme_clean
      # If you had set a title with labs(title=...), the line above removes it.
    }, error = \(e) c())
  })
}) %>%
  unlist(recursive = FALSE) %>%
  compact()

wrap_plots(fig3bParts) +
  plot_layout(
    guides = "collect",
    design = "ABCDE
              ABCDE
              ABCDE
              ABC##"
  ) &
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),  # margin around each subplot
    panel.spacing = unit(1.5, "lines") # space between panels
  )

ggsave(filename = file.path(plotDir, "fig3_b.pdf"),
       width = 16, height = 8)

# Cluster-by-condition UMAPs (no axes or strip titles)
clustByCondUmap <- DimPlot(
  allCells,
  label       = FALSE,
  repel       = TRUE,
  label.size  = 2.5,
  group.by    = "condition",
  split.by    = analysisMarkerCol,
  ncol        = 4
) + umap_theme_clean

ggsave(
  plot      = clustByCondUmap,
  filename  = file.path(plotDir, sprintf("clusterByCondition_umap_nolabel_by_%s.pdf", analysisMarkerCol)),
  width     = 18, height = 18
)

#Alluvial plot for clusters
allCells@meta.data %>%
  select(one_of(clusteringNames)) %>%
  mutate(integrated_snn_res.0.5=factor(integrated_snn_res.0.5, levels=sortClusterFactorLevels(integrated_snn_res.0.5))) %>%
  group_by(across(everything())) %>%
    summarize(Freq=n()) %>%
    ungroup() %>%
    ggplot(aes(axis1=ImmGenGroup,
               axis2=integrated_snn_res.0.5,
               y=Freq)) +
    geom_stratum(color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=2.5) +
    scale_x_discrete(limits = c("ImmGenGroup",
                                "Seurat clusters"), expand = c(.05, .05)) +
    scale_fill_brewer(palette="Paired") +
    guides(fill="none") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    labs(title="Clustering label correspondance")
ggsave(filename=file.path(plotDir, "cluster_to_cluster_alluvial.pdf"), width=10, height=12)

allCellsClusterCounts <- clusteringNames %>%
  map(\(clustering){
    sampleCounts <- count(allCells@meta.data, sample, name = "nPerSample")
    allCells@meta.data %>% 
      count(sample, condition, replicate, tissue, infect, day, .data[[clustering]]) %>%
      left_join(sampleCounts, by="sample") %>%
      mutate(cellsPerThousand=round((n/nPerSample)*1000))
  })

clusterAbundancePlots <- map(clusteringNames, \(clustering){
  plt <- allCellsClusterCounts[[clustering]] %>%
    mutate(across(one_of(clustering), \(cl)ifelse(is.na(cl), "unk", cl))) %>%
    ggplot(aes(x=factor(replicate), y=cellsPerThousand, fill=.data[[clustering]])) +
    geom_col() +
    geom_text(aes(label=trimName(.data[[clustering]])), position = position_stack(vjust = 0.5), size=2) +
    facet_grid(~tissue+day+infect, scales="free", space = "free") +
    theme_bw() +
    theme(legend.position="bottom", 
          legend.title=element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_rect(fill="white"),
          strip.text=element_text(face="bold")) +
    labs(x=NULL, y="Cells per thousand", title=clustering)
  if (clustering == "ImmGenGroup") {
    plt <- plt + scale_fill_manual(values=IMMGEN_3_COLORS)
  }
  
  ggsave(plot=plt, filename=file.path(plotDir, sprintf("rename_abundance_colPlot_%s.pdf", clustering)), width=10, height=6)
  plt
})

vlnPltConds <- c(mln.pos.d15="darkblue", mln.neg.d15="lightblue",
                mln.pos.d28="darkgreen", mln.neg.d28="lightgreen")

(VlnPlot(allCells[,allCells$ImmGenGroup %in% c("resDC2", "iMO2", "iMO3", "MC") & 
                   allCells$condition %in% names(vlnPltConds)], 
        assay="SCT", split.by = "condition",
        c("Mertk", "Cd68", "Vcam1", "Vegfa","Nos2"), 
        group.by = "ImmGenGroup", ncol = 3, raster = T) &
  theme_bw() &
    scale_fill_manual(values = vlnPltConds)) +
  plot_layout(guides="collect")
ggsave("plots/selected_vln_mc_imo_strat.pdf", width=8, height=6)

# DEGs in mtb+ vs mtb- ####
# Within tissues and timepoints

degTissues <- c("lung", "mln")
degTimes <- c(15, 28)

mtbInfDegs <- getMtbDegs(allCells,
                         mtbTissues=degTissues,
                         mtbTimes=degTimes,
                         cellTypeCol = analysisMarkerCol,
                         saveFile=sprintf("mtbInf_celltype_degs_%s.csv", analysisMarkerCol))
mtbInfDegGsea <- degGseaEnrichment(mtbInfDegs, 
                                 groupingFactors = c("tissue", "time", "celltype"),
                                 saveFile=sprintf("mtbinf_gsea_%s.rds", analysisMarkerCol))


mtbInfDegs %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC)>1 & celltype!="unk") %>%
  mutate(direction=ifelse(avg_log2FC<0, "decreased", "increased"),
         time=paste("Day", time)) %>%
  count(tissue, time, celltype, direction) %>%
  ggplot(aes(x=celltype, y=n, fill=direction)) +
    geom_col() +
    facet_grid(time~tissue) +
    scale_fill_manual(values=c(decreased="dodgerblue", increased="orangered")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          strip.background = element_rect(fill="white"),
          strip.text=element_text(face="bold")) +
    labs(title="DEG counts: Mtb+ vs Mtb- cells", x=NULL, y="N. genes") 
ggsave(sprintf("plots/mtbPosVsNeg_DEGcounts_%s.pdf", analysisMarkerCol), width=6, height=6)

mtbInfSigPathways <- mtbInfDegGsea %>%
  mutate(fdr = p.adjust(pval, method="fdr")) %>%
  filter(fdr < 0.05) %>% 
  ungroup() %>%
  slice_min(n=30, order_by=pval) %>%
  {unique(.$pathway)}


mtbInfDegGsea %>%
  filter(pathway %in% mtbInfSigPathways & !celltype %in% c("AM", "unk")) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(tissue, time, celltype) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
    geom_col(colour="black") +
    facet_grid(celltype~tissue+time, scales="free_y", space="free") +
    scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"),
          strip.text=element_text(face="bold"),
          axis.text.y=element_text(size=7)) +
  labs(title="GSEA: Mtb+ vs Mtb- cells", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/mtbPosVsNeg_GSEAres_%s.pdf", analysisMarkerCol), width=12, height=12)


mtbInfDegGsea %>%
  filter(pathway %in% selectedGoTermNames$gs_name & !celltype %in% c("AM", "unk")) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(tissue, time, celltype) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
  geom_col(colour="black") +
  facet_grid(celltype~tissue+time, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=7)) +
  labs(title="GSEA: Mtb+ vs Mtb- cells", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/mtbPosVsNeg_GSEA_selePathway_res_%s.pdf", analysisMarkerCol), width=12, height=12)




mtbPlotComps <- expand_grid(plotTissue=degTissues, 
                         plotTime=degTimes, 
                         plotCelltype=unique(allCells$lung_mouse_mtb_immune_celltypist_model)) %>%
  mutate(plotLabel=paste(plotCelltype, plotTissue, plotTime, sep="_"))

mtbInfPlots <- pmap(mtbPlotComps, \(plotTissue, plotTime, plotCelltype, plotLabel){
  plotCells <- tryCatch({
    allCells[,allCells$day==plotTime & allCells$tissue==plotTissue & allCells$lung_mouse_mtb_immune_celltypist_model==plotCelltype]
    }, error=\(e)NULL)
  markerGenes <- mtbInfDegs %>% 
    filter(time==plotTime & tissue==plotTissue & celltype==plotCelltype & genes %in% rownames(allCells)) %>%
    slice_max(order_by=avg_log2FC, n=12) %>%
    {.$genes}
  message(plotLabel)
  message(paste(markerGenes, collapse=", "))
  if (is.null(plotCells) || length(markerGenes)==0) {
    return(NA)
  }
  VlnPlot(plotCells, markerGenes, group.by="infect") &
    plot_annotation(title=sprintf("%s Mtb+ vs Mtb-: %s day %s ", plotCelltype, plotTissue, plotTime))
}) %>% set_names(mtbPlotComps$plotLabel) %>% 
  discard(\(obj)!"ggplot" %in% class(obj))

# DEGs in lung vs mln for Mo/Macs (pos and neg) ####

tissueDegs <- getTissueDegs(allCells,
                            compTissues=c("lung", "mln"),
                            times=c(15, 28),
                            infectStatuses = c("pos", "neg"),
                            cellTypeCol = analysisMarkerCol,
                            saveFile = sprintf("tissue_celltype_degs_%s.csv", analysisMarkerCol)) %>%
  mutate(plotCellType=IMMGEN_3_RENAMES[celltype])

#TODO merge infection statuses
bloodMlnDegs <- getTissueDegs(allCells,
                            compTissues=c("blood", "mln"),
                            times=c(15, 28),
                            infectStatuses = NA,
                            cellTypeCol = analysisMarkerCol,
                            saveFile = sprintf("analysis/tissue_blood_mln_celltype_degs_%s.csv", analysisMarkerCol)) %>%
  mutate(plotCellType=IMMGEN_3_RENAMES[celltype])

bloodLungDegs <- getTissueDegs(allCells,
                              compTissues=c("blood", "lung"),
                              times=c(15, 28),
                              infectStatuses = NA,
                              cellTypeCol = analysisMarkerCol,
                              saveFile = sprintf("analysis/tissue_blood_lung_celltype_degs_%s.csv", analysisMarkerCol)) %>%
  mutate(plotCellType=IMMGEN_3_RENAMES[celltype])

tissueDegGsea <- degGseaEnrichment(tissueDegs, 
                                   groupingFactors = c("time", "celltype", "infect"),
                                   saveFile=sprintf("analysis/tissue_gsea_%s.rds", analysisMarkerCol))

tissueDegs %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC)>0.5 & !celltype %in% c("AM", "Basophil", "unk")) %>%
  mutate(direction=ifelse(avg_log2FC<0, "decreased", "increased"),
         time=paste("Day", time),
         infect=ifelse(infect=="pos", "Mtb+", "Mtb-")) %>%
  count(time, plotCellType, direction, infect) %>%
  ggplot(aes(x=plotCellType, y=n, fill=direction)) +
  geom_col() +
  facet_grid(time~infect) +
  scale_fill_manual(values=c(decreased="dodgerblue", increased="orangered")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold")) +
  labs(title="DEG counts: Lung vs mLN", x=NULL, y="N. genes") 
ggsave(sprintf("plots/lungVsMLN_DEGcounts_%s.pdf", analysisMarkerCol), width=6, height=6)

##SubsetTissueDEGs
tissueDegs %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC)>0.5 & !infect %in% c("neg") & !celltype %in% c("AM", "Basophil", "unk","migDC1","migDC2","pMO","resDC1","resDC2a","resDC2b")) %>%
  mutate(direction=ifelse(avg_log2FC<0, "decreased", "increased"),
         time=paste("Day", time),
         infect=ifelse(infect=="pos", "Mtb+", "Mtb-")) %>%
  count(time, plotCellType, direction, infect) %>%
  ggplot(aes(x=plotCellType, y=n, fill=direction)) +
  geom_col() +
  facet_grid(time~infect) +
  scale_fill_manual(values=c(decreased="dodgerblue", increased="orangered")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold")) +
  labs(title="DEG counts: Lung vs medLN", x=NULL, y="Number of Genes") 
ggsave(sprintf("plots/lungVsMLN_DEGcounts_%s.pdf", analysisMarkerCol), width=6, height=6)

bloodMlnDegs %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC)>1 & !celltype %in% c("AM", "Basophil", "unk")) %>%
  mutate(direction=ifelse(avg_log2FC<0, "decreased", "increased"),
         time=paste("Day", time),
         infect=ifelse(infect=="pos", "Mtb+", "Mtb-")) %>%
  count(time, plotCellType, direction, infect) %>%
  ggplot(aes(x=plotCellType, y=n, fill=direction)) +
  geom_col() +
  facet_grid(cols = vars(time)) +
  scale_fill_manual(values=c(decreased="dodgerblue", increased="orangered")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold")) +
  labs(title="DEG counts: Blood vs mLN", x=NULL, y="N. genes") 
ggsave(sprintf("plots/bloodVsMLN_DEGcounts_%s.pdf", analysisMarkerCol), width=6, height=4)

bloodLungDegs %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC)>1 & !celltype %in% c("AM", "Basophil", "unk")) %>%
  mutate(direction=ifelse(avg_log2FC<0, "decreased", "increased"),
         time=paste("Day", time),
         infect=ifelse(infect=="pos", "Mtb+", "Mtb-")) %>%
  count(time, plotCellType, direction, infect) %>%
  ggplot(aes(x=plotCellType, y=n, fill=direction)) +
  geom_col() +
  facet_grid(cols = vars(time)) +
  scale_fill_manual(values=c(decreased="dodgerblue", increased="orangered")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold")) +
  labs(title="DEG counts: Blood vs Lung", x=NULL, y="N. genes") 
ggsave(sprintf("plots/bloodVsLung_DEGcounts_%s.pdf", analysisMarkerCol), width=6, height=4)


tissueDegs %>%
  filter(celltype %in% c("MDM", "Monocyte2a", "Monocyte2b", "DC1", "Neutrophil")) %>%
  ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(alpha=0.5) +
  geom_point(data=\(df)slice_max(group_by(df, infect, time, celltype), n=20, order_by=abs(avg_log2FC)), colour="red", shape=21) +
  geom_text_repel(data=\(df)slice_max(group_by(df, infect, time, celltype), n=20, order_by=abs(avg_log2FC)), aes(label=genes), size=2) +
  facet_grid(celltype~time+infect, scales="free_y") +
  scale_y_continuous(limits=c(0, 50), oob=scales::squish) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold")) +
  labs(title="Tissue DEGs: Lung vs mLN")
ggsave(sprintf("plots/lungVsMLN_DEGvolcano_%s.pdf", analysisMarkerCol), width=10, height=12)

topLungMlnMOdegs <- tissueDegs %>% 
  filter(tissues=="lung_mln" & plotCellType %in% c("MC", 'iMO2') & p_val_adj<0.05) %>% 
  group_by(time) %>% group_by(time) %>% slice_min(p_val, n=25) %>% 
  arrange(time, plotCellType, p_val)
write_csv(topLungMlnMOdegs, "plots/topLungMlnModegs.csv")

tissueSigPathways <- tissueDegGsea %>%
  filter(celltype %in% c("MDM", "Monocyte2a","Monocyte2b", "DC1", "Neutrophil")) %>%
  mutate(fdr = p.adjust(pval, method="fdr")) %>%
  filter(fdr < 0.05) %>% 
  ungroup() %>%
  slice_min(n=30, order_by=pval) %>%
  {unique(.$pathway)}

tissueDegGsea %>%
  filter(pathway %in% mtbInfSigPathways & celltype %in% c("MDM", "Monocyte2a", "Monocyte2b","DC1", "Neutrophil")) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(time, celltype, infect) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
  geom_col(colour="black") +
  facet_grid(celltype~infect+time, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=7)) +
  labs(title="GSEA: Lung vs mLN cells", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/lungVsMLN_GSEAres_%s.pdf", analysisMarkerCol), width=16, height=12)

tissueDegGsea %>%
  filter(pathway %in% selectedGoTermNames$gs_name & celltype %in% c("MDM", "Monocyte2a","Monocyte2b", "DC1", "Neutrophil")) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(time, celltype, infect) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
  geom_col(colour="black") +
  facet_grid(celltype~infect+time, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=7)) +
  labs(title="GSEA: Lung vs mLN cells", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/lungVsMLN_GSEAres_selePathways_%s.pdf", analysisMarkerCol), width=16, height=12)

# ggplot(plotDF, aes(x=celltype, y=pathway, colour=NES, size=-log10(padj))) +
#   geom_point() +
#   facet_grid(~tissue+time+infect, scales="free_x", space="free") +
#   scale_colour_distiller(palette="RdBu", limits=c(-3, 3), oob=scales::squish) +
#   scale_size_continuous(range=pointSizeRange) +
#   theme_minimal() +
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
#         panel.background = element_rect(fill="grey20"),
#         panel.grid = element_blank()) +
#   scale_y_discrete(limits=pathwayOrder) +
#   labs(x=NULL, y=NULL, caption=sprintf("Showing FDR < %.2f", fdrThresh))


# DEGs in D0 vs D15/28 ####

timeComps <- list(d15d0=c(t1=15, t2=0), d28d0=c(t1=28, t2=0), d28d15=c(t1=28, t2=15))
infectGroups <- list(nonpos=c("neg", "uninf", "inf"), pos=c("uninf", "pos"))

timeCompNameString <- c(d15d0="D15 vs D0", d28d0="D28 vs D0", d28d15="D28 vs D15")
infectCompNameString <- c(nonpos="Uninfected", pos="Mtb+ vs pre-inf")

getTimeAllDegs <- function(scObj, 
                           markerCol=analysisMarkerCol,
                           degTimeComps=timeComps, 
                           degInfectGroups=infectGroups,
                           degTissues=c("lung", "mln", "blood")) {
  imap_dfr(degTimeComps, \(timeComp, timeCompName){
    imap_dfr(degInfectGroups, \(infectGroup, infectGroupName){
      getTimeDegs(scObj,
                  tissues=degTissues,
                  compTimes=timeComp,
                  infectStatuses = infectGroup,
                  cellTypeCol = markerCol,
                  saveFile = sprintf("~/ondemand/scRNAseq/FD/scRNAseq/analysis/%s_celltype_%s_degs_%s_2.csv", timeCompName, infectGroupName, markerCol)) %>%
        mutate(comparison=sprintf("%s, %s", timeCompNameString[timeCompName], infectCompNameString[infectGroupName]))
    })
  })
}

timeDegsOrig <- list.files("~/ondemand/scRNAseq/FD/scRNAseq/analysis", "d28d0.*degs_immgen_refined.csv", full.names = T) %>%
  map_dfr(read_csv)

timeAllDegs <- getTimeAllDegs(allCells)

timeAllDegs_ImmGenRefined <- getTimeAllDegs(allCells, markerCol="ImmGenGroup_Refined")


timeDegGseaOrig <- read_rds("~/ondemand/scRNAseq/FD/scRNAseq/analysis/time_gsea.rds")

timeDegGsea <- degGseaEnrichment(timeAllDegs, 
                                 groupingFactors = c("tissue", "celltype", "comparison"),
                                 saveFile=sprintf("~/ondemand/scRNAseq/FD/scRNAseq/analysis/time_gsea_%s.rds", analysisMarkerCol))

timeDegGsea_ImmGenRefined <- degGseaEnrichment(timeAllDegs_ImmGenRefined, 
                                 groupingFactors = c("tissue", "celltype", "comparison"),
                                 saveFile=sprintf("~/ondemand/scRNAseq/FD/scRNAseq/analysis/time_gsea_%s.rds", "ImmGenGroup_Refined"))


extractSortedDays <- function(cmp){
  d1 <- gsub("^(D\\d*).+", "\\1", cmp)
  d2 <- gsub(".+ vs (D\\d*).+", "\\1", cmp)
  map2_chr(d1, d2, \(day1, day2){
    if(day1>day2) paste(day1, day2, sep="_") else paste(day2, day1, sep="_")
  })
}

dayDirection <- function(cmp){
  d1 <- gsub("^(D\\d*).+", "\\1", cmp)
  d2 <- gsub(".+ vs (D\\d*).+", "\\1", cmp)
  map2_int(d1, d2, \(day1, day2){
    if(day1>day2) 1 else -1
  })
}

testCor <- \(v1, v2)tryCatch(cor(v1, v2, use="complete.obs"), error=\(e)NA)

gseaSharedCols <- c("tissue", "celltype", "comparison", "pathway")
bind_rows(mutate(timeDegGseaOrig, version="orig"),
          mutate(timeDegGsea, version=analysisMarkerCol), 
          mutate(timeDegGsea_ImmGenRefined, version="ImmGenGroup_refined")) %>%
  ungroup() %>%
  mutate(status=str_split_i(comparison, ", ", 2),
         days=extractSortedDays(comparison),
         dayDir=dayDirection(comparison)) %>%
  select(tissue, celltype, pathway, NES, version, status, days) %>%
  pivot_wider(names_from="version", values_from="NES") %>%
  filter(grepl("Interferon", pathway, ignore.case=T)) %>%
  # group_by(tissue, status, days, celltype) %>%
  # summarize(orig_v_immgen=testCor(orig, ImmGenGroup_refined),
  #           orig_v_immgenv2=testCor(orig, ImmGenGroup_Refinedv2),
  #           immgen_v_immgenv2=testCor(ImmGenGroup_refined, ImmGenGroup_Refinedv2)) %>%
  filter(days=="D28_D15" & tissue=="mln")



timeDegHallmarkGsea <- timeDegGsea %>% 
  filter(grepl("^HALLMARK", pathway)) %>% 
  mutate(fdr=p.adjust(pval, "fdr"))

timeAllDegs %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC)>1 & celltype != "unk") %>%
  mutate(direction=ifelse(avg_log2FC<0, "decreased", "increased")) %>%
  count(tissue, celltype, direction, comparison) %>%
  ggplot(aes(x=celltype, y=n, fill=direction)) +
  geom_col() +
  facet_grid(comparison~tissue) +
  scale_fill_manual(values=c(decreased="dodgerblue", increased="orangered")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold")) +
  labs(title="DEG counts: Time", x=NULL, y="N. genes") 
ggsave(sprintf("~/ondemand/scRNAseq/FD/scRNAseq/plots/time_DEGcounts_%s.pdf", analysisMarkerCol), width=8, height=12)

#MDM/Monocyte volcano plot
timeAllDegs %>%
  filter(celltype %in% c("MDM", "Monocyte1", "Monocyte2a", "Neutrophil", "DC1", 'DC2') & tissue!="blood") %>%
  mutate(comparison=gsub("pre-infection", "preInf", comparison)) %>%
  ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj))) +
    geom_hex(bins=50) +
    geom_point(data=\(df)slice_max(group_by(df, celltype, tissue), n=10, order_by=abs(avg_log2FC)), colour="red", shape=21) +
    geom_text_repel(data=\(df)slice_max(group_by(df, celltype, tissue), n=10, order_by=abs(avg_log2FC)), aes(label=genes), size=2.5) +
    facet_grid(celltype~tissue+comparison, scale="free_y") +
    theme_bw() +
    scale_y_continuous(limits=c(0, 50), oob=scales::squish) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          strip.background = element_rect(fill="white"),
          strip.text=element_text(face="bold")) +
  labs(title="Time DEGs")
ggsave(sprintf("~/ondemand/scRNAseq/FD/scRNAseq/plots/time_DEGvolcano_%s.pdf", analysisMarkerCol), width=24, height=14)


timeSigPathways <- timeDegGsea %>%
  mutate(fdr = p.adjust(pval, method="fdr")) %>%
  filter(fdr < 0.05) %>%
  ungroup() %>%
  slice_min(n=40, order_by=pval) %>%
  {unique(.$pathway)}


timeDegGsea %>%
  mutate(comparison=gsub("pre-infection", "preInf", comparison)) %>%
  filter(pathway %in% timeSigPathways & 
           celltype %in% c("MDM", "Monocyte1", "Monocyte2a", "Monocyte2b","DC1", "DC2", "Neutrophil") &
           !grepl("D28 vs D15", comparison)) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(comparison, celltype, tissue) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
  geom_col(colour="black") +
  facet_grid(celltype~tissue+comparison, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=6)) +
  labs(title="GSEA: Day 15/28 vs Day 0", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/time_GSEAres_%s.pdf", analysisMarkerCol), width=24, height=10)

timeDegGsea %>%
  mutate(comparison=gsub("pre-infection", "preInf", comparison)) %>%
  filter(pathway %in% selectedGoTermNames$gs_name & 
           celltype %in% c("MDM", "Monocyte1", "Monocyte2a", "Monocyte2b", "DC1", "DC2", "Neutrophil") &
           !grepl("D28 vs D15", comparison)) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(comparison, celltype, tissue) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
  geom_col(colour="black") +
  facet_grid(celltype~tissue+comparison, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=6)) +
  labs(title="GSEA: Day 15/28 vs Day 0", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/time_GSEAres_selePathways_%s.pdf", analysisMarkerCol), width=24, height=16)


timeDegGsea %>%
  mutate(comparison=gsub("pre-infection", "preInf", comparison)) %>%
  filter(pathway %in% timeSigPathways & 
           celltype %in% c("MDM", "Monocyte1", "Monocyte2", "DC1", "DC2", "Neutrophil") &
           grepl("D28 vs D15", comparison)) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(comparison, celltype, tissue) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
  geom_col(colour="black") +
  facet_grid(celltype~tissue+comparison, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=6)) +
  labs(title="GSEA: Day 28 vs Day 15", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/time_D15_vs_28_GSEAres_%s.pdf", analysisMarkerCol), width=13, height=10)

timeDegGsea %>%
  mutate(comparison=gsub("pre-infection", "preInf", comparison)) %>%
  filter(pathway %in% selectedGoTermNames$gs_name & 
           celltype %in% c("MDM", "Monocyte1", "Monocyte2", "DC1", "DC2", "Neutrophil") &
           grepl("D28 vs D15", comparison)) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(comparison, celltype, tissue) %>%
  ggplot(aes(y=pathway, x=NES, fill=strSignif(padj))) +
  geom_col(colour="black") +
  facet_grid(celltype~tissue+comparison, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="orangered", `**`="orange", `*`="yellow", "n.s."="grey90")) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=6)) +
  labs(title="GSEA: Day 28 vs Day 15", x="GSEA Normalized enrichment score", y=NULL)
ggsave(sprintf("plots/time_D15_vs_28_GSEAres_selePathways_%s.pdf", analysisMarkerCol), width=12, height=12)

## Manual GSEA results to plot ####

manualGseaResFile <- "Meetings/GSEA_toPlot/2025_02_WG_MedLN_neg_top_GOBP.xlsx"


manualGseaRes <- map_dfr(excel_sheets(manualGseaResFile), \(sheet){
  read_excel(manualGseaResFile, sheet) %>%
    mutate(comparison=sheet)
}) %>%
  mutate(`P Value`=as.numeric(gsub("^<", "", `P Value`)),
         FDR=as.numeric(gsub("^<", "", FDR)),
         celltype=str_split_i(comparison, "-", 1),
         compTime=paste("Day", str_split_i(comparison, "-", 2)),
         goTerm=fct_reorder(paste(`Gene Set`, Description), NES))

ggplot(manualGseaRes, aes(y=goTerm, x=celltype, colour=NES, size=-log10(FDR))) +
  geom_point() +
  facet_grid(~compTime) +
  scale_colour_distiller(palette="RdBu", limits=c(-3, 3), oob=scales::squish) +
  scale_size_continuous(range=c(2, 5)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill="grey20"),
        panel.grid = element_blank(),
        strip.text = element_text(face="bold", size=12)) +
    labs(x=NULL, y=NULL, title='mLN Mtb-')
ggsave("Meetings/GSEA_toPlot/2025_02_WG_MedLN_neg_top_GOBP.pdf", width=7, height=12)

# ggplot(plotDF, aes(x=celltype, y=pathway, colour=NES, size=-log10(padj))) +
#   geom_point() +
#   facet_grid(~tissue+time+infect, scales="free_x", space="free") +
#   scale_colour_distiller(palette="RdBu", limits=c(-3, 3), oob=scales::squish) +
#   scale_size_continuous(range=pointSizeRange) +
#   theme_minimal() +
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
#         panel.background = element_rect(fill="grey20"),
#         panel.grid = element_blank()) +
#   scale_y_discrete(limits=pathwayOrder) +
#   labs(x=NULL, y=NULL, caption=sprintf("Showing FDR < %.2f", fdrThresh))

# Monocyte DEGs ####


# For D15 and D28 Mtb+ medLN and Lung:
#   -iMO2 vs MC cell DEGs
#   -iMO1 vs MC cell DEGs
# volcano plots with the top 50 DEGs labeled (please show all DEGs on plot),
# GO biological processes pathway analysis on DEGs with cutoffs for log2FC greater than 0.5 or less than -0.5, and adj-p val<0.05?

monoComps <- expand.grid(tissue=c("mln", "lung"), 
                        timepoint=c(15, 28), 
                        infect="pos",
                        ct1=c("iMO2"),
                        ct2=c("MC"))

monoDegSaveFile <- "~/ondemand/scRNAseq/FD/scRNAseq/mc_imo_degs_mtbpos.csv"

monoDegs <- tryCatch({
  read_csv(monoDegSaveFile)
}, error=\(e){
  message("Could not read file: ", monoDegSaveFile)
  message(as.character(e))
  message("Recalculating...")
  degs <- pmap_dfr(monoComps, \(tissue, timepoint, infect, ct1, ct2){
    filtCells <- allCells %>%
      {.[,.$tissue==tissue & .$infect==infect & .$day==timepoint & .$ImmGenGroup_Refinedv2_rename %in% c(ct1, ct2)]}
    Idents(filtCells) <- "ImmGenGroup_Refinedv2_rename"
    FindMarkers(filtCells,
                ident.1=ct1,
                ident.2=ct2,
                assay="SCT",
                test.use="negbinom",
                min.pct=0.2,
                logfc.threshold = 0.5,
                recorrect_umi=FALSE) %>%
      data.frame(check.names = F) %>%
      rownames_to_column("genes") %>%
      mutate(ct1=ct1,
             ct2=ct2,
             tissue=tissue, 
             day=timepoint,
             infect=paste(infect, collapse="/"))
  })
  write_csv(degs, monoDegSaveFile)
  read_csv(monoDegSaveFile)
}) %>% mutate(comparison = paste(ct1, ct2, sep=" vs ")) 

sigDegs <- monoDegs %>%
  filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)

topDegs <- sigDegs %>%
  slice_min(by=c(comparison, tissue, day), order_by=p_val_adj, n=50)



ggplot(monoDegs, aes(x=avg_log2FC, y=-log10(p_val_adj), label=genes)) +
  geom_point() +
  geom_point(data=topDegs, colour="red") +
  geom_text_repel(data=topDegs, colour="red", max.overlaps=100, size=3) +
  scale_y_continuous(limits=c(0, 20), oob=scales::squish) +
  facet_grid(day+tissue~comparison) +
  theme_minimal() +
  theme(strip.text=element_text(size=12, face="bold"),
        plot.title=element_text(size=14, face="bold")) +
  labs(title="Monocyte subset DEGS: Mtb+")
ggsave("~/ondemand/scRNAseq/FD/scRNAseq/analysis/mc_imo_degs_mtbpos_volcano.pdf", width=10, height=20)

d28mlnMonoDegs <- monoDegs %>%
  filter(day %in% c(28, 15) & tissue=="mln" & comparison == "iMO2 vs MC") %>%
  group_by(day, comparison) %>%
  slice_min(p_val, n=40) %>%
  ungroup() %>%
  mutate(genes=fct_reorder(genes, avg_log2FC), 
         infect=paste("Mtb", infect),
         day=paste("Day", day))

ggplot(d28mlnMonoDegs, aes(x=comparison, y=genes, size=abs(pct.1-pct.2), colour=avg_log2FC)) +
  geom_point() +
  facet_grid(~tissue+day+infect, space="free", scales="free_x") +
  scale_colour_distiller(palette="RdBu", limits=c(-5, 5), oob=scales::squish) +
  #scale_size_continuous(range=pointSizeRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill="grey20"),
        panel.grid = element_blank()) +
  labs(x=NULL, y=NULL)
ggsave("~/ondemand/scRNAseq/FD/scRNAseq/plots/mc_imo2_degs_d15_d28_mtbpos_bubble.pdf", width=5, height=12)

# #Bubble plots: top 50 DEGs for D28 medLN Mtb-infected MCs vs iMO2s TODO

# elyaHmapd38Mono <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("28_15")),
#                                   unique(d28_15_mln_degs$gene),
#                                   c("Monocyte2b", "MDM"),
#                                   c("Monocyte2b", "MDM"), pointSizeRange = c(1,4))
# ggsave("plots/elya-hmaps/hmap2a_imo2_mc_d28_25_pos.pdf",  plot=elyaHmap2a, width=8, height=13)



monoDegGsea <- degGseaEnrichment(monoDegs, 
                                 groupingFactors = c("day", "tissue", "comparison"),
                                 saveFile=gsub(".csv$", "_gsea.rds", monoDegSaveFile))

write_csv(monoDegGsea, gsub(".csv$", "_gsea.csv", monoDegSaveFile))

sigMonoGseaGoPathways <- filter(monoDegGsea, startsWith(pathway, "GOBP") & padj<0.01 & abs(NES)>2) %>%
  {unique(.$pathway)}

monoDegGsea %>%
  filter(pathway %in% sigMonoGseaGoPathways) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, NES)) %>%
  group_by(day, tissue, comparison) %>%
  ggplot(aes(y=pathway, x=comparison, colour=NES, size=-log10(padj))) +
    geom_point() +
    facet_grid(~day+tissue, scales="free_y", space="free") +
    scale_colour_distiller(palette="RdBu", limits=c(-3, 3), oob=scales::squish) +
    scale_size_continuous(range=c(1,7)) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.background = element_rect(fill="grey20"),
          panel.grid = element_blank()) +
    labs(title="GSEA: Mtb+ monocytes", x="GSEA Normalized enrichment score", y=NULL)
ggsave("~/ondemand/scRNAseq/FD/scRNAseq/analysis/mc_imo_degs_mtbpos_gobp_gsea.pdf", width=10, height=12)


monoDegGsea %>%
  filter(pathway %in% sigMonoGseaGoPathways) %>%
  ungroup() %>%
  mutate(pathway=fct_reorder(pathway, padj, max),
         sig=case_when(padj<0.001~'***',
                       padj<0.01~"**",
                       padj<0.05~"*",
                       padj<0.2~".",
                       TRUE~'')) %>%
  group_by(day, tissue, comparison) %>%
  ggplot(aes(y=pathway, x=NES, fill=sig)) +
  geom_col() +
  facet_grid(~day+tissue+comparison, scales="free_y", space="free") +
  scale_fill_manual(values=c(`***`="darkred", `**`="red", `*`="orange", `.`="yellow")) +
  scale_size_continuous(range=c(1,7)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.grid = element_blank()) +
  labs(title="GSEA: Mtb+ monocytes", x="GSEA Normalized enrichment score", y=NULL)
ggsave("~/ondemand/scRNAseq/FD/scRNAseq/analysis/mc_imo_degs_mtbpos_gobp_gsea_barplot.pdf", width=10, height=12)

#Use clusterProfiler to reannotate GSEA of these samples

keytypes(org.Mm.eg.db)


monoDegGseaCP <- enrichGO(gene=monoDegs,
                          OrgDb= org.Mm.eg.db,
                          keyType = "SYMBOL",
                          ont="BP",
                          universe=topDegs$genes)


s_monoDegGsea <- clusterProfiler::simplify(monoDegGsea)

# Save for REVIGO

mouseGOmap <- msigdbr(species="Mus musculus", category="C5", subcategory="GO:BP") %>%
  filter(gs_name %in% sigMonoGseaGoPathways) %>%
  select(gs_name, gs_exact_source) %>%
  distinct() %>%
  {set_names(.$gs_exact_source, .$gs_name)}

dir.create("analysis/REVIGO")
revigoMO1vsMC <- monoDegGsea %>%
  filter(pathway %in% sigMonoGseaGoPathways) %>%
  mutate(group=paste(tissue, "day", day)) %>%
  ungroup() %>%
  select(group, pathway, padj, NES, comparison) %>%
  mutate(goID=mouseGOmap[pathway]) %>%
  select(goID, padj, group) %>%
  split(.$group)

imap(revigoMO1vsMC, \(df, grp)write_csv(select(df, -group), sprintf("analysis/REVIGO/%s.csv", grp)))

# Save Seurat output ####

write_rds(allCells, file="analysis/annotated_seurat.rds")

# Heatmaps ####

elyaDegSpreadsheet <- "2024_12 DEG lists for heatmaps.xlsx"
elyaDegSheets <- excel_sheets(elyaDegSpreadsheet)
elyaDegs <- map_dfr(elyaDegSheets, \(s){
  read_excel(elyaDegSpreadsheet, sheet=s) %>% 
    pivot_longer(everything(), names_to="group", values_to="gene") %>%
    filter(!is.na(gene)) %>%
    mutate(celltype=s)
})

imo2D28_15_degs <- read_excel("F4I_medLN_Mtbinf_iMO2_MC_DEGs.xlsx", range = "A2:C36") %>%
  pivot_longer(everything(), names_to="group", values_to="gene") %>%
  filter(!is.na(gene)) %>%
  mutate(celltype='iMO2')
  
mcD28_15_degs <- read_excel("F4I_medLN_Mtbinf_iMO2_MC_DEGs.xlsx", range = "E2:G28") %>%
  pivot_longer(everything(), names_to="group", values_to="gene") %>%
  filter(!is.na(gene)) %>%
  mutate(celltype='MC')

d28_15_mln_degs <- bind_rows(imo2D28_15_degs, mcD28_15_degs)

#Elya hmap (1)
#A heatmap of only monocyte subset cells MO1, MO2a, MO2b, and MDM with 
#the genes listed on the MO2b- and MDM- sheets. 
#This is a list of top DEGs for D28 and D15 (vs D0) for MO2bs and for MDMs in the 
#medLN from mtb- fractions.
dir.create("plots/elya-hmaps")
elyaHmap1 <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("15_0", "28_0")), 
                                unique(filter(elyaDegs, celltype %in% c("MO2b-", "MDM-"))$gene), 
                                 c("Monocyte1", "Monocyte2a", "Monocyte2b", "MDM"), 
                                c("Monocyte1", "Monocyte2a", "Monocyte2b", "MDM"), pointSizeRange = c(1,4))
ggsave("plots/elya-hmaps/hmap1_mono_neg.pdf", plot=elyaHmap1, width=15, height=13)
#Elya hmap (2)
#A similar heatmap as (1) but with the MO2b+ and MDM+ gene lists in that same 
#file, only including monocyte subset cells in the heatmap.
elyaHmap2 <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("15_0", "28_0")), 
                                 unique(filter(elyaDegs, celltype %in% c("MO2b+", "MDM+"))$gene), 
                                 c("Monocyte1", "Monocyte2a", "Monocyte2b", "MDM"), 
                                 c("Monocyte1", "Monocyte2a", "Monocyte2b", "MDM"), pointSizeRange = c(1,4))
ggsave("plots/elya-hmaps/hmap2_mono_pos.pdf",  plot=elyaHmap2, width=15, height=13)


#Bubble plots: d28_vs_d25 TODO
#Bubble plots: top 50 DEGs for D28 medLN Mtb-infected MCs vs iMO2s TODO
# medLN vs lung analysis of Mtb-infected iMO2s at D15 and Mtb-infected MCs at D28?
elyaHmap2a <- clustDotGeneHeatmap(mutate(filter(tissueDegs, tissues %in% c("lung_mln") & infect=="pos"), times=time, tissue=tissues), 
                                 unique(d28_15_mln_degs$gene), 
                                 c("Monocyte2b", "MDM"), 
                                 c("Monocyte2b", "MDM"), pointSizeRange = c(1,4),
                                 inflabels=NA, timelabels = NA) +
  labs(caption=NULL)
ggsave("plots/elya-hmaps/hmap2a_imo2_mc_d28_25_pos.pdf",  plot=elyaHmap2a, width=4, height=13)

#Elya hmap (3)
#Aheatmap of only DC subset cells resDC1, resDC2a, resDC2b, migDC1, migDC2, 
#with the genes listed on the DCa- and DCb-. 
#This is a list of top DEGs between D28 and D15 (vs D0) for migDC and resDC2a subsets.
elyaHmap3 <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("15_0", "28_0")), 
                                 unique(filter(elyaDegs, celltype %in% c("DCa-", "DC2-"))$gene), 
                                 c("resDC1", "resDCa", "resDCb", "migDC1", "migDC2"), 
                                 c("resDC1", "resDCa", "resDCb", "migDC1", "migDC2"), pointSizeRange = c(1,4))
ggsave("plots/elya-hmaps/hmap3_DC_neg.pdf", plot=elyaHmap3, width=16, height=14)

#Elya hmap (4)
#A similar heatmap as (3) but with the DCa+ and DC2+ gene lists in that same file, 
# only including DC subset cells in the heatmap.
elyaHmap4 <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("15_0", "28_0")), 
                                 unique(filter(elyaDegs, celltype %in% c( "DC2+"))$gene), 
                                 c("resDC1", "resDCa", "resDCb", "migDC1", "migDC2"), 
                                 c("resDC1", "resDCa", "resDCb", "migDC1", "migDC2"), pointSizeRange = c(1,4))
ggsave("plots/elya-hmaps/hmap4_DC_pos.pdf", plot=elyaHmap4, width=16, height=12)

#Heatmaps
analysis1genes <- c("Tnf", "Vegfa", "Il12a", "Il12b", "Cxcl9", "Cxcl10", 
                    "Ly6c2", "Ccr2", "Fcgr1", "Tgfb1", "Tgfbr2", "Acod1", 
                    "Il1rn", "Il10", "Cd40", "Sirrpa", "Dpp4", "C5ar1", 
                    "Itgam", "Itgax", "Il1b", "Hif1a", "Apoe", "Mertk", 
                    "Ccl5", "Ccr5", "Ccr7", "Lyz2", "Nos2", "Il2ra", "Il2rb", 
                    "Pparg", "Ptges2", "Ifng", "Ifnb1", "Ifngr2", "Fabp5", "Cd1d1", 
                    "Fpr2", "Slc7a11", "Camk2a", "Alox5", "Alox5ap", "Ccl19", 
                    "Ceacam1", "Ciita", "Rfx5", "Rfxap", "Rfxank", 
                    "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-Eb2", "Cd74", "H2-DMa", 
                    "H2-DMb1", "H2-Oa", "H2-Ob")

# analysis1cellTypes <- c("Monocyte1", "Monocyte2a", "Monocyte2b", 
#                         "MDM", "DC1", "DC2", "DCa", "DCb")
analysis1cellTypes <- c("Monocyte1", "Monocyte2a", "Monocyte2b",
                        "MDM", "resDC1", "migDC1", "migDC2", "resDCa", "resDCb")

analysis1MtbposcellTypes <- c("Monocyte2a", "Monocyte2b", 
                              "MDM", "migDC2")


analysis2goPathways <- read_delim("hmap_go_terms.txt", 
                                  delim=": ",
                                  col_names = c("GOname", "GOid"))

analysis2goGseaNames <- msigdbr(species="Mus musculus", category="C5", subcategory="GO:BP") %>% 
  select(gs_name, gs_exact_source) %>% 
  distinct() %>% 
  filter(gs_exact_source %in% analysis2goPathways$GOid) 

analysis2goPathwaysAnnot <- analysis2goPathways %>%
  left_join(analysis2goGseaNames, by=c("GOid"="gs_exact_source"))

examinedGenes <- Reduce(union, 
                        list(analysis1genes,
                             unique(unlist(getMousePathways()[analysis2goGseaNames$gs_name])),
                             unique(unlist(getMousePathways()[grepl("^HALLMARK", names(getMousePathways()))]))))


unexaminedDEGs <- timeAllDegs %>% 
  filter(!genes %in% examinedGenes & 
           times %in% c("15_0", "28_0") &
           abs(avg_log2FC)>3 & abs(pct.1-pct.2)>0.33 & 
           p_val_adj<1e-3 & !grepl("^Gm", genes) & 
           !grepl("Rik\\d?$", genes)) %>%
  arrange(desc(abs(avg_log2FC)))

#D15/D28 vs D0 DEGs
analysis1hmap <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("15_0", "28_0")), 
                                     analysis1genes, 
                                     analysis1cellTypes, 
                                     analysis1MtbposcellTypes)

ggsave(plot=analysis1hmap, 
       filename=sprintf("analysis/timeDEG_hmap_%s.pdf", analysisMarkerCol), 
       width=14, height=10)

analysis1UnexaminedHmap <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("15_0", "28_0")), 
                                               unique(unexaminedDEGs$genes), 
                                               analysis1cellTypes, 
                                               analysis1MtbposcellTypes)
ggsave(plot=analysis1UnexaminedHmap, 
       filename=sprintf("analysis/timeDEG_hmap_noPathwayGenes_%s.pdf", analysisMarkerCol), 
       width=14, height=14)


analysis1d15d28hmap <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("28_15")), 
                                           analysis1genes, 
                                           analysis1cellTypes, 
                                           analysis1MtbposcellTypes)

ggsave(plot=analysis1d15d28hmap, 
       filename=sprintf("analysis/timeDEG_d28vs15_hmap_%s.pdf", analysisMarkerCol), 
       width=8, height=10)

analysis1UnexaminedD15D28Hmap <- clustDotGeneHeatmap(filter(timeAllDegs, times %in% c("28_15")), 
                                               unique(unexaminedDEGs$genes), 
                                               analysis1cellTypes, 
                                               analysis1MtbposcellTypes)
ggsave(plot=analysis1UnexaminedD15D28Hmap, 
       filename=sprintf("analysis/timeDEG_hmap_noPathwayGenes_d28d15_%s.pdf", analysisMarkerCol), 
       width=10, height=14)



analysis2hmap <- clustDotplotGseaHeatmap(filter(timeDegGsea, pathway %in% analysis2goPathwaysAnnot$gs_name & !grepl("D28 vs D15", comparison)),
                                         analysis1cellTypes,
                                         analysis1MtbposcellTypes)

ggsave(plot=analysis2hmap, 
       filename=sprintf("analysis/timeGSEA_hmap_%s.pdf", analysisMarkerCol), 
       width=15, height=8)


analysis2hmapD15_v_28 <- clustDotplotGseaHeatmap(filter(timeDegGsea, pathway %in% analysis2goPathwaysAnnot$gs_name & grepl("D28 vs D15", comparison)),
                                                 analysis1cellTypes,
                                                 analysis1MtbposcellTypes)

ggsave(plot=analysis2hmapD15_v_28, 
       filename=sprintf("analysis/timeGSEA_d15vs28_hmap_%s.pdf", analysisMarkerCol), 
       width=12, height=8)



analysis2HallmarkHmapOrig <- clustDotplotGseaHeatmap(filter(timeDegGseaOrig, grepl("^HALLMARK", pathway) & !grepl("D28 vs D15", comparison)),
                                                 analysis1cellTypes,
                                                 analysis1MtbposcellTypes, fdrThresh = 1, 
                                                 pointSizeRange = c(1,6))


analysis2HallmarkHmap_rename <- filter(timeDegGsea, grepl("^HALLMARK", pathway) & !grepl("D28 vs D15", comparison)) %>%
  mutate(celltype=IMMGEN_3_RENAMES[celltype]) %>%
  clustDotplotGseaHeatmap(IMMGEN_3_RENAMES[analysis1cellTypes],
                          IMMGEN_3_RENAMES[analysis1MtbposcellTypes], 
                          fdrThresh = 0.05, pointSizeRange=c(2,5))

analysis2HallmarkHmap <- clustDotplotGseaHeatmap(filter(timeDegGsea, grepl("^HALLMARK", pathway)),
                                                 analysis1cellTypes,
                                                 analysis1MtbposcellTypes)
ggsave(plot=analysis2HallmarkHmap, 
       filename=sprintf("analysis/timeGSEA_hallmark_hmap_%s.pdf", analysisMarkerCol), 
       width=20, height=8)

ggsave(plot=analysis2HallmarkHmap_rename, 
       filename=sprintf("analysis/timeGSEA_hallmark_hmap_%s_rename.pdf", analysisMarkerCol), 
       width=20, height=8)


analysis2GoHallmarkHmapOrig <- clustDotplotGseaHeatmap(filter(timeDegGseaOrig, (grepl("^HALLMARK", pathway) | pathway %in% analysis2goPathwaysAnnot$gs_name) & 
                                                            !grepl("D28 vs D15", comparison)),
                                                   analysis1cellTypes,
                                                   analysis1MtbposcellTypes, 
                                                   trimPathNames = F)

analysis2GoHallmarkHmap <- clustDotplotGseaHeatmap(filter(timeDegGsea, (grepl("^HALLMARK", pathway) | pathway %in% analysis2goPathwaysAnnot$gs_name)),
                                                   analysis1cellTypes,
                                                   analysis1MtbposcellTypes, 
                                                   trimPathNames = F)
ggsave(plot=analysis2GoHallmarkHmap, 
       filename=sprintf("analysis/timeGSEA_go_hallmark_hmap_%s.pdf", analysisMarkerCol), 
       width=24, height=12)

analysis2HallmarkD28vsD15Hmap <- clustDotplotGseaHeatmap(filter(timeDegGsea, grepl("^HALLMARK", pathway) & grepl("D28 vs D15", comparison)),
                                                         analysis1cellTypes,
                                                         analysis1MtbposcellTypes)

ggsave(plot=analysis2HallmarkD28vsD15Hmap, 
       filename=sprintf("analysis/timeGSEA_D15vsd28_hallmark_hmap_%s.pdf", analysisMarkerCol), 
       width=12, height=8)

analysis2GoHallmarkD28vs15Hmap <- clustDotplotGseaHeatmap(filter(timeDegGsea, (grepl("^HALLMARK", pathway) | pathway %in% analysis2goPathwaysAnnot$gs_name) & 
                                                            grepl("D28 vs D15", comparison)),
                                                   analysis1cellTypes,
                                                   analysis1MtbposcellTypes, 
                                                   trimPathNames = F)
ggsave(plot=analysis2GoHallmarkD28vs15Hmap, 
       filename=sprintf("analysis/timeGSEA_D15vsd28_go_hallmark_hmap_%s.pdf", analysisMarkerCol), 
       width=12, height=12)


#Signature analysis ####

#TODO update to analysisMarkerCol

#Also, must translate to mouse
cancerMonocyteSigs <- read_excel("Gene signature lists_TAM_M1_M2.xlsx",range = "Sheet1!A2:C56") %>%
  as.list() %>%
  map(\(x)humanGenesToMouse(as.character(na.omit(x)))) %>%
  set_names(gsub("[ -]", "_", names(.)))


cytokineDrivenCellTypeSigs <- read_excel("cell_cytokine_driver_dict.xlsx") %>%
  mutate(cellType=make.names(gsub(",$", "", `Polarization,`))) %>%
  split(.$cellType) %>%
  map(\(df)str_split_1(df$`Top marker genes`, pattern=fixed(":")))

moDccytokineDrivenCellTypeSigs <- cytokineDrivenCellTypeSigs %>%
  keep_at(\(nm)grepl("cDC|Mac|Mono", nm))

monoMacSigs <- c(cancerMonocyteSigs, moDccytokineDrivenCellTypeSigs)




ucellScoreFile <- "analysis/cancerMonocyteUcellScoreFile.csv"

allCells <- getUCellSigScores(allCells, monoMacSigs, ucellScoreFile)
signatureColumns <- c(paste0(names(monoMacSigs), "_UCell"))

oldIdents <- Idents(allCells)
Idents(allCells) <- "ImmGenGroup_Refined"
FeaturePlot(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="0"], 
            features =signatureColumns, 
            label=T, cols=c("grey", "red", "blue"), raster = T, split.by = "tissue")
ggsave("plots/m1_m2_monocyte_signatures_d0_UMAP.pdf", width=20, height=16)

FeaturePlot(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="15"], 
            features =signatureColumns, 
            label=T, cols=c("grey", "red", "blue"), raster = T, split.by = "tissue")
ggsave("plots/m1_m2_monocyte_signatures_d15_UMAP.pdf", width=20, height=16)

FeaturePlot(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="28"], 
            features =signatureColumns, 
            label=T, cols=c("grey", "red", "blue"), raster = T, split.by = "tissue")
ggsave("plots/m1_m2_monocyte_signatures_d28_UMAP.pdf", width=20, height=16)

# 
# FeaturePlot(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined)], 
#             features =c( "TAM_genes_UCell"), 
#             label=T, cols=c("grey90", "forestgreen"), raster = T)
# ggsave("plots/tam_monocyte_signatures_UMAP.pdf", width=8, height=6)
# 

plotSigViolinByTissue <- function(seuratObj, sigNames, pltSubtitle=NA, invertPattern=F) {
  (VlnPlot(seuratObj, features=sigNames, pt.size = 0, split.by = "tissue") & 
     theme_bw() & 
     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) &
     scale_y_continuous(limits=c(0, 0.4), oob=scales::squish) &
     labs(x=NULL, y="Signature score", subtitle=pltSubtitle)) +
     plot_layout(guides="collect", ncol=1) 
}


day0 <- plotSigViolinByTissue(allCells[,!grepl("^unk", allCells$ImmGenGroup_Refined) & allCells$day=="0"],
                                       sigNames=paste0(names(cancerMonocyteSigs), "_UCell"), 
                              pltSubtitle ="Day 0")

day15 <- plotSigViolinByTissue(allCells[,!grepl("^unk", allCells$ImmGenGroup_Refined) & allCells$day=="15"],
                               sigNames=paste0(names(cancerMonocyteSigs), "_UCell"), 
                               pltSubtitle ="Day 15")

day28 <- plotSigViolinByTissue(allCells[,!grepl("^unk", allCells$ImmGenGroup_Refined) & allCells$day=="28"],
                               sigNames=paste0(names(cancerMonocyteSigs), "_UCell"), 
                               pltSubtitle ="Day 28")

(day0 | day15 | day28) 
ggsave("plots/ucell_signature_violinplots.pdf", width=18, height=10)


day0Mono <- plotSigViolinByTissue(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="0"],
                              sigNames=paste0(names(cancerMonocyteSigs), "_UCell"), 
                              pltSubtitle ="Day 0")

day15Mono <- plotSigViolinByTissue(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="15"],
                               sigNames=paste0(names(cancerMonocyteSigs), "_UCell"), 
                               pltSubtitle ="Day 15")

day28Mono <- plotSigViolinByTissue(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="28"],
                               sigNames=paste0(names(cancerMonocyteSigs), "_UCell"), 
                               pltSubtitle ="Day 28")

(day0Mono | day15Mono | day28Mono) 
ggsave("plots/ucell_signature_violinplots_monoMdm.pdf", width=14, height=10)


groupedMacMonoDCSigs <- data.frame(sigs=names(moDccytokineDrivenCellTypeSigs)) %>% 
  mutate(group=str_split_i(sigs, '\\.', 1)) %>% 
  split(.$group) %>% 
  map(\(df)df$sigs)

dictSigPlots <- imap(groupedMacMonoDCSigs[3:4], \(sigs, grp){
  d0 <- plotSigViolinByTissue(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="0"],
                                        sigNames=paste0(sigs, "_UCell"), 
                                        pltSubtitle ="Day 0")
  
  d15 <- plotSigViolinByTissue(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="15"],
                                         sigNames=paste0(sigs, "_UCell"), 
                                         pltSubtitle ="Day 15")
  
  d28 <- plotSigViolinByTissue(allCells[,grepl("^(Mono|MDM)", allCells$ImmGenGroup_Refined) & allCells$day=="28"],
                                         sigNames=paste0(sigs, "_UCell"), 
                                         pltSubtitle ="Day 28")
  d0 | d15 | d28
  ggsave(sprintf("plots/ucell_signature_violinplots_cytokineResp_%s.pdf", grp), width=14, height=2.5*length(sigs))
})

dictDcSigPlots <- imap(groupedMacMonoDCSigs[1:2], \(sigs, grp){
  d0 <- plotSigViolinByTissue(allCells[,grepl("^(p?DC)", allCells$ImmGenGroup_Refined) & allCells$day=="0"],
                              sigNames=paste0(sigs, "_UCell"), 
                              pltSubtitle ="Day 0")
  
  d15 <- plotSigViolinByTissue(allCells[,grepl("^(p?DC)", allCells$ImmGenGroup_Refined) & allCells$day=="15"],
                               sigNames=paste0(sigs, "_UCell"), 
                               pltSubtitle ="Day 15")
  
  d28 <- plotSigViolinByTissue(allCells[,grepl("^(p?DC)", allCells$ImmGenGroup_Refined) & allCells$day=="28"],
                               sigNames=paste0(sigs, "_UCell"), 
                               pltSubtitle ="Day 28")
  d0 | d15 | d28
  ggsave(sprintf("plots/ucell_signature_violinplots_cytokineResp_%s.pdf", grp), width=14, height=2.5*length(sigs))
})




Idents(allCells) <- oldIdents



