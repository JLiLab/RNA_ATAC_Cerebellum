source("/home/CAM/jali/TOOLS/helper_scripts/findUniquePeaks.R")
library(ArchR);library(Signac);library(cowplot)   

# Setting up
set.seed(123)
addArchRThreads(threads = 20) 
proj <- loadArchRProject(path = "E12_14")

# # Manully annotate the cell clusters in excel and add the cell type information to the ArchR object
# sig_genes <- openxlsx::read.xlsx("results/marker_annotated.xlsx")
# tab <- sig_genes[!duplicated(sig_genes$cluster),]
# clusters <- tab$cluster
# cellType <- tab$cellType
# proj$cellType <- proj$Clusters
# proj$cellType <- mapLabels(proj$Clusters, newLabels = cellType, oldLabels = clusters)
# 
# p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
# p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cellType", embedding = "UMAP")
# 
# pdf("E12_14/Plots/UMAP_cellType.pdf", h =5, w =10)
# plot_grid(p1,p2)
# dev.off()
# 
# # Save object
# saveArchRProject(ArchRProj = proj, outputDirectory = "E12_14", load = F)
# 
# Making Pseudo-bulk Replicates
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "cellType", force = FALSE)

# 10 Calling Peaks w/ Macs2
pathToMacs2 <- "/home/CAM/jali/anaconda3/envs/python27/bin/macs2"
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "cellType",
  pathToMacs2 = pathToMacs2,
  force = FALSE
)

proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

# Save object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "E12_14",
  overwrite = TRUE,
  load = FALSE
)

# # Find cluster-specific peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "cellType",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
daPeakList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

da <- daPeakList[[3]]
da$cellType <- names(daPeakList)[3]
for (i in 4:length(names(daPeakList))){
  res <- daPeakList[[i]]
  if(nrow(res) > 0){
    res$cellType <- names(daPeakList)[i]
    da <- rbind(da,res)
  }
}
table(da$cellType)
# GABA.Pre  GABA.Pro        GC       GCP        IN    IN.Pre Isth.Lhx9  Isth.Otp Isth.Tlx3       lCN  lCN-like       mCN  Mes.Isl1      MidO       NPC     NPC.e 
# 12443     12246     20281     14071      3917      4652      6825      4107     10220     18940        19     16363     11186      3344     18659     15192 
# NPC.g    PC.En1       PC1       PC3       PC4       PC5       PCx      PCx1  Pericyte       PTZ       PVM      VECA 
# 22616      8597      8222      6876      9943      8244      2271      5956     14443      9188      1346       155 

median(table(da$cellType))
# [1] 8802

da$peak <- paste0(da$seqnames,":",paste0(da$start,"-",paste0(da$end)))
nrow(da[!duplicated(da$peak),])
# [1] 162309
length(proj@peakSet)
# [1] 401835
# 162309/401835*100 = 40.39195


# saveRDS(markersPeaks, "results/markersPeaks.rds")
# 
# # 11.2 Plotting Marker Peaks in ArchR
# heatmapPeaks <- plotMarkerHeatmap(
#   seMarker = markersPeaks, 
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
#   transpose = TRUE
# )
# 
# # draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)
# 
# # using HOMER motif
# proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "homer", name = "Motif", force = TRUE)
# 
# # 12.2 Motif Enrichment in Marker Peaks
# enrichMotifs <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
# )
# 
# heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = F, cutOff = 5)
# plotPDF(heatmapEM, name = "Motifs-HOMER-Heatmap", width = 12, height = 10, ArchRProj = proj, addDOC = FALSE)
# 
# # using JASPAR
# proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", name = "Motif", force = TRUE)
# 
# # 12.2 Motif Enrichment in Marker Peaks
# enrichMotifs <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
# )
# 
# heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = F, cutOff = 5)
# plotPDF(heatmapEM, name = "Motifs-JASPAR2020-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)
# 
# proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = T)
# 
# # 12.2 Motif Enrichment in Marker Peaks
# enrichMotifs <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
# )
# enrichMotifs
# 
# heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = F, cutOff = 10)
# plotPDF(heatmapEM, name = "Motifs-cisbp-Heatmap", width = 6, height = 15, ArchRProj = proj, addDOC = FALSE)


# 13.1 Motif Deviations
if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

motifs <- c("Lhx1", "En1", "Foxp1", "Etv1", "Tfap2b", "Atoh1")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

proj <- addImputeWeights(proj)
p <- plotGroups(ArchRProj = proj, 
                groupBy = "cellType", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})

pdf("E12_14/Plots/Motiff_ChromVar.pdf")
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
dev.off()

proj <- addArchRAnnotations(ArchRProj = proj, collection = "CistromeTFBS")

# Save object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "E12_14",
  overwrite = TRUE,
  load = FALSE
)
