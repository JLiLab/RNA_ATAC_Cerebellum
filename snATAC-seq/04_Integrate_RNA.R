library(ArchR);library(Seurat);library(cowplot)
set.seed(1)
threadNum = 5
ArchR_proj = "E12_14"

addArchRThreads(threads = threadNum) 
proj <- loadArchRProject(path = ArchR_proj)

# Cross-platform linkage of scATAC-seq cells with scRNA-seq cells
cdS <- readRDS("E12_14/RNA/data/Seurat_mnn.rds")

seRNA <- as.SingleCellExperiment(cdS)
table(colData(seRNA)$cellType1)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cellType", embedding = "UMAP")

# unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "cellType1",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(proj$cellType, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
#      preClust              
# [1,] "NPC.g"    "NPC"      
# [2,] "PC.Pcp4"  "PC4"      
# [3,] "GCP"      "GCP"      
# [4,] "Endo"     "Pericyte" 
# [5,] "X.Tlx3"   "Isth.Tlx3"
# [6,] "NPC"      "NPC.e"    
# [7,] "PC1"      "PC1"      
# [8,] "CN.gaba"  "CN.gaba"  
# [9,] "PC5"      "PC5"      
# [10,] "NPC.g"    "NPC.g"    
# [11,] "PC.Pre"   "GABA.Pre" 
# [12,] "X.Lhx9"   "lCN"      
# [13,] "Isth"     "Mes.Isl1" 
# [14,] "lCN"      "Isth.Lhx9"
# [15,] "mCN"      "mCN"      
# [16,] "GC"       "GC"       
# [17,] "NPC.S"    "PTZ"      
# [18,] "PC3"      "PC3"      
# [19,] "PC5"      "PC.En1"   
# [20,] "Isth.Otp" "Isth.Otp" 
# [21,] "GABA.Pro" "GABA.Pro" 
# [22,] "X.Lhx9"   "lCN-like" 
# [23,] "Endo"     "PVM"      
# [24,] "IN"       "IN"       
# [25,] "Endo"     "VECA"     
# [26,] "PC5"      "PCx"      
# [27,] "PC3"      "PCx1"     
# [28,] "IN"       "IN.Pre"   
# [29,] "Endo"     "C2"       
# [30,] "Mes"      "MidO"  

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
p1+p2

plotPDF(plot_grid(p1,p2), name = "RNA-Integration_unconstrained.pdf", ArchRProj = proj, addDOC = FALSE, width = 8, height = 8)


# Constrained Integration
# group <- c("PC1","PC3","PC4","PC5","PC.En1","PCx","PCx1")
group <- c("CN.gaba")
clustNon_group <- setdiff(unique(proj$cellType),group)

# rna_group <- colnames(seRNA)[colData(seRNA)$cellType %in% c("PC.Pcp4","PC1","PC5","PC3","PC1","PC.Cdh9","PCx")]
rna_group <- colnames(seRNA)[colData(seRNA)$cellType %in% c("CN.gaba")]
rnaNon_group <- setdiff(colnames(seRNA),rna_group)

groupList <- SimpleList(
  Group = SimpleList(
    ATAC = proj$cellNames[proj$cellType %in% group],
    RNA = rna_group
  ),
  Non_Group = SimpleList(
    ATAC = proj$cellNames[proj$cellType %in% clustNon_group],
    RNA = rnaNon_group
  )    
)

# We pass this list to the `groupList` parameter of the `addGeneIntegrationMatrix()` function to constrain our integration. Note that, in this case, we are still not adding these results to the Arrow files (`addToArrow = FALSE`). We recommend checking the results of the integration thoroughly against your expectations prior to saving the results in the Arrow files. We illustrate this process in the next section of the book.
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "cellType1",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

# compare results
pal <- paletteDiscrete(values = proj$predictedGroup_Un)
p2 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal
)

p3 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "predictedGroup", 
  pal = pal
)

plotPDF(plot_grid(p1,p2,p3, ncol = 2), name = "RNA-Integration_comparison.pdf", ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)

getAvailableMatrices(proj)
  
markerGenes  <- c(
  "Foxp1","Foxp2","Sox14",
  "Atoh1", "Tlx3","Eomes","Isl1","Otx2",
  "Id3", "Lmx1a","Olig2",
  "Tfap2b", "Pax2", 
  "Ptf1a", "Wnt1"
)

p1 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p2 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

pdf(paste0(ArchR_proj,"/Plots/Markers_IntegratedMatrix.pdf"), w =9, h =16)
do.call(plot_grid, c(list(ncol = 3), p1c))
dev.off()

pdf(paste0(ArchR_proj,"/Plots/Markers_geneScoreMatrix.pdf"), w =9, h =16)
do.call(plot_grid, c(list(ncol = 3), p2c))
dev.off()

markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneIntegrationMatrix", 
  groupBy = "cellType",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

sig_genes <- data.frame(unlist(markerList))
sig_genes$cluster <- rownames(sig_genes)
openxlsx::write.xlsx(sig_genes, "results/marker_integrated.xlsx")

markerGenes  <- c(
  "Tlx3", 
  "Meis2", "Atoh1","Lhx9","Lhx2",
  "Pax2", "Kdr", "EBF1", "Tfap2b", #B-Cell Trajectory
  "Foxp1", "Foxp2", "Sox14",
  "Ascl1","Ptf1a", 
  "Hes5", "Notch1", "Nes", "Wls" #TCells
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

pdf(paste0(ArchR_proj,"/Plots/heatmap_IGS.pdf"), h =5, w =8)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveArchRProject(ArchRProj = proj, outputDirectory = ArchR_proj, load = FALSE)

