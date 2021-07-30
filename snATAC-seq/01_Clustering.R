library(ArchR);library(cowplot);
set.seed(123)
#set threads to one half of cores requested for job
addArchRThreads(threads = 20)
addArchRGenome("mm10")
inputFiles <- c("../JL20001/E12_fragments.tsv.gz","../JL20002/E13_fragments.tsv.gz","../JL20003/E14_fragments.tsv.gz")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = c("E12", "E13", "E14"),
  filterTSS = 5, #Dont set this too high because you can always increase later
  filterFrags = 1250,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = FALSE
)

# Inferring Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "E12_14",
  copyArrows = TRUE 
)

length(proj$cellNames)

# Filtering Doublets from an ArchRProject
proj <- filterDoublets(ArchRProj = proj)
length(proj$cellNames)

# Dimensionality Reduction and Clustering
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:35,
  force = TRUE
)  

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

# Visualizing in a 2D UMAP Embedding
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.3, 
  metric = "cosine",
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

plotPDF(p1,p2, name = "Plot-UMAP-Clusters_new.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 5)

# Find feature genes using gene score matrix
markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

sig_genes <- data.frame(unlist(markerList))
sig_genes$cluster <- rownames(sig_genes)
openxlsx::write.xlsx(sig_genes, "results/marker.xlsx")

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

pdf("E12_14/Plots/heatmapGS.pdf", h =5, w =8)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

proj <- addImputeWeights(proj)

markerGenes  <- c(
  "Foxp1","Foxp2","Sox14","Etv1",
  "Atoh1", "Tlx3","Wnt1","En1","Otx2",
  "Rorb", "Lmx1a","Olig2",
  "Tfap2b", "Pax2", 
  "Nr2f2", "Bcl11b"
)

p2 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

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

pdf("E12_14/Plots/Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", w =16, h =16)
do.call(plot_grid, c(list(ncol = 4), p2c))
dev.off()

# Save object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "E12_14",
  overwrite = TRUE,
  load = FALSE
)






