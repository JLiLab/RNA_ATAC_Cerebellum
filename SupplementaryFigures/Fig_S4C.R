library(ArchR);library(Signac);library(Seurat)
library(cowplot)
library(patchwork)
set.seed(1234)

# Load ArchR object
threadNum = 5
ArchR_proj = "snATACseq/JL20001/ArchR/E12_ATAC"
proj <- loadArchRProject(path = ArchR_proj)

# keep Cb-only brain cells
id <- proj$cellNames[!proj$cellType %in% c("Endothelium","Pericytes","Microglia","Mes")]
proj1=proj[id, ]

# rotate UMAP for better visualization
proj1@embeddings$UMAP$df=proj1@embeddings$UMAP$df*-1
proj1 <- addImputeWeights(proj1) # repeat imputation

# select markers
markerGenes  <- c(
  "Atoh1","Cdon","Dkk3","Hes1","Lmx1a","Msx1","Reln","Rspo1","Wnt1","Wnt9a"
)

p2 <- plotEmbedding(
  ArchRProj = proj1, 
  colorBy = "GeneIntegrationMatrix", 
  continuousSet = "greyMagma",
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj1)
)

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    labs(x="",y="")+
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.line=element_blank(),
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),
      axis.text.x=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank()
    )
})

pdf("figures/UMAP_PTZ_markers.pdf", w =15, h =6.4)
do.call(plot_grid, c(list(ncol = 5), p2c))
dev.off()

pdf("figures/UMAP_PTZ_Legend.pdf", w =5, h =6.4)
plotEmbedding(
  ArchRProj = proj1, 
  colorBy = "GeneIntegrationMatrix", 
  continuousSet = "greyMagma",
  name = "Atoh1", 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj1)
)
dev.off()

# Plot RNA
# Load E12.5 scRNA-seq seurat object
cdS <- readRDS("E12/Cb_only.rds")

# Keep Cb-only brain cells
cdS1=subset(cdS, idents= c(23,25,26,18,27), invert = T)
DimPlot(cdS1, reduction = "ref.umap", label = T)+NoLegend()

library(SeuratWrappers)
cdS1=RunALRA(cdS1)

pn <- FeaturePlot(cdS1, markerGenes, combine = FALSE, reduction = "ref.umap",
                  cols = c("grey","#FB8861FF","#B63679FF","#51127CFF","#000004FF"), raster = T)

for(i in 1:length(pn)) {
  pn[[i]] <- pn[[i]] + NoLegend()+ NoAxes()
}

pdf("figures/UMAP_PTZ_markers_alra.pdf", w =15, h =6.4)
do.call(plot_grid, c(list(ncol = 5), pn))
dev.off()

