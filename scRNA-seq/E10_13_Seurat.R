library(Seurat);library(tidyverse);library(SeuratWrappers)
cdS <- readRDS('Seurat.rds')
cdS <- NormalizeData(cdS,verbose = F) %>% FindVariableFeatures(verbose = F)

cdS.mnn <- RunFastMNN(object.list = SplitObject(cdS, split.by = "dataset"))
cdS[["mnn"]] <- CreateDimReducObject(Embeddings(cdS.mnn, reduction = "mnn"), assay = "RNA")

dir.create("figures/UMAPtest")
for (n in seq(from =20, to = 90, by = 5)){
  cdS <- RunUMAP(cdS, reduction = "mnn", dims = 1:30, return.model = T, umap.method = 'uwot', n.neighbors = n, ,min.dist = 0.2)
  pdf(paste0("figures/UMAPtest/UMAP_", n, "L.pdf"))
  print(DimPlot(cdS, label = TRUE, group.by = "cellType")+NoLegend()+NoAxes())
  dev.off()
}

# inspection of clusters
dir.create("figures/clusters")
for (j in levels(Idents(cdS))){
  id <- WhichCells(cdS, idents = j)
  pdf(paste0("figures/clusters/UMAP_Ident", j, ".pdf"))
  print(DimPlot(cdS, cells.highlight = id))
  dev.off()
}

# Remove putative LQC, cluster 19
id <- WhichCells(cdS, idents = 19)

pdf("figures/LQC.pdf")
DimPlot(cdS, cells.highlight = id)+NoLegend()
dev.off()

cdS1 <- subset(cdS, idents = 19, invert = T)

dir.create("figures/UMAPtest1")
for (n in seq(from =20, to = 90, by = 5)){
  cdS1 <- RunUMAP(cdS1, reduction = "mnn", dims = 1:30, return.model = T, umap.method = 'uwot', n.neighbors = n, ,min.dist = 0.25)
  pdf(paste0("figures/UMAPtest1/UMAP_", n, "L.pdf"))
  print(DimPlot(cdS1, label = TRUE, group.by = "cellType1")+NoLegend()+NoAxes())
  dev.off()
}

dir.create("figures/UMAPtest2")
for (n in seq(from =20, to = 90, by = 5)){
  cdS1 <- RunUMAP(cdS1, reduction = "mnn", dims = 1:35, return.model = T, umap.method = 'uwot', n.neighbors = n, ,min.dist = 0.30)
  pdf(paste0("figures/UMAPtest2/UMAP_", n, "L.pdf"))
  print(DimPlot(cdS1, label = TRUE, group.by = "cellType1")+NoLegend()+NoAxes())
  dev.off()
}

dir.create("figures/UMAPtest3")
for (n in seq(from =45, to = 55, by = 2)){
  cdS1 <- RunUMAP(cdS1, reduction = "mnn", dims = 1:35, return.model = T, umap.method = 'uwot', n.neighbors = n, ,min.dist = 0.20)
  pdf(paste0("figures/UMAPtest3/UMAP_", n, "L.pdf"))
  print(DimPlot(cdS1, label = TRUE, group.by = "cellType1")+NoLegend()+NoAxes())
  dev.off()
}

dir.create("figures/UMAPtest4")
for (n in seq(from =30, to = 40, by = 2)){
  cdS1 <- RunUMAP(cdS1, reduction = "mnn", dims = 1:35, return.model = T, umap.method = 'uwot', n.neighbors = n, ,min.dist = 0.35)
  pdf(paste0("figures/UMAPtest4/UMAP_", n, "L.pdf"))
  print(DimPlot(cdS1, label = TRUE, group.by = "cellType1")+NoLegend()+NoAxes())
  dev.off()
}

pdf("./figures/UMP_stages.pdf", w =16, h =4.7)
DimPlot(cdS, group.by = "ident", split.by = 'stage', label = T)+NoLegend()
dev.off()

cdS1 <- RunUMAP(cdS1, reduction = "mnn", dims = 1:35, return.model = T, umap.method = 'uwot', n.neighbors = 35, ,min.dist = 0.35)
cdS1 <- FindNeighbors(cdS1, reduction = "mnn", dims = 1:35)
cdS1 <- FindClusters(cdS1, resolution = 1.4)

DimPlot(cdS1, label = T)+NoLegend()
DimPlot(cdS1, label = T, group.by = "cellType1")+NoLegend()

saveRDS(cdS1, "data/E10_13.rds")

test <- readRDS("/Volumes/Backup_JamesLi/Experiments1/scRNAseq/CerebellarAnlage/E12/Seurat4/data/Seurat1.rds")
DimPlot(test, label = T, group.by = "cellType")+NoLegend()
test@assays$RNA

anchors <- FindTransferAnchors(
  reference = cdS,
  query = test,
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  reference.reduction = "mnn",
  dims = 1:50
)

test.m <- MapQuery(
  anchorset = anchors,
  query = test,
  reference = cdS1,
  refdata = "seurat_clusters",
  reference.reduction = "pca", 
  reduction.model = "umap"
)

p1 <- DimPlot(test.m, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE)+NoAxes()+NoLegend()
p2 <- DimPlot(test.m, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)+NoAxes()+NoLegend()
p3 <- DimPlot(test.m, reduction = "ref.umap", group.by = "cellType", label = TRUE, repel = TRUE)+NoAxes()+NoLegend()

pdf("figures/UMAP_inhouse_E12.pdf", w =15, h =5)
p1 + p2 + p3
dev.off()


library(Seurat);
cdS <- readRDS("data/E10_13n.rds")

# Determine the cell cycle phases among the presumptive NPCs
cc.genes <- read.csv("/Volumes/jali/Genome/Cell Cycle Genes.csv", stringsAsFactors = F)
s.genes <- cc.genes$mGene[cc.genes$cellCycle == "G1/S"]
m.genes <- cc.genes$mGene[cc.genes$cellCycle == "G2/M"]
cdS <- CellCycleScoring(cdS, s.features = s.genes, g2m.features = m.genes)


pdf("figures/UMAP_cellType1.pdf", w =6, h=6.2)
DimPlot(cdS, group.by = "cellType1", label = T)+NoLegend()+NoAxes()
DimPlot(cdS, group.by = "stage")+NoLegend()+NoAxes()
DimPlot(cdS, group.by = "Phase")+NoLegend()+NoAxes()
dev.off()

cdS1 <- subset(cdS, idents = c(17,18,20,21,22), invert=T)
pdf("figures/UMAP_stages1.pdf", w =9, h=9.2)
DimPlot(cdS1, group.by = "cellType1", split.by = "stage", ncol = 2)+NoLegend()+NoAxes()
dev.off()


library(ggplot2)
DefaultAssay(cdS1) <- "SCT"

genes <- c("Ascl1","Neurog1","Ptf1a","Atoh1",
           "Lhx1","Tlx3","Kirrel2","Tfap2b")
plots <- list()
for (g in genes){
  plots[[g]] <- FeaturePlot(cdS1, g ,pt.size = 1)+
    NoAxes()+NoLegend()+
    ylim(-4,3)+xlim(0,6)
}

pdf("figures/markers_sct.pdf", w = 12, h = 6)
do.call(cowplot::plot_grid, c(list(ncol = 4),plots))
dev.off()

FeaturePlot(cdS, "Aldoc",pt.size = 1)+
  NoAxes()+NoLegend()+
  ylim(-4,3)+xlim(0,6)

cdS <- readRDS("data/Cb_only.rds")
FeaturePlot(cdS, "Aldoc",pt.size = 1, reduction = "ref.umap")
