library("Seurat");library(SeuratWrappers);library(tidyverse)
e10 <- readRDS("../E10/Seurat_E10.rds") 
e11 <- readRDS("../E11/Seurat_E11a.rds")     
e12 <- readRDS("../E12/Seurat_E12.rds")
e13 <- readRDS("../E13/Seurat_E13.rds")
e14 <- readRDS("../E14/Seurat_E14.rds")
e15 <- readRDS("../E15/Seurat_E15.rds")
e16 <- readRDS("../E16/Seurat_E16.rds")
e17 <- readRDS("../E17/Seurat_E17.rds")

cdS <- merge(e10, c(e11,e12,e13,e14,e15,e16,e17))
cdS
# 18481 features across 44681 samples within 1 assay 

cdS$dataset <- cdS$orig.ident
cdS$stage <- substr(cdS$dataset, 1,3)
# cdS <- NormalizeData(cdS,verbose = T) %>% FindVariableFeatures(verbose = T) %>% ScaleData(verbose = T) %>% RunPCA(verbose = T)
cdS <- SCTransform(cdS,verbose = T) %>% RunPCA(verbose = T)
cdS.mnn <- RunFastMNN(object.list = SplitObject(cdS, split.by = "dataset"))

# cdS[["mnn"]] <- CreateDimReducObject(Embeddings(cdS.mnn, reduction = "mnn"), assay = "RNA")
cdS[["mnn"]] <- CreateDimReducObject(Embeddings(cdS.mnn, reduction = "mnn"), assay = "SCT")
cdS <- RunUMAP(cdS, reduction = "mnn", dims = 1:30, return.model = T, umap.method = 'uwot')
cdS <- FindNeighbors(cdS, reduction = "mnn", dims = 1:30)
cdS <- FindClusters(cdS, resolution = 1.4)

p1 <- DimPlot(cdS, label = T,group.by = "seurat_clusters")+NoAxes()+NoLegend()
p2 <- DimPlot(cdS, group.by = "stage")+NoAxes()

pdf("figures/UMAP_original.pdf", w = 10, h =5)
p1+p2
dev.off()

dir.create("Seurat4_SCT");
setwd("Seurat4_SCT")
dir.create("figures");dir.create("tables");dir.create("data")

# inspection of clusters
dir.create("figures/clusters")
for (j in levels(Idents(cdS))){
  id <- WhichCells(cdS, idents = j)
  pdf(paste0("figures/clusters/UMAP_Ident", j, ".pdf"))
  print(DimPlot(cdS, cells.highlight = id))
  dev.off()
}

# Remove putative LQC, cluster 25
id <- WhichCells(cdS, idents = 25)

pdf("figures/LQC.pdf")
DimPlot(cdS, cells.highlight = id)+NoLegend()
dev.off()

cdS1 <- subset(cdS, idents = 25, invert = T)

dir.create("figures/UMAPtest")
for (n in seq(from =20, to = 90, by = 5)){
  cdS1 <- RunUMAP(cdS1, reduction = "mnn", dims = 1:40, return.model = T, umap.method = 'uwot', n.neighbors = n, ,min.dist = 0.3)
  pdf(paste0("figures/UMAPtest/UMAP_", n, "L.pdf"))
  print(DimPlot(cdS1, label = TRUE, group.by = "cellType")+NoLegend()+NoAxes())
  dev.off()
}

cdS1 <- RunUMAP(cdS1, reduction = "mnn", dims = 1:30, return.model = T, umap.method = 'uwot', n.neighbors = 30, ,min.dist = 0.3)
cdS1 <- FindNeighbors(cdS1, reduction = "mnn", dims = 1:30)
cdS1 <- FindClusters(cdS1, resolution = 3)
DimPlot(cdS1, label = T)+NoAxes()+NoLegend()

Idents(cdS1) <- cdS1$SCT_snn_res.3
id.isl1 <- WhichCells(cdS1, idents = 43)

plot <- FeaturePlot(cdS1, "Eomes")
ids.ubs <- CellSelector(plot)

Idents(cdS1) <- cdS1$SCT_snn_res.1.4

Idents(cdS1, cells = id.isl1) <- 32
Idents(cdS1, cells = ids.ubs) <- 33

DimPlot(cdS1, label = T)+NoAxes()+NoLegend()

library(dplyr);library(presto)
mGenes <- readRDS(file="~/Desktop/Annotation1.rds")
res <- wilcoxauc(cdS1)
res <- mutate(res,pct.diff=pct_in-pct_out)

res$TF <- "no"
id <- intersect(mGenes$mgi_symbol,res$feature)

x <- res[res$feature %in% id,]
res$TF[res$feature %in% id] <- mGenes$Type[match(x$feature,mGenes$mgi_symbol)]
res$description <- mGenes$description[match(as.character(res$feature),mGenes$mgi_symbol)]
head(res)

res <- res[, c("feature","TF","description","statistic","logFC","auc","pval","padj","pct_in","pct_out","pct.diff","group")]
sig_genes <- res %>% filter(padj < 0.05 & logFC > 0.25 & pct.diff > 25 & auc > 0.7)
table(sig_genes$TF)
openxlsx::write.xlsx(sig_genes,"./tables/sigGenes_all.xlsx")

top <- top_markers(res, n = 10, auc_min = 0.8, pval_max = 1, padj_max = 0.05,
                   pct_in_min = 50, pct_out_max = 100)
openxlsx::write.xlsx(top, "tables/top_markers.xlsx") # annotate in excel

# Manually annotate the cell clusters in Excel and read in the modified table
top1 <- openxlsx::read.xlsx("tables/top_markers.xlsx", sheet = 2)
length(unique(top1$cellType));length(unique(Idents(cdS1)))

cdS1$cellType1 <- plyr::mapvalues(Idents(cdS1),from = top1$rank, to = top1$cellType)
DimPlot(cdS1, label = TRUE, group.by = "cellType1")+NoLegend()+NoAxes()
DimPlot(cdS1, label = TRUE)+NoLegend()+NoAxes()

# Determine the cell cycle phases among the presumptive NPCs
cc.genes <- read.csv("/Volumes/jali/Genome/Cell Cycle Genes.csv", stringsAsFactors = F)
s.genes <- cc.genes$mGene[cc.genes$cellCycle == "G1/S"]
m.genes <- cc.genes$mGene[cc.genes$cellCycle == "G2/M"]
cdS1 <- CellCycleScoring(cdS1, s.features = s.genes, g2m.features = m.genes)

pdf("./figures/UMP_cellCycle.pdf", w =5, h =4)
DimPlot(cdS1, group.by = "Phase")
dev.off()

plot.cc <- DimPlot(cdS1, group.by = "Phase")+NoAxes()

pdf("./figures/UMP_cellCycle1.pdf", w =4.8, h =5)
AugmentPlot(plot.cc,dpi = 300)
dev.off()

p1 <- DimPlot(cdS1, group.by = "cellType1")+NoAxes()+NoLegend()
p2 <- DimPlot(cdS1, group.by = "stage")+NoAxes()
p1a=AugmentPlot(p1,dpi = 300)
p2a=AugmentPlot(p2,dpi = 300)

p3 <- p1a+p2a
p4 <- DimPlot(cdS1, split.by = "stage", ncol = 4, group.by = "cellType1")+NoAxes()+NoLegend()
p4a <- AugmentPlot(p4,dpi = 300)

pdf("figures/UMAP_cellType1.pdf", w =10, h =10)
cowplot::plot_grid(p3,p4a, rel_heights = c(1,1), ncol = 1)
dev.off()

p1 <- DimPlot(cdS1, label = T, group.by = "cellType1")+NoAxes()+NoLegend()
p2 <- DimPlot(cdS1, group.by = "stage")+NoAxes()
p3 <- p1+p2
p4 <- DimPlot(cdS1, split.by = "stage", ncol = 4, group.by = "cellType1")+NoAxes()+NoLegend()

pdf("figures/UMAP_cellType.pdf", w =10, h =10)
cowplot::plot_grid(p3,p4, rel_heights = c(1,1), ncol = 1)
dev.off()

saveRDS(cdS1, "data/E10_17.rds")

cdS
# 36282 features across 44202 samples within 2 assays 
table(cdS$stage)

levels(Idents(cdS))


## Calculate the composition based on stages
comp <- prop.table(table(cdS1$cellType1, cdS1$stage),1)
comp <- data.frame(comp)
colnames(comp) <- c("cellType","Stage","Ratio")

comp1 <- table(cdS1$cellType1, cdS1$stage)
comp1 <- data.frame(comp1)
colnames(comp1) <- c("cellType","Stage","Number")

# Re-order cell types based on the contribution of E12.5 cells
order <- comp %>% 
  spread(Stage, Ratio) %>% 
  mutate(early = (E10+E11)/2, late = (E16+E17)/2, Sort.order = (late - early)) %>% 
  arrange(Sort.order)

comp$cellType <- factor(comp$cellType, levels = as.character(order$cellType))
comp1$cellType <- factor(comp1$cellType, levels = as.character(order$cellType))

library(ggplot2)
pdf("figures/Stage_composition.pdf", h =5, w = 6)
ggplot(comp, aes(x = cellType, y = Ratio, fill = Stage))+
  geom_bar(stat = "identity")+
  cowplot::theme_cowplot(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(comp1, aes(x = cellType, y = Number, fill = Stage))+
  geom_bar(stat = "identity")+
  cowplot::theme_cowplot(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

Idents(cdS1) <- cdS1$SCT_snn_res.2
DimPlot(cdS1, label = T)+NoLegend()

res=FindMarkers(cdS, ident.1 = 18, ident.2 = 3)
tail(res, 20)
