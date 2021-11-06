library(Seurat);library(dplyr)
count <- Read10X_h5("../filtered_gene_bc_matrices_h5.h5") 

cdS <- CreateSeuratObject(count, project = "E12") 
cdS[["percent.mt"]] <- PercentageFeatureSet(cdS, pattern = "^mt-")
cdS
# 28692 features across 6973 samples within 1 assay 

dir.create("figures");dir.create("tables");dir.create("data")
pdf("./figures/VlnPlot_QC.pdf",h=8,w=12)
VlnPlot(object = cdS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1 <- FeatureScatter(object = cdS, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(object = cdS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cowplot::plot_grid(p1,p2)
dev.off()

library(ggplot2)
VlnPlot(object = cdS, features = "nFeature_RNA")+
  geom_hline(yintercept = c(1000, 6000), linetype="dotted", 
             color = c("blue","red"), size=1)
VlnPlot(object = cdS, features = "percent.mt")+
  geom_hline(yintercept = 8, linetype="dotted", 
             color = "red", size=1)

# filter cells
cdS1 <- subset(cdS, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 8 & nCount_RNA < 30000)

pdf("./figures/VlnPlot_QC1.pdf",h=8,w=12)
VlnPlot(object = cdS1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1 <- FeatureScatter(object = cdS1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(object = cdS1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cowplot::plot_grid(p1,p2)
dev.off()

# number of cells removed
ncol(cdS) - ncol(cdS1)
# [1] 158

cdS=cdS1
rm(cdS1)

# ==== normalization of data ========
cdS <- SCTransform(cdS, verbose = F) %>% RunPCA() %>% RunUMAP(dims = 1:40,return.model = T,umap.method = "uwot",n.neighbors = 30L)
cdS <- FindNeighbors(object = cdS, dims = 1:40, verbose = FALSE)
cdS <- FindClusters(object = cdS, verbose = FALSE)
cdS <- RunUMAP(cdS, dims = 1:30, return.model = T, umap.method = 'uwot', n.neighbors = 30,min.dist = 0.3)

pdf("./figures/UMAP_01.pdf",h=5,w=5)
DimPlot(cdS, label = TRUE)+NoLegend()+NoAxes()
dev.off()

# Inspect the distributio reads in each cluster. No obvious outliers
pdf("./figures/VlnPlot_clusters_01.pdf",h=8,w=12)
VlnPlot(object = cdS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 1)
FeaturePlot(cdS, c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
dev.off()

# Mannually select cell clusters
p1 <- DimPlot(cdS, label = TRUE)+NoLegend()+NoAxes()
id1 <- CellSelector(p1)
id2 <- CellSelector(p1)
Idents(cdS, cells = id1) <- 26
Idents(cdS, cells = id2) <- 27

DimPlot(cdS, label = TRUE)+NoLegend()+NoAxes()
cdS$seurat_clusters1 <- Idents(cdS)

# find cluster-specific markers for cell-type annotation
library(presto)
res <- wilcoxauc(cdS, "seurat_clusters1")
top <- top_markers(res, n = 10)
openxlsx::write.xlsx(top, "tables/top_markers.xlsx") # annotate in excel

# Annotate cell types in excel and add results to the Seurat object
top1 <- openxlsx::read.xlsx("tables/top_markers.xlsx", sheet = 2)

cdS$cellType <- plyr::mapvalues(Idents(cdS),from = top1$rank, to = top1$cellType)
DimPlot(cdS, label = TRUE, group.by = "cellType")+NoLegend()+NoAxes()

p2 <- DimPlot(cdS, label = TRUE, group.by = "cellType", repel = T)+NoLegend()+NoAxes()
p1 <- DimPlot(cdS, label = TRUE, group.by = "seurat_clusters1", repel = T)+NoLegend()+NoAxes()
pdf("./figures/UMAP_03_cellType.pdf",h=5,w=10)
p1+p2
dev.off()

# save results
saveRDS(cdS, "data/E12.rds")

# find cell-type specific markers
res <- wilcoxauc(cdS, "cellType")
mGenes <- readRDS(file="~/Desktop/Annotation1.rds")
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


# Project cells according to the Carter's E10.5-13.5 data
## select cerebellum-cells
cdS.cb <- subset(cdS, idents = c(8,22,24,13,17), invert =T)

ref <- readRDS("data/E10_13.rds")
DimPlot(ref, label = T, group.by = "cellType1")+NoLegend()

anchors <- FindTransferAnchors(
  reference = ref,
  query = cdS.cb,
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

test <- MapQuery(
  anchorset = anchors,
  query = cdS.cb,
  reference = ref,
  refdata = "cellType1",
  reference.reduction = "pca", 
  reduction.model = "umap"
)

p1 <- DimPlot(test, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE)+NoAxes()+NoLegend()
p2 <- DimPlot(test, reduction = "ref.umap", group.by = "cellType", label = TRUE, repel = TRUE)+NoAxes()+NoLegend()
p3 <- DimPlot(test, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE)+NoAxes()+NoLegend()
p4 <- DimPlot(test, reduction = "umap", group.by = "cellType", label = TRUE, repel = TRUE)+NoAxes()+NoLegend()

pdf("figures/UMAP_Carter_projection1.pdf", w =12, h=12)
(p1+p2)/(p3+p4)
dev.off()

saveRDS(test, "data/Cb_only.rds")

