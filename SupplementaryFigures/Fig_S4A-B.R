library(Seurat)

# Load E12.5 scRNA-seq seurat object
cdS <- readRDS("E12/Cb_only.rds")

# Load meta data created with Scanpy
meta <- read.csv("tables/meta.csv", row.names = 1)
rownames(meta) <- paste0(rownames(meta), "-1")

# select cells that were included in scnapy analysi and add pseudotime information
cdS_rna <- subset(cdS, cells = rownames(meta))
meta <- meta[rownames(cdS_rna@meta.data),]
cdS_rna@meta.data <- cbind(cdS_rna@meta.data, meta[,c("velocity_pseudotime","latent_time","dpt_pseudotime")])

# Load ATAC data
library(ArchR);library(Signac)
set.seed(1234)
threadNum = 5
ArchR_proj = "snATACseq/JL20001/ArchR/E12_ATAC"
proj <- loadArchRProject(path = ArchR_proj)

# convert the ArchR object to a seurat object
se <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

rowRanges(se) <-  proj@peakSet
rownames(se) <- GRangesToString(grange = proj@peakSet, sep = c(":", "-"))

cdS_atac <- CreateSeuratObject(counts = assay(se), assay = 'peaks', min.cells=0, min.features=0)
cdS_atac$dataset <- proj$Sample

cdS_atac@meta.data <- as.data.frame(cbind(cdS_atac@meta.data, proj@cellColData))
cdS_atac@meta.data$velocity_pseudotime=meta$velocity_pseudotime[match(cdS_atac@meta.data$predictedCell, rownames(meta))]
cdS_atac@meta.data$latent_time=meta$latent_time[match(cdS_atac@meta.data$predictedCell, rownames(meta))]
cdS_atac@meta.data$dpt_pseudotime=meta$dpt_pseudotime[match(cdS_atac@meta.data$predictedCell, rownames(meta))]

embeddings <- (proj@embeddings$UMAP$df)*-1
colnames(embeddings) <- c("UMAP_1","UMAP_2")
cdS_atac[["umap"]] <- CreateDimReducObject(as.matrix(embeddings), assay = "peaks", key = "umap_")

# remove non-Cb cells
id <- rownames(cdS_atac@meta.data[cdS_atac$cellType %in% c("Endothelium","Pericytes","Microglia","Mes"), ])

cdS_atac1 <- subset(cdS_atac, cells = id, invert = T)  

# plot figure 4SA
p2 <- DimPlot(cdS_atac1, label = T, repel = T, group.by = "cellType")+NoLegend()+NoAxes()+ggtitle("cellType")
p1 <- DimPlot(cdS_rna, reduction = "ref.umap", label = T, group.by = "cellType")+NoLegend()+NoAxes()

pdf(paste0(ArchR_proj,"/Plots/RNA_ATAC.pdf"), w =9.4, h =5)
p1+p2
dev.off()

# plot figure 4SB
p1a=FeaturePlot(cdS_rna, reduction = "ref.umap", label = F, features = "latent_time")+NoLegend()+NoAxes()
p2a=FeaturePlot(cdS_atac1, label = F, features = "latent_time")+NoLegend()+NoAxes()

pdf(paste0(ArchR_proj,"/Plots/RNA_ATAC_pseudotime.pdf"), w =9.6, h =5)
cowplot::plot_grid(p1a,p2a)
dev.off()

