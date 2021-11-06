library(Seurat);library(dplyr);library(openxlsx);library(ggplot2);
library(pheatmap);library(RColorBrewer)

# load disease candidate gene lists that download from DisGeNet
tab=openxlsx::read.xlsx("tables/DiseaseGenes.xlsx")

# Reload seurat object and select Cb-only cells
cdS_carter <- readRDS("E10-E13/data/E10_13_final.rds")
DimPlot(cdS_carter, label = T)+NoLegend()

# Keep Cb-only brain cells
cdS_carter_sub <- subset(cdS_carter, idents = c(17,18,20:22), invert=T)

Idents(cdS_carter_sub)=cdS_carter_sub$cellType
DimPlot(cdS_carter_sub, label = T)+NoLegend()

DWM=tab$symbol[tab$Disease1=="DWM"]
ids <- intersect(DWM, rownames(cdS_carter))

# Calculate cell-type average expression matrix
cluster.averages <- AverageExpression(object = cdS_carter_sub, assays = "SCT")
mat=cluster.averages$SCT
mat=data.frame(mat)
mat <- as.data.frame(mat)[ids,]

z.mat <- t(scale(t(mat), center=TRUE, scale=TRUE)) 
z.mat <- z.mat[!is.na(rowSums(z.mat)),]

pdf("figures/DWM_genes_carter.pdf", h = 6, w = 5)
pheatmap(z.mat, labels_row = NULL,border_color =NA,show_rownames = T,
         col=colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), 
         cutree_rows = 4,cutree_cols = 1,
         treeheight_row = 0,
         main = "Cell-type average expression of DWM genes \n E10-E13 (Carter)")
dev.off()

# Repeat the analysis with in-house E12.5 data
cdS=readRDS("E12.5/data/Cb_final.rds")
DimPlot(cdS, label = T, group.by = "cellType", reduction = "ref.umap", repel = T)+NoLegend()+NoAxes()

# subset cb-only cells
cellTypeCB <- c("NPC.cyc","PC","NPC","CN.late","GABA.Pre","GABA.Pro","Glu.Pro","C4.Otp","CN.early2","CN.early1","CN.gaba","C4.Otp.Pre","Isth.N","NPCa","PTZ","IN")
Idents(cdS)=cdS$cellType
cdSsub = subset(cdS, idents= cellTypeCB)
DimPlot(cdSsub, label = T, reduction = "ref.umap", repel = T)+NoLegend()+NoAxes()

ids <- intersect(DWM, rownames(cdSsub))

# Calculate cell-type average expression matrix
cluster.averages <- AverageExpression(object = cdSsub, assays = "SCT")
mat=cluster.averages$SCT
mat=data.frame(mat)
mat <- as.data.frame(mat)[ids,]

z.mat <- t(scale(t(mat), center=TRUE, scale=TRUE)) 
z.mat <- z.mat[!is.na(rowSums(z.mat)),]

pdf("figures/DWM_genes_E12.pdf", h = 6, w = 5)
pheatmap(z.mat, labels_row = NULL,border_color =NA,show_rownames = T,
         col=colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), 
         cutree_rows = 4,cutree_cols = 1,
         treeheight_row = 0,
         main = "Cell-type average expression of DWM genes \n E2.5 (in-house)")
dev.off()

# Perform enrichment analysis
## Find and filter differential marker genes 
res <- FindAllMarkers(object = cdSsub, logfc.threshold = 0.3,
                            test.use = "wilcox", min.pct = 0.25, min.diff.pct = 0.1,
                            print.bar = TRUE, only.pos = TRUE, max.cells.per.ident = Inf,
                            return.thresh = 0.001, do.print = FALSE, random.seed = 1,
                            min.cells.gene = 3, min.cells.group = 3, latent.vars = NULL,
                            assay.type = "SCT")
res <- res %>% mutate(pct.diff=pct.1-pct.2)
sig_genes <-  res %>% filter(pct.2<0.3, p_val_adj<0.0001, pct.diff>0.2) 
table(sig_genes$cluster)

all=rownames(cdSsub)
disease_genes_filtered=intersect(all, DWM)
cellType=as.character(unique(res$cluster))

# perform Fisher's exact test for markers genes of all cb cell types
p.value=list();odds_ratio=list();intersect=list();notInList=list()

for (j in 1:length(cellType)){
  name=cellType[j]
  markers=sig_genes$gene[sig_genes$cluster ==name]
  m = length(intersect(markers, disease_genes_filtered));m
  n=length(setdiff(disease_genes_filtered,markers));n
  k=length(setdiff(markers,disease_genes_filtered));k
  i=length(all)-m-n-k;i
  contingency.table <- matrix(c(m,k,n,i), nrow = 2)
  ft=fisher.test(contingency.table, alternative = 'greater')
  p.value[[name]]=ft$p.value
  odds_ratio[[name]]=ft$estimate
  intersect[[name]]=m
  notInList[[name]]=n
}

# output histogram of enrichment results
df=data.frame(p.value=unlist(p.value), enrichment=unlist(odds_ratio), intersect=unlist(intersect), notInList = unlist(notInList))
df=df %>% arrange(p.value)
df$cellType <- rownames(df)
df$p_val <- -log10(df$p.value)

txt=paste0(df$intersect, "/", df$notInList)

library(cowplot)
p <- ggplot(df)  + 
  geom_bar(aes(x=reorder(cellType, p_val), y=enrichment),stat="identity", fill="cyan")+ 
  xlab("cell type")+
  geom_line(aes(x=reorder(cellType, p_val), y=p_val*4, group = 1),stat="identity",colour="red")+
  geom_point(aes(x=reorder(cellType, p_val), y=p_val*4, group = 1),stat="identity",colour="red")+
  geom_text(aes(label=txt, x=cellType, y=0.85*enrichment), colour="black")+
  scale_y_continuous(sec.axis = sec_axis(~./4,name="-log10(p.value)"))+coord_flip()+
  geom_hline(yintercept=2, linetype="dashed", 
               color = "black", size=0.5)+
  theme_cowplot(12)

pdf("figures/enrichement.pdf", h =6, w=4)
p
dev.off()

# extract PTZ-specific DWM genes
markers=sig_genes$gene[sig_genes$cluster =="PTZ"]
genes <- intersect(markers, disease_genes_filtered)
#"Zic4"  "Lamc1" "Nphp1" "Fgf17" "Flna"  "Zic5" 

# Plot DWM Module on UMAP
cdS <- AddModuleScore(cdS, genes, name="DWM")
FeaturePlot(cdS, "DWM1", min.cutoff = 1, reduction = "ref.umap",raster = T, pt.size =2, cols = c("grey92","blue")) + NoAxes()+ggtitle("DW Module")+guides(colour = guide_legend(title.hjust =0.5))




