## ==== annotate cluster-unique peaks ======
library(ArchR);library(Signac)
library(ChIPseeker);library(clusterProfiler)
library(tidyverse);library(genomation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

threadNum = 5
ArchR_proj = "E12_14"
proj <- loadArchRProject(path = ArchR_proj)

# # Obtain signatures (if it has not been run before)
# path_to_signatures <- '/Volumes/Backup_JamesLi/GSE60731_RAW/Consensus_mm10/'
# fileNames <- paste(path_to_signatures, list.files(path_to_signatures), sep='')
# labels  <- c("ASCL1","ATOH1","Cb_8W_DHS","Cb_DHS","CTCF_E14Cb","Foxp2",'MEF2A','MEF2D',"H3K27me3_Ezh2KO","H3K27me3","Hb_E14_DHS",'P0_DNaseHS', 'P14_DAaseHS', 'P3_DNaseHS','P60_DNaseHS','P60_Zic_ChIP','P7_DNaseHS','P7_H3K27ac','P7_Zic_ChIP',"PTF1A","VISTA_Cb")
# 
# signatures <- readBed(fileNames[1]) 
# signatures$name <- labels[1]
# 
# for (i in 2:length(fileNames)) {
#   file <- fileNames[i]
#   temp <- readBed(file) 
#   temp$name <- labels[i]
#   signatures <- append(signatures, temp)
# }
# 
# load('/Volumes/Backup_JamesLi/Experiments1/scRNAseq/CerebellarAnlage/ATACseq/MouseAtlas/Sig_peaks.Rdata')
# PC.mm10$name <- "PC_mouseAtlas"
# GC.mm10$name <- "GC_mouseAtlas"
# signatures <- append(signatures, PC.mm10)
# signatures <- append(signatures, GC.mm10)
# 
# signatures.list=GenomicRanges::split(signatures, signatures$name,drop=TRUE)
# 
# path_to_ChromHMM <- '/Volumes/Backup_JamesLi/GSE60731_RAW/ChromHMM/'
# setwd(path_to_ChromHMM)
# fileNames <- list.files(path = path_to_ChromHMM, pattern = "rep1")
# chrHMM <- readBed(fileNames[1]) 
# names_short <- gsub("_day_hindbrain_rep1_18state_dense.sorted.bed","", fileNames[1])
# names_short <- gsub("embryonic_","E", names_short)
# names_short <- gsub("postnatal_","P", names_short)
# chrHMM$stage <- names_short
# 
# for (i in 2:length(fileNames)) {
#   file <- fileNames[i]
#   names_short <- gsub("_day_hindbrain_rep1_18state_dense.sorted.bed","", file)
#   names_short <- gsub("_day_hindbrain_rep2_18state_dense.sorted.bed","_rep2", names_short)
#   names_short <- gsub("embryonic_","E", names_short)
#   names_short <- gsub("postnatal_","P", names_short)
#   temp <- readBed(file) 
#   temp$stage <- names_short
#   chrHMM <- append(chrHMM, temp)
# }
# chrHMM$group <- paste0(chrHMM$stage, "_",chrHMM$name)
# chrHMM.list=GenomicRanges::split(chrHMM, chrHMM$group,drop=TRUE)
# save(signatures.list,chrHMM.list, file = "/Volumes/Backup_JamesLi/GSE60731_RAW/Consensus_mm10/forEnrichment.Rdata")

# load previously save signature files
load("/Volumes/Backup_JamesLi/GSE60731_RAW/Consensus_mm10/forEnrichment.Rdata")

proj <- addPeakAnnotations(ArchRProj = proj, regions = signatures.list, name = "Signature", force = TRUE)
proj <- addPeakAnnotations(ArchRProj = proj, regions = chrHMM.list, name = "ChromHMM", force = TRUE)
saveArchRProject(ArchRProj = proj, outputDirectory = ArchR_proj, load = FALSE)

se_da <- readRDS("E12_14/results/UniquePeak_se.rds")

## test enrichment for signature.list
enrichRegions <- peakAnnoEnrichment(
  seMarker = se_da,
  ArchRProj = proj,
  peakAnnotation = "Signature",
  cutOff = "binary > 0"
)

heatmapRegions <- plotEnrichHeatmap(enrichRegions, cutOff = 0,transpose = TRUE)

pdf(paste0(ArchR_proj,"/Plots/ChIP_Signature.pdf"), width=5, height=5)
ComplexHeatmap::draw(heatmapRegions, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

enrichRegions.chromHMM <- peakAnnoEnrichment(
  seMarker = se_da,
  ArchRProj = proj,
  peakAnnotation = "ChromHMM",
  cutOff = "binary > 0"
)
heatmapRegions <- plotEnrichHeatmap(enrichRegions.chromHMM, cutOff = 1,transpose = F)

pdf(paste0(ArchR_proj,"/Plots/ChromHMM_Enrichment.pdf"), width=6, height=9)
ComplexHeatmap::draw(heatmapRegions, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

library(regioneR)
vista <- readBed("/Volumes/Backup_JamesLi/GSE60731_RAW/Consensus_mm10/VISTA_Cb.tsv")
cCRE <- StringToGRanges(unique(cCREs), sep = c(":", "-"))
numOverlaps(vista, cCRE, count.once=TRUE)
# 72
numOverlaps(randomizeRegions(vista, genome = "mm10"), cCRE, count.once=TRUE)

pt <- overlapPermTest(A=vista, B=cCRE, ntimes=500, genome='mm10', non.overlapping=FALSE, verbose=TRUE)
pt
# $numOverlaps
# P-value: 0.00199600798403194
# Z-score: 33.0859
# Number of iterations: 500
# Alternative: greater
# Evaluation of the original region set: 126
# Evaluation function: numOverlaps
# Randomization function: randomizeRegions

pt <- overlapPermTest(A=vista, B=cCRE, ntimes=500, genome='mm10', non.overlapping=FALSE, verbose=TRUE)
pt
plot(pt)

dhs <- readBed("/Volumes/Backup_JamesLi/GSE60731_RAW/Consensus_mm10/Hb_E14_DHS.bed")
dhs <- dhs[1:3450,]
peaks=GRangesToString(grange = proj@peakSet, sep = c(":", "-"))
peaks.all <- StringToGRanges(peaks, sep = c(":", "-"))
numOverlaps(dhs, peaks.all, count.once=TRUE)
# [1] 2585
numOverlaps(randomizeRegions(dhs, genome = "mm10"), peaks.all, count.once=TRUE)
# [1] 523

pt <- overlapPermTest(A=dhs, B=peaks.all, ntimes=500, genome='mm10', non.overlapping=FALSE, verbose=TRUE)
pt
# $numOverlaps
# P-value: 0.00199600798403194
# Z-score: 123.5035
# Number of iterations: 500
# Alternative: greater
# Evaluation of the original region set: 3765
# Evaluation function: numOverlaps
# Randomization function: randomizeRegions
plot(pt)

plot(pt)

sum <- data.frame(pval = pt$numOverlaps[[1]],Zscore=pt$numOverlaps[[6]],Overlaps=pt$numOverlaps[[4]],N = pt$numOverlaps[[2]],bed = file)

sum <- rbind(sum,temp)
data.frame(sum)

peaks=GRangesToString(grange = proj@peakSet, sep = c(":", "-"))
p2g.peaks <- readRDS(paste0(ArchR_proj,"/results/peak2gene.rds"))
cCREs <- p2g.peaks$peaks




# reload cluster-unique peaks results
uf <- readRDS(paste0(ArchR_proj,"/results/Unique_Peaks.rds"))
dim(uf$binaryMat)
# [1] 46772    30

da <- as.data.frame(uf$binaryMat)
n=ncol(da)
da$gene <- rownames(da)
da$padj <- uf$padj
da <- gather(da, cluster, binary, 1:n)
da <- da[da$binary == 1,]

# Generate a list of unique peaks of each cluster
peakList=list()
for(i in 1:length(unique(da$cluster))){
  cluster_i <- unique(da$cluster)[i]
  peakList[i] <- StringToGRanges(da$gene[da$cluster == cluster_i], sep = c(":", "-"))
}
names(peakList) <- unique(da$cluster)

peakList[["cCRE"]] <- StringToGRanges(unique(cCREs), sep = c(":", "-"))
peakList[["all"]] <- StringToGRanges(peaks, sep = c(":", "-"))

# peakList <- peakList[paste0("Cluster",0:23)]
# peakList.unique <- peakList
# names(peakList.unique) <- 0:23

# annotate peaks
peakAnnoList <- lapply(peakList, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

# output plots
pdf(paste0(ArchR_proj,"/Plots/AnnoBar_unique.pdf"), h =6, w = 6)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()



upsetplot(peakAnnoList[[1]])
upsetplot(peakAnnoList[[1]], vennpie = T)

# annoate peaks identified by Signac
library(dplyr)
da_peaks_all <- readRDS("data/da_peaks_all.rds")

# Generate a list of unique peaks of each cluster
peakList=list()
for(i in 1:length(unique(da_peaks_all$cluster))){
  cluster_i <- unique(da_peaks_all$cluster)[i]
  peakList[i] <- StringToGRanges(da_peaks_all$gene[da_peaks_all$cluster == cluster_i], sep = c(":", "-"))
}
names(peakList) <- unique(da_peaks_all$cluster)
peakList.da <- peakList

# annotate peaks
peakAnnoList <- lapply(peakList, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

# output plots
pdf('figures/AnnoBar_da_all.pdf', h =6, w = 6)
plotAnnoBar(peakAnnoList)
dev.off()

pdf('figures/Dist2TSS_da_all.pdf', h =6, w = 6)
plotDistToTSS(peakAnnoList,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()

# Find the intersected peaks
peakList=list()
for(i in 0:23){
  peaks.unique <- da$gene[da$cluster == paste0("Cluster",i)]
  peaks.da <- da_peaks_all$gene[da_peaks_all$cluster == i]
  print(paste0("number of peaks of cluster", i, ": ", length(intersect(peaks.unique, peaks.da))))
  peakList[i] <- StringToGRanges(intersect(peaks.unique, peaks.da), sep = c(":", "-"))
}

# annotate peaks
peakAnnoList <- lapply(peakList, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

# output plots
pdf('figures/AnnoBar_intersect.pdf', h =6, w = 6)
plotAnnoBar(peakAnnoList)
dev.off()

pdf('figures/Dist2TSS_intersect.pdf', h =6, w = 6)
plotDistToTSS(peakAnnoList,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()


# compare functional profiles of each cluster
names(peakAnnoList) <- 0:23
genes = lapply(peakAnnoList, function(x) as.data.frame(x)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           organism="mouse")
write.csv(compKEGG, 'tables/comKEGG.csv')

pdf('figures/compKEGG.pdf', h = 12, w =15)
dotplot(compKEGG, showCategory = 5, title = "KEGG Pathway Enrichment Analysis")
dev.off()

