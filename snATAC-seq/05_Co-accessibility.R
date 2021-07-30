library(ArchR)
set.seed(1234)

threadNum = 12
ArchR_proj = "E12_14"
proj <- loadArchRProject(path = ArchR_proj)

# 15.1 Creating Low-Overlapping Aggregates of Cells
# 15.2 Co-accessibility with ArchR
proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = proj,
  corCutOff = 0.5,
  resolution = 100,
  returnLoops = TRUE
)

markerGenes  <- c("Shh","Calb1","Tlx3","Fgf3","Pcdh10","Foxp2")

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(proj)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 10, height = 8)

# 15.3 Peak2GeneLinkage with ArchR
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

saveArchRProject(ArchRProj = proj, outputDirectory = ArchR_proj, load = FALSE)

matrix_p2g <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "cellType", returnMatrices = T)
dim(matrix_p2g$Peak2GeneLinks)
dim(matrix_p2g$RNA$matrix)

# Plot Peaks to genes association heatmap
plot_p2g <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "cellType")
plotPDF(plot_p2g, name = "Peak2GeneHeatmap2", width = 10, height = 8, ArchRProj = proj, addDOC = FALSE)

plot_p2g <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "cellType", k = 30)
plotPDF(plot_p2g, name = "Peak2GeneHeatmap", width = 10, height = 10, ArchRProj = proj, addDOC = FALSE)

# get the entire peak2gene list
p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
head(p2geneDF)

library(dplyr)
# positive correlated peak2gene links
p2gene.sig <- as.data.frame(p2geneDF) %>% filter(FDR < 1e-04 & VarQATAC > 0.25 & VarQRNA > 0.25 & Correlation > 0.45)
head(p2gene.sig)
length(unique(p2gene.sig$peakName))

# negative correlated peak2gene links
p2gene.sig.n <- as.data.frame(p2geneDF) %>% filter(FDR < 1e-04 & VarQATAC > 0.25 & VarQRNA > 0.25 & Correlation < -0.45) %>% arrange(Correlation)
saveRDS(p2gene.sig.n, paste0(ArchR_proj,"/results/peak2gene_negative.rds"))

p2gene.sig.n[601:800,]
length(unique(p2gene.sig.n$peakName))
length(unique(p2gene.sig.n$geneName))
p2gene.sig.n[p2gene.sig.n$geneName == "Chd3",]
intersect("Gdf10", unique(p2gene.sig.n$geneName))



library(TxDb.Mmusculus.UCSC.mm10.knownGene)
mygenes = subset(transcripts(TxDb.Mmusculus.UCSC.mm10.ensGene, columns=c("tx_id", "tx_name","gene_id")), gene_id %in% p2g.peaks$linkedGene)
mygenes = subset(transcripts(TxDb.Mmusculus.UCSC.mm10.ensGene, gene_id %in% p2g.peaks$p2g.peaks$linkedGene))
mygenes.transcripts = subset(GenomicFeatures::transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, columns=c("tx_id", "tx_name","gene_id")), gene_id %in% p2g.peaks$linkedGene)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g.peaks <- proj@peakSet[p2g$idxATAC,]
p2g$peaks <- Signac::GRangesToString(grange = p2g.peaks, sep = c(":", "-"))
p2g$gene <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2g.db$gene <- getFeatures(proj,useMatrix = "GeneIntegrationMatrix")[p2g.db$idxRNA]
head(p2g.db)




p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 10000,
  returnLoops = TRUE
)
p2g[[1]]

markerGenes  <- c("Shh","Calb1","Bcl11b","Fgf3","Pcdh10","Foxp2")
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = markerGenes, 
  upstream = 100000,
  downstream = 100000,
  loops = getPeak2GeneLinks(proj)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 12, height = 8)

markerGenes  <- c("Nrk","Sphk1","Htr2a","Gabra1","Gabra2","Foxp1","Cck")

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = markerGenes, 
  upstream = 100000,
  downstream = 100000,
  loops = getPeak2GeneLinks(proj)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks1.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 10, height = 5)

markerGenes  <- c("Plcb4","Bcl11b","Ndnf","Rorb","Pcp2","Nrgn","Eya1","Eya2","Etv4")

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = markerGenes, 
  upstream = 100000,
  downstream = 100000,
  loops = getPeak2GeneLinks(proj)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks2.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 10, height = 5)

p1 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = "Bcl11b", 
  upstream = 105000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj)
)

p2 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = "Bcl11b", 
  upstream = 105000,
  downstream = 50000,
  loops = getCoAccessibility(proj)
)

plotPDF(plotList = c(p1,p2), 
        name = "Plot-Tracks-Peak2_Bcl11b.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 12, height = 9)

p1 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = "Dab1", 
  upstream = 55000,
  downstream = 250000,
  loops = getPeak2GeneLinks(proj)
)

p2 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = "Dab1", 
  upstream = 55000,
  downstream = 550000,
  loops = getCoAccessibility(proj)
)

plotPDF(plotList = c(p1,p2), 
        name = "Plot-Tracks-Peak2_Dab1.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 12, height = 9)


p1 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = "Calb1", 
  upstream = 250000,
  downstream = 30000,
  loops = getPeak2GeneLinks(proj)
)

p2 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cellType", 
  geneSymbol = "Calb1", 
  upstream = 250000,
  downstream = 30000,
  loops = getCoAccessibility(proj)
)

plotPDF(plotList = c(p1,p2), 
        name = "Plot-Tracks-Peak2_Calb1.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 12, height = 9)


grid::grid.newpage()
grid::grid.draw(p1$Calb1)

ArchRBrowser(proj)


