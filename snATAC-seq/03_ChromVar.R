library(ArchR);library(Signac);library(cowplot)   

# Setting up
set.seed(123)
addArchRThreads(threads = 5) 
proj <- loadArchRProject(path = "E12_14")

# 13.1 Motif Deviations
if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

motifs <- c("Lhx1", "En1", "Foxp1", "Etv1", "Tfap2b", "Atoh1")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

proj <- addImputeWeights(proj)
p <- plotGroups(ArchRProj = proj, 
                groupBy = "cellType", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})

pdf("E12_14/Plots/Motiff_ChromVar.pdf")
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
dev.off()

# Save object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "E12_14",
  overwrite = TRUE,
  load = FALSE
)

system("wget https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra-Human-Motifs.rds")
motifPWMs <- readRDS("Vierstra-Human-Motifs.rds")
proj <- addMotifAnnotations(proj, motifPWMs = motifPWMs, name = "Vierstra")
# proj <- addDeviationsMatrix(proj, peakAnnotation = "Vierstra", force = TRUE)

names(motifPWMs)

