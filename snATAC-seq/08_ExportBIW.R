library(ArchR);
set.seed(1)
threadNum = 5
ArchR_proj = "E12_14"

addArchRThreads(threads = threadNum) 
proj <- loadArchRProject(path = ArchR_proj)

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cellType", embedding = "UMAP")

getGroupBW(
  ArchRProj = proj,
  groupBy = "cellType",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

peaks <- data.frame(proj@peakSet);dim(peaks)
saveRDS(peaks, "E12_14/GroupBigWigs/allPeaks.rds")
write.table(peaks[,1:3],file="tables/ArchR_peaks.txt",quote=FALSE,
            row.names=FALSE,col.names=FALSE,sep="\t")

