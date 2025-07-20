# FastCAR
library(Matrix)
# seurat v4
library(Seurat)
library(qlcMatrix)
library(FastCAR)
library(ggplot2)
library(gridExtra)

#Specify the locations of the expression matrices
message("load data...\n")
cellExpressionFolder  = c("F18W/outs/filtered_feature_bc_matrix/")
fullMatrixFolder      = c("F18W/outs/raw_feature_bc_matrix/")

#Load both the cell matrix and the full matrix
# generate a list
cellMatrix     = read.cell.matrix(cellExpressionFolder)
fullMatrix     = read.full.matrix(fullMatrixFolder)

# all dgCMatrix
#32285 genes x 20000 cells, filtered gene expression matrix
filtered_GEX <- cellMatrix$'Gene Expression' 
#107476 peaks x 20000 cells, filtered atac matrix
filtered_atac <- cellMatrix$Peaks

#32285 genes x 648028 cells, raw gene expression matrix
raw_GEX <- fullMatrix$'Gene Expression' 
#107476 x 648028 cells, raw atac matrix
raw_atac <- fullMatrix$Peaks
message("done.\n")

message("run FastCAR...\n")
#======x require a dgCMatrix====
ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = raw_GEX, 
                                           start = 10, 
                                           stop = 500, 
                                           by = 10, 
                                           contaminationChanceCutoff = 0.05)

pdf("ambient_profile.pdf", width = 6, height = 5)
plot.ambient.profile(ambProfile)
dev.off()

emptyDropletCutoff = recommend.empty.cutoff(ambProfile) # recommand 490
#emptyDropletCutoff        = 100 
contaminationChanceCutoff = 0.05

ambientProfile = FastCAR::determine.background.to.remove(raw_GEX, emptyDropletCutoff, contaminationChanceCutoff)
filtered_GEX     = remove.background(filtered_GEX, ambientProfile)
#32285genes x 20000cells
message("filtered gene expression dim: ", nrow(filtered_GEX),"genes x ",ncol(filtered_GEX),"cells \n")

FastCAR_seuratObject = CreateSeuratObject(filtered_GEX) 
message("save FastCAR seurat object...\n ")
saveRDS(FastCAR_seuratObject, file = "seuratobj_rds/FastCAR_seuratObject.rds")

# 原始计数矩阵（稀疏矩阵）
FastCAR_seuratObject[["RNA"]]@counts[1:5, 1:5]
# 归一化数据矩阵（默认是 log-normalized）
FastCAR_seuratObject[["RNA"]]@data[1:5, 1:5]
# meta information
head(FastCAR_seuratObject@meta.data)
message("done.\n")

#===========================================================
message("run EmptyNN...\n")
#EmptyNN
library(EmptyNN)
# Python tensorflow 2.15.0 // Python Keras 2.15.0 // R keras 2.15.0 // cudnn 8.9.4 // cuda 12
# dependencies('caret', python 'tensorflow')

# Transpose the count matrix, so rows are cells and columns are genes
EmptyNN_counts <- t(raw_GEX)

#按需分配显存
library(reticulate)
tf <- import("tensorflow")
gpus <- tf$config$experimental$list_physical_devices('GPU')
if (length(gpus) > 0) {
  for (gpu in gpus) {
    tf$config$experimental$set_memory_growth(gpu, TRUE)
  }
}
# Run emptynn()
nn.res <- emptynn(EmptyNN_counts, threshold = 100, k = 10, iteration = 10, verbose = TRUE)
#32285genes x 21206cells 
message("retained counts dims: ", nrow(EmptyNN_counts[nn.res$nn.keep,]),"genes x ",ncol(EmptyNN_counts[nn.res$nn.keep,]), "cells \n")

# Downstream analysis
# input Seurat should: genes x cells
EmptyNN_seuratObject <- runSeurat(counts = t(EmptyNN_counts[nn.res$nn.keep,]), resolution = 0.2)
message("save EmptyNN seurat object...\n ")
saveRDS(EmptyNN_seuratObject, file = "seuratobj_rds/EmptyNN_seuratObject.rds")

pdf("EmptyNN_tsne.pdf", width = 10, height = 8)
DimPlot(EmptyNN_seuratObject, reduction = "tsne") + ggtitle("EmptyNN") + NoLegend()
dev.off()
message("done.\n")

#===========================================================
message("run Emptydrops...\n")
#Emptydrops
library(DropletUtils)
library(Matrix)
library(UpSetR)
library(SingleCellExperiment)

saveRaster <- function(...) {
  tmp <- tempfile(fileext=".png")
  png(tmp, width=7, height=7, units="in", res=300, pointsize=12)
  par(mar=c(0,0,0,0))
  options(bitmapType="cairo")
  plot(..., xlab="", ylab="", axes=FALSE)
  dev.off()
  tmp
}

loadRaster <- function(fpath, log=FALSE) {
  limits <- par("usr")
  img <- png::readPNG(fpath)
  if (log) { limits <- 10^limits }
  rasterImage(img, limits[1], limits[3], limits[2], limits[4])
  invisible(limits)
}

plotHistogramOutline <- function(breaks, heights, ...) {
  x <- rep(breaks, each=2)
  y <- c(0, rep(heights, each=2), 0)
  lines(x, y, ...)
}

Emptydrops_counts <- raw_GEX
stats <- barcodeRanks(Emptydrops_counts)
expected_cells_num <- 20000

# Testing emptyDrops.
e.out <- emptyDrops(Emptydrops_counts)
e.keep <- e.out$FDR <= 0.001
e.keep[is.na(e.keep)] <- FALSE

#32285genes x 22258cells
EmptyDrops_retained <- Emptydrops_counts[,e.keep]
message("filtered gene expression dim: ", nrow(EmptyDrops_retained),"genes x ",ncol(EmptyDrops_retained),"cells \n")
EmptyDrops_seuratObject = CreateSeuratObject(EmptyDrops_retained) 
message("save EmptyDrops seurat object...\n ")
saveRDS(EmptyDrops_seuratObject, file = "seuratobj_rds/EmptyDrops_seuratObject.rds")
message("done.\n")


# 提取拐点信息（metadata）
inflection <- metadata(stats)$inflection
knee <- metadata(stats)$knee
# Keeping everything above the knee point.
k.keep <- knee <= stats$total

# Using the CellRanger approach.
c.keep <- defaultDrops(Emptydrops_counts, expected=expected_cells_num)

############################
# Quantifying the intersections.
at.least.one <- k.keep | e.keep | c.keep
collected <- data.frame(EmptyDrops=as.integer(e.keep), 
                        `Knee point`=as.integer(k.keep), 
                        CellRanger=as.integer(c.keep),
                        check.names=FALSE)[at.least.one,]

pdf("EmptyDrops_intersect.pdf", onefile=FALSE)
upset(collected, sets.bar.color = "#56B4E9", point.size=5, order.by = "freq",
      text.scale=c(2, 1.5, 1, 1, 1.5, 1.5))
dev.off()

############################
# Having a look at the distribution of retained cells in more detail.
limits <- log10(range(stats$total[at.least.one]))
breaks <- seq(limits[1], limits[2], length.out=21)
modes <- list("EmptyDrops"=e.keep, "Knee point"=k.keep, "CellRanger"=c.keep)

collected.x <- collected.y <- list()
for (mode in names(modes)) {
  #EmptyDrops:22258cells, Knee point:14952 cells, CellRanger:11774 cells
  keep <- modes[[mode]]
  d <- hist(log10(stats$total[keep]), breaks=breaks, plot=FALSE)
  collected.x[[mode]] <-d$mid
  collected.y[[mode]] <-d$count
}
xrange <- range(sapply(collected.x, range))
yrange <- range(sapply(collected.y, range))*1.2

pdf("EmptyDrops_kept.pdf")
plot(0,0,type="n", xlim=xrange, ylim=yrange, xlab=expression(Log[10]~"Total count"), 
     ylab="Number of cells", cex.axis=1.2, cex.lab=1.4, main="EmptyDrops", cex.main=1.4)
shift <- 0
my_colors <- c("EmptyDrops" = "blue", "Knee point" = "green", "CellRanger" = "red")
for (mode in names(modes)) { 
  plotHistogramOutline(breaks+shift, collected.y[[mode]], col=my_colors[mode], lwd=2)
  shift <- shift + 0.01
}
legend("topleft", col=my_colors[names(modes)], legend=names(modes), lwd=2, cex=1.2)
dev.off()

############################
# Examining the distribution of deviances.
X <- e.out$Total/1000
Y <- -e.out$LogProb/1000
xlim <- range(X[!is.na(e.out$LogProb)]) 
tmp <- saveRaster(X, Y, pch=16, cex=0.5, col=ifelse(e.keep, "red", "grey80"), log="xy", xlim=xlim)

pdf("EmptyDrops_deviances.pdf")
par(mar=c(5.1, 5.1, 4.1, 1.1))

plot(X, Y, log="xy", 
     xlab=expression("Total count (x "*10^3*")"), 
     ylab=expression("-Log likelihood (x "*10^3*")"), 
     xlim=xlim, type="n", cex.axis=1.2, cex.lab=1.4,
     main="EmptyDrops", cex.main=1.4)
loadRaster(tmp, log=TRUE)
box()

legend("bottomright", col=c("red", "grey80"), pch=16, cex=1.2,
       legend=c("Putative cell", "Empty droplet"))
dev.off()

############################
# Creating a MA plot.
ambient.prof <- rowSums(Emptydrops_counts[,colSums(Emptydrops_counts) <= 100])
cell.prof <- rowSums(Emptydrops_counts[,e.keep])
summed <- cbind(ambient.prof, cell.prof)
vals <- edgeR::cpm(summed, log=TRUE, prior.count=5)
M <- vals[,1] - vals[,2]
A <- (vals[,1] + vals[,2])/2

tmp <- saveRaster(A, M, pch=16, cex=0.5)

pdf("ma_EmptyDrops.pdf")
par(mar=c(5.1, 4.1, 2.1, 5.1))
plot(A, M, xlab="A", ylab="M (ambient - cells)", cex.main=1.4, main="EmptyDrops",
     cex.axis=1.2, cex.lab=1.4, pch=16, cex=0.2, type="n")
limits <- loadRaster(tmp)

# Choosing the top 10 DE genes to show, by fitting a Poisson GLM.
library(edgeR)
y <- DGEList(summed)
retained <- which(aveLogCPM(y) > 1)
y <- y[retained,]
y <- calcNormFactors(y)
design <- cbind(1, c(0, 1))
fit <- glmFit(y, design, dispersion=0) 
res <- glmLRT(fit)
o <- order(res$table$PValue)
#差异表达最显著的基因（索引）
extreme <- retained[head(o, 10)]

ypos <- rank(M[extreme]) / length(extreme) * diff(range(M)) + min(M)
xpos <- limits[2] + 0.05 * diff(limits)
#text(xpos, ypos, col="red", pos=4, xpd=TRUE, offset=0.1,rowData(raw_GEX)$Symbol[match(names(M)[extreme], rownames(raw_GEX))])
matched_symbols <- names(M)[extreme]
text(xpos, ypos, col="red", pos=4, xpd=TRUE, offset=0.1, matched_symbols)
segments(xpos, ypos, A[extreme], M[extreme], col="red", xpd=TRUE)
box()
dev.off()











