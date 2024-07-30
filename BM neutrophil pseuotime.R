library(slingshot)
library(tradeSeq)
library(BiocParallel)
sce <- as.SingleCellExperiment(obj_big, assay = "RNA") # SCT
sce_slingshot1 <- slingshot(sce, 
                            reducedDim = 'UMAP', 
                            clusterLabels = sce$celltype) 
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_slingshot1$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce_slingshot1)$UMAP, col = plotcol, pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sce_slingshot1), lwd=2, col='black')

###基因变化
slingsce<-SlingshotDataSet(sce_slingshot1)
pseudotimeED <- slingPseudotime(slingsce, na = FALSE)
cellWeightsED <- slingCurveWeights(slingsce)
counts<-sce_slingshot1@assays@data@listData$counts
# 拟合负二项式模型确定nknots
icMat <- evaluateK(counts = counts, 
                   sds = slingsce, 
                   nGenes = 50,
                   k = 3:10,
                   verbose = T, 
                   plot = TRUE)
sce_slinghot <- fitGAM(counts = counts, pseudotime = pseudotimeED, cellWeights = cellWeightsED, nknots = 5, verbose = T
                       ,BPPARAM = MulticoreParam(20),parallel=T)
# 与伪时间存在显著关联的基因
rowData(sce_slinghot)$assocRes <- associationTest(sce_slinghot, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce_slinghot)$assocRes
oAsso <- order(assocRes$waldStat, decreasing = TRUE)
Geneasso <- names(sce_slinghot)[oStart[1]]
plotSmoothers(sce_slinghot, assays(sce_slinghot)$counts,
              gene =Geneasso, alpha = 0.6, border = T, lwd = 2)+
  ggtitle(Geneasso)
# 谱系内比较
startRes <- startVsEndTest(sce_slinghot,lineages=T)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce_slinghot)[oStart[1]]
plotSmoothers(sce_slinghot, assays(sce_slinghot)$counts, gene = sigGeneStart)+
  ggtitle(sigGeneStart)
plotGeneCount(slingsce, assays(sce_slinghot)$counts, gene = sigGeneStart)+
  ggtitle(sigGeneStart)
# 谱系间比较
endRes <- diffEndTest(sce_slinghot,pairwise=T)
oEnd <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce_slinghot)[o[1]]
plotSmoothers(sce_slinghot, assays(sce_slinghot)$counts, gene = sigGene)+
  ggtitle(sigGene)
plotGeneCount(slingsce, assays(sce_slinghot)$counts, gene = sigGene)+
  ggtitle(sigGene)

