install.packages("spatstat.utils")
install.packages("fastmap")
install.packages("Seurat")
install.packages("ggsci")
library(tidyverse)
library(Matrix)
library(Seurat)
library(future)
library(dplyr)
library(knitr)
library(ggsci)

setwd("C:/Users/Administrator/Desktop/SAMP6/")

#数据地址
data_path_SAMP6 =  "C:/Users/Administrator/Desktop/SAMP6/"
data_path_SAMR1 =  "C:/Users/Administrator/Desktop/SAMR1/"

#输出图像地址
output_path = "C:/Users/Administrator/Desktop/SAMP6/"

gene = readr::read_tsv(paste0(data_path_SAMR1,"genes.tsv"), col_names = FALSE)

data_6 = Matrix::readMM(paste0(data_path_SAMP6,"matrix.mtx") )
data_1 = Matrix::readMM(paste0(data_path_SAMR1,"matrix.mtx") )

data = cbind(data_1,data_6)

cellnames = c( paste0( "cell",1:ncol(data_1),"-","m1"),paste0( "cell",1:ncol(data_6),"-","m6") )
genenames = paste0(gene$X2,"-",gene$X1)

rownames(data) = genenames
colnames(data) = cellnames

Sobj = CreateSeuratObject(counts = data)
Sobjr1 = CreateSeuratObject(counts = data_1)

groups <- c(rep("r1", ncol(data_1)), rep("p6", ncol(data_6)))

# 将 group 标签添加到 Seurat 对象的元数据中
Sobj@meta.data$group <- groups
#Sobj@meta.data$groups =  as.factor( c(rep("g1",ncol(data_1)),rep("g2",ncol(data_6)) ) )

metadata <- c("groups")
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

table_html <-table(Idents(Sobj)) %>%  kable("html",table.attr = "style='width:100%;'")

table_p6r1 <-head(Sobj@meta.data) %>%  kable("html",table.attr = "style='width:100%;'")

getwd()
writeLines(table_html, "table.html")
writeLines(table_p6r1, "table2.html")




library(GGally)
## Registered S3 method overwritten by 'GGally':
##   method from
##   +.gg   ggplot2
# 'seq_folder', 'nUMI', 'nGene', 'groups', 'percent.mt', 'percent.HB', 'log10GenesPerUMI', 'mitoRatio', 'cells', 'sample'
f1 = ggpairs(Sobj@meta.data, mapping = aes(color= groups),columns = c("nFeature_RNA", "nCount_RNA") )
# 'seq_folder', 'nUMI', 'nGene', 'groups', 'percent.mt', 'percent.HB', 'log10GenesPerUMI', 'mitoRatio', 'cells', 'sample'
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

print(f1)


#log1p( (Sobj@assays$RNA@data[,1]/Sobj$nCount_RNA[1])*10000 ))
Sobj <- NormalizeData(Sobj, normalization.method = "LogNormalize", scale.factor = 10000)
#Identification of highly variable features (feature selection)
Sobj <- FindVariableFeatures(Sobj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Sobj), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Sobj)
plot(plot1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)
CombinePlots(plots = list(plot1, plot2),legend="bottom")

?PercentageFeatureSet
Sobj[[]]
Sobj[["percent.mt"]] <- PercentageFeatureSet(Sobj, pattern = "^mt-")
Hb.genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz")
Hb_m <- match(Hb.genes, rownames(Sobj@assays$RNA)) 
Hb.genes <- rownames(Sobj@assays$RNA)[Hb_m] 
Hb.genes <- Hb.genes[!is.na(Hb.genes)] 
Sobj[["percent.Hb"]]<-PercentageFeatureSet(Sobj, features=Hb.genes) 
#head(scRNA@meta.data)
col.num <- length(levels(Sobj@active.ident))
violin <- VlnPlot(Sobj,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
###把图片画到画板上面
plot(violin)
#####以后保存图片都手动保存 不要用代码保存了。
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
###这几个指标之间的相关性。 把图画到画板上，然后手动保存
plot1=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
plot1
plot2
plot3
pearplot

# Scaling the data
#all.genes <- rownames(Sobj)
# Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito.
Sobj <- ScaleData(Sobj, features = VariableFeatures(Sobj),vars.to.regress=c('nCount_RNA'))
Sobj <- RunPCA(Sobj, features = VariableFeatures(object = Sobj ))


table_ordinatemerge <-head(Sobj@reductions$pca@cell.embeddings) %>% knitr::kable(caption = " pca  cell.embeddings")
writeLines(table_ordinatemerge, "table_ordinatemerge.excel")

head(Sobj@reductions$pca@feature.loadings)  %>% knitr::kable(caption = " pca  feature.loadings")
head(Sobj@reductions$pca@assay.used)
head(Sobj@reductions$pca@stdev)
head(Sobj@reductions$pca@key)
# Get the feature loadings for a given DimReduct(Loadings(object = seuo[["pca"]])[1:5,1:5])
new.loadings <- Loadings(object = Sobj[["pca"]])
new.loadings <- new.loadings + 0.01
Loadings(object = Sobj[["pca"]]) <- new.loadings
VizDimLoadings(Sobj,dims = 1:4)

# Examine and visualize PCA results a few different ways
print(Sobj[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(Sobj, reduction = "pca")
#DimHeatmap(seuo, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Sobj, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Sobj <- JackStraw(Sobj, num.replicate = 100)
Sobj <- ScoreJackStraw(Sobj, dims = 1:20)
plot1<-JackStrawPlot(Sobj, dims = 1:15)
plot2<-ElbowPlot(Sobj)
CombinePlots(plots = list(plot1, plot2),legend="bottom")

#Returns a set of genes, based on the JackStraw analysis, that have statistically significant associations with a set of PCs.
#?PCASigGenes
head(PCASigGenes(Sobj,pcs.use=2,pval.cut = 0.7))


#聚类分析
# Cluster the cells
#Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
Sobj <- FindNeighbors(Sobj, dims = 1:18)

Sobj <- FindClusters(
  object = Sobj,
  resolution = c(seq(0,2,.2))
)


library(clustree)
clustree(Sobj@meta.data, prefix = "RNA_snn_res.")


Idents(Sobj)  <- 'RNA_snn_res.0.8'  # 指定聚类 resolution。
table(Sobj@active.ident)   %>% knitr::kable(caption = "  cluster  cell number ") # 查看每一类有多少个细胞

qcvl<- c("nCount_RNA")

for( i in qcvl){
  cat(sprintf('\n###  VlnPlot   %s    \n', i ))
  print(VlnPlot(   Sobj , features = i ,pt.size= 0 , ncol = 1))
  cat(sprintf('\n\n'))
}

filename = paste0(output_path,"ftmp.png")
ggsave(filename, width = 10, height = 20,  units = "cm")


library(kableExtra)
cat(sprintf('\n\n### cluster    umap embedding   \n\n' ))

print(DimPlot(Sobj,label= T) +NoLegend())
table_cluster<- table(Idents(Sobj))  %>%  kable("html",table.attr = "style='width:100%;'")
writeLines(table_cluster, "table_cluster.html")

table_cluster2 <- prop.table(table(Idents(Sobj))) %>% kable("html", table.attr = "style='width:100%;'")
writeLines(table_cluster2, "table_cluster2.html")

pdf <- prop.table(table(Idents(Sobj)))
td <- as.data.frame(pdf)

allcolour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
Sobj@meta.data$groups =  as.factor( c(rep("samr1",ncol(data_1)),rep("samp6",ncol(data_6)) ) )

#在计算samp6和samr1里27个群分别的占比
as.data.frame(prop.table(table(Idents(Sobj), Sobj@meta.data[,"groups"]), margin = 2))-> pdf -> td
plt <- ggplot(td, aes(x = td[,2], y = Freq, fill = td[,1])) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "sample", y = "Cells Ratio") +
  theme(panel.background = element_rect(fill = "transparent", color = "black"),
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        axis.text = element_text(color = "black")) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_manual(values = allcolour) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol = 1, title = "Cluster"))

plt + coord_flip()


# 统计表格
table(subsobj$target, subsobj$groups) %>% kable("html", table.attr = "style='width:100%;'")
prop.table(table(subsobj$target, subsobj$groups), margin = 2) %>% kable("html", table.attr = "style='width:100%;'")

# 计算百分比并转换为数据框
td <- as.data.frame(prop.table(table(subsobj$target, subsobj@meta.data[,"groups"]), margin = 2))

# 生成柱状图并添加百分比标签
plt <- ggplot(td, aes(x = td[,2], y = Freq, fill = td[,1])) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = scales::percent(Freq)), position = position_stack(vjust = 0.5), size = 3) + # 添加百分比标签
  labs(x = "sample", y = "Cells Ratio") +
  theme(
    panel.background = element_rect(fill = "transparent", color = "black"),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 14) # 修改这里来增大SAMP6, SAMR1的字号
  ) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_manual(values = allcolour) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol = 1, title = "Cluster"))
plt + coord_flip()

library(ggforce)
p1 <- Sobj@meta.data %>%
  count(RNA_snn_res.0.8) %>%
  mutate(focus = ifelse(RNA_snn_res.0.8 == "0", 0.2, 0)) %>%
  ggplot()+
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1, amount = n,
                   fill = RNA_snn_res.0.8, explode = focus),
               alpha = 1, stat = "pie") +
  scale_fill_manual(values = allcolour) +
  theme_bw()

p1

install.packages('ape')
#Constructs a phylogenetic tree relating the 'average' cell from each identity class.
# Tree is estimated based on a distance matrix constructed in either gene expression space or PCA spac

# ?BuildClusterTree
Sobj<-BuildClusterTree(Sobj)
Tool(object = Sobj, slot = 'BuildClusterTree')
PlotClusterTree(Sobj)

#Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
#绘制UMAP图
head(Sobj)
Sobj <- RunUMAP(Sobj, dims = 1:18)

head(Sobj@reductions$umap@cell.embeddings) # 提取UMAP坐标值。


plot1<-DimPlot(Sobj, reduction = "umap",label = TRUE, label.size = 6.5)
plot1 <- plot1 + theme(legend.text = element_text(size = 18))  # 修改 size 为你想要的大小

plot(plot1)

plot1 <- plot1 + theme(
  plot.title = element_text(size = 20), # 调整标题的字体大小
  axis.title = element_text(size = 20), # 调整轴标题的字体大小
  axis.text = element_text(size = 20),   # 调整轴标签的字体大小
  legend.title = element_text(size = 20), # 调整图例标题的字体大小
  legend.text = element_text(size = 25)   # 调整图例文本的字体大小
)
plot(plot1)
anno <- read.csv("anno0724.csv", header= F)


head(anno)
dim(anno)

cluster.ids = as.array(anno[,2])
view(cluster.ids)
names(cluster.ids) <- levels(Sobj)
levels(Sobj)
Sobj <- RenameIdents(Sobj, cluster.ids)
Idents(Sobj)
head(Sobj)
Sobj[['celltype']] = Sobj@active.ident #添加'cell_types'信息

Sobj@misc$markers <- FindAllMarkers(Sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 将marker基因数据输出为CSV文件
write.csv(Sobj@misc$markers, file = "Sobj_markers.csv", row.names = FALSE)


plot3<-FeaturePlot(Sobj, features = c("Ly6g-ENSMUSG00000022582"),min.cutoff = 0, max.cutoff = 4)
plot(plot3)
plot4<-FeaturePlot(Sobj, features = c("Cd79a-ENSMUSG00000003379"),min.cutoff = 0, max.cutoff = 4)
plot(plot4)
plot4<-FeaturePlot(Sobj, features = c("Cd19-ENSMUSG00000030724"),min.cutoff = 0, max.cutoff = 4)
plot(plot4)
plot4<-FeaturePlot(Sobj, features = c("Tfrc-ENSMUSG00000022797"),min.cutoff = 0, max.cutoff = 4)
plot(plot4)
plot5<-FeaturePlot(Sobj, features = c("Csf1r-ENSMUSG00000024621"),min.cutoff = 0, max.cutoff = 4)
plot(plot5)
plot6<-FeaturePlot(Sobj, features = c("Kit-ENSMUSG00000005672"),min.cutoff = 0, max.cutoff = 4)
plot(plot6)
plot7<-FeaturePlot(Sobj, features = c("Cd3e-ENSMUSG00000032093"),min.cutoff = 0, max.cutoff = 4)
plot(plot7)
plot8<-FeaturePlot(Sobj, features = c("Mpo-ENSMUSG00000009350"),min.cutoff = 0, max.cutoff = 4)
plot(plot8)
plot9<-FeaturePlot(Sobj, features = c("Klrb1c-ENSMUSG00000030325"),min.cutoff = 0, max.cutoff = 4)
plot(plot9)
plot10<-FeaturePlot(Sobj, features = c("Siglech-ENSMUSG00000051504"),min.cutoff = 0, max.cutoff = 4)
plot(plot10)
plot11<-FeaturePlot(Sobj, features = c("Hbb-bt-ENSMUSG00000073940"),min.cutoff = 0, max.cutoff = 4)
plot(plot11)
plot12<-FeaturePlot(Sobj, features = c("Sparc-ENSMUSG00000018593"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Pf4-ENSMUSG00000029373"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Jchain-ENSMUSG00000067149"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Prtn3-ENSMUSG00000057729"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Vcam1-ENSMUSG00000027962"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Cd68-ENSMUSG00000018774"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Ms4a2-ENSMUSG00000024680"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Cd34-ENSMUSG00000016494"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Ms4a2-ENSMUSG00000024680"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Cd79b-ENSMUSG00000040592"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)
plot12<-FeaturePlot(Sobj, features = c("Col1a2-ENSMUSG00000029661"),min.cutoff = 0, max.cutoff = 4)
plot(plot12)


#基因表达点图
library(Seurat)
DotPlot(Sobj, features= c("Cd79a-ENSMUSG00000003379", "Cd19-ENSMUSG00000030724", "Ly6g-ENSMUSG00000022582", "Cd177-ENSMUSG00000052212", "Siglech-ENSMUSG00000051504","Mpo-ENSMUSG00000009350", "Prtn3-ENSMUSG00000057729","Csf1r-ENSMUSG00000024621","Cd68-ENSMUSG00000018774","Kit-ENSMUSG00000005672","Cd3e-ENSMUSG00000032093","Klrb1c-ENSMUSG00000030325","Hbb-bt-ENSMUSG00000073940","Cd79b-ENSMUSG00000040592","Jchain-ENSMUSG00000067149","Pf4-ENSMUSG00000029373"))

# 加载Seurat包
library(Seurat)

# 定义特征名称和显示名称的向量
features <- c("Hbb-bt-ENSMUSG00000073940","Klrb1c-ENSMUSG00000030325",
              "Ms4a2-ENSMUSG00000024680", "Siglech-ENSMUSG00000051504","Cd79b-ENSMUSG00000040592",
              "Cd3e-ENSMUSG00000032093","Kit-ENSMUSG00000005672", 
              "Mpo-ENSMUSG00000009350","Prtn3-ENSMUSG00000057729","Cd79a-ENSMUSG00000003379", "Cd19-ENSMUSG00000030724", 
              "Ly6g-ENSMUSG00000022582", "Cd177-ENSMUSG00000052212", 
              "Csf1r-ENSMUSG00000024621" )

# 创建一个映射显示名称的向量
display_names <- c("Hbb-bt","Klrb1c",
                   "Ms4a2", "Siglech","Cd79b",
                   "Cd3e","Kit", 
                   "Mpo","Prtn3","Cd79a", "Cd19", 
                   "Ly6g", "Cd177", 
                   "Csf1r" )


# 绘制DotPlot并替换标签
# 使用更深的绿色 (自定义颜色代码)
DotPlot(Sobj, features = features) + 
  scale_x_discrete(labels = setNames(display_names, features)) + 
  RotatedAxis() +
  scale_color_gradient(low = "white", high = "#21a675") +  # 使用自定义颜色代码
  theme(legend.position = "right")  # 可选：调整图例的位置

top6 <- Sobj@misc$markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)

Sobjheat <- ScaleData(Sobj,features = top6$gene)
head(Sobjheat)

library(scales)
plot15 <- DoHeatmap(subset(Sobjheat, downsample = 500), features = top6$gene) +
  scale_fill_gradientn(colors = c("blue","white", "red")) + NoLegend() +
  theme(axis.text.x = element_blank(), # 去掉横轴标签
        axis.title.x = element_blank(), # 去掉横轴标题
        axis.text.y = element_text(size = 30)) # 调整纵轴标签的字体大小
plot15 <- DoHeatmap(subset(Sobjheat, downsample = 500), features = top7$gene) +
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    limits = c(-3, 3), # 设置显示范围
    values = rescale(c(-3, 0, 3)) # 确保0是白色
  ) +
  NoLegend() +
  theme(
    axis.text.x = element_blank(), # 去掉横轴标签
    axis.title.x = element_blank(), # 去掉横轴标题
    axis.text.y = element_text(size = 20) # 调整纵轴标签的字体大小
  )

# 添加图例作为横轴标签，并调整图例文字和点的大小
plot15 <- plot15 + 
  theme(
    legend.position = "right", # 图例放在右侧
    legend.title = element_blank(), # 不显示图例标题
    legend.text = element_text(size = 23) # 调整图例文字的字体大小
  ) 
plot(plot15)

