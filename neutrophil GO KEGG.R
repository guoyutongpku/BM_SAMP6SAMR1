

library("AnnotationDbi")
library("org.Mm.eg.db")
install.packages("clusterProfiler")
library(clusterProfiler)
# 提取 gene 列中 "-" 之前的部分作为行名
DEGs$genename <- gsub("-.*", "", DEGs$gene)
# BP, CC和MF三种通路都一起富集
subsobjgo <-   enrichGO(gene  = DEGs$genename,
                        #universe     = row.names(dge.celltype),
                        OrgDb         = 'org.Mm.eg.db',
                        keyType       = 'SYMBOL',
                        ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)

cellcluster_msc = unique(DEGs$group1)

write.csv(subsobjgo, file = "./GO/go.csv", row.names = T)
dotplot(subsobjgo, showCategory = 20)  

# 假设 subsobjgo 已经通过 enrichGO 函数得到了 GO 富集分析的结果

# 从 DEGs 中筛选出 group1 列为 "C10" 的基因
genes_of_interest <- DEGs$genename[DEGs$group1 == "C10"]

# 进行 GO 富集分析
subsobjgo_c10 <- enrichGO(gene = genes_of_interest,
                          OrgDb = 'org.Mm.eg.db',
                          keyType = 'SYMBOL',
                          ont = "ALL",            # 设置为ALL时BP, CC, MF都计算
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05)
subsobjgo_c10bp <- enrichGO(gene = genes_of_interest,
                            OrgDb = 'org.Mm.eg.db',
                            keyType = 'SYMBOL',
                            ont = "BP",            # 设置为ALL时BP, CC, MF都计算
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.05)

top_go_datac10 <- subsobjgo_c10[order(-go_data$Count), ][1:15, ]

# 创建dotplot图
dotplot <- ggplot(top_go_datac10, aes(x = Count, y = reorder(Description, Count))) +
  geom_point(aes(size = -log10(pvalue)), color = 'red', alpha = 0.9) +
  scale_size_area(max_size = 10) +
  labs(x = 'Gene Count', y = 'GO Term', title = 'Top 15 GO Analysis Dotplot for C10', size = 'p-value (-log10)') +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 16))

# 反转y轴顺序以使基因数最大的GO Term显示在顶部
dotplot <- dotplot + scale_y_discrete(limits = rev(levels(reorder(top_go_datac10$Description, top_go_datac10$Count))))

# 显示图表
print(dotplot)


# 可以进一步对 subsobjgo_c10 进行可视化或其他分析
top_go_datac10 <- subsobjgo_c10[order(-go_data$Count), ][1:20, ]

# 创建dotplot图
dotplot <- ggplot(subsobjgo_c10, aes(x = Count, y = reorder(Description, Count))) +
  geom_point(aes(size = -log10(pvalue)), color = 'red', alpha = 0.9) +
  scale_size_area(max_size = 10) +
  labs(x = 'Gene Count', y = 'GO Term', title = 'Top 20 GO Analysis Dotplot in C10', size = 'p-value (-log10)') +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 18))

# 反转y轴顺序以使基因数最大的GO Term显示在顶部
dotplot <- dotplot + scale_y_discrete(limits = rev(levels(reorder(subsobjgo_c10$Description, top_go_data$Count))))

# 显示图表
print(dotplot)

##KEGG分析
genes_of_interest <- DEGs$genename[DEGs$group1 == "C10"]

subsobjkegg_c10 <- lapply(DEGs, function(x){
  bitr(gene = genes_of_interest, fromType="SYMBOL",
       toType="ENTREZID", OrgDb='org.Mm.eg.db') %>% pull(ENTREZID) %>% enrichKEGG(organism = 'mmu',pAdjustMethod = 'fdr') #hsa人；mmu小鼠
})

write.csv(subsobjkegg_c10, file = "./GO/subsobjkeggc10.csv", row.names = T)
plots <- lapply(subsobjkegg_c10, function(result){
  if (class(result) == "enrichResult") 
    dotplot(result, showCategory=10, title="KEGG Pathway Enrichment")
  
})

# 显示所有生成的图
for (plot in plots) {
  print(plot)
}



