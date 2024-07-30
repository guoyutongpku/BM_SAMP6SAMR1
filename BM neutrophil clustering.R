# How does cluster membership vary by replicate?
table(Idents(subsobj), subsobj$groups)  %>%  kable("html",table.attr = "style='width:100%;'")
prop.table(table(Idents(subsobj),  subsobj$groups), margin = 2) %>%  kable("html",table.attr = "style='width:100%;'")
#把上面的表格频率改成柱状图
as.data.frame(prop.table(table(Idents(subsobj), subsobj@meta.data[,"groups"]), margin = 2))-> pdf -> td
plt <- ggplot(td, aes(x = td[,2], y = Freq, fill = td[,1])) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "sample", y = "Cells Ratio") +
  theme(panel.background = element_rect(fill = "transparent", color = "black"),
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        axis.text = element_text(color = "black")) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_manual(values = allcolour) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol = 1, title = "Cluster"))
plt +   coord_flip()

library(ggplot2)

# Assuming `td` is a data frame created from the previous steps
as.data.frame(prop.table(table(Idents(subsobj), subsobj@meta.data[,"groups"]), margin = 2))-> pdf -> td
plt <- ggplot(td, aes(x = td[,2], y = Freq, fill = td[,1])) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "sample", y = "Cells Ratio") +
  theme(
    panel.background = element_rect(fill = "transparent", color = "black"),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(color = "black", size = 14),   # Change axis text size
    axis.title = element_text(size = 30),                  # Change axis title size
    legend.title = element_text(size = 20),                # Change legend title size
    legend.text = element_text(size = 20),                 # Change legend text size
    plot.title = element_text(size = 20, face = "bold")    # Optional: Change plot title size
  ) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_manual(values = allcolour) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol = 1, title = "Cluster"))

plt + coord_flip()

#画subneu的热图
# find markers for every cluster compared to all remaining cells, report only the positive ones
subsobj@misc$markers <- FindAllMarkers(subsobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(subsobj@misc$markers)
subsobj@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  %>%  kable() 
table_subneu <- subsobj@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  %>%  kable() 
write_xlsx(results_neutop10dsg, "neutop10dsg.xlsx")
writeLines(table_subneu, "table_subneu.html")


top10sub <- subsobj@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
valid_genes <- top10sub$gene %in% rownames(subsobj)

# 过滤掉无效基因
filtered_genes <- top10sub$gene[valid_genes]

# 打印无效基因
if (any(!valid_genes)) {
  cat("以下基因不存在于数据集中:\n")
  print(top10sub$gene[!valid_genes])
}


# 缩放数据
sobjneu <- ScaleData(subsobj, features = filtered_genes)
head(sobjneu)


sobjneu <- ScaleData(subsobj,features = top10sub$gene)
plot16 <-DoHeatmap(subset(sobjneu, downsample = 500), features = top10$gene) + NoLegend()
+
  NoLegend() +
  theme(
    axis.text.x = element_blank(), # 去掉横轴标签
    axis.title.x = element_blank(), # 去掉横轴标题
    axis.text.y = element_text(size = 20) # 调整纵轴标签的字体大小
  )

# 添加图例作为横轴标签，并调整图例文字和点的大小
plot15 <- plot16 + 
  theme(
    legend.position = "right", # 图例放在右侧
    legend.title = element_blank(), # 不显示图例标题
    legend.text = element_text(size = 23) # 调整图例文字的字体大小
  ) 

VlnPlot(subset(subsobj), features = c('Cd55'),group.by = "target",slot = "data")


plot1 <- FeaturePlot(subsobj_samp6, features = c('Cd55'), min.cutoff = 0, max.cutoff = 4)
plot(plot1)


