
# F7: Activation and cytotoxic function of hepatic NK cells are impaired during clinical ketosis.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 NK细胞亚群数据处理
# =============================================

# 加载R包
library(Seurat)
library(Rcpp)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
library(BPCells)
library(presto)
library(glmGamPoi)
library(data.table)
library(tidyverse)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(ggsci)
library(cowplot)
library(ggplotify)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(corrplot)
library(rio)
library(dplyr)
library(ggsignif)
library(DESeq2)
library(RColorBrewer)

# 加载数据
load("CD45去双细胞后.Rdata")  
  
# 提取出里面的NK细胞
NK <- subset(seurat_object, celltype %in% c("NK"))

# 标准化 
NK <- NormalizeData(NK, normalization.method = "LogNormalize", scale.factor = 10000) 

# 高变基因
NK <- FindVariableFeatures(NK, selection.method = "vst", nfeatures = 2000) 

# 周期打分
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene) 
X = homologene(s.genes,inTax = 9606,outTax = 9913) 
Y = homologene(g2m.genes,inTax = 9606,outTax = 9913) 
s.genes = X$"9913"
g2m.genes = Y$"9913" 
NK <- CellCycleScoring(NK, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) 

# 归一化缩放去除周期影响
NK <- ScaleData(NK, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(NK))

# 线性降维PCA 默认用高变基因集
NK <- RunPCA(NK, features = VariableFeatures(object = NK))

# 肘部图
ElbowPlot(NK, 50)

# 计算KNN和SNN
NK = FindNeighbors(NK, dims = 1:15) 

# 分群数量
NK = FindClusters(NK, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(NK, prefix = "RNA_snn_res.") 
Idents(NK) <- "RNA_snn_res.0.1"
NK$seurat_clusters <- NK@active.ident

# UMAP非线性降维
NK <- RunUMAP(NK, dims = 1:25) 
save(NK, file = "NK降维聚类.Rdata") 
rm(list = ls())
load("NK降维聚类.Rdata")  
NK = RunHarmony(NK, "orig.ident", plot_convergence = TRUE)   
NK = FindNeighbors(NK, reduction = "harmony", dims = 1:15)
NK <- RunUMAP(NK, reduction = "harmony", dims = 1:30)   
NK = FindClusters(NK, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(NK, prefix = "RNA_snn_res.")
Idents(NK) <- "RNA_snn_res.1"
NK$seurat_clusters <- NK@active.ident  
NK <- subset(NK, seurat_clusters %in% c("12"), invert = TRUE) 
  
# UMAP图
Seurat::DimPlot(NK, reduction = "umap", pt.size = 0.1, label = T)  
  
# 需要重新走降维聚类流程
NK = FindNeighbors(NK, reduction = "harmony", dims = 1:15)  

# 批量设置分辨率
NK = FindClusters(NK, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(NK, prefix = "RNA_snn_res.") 
Idents(NK) <- "RNA_snn_res.0.4"
NK$seurat_clusters <- NK@active.ident  
  
# UMAP非线性降维
NK <- RunUMAP(NK, reduction = "harmony", dims = 1:25)  
NK <- subset(NK, seurat_clusters %in% c("5", "6"), invert = TRUE) 
Seurat::DimPlot(NK, reduction = "umap", pt.size = 0.1, label = T)  
NK = FindNeighbors(NK, reduction = "harmony", dims = 1:15)  
NK = FindClusters(NK, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(NK, prefix = "RNA_snn_res.") 
Idents(NK) <- "RNA_snn_res.0.3"
NK$seurat_clusters <- NK@active.ident  
NK <- RunUMAP(NK, reduction = "harmony", dims = 1:30)  
  
# UMAP图
Seurat::DimPlot(NK, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(NK, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(NK, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(NK, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)  

# 未注释前细胞cluster进行差异基因分析 用来鉴定细胞
dif<-FindAllMarkers(NK, 
                    group.by = NK@meta.data$seurat_clusters, 
                    logfc.threshold = log2(1.2),                        
                    min.pct = 0.2,                                      
                    only.pos = T                                        
                    )        
dif$pct_diff <- dif$pct.1 - dif$pct.2 
table(dif$cluster)                    
dif<-dif %>%
  group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
save(dif, file = "NK整体cluster差异基因.Rdata")

# 为分群重新指定细胞类型 
new.cluster.ids <- c("GNLY_mNK",
                     "EOMES_mNK",
                     "ZNF683_NK",
                     "IL7R_iNK") 
new.cluster.ids
names(new.cluster.ids) 
levels(NK)
names(new.cluster.ids) <- levels(NK) 
names(new.cluster.ids)
new.cluster.ids
NK <- RenameIdents(NK, new.cluster.ids) 

NK[["celltype"]] <- Idents(NK) 
table(NK@meta.data[["orig.ident"]])
unique(NK@meta.data[["orig.ident"]])
NK[["group"]]<- c(rep("Health", 1501), rep("Ketosis", 559))
NK[["group.celltype"]]<-paste(NK$group, Idents(NK), sep = '_') 
table(NK@meta.data[["orig.ident"]])
table(NK@meta.data[["group"]])
table(NK@meta.data[["celltype"]])
table(NK@meta.data[["group.celltype"]])
table(NK@meta.data[["seurat_clusters"]])

# 细胞水平信息
Idents(NK) <- factor(Idents(NK),
                       levels = c("GNLY_mNK", "EOMES_mNK", "ZNF683_NK", "IL7R_iNK"))
NK[["celltype"]] <- Idents(NK) 

# 绘制总umap图
Seurat::DimPlot(NK, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(NK, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(NK, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(NK, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 保存工作空间
save(NK,file = "NK细胞鉴定.Rdata")
rm(list = ls())
load("NK细胞鉴定.Rdata")


# =============================================
# 🎨 Fig. 7f
# =============================================

# 样本点图 以CD45总细胞为背景进行计算
# 计算每个细胞群体在总细胞中的比例
load("NK细胞鉴定.Rdata")
load("CD45去双细胞后.Rdata")
table(seurat_object$orig.ident)
table(Idents(NK), NK$orig.ident)

# 获取每个亚群在每个样本中的细胞数
NK_subgroup_counts <- table(Idents(NK), NK$orig.ident)

# 获取每个样本中的细胞数
total_cells <- table(seurat_object$orig.ident)
total_cells <- total_cells[match(colnames(NK_subgroup_counts), names(total_cells))]
total_cells_matrix <- matrix(total_cells, nrow = nrow(NK_subgroup_counts), ncol = ncol(NK_subgroup_counts), byrow = TRUE)

# 计算每个亚群在每个样本中的占比
NK_subgroup_percentage <- NK_subgroup_counts / total_cells_matrix * 100

# 查看比例
NK_subgroup_percentage

# 转成数据框
Cellratio <- as.data.frame(as.table(NK_subgroup_percentage))

# 修改列名
colnames(Cellratio) <- c("Celltype", "orig.ident", "Freq")

# 创建分组变量 H / K
Cellratio$group <- Cellratio$orig.ident

# 去掉最后一个字符，只保留 H 或 K
Cellratio$group <- gsub(".{1}$", "", Cellratio$group)

# 设置因子顺序
Cellratio$group <- factor(Cellratio$group, levels = c("H", "K"))

# 查看结果
Cellratio

# GNLY_mNK 频率分布图
set.seed(123)
ggplot(data=Cellratio[Cellratio$Celltype=="GNLY_mNK",], aes(group, Freq)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(0, 4.8) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
        legend.position = "none") +
  stat_summary(fun = mean,
               geom = "errorbar",
               size = 0.5,
               width = 0.15,
               color = "black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) + 
  # 使用平均值mean
  stat_summary(fun = mean,         
               geom = "crossbar",
               size = 0.2,        
               color = "black",  
               width = 0.2)        

Cellratio[Cellratio$Celltype == "GNLY_mNK",]
shapiro.test(Cellratio[Cellratio$Celltype=="GNLY_mNK",]$Freq[1:3]) 
shapiro.test(Cellratio[Cellratio$Celltype=="GNLY_mNK",]$Freq[4:8]) 
leveneTest(Freq ~ group, data = Cellratio[Cellratio$Celltype=="GNLY_mNK",], center = "mean") 
t.test(Freq ~ group, data = Cellratio[Cellratio$Celltype=="GNLY_mNK",], var.equal = TRUE) # 0.0001214 独立两组-正态-方差齐性

# EOMES_mNK 频率分布图
set.seed(124)
ggplot(data=Cellratio[Cellratio$Celltype=="EOMES_mNK",], aes(group, Freq)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(-0.08, 2.25) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
        legend.position = "none") +
  stat_summary(fun = mean,
               geom = "errorbar",
               size = 0.5,
               width = 0.15,
               color = "black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean,          
               geom = "crossbar",   
               size = 0.2,       
               color = "black",  
               width = 0.2)      

Cellratio[Cellratio$Celltype == "EOMES_mNK",]
shapiro.test(Cellratio[Cellratio$Celltype=="EOMES_mNK",]$Freq[1:3]) 
shapiro.test(Cellratio[Cellratio$Celltype=="EOMES_mNK",]$Freq[4:8]) 
leveneTest(Freq ~ group, data = Cellratio[Cellratio$Celltype=="EOMES_mNK",], center = "mean") 
t.test(Freq ~ group, data = Cellratio[Cellratio$Celltype=="EOMES_mNK",], var.equal = TRUE) # 0.08687 独立两组-正态-方差齐性

# EOMES_mNK 频率分布图
set.seed(124)
ggplot(data=Cellratio[Cellratio$Celltype=="ZNF683_NK",], aes(group, Freq)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(0, 0.28) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
        legend.position = "none") +
  stat_summary(fun = mean,
               geom = "errorbar",
               size = 0.5,
               width = 0.15,
               color = "black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) + 
  stat_summary(fun = mean,         
               geom = "crossbar",  
               size = 0.2,         
               color = "black",  
               width = 0.2)     

Cellratio[Cellratio$Celltype == "ZNF683_NK",]
shapiro.test(Cellratio[Cellratio$Celltype=="ZNF683_NK",]$Freq[1:3]) 
shapiro.test(Cellratio[Cellratio$Celltype=="ZNF683_NK",]$Freq[4:8]) 
leveneTest(Freq ~ group, data = Cellratio[Cellratio$Celltype=="ZNF683_NK",], center = "mean")
t.test(Freq ~ group, data = Cellratio[Cellratio$Celltype=="ZNF683_NK",], var.equal = TRUE) # 0.01957 独立两组-正态-方差齐性

# IL7R_iNK 频率分布图
set.seed(124)
ggplot(data=Cellratio[Cellratio$Celltype=="IL7R_iNK",], aes(group, Freq)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(0, 0.353) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
        legend.position = "none") +
  stat_summary(fun = mean,
               geom = "errorbar",
               size = 0.5,
               width = 0.15,
               color = "black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean,          
               geom = "crossbar",   
               size = 0.2,         
               color = "black",    
               width = 0.2)       

Cellratio[Cellratio$Celltype == "IL7R_iNK",]
shapiro.test(Cellratio[Cellratio$Celltype=="IL7R_iNK",]$Freq[1:3]) 
shapiro.test(Cellratio[Cellratio$Celltype=="IL7R_iNK",]$Freq[4:8]) 
leveneTest(Freq ~ group, data = Cellratio[Cellratio$Celltype=="IL7R_iNK",], center = "mean") 
t.test(Freq ~ group, data = Cellratio[Cellratio$Celltype=="IL7R_iNK",], var.equal = TRUE) # 0.2278 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 7g
# =============================================
  
# 分组的细胞比例图 
load("NK细胞鉴定.Rdata")  
library(dplyr)
library(ggplot2)
library(gtools)
library(ggalluvial)

# 准备细胞比例输入数据
prop_df <- NK@meta.data %>%
  dplyr::select(Sample = group, Celltype = celltype) %>%
  group_by(Sample, Celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  dplyr::select(Sample, Celltype, Proportion)
prop_df

# 设置数据格式
prop_df$Celltype = factor(prop_df$Celltype, levels = unique(prop_df$Celltype))
prop_df$Sample = factor(prop_df$Sample, levels = mixedsort(unique(prop_df$Sample)))

# 开始画图
ggplot(as.data.frame(prop_df), aes(x = Sample, y = Proportion, fill = Celltype, stratum = Celltype, alluvium = Celltype)) +
  geom_flow(width = 0.6, alpha = 0.3, knot.pos = 0.1) +  
  geom_col(width = 0.6) +  
  scale_y_continuous(expand = c(0, 0)) +   
  scale_fill_manual(values = c("#ae7eb5","#f98177",
                               "#54c6a8","#568389")) +  
  xlab("") + 
  ylab("Cell proportion") +  
  theme_classic() +  
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),    
        axis.text.y = element_text(size = 16),    
        axis.title.y = element_text(size = 16),    
        legend.text=element_text(size = 16),    
        legend.title=element_text(size = 16) 
        ) + geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0,
                          fill = NA, color = 'white', size = 0.5)
  

# =============================================
# 🎨 Fig. 7a
# =============================================

load("NK细胞鉴定.Rdata")

Seurat::DimPlot(NK,
                group.by = "celltype",
                cols = c("#ae7eb5",
                         "#f98177",
                         "#54c6a8",
                         "#568389"), 
                pt.size = 0.5,
                label = T) +
  NoLegend()+ 
  labs(title = NULL)  
dev.off()  


# =============================================
# 🎨 Fig. 7b
# =============================================

load("NK整体细胞类型差异基因.Rdata")
head(diff[diff$cluster == unique(diff$cluster)[1],], 50) $ gene 
library(scplotter)
library(plotthis)
plotthis::show_palettes(type = "continuous", index = 1:30)
Seurat::DotPlot(NK, 
                features = c("GNLY","CX3CR1", "CD81", "PSTPIP2", 
                             "EOMES", "CCR5", "CD7", "CST7", 
                             "ZNF683","CCR8", "CD300LF", "XCL1", 
                             "IL7R", "CRABP1", "KIT","TNFRSF18"),
                group.by = "celltype") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid = element_blank(), 
        axis.line = element_blank()) + 
  scale_size(range = c(-0.2, 8)) +
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(11)))


# =============================================
# 🎨 Fig. 7c
# =============================================

plot_density(NK, 
             reduction = "umap",
             features = c("NCR1"),
             adjust = 1,  
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 
 
plot_density(NK, 
             reduction = "umap",
             features = c("KLRD1"),
             adjust = 1, 
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(NK, 
             reduction = "umap",
             features = c("GNLY"),
             adjust = 1,  
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(NK, 
             reduction = "umap",
             features = c("EOMES"),
             adjust = 1,  
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 
dev.off() 

plot_density(NK, 
             reduction = "umap",
             features = c("ZNF683"),
             adjust = 1, 
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(NK, 
             reduction = "umap",
             features = c("IL7R"),
             adjust = 1,  
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(NK, 
             reduction = "umap",
             features = c("CTSW"),
             adjust = 1, 
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

plot_density(NK, 
             reduction = "umap",
             features = c("NKG7"),
             adjust = 1,  
             raster = T, 
             size = 2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 


# =============================================
# 🎨 Fig. 7d
# =============================================

# 加载数据
load("NK细胞鉴定.Rdata")

# 随机抽样
table(NK@meta.data[["celltype"]])
cell_types <- unique(NK$celltype)
subset_cells <- list()  

# 对每个细胞类型进行处理
for (cell_type in cell_types) {
  cells_of_type <- WhichCells(NK, expression = celltype == cell_type) 
  if (length(cells_of_type) >= 115) {
    selected_cells <- sample(cells_of_type, 115)                    
  } else {
    selected_cells <- cells_of_type                                
  }
  subset_cells[[cell_type]] <- selected_cells                    
}  

# 合并所有子集细胞
subset_all_cells <- unlist(subset_cells)

# 根据选中的细胞创建新的 Seurat 对象
new_seurat_object <- subset(NK, cells = subset_all_cells) 

# 查看新的 Seurat 对象
new_seurat_object

# 基于上调基因分析挑选用于绘图的基因
dif<-FindAllMarkers(NK,
                    group.by = NK@active.ident,
                    logfc.threshold = log2(1.2),              
                    min.pct = 0.2,                           
                    only.pos = T)
sig.dif<-dif%>%
  group_by(cluster)%>%
  top_n(n=10,wt=avg_log2FC) 
genes<-unique(sig.dif$gene) 
genes[1:5]  

# 获取绘图使用的表达量信息 即标准化后的基因数据
data <- new_seurat_object@assays$RNA$data[genes,] 

# 加上细胞类型排序
celltype<-new_seurat_object$celltype
levels(celltype)                     
celltype[1:20]
celltype <- celltype[order(celltype)] 
celltype[1:20]
celltype <- data.frame(celltype)  
data <- data[,rownames(celltype)] 
identical(rownames(celltype),colnames(data)) 
sig.dif<-sig.dif[!duplicated(sig.dif$gene),] 
gene.anno<-data.frame(gene.anno=sig.dif$cluster,row.names = sig.dif$gene) 

# 画图
pheatmap(data1,scale = "none",cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,show_rownames = FALSE,
         annotation_col = celltype,    
         annotation_row = gene.anno,   
         annotation_names_row = FALSE, 
         color = colorRampPalette(c("#040509","#608fe4", "#ffd700"))(100)
         )


# =============================================
# 🎨 Fig. 7e
# =============================================

# 首先进行亚群间的差异分析
dif <- FindAllMarkers(NK, 
                      group.by = NK@active.ident, 
                      logfc.threshold = 0, 
                      min.pct = 0.2, 
                      only.pos = T 
                      ) 
dif$pct_diff <- dif$pct.1 - dif$pct.2
table(dif$cluster)                   
dif<-dif %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)

sigposDEG.all <- subset(dif, p_val < 0.05 & avg_log2FC > log2(1.2))
table(sigposDEG.all$cluster)

table(NK@active.ident)
GNLY_mNK  <- subset(sigposDEG.all, cluster=='GNLY_mNK') 
EOMES_mNK <- subset(sigposDEG.all, cluster=='EOMES_mNK')
ZNF683_NK <- subset(sigposDEG.all, cluster=='ZNF683_NK')
IL7R_iNK  <- subset(sigposDEG.all, cluster=='IL7R_iNK')  

# 将不同细胞群体的上调基因保存为列表
list <- list(GNLY_mNK, EOMES_mNK, ZNF683_NK, IL7R_iNK)
names(list)[1:4] <- c("GNLY_mNK", "EOMES_mNK", "ZNF683_NK", "IL7R_iNK")
names(list)  

# GO分析
library(clusterProfiler)
for(i in 1:length(list)){GO <- enrichGO(gene =list[[i]]$gene,
                                        OrgDb = 'org.Bt.eg.db',
                                        keyType = 'SYMBOL',
                                        ont = "BP", 
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 1, 
                                        qvalueCutoff = 0.2)
GO<- data.frame(GO)
write.csv(GO,paste0("GO_new_", names(list[i]),".CSV"))
}  

# 读取GO富集分析结果
GNLY_mNK <- read.csv("GO_new_GNLY_mNK.CSV", row.names = 1) 
EOMES_mNK <- read.csv("GO_new_EOMES_mNK.CSV", row.names = 1)       
ZNF683_NK <- read.csv("GO_new_ZNF683_NK.CSV", row.names = 1)  
IL7R_iNK <- read.csv("GO_new_IL7R_iNK.CSV", row.names = 1)   

# 为每个细胞群体添加标签
GNLY_mNK$group <- "GNLY_mNK"
EOMES_mNK$group <- "EOMES_mNK"
ZNF683_NK$group <- "ZNF683_NK"
IL7R_iNK$group <- "IL7R_iNK"

# 选择TOP通路
# GNLY_mNK 
select_GNLY_mNK = c("cell killing",  
                 "leukocyte mediated cytotoxicity", 
                 "natural killer cell mediated immunity",
                 "natural killer cell mediated cytotoxicity")
# EOMES_mNK
select_EOMES_mNK = c("regulation of intracellular signal transduction",
                     "regulation of ERK1 and ERK2 cascade",
                    "immune response-activating signaling pathway",
                    "positive regulation of immune system process")
# ZNF683_NK 
select_ZNF683_NK = c("activation of immune response",
                     "T cell migration",
                     "T cell receptor signaling pathway",
                     "regulation of leukocyte migration")
# IL7R_iNK 
select_IL7R_iNK = c("positive regulation of mitotic cell cycle",
                     "positive regulation of mitotic cell cycle phase transition",
                     "protein maturation",
                     "protein folding")

# 选择每个亚群的通路
GNLY_mNK <- GNLY_mNK[GNLY_mNK$Description %in% select_GNLY_mNK,]
EOMES_mNK <- EOMES_mNK[EOMES_mNK$Description %in% select_EOMES_mNK,]
ZNF683_NK <- ZNF683_NK[ZNF683_NK$Description %in% select_ZNF683_NK,]
IL7R_iNK <- IL7R_iNK[IL7R_iNK$Description %in% select_IL7R_iNK,] 

# 生成新的P值列
GNLY_mNK$`-log10pvalue` <- -log10(GNLY_mNK$pvalue)
EOMES_mNK$`-log10pvalue` <- -log10(EOMES_mNK$pvalue)
ZNF683_NK$`-log10pvalue` <- -log10(ZNF683_NK$pvalue)
IL7R_iNK$`-log10pvalue` <- -log10(IL7R_iNK$pvalue) 

# 合并所有数据
all <- rbind(GNLY_mNK, EOMES_mNK, ZNF683_NK, IL7R_iNK)
all$Description <- gsub("-", " ", all$Description)  

# 为Description列设置因子顺序
library(forcats)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)  

# 分割GeneRatio列并计算生成比率
dat2 <- str_split(all$GeneRatio, "/", simplify = TRUE)[,1]
all$Num <- as.numeric(dat2)
dat3 <- str_split(all$GeneRatio, "/", simplify = TRUE)[,2]
all$Tot <- as.numeric(dat3)
all$Generatio <- as.numeric(all$Num / all$Tot) 

# 指定分组顺序
My_levels <- c("GNLY_mNK", "EOMES_mNK", "ZNF683_NK", "IL7R_iNK")
all$group <- factor(all$group, levels= My_levels)

# 绘制GO气泡图
ggplot(all, aes(group, Description)) +
  theme_bw() +
  geom_point(aes(fill = `-log10pvalue`, size = FoldEnrichment), shape = 21, colour = "black", alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_line(color = "#ececec", size = 0.5),
        panel.grid.minor = element_line(color = "#ececec", size = 0.5),
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 1)) +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_size_continuous(range = c(3, 6.5))


# =============================================
# 🎨 Fig. 7h, i
# =============================================

# 加载R包
library(slingshot)
library(Seurat)
library(Seurat)
library(cowplot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tradeSeq)
library(RColorBrewer)
library(scales)  

# 加载数据  
load("NK细胞鉴定.Rdata")

# 将Seurat对象转为SingleCellExperiment对象
sce <- as.SingleCellExperiment(NK, assay = "RNA") 

# 运行主函数
sce_slingshot1 <- slingshot(data = sce,                     
                            reducedDim = 'UMAP',            
                            clusterLabels = sce$celltype, 
                            start.clus = NULL,      
                            end.clus = NULL,           
                            approx_points = 150         
                            )

# 保存slingshot分析结果
save(sce_slingshot1, file = "NK拟时之Slingshot分析.Rdata")  
load("NK拟时之Slingshot分析.Rdata")

# SlingshotDataSet函数查看轨迹信息
SlingshotDataSet(sce_slingshot1)
dim(slingPseudotime(sce_slingshot1)) 
slingPseudotime(sce_slingshot1)[1:2,1:2]  

# 以细胞类型为基础进行可视化
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$celltype <- sce_slingshot1$celltype
umap_df$group <- sce_slingshot1$group
cell_colors <- c("GNLY_mNK" = "#ae7eb5",  
                 "EOMES_mNK" = "#f98177",  
                 "ZNF683_NK" = "#54c6a8",
                 "IL7R_iNK" = "#568389")
curves <- slingCurves(sce_slingshot1)
curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  crv <- as.data.frame(curves[[i]]$s)
  crv$lineage <- paste0("Lineage", i)
  colnames(crv)[1:2] <- c("UMAP_1", "UMAP_2")
  return(crv)
}))

# 画总的UMAP轨迹图
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 1, alpha = 0.8, show.legend = FALSE) + 
  scale_color_manual(values = cell_colors) +
  geom_path(data = curve_df, 
            aes(x = UMAP_1, y = UMAP_2, group = lineage), 
            inherit.aes = FALSE,
            color = "black", size = 1) +
  theme_classic2() + 
  labs(x = "UMAP_1", y = "UMAP_2", color = "Cell Type") +
  theme(legend.position = "right")

# 第1条轨迹
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$pseudotime <- sce_slingshot1$slingPseudotime_1
curve_obj <- slingCurves(sce_slingshot1)[[1]] 
curve_df <- as.data.frame(curve_obj$s)
colnames(curve_df) <- c("UMAP_1", "UMAP_2")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1.5, alpha = 0.9, show.legend = FALSE) +
  scale_color_viridis_c(option = "D", na.value = "grey80") +
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
            inherit.aes = FALSE, color = "black", size = 1) +
  theme_void() + 
  labs(x = "UMAP_1", y = "UMAP_2", color = "Pseudotime")
dev.off()

# 第2条轨迹
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$pseudotime <- sce_slingshot1$slingPseudotime_2
curve_obj <- slingCurves(sce_slingshot1)[[2]]  
curve_df <- as.data.frame(curve_obj$s)
colnames(curve_df) <- c("UMAP_1", "UMAP_2")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1.5, alpha = 0.9, show.legend = FALSE) +
  scale_color_viridis_c(option = "D", na.value = "grey80") +
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
            inherit.aes = FALSE, color = "black", size = 1) +
  theme_void() + 
  labs(x = "UMAP_1", y = "UMAP_2", color = "Pseudotime")

# 提取拟时间值加到元数据列即可
slingPseudotime = as.data.frame(slingPseudotime(sce_slingshot1))
identical(rownames(slingPseudotime),colnames(NK))
NK$Lineage1 = slingPseudotime$Lineage1
NK$Lineage2 = slingPseudotime$Lineage2 

# 轨迹1 
merge <- NK@meta.data[, c("Lineage1", "celltype")]
merge$cell <- rownames(merge)
merge <- merge[!is.na(merge$Lineage1), ]
table(merge$celltype)
merge <- merge[merge$celltype != "ZNF683_NK", ]
table(merge$celltype)
merge <- merge[order(merge$Lineage1), ]
merge$celltype <- factor(merge$celltype, levels = c("IL7R_iNK", "EOMES_mNK", "GNLY_mNK"))
ggplot(merge, aes(x=Lineage1,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#568389", 
                               "#f98177", 
                               "#ae7eb5")) 

# 轨迹2 
merge <- NK@meta.data[, c("Lineage2", "celltype")]
merge$cell <- rownames(merge)
merge <- merge[!is.na(merge$Lineage2), ]
table(merge$celltype)
merge <- merge[merge$celltype != "GNLY_mNK", ]
merge <- merge[merge$celltype != "EOMES_mNK", ]
table(merge$celltype)
merge <- merge[order(merge$Lineage2), ]
merge$celltype <- factor(merge$celltype, levels = c("IL7R_iNK", "ZNF683_NK")) 
ggplot(merge, aes(x=Lineage2,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#568389", "#54c6a8")) 
dev.off()


# =============================================
# 🎨 Fig. 7j
# =============================================

# 加载数据  
load("NK细胞鉴定.Rdata")

# 加载R包
library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr) 

# 构建CDS对象
expression_matrix <- as(as.matrix(NK@assays$RNA$counts), 'sparseMatrix') 
cell_metadata <- NK@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)    
  
# 构建Monocle3的专属CDS对象
cds <- new_cell_data_set(expression_matrix,             
                         cell_metadata = cell_metadata,  
                         gene_metadata = gene_annotation  
                         ) 
  
# CDS预处理和降维聚类
cds <- preprocess_cds(cds,
                      norm_method = "log", 
                      method = "PCA",      
                      num_dim = 50         
                      )   
  
# 绘制PCA方差解释图 显示前几个主成分解释的数据方差比例
plot_pc_variance_explained(cds)

# 降维聚类
cds <- reduce_dimension(cds,
                        reduction_method = 'UMAP',
                        preprocess_method = 'PCA'
                        )  

colnames(colData(cds)) 
plot_cells(cds, color_cells_by = "celltype", cell_size = 0.5, group_label_size = 5)
plot_cells(cds, color_cells_by = "orig.ident", cell_size = 0.5, group_label_size = 5)
Seurat::DimPlot(NK, reduction = "umap", pt.size = 0.2, label = T)   
cds <- cluster_cells(cds)
plot_cells(cds, cell_size = 0.5, group_label_size = 5) 
cds.embed <- cds@int_colData$reducedDims$UMAP           
int.embed <- Embeddings(NK, reduction = "umap")      
int.embed <- int.embed[rownames(cds.embed),]          
cds@int_colData$reducedDims$UMAP <- int.embed       
plot_cells(cds,
           color_cells_by = "celltype",
           cell_size = 0.5,
           group_label_size = 4)   

# 保存
mycds <- cds
save(mycds,file = "NK_mycds.Rdata")
rm(list = ls())
load("NK_mycds.Rdata") 

# 轨迹推断
mycds <- learn_graph(mycds,
                     verbose = T,
                     use_partition = F,         
                     close_loop = T,           
                     learn_graph_control = NULL 
                     )

# 可视化轨迹树
plot_cells(mycds3, 
           label_cell_groups = F,
           color_cells_by = 'celltype',
           label_groups_by_cluster = FALSE,
           label_leaves = F,
           label_branch_points = F,
           label_roots = F,
           graph_label_size = 2, 
           trajectory_graph_segment_size = 1.5,
           cell_size = 1,
           trajectory_graph_color = "#2c2828", 
           group_label_size = 4) +
  scale_color_manual(values = c("GNLY_mNK" = "#ae7eb5",  
                                "EOMES_mNK" = "#f98177",  
                                "ZNF683_NK" = "#54c6a8",
                                "IL7R_iNK" = "#568389")) +
  theme_void() +
  theme(legend.position = "none")

# 定义root cell
mycds3 <- mycds  
mycds3 <- order_cells(mycds3,
                      root_cells = colnames(mycds3[, mycds3@colData@listData[["celltype"]] == c("IL7R_iNK")]))
# 可视化拟时图
plot_cells(mycds3,
           label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           label_roots = F,
           graph_label_size = 5,
           cell_size = 1.5, 
           trajectory_graph_segment_size = 1.5,
           trajectory_graph_color = "#14ffb1"
           ) +
  scale_color_viridis(option = "A") + 
  theme_void() +
  theme(legend.position = "none")


# =============================================
# 🎨 Fig. 7k
# =============================================

load("NK细胞鉴定.Rdata")
load("NK_mycds3.RData")
mycds1 = mycds3
mycds1_L1 <- mycds1[, mycds1@colData@listData[["celltype"]] %in% c("IL7R_iNK", "EOMES_mNK", "GNLY_mNK") ] 

# 计算基因按照轨迹的显著性变化
mycds1_L1_res <- graph_test(mycds1_L1, 
                            neighbor_graph = "principal_graph",
                            cores = 1
                            )

# 保存
save(mycds1_L1_res,file = "./NK_mycds1_L1_res.Rdata")
load("NK_mycds1_L1_res.Rdata")

# 提取显著性变化基因的行
mycds2_res <- subset(mycds1_L1_res, q_value < 0.01)

# 按Moran’s I从高到低排序
mycds2_res <- mycds2_res[order(-mycds2_res$morans_I), ]

# 沿伪时间轨迹变化最显著的前1000个差异表达基因
genes = rownames(mycds2_res)[1:1000] 

# 提取差异基因的表达矩阵
library(ClusterGVis)
mat <- pre_pseudotime_matrix(cds_obj = mycds1_L1, 
                             gene_list = genes)
head(mat[1:5,1:5])
class(mat)

# kmeans 聚类
ck <- clusterData(obj = as.data.frame(mat),
                  cluster.method = "kmeans",
                  cluster.num = 4
                  )

# 保存
save(ck, file = "NK_L1_monocle3_ck.Rdata")

# 热图
visCluster(object = ck, 
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 11),
             col_color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11)
           ),
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = rownames(mat)[1:10]
           )


# =============================================
# 🎨 Fig. 7l
# =============================================

# GZMA
FeaturePlot(subset(NK, group == "Health"),
            features = "GZMA",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "blue"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) # 添加四周的框线
dev.off()

FeaturePlot(subset(NK, group == "Ketosis"),
            features = "GZMA",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "blue"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# PRF1
FeaturePlot(subset(NK, group == "Health"),
            features = "PRF1",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "blue"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

FeaturePlot(subset(NK, group == "Ketosis"),
            features = "PRF1",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "blue"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# NCR1
FeaturePlot(subset(NK, group == "Health"),
            features = "NCR1",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

FeaturePlot(subset(NK, group == "Ketosis"),
            features = "NCR1",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# NCR3
FeaturePlot(subset(NK, group == "Health"),
            features = "NCR3",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

FeaturePlot(subset(NK, group == "Ketosis"),
            features = "NCR3",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# FCGR3A
FeaturePlot(subset(NK, group == "Health"),
            features = "FCGR3A",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

FeaturePlot(subset(NK, group == "Ketosis"),
            features = "FCGR3A",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# KLRK1
FeaturePlot(subset(NK, group == "Health"),
            features = "KLRK1",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

FeaturePlot(subset(NK, group == "Ketosis"),
            features = "KLRK1",
            reduction = "umap",
            pt.size = 0.9,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

