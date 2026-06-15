
# F5/S9: Subclustering analysis of hepatic CD4⁺ and CD8⁺ T cells.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 CD4 T细胞亚群数据处理
# =============================================

# 对CD4细胞进行再次降维聚类分群
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
  
load("TT细胞鉴定.Rdata")
CD4T <- subset(TT, celltype %in% c("CD4 T"))
CD4T <- NormalizeData(CD4T, normalization.method = "LogNormalize", scale.factor = 10000) 
CD4T <- FindVariableFeatures(CD4T, selection.method = "vst", nfeatures = 2000) 

# 周期打分
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene) 
X = homologene(s.genes,inTax = 9606,outTax = 9913) 
Y = homologene(g2m.genes,inTax = 9606,outTax = 9913) 
s.genes = X$"9913"
g2m.genes = Y$"9913" 
CD4T <- CellCycleScoring(CD4T, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# 归一化缩放去除周期影响
CD4T <- ScaleData(CD4T, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(CD4T))

# 线性降维PCA 默认用高变基因集
CD4T <- RunPCA(CD4T, features = VariableFeatures(object = CD4T))

# 肘部图
ElbowPlot(CD4T, 50)
CD4T = FindNeighbors(CD4T, dims = 1:20) 
CD4T = FindClusters(CD4T, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(CD4T, prefix = "RNA_snn_res.") 
Idents(CD4T) <- "RNA_snn_res.0.5"
CD4T$seurat_clusters <- CD4T@active.ident
CD4T <- RunUMAP(CD4T, dims = 1:25)

# UMAP图
Seurat::DimPlot(CD4T, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(CD4T, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(CD4T, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(CD4T, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

# 保存数据
save(CD4T, file = "CD4T降维聚类.Rdata") 
rm(list = ls())
load("CD4T降维聚类.Rdata")

# 删掉未知身份的零碎小亚群
CD4T <- subset(CD4T, seurat_clusters %in% c("5"), invert = TRUE) 

# 需要重新走降维聚类流程
CD4T <- FindNeighbors(CD4T,
                      dims = 1:20 
                      )  

CD4T <- FindClusters(object = CD4T,
                   resolution = c(seq(0.1, 1, 0.1))
                   )
library(clustree)
clustree(CD4T, prefix = "RNA_snn_res.") 
Idents(CD4T) <- "RNA_snn_res.0.3"
CD4T$seurat_clusters <- CD4T@active.ident

head(Idents(CD4T), 10) 

CD4T <- RunUMAP(CD4T, dims = 1:20) 

# 保存数据
save(CD4T, file = "CD4T原始数据.Rdata") 
rm(list = ls())
load("CD4T原始数据.Rdata") 

CD4T = RunHarmony(CD4T, "orig.ident", plot_convergence = TRUE)   
CD4T = FindNeighbors(CD4T, reduction = "harmony", dims = 1:20)
CD4T <- RunUMAP(CD4T, reduction = "harmony", dims = 1:25) 
library(clustree)
clustree(CD4T, prefix = "RNA_snn_res.")
Idents(CD4T) <- "RNA_snn_res.0.2"
CD4T$seurat_clusters <- CD4T@active.ident

dif<-FindAllMarkers(CD4T, 
                    group.by = CD4T@meta.data$seurat_clusters, 
                    logfc.threshold = log2(1.2),                        
                    min.pct = 0.2,                                      
                    only.pos = T                                        
                    )        
dif$pct_diff <- dif$pct.1 - dif$pct.2 
table(dif$cluster)                    
dif<-dif %>%
  group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
save(dif, file = "CD4T整体cluster差异基因.Rdata")
  
# 鉴定细胞
load("CD4T整体cluster差异基因.Rdata")
    
# Marker  
# 炎症驱动 + 杀伤型 CD4 
FeaturePlot(CD4T, features = "XCL2", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "CST7", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "CXCR6", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "CCL4", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "CCL5", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "TNF", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "CTSW", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD4T, 
                features = c("XCL2", "CCL4", "CCL5", "CST7", "CTSW", "FCGR3A", "CD96", 
                             "CD69", "CD83", "TNFSF9",  
                             "IFNG", "TNF"
                             ), 
                group.by = "seurat_clusters") + RotatedAxis()

# naive-like
FeaturePlot(CD4T, features = "PHGDH", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "S1PR1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "KLF2", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "ITGB1", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD4T, 
                features = c("PHGDH", 
                             "S1PR1", "KLF2", 
                             "ITGA4", "SELPLG", "ITGB1", 
                             "TNFSF10", "LGALS3" 
                             ), 
                group.by = "seurat_clusters") + RotatedAxis() 

# Naive T
FeaturePlot(CD4T, features = "TCF7", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "LEF1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "CCR7", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "SELL", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "FOXP1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "BACH2", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "ID3", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD4T, 
                features = c("TCF7", "LEF1", "CCR7", "SELL", "FOXP1", "BACH2", "ID3"), 
                group.by = "seurat_clusters") + RotatedAxis()

# Treg 
FeaturePlot(CD4T, features = "FOXP3", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "IL2RA", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "CTLA4", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "LAG3", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "PDCD1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD4T, features = "TOX", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD4T, 
                features = c("FOXP3", "IL2RA", "CTLA4",
                             "LAG3", "PDCD1", "TOX"), 
                group.by = "seurat_clusters") + RotatedAxis()   
  
# 为分群重新指定细胞类型 
new.cluster.ids <- c("XCL2 CD4",
                     "S1PR1 CD4",
                     "Naive CD4",
                     "CD4 Treg") 
new.cluster.ids
names(new.cluster.ids) 
levels(CD4T)
names(new.cluster.ids) <- levels(CD4T) 
names(new.cluster.ids)
new.cluster.ids
CD4T <- RenameIdents(CD4T, new.cluster.ids)

CD4T[["celltype"]] <- Idents(CD4T) 
table(CD4T@meta.data[["orig.ident"]])
unique(CD4T@meta.data[["orig.ident"]])
CD4T[["group"]]<- c(rep("Health", 4695), rep("Ketosis", 3567))
CD4T[["group.celltype"]]<-paste(CD4T$group, Idents(CD4T), sep = '_') 
table(CD4T@meta.data[["orig.ident"]])
table(CD4T@meta.data[["group"]])
table(CD4T@meta.data[["celltype"]])
table(CD4T@meta.data[["group.celltype"]])
table(CD4T@meta.data[["seurat_clusters"]])

# 细胞水平信息
Idents(CD4T) <- factor(Idents(CD4T),
                      levels = c("Naive CD4", "S1PR1 CD4", "CD4 Treg", "XCL2 CD4"))
CD4T[["celltype"]] <- Idents(CD4T) 

# UMAP
Seurat::DimPlot(CD4T, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(CD4T, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(CD4T, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(CD4T, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 保存工作空间
save(CD4T,file = "CD4T细胞鉴定.Rdata") 
rm(list = ls())
load("CD4T细胞鉴定.Rdata")
  
# 对注释后细胞亚群 进行差异基因分析
diff <- FindAllMarkers(CD4T, 
                      group.by = CD4T@active.ident, 
                      logfc.threshold = log2(1.2), 
                      min.pct = 0.2, 
                      only.pos = T 
                      )
diff$pct_diff <- diff$pct.1 - diff$pct.2 
table(diff$cluster) 
head(diff[diff$cluster == unique(diff$cluster)[1],], 50) $ gene 
diff<-diff %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)
save(diff, file = "CD4T整体细胞类型差异基因.Rdata")    
  

# =============================================
# 🎨 Fig. 5a
# ============================================= 
  
load("CD4T细胞鉴定.Rdata")
Seurat::DimPlot(CD4T,
                group.by = "RNA_snn_res.0.3",
                cols = c(
                  "#de9f99",  
                  "#fcb4dc",  
                  "#f7bf92",  
                  "#99badf",
                  "#f29897",
                  "#96d0cb"), 
                pt.size = 0.1,
                label = T) +
  NoLegend()+ 
  labs(title = NULL)  
dev.off()  
  
  
# =============================================
# 🎨 Fig. 5b
# ============================================= 

# Marker基因气泡图
load("CD4T整体细胞类型差异基因.Rdata")
Seurat::DotPlot(CD4T, 
                features = c("SELL","LEF1", "TCF7","BACH2", "CCR7", 
                             "LGALS3", "ITGB1", "PHGDH", "S1PR1", "SELPLG",
                             "CTLA4","PDCD1", "LAG3", "FOXP3", "IL2RA", 
                             "CD69","TNF","XCL2","IFNG", "CTSW"),
                group.by = "celltype") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid = element_blank(), 
        axis.line = element_blank()) +   
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(11)))


# =============================================
# 🎨 Fig. 5c
# ============================================= 
  
plot_density(CD4T, 
             reduction = "umap",
             features = c("SELL"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 
  
plot_density(CD4T, 
             reduction = "umap",
             features = c("CCR7"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 
   
plot_density(CD4T, 
             reduction = "umap",
             features = c("S1PR1"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(CD4T, 
             reduction = "umap",
             features = c("PHGDH"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

plot_density(CD4T, 
             reduction = "umap",
             features = c("XCL2"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

plot_density(CD4T, 
             reduction = "umap",
             features = c("IFNG"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

plot_density(CD4T, 
             reduction = "umap",
             features = c("CST7"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(CD4T, 
             reduction = "umap",
             features = c("FOXP3"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

plot_density(CD4T, 
             reduction = "umap",
             features = c("CTLA4"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)


# =============================================
# 🎨 Fig. 5d
# ============================================= 

# 首先进行亚群间的差异分析
dif <- FindAllMarkers(CD4T, 
                      group.by = CD4T@active.ident, 
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

table(CD4T@active.ident)
Naive <- subset(sigposDEG.all, cluster=='Naive CD4')
S1PR1 <- subset(sigposDEG.all, cluster=='S1PR1 CD4')
Treg <- subset(sigposDEG.all, cluster=='CD4 Treg')
XCL2 <- subset(sigposDEG.all, cluster=='XCL2 CD4')  
  
# 将不同细胞群体的上调基因保存为列表
list <- list(Naive, S1PR1, Treg, XCL2)
names(list)[1:4] <- c("Naive", "S1PR1", "Treg", "XCL2")
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
Naive <- read.csv("GO_new_Naive.CSV", row.names = 1)
S1PR1 <- read.csv("GO_new_S1PR1.CSV", row.names = 1)
Treg <- read.csv("GO_new_Treg.CSV", row.names = 1)
XCL2 <- read.csv("GO_new_XCL2.CSV", row.names = 1) 
  
# 为每个细胞群体添加标签
Naive$group <- "Naive"
S1PR1$group <- "S1PR1"
Treg$group <- "Treg"
XCL2$group <- "XCL2"  
  
# 选择TOP通路
# Naive
select_Naive = c("autophagy",   # 幼稚CD8 T细胞具有较高的自噬通量 https://www.nature.com/articles/s41590-025-02090-1
               "cellular response to starvation",
               "cellular response to nutrient levels",
               "response to amino acid starvation")
# S1PR1
select_S1PR1 = c("cellular extravasation",
                 "Ras protein signal transduction",
                 "small GTPase-mediated signal transduction",
               "ribosome assembly")
# Treg
select_Treg = c("negative regulation of immune system process",
              "negative regulation of leukocyte activation",
              "negative regulation of cell activation",
              "regulation of T cell activation")
# XCL2
select_XCL2 = c("inflammatory response",
               "cellular response to type II interferon",
               "positive regulation of T cell activation",
               "positive regulation of leukocyte chemotaxis")  
  
# 选择每个亚群的通路
Naive <- Naive[Naive$Description %in% select_Naive,]
S1PR1 <- S1PR1[S1PR1$Description %in% select_S1PR1,]
Treg <- Treg[Treg$Description %in% select_Treg,]
XCL2 <- XCL2[XCL2$Description %in% select_XCL2,]  
  
# 生成新的P值列
Naive$`-log10pvalue` <- -log10(Naive$pvalue)
S1PR1$`-log10pvalue` <- -log10(S1PR1$pvalue)
Treg$`-log10pvalue` <- -log10(Treg$pvalue)
XCL2$`-log10pvalue` <- -log10(XCL2$pvalue)  
  
# 合并所有数据
all <- rbind(Naive, S1PR1, Treg, XCL2)
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
My_levels <- c("Naive", "S1PR1", "Treg", "XCL2")
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
  scale_fill_viridis(option = "A", direction = -1, limits = c(2, 6), oob = scales::squish) +
  scale_size_continuous(range = c(2, 6))
dev.off()


# =============================================
# 🎨 Fig. 5e
# ============================================= 

CellStatPlot(CD4T, group_by = "group", frac = "group", ident = "celltype", x_text_angle = 60) 


# =============================================
# 🎨 Fig. 5f, g
# =============================================

# 读取单细胞对象
load("CD4T细胞鉴定.Rdata")
# 修改细胞名为: 样本 + barcodes
colnames(CD4T)[1:5] 
meta.data <- CD4T@meta.data
meta.data$cellname <- rownames(meta.data)
rownames(meta.data) <- gsub('-1', '', rownames(meta.data))
rownames(meta.data) <- gsub('H1', 'Cont_3', rownames(meta.data))
rownames(meta.data) <- gsub('H2', 'FL_3', rownames(meta.data))
rownames(meta.data) <- gsub('H3', 'FL_4', rownames(meta.data))
rownames(meta.data) <- gsub('K1', 'Cont_1', rownames(meta.data))
rownames(meta.data) <- gsub('K2', 'FL_1', rownames(meta.data))
rownames(meta.data) <- gsub('K3', 'FL_5', rownames(meta.data))
rownames(meta.data) <- gsub('K4', 'FL_2', rownames(meta.data))
rownames(meta.data) <- gsub('K5', 'FL_6', rownames(meta.data))
rownames(meta.data)[1:5] 

# 将修改好的细胞名替换至Seurat对象中
CD4T <- Seurat::RenameCells(CD4T, new.names = rownames(meta.data))
colnames(CD4T)[1:5]

# 接下来将Seurat对象转为python中的数据格式 h5ad
# 加载R包
library(SCNT)
library(reticulate)
reticulate::use_python("/mnt/data/tool/miniconda3/envs/scanpy/bin/python", required = TRUE)
reticulate::py_config()
GetH5ad(CD4T,
        output_path = "CD4T.h5ad",
        mode = "sc",
        assay = "RNA")

# 接下来在JupyterLab中的scvelo小环境内核中操作
pwd
import scanpy as sc 
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scvelo as scv
sce=anndata.read_h5ad("CD4T.h5ad")
sce 
sce.obs.head() 
sce.obs['celltype'] = sce.obs['celltype'].astype('category')
sc.pl.umap(sce, color=['celltype'], frameon=False) 
ldat = anndata.read_loom('NCBI.loom')
ldat
ldat.obs.head() 
barcodes = [bc.replace('x',"") for bc in ldat.obs.index.tolist()] 
barcodes = [bc.replace(':',"_") for bc in barcodes]
ldat.obs.index = barcodes
ldat.var_names_make_unique() 
adata = scv.utils.merge(sce, ldat)
adata
adata.obs.head() 
adata
scv.pl.proportions(adata, groupby='celltype',layers=['spliced','unspliced','ambiguous'], save='RNA速率1.pdf')
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_grid(adata, 
                               basis='umap',
                               color='celltype',
                               save='RNA速率2.pdf', 
                               figsize=(6, 6), 
                               legend_loc="on data")
scv.tl.velocity_pseudotime(adata)
adata.uns['neighbors']['distances'] = adata.obsp['distances'] 
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities'] 
scv.tl.paga(adata, groups='celltype')  
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')  
scv.pl.paga(adata, 
            basis='umap',       
            size=50,             
            alpha=.3,           
            min_edge_width=2,   
            node_size_scale=1.5, 
            save='RNA速率3.pdf',
            figsize=(6, 6)
            ) 


# =============================================
# 🎨 Fig. 5h
# =============================================

load("CD4T细胞鉴定.Rdata")

library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)    

expression_matrix <- as(as.matrix(CD4T@assays$RNA$counts), 'sparseMatrix') 
cell_metadata <- CD4T@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix)) 
rownames(gene_annotation) <- rownames(expression_matrix) 

cds <- new_cell_data_set(expression_matrix,               
                         cell_metadata = cell_metadata,  
                         gene_metadata = gene_annotation  
                         )

cds <- preprocess_cds(cds,
                      norm_method = "log",
                      method = "PCA",     
                      num_dim = 50        
                      ) 

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,
                        reduction_method = 'UMAP',
                        preprocess_method = 'PCA'
                        )  

colnames(colData(cds))
plot_cells(cds, color_cells_by = "celltype", cell_size = 0.5, group_label_size = 5)
plot_cells(cds, color_cells_by = "orig.ident", cell_size = 0.5, group_label_size = 5)
Seurat::DimPlot(CD4T, reduction = "umap", pt.size = 0.2, label = T)  

cds <- cluster_cells(cds)
plot_cells(cds, cell_size = 0.5, group_label_size = 5)

cds.embed <- cds@int_colData$reducedDims$UMAP           
int.embed <- Embeddings(CD4T, reduction = "umap")      
int.embed <- int.embed[rownames(cds.embed),]            
cds@int_colData$reducedDims$UMAP <- int.embed    

plot_cells(cds,
           color_cells_by = "celltype",
           cell_size = 0.5,
           group_label_size = 4)   

mycds <- cds
save(mycds,file = "CD4T_mycds.Rdata")
rm(list = ls())
load("CD4T_mycds.Rdata") 

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
  scale_color_manual(values = c("Naive CD4" = "#f08e8d",  
                                "S1PR1 CD4" = "#89c4ed",  
                                "CD4 Treg" = "#8eccc7",
                                "XCL2 CD4" = "#a399cc")) +
  theme(legend.position = "none")

# 定义root cell
mycds3 <- mycds  
mycds3 <- order_cells(mycds3,
                      root_cells = colnames(mycds3[, mycds3@colData@listData[["celltype"]] == c("Naive CD4")]))
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
           trajectory_graph_color = "white"
           ) +
  scale_color_viridis(option = "D") + 
  theme(legend.position = "none")


# =============================================
# 🎨 Fig. 5i
# =============================================

# Monocle3拟时结果可视化
load("CD4T细胞鉴定.Rdata")
load("CD4T_mycds3.RData")

# 提取细胞的伪时间
pd <- pseudotime(mycds3, reduction_method = 'UMAP') 
head(pd)

# 将拟时结果添加到Seurat对象的元数据中 新增列命名为pseudotime
CD4T <- AddMetaData(CD4T, 
                    metadata = pd,
                    col.name = 'Monocle3_pseudotime'
                    ) 
colnames(CD4T@meta.data)

# 保存
save(CD4T, file = 'CD4T细胞鉴定.Rdata')

# 计算基因按照轨迹的显著性变化 
mycds3_res <- graph_test(mycds3, 
                         neighbor_graph = "principal_graph",
                         cores = 1
                         )

# 运行浪费时间 这里直接保存文件
save(mycds3_res,file = "CD4T_mycds3_res.Rdata")

# 提取显著性变化基因的行
mycds2_res <- subset(mycds3_res, q_value < 0.01)

# 按Moran’s I从高到低排序
mycds2_res <- mycds2_res[order(-mycds2_res$morans_I), ]

# 沿伪时间轨迹变化最显著的前1000个差异表达基因
genes = rownames(mycds2_res)[1:1000]

# 提取差异基因的表达矩阵
library(ClusterGVis)
mat <- pre_pseudotime_matrix(cds_obj = mycds3,
                             gene_list = genes)
head(mat[1:5,1:5])
class(mat)

# kmeans 聚类
ck <- clusterData(obj = as.data.frame(mat),
                  cluster.method = "kmeans",
                  cluster.num = 4 
                  )

# 保存一下这个ck 每次都不太一样
save(ck, file = "CD4T_monocle3_ck.Rdata")

# 热图 
visCluster(object = ck, 
           ht.col.list = list(col_range = seq(from = -2, to = 2, length.out = 25), 
                              col_color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(25))),
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = rownames(mat)[1:25]
           )

# 添加富集信息 KEGG
library(clusterProfiler)
library(org.Bt.eg.db)
enrichKEGG <- enrichCluster(object = ck,
                            type = "KEGG",
                            OrgDb = org.Bt.eg.db,
                            id.trans = TRUE,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            organism = "bta",
                            pvalueCutoff = 0.05,
                            topn = 20,
                            readable = TRUE,
                            seed = 5201314)

# 提取细胞的拟时间值
pseudotime <- pseudotime(mycds3) %>% as.data.frame() 
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "peu"

# 提取细胞注释信息
celltype <- data.frame(
  celltype = mycds3@colData@listData[["celltype"]],
  cell = rownames(mycds3@colData))

# 合并拟时间值和细胞注释信息 就可以得到细胞类型在拟时轴上的分布信息
merge <-merge(pseudotime, celltype, by = 'cell')

# 数据按照拟时序从小到大排列
merge <- merge[order(merge$peu), ]

# 分开展示
merge$celltype <- factor(merge$celltype, levels = levels(mycds3@colData@listData[["celltype"]]))
library(ggridges)
ggplot(merge, aes(x=peu,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#f08e8d", 
                               "#89c4ed", 
                               "#8eccc7", 
                               "#a399cc")) 


# =============================================
# 🎨 Supplementary Figure 9a
# =============================================
  
load("CD4T细胞鉴定.Rdata")
library(clustree)
clustree(CD4T, prefix = "RNA_snn_res.") 
Idents(CD4T) <- "RNA_snn_res.0.3"
CD4T$seurat_clusters <- CD4T@active.ident

# 画seurat_clusters分群的结果UMAP图
Seurat::DimPlot(CD4T, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)


# =============================================
# 🎨 Supplementary Figure 9b
# =============================================

library(Seurat)
library(pheatmap)
library(ggplot2)
library(dplyr) 
load("CD4T细胞鉴定.Rdata")
table(CD4T@meta.data[["celltype"]])
  
# 随机抽样
table(CD4T@meta.data[["celltype"]])
cell_types <- unique(CD4T$celltype)
subset_cells <- list()  
  
# 对每个细胞类型进行处理
for (cell_type in cell_types) {
  cells_of_type <- WhichCells(CD4T, expression = celltype == cell_type) 
  if (length(cells_of_type) >= 200) {
    selected_cells <- sample(cells_of_type, 200)                     
  } else {
    selected_cells <- cells_of_type                                 
  }
  subset_cells[[cell_type]] <- selected_cells                       
}  
  
# 合并所有子集细胞
subset_all_cells <- unlist(subset_cells)

# 根据选中的细胞创建新的 Seurat 对象
new_seurat_object <- subset(CD4T, cells = subset_all_cells) 
  
# 查看新的 Seurat 对象
new_seurat_object
  
# 基于上调基因分析挑选用于绘图的基因
dif<-FindAllMarkers(CD4T,
                    group.by = CD4T@active.ident, 
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


pheatmap(data1,scale = "none",cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,show_rownames = FALSE,
         annotation_col = celltype,    
         annotation_row = gene.anno,   
         annotation_names_row = FALSE,
         color = colorRampPalette(c("#040509","#608fe4", "#ffd700"))(100)
         )


# =============================================
# 🎨 CD8 T细胞亚群数据处理
# =============================================

load("~/奶牛肝脏解离单细胞3版/TT细胞鉴定.Rdata")
CD8T <- subset(TT, celltype %in% c("CD8 T"))
CD8T <- NormalizeData(CD8T, normalization.method = "LogNormalize", scale.factor = 10000) 
CD8T <- FindVariableFeatures(CD8T, selection.method = "vst", nfeatures = 2000) 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene) 
X = homologene(s.genes,inTax = 9606,outTax = 9913) 
Y = homologene(g2m.genes,inTax = 9606,outTax = 9913) 
s.genes = X$"9913"
g2m.genes = Y$"9913" 
CD8T <- CellCycleScoring(CD8T, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
CD8T <- ScaleData(CD8T, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(CD8T))
CD8T <- RunPCA(CD8T, features = VariableFeatures(object = CD8T))
ElbowPlot(CD8T, 50)
CD8T = FindNeighbors(CD8T, dims = 1:20)
CD8T = FindClusters(CD8T, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(CD8T, prefix = "RNA_snn_res.") 
Idents(CD8T) <- "RNA_snn_res.0.7"
CD8T$seurat_clusters <- CD8T@active.ident
CD8T <- RunUMAP(CD8T, dims = 1:25)  
  
# UMAP图
Seurat::DimPlot(CD8T, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(CD8T, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(CD8T, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(CD8T, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

# 保存数据
save(CD8T, file = "CD8T降维聚类.Rdata") 
rm(list = ls())
load("CD8T降维聚类.Rdata")  
CD8T = RunHarmony(CD8T, "orig.ident", plot_convergence = TRUE)   
CD8T = FindNeighbors(CD8T, reduction = "harmony", dims = 1:20)
CD8T <- RunUMAP(CD8T, reduction = "harmony", dims = 1:20) 
CD8T = FindClusters(CD8T, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(CD8T, prefix = "RNA_snn_res.") 
Idents(CD8T) <- "RNA_snn_res.0.5"
CD8T$seurat_clusters <- CD8T@active.ident

dif<-FindAllMarkers(CD8T, 
                    group.by = CD8T@meta.data$seurat_clusters, 
                    logfc.threshold = log2(1.2),                        
                    min.pct = 0.2,                                      
                    only.pos = T                                        
                    )        
dif$pct_diff <- dif$pct.1 - dif$pct.2 
table(dif$cluster)                    
dif<-dif %>%
  group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
head(dif[dif$cluster == unique(dif$cluster)[4],]$gene, 50)   
# 保存结果
save(dif, file = "CD8T整体cluster差异基因.Rdata")

# Marker  
# 0 CX3CR1
FeaturePlot(CD8T, features = "CX3CR1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "S1PR5", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "GNLY", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "ZEB2", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "KLF2", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD8T, 
                features = c("CX3CR1", "S1PR5", "GNLY", "ZEB2", "KLF2"), 
                group.by = "seurat_clusters") + RotatedAxis()

# 7 Exhausted CD8 T  
FeaturePlot(CD8T, features = "LAG3", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "PDCD1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "ENTPD1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "NR4A3", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD8T, 
                features = c("LAG3", "PDCD1", "ENTPD1", "NR4A3"), 
                group.by = "seurat_clusters") + RotatedAxis()

# 2和5 Tissue-resident memory CD8 T   
FeaturePlot(CD8T, features = "IKZF2", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "CD101", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "CD244", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "KLRD1", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "CD7", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "TRGC6", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "ZNF683", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "CD69", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "CD63", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD8T, 
                features = c("IKZF2", "CD101", "CD244", "KLRD1",
                             "CD7", "TRGC6", "ZNF683", "CD69", "CD63"), 
                group.by = "seurat_clusters") + RotatedAxis()

# 4/6   GZMK/EOMES/KLRG1   Effector/Terminal Effector CD8⁺T 
FeaturePlot(CD8T, features = "GZMK", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "EOMES", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "KLRG1", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD8T, 
                features = c("GZMK", "EOMES", "KLRG1"), 
                group.by = "seurat_clusters") + RotatedAxis()

# 1  GZMK/ZNF831   Effector/Memory-like CD8⁺T  
FeaturePlot(CD8T, features = "GZMK", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "ZNF831", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD8T, 
                features = c("GZMK", "ZNF831"), 
                group.by = "seurat_clusters") + RotatedAxis()

# 3  GZMK/IL7R  Early Activated/Memory CD8⁺T 
FeaturePlot(CD8T, features = "GZMK", reduction = "umap", pt.size = 0.1)
FeaturePlot(CD8T, features = "IL7R", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(CD8T, 
                features = c("GZMK", "IL7R"), 
                group.by = "seurat_clusters") + RotatedAxis()

# 为分群重新指定细胞类型 
new.cluster.ids <- c("CX3CR1_eff",
                     "GZMK_em",
                     "ZNF683_rm",
                     "GZMK_em",
                     "GZMK_em",
                     "ZNF683_rm",
                     "GZMK_em",
                     "LAG3_ex")
new.cluster.ids
names(new.cluster.ids) 
levels(CD8T)
names(new.cluster.ids) <- levels(CD8T) 
names(new.cluster.ids)
new.cluster.ids
CD8T <- RenameIdents(CD8T, new.cluster.ids) 

CD8T[["celltype"]] <- Idents(CD8T) 
table(CD8T@meta.data[["orig.ident"]])
unique(CD8T@meta.data[["orig.ident"]])
CD8T[["group"]]<- c(rep("Health", 7277), rep("Ketosis", 6367))
CD8T[["group.celltype"]]<-paste(CD8T$group, Idents(CD8T), sep = '_') 
table(CD8T@meta.data[["orig.ident"]])
table(CD8T@meta.data[["group"]])
table(CD8T@meta.data[["celltype"]])
table(CD8T@meta.data[["group.celltype"]])
table(CD8T@meta.data[["seurat_clusters"]])

# 细胞水平信息
Idents(CD8T) <- factor(Idents(CD8T),
                       levels = c("CX3CR1_eff", 
                                  "GZMK_em",
                                  "ZNF683_rm", 
                                  "LAG3_ex"))
CD8T[["celltype"]] <- Idents(CD8T)

Seurat::DimPlot(CD8T, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(CD8T, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(CD8T, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(CD8T, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 保存
save(CD8T,file = "CD8T细胞鉴定.Rdata")
rm(list = ls())
load("CD8T细胞鉴定.Rdata")

# 对注释后细胞亚群 进行差异基因分析
diff <- FindAllMarkers(CD8T, 
                      group.by = CD8T@active.ident, 
                      logfc.threshold = log2(1.2), 
                      min.pct = 0.2, 
                      only.pos = T 
                      )
diff$pct_diff <- diff$pct.1 - diff$pct.2 
table(diff$cluster) 
head(diff[diff$cluster == unique(diff$cluster)[4],], 50) $ gene
diff<-diff %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)
save(diff, file = "CD8T整体细胞类型差异基因.Rdata")    


# =============================================
# 🎨 Fig. 5j
# =============================================
  
load("~/奶牛肝脏解离单细胞3版/CD8T细胞鉴定.Rdata")
Seurat::DimPlot(CD8T,
                group.by = "seurat_clusters",
                cols = c(
                  "#c9a1b2",  
                  "#9bc0dc",  
                  "#f5b394",  
                  "#c6e2cd",
                  "#98d2d7",
                  "#e59589",
                  "#98a9cc",
                  "#b9b9b9"), 
                pt.size = 0.1,
                label = T) +
  NoLegend()+ 
  labs(title = NULL)  


# =============================================
# 🎨 Fig. 5k
# =============================================

# Marker基因的气泡图
load("CD8T整体细胞类型差异基因.Rdata")
library(scplotter)
library(plotthis)
plotthis::show_palettes(type = "continuous", index = 1:30)
Seurat::DotPlot(CD8T, 
                features = c("GNLY", "KLF2","CX3CR1", "ZEB2", "S1PR5",   # 循环效应型 外周巡游 高杀伤力 归巢相关
                             
                             "CD69", "GZMK", "CCR5", "IL7R", "GZMM",  # 活化早期记忆型 部分组织驻留潜力
                             
                             "IL2RB", "CD244", "CD7", "FCER1G", "ZNF683", # 组织驻留型 具有NK样特征 可能拥有增强杀伤力或快速应答能力 类似NK细胞 这是TRMCD8T的常见现象
                             
                             "LAG3", "GZMA", "NFATC1", "PDCD1", "ENTPD1" # 典型耗竭 免疫检查点 残余杀伤力
                             ),
                group.by = "celltype") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank()) + 
  scale_size(range = c(-1, 7.5)) +
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(11)))
dev.off()


# =============================================
# 🎨 Fig. 5l
# =============================================

plot_density(CD8T, 
             reduction = "umap",
             features = c("CX3CR1"),
             adjust = 1,  
             raster = T, 
             size = 1.2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)
dev.off()  

plot_density(CD8T, 
             reduction = "umap",
             features = c("S1PR5"),
             adjust = 1, 
             raster = T, 
             size = 1.2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

plot_density(CD8T, 
             reduction = "umap",
             features = c("GZMK"),
             adjust = 1,  
             raster = T, 
             size = 1.2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(CD8T, 
             reduction = "umap",
             features = c("ZNF683"),
             adjust = 1,  
             raster = T, 
             size = 1.2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(CD8T, 
             reduction = "umap",
             features = c("FCER1G"),
             adjust = 1,  
             raster = T, 
             size = 1.2) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

plot_density(CD8T, 
             reduction = "umap",
             features = c("LAG3"),
             adjust = 1,  
             raster = T, 
             size = 1.2) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 


# =============================================
# 🎨 Fig. 5m
# =============================================

# 首先进行亚群间的差异分析
dif <- FindAllMarkers(CD8T, 
                      group.by = CD8T@active.ident, 
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

table(CD8T@active.ident)
CX3CR1_eff <- subset(sigposDEG.all, cluster=='CX3CR1_eff') 
GZMK_em <- subset(sigposDEG.all, cluster=='GZMK_em')
ZNF683_rm <- subset(sigposDEG.all, cluster=='ZNF683_rm')
LAG3_ex <- subset(sigposDEG.all, cluster=='LAG3_ex')  

# 将不同细胞群体的上调基因保存为列表
list <- list(CX3CR1_eff, GZMK_em, ZNF683_rm, LAG3_ex)
names(list)[1:4] <- c("CX3CR1_eff", "GZMK_em", "ZNF683_rm", "LAG3_ex")
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
CX3CR1_eff <- read.csv("GO_new_CX3CR1_eff.CSV", row.names = 1) 
GZMK_em <- read.csv("GO_new_GZMK_em.CSV", row.names = 1)      
ZNF683_rm <- read.csv("GO_new_ZNF683_rm.CSV", row.names = 1)  
LAG3_ex <- read.csv("GO_new_LAG3_ex.CSV", row.names = 1)     

# 为每个细胞群体添加标签
CX3CR1_eff$group <- "CX3CR1_eff"
GZMK_em$group <- "GZMK_em"
ZNF683_rm$group <- "ZNF683_rm"
LAG3_ex$group <- "LAG3_ex"

# CX3CR1_eff
select_CX3CR1_eff = c("cell killing",  
                 "actin cytoskeleton organization", 
                 "small GTPase-mediated signal transduction",
                 "regulation of immune response")
# GZMK_em
select_GZMK_em = c("regulation of immune system process",
                 "T cell receptor signaling pathway",
                 "T cell differentiation", 
                 "leukocyte migration") 

select_ZNF683_rm = c("immune response-activating signaling pathway",
                "natural killer cell mediated cytotoxicity",
                "immune effector process",
                "ERK1 and ERK2 cascade")
# LAG3_ex 
select_LAG3_ex = c("positive regulation of MAPK cascade",
                "canonical NF-kappaB signal transduction",
                "response to tumor necrosis factor",
                "protein folding")  

CX3CR1_eff <- CX3CR1_eff[CX3CR1_eff$Description %in% select_CX3CR1_eff,]
GZMK_em <- GZMK_em[GZMK_em$Description %in% select_GZMK_em,]
ZNF683_rm <- ZNF683_rm[ZNF683_rm$Description %in% select_ZNF683_rm,]
LAG3_ex <- LAG3_ex[LAG3_ex$Description %in% select_LAG3_ex,]  

# 生成新的P值列
CX3CR1_eff$`-log10pvalue` <- -log10(CX3CR1_eff$pvalue)
GZMK_em$`-log10pvalue` <- -log10(GZMK_em$pvalue)
ZNF683_rm$`-log10pvalue` <- -log10(ZNF683_rm$pvalue)
LAG3_ex$`-log10pvalue` <- -log10(LAG3_ex$pvalue)  

# 合并所有数据
all <- rbind(CX3CR1_eff, GZMK_em, ZNF683_rm, LAG3_ex)
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
My_levels <- c("CX3CR1_eff", "GZMK_em", "ZNF683_rm", "LAG3_ex")
all$group <- factor(all$group, levels= My_levels)

# 绘制GO气泡图
all$FoldEnrichment <- pmin(all$FoldEnrichment, 12.5)
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
  scale_size_continuous(range = c(3, 6))


# =============================================
# 🎨 Fig. 5n
# =============================================

scplotter::CellStatPlot(CD8T, group_by = "group", frac = "group", ident = "celltype", x_text_angle = 60) 


# =============================================
# 🎨 Fig. 5o, p
# =============================================

# 读取单细胞对象
load("CD8T细胞鉴定.Rdata")

# 修改细胞名为: 样本 + barcodes
colnames(CD8T)[1:5] 
meta.data <- CD8T@meta.data
meta.data$cellname <- rownames(meta.data)
rownames(meta.data) <- gsub('-1', '', rownames(meta.data))
rownames(meta.data) <- gsub('H1', 'Cont_3', rownames(meta.data))
rownames(meta.data) <- gsub('H2', 'FL_3', rownames(meta.data))
rownames(meta.data) <- gsub('H3', 'FL_4', rownames(meta.data))
rownames(meta.data) <- gsub('K1', 'Cont_1', rownames(meta.data))
rownames(meta.data) <- gsub('K2', 'FL_1', rownames(meta.data))
rownames(meta.data) <- gsub('K3', 'FL_5', rownames(meta.data))
rownames(meta.data) <- gsub('K4', 'FL_2', rownames(meta.data))
rownames(meta.data) <- gsub('K5', 'FL_6', rownames(meta.data))
rownames(meta.data)[1:5] 

# 将修改好的细胞名替换至Seurat对象中
CD8T <- Seurat::RenameCells(CD8T, new.names = rownames(meta.data))
colnames(CD8T)[1:5]

# 接下来将Seurat对象转为python中的数据格式 h5ad
library(SCNT)
library(reticulate)
reticulate::use_python("/mnt/data/tool/miniconda3/envs/scanpy/bin/python", required = TRUE)
reticulate::py_config() 
GetH5ad(CD8T,
        output_path = "CD8T.h5ad",
        mode = "sc",
        assay = "RNA")

# 接下来在JupyterLab中的scvelo小环境内核中操作
pwd
import scanpy as sc 
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scvelo as scv
sce=anndata.read_h5ad("CD8T.h5ad")
sce 
sce.obs.head() 
sce.obs['celltype'] = sce.obs['celltype'].astype('category')
sc.pl.umap(sce, color=['celltype'], frameon=False) 
ldat = anndata.read_loom('NCBI.loom')
ldat
ldat.obs.head() 
barcodes = [bc.replace('x',"") for bc in ldat.obs.index.tolist()] 
barcodes = [bc.replace(':',"_") for bc in barcodes]
ldat.obs.index = barcodes
ldat.var_names_make_unique() 
adata = scv.utils.merge(sce, ldat)
adata
adata.obs.head() 
adata
scv.pl.proportions(adata, groupby='celltype',layers=['spliced','unspliced','ambiguous'], save='RNA速率4.pdf')
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_grid(adata, 
                               basis='umap',
                               color='celltype',
                               save='RNA速率5.pdf', 
                               figsize=(6, 6), 
                               legend_loc="on data")
scv.tl.velocity_pseudotime(adata)
adata.uns['neighbors']['distances'] = adata.obsp['distances'] 
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities'] 
scv.tl.paga(adata, groups='celltype')  
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')  
scv.pl.paga(adata, 
            basis='umap',       
            size=50,             
            alpha=.3,           
            min_edge_width=2,   
            node_size_scale=1.5, 
            save='RNA速率6.pdf',
            figsize=(6, 6)
            ) 


# =============================================
# 🎨 Fig. 5q
# =============================================

# 加载数据  
load("CD8T细胞鉴定.Rdata")

# 加载R包
library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)      
  
# 构建CDS对象
expression_matrix <- as(as.matrix(CD8T@assays$RNA$counts), 'sparseMatrix') 
cell_metadata <- CD8T@meta.data 
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
Seurat::DimPlot(CD8T, reduction = "umap", pt.size = 0.2, label = T)   
cds <- cluster_cells(cds)
plot_cells(cds, cell_size = 0.5, group_label_size = 5) 

cds.embed <- cds@int_colData$reducedDims$UMAP           
int.embed <- Embeddings(CD8T, reduction = "umap")      
int.embed <- int.embed[rownames(cds.embed),]           
cds@int_colData$reducedDims$UMAP <- int.embed      

plot_cells(cds,
           color_cells_by = "celltype",
           cell_size = 0.5,
           group_label_size = 4)   

# 保存
mycds <- cds
save(mycds,file = "./CD8T_mycds.Rdata")
rm(list = ls())
load("./CD8T_mycds.Rdata") 

# 轨迹推断
mycds <- learn_graph(mycds,
                     verbose = T,
                     use_partition = F,        
                     close_loop = T,            
                     learn_graph_control = NULL 
                     )

# 可视化轨迹树
plot_cells(mycds1,
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
  scale_color_manual(values = c("CX3CR1_eff" = "#c69daf",  
                                "GZMK_em" = "#8fb6b4",  
                                "ZNF683_rm" = "#f6cd96",
                                "LAG3_ex" = "#b3b3b3")) +
  theme(legend.position = "none")
dev.off()

# 运行order_cells 在交互界面中手动选择起点细胞
mycds1 <- mycds
mycds1 <- order_cells(mycds1)

# 可视化拟时图
plot_cells(mycds1,
           label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           label_roots = F,
           graph_label_size = 5,
           cell_size = 1.5, 
           trajectory_graph_segment_size = 1.5,
           trajectory_graph_color = "white"
           ) +
  scale_color_viridis(option = "D") + 
  theme(legend.position = "none")
dev.off()


# =============================================
# 🎨 Fig. 5r
# =============================================

# 加载数据 删掉数据中除了C6和CX3CR1_eff之外的亚群
load("CD8T细胞鉴定.Rdata")
load("CD8T_mycds1.RData")
mycds1_C6_CX3CR1 <- mycds1[, mycds1@colData@listData[["seurat_clusters"]] %in% c("0","6") ] 

# 计算基因按照轨迹的显著性变化 该函数返回数据框 包含了每个轨迹分支的显著性测试结果
mycds1_C6_CX3CR1_res <- graph_test(mycds1_C6_CX3CR1, 
                                  neighbor_graph = "principal_graph",
                                  cores = 1
                                  )

# 保存文件
save(mycds1_C6_CX3CR1_res,file = "CD8T_mycds1_C6_CX3CR1_res.Rdata")
load("CD8T_mycds1_C6_CX3CR1_res.Rdata")

# 提取显著性变化基因的行
mycds2_res <- subset(mycds1_C6_CX3CR1_res, q_value < 0.01)

# 按Moran’s I从高到低排序
mycds2_res <- mycds2_res[order(-mycds2_res$morans_I), ]
mycds2_res[rownames(mycds2_res) == "CX3CR1",]

# 沿伪时间轨迹变化最显著的前1000个差异表达基因
genes = rownames(mycds2_res)[1:1000] 

# 提取差异基因的表达矩阵
library(ClusterGVis)
mat <- pre_pseudotime_matrix(cds_obj = mycds1_C6_CX3CR1, 
                             gene_list = genes)
head(mat[1:5,1:5])
class(mat)

# kmeans 聚类
ck <- clusterData(obj = as.data.frame(mat),
                  cluster.method = "kmeans",
                  cluster.num = 4 
                  )

# 保存
save(ck, file = "CD8T_C6_CX3CR1_monocle3_ck.Rdata") 

# 热图 
visCluster(object = ck, 
           ht.col.list = list(col_range = seq(from = -2, to = 2, length.out = 25), 
                              col_color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(25))),
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = c(rownames(mat)[c(1:3, 6, 9:11, 13, 16:27, 29)], "RBPJ")
           )


# =============================================
# 🎨 Fig. 5s
# =============================================

# 加载数据 删掉数据中除了C6和ZNF683_rm之外的亚群
load("CD8T细胞鉴定.Rdata")
load("CD8T_mycds1.RData")
mycds1_C6_ZNF683 <- mycds1[, mycds1@colData@listData[["seurat_clusters"]] %in% c("2", "5", "6") ] 

# 计算基因按照轨迹的显著性变化 该函数返回数据框 包含了每个轨迹分支的显著性测试结果
mycds1_C6_ZNF683_res <- graph_test(mycds1_C6_ZNF683, 
                                  neighbor_graph = "principal_graph",
                                  cores = 1
                                  )

# 保存文件
save(mycds1_C6_ZNF683_res,file = "CD8T_mycds1_C6_ZNF683_res.Rdata")
load("CD8T_mycds1_C6_ZNF683_res.Rdata")

# 提取显著性变化基因的行
mycds2_res <- subset(mycds1_C6_ZNF683_res, q_value < 0.01)

# 按Moran’s I从高到低排序
mycds2_res <- mycds2_res[order(-mycds2_res$morans_I), ]
mycds2_res[rownames(mycds2_res) == "CX3CR1",]

# 沿伪时间轨迹变化最显著的前1000个差异表达基因
genes = rownames(mycds2_res)[1:1000]

# 提取差异基因的表达矩阵
library(ClusterGVis)
mat <- pre_pseudotime_matrix(cds_obj = mycds1_C6_ZNF683, 
                             gene_list = genes)
head(mat[1:5,1:5])
class(mat)

# kmeans 聚类
ck <- clusterData(obj = as.data.frame(mat),
                  cluster.method = "kmeans",
                  cluster.num = 4 
                  )

# 保存
save(ck, file = "CD8T_C6_ZNF683_monocle3_ck.Rdata") 

# 热图
visCluster(object = ck, 
           ht.col.list = list(col_range = seq(from = -2, to = 2, length.out = 25), 
                              col_color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(25))),
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = rownames(mat)[c(1, 3, 5:12, 14:18, 20:25)]
           )


# =============================================
# 🎨 Fig. 5t
# =============================================

# 交集基因
load("~/奶牛肝脏解离单细胞3版/CD8T_C6_ZNF683_monocle3_ck.Rdata") 
ZNF683_C4 = ck[["cluster.list"]][["C4"]] 

load("~/奶牛肝脏解离单细胞3版/CD8T_C6_CX3CR1_monocle3_ck.Rdata") 
CX3CR1_C1 = ck[["cluster.list"]][["C1"]] 

intersect(ZNF683_C4, CX3CR1_C1) 


# =============================================
# 🎨 Fig. 5u
# =============================================
  
load("CD8T细胞鉴定.Rdata")
CD8T_subset <- subset(CD8T, seurat_clusters %in% c("0", "2", "5", "6"))

Seurat::DotPlot(CD8T_subset, 
        features = c("RBPJ"),
        group.by = "celltype",
        scale = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_size_continuous(range = c(2, 6)) + 
  scale_color_distiller(palette = "YlOrRd", direction = 1)  +
  theme(panel.grid.major.x = element_line(color = "#ececec"),  
        panel.grid.major.y = element_line(color = "#ececec")) +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank()) + 
  coord_flip()

Seurat::DotPlot(CD8T_subset, 
                features = c("RBPJ"),
                group.by = "group",
                scale = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_size_continuous(range = c(2, 6)) +               
  scale_color_distiller(palette = "YlGnBu", direction = 1)  +
  theme(panel.grid.major.x = element_line(color = "#ececec"), 
        panel.grid.major.y = element_line(color = "#ececec")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank()) + 
  coord_flip()


# =============================================
# 🎨 Supplementary Figure 9c
# =============================================

# 加载对象
load("CD8T细胞鉴定.Rdata")
  
# 分群不变
library(clustree)
clustree(CD8T, prefix = "RNA_snn_res.") 
Idents(CD8T) <- "RNA_snn_res.0.3"
CD8T$seurat_clusters <- CD8T@active.ident

# 画seurat_clusters分群的结果UMAP图
Seurat::DimPlot(CD8T, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)


# =============================================
# 🎨 Supplementary Figure 9d
# =============================================

# 加载数据
load("CD8T细胞鉴定.Rdata")
table(CD8T@meta.data[["celltype"]]) 
  
# 随机抽样
table(CD8T@meta.data[["celltype"]])
cell_types <- unique(CD8T$celltype)
subset_cells <- list()  
  
# 对每个细胞类型进行处理
for (cell_type in cell_types) {
  cells_of_type <- WhichCells(CD8T, expression = celltype == cell_type) 
  if (length(cells_of_type) >= 200) {
    selected_cells <- sample(cells_of_type, 200)                   
  } else {
    selected_cells <- cells_of_type                                  
  }
  subset_cells[[cell_type]] <- selected_cells                       
}  
  
# 合并所有子集细胞
subset_all_cells <- unlist(subset_cells)

# 根据选中的细胞创建新的 Seurat 对象
new_seurat_object <- subset(CD8T, cells = subset_all_cells) 
  
# 查看新的 Seurat 对象
new_seurat_object
  
# 基于上调基因分析挑选用于绘图的基因
dif<-FindAllMarkers(CD8T,
                    group.by = CD8T@active.ident, 
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
# 🎨 Supplementary Figure 9e
# =============================================

# 加载数据
load("CD8T细胞鉴定.Rdata")  

# 画图
Seurat::DimPlot(CD8T,
                reduction = "umap", 
                group.by = "seurat_clusters",
                cols = c("#c9a1b2", 
                         "#eeeeee", 
                         "#f5d09b", 
                         "#eeeeee",
                         "#eeeeee",
                         "#f5d09b",
                         "#a4add0",
                         "#eeeeee"
                         ), 
                pt.size = 0.1,
                label = F) 


# =============================================
# 🎨 Supplementary Figure 9f
# =============================================

# 加载数据 删掉数据中除了C6和CX3CR1_eff之外的亚群
load("CD8T细胞鉴定.Rdata")
load("CD8T_mycds1.RData")
mycds1 <- mycds1[, mycds1@colData@listData[["seurat_clusters"]] %in% c("0","6") ] 

# 提取细胞的拟时间值
pseudotime <- pseudotime(mycds1) %>% as.data.frame() 
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "peu"

# 提取细胞注释信息
celltype <- data.frame(
  celltype = mycds1@colData@listData[["celltype"]],
  cell = rownames(mycds1@colData))

# 合并拟时间值和细胞注释信息 就可以得到细胞类型在拟时轴上的分布信息
merge <-merge(pseudotime, celltype, by = 'cell')

# 数据按照拟时序从小到大排列
merge <- merge[order(merge$peu), ]
merge$celltype <- factor(merge$celltype, levels = levels(mycds1@colData@listData[["celltype"]]))

# 作图
ggplot(merge, aes(x=peu,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#c69daf", 
                               "#a4add0")) 

# 加载数据 删掉数据中除了C6和CX3CR1_eff之外的亚群
load("CD8T细胞鉴定.Rdata")
load("CD8T_mycds1.RData")
mycds1 <- mycds1[, mycds1@colData@listData[["seurat_clusters"]] %in% c("2","5","6") ] 

# 提取细胞的拟时间值
pseudotime <- pseudotime(mycds1) %>% as.data.frame() 
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "peu"

# 提取细胞注释信息
celltype <- data.frame(
  celltype = mycds1@colData@listData[["celltype"]],
  cell = rownames(mycds1@colData))

# 合并拟时间值和细胞注释信息 就可以得到细胞类型在拟时轴上的分布信息
merge <-merge(pseudotime, celltype, by = 'cell')

# 数据按照拟时序从小到大排列
merge <- merge[order(merge$peu), ]
merge$celltype <- factor(merge$celltype, levels = levels(mycds1@colData@listData[["celltype"]]))

# 作图
ggplot(merge, aes(x=peu,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#a4add0", 
                               "#f6cd96")) 

