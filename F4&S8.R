
# F4/S8: Clinical ketosis is associated with functional dysregulation of hepatic T cell subsets.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 T细胞亚群数据处理
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
  

load("CD45去双细胞后.Rdata")
TT <- subset(seurat_object, celltype %in% c("T"))
TT <- NormalizeData(TT, normalization.method = "LogNormalize", scale.factor = 10000) 
TT <- FindVariableFeatures(TT, selection.method = "vst", nfeatures = 2000) 

# 周期打分
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene) 
X = homologene(s.genes,inTax = 9606,outTax = 9913) 
Y = homologene(g2m.genes,inTax = 9606,outTax = 9913) 
s.genes = X$"9913"
g2m.genes = Y$"9913" 
TT <- CellCycleScoring(TT, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) 

# 归一化缩放去除周期影响
TT <- ScaleData(TT, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(TT))

# 线性降维PCA 默认用高变基因集
TT <- RunPCA(TT, features = VariableFeatures(object = TT))

# 肘部图
ElbowPlot(TT, 50)
TT = FindNeighbors(TT, dims = 1:30)

# 分群数量
TT = FindClusters(TT, resolution = c(seq(0.1, 2, 0.1)))
library(clustree)
clustree(TT, prefix = "RNA_snn_res.")
Idents(TT) <- "RNA_snn_res.0.2"
TT$seurat_clusters <- TT@active.ident
TT <- RunUMAP(TT, dims = 1:25)

# UMAP图
Seurat::DimPlot(TT, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(TT, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(TT, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(TT, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

Seurat::DotPlot(TT, features = c("CD3D","CD3E",
                                 "CD4", 
                                 "CD8A", "CD8B", 
                                 "GNLY", "NKG7", "KLRB1", "KLRF1", "NCR1", 
                                 "LOC101908015", "SOX13",
                                 "LOC112442408", "SYNE1"), group.by = "seurat_clusters",scale = T) + RotatedAxis()

# 保存数据
save(TT, file = "TT降维聚类.Rdata") 
rm(list = ls())
load("TT降维聚类.Rdata")
  
# 删掉低质量零碎小群
TT <- subset(TT, seurat_clusters %in% c("2","7"), invert = TRUE)  

# 需要重新走降维聚类流程
TT <- FindNeighbors(TT,
                    dims = 1:30 
                    )  

# resolution参数决定下游聚类分析得到的分群数
TT <- FindClusters(object = TT,
                   resolution = c(seq(0.1, 1, 0.1))
                   )
library(clustree)
clustree(TT, prefix = "RNA_snn_res.") 
Idents(TT) <- "RNA_snn_res.0.1"
TT$seurat_clusters <- TT@active.ident

# UMAP非线性降维
TT <- RunUMAP(TT, dims = 1:25)

# UMAP图
Seurat::DimPlot(TT, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(TT, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(TT, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(TT, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)

# 保存数据
save(TT, file = "TT原始数据.Rdata") 
rm(list = ls())
load("TT原始数据.Rdata")   
  
# 整体差异分析
dif<-FindAllMarkers(TT, 
                    group.by = TT@meta.data$seurat_clusters, 
                    logfc.threshold = log2(1.2),                        
                    min.pct = 0.2,                                      
                    only.pos = T                                        
                    )        
dif$pct_diff <- dif$pct.1 - dif$pct.2 
table(dif$cluster)                    
dif<-dif %>%
  group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
save(dif, file = "TT整体cluster差异基因.Rdata")

# Marker
FeaturePlot(TT, features = "CD3D", reduction = "umap", pt.size = 0.1) 
FeaturePlot(TT, features = "CD3E", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(TT, features = c("CD3D","CD3E"), group.by = "seurat_clusters")

FeaturePlot(TT, features = "CD4", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(TT, features = c("CD4"), group.by = "seurat_clusters")

FeaturePlot(TT, features = "CD8A", reduction = "umap", pt.size = 0.1)
FeaturePlot(TT, features = "CD8B", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(TT, features = c("CD8A", "CD8B"), group.by = "seurat_clusters")

FeaturePlot(TT, features = "GNLY", reduction = "umap", pt.size = 0.1)
FeaturePlot(TT, features = "NKG7", reduction = "umap", pt.size = 0.1)
FeaturePlot(TT, features = "KLRB1", reduction = "umap", pt.size = 0.1)
FeaturePlot(TT, features = "KLRF1", reduction = "umap", pt.size = 0.1)
FeaturePlot(TT, features = "NCR1", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(TT, features = c("GNLY", "NKG7", "KLRB1", "KLRF1", "NCR1"), group.by = "seurat_clusters")

FeaturePlot(TT, features = "LOC101908015", reduction = "umap", pt.size = 0.1) # 又称TRDC
FeaturePlot(TT, features = "TRGC3", reduction = "umap", pt.size = 0.1)
FeaturePlot(TT, features = "TRGC4", reduction = "umap", pt.size = 0.1)
FeaturePlot(TT, features = "TRGC6", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(TT, features = c("LOC101908015", "TRGC3", "TRGC4", "TRGC6"), group.by = "seurat_clusters")

Seurat::DotPlot(TT, features = c("CD3D","CD3E",
                                  "CD4", 
                                  "CD8A", "CD8B", 
                                  "GNLY", "NKG7", "KLRB1", "ZBTB16", "KLRD1", "IFNG",
                                  "LOC101908015"), group.by = "seurat_clusters",scale = T) + RotatedAxis()

# 为分群重新指定细胞类型 
new.cluster.ids <- c("CD4 T",
                     "CD8 T",
                     "CD8 T",
                     "CD8 T",
                     "γδT",
                     "NKT") 
new.cluster.ids
names(new.cluster.ids) 
levels(TT)
names(new.cluster.ids) <- levels(TT)
names(new.cluster.ids)
new.cluster.ids
TT <- RenameIdents(TT, new.cluster.ids)

TT[["celltype"]] <- Idents(TT) 
table(TT@meta.data[["orig.ident"]])
unique(TT@meta.data[["orig.ident"]])
TT[["group"]]<- c(rep("Health", 15027), rep("Ketosis", 12858))
TT[["group.celltype"]]<-paste(TT$group, Idents(TT), sep = '_') 

# 细胞水平信息
Idents(TT) <- factor(Idents(TT),
                     levels = c("CD4 T", "CD8 T", "NKT", "γδT"))
TT[["celltype"]] <- Idents(TT) 

# UMAP
Seurat::DimPlot(TT, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(TT, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(TT, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(TT, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

Seurat::DotPlot(TT, features = c("CD3D","CD3E",
                                 "CD4", 
                                 "CD8A", "CD8B", 
                                 "GNLY", "NKG7", "KLRD1", "IFNG", "KLRB1", "ZBTB16", 
                                 "LOC101908015","SOX13"), group.by = "celltype",scale = T) + RotatedAxis()

# 保存
save(TT,file = "TT细胞鉴定.Rdata")
rm(list = ls())
load("TT细胞鉴定.Rdata")

# 对注释后细胞亚群 进行差异基因分析
diff <- FindAllMarkers(TT, 
                      group.by = TT@active.ident, 
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
save(diff, file = "TT整体细胞类型差异基因.Rdata")  


# =============================================
# 🎨 Fig. 4a
# =============================================

load("TT细胞鉴定.Rdata")
Seurat::DimPlot(TT,
                group.by = "celltype",
                cols = c("#8dd3c7", 
                         "#80b1d3",
                         "#fa9fb5", 
                         "#beaed4"), 
                pt.size = 0.1,
                label = T) +
  NoLegend()+ 
  labs(title = NULL)  
dev.off()


# =============================================
# 🎨 Fig. 4b
# =============================================

load("TT整体细胞类型差异基因.Rdata")
library(scplotter)
library(plotthis)
plotthis::show_palettes(type = "continuous", index = 1:30)
Seurat::DotPlot(TT, 
                features = c("CD3D","CD3E",
                             "CD4", "LEF1", "CD5",
                             "CD8A", "CD8B", 
                             "GNLY", "NKG7", "KLRB1", "ZBTB16", 
                             "LOC101908015","RHEX","SOX13"),
                group.by = "celltype") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_distiller(palette = "PuRd",direction = 1)  +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid = element_blank(), 
        axis.line = element_blank()) +  
  scale_size(range = c(-0.1, 8)) 
dev.off()


# =============================================
# 🎨 Fig. 4c
# =============================================

# 经典Marker基因的UMAP图
load("TT细胞鉴定.Rdata")
library(Seurat)
library(Nebulosa)
library(ggnetwork)
library(dplyr)

# CD3E
plot_density(TT, 
             reduction = "umap",
             features = c("CD3E"),
             pal = 'magma', 
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

# CD4
plot_density(TT, 
             reduction = "umap",
             features = c("CD4"),
             pal = 'magma', 
             adjust = 1,  # 调整带宽
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

# CD8B
plot_density(TT, 
             reduction = "umap",
             features = c("CD8B"),
             pal = 'magma', 
             adjust = 1,  # 调整带宽
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# GNLY
plot_density(TT, 
             reduction = "umap",
             features = c("GNLY"),
             pal = 'magma', 
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# KLRB1
plot_density(TT, 
             reduction = "umap",
             features = c("KLRB1"),
             pal = 'magma', 
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

# NKG7
plot_density(TT, 
             reduction = "umap",
             features = c("NKG7"),
             pal = 'magma', 
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

# LOC101908015
plot_density(TT, 
             reduction = "umap",
             features = c("LOC101908015"),
             pal = 'magma', 
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma", 
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)

# SOX13
plot_density(TT, 
             reduction = "umap",
             features = c("SOX13"),
             pal = 'magma', 
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "magma",
                        oob = scales::squish) +
  theme_blank() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1)


# =============================================
# 🎨 Fig. 4d
# =============================================

# 加载数据
load("TT细胞鉴定.Rdata")

# 随机抽样
table(TT@meta.data[["celltype"]])
cell_types <- unique(TT$celltype)
subset_cells <- list()

# 对每个细胞类型进行处理
for (cell_type in cell_types) {
  cells_of_type <- WhichCells(TT, expression = celltype == cell_type) 
  if (length(cells_of_type) >= 300) {
    selected_cells <- sample(cells_of_type, 300)                                 
  } else {
    selected_cells <- cells_of_type                                              
  }
  subset_cells[[cell_type]] <- selected_cells                                   
}

# 合并所有子集细胞
subset_all_cells <- unlist(subset_cells)

# 根据选中的细胞创建新的 Seurat 对象
new_seurat_object <- subset(TT, cells = subset_all_cells)

# 查看新的Seurat对象
new_seurat_object

# 基于上调基因分析挑选用于绘图的基因
dif<-FindAllMarkers(TT,
                    group.by = TT@active.ident, 
                    logfc.threshold = log2(1.2),            
                    min.pct = 0.2,                        
                    only.pos = T)
sig.dif<-dif%>%
  group_by(cluster)%>%
  top_n(n=10,wt=avg_log2FC) 
genes<-unique(sig.dif$gene)

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
# 🎨 Fig. 4e
# =============================================

dif <- FindAllMarkers(TT, 
                      group.by = TT@active.ident, 
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

# 将不同细胞群体的上调基因保存为列表
list <- list(CD4, CD8, NKT, y6T)
names(list)[1:4] <- c("CD4", "CD8", "NKT", "y6T")
names(list)

# GO分析
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
CD4 <- read.csv("GO_new_CD4.CSV", row.names = 1)
CD8 <- read.csv("GO_new_CD8.CSV", row.names = 1)
NKT <- read.csv("GO_new_NKT.CSV", row.names = 1)
y6T <- read.csv("GO_new_y6T.CSV", row.names = 1)

# 为每个细胞群体添加标签
CD4$group <- "CD4"
CD8$group <- "CD8"
NKT$group <- "NKT"
y6T$group <- "y6T"

# 选择TOP通路
# CD4
select_CD4 = c("leukocyte cell-cell adhesion",
               "regulation of T cell activation",
               "T cell costimulation",
               "lymphocyte costimulation")
# CD8
select_CD8 = c("cell killing",
                 "leukocyte mediated cytotoxicity",
                 "immune effector process",
               "immune response-activating cell surface receptor signaling pathway")
# NKT
select_NKT = c("T cell receptor signaling pathway",
              "antigen receptor-mediated signaling pathway",
              "natural killer cell mediated cytotoxicity",
              "natural killer cell mediated immunity")
# y6T
select_y6T = c("calcium-mediated signaling",
               "interferon-mediated signaling pathway",
               "cellular extravasation",
               "intracellular calcium ion homeostasis")

# 选择每个亚群的通路
CD4 <- CD4[CD4$Description %in% select_CD4,]
CD8 <- CD8[CD8$Description %in% select_CD8,]
NKT <- NKT[NKT$Description %in% select_NKT,]
y6T <- y6T[y6T$Description %in% select_y6T,]

# 生成新的P值列
CD4$`-log10pvalue` <- -log10(CD4$pvalue)
CD8$`-log10pvalue` <- -log10(CD8$pvalue)
NKT$`-log10pvalue` <- -log10(NKT$pvalue)
y6T$`-log10pvalue` <- -log10(y6T$pvalue)

# 合并所有数据
all <- rbind(CD4, CD8, NKT, y6T)
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
My_levels <- c("CD4", "CD8", "NKT", "y6T")
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
  scale_fill_viridis(option = "C", direction = -1, limits = c(0, 8), oob = scales::squish) +
  scale_size_continuous(range = c(2, 6))
dev.off()


# =============================================
# 🎨 Fig. 4f
# =============================================

load("~/奶牛肝脏解离单细胞3版/TT细胞鉴定.Rdata")  
load("~/奶牛肝脏解离单细胞3版/CD45去双细胞后.Rdata")
table(TT@active.ident, TT$orig.ident)
table(seurat_object$orig.ident)

# 手动计算
# H1   H2   H3   K1   K2   K3   K4   K5
# CD4  1509/6478  1587/10360  1771/11676   833/10960   828/8113   538/9263   1133/9334   349/7538 
# CD8  928/6478  3107/10360  3242/11676  818/10960  1922/8113   980/9263  2090/9334   557/7538
# NKT  278/6478  511/10360  562/11676  137/10960  330/8113  186/9263  436/9334  205/7538
# γδT  318/6478  609/10360  605/11676  570/10960  176/8113  336/9263  215/9334  219/7538 
  
# CD4  0.2329  0.1532  0.1517   0.076   0.1021   0.05808   0.1214   0.0463 
# CD8  0.1433  0.2999  0.2777  0.07464  0.2369   0.1058  0.2239   0.07389
# NKT  0.04291  0.04932  0.04813  0.0125  0.04068  0.02008  0.04671  0.0272
# γδT  0.04909  0.05878  0.0518  0.05201  0.02169  0.03627  0.02303  0.02905  
  
# 整理数据
CellRa <- data.frame(
  CD4 = c(0.2329, 0.1532, 0.1517, 0.076, 0.1021, 0.05808, 0.1214, 0.0463),
  CD8 = c(0.1433, 0.2999, 0.2777, 0.07464, 0.2369, 0.1058, 0.2239, 0.07389),
  NKT = c(0.04291, 0.04932, 0.04813, 0.0125, 0.04068, 0.02008, 0.04671, 0.0272),
  y6T = c(0.04909, 0.05878, 0.0518, 0.05201, 0.02169, 0.03627, 0.02303, 0.02905)
  )
CellRa$group <- c("H", "H", "H", "K", "K", "K", "K", "K")

library(tidyr)
CellRa_long <- pivot_longer(CellRa,
                            cols = c(CD4, CD8, NKT, y6T),
                            names_to = "Var1",
                            values_to = "Freq")

# CD4频率分布图
set.seed(123)
ggplot(data = CellRa_long[CellRa_long$Var1=="CD4",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 29) +  # 设置Y轴范围
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 13, color = "black"),  
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
  
shapiro.test(CellRa_long[CellRa_long$Var1=="CD4",]$Freq[1:3])
shapiro.test(CellRa_long[CellRa_long$Var1=="CD4",]$Freq[4:8])
leveneTest(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="CD4",], center = "mean") 
wilcox.test(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="CD4",]) # 0.03571 独立两组-非正态
  
# CD8频率分布图
set.seed(123)
ggplot(data = CellRa_long[CellRa_long$Var1=="CD8",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 38) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 13, color = "black"),  
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

shapiro.test(CellRa_long[CellRa_long$Var1=="CD8",]$Freq[1:3]) 
shapiro.test(CellRa_long[CellRa_long$Var1=="CD8",]$Freq[4:8]) 
leveneTest(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="CD8",], center = "mean")
t.test(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="CD8",], var.equal = TRUE)$p.value # 0.1563 独立两组-正态-方差齐性
  
# NKT频率分布图
dev.off()
pdf("P148.pdf", width = 1.65, height = 2)
set.seed(123)
ggplot(data = CellRa_long[CellRa_long$Var1=="NKT",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 6.1) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 13, color = "black"),  
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
dev.off()  
shapiro.test(CellRa_long[CellRa_long$Var1=="NKT",]$Freq[1:3])
shapiro.test(CellRa_long[CellRa_long$Var1=="NKT",]$Freq[4:8])
leveneTest(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="NKT",], center = "mean")
t.test(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="NKT",], var.equal = TRUE)$p.value # 0.08948 独立两组-正态-方差齐性

# γδT频率分布图
set.seed(123)
ggplot(data = CellRa_long[CellRa_long$Var1=="y6T",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 7.2) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 13, color = "black"),  
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
dev.off()  
shapiro.test(CellRa_long[CellRa_long$Var1=="y6T",]$Freq[1:3])
shapiro.test(CellRa_long[CellRa_long$Var1=="y6T",]$Freq[4:8])
leveneTest(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="y6T",], center = "mean")
t.test(Freq ~ group, data = CellRa_long[CellRa_long$Var1=="y6T",], var.equal = TRUE)$p.value # 0.03507 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 4g
# =============================================

load("TT细胞鉴定.Rdata")  
library(dplyr)
library(ggplot2)
library(gtools)
library(ggalluvial)

# 准备细胞比例输入数据
prop_df <- TT@meta.data %>%
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
  scale_fill_manual(values = c("#90d5c9",
                               "#86b3d5",
                               "#faa3b8",
                               "#c2b1d6")) +  
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
# 🎨 Fig. 4j
# =============================================

# 对所有的组间细胞亚群进行差异基因分析
load("TT细胞鉴定.Rdata")

# CD4
group.dif.CD4 <- FindMarkers(TT, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_CD4 T",   
                            ident.2 = "Health_CD4 T",   
                            logfc.threshold = 0,         
                            min.pct = 0.1,              
                            only.pos = FALSE           
                            ) 
group.dif.CD4$pct_diff <- group.dif.CD4$pct.1 - group.dif.CD4$pct.2 
group.dif.CD4<-group.dif.CD4 %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.CD4, file = 'group.dif.CD4.Rdata')
write.table(group.dif.CD4,"group.dif.CD4.xls",row.names = T,col.names = NA,quote = F,sep = "\t") 
table(group.dif.CD4$p_val < 0.05 & group.dif.CD4$avg_log2FC > log2(1.2))
table(group.dif.CD4$p_val < 0.05 & abs(group.dif.CD4$avg_log2FC) > log2(1.2))
table(group.dif.CD4$p_val < 0.05 & group.dif.CD4$avg_log2FC > 0)
table(group.dif.CD4$p_val < 0.05 & abs(group.dif.CD4$avg_log2FC) > 0) # 2153差异基因 541上调 1612下调

# CD8
group.dif.CD8 <- FindMarkers(TT, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_CD8 T",    
                            ident.2 = "Health_CD8 T",    
                            logfc.threshold = 0,         
                            min.pct = 0.1,             
                            only.pos = FALSE        
                            ) 
group.dif.CD8$pct_diff <- group.dif.CD8$pct.1 - group.dif.CD8$pct.2 
group.dif.CD8<-group.dif.CD8 %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.CD8, file = 'group.dif.CD8.Rdata')
write.table(group.dif.CD8,"group.dif.CD8.xls",row.names = T,col.names = NA,quote = F,sep = "\t") 
table(group.dif.CD8$p_val < 0.05 & group.dif.CD8$avg_log2FC > log2(1.2))
table(group.dif.CD8$p_val < 0.05 & abs(group.dif.CD8$avg_log2FC) > log2(1.2))
table(group.dif.CD8$p_val < 0.05 & group.dif.CD8$avg_log2FC > 0)
table(group.dif.CD8$p_val < 0.05 & abs(group.dif.CD8$avg_log2FC) > 0) # 2271差异基因 1245上调 1026下调

# NKT
group.dif.NKT <- FindMarkers(TT, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_NKT",  
                            ident.2 = "Health_NKT",    
                            logfc.threshold = 0,       
                            min.pct = 0.1,          
                            only.pos = FALSE       
                            ) 
group.dif.NKT$pct_diff <- group.dif.NKT$pct.1 - group.dif.NKT$pct.2 
group.dif.NKT<-group.dif.NKT %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.NKT, file = 'group.dif.NKT.Rdata')
write.table(group.dif.NKT,"group.dif.NKT.xls",row.names = T,col.names = NA,quote = F,sep = "\t") 
table(group.dif.NKT$p_val < 0.05 & group.dif.NKT$avg_log2FC > log2(1.2))
table(group.dif.NKT$p_val < 0.05 & abs(group.dif.NKT$avg_log2FC) > log2(1.2))
table(group.dif.NKT$p_val < 0.05 & group.dif.NKT$avg_log2FC > 0)
table(group.dif.NKT$p_val < 0.05 & abs(group.dif.NKT$avg_log2FC) > 0) # 2774差异基因 374上调 2400下调

# γδT
group.dif.γδT <- FindMarkers(TT, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_γδT",   
                            ident.2 = "Health_γδT",    
                            logfc.threshold = 0,      
                            min.pct = 0.1,        
                            only.pos = FALSE   
                            ) 
group.dif.γδT$pct_diff <- group.dif.γδT$pct.1 - group.dif.γδT$pct.2 
group.dif.γδT<-group.dif.γδT %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.γδT, file = 'group.dif.γδT.Rdata')
write.table(group.dif.γδT,"group.dif.γδT.xls",row.names = T,col.names = NA,quote = F,sep = "\t") 
table(group.dif.γδT$p_val < 0.05 & group.dif.γδT$avg_log2FC > log2(1.2))
table(group.dif.γδT$p_val < 0.05 & abs(group.dif.γδT$avg_log2FC) > log2(1.2))
table(group.dif.γδT$p_val < 0.05 & group.dif.γδT$avg_log2FC > 0)
table(group.dif.γδT$p_val < 0.05 & abs(group.dif.γδT$avg_log2FC) > 0) # 2462差异基因 312上调 2150下调

# 逐个修改差异基因列表的格式
# CD4
head(group.dif.CD4)
group.dif.CD4$gene <- rownames(group.dif.CD4)
group.dif.CD4$Group <- "CD4"
group.dif.CD4$Regulated <- ifelse(group.dif.CD4$p_val < 0.05 & group.dif.CD4$avg_log2FC > 0, "Up",
                                ifelse(group.dif.CD4$p_val < 0.05 & group.dif.CD4$avg_log2FC < 0, "Down", "Stable"))
group.dif.CD4$log10fdr <- -log10(group.dif.CD4$p_val)
head(group.dif.CD4)
# CD8
head(group.dif.CD8)
group.dif.CD8$gene <- rownames(group.dif.CD8)
group.dif.CD8$Group <- "CD8"
group.dif.CD8$Regulated <- ifelse(group.dif.CD8$p_val < 0.05 & group.dif.CD8$avg_log2FC > 0, "Up",
                                 ifelse(group.dif.CD8$p_val < 0.05 & group.dif.CD8$avg_log2FC < 0, "Down", "Stable"))
group.dif.CD8$log10fdr <- -log10(group.dif.CD8$p_val)
head(group.dif.CD8)
# NKT
head(group.dif.NKT)
group.dif.NKT$gene <- rownames(group.dif.NKT)
group.dif.NKT$Group <- "NKT"
group.dif.NKT$Regulated <- ifelse(group.dif.NKT$p_val < 0.05 & group.dif.NKT$avg_log2FC > 0, "Up",
                                ifelse(group.dif.NKT$p_val < 0.05 & group.dif.NKT$avg_log2FC < 0, "Down", "Stable"))
group.dif.NKT$log10fdr <- -log10(group.dif.NKT$p_val)
head(group.dif.NKT)
# PC
head(group.dif.γδT)
group.dif.γδT$gene <- rownames(group.dif.γδT)
group.dif.γδT$Group <- "γδT"
group.dif.γδT$Regulated <- ifelse(group.dif.γδT$p_val < 0.05 & group.dif.γδT$avg_log2FC > 0, "Up",
                                     ifelse(group.dif.γδT$p_val < 0.05 & group.dif.γδT$avg_log2FC < 0, "Down", "Stable"))
group.dif.γδT$log10fdr <- -log10(group.dif.γδT$p_val)
head(group.dif.γδT)

# 上下合并
diff_sc <- rbind(group.dif.CD4,
                 group.dif.CD8,
                 group.dif.NKT,
                 group.dif.γδT)

# 设置水平信息
levels(diff_sc$Group)
diff_sc$Group <- factor(diff_sc$Group, levels = c("CD4", "CD8", "NKT", "γδT"))
levels(diff_sc$Group)

# 展示上调TOP基因
up_genes <- data.frame(Group = c(rep("CD4", 5), rep("CD8", 5), rep("NKT", 5), rep("γδT", 5)),
                       gene = c("CXCR4", "CD69", "IL2", "ISG15", "FASLG",  # CD4
                                "XCL1", "XCL2", "CXCR4", "ZNF683", "IL2RB", # CD8
                                "IFNG", "TNF", "FASLG", "CCL4", "CD69", # NKT
                                "CCL4", "ISG15", "IL6R", "CXCL12", "GPR35" # γδT
                                )
                       )

# 展示下调TOP基因 
down_genes <- data.frame(Group = c(rep("CD4", 5), rep("CD8", 5), rep("NKT", 5), rep("γδT", 5)),
                         gene = c("TCF7", "LEF1", "FOXP1", "SELL", "BACH2",  # CD4
                                "CX3CR1", "GNLY", "GZMB", "KLRG1", "ZEB2", # CD8
                                "KLRD1", "KLRG1", "CPT1A", "IL6ST", "IKZF3", # NKT
                                "PRF1", "KLRD1", "ZBTB16", "AHR", "PRDM1" # γδT
                                )
                         )

# 从diff_sc中提取上下调挑选基因 
selected_up <- diff_sc %>%
  inner_join(up_genes, by = c("Group", "gene")) %>%
  arrange(Group, desc(avg_log2FC)) 
selected_down <- diff_sc %>%
  inner_join(down_genes, by = c("Group", "gene")) %>%
  arrange(Group, avg_log2FC)     

# 查看前几行确认是否正确
head(selected_up, 10)
head(selected_down, 10)

# 细胞类型的颜色
mycol <- c("#90d5c9",
           "#86b3d5",
           "#faa3b8",
           "#c2b1d6")

# 明确细胞顺序
cell_order <- c("CD4", "CD8", "NKT", "γδT")
diff_sc$Group <- factor(diff_sc$Group, levels = cell_order)
selected_up$Group <- factor(selected_up$Group, levels = cell_order)
selected_down$Group <- factor(selected_down$Group, levels = cell_order)

# 基础火山图绘制
p <- ggplot() +
  geom_point(data = diff_sc %>% dplyr::filter(Regulated == "Stable"),
             aes(x = avg_log2FC, y = log10fdr),
             size = 0.8,
             color = "grey20",
             alpha = 0.8) +
  geom_point(data = diff_sc %>% dplyr::filter(Regulated != "Stable"),
             aes(x = avg_log2FC, y = log10fdr),
             size = 0.8,
             color = "grey",
             alpha = 0.6) +
  coord_flip() +
  facet_grid(. ~ Group, scales = "free") +
  geom_point(data = selected_up,
             aes(x = avg_log2FC, y = log10fdr, color = Group),
             size = 1.5) +
  geom_point(data = selected_down,
             aes(x = avg_log2FC, y = log10fdr, color = Group),
             size = 1.5) +
  geom_vline(xintercept = 0, size = 0.5, color = "grey50", lty = 'dashed') +
  scale_color_manual(values = mycol) +
  xlab("avg_log2FC (Ketosis vs. Health)") +
  ylab("-Log10(adjusted P)") +
  scale_x_continuous(limits = c(-4, 4)) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 0.8),
        strip.text.x = element_text(size = 10, face = 'bold'))

# 添加TOP基因
library(ggrepel)
p + geom_text_repel(data = selected_up,
                    aes(x = avg_log2FC, y = log10fdr, label = gene, color = Group),
                    fontface = 'italic',
                    seed = 233,
                    size = 3,
                    box.padding = 0.4,       
                    point.padding = 0.3,     
                    force = 3,               
                    force_pull = 1,          
                    max.overlaps = 20,       
                    segment.alpha = 0.6,
                    segment.size = 0.3,
                    direction = "both"      
                    ) +
  geom_text_repel(data = selected_down,
                  aes(x = avg_log2FC, y = log10fdr, label = gene, color = Group),
                  fontface = 'italic',
                  seed = 233,
                  size = 3,
                  box.padding = 0.4,       
                  point.padding = 0.3,   
                  force = 3,              
                  force_pull = 1,       
                  max.overlaps = 20,     
                  segment.alpha = 0.6,
                  segment.size = 0.3,
                  direction = "both"     
                  )


# =============================================
# 🎨 Supplementary Figure 8
# =============================================

# 使用TCellSI包进行T细胞功能状态评分 首个全面评估T细胞状态的创新工具
# 目前准确评估T细胞的八种不同状态 这些状态包括：
  # ● 静止（Quiescence）
  # ● 调节（Regulating）
  # ● 增殖（Proliferation）
  # ● 辅助（Helper） 对应CD4
  # ● 细胞毒性（Cytotoxicity） 对应CD8
  # ● 初始耗竭（Progenitor exhaustion）
  # ● 终末耗竭（Terminal exhaustion）
  # ● 衰老（Senescence）  
  
# 加载单细胞数据和R包
library(TCellSI)
load("TT细胞鉴定.Rdata")  

# 作者推荐使用使用伪bulk进行打分 计算每个样本的状态分数 这可以减少单细胞数据中的dropout问题
expr_mat <- GetAssayData(TT, assay = "RNA", layer = "data")
expr_mat <- as.matrix(expr_mat)

# 从meta.data里整理出细胞ID和注释
annotation_data <- data.frame(UniqueCell_ID = colnames(TT),
                              annotation = paste(TT@meta.data$celltype, TT@meta.data$group, sep = "_"))
head(annotation_data)

# 创建pseudo-bulk
set.seed(123)
pseudo_bulk <- TCellSI::create_pseudo_bulk(annotation_data = annotation_data,
                                           expression_data = expr_mat,
                                           cluster_col = "annotation",
                                           cell_id_col = "UniqueCell_ID",
                                           n_clusters = length(unique(annotation_data$annotation)),
                                           factor = 5,
                                           sampling_rate = 0.6
                                           )
dim(pseudo_bulk)
pseudo_bulk[1:5, 1:5]
save(pseudo_bulk,file = "单细胞TCellSI评分pseudo_bulk.Rdata") 

# 对pseudo-bulk计算TCellSI
pb_scores <- TCellSI::TCSS_Calculate(pseudo_bulk)
save(pb_scores,file = "单细胞TCellSI评分pb_scores.Rdata") 
load("单细胞TCellSI评分pb_scores.Rdata")
pb_scores <- as.data.frame(t(pb_scores))
dim(pb_scores)
head(pb_scores)

# 处理数据框
pb_scores$sample <- rownames(pb_scores)
pb_scores$celltype <- sub("_.*", "", pb_scores$sample)
pb_scores$group <- sub(".*_(Health|Ketosis).*", "\\1", pb_scores$sample)
head(pb_scores)

# 转成长格式
library(reshape2)
df_long <- melt(
  pb_scores,
  id.vars = c("sample", "celltype", "group"))
head(df_long)

# Quiescence
df_cyto <- subset(df_long, variable == "Quiescence")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,   
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.325) +
  labs(x = NULL, y = "Quiescence score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")

# Regulating
df_cyto <- subset(df_long, variable == "Regulating")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,  
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.0027) +
  labs(x = NULL, y = "Regulating score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")

# Proliferation
df_cyto <- subset(df_long, variable == "Proliferation")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,  
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.16) +
  labs(x = NULL, y = "Proliferation score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")

# Helper
df_cyto <- subset(df_long, variable == "Helper")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4, 
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.041) +
  labs(x = NULL, y = "Helper score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")

# Cytotoxicity
df_cyto <- subset(df_long, variable == "Cytotoxicity")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,  
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.202) +
  labs(x = NULL, y = "Cytotoxicity score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")

# Progenitor_exhaustion
df_cyto <- subset(df_long, variable == "Progenitor_exhaustion")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,  
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.245) +
  labs(x = NULL, y = "Progenitor_exhaustion score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")

# Terminal_exhaustion
df_cyto <- subset(df_long, variable == "Terminal_exhaustion")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,   
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.091) +
  labs(x = NULL, y = "Terminal_exhaustion score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")

# Senescence
df_cyto <- subset(df_long, variable == "Senescence")
ggplot(df_cyto, aes(x = celltype, y = value, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,   
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.159) +
  labs(x = NULL, y = "Senescence score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif")


# =============================================
# 🎨 Fig. 4l
# =============================================

# 基因集打分 T细胞6种功能状态相关基因表达情况
Quiescence = c("BTG2", "KLF2", "FOXP1", "RUNX1", "MYC")
Helper1_2 = c("LTA", "TNF", "IL2", "IL18", "TBX21", "LTBR", "CXCR3", 
              "IFNG", "KLRD1", "HAVCR2", "DPP4", 
              "CCR5", "IFNGR1", "STAT1", "STAT4", "CCR4", "CCR6", 
              "GZMK", "IFNGR1", "PARP8", "PVRIG", "SLC4A10") 
Helper17 = c("TGFB1", "STAT3", "RORA", "RORC", "IL22", "IL17F", 
             "IL17A", "CCR6", "GZMK", "IFNGR1", "PARP8", "PVRIG", "SLC4A10") 
Fhelper = c("STAT3", "ICOS", "CXCR5", "BCL6", "CXCL13", "PDCD1", "BATF", "CXCR3", "KLRB1", 
            "LTA", "MAP4K1", "CCR4","CCR6", "GZMK", "IFNGR1", "PARP8", "PVRIG", "SLC4A10") 
Cytotoxicity = c("GZMB", "GZMK", "GZMA", "TIA1", "PRF1", "LAMP1",
                 "GNLY", "FASLG", "SLAMF7", "ZAP70", "CD69", "TNF")
Senescence = c("CD3G", "CD3D", "CD3E", "CD247", "KLRG1", "B3GAT1", "TP53", 
               "MAPK14", "MAPK8", "MAPK1", "CDKN2A")

load("TT细胞鉴定.Rdata")

Quiescence %in% rownames(TT)
Quiescence <- intersect(Quiescence, rownames(TT))
Quiescence <- list(Quiescence)  
  
Helper1_2 %in% rownames(TT) 
Helper1_2 <- intersect(Helper1_2, rownames(TT))
Helper1_2 <- list(Helper1_2) 

Helper17 %in% rownames(TT) 
Helper17 <- intersect(Helper17, rownames(TT))
Helper17 <- list(Helper17) 

Fhelper %in% rownames(TT) 
Fhelper <- intersect(Fhelper, rownames(TT))
Fhelper <- list(Fhelper) 

Cytotoxicity %in% rownames(TT) 
Cytotoxicity <- intersect(Cytotoxicity, rownames(TT))
Cytotoxicity <- list(Cytotoxicity)

Senescence %in% rownames(TT) 
Senescence <- intersect(Senescence, rownames(TT))
Senescence <- list(Senescence)
  
# 进行基因集打分
TT <- AddModuleScore(TT, features = Quiescence, name = "Quiescence")
TT <- AddModuleScore(TT, features = Helper1_2, name = "Helper1_2")
TT <- AddModuleScore(TT, features = Helper17, name = "Helper17")
TT <- AddModuleScore(TT, features = Fhelper, name = "Fhelper")
TT <- AddModuleScore(TT, features = Cytotoxicity, name = "Cytotoxicity")
TT <- AddModuleScore(TT, features = Senescence, name = "Senescence")
  
# 基因集打分存储在metadata信息中
colnames(TT@meta.data)   
  
# 提取需要修改的列名
col_names <- colnames(TT@meta.data)[-6:-1]
for (name in col_names) {
  new_name <- gsub("1", "", name) 
  colnames(TT@meta.data) <- gsub(name, new_name, colnames(TT@meta.data))
  }
colnames(TT@meta.data)  
  
# 保存结果
save(TT,file = "TT细胞鉴定8大功能打分版.Rdata") 

# Quiescence
ggplot(TT@meta.data, aes(x = celltype, y = Quiescence, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,  
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 1.6) +
  labs(x = NULL, y = "Quiescence score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) 

# Helper1_2
ggplot(TT@meta.data, aes(x = celltype, y = Helper_2, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4, 
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.5) +
  labs(x = NULL, y = "Helper1_2 score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) 

# Helper17
ggplot(TT@meta.data, aes(x = celltype, y = Helper7, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,   
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.52) +
  labs(x = NULL, y = "Helper17 score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) 

# Fhelper
ggplot(TT@meta.data, aes(x = celltype, y = Fhelper, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4, 
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.35) +
  labs(x = NULL, y = "Fhelper score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) 

# Cytotoxicity
ggplot(TT@meta.data, aes(x = celltype, y = Cytotoxicity, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,  
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 1.0) +
  labs(x = NULL, y = "Cytotoxicity score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) 
dev.off()

# Senescence
ggplot(TT@meta.data, aes(x = celltype, y = Senescence, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4, 
               position = position_dodge(width = 0.75),
               linewidth = 0.5) + 
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.65,
               linewidth = 0.4,
               outlier.shape = NA) + 
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 0.9) +
  labs(x = NULL, y = "Senescence score") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.border = element_rect(linewidth = 0.8, colour = "black")) 


# =============================================
# 🎨 Fig. 4m
# =============================================
  
# 基因集中关键基因表达热图或者气泡图 T细胞6种功能状态相关基因表达情况
load("TT细胞鉴定.Rdata")

# 设置对象的水平信息
TT@meta.data$group.celltype <- factor(
  TT@meta.data$group.celltype,
  levels = c(
    "Health_CD4 T","Ketosis_CD4 T",
    "Health_CD8 T","Ketosis_CD8 T",
    "Health_NKT","Ketosis_NKT",
    "Health_γδT","Ketosis_γδT"))

levels(TT@meta.data$group.celltype)
  
# 想要展示的基因
genes_use <- c("IL2","CD69",  # Activation
               
               "BACH2", "LEF1", "TCF7", "SELL", "CCR7", "FOXP1", # Quiescence
               
               "IL18",  "GATA3",  "IFNG", "STAT4", "CCR4", # Th1/2
               
               "GZMB", "CCL5", "NKG7", "TIA1", "GNLY", "TNF","CST7", "FASLG",    # Cytotoxicity 
               
               "RORC","CCR6","STAT3","AHR","IFNGR1"  # Th17
               ) 
genes_use <- intersect(genes_use, rownames(TT))
genes_use
  
# 画图
Seurat::DotPlot(TT, 
                features = genes_use,           
                scale = T,
                group.by = "group.celltype") +
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11, "PiYG"))(11))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10)) +
  theme(axis.text.y = element_text(face = "italic", size = 10)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid = element_blank(), 
        axis.line = element_blank()) + 
  coord_flip() +
  scale_size_continuous(range = c(0, 6))
  
  