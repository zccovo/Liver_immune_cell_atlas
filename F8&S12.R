
# F8/S12: Clinical ketosis drives lipid-associated macrophage expansion and MS4A7-mediated inflammatory activation.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 Mono/Mac细胞亚群数据处理
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

# 提取出里面的单核巨噬细胞
Mac <- subset(seurat_object, celltype %in% c("Mono", "Mac"))

library(clustree)
clustree(Mac, prefix = "RNA_snn_res.")
Idents(Mac) <- "RNA_snn_res.0.9"
Mac$seurat_clusters <- Mac@active.ident

# 进行UMAP非线性降维
Mac <- RunUMAP(Mac, dims = 1:30) 

# UMAP图
Seurat::DimPlot(Mac, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(Mac, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

# 保存数据
save(Mac, file = "Mac原始数据.Rdata") 
rm(list = ls())
load("Mac原始数据.Rdata") 

# 加载R包
library(decontX)
load("Mac原始数据.Rdata")
seurat_object = Mac
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])

# 运行自定义封装函数
ks_runDecontX <- function(seurat_obj,  
                          idents,     
                          seed = 1){  
  counts <- GetAssayData(object = seurat_obj, assay = 'RNA', layer = 'counts')
  Idents(seurat_obj) <- idents
  clusters <- Idents(seurat_obj)
  x <- counts[rowSums(counts)>0,]
  decon <- decontX(x,              
                   z = clusters,  
                   verbose = TRUE, 
                   seed = seed)   
  seurat_obj[["origCounts"]] <- CreateAssayObject(counts = counts)
  newCounts <- decon$decontXcounts
  newCounts <- rbind(newCounts, counts[rowSums(counts)==0,])[rownames(counts),]
  seurat_obj[["newCounts"]] <- CreateAssayObject(counts = as(round(newCounts), "sparseMatrix"))
  seurat_obj$estConp <- decon$contamination 
  return(seurat_obj)
  }

adj_scRNA <- ks_runDecontX(seurat_obj = seurat_object,
                           idents = "seurat_clusters")

# 保存结果
save(adj_scRNA,file = "./Mac之decontX分析结果.Rdata")

# 可视化UMAP降维图中的污染分布
library(ggplot2)
FeaturePlot(adj_scRNA, 
            features = 'estConp') +       
  scale_color_viridis_c(option = "B") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  xlab('UMAP_1')+
  ylab('UMAP_2')

# 根据污染分数删除细胞
adj_scRNA1 = adj_scRNA[,adj_scRNA$estConp < 0.2]
Seurat::DotPlot(adj_scRNA1, features = c("GPNMB", "TREM2", "SPP1", "MS4A7", "APOE", "ALOX5"), group.by = "group") + RotatedAxis() 

# 对比去除污染前后的结果
DimPlot(adj_scRNA1, label = T) + DimPlot(seurat_object, label = T)

# 重命名
Mac = adj_scRNA1

# UMAP图
Seurat::DimPlot(Mac, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(Mac, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)
library(clustree)
clustree(Mac, prefix = "RNA_snn_res.") 
Idents(Mac) <- "RNA_snn_res.0.6"
Mac$seurat_clusters <- Mac@active.ident

# 保存去环境RNA后的数据
save(Mac, file = "Mac去环境RNA后.Rdata") 
rm(list = ls())
load("Mac去环境RNA后.Rdata")

# 分群数量
Mac = FindClusters(Mac, resolution = c(seq(0.1, 2, 0.1)))
library(clustree)
clustree(Mac, prefix = "RNA_snn_res.") 
Idents(Mac) <- "RNA_snn_res.0.6"
Mac$seurat_clusters <- Mac@active.ident

# UMAP图
Seurat::DimPlot(Mac, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(Mac, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

# 保存结果
save(Mac, file = "Mac最终版.Rdata") 

# 细胞鉴定之前 先用FindAllMarkers()函数寻找Cluster之间的差异基因
dif<-FindAllMarkers(Mac, 
                    group.by = Mac@meta.data$seurat_clusters, 
                    logfc.threshold = log2(1.2),                        
                    min.pct = 0.2,                                      
                    only.pos = T                                        
                    )        
dif$pct_diff <- dif$pct.1 - dif$pct.2
table(dif$cluster)                    
dif<-dif %>%
  group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
save(dif, file = "Mac整体cluster差异基因.Rdata")

# 0/1 GPNMB_Mac
FeaturePlot(Mac, features = "CLEC4F", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "GPNMB", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "APOE", reduction = "umap", pt.size = 0.1)  
Seurat::DotPlot(Mac, features = c("GPNMB", "CLEC4F", "APOE"), group.by = "seurat_clusters") + RotatedAxis()
# 7/10 KANK1_Mac
FeaturePlot(Mac, features = "KANK1", reduction = "umap", pt.size = 0.1) # 细胞骨架/迁移调控
FeaturePlot(Mac, features = "TNF", reduction = "umap", pt.size = 0.1) # 炎症
FeaturePlot(Mac, features = "RGS1", reduction = "umap", pt.size = 0.1) # 白细胞迁移调控
FeaturePlot(Mac, features = "PTGS2", reduction = "umap", pt.size = 0.1) # COX-2 急性炎症标志
Seurat::DotPlot(Mac, features = c("KANK1", "TNF", "RGS1", "PTGS2"), group.by = "seurat_clusters") + RotatedAxis()
# 9 CXCL10_Mac
FeaturePlot(Mac, features = "CXCL11", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "CXCL10", reduction = "umap", pt.size = 0.1)
FeaturePlot(Mac, features = "CXCL9", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(Mac, features = c("CXCL11", "CXCL10", "CXCL9"), group.by = "seurat_clusters") + RotatedAxis()
# 2/3/4 MARCO_Mac
FeaturePlot(Mac, features = "MARCO", reduction = "umap", pt.size = 0.1)  
Seurat::DotPlot(Mac, features = c("MARCO"), group.by = "seurat_clusters") + RotatedAxis()
# 11 LTB4R_Mac
FeaturePlot(Mac, features = "ALOX5", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "CX3CR1", reduction = "umap", pt.size = 0.1) # 单核/巨噬迁入特征 提示这个群可能是单核来源或巡逻型巨噬
FeaturePlot(Mac, features = "LTB4R", reduction = "umap", pt.size = 0.1) # 白三烯B4受体 炎症趋化因子 调控免疫细胞迁移
FeaturePlot(Mac, features = "CD1E", reduction = "umap", pt.size = 0.1) # 脂质抗原呈递给T细胞
FeaturePlot(Mac, features = "TREM2", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "FABP5", reduction = "umap", pt.size = 0.1)
FeaturePlot(Mac, features = "B4GALT6", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "CREB5", reduction = "umap", pt.size = 0.1)
Seurat::DotPlot(Mac, features = c("ALOX5", "CX3CR1", "LTB4R"), group.by = "seurat_clusters") + RotatedAxis()
# 8 OSM_Mono 促炎/趋化 招募中性粒细胞
FeaturePlot(Mac, features = "OSM", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "THBS1", reduction = "umap", pt.size = 0.1) 
Seurat::DotPlot(Mac, features = c("OSM", "THBS1", "CXCL2", "CXCL3", "CXCL5", "CXCL8"), group.by = "seurat_clusters") + RotatedAxis()
# 5 CCR2_Mono
FeaturePlot(Mac, features = "CCR2", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "ABCA10", reduction = "umap", pt.size = 0.1) 
FeaturePlot(Mac, features = "RGS14", reduction = "umap", pt.size = 0.1) 
Seurat::DotPlot(Mac, features = c("CCR2", "ABCA10", "RGS14"), group.by = "seurat_clusters") + RotatedAxis()
# 6 HES4_Mac
FeaturePlot(Mac, features = "HES4", reduction = "umap", pt.size = 0.1)  
FeaturePlot(Mac, features = "CX3CR1", reduction = "umap", pt.size = 0.1)  
Seurat::DotPlot(Mac, features = c("HES4","CX3CR1"), group.by = "seurat_clusters") + RotatedAxis()

# 为分群重新指定细胞类型 
new.cluster.ids <- c("GPNMB_Mac",
                     "GPNMB_Mac",
                     "MARCO_Mac",
                     "MARCO_Mac",
                     "MARCO_Mac",
                     "CCR2_Mono",
                     "HES4_Mac",
                     "KANK1_Mac",
                     "OSM_Mono",
                     "CXCL10_Mac",
                     "KANK1_Mac",
                     "LTB4R_Mac") 
new.cluster.ids
names(new.cluster.ids) 
levels(Mac)
names(new.cluster.ids) <- levels(Mac) 
names(new.cluster.ids)
new.cluster.ids
Mac <- RenameIdents(Mac, new.cluster.ids)

Mac[["celltype"]] <- Idents(Mac) 
table(Mac@meta.data[["orig.ident"]])
unique(Mac@meta.data[["orig.ident"]])
Mac[["group"]]<- c(rep("Health", 7252), rep("Ketosis", 22678))
Mac[["group.celltype"]]<-paste(Mac$group, Idents(Mac), sep = '_') 
table(Mac@meta.data[["orig.ident"]])
table(Mac@meta.data[["group"]])
table(Mac@meta.data[["celltype"]])
table(Mac@meta.data[["group.celltype"]])
table(Mac@meta.data[["seurat_clusters"]])

# 细胞水平信息
Idents(Mac) <- factor(Idents(Mac),
                     levels = c("CCR2_Mono", "OSM_Mono", 
                                "HES4_Mac", "LTB4R_Mac", "KANK1_Mac",
                                "CXCL10_Mac", "GPNMB_Mac", "MARCO_Mac"))
Mac[["celltype"]] <- Idents(Mac) 

# 绘制总umap图
Seurat::DimPlot(Mac, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(Mac, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 保存工作空间
save(Mac,file = "Mac细胞鉴定.Rdata") 
rm(list = ls())
load("Mac细胞鉴定.Rdata")

# 对注释后细胞亚群 进行差异基因分析
diff <- FindAllMarkers(Mac, 
                      group.by = Mac@active.ident, 
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
save(diff, file = "Mac整体细胞类型差异基因.Rdata")    


# =============================================
# 🎨 Fig. 8a
# =============================================

load("Mac细胞鉴定.Rdata")
Seurat::DimPlot(Mac,
                group.by = "celltype",
                cols = c("#e64b35",
                         "#4dbbd5",
                         "#00a087",
                         "#3e7cbd",
                         "#91d1c2",
                         "#b09c85",
                         "#8491b4",
                         "#f39b7f"), 
                pt.size = 0.5,
                label = T) +
  NoLegend()+ 
  labs(title = NULL)  


# =============================================
# 🎨 Fig. 8b
# =============================================

# 加载数据
load("Mac细胞鉴定.Rdata")
plot_density(Mac, 
             reduction = "umap",
             features = c("C1QA"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("FCN1"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        limits = c(0, 0.075),
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("CCR2"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("OSM"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A",
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("HES4"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("LTB4R"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("KANK1"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("CXCL10"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("GPNMB"),
             adjust = 1, 
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())

plot_density(Mac, 
             reduction = "umap",
             features = c("MARCO"),
             adjust = 1,  
             raster = T, 
             size = 1.5) +
  scale_color_viridis_c(option = "A", 
                        oob = scales::squish) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())


# =============================================
# 🎨 Fig. 8c
# =============================================

# 加载数据
load("Mac细胞鉴定.Rdata")
table(Mac@meta.data[["celltype"]])

# 随机抽样
table(Mac@meta.data[["celltype"]])
cell_types <- unique(Mac$celltype)

# 创建一个空的 list 来存储每个细胞类型的子集
subset_cells <- list() 

# 对每个细胞类型进行处理
for (cell_type in cell_types) {
  cells_of_type <- WhichCells(Mac, expression = celltype == cell_type)
  if (length(cells_of_type) >= 150) {
    selected_cells <- sample(cells_of_type, 150)                     
  } else {
    selected_cells <- cells_of_type                   
  }
  subset_cells[[cell_type]] <- selected_cells   
}  

# 合并所有子集细胞
subset_all_cells <- unlist(subset_cells)

# 根据选中的细胞创建新的 Seurat 对象
new_seurat_object <- subset(Mac, cells = subset_all_cells) 

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
pheatmap(data,scale = "none",cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,show_rownames = FALSE,
         annotation_col = celltype,     
         annotation_row = gene.anno,   
         annotation_names_row = FALSE, 
         color = colorRampPalette(c("#040509","#608fe4", "#ffd700"))(100)
         )


# =============================================
# 🎨 Fig. 8d
# =============================================

# 首先进行亚群间的差异分析
dif <- FindAllMarkers(Mac, 
                      group.by = Mac@active.ident, 
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

table(Mac@active.ident)
CCR2_Mono  <- subset(sigposDEG.all, cluster=='CCR2_Mono') 
OSM_Mono <- subset(sigposDEG.all, cluster=='OSM_Mono')
HES4_Mac <- subset(sigposDEG.all, cluster=='HES4_Mac')
LTB4R_Mac  <- subset(sigposDEG.all, cluster=='LTB4R_Mac')  
KANK1_Mac  <- subset(sigposDEG.all, cluster=='KANK1_Mac') 
CXCL10_Mac <- subset(sigposDEG.all, cluster=='CXCL10_Mac')
GPNMB_Mac <- subset(sigposDEG.all, cluster=='GPNMB_Mac')
MARCO_Mac  <- subset(sigposDEG.all, cluster=='MARCO_Mac')  

list <- list(CCR2_Mono, OSM_Mono, HES4_Mac, LTB4R_Mac, KANK1_Mac, CXCL10_Mac, GPNMB_Mac, MARCO_Mac)
names(list)[1:8] <- c("CCR2_Mono", "OSM_Mono", "HES4_Mac", "LTB4R_Mac","KANK1_Mac", "CXCL10_Mac", "GPNMB_Mac", "MARCO_Mac")
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
CCR2_Mono <- read.csv("GO_new_CCR2_Mono.CSV", row.names = 1) 
OSM_Mono <- read.csv("GO_new_OSM_Mono.CSV", row.names = 1)       
HES4_Mac <- read.csv("GO_new_HES4_Mac.CSV", row.names = 1)  
LTB4R_Mac <- read.csv("GO_new_LTB4R_Mac.CSV", row.names = 1)   
KANK1_Mac <- read.csv("GO_new_KANK1_Mac.CSV", row.names = 1)  
CXCL10_Mac <- read.csv("GO_new_CXCL10_Mac.CSV", row.names = 1)   
GPNMB_Mac <- read.csv("GO_new_GPNMB_Mac.CSV", row.names = 1)  
MARCO_Mac <- read.csv("GO_new_MARCO_Mac.CSV", row.names = 1)   

# 为每个细胞群体添加标签
CCR2_Mono$group <- "CCR2_Mono"
OSM_Mono$group <- "OSM_Mono"
HES4_Mac$group <- "HES4_Mac"
LTB4R_Mac$group <- "LTB4R_Mac"
KANK1_Mac$group <- "KANK1_Mac"
CXCL10_Mac$group <- "CXCL10_Mac"
GPNMB_Mac$group <- "GPNMB_Mac"
MARCO_Mac$group <- "MARCO_Mac"

# 选择TOP通路
# CCR2_Mono
select_CCR2_Mono = c("myeloid leukocyte migration",  
                     "canonical NF-kappaB signal transduction", 
                     "toll-like receptor signaling pathway",
                     "activation of innate immune response")
# OSM_Mono
select_OSM_Mono = c("inflammatory response",  
                    "leukocyte chemotaxis", 
                    "cytokine-mediated signaling pathway",
                    "cellular response to lipid")
# HES4_Mac
select_HES4_Mac = c("actin cytoskeleton organization", 
                    "positive regulation of RNA metabolic process",
                    "regulation of RNA splicing",
                    "regulation of cell cycle")
# LTB4R_Mac
select_LTB4R_Mac = c("MHC class II protein complex assembly",  
                     "antigen processing and presentation of peptide antigen via MHC class II", 
                     "antigen processing and presentation of exogenous antigen",
                     "immune response-activating signaling pathway")
# KANK1_Mac
select_KANK1_Mac = c("phosphatidylinositol 3-kinase/protein kinase B signal transduction",  
                     "regulation of intracellular signal transduction", 
                     "regulation of GTPase activity",
                     "small GTPase-mediated signal transduction")
# CXCL10_Mac
select_CXCL10_Mac = c("chemokine-mediated signaling pathway",  
                     "positive regulation of lymphocyte activation",
                     "response to type II interferon",
                     "regulation of leukocyte activation")
# GPNMB_Mac
select_GPNMB_Mac = c("fatty acid metabolic process", 
                     "cellular lipid metabolic process", 
                     "monocarboxylic acid metabolic process",
                     "lipid catabolic process")
# MARCO_Mac
select_MARCO_Mac = c("positive regulation of endocytosis",  
                     "mitochondrion organization", 
                     "regulation of cell-cell adhesion",
                     "actin filament organization")

# 选择每个亚群的通路
CCR2_Mono <- CCR2_Mono[CCR2_Mono$Description %in% select_CCR2_Mono,]
OSM_Mono <- OSM_Mono[OSM_Mono$Description %in% select_OSM_Mono,]
HES4_Mac <- HES4_Mac[HES4_Mac$Description %in% select_HES4_Mac,]
LTB4R_Mac <- LTB4R_Mac[LTB4R_Mac$Description %in% select_LTB4R_Mac,] 
KANK1_Mac <- KANK1_Mac[KANK1_Mac$Description %in% select_KANK1_Mac,]
CXCL10_Mac <- CXCL10_Mac[CXCL10_Mac$Description %in% select_CXCL10_Mac,]
GPNMB_Mac <- GPNMB_Mac[GPNMB_Mac$Description %in% select_GPNMB_Mac,]
MARCO_Mac <- MARCO_Mac[MARCO_Mac$Description %in% select_MARCO_Mac,] 

# 生成新的P值列
CCR2_Mono$`-log10pvalue` <- -log10(CCR2_Mono$pvalue)
OSM_Mono$`-log10pvalue` <- -log10(OSM_Mono$pvalue)
HES4_Mac$`-log10pvalue` <- -log10(HES4_Mac$pvalue)
LTB4R_Mac$`-log10pvalue` <- -log10(LTB4R_Mac$pvalue) 
KANK1_Mac$`-log10pvalue` <- -log10(KANK1_Mac$pvalue)
CXCL10_Mac$`-log10pvalue` <- -log10(CXCL10_Mac$pvalue)
GPNMB_Mac$`-log10pvalue` <- -log10(GPNMB_Mac$pvalue)
MARCO_Mac$`-log10pvalue` <- -log10(MARCO_Mac$pvalue)

# 合并所有数据
all <- rbind(CCR2_Mono, OSM_Mono, HES4_Mac, LTB4R_Mac, KANK1_Mac, CXCL10_Mac, GPNMB_Mac, MARCO_Mac)
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
My_levels <- c("CCR2_Mono", "OSM_Mono", "HES4_Mac", "LTB4R_Mac", "KANK1_Mac", "CXCL10_Mac", "GPNMB_Mac", "MARCO_Mac")
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
  scale_size_continuous(range = c(3, 6))


# =============================================
# 🎨 Fig. 8e
# =============================================

# 分组的细胞比例图
load("Mac细胞鉴定.Rdata")  
library(dplyr)
library(ggplot2)
library(gtools)
library(ggalluvial)

# 准备细胞比例输入数据
prop_df <- Mac@meta.data %>%
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
  geom_flow(width = 0.6, alpha = 0.3, knot.pos = 0.3) +  
  geom_col(width = 0.6,) +  
  scale_y_continuous(expand = c(0, 0)) +   
  scale_fill_manual(values = c("#e64b35",
                               "#4dbbd5",
                               "#00a087",
                               "#3e7cbd",
                               "#91d1c2",
                               "#b09c85",
                               "#8491b4",
                               "#f39b7f")) +  
  xlab("") + 
  ylab("Cell proportion") +  
  theme_classic() +  
  theme(axis.text.x = element_text(size = 16),    
        axis.text.y = element_text(size = 16),    
        axis.title.y = element_text(size = 16),    
        legend.text=element_text(size = 16),    
        legend.title=element_text(size = 16)
        )


# =============================================
# 🎨 Fig. 8f
# =============================================

# 计算每个细胞群体在总细胞中的比例
load("Mac细胞鉴定.Rdata")
load("CD45去双细胞后.Rdata")
table(seurat_object$group)
table(Idents(Mac), Mac$group)
# 获取每个亚群在每个样本中的细胞数
Mac_subgroup_counts <- table(Idents(Mac), Mac$group)
# 获取每个样本中的细胞数
total_cells <- table(seurat_object$group)
total_cells <- total_cells[match(colnames(Mac_subgroup_counts), names(total_cells))]
total_cells_matrix <- matrix(total_cells, nrow = nrow(Mac_subgroup_counts), ncol = ncol(Mac_subgroup_counts), byrow = TRUE)
# 计算每个亚群在每个样本中的占比
Mac_subgroup_percentage <- Mac_subgroup_counts / total_cells_matrix * 100
# 查看比例
Mac_subgroup_percentage
# 转换为 matrix
mac_matrix <- as.matrix(Mac_subgroup_percentage)
# 计算 Ketosis / Health 倍数
fold_change <- mac_matrix[, "Ketosis"] / mac_matrix[, "Health"]
fold_change
# table转data.frame
df <- as.data.frame.matrix(Mac_subgroup_percentage)
# 计算FC和log2FC
df$Subgroup <- rownames(df)
df$FoldChange <- df$Ketosis / df$Health
df$log2FC <- log2(df$FoldChange)
# 保持原始亚群顺序
df$Subgroup <- factor(df$Subgroup, levels = rownames(Mac_subgroup_percentage))

# 自定义颜色
my_colors <- c("#e64b35","#4dbbd5","#00a087","#3e7cbd",
               "#91d1c2","#b09c85","#8491b4","#f39b7f")

# 画图
ggplot(df, aes(x = Subgroup, y = log2FC, fill = Subgroup)) +
  geom_col(color = "black",width = 0.7,linewidth = 0.8) +
  geom_hline(yintercept = 0,color = "black",linewidth = 0.8) +
  scale_fill_manual(values = my_colors) +
  labs(x = NULL,y = expression(Log[2]~"(Ketosis/Health)")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 13),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(linewidth = 0.7))


# =============================================
# 🎨 Fig. 8g
# =============================================

# GPNMB_Mac 频率分布图
set.seed(123)
ggplot(data=Cellratio[Cellratio$Celltype=="GPNMB_Mac",], aes(group, Freq)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(NA, 25) + 
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

Cellratio[Cellratio$Celltype == "GPNMB_Mac",]
shapiro.test(Cellratio[Cellratio$Celltype=="GPNMB_Mac",]$Freq[1:3])
shapiro.test(Cellratio[Cellratio$Celltype=="GPNMB_Mac",]$Freq[4:8]) 
leveneTest(Freq ~ group, data = Cellratio[Cellratio$Celltype=="GPNMB_Mac",], center = "mean") 
t.test(Freq ~ group, data = Cellratio[Cellratio$Celltype=="GPNMB_Mac",], var.equal = TRUE) # 0.02933 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 8h
# =============================================

# GPNMB+/DAPI (Fold)
group <- c(rep("Health", 4), rep("Ketosis", 4))
value <- c(0.711160965,
           0.867077958,
           0.927650028,
           1.494111049,
           2.765002804,
           2.969153113,
           4.325294448,
           4.66741447)
df <- data.frame(Group = group, Value = value)
set.seed(123)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(0, 5) + 
  ylab(NULL) +
  xlab(NULL)  +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black"),  
        axis.title.y = element_text(size = 19, color = "black"),  
        legend.position = "none") +
  stat_summary(fun = mean,
               geom = "errorbar",
               size = 0.5,
               width = 0.1,
               color = "black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x)) +
  stat_summary(fun = mean,         
               geom = "crossbar",  
               size = 0.2,        
               color = "black", 
               width = 0.15)     

shapiro.test(df$Value[df$Group == "Health"])
shapiro.test(df$Value[df$Group == "Ketosis"])
leveneTest(Value ~ Group, data = df, center = "mean") 
t.test(Value ~ Group, data = df) # 0.007291 独立两组-正态-方差非齐


# =============================================
# 🎨 Supplementary Figure 12f
# =============================================

load("Mac细胞鉴定.Rdata")

# 整体UMAP图显示MS4A7主要在巨噬细胞中表达
FeaturePlot(seurat_object,
            features = "MS4A7",
            reduction = "umap",
            pt.size = 0.5) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())


# =============================================
# 🎨 Supplementary Figure 12g
# =============================================

# 两组分开 用红色点
FeaturePlot(subset(subset(Mac, celltype %in% c("CCR2_Mono", "OSM_Mono"), invert = TRUE), group == "Health"),
            features = "MS4A7",
            reduction = "umap",
            pt.size = 0.4,
            cols = c("#d3d3d3", "red"),
            order = TRUE) +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 

FeaturePlot(subset(subset(Mac, celltype %in% c("CCR2_Mono", "OSM_Mono"), invert = TRUE), group == "Ketosis"),
            features = "MS4A7",
            reduction = "umap",
            pt.size = 0.4,
            cols = c("#d3d3d3", "red"),
            order = TRUE,
            max.cutoff = "q99") +
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank())+  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            color = "black", fill = NA, size = 1) 


# =============================================
# 🎨 Supplementary Figure 12h
# =============================================

# 伪BulkRNA分析显示 MS4A7在CD45免疫细胞中显著上调
load("~/奶牛肝脏解离单细胞3版/Bulk原始矩阵.rdata")     
load("~/奶牛肝脏解离单细胞3版/Bulk标准化矩阵.Rdata") 
load("~/奶牛肝脏解离单细胞3版/DEG_deseq2.Rdata")
DEG_deseq2[rownames(DEG_deseq2) == "MS4A7",] # 3.482202倍

# 看看CPM标准化基因表达
count_data[rownames(count_data) == "MS4A7",]
dat[rownames(dat) == "MS4A7",]

# 比较目标基因表达
boxplot(dat,las = 2) 
dat = log2(dat + 1) 
dat[1:4, 1:4]
boxplot(dat,las = 2) 
dat = t(dat) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(.) %>% 
  mutate(.,Group = Group) 

# 基因表达散点图
set.seed(123)
ggplot(dat, aes(x = Group, y = MS4A7)) +
  geom_jitter(aes(fill = Group), position = position_jitter(0), shape = 21, size = 4.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),  # 👈这里
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

# 计算P值
shapiro.test(dat[,colnames(dat)=="MS4A7"][1:3]) 
shapiro.test(dat[,colnames(dat)=="MS4A7"][4:8]) 
leveneTest(MS4A7 ~ Group, data = dat, center = "mean") 
t.test(MS4A7 ~ Group, data = dat, var.equal = TRUE) # 0.004835 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 8j
# =============================================

# 关键基因的表达情况
Seurat::DotPlot(Mac, features = c("GPNMB", "TREM2", "SPP1", "LGALS3",       # LAM
                                  "APOE", "FABP4", "FABP5", "LPL", "PLIN2", # 脂代谢
                                  "CD63","CTSB",                            # 溶酶体
                                  "IL1B", "TNF", "IL6", "IL18", "IL10",     # 炎症相关细胞因子
                                  "MS4A7"), 
                group.by = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 13)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid = element_blank(), 
        axis.line = element_blank()) +  
  scale_size(range = c(1, 6.5)) 


# =============================================
# 🎨 Fig. 8l
# =============================================

# 证明MS4A7是酮病期间高度诱导的巨噬细胞富集因子 只看巨噬细胞的小提琴图
VlnPlot(seurat_object, features = c("MS4A7"), group.by = "group", pt.size = 0)
VlnPlot(Mac, features = c("MS4A7"), group.by = "group", pt.size = 0)
VlnPlot(subset(Mac, celltype %in% c("CCR2_Mono", "OSM_Mono"), invert = TRUE), 
        features = c("MS4A7"), 
        group.by = "group", pt.size = 0.01) +  
  scale_fill_manual(values = c("Health" = "#4cdafe", "Ketosis" = "#ff6362"))+
  stat_summary(fun = mean,            
               geom = "point",        
               shape = 18,          
               size = 5,             
               color = "#f8f65c")     
dev.off()
expr_df <- FetchData(subset(Mac, celltype %in% c("CCR2_Mono", "OSM_Mono", invert = TRUE)), vars = c("MS4A7", "group"))
wilcox.test(MS4A7 ~ group, data = expr_df)$p.value # 5.32548e-18


# =============================================
# 🎨 Fig. 8k
# =============================================

# 基因集打分显示酮病组巨噬细胞炎症反应增强 
load("CD45去双细胞后打分版.Rdata")
expr_df <- FetchData(subset(seurat_object, celltype %in% c("Mac")),
                     vars = c("inflammatory_response", "group"))

ggplot(expr_df, aes(x = group, y = inflammatory_response, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0)) +  
  scale_fill_manual(values = c("Health" = "#4cdafe", "Ketosis" = "#ff6362")) + 
  scale_y_continuous(limits = c(NA, 0.2)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"),
        legend.position = "none") 例

df <- FetchData(subset(seurat_object, celltype %in% c("Mac")),vars = c("inflammatory_response", "group"))
head(df)
wilcox.test(inflammatory_response ~ group, data = df)$p.value # P-value = 3.500089e-21


# =============================================
# 🎨 Fig. 8m
# =============================================

# 在总体巨噬细胞 或者 仅在酮病巨噬细胞中 表达MS4A7的细胞 炎症评分更强
Mac_subset <- subset(seurat_object, subset = celltype == "Mac")
Mac_subset$MS4A7_status <- ifelse(FetchData(Mac_subset, vars = "MS4A7")$MS4A7 > 0, 
                                  "MS4A7_pos", "MS4A7_neg")
table(Mac_subset@meta.data[["MS4A7_status"]])

VlnPlot(Mac_subset,
        features = "inflammatory_response",
        group.by = "MS4A7_status",
        pt.size = 0)+
  scale_fill_manual(values = c("MS4A7_neg" = "#8dd3c7", "MS4A7_pos" = "#b3de69")) +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.x =  element_blank(),
        axis.text.x =  element_blank(), 
        legend.position = "none") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "red", color = "red") 

expr_df <- FetchData(Mac_subset, vars = c("inflammatory_response", "MS4A7_status"))
wilcox.test(inflammatory_response ~ MS4A7_status, data = expr_df)$p.value # 4.75999e-45

Mac_subset$group_MS4A7 <- paste0(Mac_subset$group, "_", Mac_subset$MS4A7_status)

VlnPlot(Mac_subset,
        features = "inflammatory_response",
        group.by = "group_MS4A7",
        pt.size = 0.0) +
  scale_fill_manual(values = c("Health_MS4A7_neg" = "#bc80bd", "Health_MS4A7_pos" = "#80b1d3", "Ketosis_MS4A7_neg" = "#fb8072", "Ketosis_MS4A7_pos" = "#fdb462")) +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.x =  element_blank(),
        axis.text.x =  element_blank(),  # 横坐标字体显示并倾斜45度
        legend.position = "none") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "red", color = "red") 

expr_df <- FetchData(Mac_subset, vars = c("inflammatory_response", "group_MS4A7"))
kruskal.test(inflammatory_response ~ group_MS4A7, data = expr_df) 
dunn.test(x = expr_df$inflammatory_response, g = expr_df$group_MS4A7, method = "bh", kw = TRUE) 


# =============================================
# 🎨 Supplementary Figure 12e
# =============================================

load("~/可分析奶牛公共数据/33189277/标准化后的表达矩阵.Rdata")  
load("~/可分析奶牛公共数据/33189277/Deseq2差异分析结果.Rdata") 
colnames(标准化后的表达矩阵)
mat = 标准化后的表达矩阵[,-c(22,36)]
colnames(mat)
Group <- c(rep("After", 17), rep("Before", 17))
Group <- factor(Group, levels = c("Before", "After")) 
Group

# 存储存储整理好的数据
save(mat, Group, file = "产前产后配对标准化表达矩阵.Rdata") 
load("产前产后配对标准化表达矩阵.Rdata")

# 用标准化后的TPM数据绘制样本表达总体分布图
boxplot(mat, las = 2)  
dat = log2(mat + 1) 
boxplot(dat, las = 2)  
dat[1:4, 1:4]

# 转置成绘图格式
dat = t(dat) %>%               
  as.data.frame(.) %>%          
  rownames_to_column(.) %>%     
  mutate(.,Group = Group)      

# GPNMB
df <- dat %>%
  dplyr::select(Sample = rowname, GPNMB) %>%
  mutate(Group = Group) 
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 
set.seed(123)
ggplot(df, aes(x = Group, y = GPNMB, fill = Group)) +
  geom_boxplot(width = 0.5, fill = NA, color = "black", outlier.shape = NA)+
  geom_point(aes(fill = Group), shape = 21, color = "grey", size = 3, alpha = 0.85) +
  geom_line(aes(group = ID), color = "grey", alpha = 0.6) +
  scale_fill_manual(values = c("Before" = "#5b189c", "After" = "#ffb702")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5))

shapiro.test(df$GPNMB[df$Group == "Before"])
shapiro.test(df$GPNMB[df$Group == "After"])
wilcox.test(df$GPNMB[df$Group == "Before"], df$GPNMB[df$Group == "After"], paired = TRUE) # 1.526e-05

# TREM2
df <- dat %>%
  dplyr::select(Sample = rowname, TREM2) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))  
set.seed(123)
ggplot(df, aes(x = Group, y = TREM2, fill = Group)) +
  geom_boxplot(width = 0.5, fill = NA, color = "black", outlier.shape = NA)+
  geom_point(aes(fill = Group), shape = 21, color = "grey", size = 3, alpha = 0.85) +
  geom_line(aes(group = ID), color = "grey", alpha = 0.6) +
  scale_fill_manual(values = c("Before" = "#5b189c", "After" = "#ffb702")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

shapiro.test(df$TREM2[df$Group == "Before"])
shapiro.test(df$TREM2[df$Group == "After"])
t.test(df$TREM2[df$Group == "Before"], df$TREM2[df$Group == "After"], paired = TRUE) # 5.511e-05

# SPP1
df <- dat %>%
  dplyr::select(Sample = rowname, SPP1) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))  
set.seed(123)
ggplot(df, aes(x = Group, y = SPP1, fill = Group)) +
  geom_boxplot(width = 0.5, fill = NA, color = "black", outlier.shape = NA)+
  geom_point(aes(fill = Group), shape = 21, color = "grey", size = 3, alpha = 0.85) +
  geom_line(aes(group = ID), color = "grey", alpha = 0.6) +
  scale_fill_manual(values = c("Before" = "#5b189c", "After" = "#ffb702")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5))

shapiro.test(df$SPP1[df$Group == "Before"])
shapiro.test(df$SPP1[df$Group == "After"])
t.test(df$SPP1[df$Group == "Before"], df$SPP1[df$Group == "After"], paired = TRUE) # 0.000169

# LGALS3
df <- dat %>%
  dplyr::select(Sample = rowname, LGALS3) %>%
  mutate(Group = Group) 
df
df <- df %>%
  mutate(ID = rep(1:17, 2))  
set.seed(123)
ggplot(df, aes(x = Group, y = LGALS3, fill = Group)) +
  geom_boxplot(width = 0.5, fill = NA, color = "black", outlier.shape = NA)+
  geom_point(aes(fill = Group), shape = 21, color = "grey", size = 3, alpha = 0.85) +
  geom_line(aes(group = ID), color = "grey", alpha = 0.6) +
  scale_fill_manual(values = c("Before" = "#5b189c", "After" = "#ffb702")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5))

shapiro.test(df$LGALS3[df$Group == "Before"])
shapiro.test(df$LGALS3[df$Group == "After"])
wilcox.test(df$LGALS3[df$Group == "Before"], df$LGALS3[df$Group == "After"], paired = TRUE) # 0.0003815


# =============================================
# 🎨 Supplementary Figure 12b, c
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
load("Mac细胞鉴定.Rdata")

# 将Seurat对象转为SingleCellExperiment对象
sce <- as.SingleCellExperiment(Mac, assay = "RNA") 

# 运行主函数
sce_slingshot1 <- slingshot(data = sce,                   
                            reducedDim = 'UMAP',           
                            clusterLabels = sce$celltype,   
                            start.clus = NULL,      
                            end.clus = NULL,                
                            approx_points = 150           
                            )

# 保存slingshot分析结果
save(sce_slingshot1, file = "Mac拟时之Slingshot分析.Rdata")  
load("Mac拟时之Slingshot分析.Rdata")

# SlingshotDataSet函数查看轨迹信息
SlingshotDataSet(sce_slingshot1)
dim(slingPseudotime(sce_slingshot1)) 
slingPseudotime(sce_slingshot1)[1:2,1:2] 

# 细胞类型为基础进行可视化
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$celltype <- sce_slingshot1$celltype
umap_df$group <- sce_slingshot1$group

# 自定义颜色
cell_colors <- c("CCR2_Mono" = "#e64b35",  
                 "OSM_Mono" = "#4dbbd5",  
                 "HES4_Mac" = "#00a087",
                 "LTB4R_Mac" = "#3e7cbd",
                 "KANK1_Mac" = "#91d1c2",
                 "CXCL10_Mac" = "#b09c85",
                 "GPNMB_Mac" = "#8491b4",
                 "MARCO_Mac" = "#f39b7f")

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
dev.off()  

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

# 第3条轨迹
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$pseudotime <- sce_slingshot1$slingPseudotime_3
curve_obj <- slingCurves(sce_slingshot1)[[3]]  
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

# 画拟时间细胞类型密度图
load("Mac细胞鉴定.Rdata")
load("Mac拟时之Slingshot分析.Rdata")

# 直接提取拟时间值加到元数据列即可
slingPseudotime = as.data.frame(slingPseudotime(sce_slingshot1))
identical(rownames(slingPseudotime),colnames(Mac))
Mac$Lineage1 = slingPseudotime$Lineage1
Mac$Lineage2 = slingPseudotime$Lineage2 
Mac$Lineage3 = slingPseudotime$Lineage3 

# 轨迹1 
merge <- Mac@meta.data[, c("Lineage1", "celltype")]
merge$cell <- rownames(merge)
merge <- merge[!is.na(merge$Lineage1), ]
table(merge$celltype)
merge <- merge[merge$celltype != "OSM_Mono", ]
merge <- merge[merge$celltype != "LTB4R_Mac", ]
table(merge$celltype)
merge <- merge[order(merge$Lineage1), ]
merge$celltype <- factor(merge$celltype, levels = c("CCR2_Mono", "HES4_Mac", "CXCL10_Mac", "GPNMB_Mac", "MARCO_Mac")) # 轨迹1的5个亚群
ggplot(merge, aes(x=Lineage1,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#e64b35",          "#00a087",
                               "#b09c85","#8491b4","#f39b7f")) 

# 轨迹2 
merge <- Mac@meta.data[, c("Lineage2", "celltype")]
merge$cell <- rownames(merge)
merge <- merge[!is.na(merge$Lineage2), ]
table(merge$celltype)
merge <- merge[merge$celltype != "GPNMB_Mac", ]
merge <- merge[merge$celltype != "OSM_Mono", ]
table(merge$celltype)
merge <- merge[order(merge$Lineage2), ]
merge$celltype <- factor(merge$celltype, levels = c("CCR2_Mono","HES4_Mac", "LTB4R_Mac", "KANK1_Mac")) # 轨迹1的5个亚群
ggplot(merge, aes(x=Lineage2,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#e64b35","#00a087","#3e7cbd",
                               "#91d1c2")) 

# 轨迹3 
merge <- Mac@meta.data[, c("Lineage3", "celltype")]
merge$cell <- rownames(merge)
merge <- merge[!is.na(merge$Lineage3), ]
table(merge$celltype)
merge <- merge[merge$celltype != "HES4_Mac", ]
merge <- merge[merge$celltype != "LTB4R_Mac", ]
table(merge$celltype)
merge <- merge[order(merge$Lineage3), ]
merge$celltype <- factor(merge$celltype, levels = c("CCR2_Mono","OSM_Mono")) 
library(ggridges)
ggplot(merge, aes(x=Lineage3,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#e64b35","#4dbbd5")) 


# =============================================
# 🎨 Supplementary Figure 12a
# =============================================

# 读取单细胞对象
load("Mac细胞鉴定.Rdata")

# 修改细胞名为: 样本 + barcodes
colnames(Mac)[1:5] 
meta.data <- Mac@meta.data
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
Mac <- Seurat::RenameCells(Mac, new.names = rownames(meta.data))
colnames(Mac)[1:5]  

# 接下来将Seurat对象转为python中的数据格式 h5ad
library(SCNT)
library(reticulate)
reticulate::use_python("/mnt/data/tool/miniconda3/envs/scanpy/bin/python", required = TRUE)
reticulate::py_config() 
GetH5ad(Mac,
        output_path = "Mac.h5ad",
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
sce=anndata.read_h5ad("Mac.h5ad")
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
scv.pl.proportions(adata, groupby='celltype',layers=['spliced','unspliced','ambiguous'], save='RNA速率7.pdf')
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, 
                                 basis='umap',
                                 color='celltype',
                                 figsize=(5, 5),
                                 save='RNA速率8.pdf',
                                 arrow_size=1.25)
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
            save='RNA速率9.pdf',
            figsize=(6, 6)
            ) 


# =============================================
# 🎨 Supplementary Figure 12d
# =============================================

# 加载数据
load("Mac细胞鉴定.Rdata")

# 加载R包
library(CytoTRACE) 
  
# 提取count矩阵
exp1 <- as.matrix(Mac@assays$RNA$counts)
dim(exp1)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,] 
dim(exp1) 
  
# 运行CytoTRACE函数
results <- CytoTRACE(exp1,                   
                     batch = NULL,          
                     enableFast = TRUE,     
                     ncores = 8,             
                     subsamplesize = 1000   
                     ) 

# 保存
save(results, file = "Mac拟时之CytoTRACE分析.Rdata")   
load("Mac拟时之CytoTRACE分析.Rdata")
class(results)

# 提取出每个细胞的CytoTRACE结果   
CytoTRACE = as.data.frame(results[[1]]) 

# 确定行名和对象元数据行名一致
identical(rownames(CytoTRACE),rownames(Mac@meta.data)) 

# 为CytoTRACE数据框加上细胞类型分组列
CytoTRACE$celltype <- Mac@meta.data$celltype
colnames(CytoTRACE)[1] = "value"
colnames(CytoTRACE)

# 重新设置水平信息
unique(CytoTRACE$celltype)
CytoTRACE$celltype <- factor(CytoTRACE$celltype, 
                             levels = c("CCR2_Mono", "OSM_Mono", "HES4_Mac", "LTB4R_Mac", "KANK1_Mac" , "CXCL10_Mac", "GPNMB_Mac", "MARCO_Mac"))

# 画箱线图
ggplot(CytoTRACE, aes(x = celltype,y = value)) + 
  geom_jitter(col="#00000020", 
              pch=19,         
              cex=0.8,          
              position = position_jitter(0.2)) + 
  geom_boxplot(position=position_dodge(0),  
               aes(color = factor(celltype),fill = factor(celltype)), 
               outlier.shape = NA, 
               alpha = 0.75) +
  scale_fill_manual(values = c("#e64b35","#4dbbd5","#00a087","#3e7cbd",
                               "#91d1c2","#b09c85","#8491b4","#f39b7f"))+
  scale_color_manual(values = c("#9a9a9a", "#9a9a9a", "#9a9a9a", "#9a9a9a",
                                "#9a9a9a", "#9a9a9a", "#9a9a9a", "#9a9a9a"), guide = "none") +
  labs(x=NULL,
       y="Predicted ordering\nby CytoTRACE",
       color = "celltype") +
  coord_cartesian(ylim = c(0, 1.16)) + 
  theme_bw() + 
  theme(text = element_text(family = "sans"),                        
        plot.title = element_text(size = rel(1.16),hjust = 0.5),     
        axis.title = element_text(size = rel(1.16),colour = "black"),
        axis.text = element_text(size=rel(1),colour = "black"),     
        axis.text.x = element_text(angle = 45, hjust = 1),        
        legend.position = "none"                                  
        ) 
dev.off()

# 差异分析
kruskal.test(value ~ celltype, data = CytoTRACE) # p-value < 2.2e-16

