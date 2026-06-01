
# F6/S10: Clinical ketosis is associated with impaired IgG-associated plasma cell responses.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 B/PC细胞亚群数据处理
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
  
# 提取出里面的B/PC细胞
BPC <- subset(seurat_object, celltype %in% c("B", "PC"))
  
# 标准化 
BPC <- NormalizeData(BPC, normalization.method = "LogNormalize", scale.factor = 10000) 

# 高变基因
BPC <- FindVariableFeatures(BPC, selection.method = "vst", nfeatures = 2000) 

# 周期打分
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(homologene) 
X = homologene(s.genes,inTax = 9606,outTax = 9913) 
Y = homologene(g2m.genes,inTax = 9606,outTax = 9913) 
s.genes = X$"9913"
g2m.genes = Y$"9913" 
BPC <- CellCycleScoring(BPC, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) 

# 归一化缩放去除周期影响
BPC <- ScaleData(BPC, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(BPC))

# 线性降维PCA 默认用高变基因集
BPC <- RunPCA(BPC, features = VariableFeatures(object = BPC))

# 肘部图
ElbowPlot(BPC, 50)

# 计算KNN和SNN 
BPC = FindNeighbors(BPC, dims = 1:20)

# 分群数量
BPC = FindClusters(BPC, resolution = c(seq(0.1, 1, 0.1)))
library(clustree)
clustree(BPC, prefix = "RNA_snn_res.")
Idents(BPC) <- "RNA_snn_res.0.4"
BPC$seurat_clusters <- BPC@active.ident

# UMAP非线性降维
BPC <- RunUMAP(BPC, dims = 1:25) 
  
# UMAP图
Seurat::DimPlot(BPC, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(BPC, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

# 保存数据
save(BPC, file = "BPC降维聚类.Rdata") 
rm(list = ls())
load("BPC降维聚类.Rdata")  
  
# 删掉细胞数量很少的7/8
BPC <- subset(BPC, seurat_clusters %in% c("7","8"), invert = TRUE)  
  
# UMAP图
Seurat::DimPlot(BPC, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(BPC, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)  

# 重新走降维聚类流程
BPC <- FindNeighbors(BPC,
                    dims = 1:20 
                    )  

# 批量设置分辨率
BPC <- FindClusters(object = BPC,
                    resolution = c(seq(0.1, 1, 0.1))
                   )
library(clustree)
clustree(BPC, prefix = "RNA_snn_res.") 
Idents(BPC) <- "RNA_snn_res.0.3"
BPC$seurat_clusters <- BPC@active.ident

# UMAP非线性降维
BPC <- RunUMAP(BPC, dims = 1:25) 

# UMAP图
Seurat::DimPlot(BPC, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(BPC, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)

# 保存数据
save(BPC, file = "BPC原始数据.Rdata") 
rm(list = ls())
load("BPC原始数据.Rdata")   
  
# 未注释前细胞cluster进行差异基因分析 用来鉴定细胞
dif<-FindAllMarkers(BPC, 
                    group.by = BPC@meta.data$seurat_clusters, 
                    logfc.threshold = 0.25,                        
                    min.pct = 0.2,                                      
                    only.pos = T                                        
                    )        
dif$pct_diff <- dif$pct.1 - dif$pct.2 
table(dif$cluster)                    
dif<-dif %>%
  group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
head(dif[dif$cluster == unique(dif$cluster)[1],]$gene, 50)   
save(dif, file = "BPC整体cluster差异基因.Rdata")

# 使用BPC整体cluster差异基因进行辅助 来鉴定细胞
load("BPC整体cluster差异基因.Rdata")

# TOP Marker
Seurat::DotPlot(BPC, features = c("BACH2",
                                  "TSPAN2",
                                  "COPZ2",
                                  "CD99", 
                                  "ERN1", 
                                  "IGCGAMMA"), group.by = "seurat_clusters") + RotatedAxis()

# 为分群重新指定细胞类型 
new.cluster.ids <- c("BACH2",
                     "TSPAN2",
                     "COPZ2",
                     "CD99",
                     "ERN1",
                     "IGCGAMMA") 
new.cluster.ids
names(new.cluster.ids) 
levels(BPC)
names(new.cluster.ids) <- levels(BPC) 
names(new.cluster.ids)
new.cluster.ids
BPC <- RenameIdents(BPC, new.cluster.ids) 

BPC[["celltype"]] <- Idents(BPC) 
table(BPC@meta.data[["orig.ident"]])
unique(BPC@meta.data[["orig.ident"]])
BPC[["group"]]<- c(rep("Health", 905), rep("Ketosis", 2763))
BPC[["group.celltype"]]<-paste(BPC$group, Idents(BPC), sep = '_') 
table(BPC@meta.data[["orig.ident"]])
table(BPC@meta.data[["group"]])
table(BPC@meta.data[["celltype"]])
table(BPC@meta.data[["group.celltype"]])
table(BPC@meta.data[["seurat_clusters"]])

# 细胞水平信息
Idents(BPC) <- factor(Idents(BPC),
                       levels = c("BACH2", 
                                  "TSPAN2",
                                  "CD99", 
                                  "COPZ2",
                                  "ERN1",
                                  "IGCGAMMA"))
BPC[["celltype"]] <- Idents(BPC)

# 绘制总umap图
Seurat::DimPlot(BPC, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(BPC, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 保存工作空间
save(BPC,file = "BPC细胞鉴定.Rdata") 
rm(list = ls())
load("BPC细胞鉴定.Rdata")
  
# 对注释后细胞亚群 进行差异基因分析
diff <- FindAllMarkers(BPC, 
                      group.by = BPC@active.ident, 
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
head(diff[diff$cluster == unique(diff$cluster)[1],]$gene, 50)
save(diff, file = "BPC整体细胞类型差异基因.Rdata")    


# =============================================
# 🎨 Fig. 6a
# =============================================

load("BPC细胞鉴定.Rdata")

Seurat::DimPlot(BPC,
                group.by = "celltype",
                cols = c(
                  "#e64b35", 
                  "#4dbbd5", 
                  "#00a087",  
                  "#bc80bd",
                  "#f39b7f",
                  "#7384aa"
                  ), 
                pt.size = 0.5,
                label = F) +
  NoLegend()+ 
  labs(title = NULL)  
dev.off()  


# =============================================
# 🎨 Fig. 6b
# =============================================

plot_density(BPC, 
             reduction = "umap",
             features = c("MS4A1"),
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
  
plot_density(BPC, 
             reduction = "umap",
             features = c("MZB1"),
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

plot_density(BPC, 
             reduction = "umap",
             features = c("BACH2"),
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

plot_density(BPC, 
             reduction = "umap",
             features = c("TSPAN2"),
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

plot_density(BPC, 
             reduction = "umap",
             features = c("CD99"),
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

plot_density(BPC, 
             reduction = "umap",
             features = c("COPZ2"),
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

plot_density(BPC, 
             reduction = "umap",
             features = c("ERN1"),
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

plot_density(BPC, 
             reduction = "umap",
             features = c("IGCGAMMA"),
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
# 🎨 Fig. 6c
# =============================================

# 加载数据
load("BPC细胞鉴定.Rdata")
table(BPC@meta.data[["celltype"]])
  
# 随机抽样 
table(BPC@meta.data[["celltype"]])
cell_types <- unique(BPC$celltype)
subset_cells <- list()  
  
# 对每个细胞类型进行处理
for (cell_type in cell_types) {
  cells_of_type <- WhichCells(BPC, expression = celltype == cell_type) 
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
new_seurat_object <- subset(BPC, cells = subset_all_cells) 
  
# 查看新的 Seurat 对象
new_seurat_object
  
# 基于上调基因分析挑选用于绘图的基因
dif<-FindAllMarkers(BPC,
                    group.by = BPC@active.ident,
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
# 🎨 Fig. 6d
# =============================================

# 首先进行亚群间的差异分析
dif <- FindAllMarkers(BPC, 
                      group.by = BPC@active.ident, 
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

# 单拎出来每个群体的显著基因
table(BPC@active.ident)
BACH2 <- subset(sigposDEG.all, cluster=='BACH2') 
TSPAN2 <- subset(sigposDEG.all, cluster=='TSPAN2')
CD99 <- subset(sigposDEG.all, cluster=='CD99')
COPZ2 <- subset(sigposDEG.all, cluster=='COPZ2')  
ERN1 <- subset(sigposDEG.all, cluster=='ERN1')
IGCGAMMA <- subset(sigposDEG.all, cluster=='IGCGAMMA') 

# 将不同细胞群体的上调基因保存为列表
list <- list(BACH2, TSPAN2, CD99, COPZ2, ERN1, IGCGAMMA)
names(list)[1:6] <- c("BACH2", "TSPAN2", "CD99", "COPZ2", "ERN1", "IGCGAMMA")
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
BACH2 <- read.csv("GO_new_BACH2.CSV", row.names = 1) 
TSPAN2 <- read.csv("GO_new_TSPAN2.CSV", row.names = 1)       
CD99 <- read.csv("GO_new_CD99.CSV", row.names = 1)  
COPZ2 <- read.csv("GO_new_COPZ2.CSV", row.names = 1)   
ERN1 <- read.csv("GO_new_ERN1.CSV", row.names = 1)  
IGCGAMMA <- read.csv("GO_new_IGCGAMMA.CSV", row.names = 1)   

# 为每个细胞群体添加标签
BACH2$group <- "BACH2"
TSPAN2$group <- "TSPAN2"
CD99$group <- "CD99"
COPZ2$group <- "COPZ2"
ERN1$group <- "ERN1"
IGCGAMMA$group <- "IGCGAMMA"

# 选择TOP通路
# BACH2
select_BACH2 = c("regulation of lymphocyte activation",  
                 "regulation of B cell activation", 
                 "MHC class II protein complex assembly",
                 "B cell receptor signaling pathway")
# TSPAN2
select_TSPAN2 = c("actin cytoskeleton organization",  
                  "actin filament organization",
                  "small GTPase-mediated signal transduction",
                  "regulation of cell-cell adhesion")
# CD99
select_CD99 = c("antigen processing and presentation of peptide antigen via MHC class II",
                "immunoglobulin production involved in immunoglobulin-mediated immune response",
                "positive regulation of T cell activation",
                "regulation of protein-containing complex assembly")
# COPZ2
select_COPZ2 = c("Golgi vesicle transport",
                "endoplasmic reticulum to Golgi vesicle-mediated transport",
                "protein targeting to ER",
                "establishment of protein localization to endoplasmic reticulum")
# ERN1
select_ERN1 = c("response to endoplasmic reticulum stress",
                 "endoplasmic reticulum unfolded protein response",
                 "protein folding",
                 "ERAD pathway")
# IGCGAMMA
select_IGCGAMMA = c("protein glycosylation",
                "macromolecule glycosylation",
                "protein maturation",
                "intracellular protein transport")

# 选择每个亚群的通路
BACH2 <- BACH2[BACH2$Description %in% select_BACH2,]
TSPAN2 <- TSPAN2[TSPAN2$Description %in% select_TSPAN2,]
CD99 <- CD99[CD99$Description %in% select_CD99,]
COPZ2 <- COPZ2[COPZ2$Description %in% select_COPZ2,] 
ERN1 <- ERN1[ERN1$Description %in% select_ERN1,]
IGCGAMMA <- IGCGAMMA[IGCGAMMA$Description %in% select_IGCGAMMA,] 

# 生成新的P值列
BACH2$`-log10pvalue` <- -log10(BACH2$pvalue)
TSPAN2$`-log10pvalue` <- -log10(TSPAN2$pvalue)
CD99$`-log10pvalue` <- -log10(CD99$pvalue)
COPZ2$`-log10pvalue` <- -log10(COPZ2$pvalue) 
ERN1$`-log10pvalue` <- -log10(ERN1$pvalue)
IGCGAMMA$`-log10pvalue` <- -log10(IGCGAMMA$pvalue) 

# 合并所有数据
all <- rbind(BACH2, TSPAN2, CD99, COPZ2, ERN1, IGCGAMMA)
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
My_levels <- c("BACH2", "TSPAN2", "CD99", "COPZ2", "ERN1", "IGCGAMMA")
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
  scale_size_continuous(range = c(2, 5))


# =============================================
# 🎨 Fig. 6e, f and Supplementary Figure 10
# =============================================

# 加载数据
load("BPC细胞鉴定.Rdata")

# 在R中输出原始表达矩阵
head(as.matrix(BPC@assays$RNA$counts)[1:3,1:3])
colnames(as.matrix(BPC@assays$RNA$counts))[1:10]
rownames(as.matrix(BPC@assays$RNA$counts))[1:10]
write.csv(as.matrix(BPC@assays$RNA$counts), file = "counts.csv")

# 在JupyterLab(pyscenic)中将原始表达矩阵文件转为loom文件
import os,sys                                                   
os.getcwd()                                                     
os.listdir(os.getcwd())                                        
import loompy as lp                                            
import numpy as np                                              
import pandas as pd                                             
counts = pd.read_csv("counts.csv", index_col = 0, header = 0)  
features = counts.index.tolist()                               
barcodes = counts.columns.tolist()                            
row_attrs = {'Gene':np.array(features)}                       
col_attrs = {'CellID':np.array(barcodes)}                     
lp.create("counts.loom", counts.values, row_attrs, col_attrs)   

# Linux终端pyscenic小环境中运行 三步走(grn/ctx/AUCell) 最终输出的aucell.loom是我们需要的
conda activate pyscenic

# 第一步 grn
nohup pyscenic grn \
  --num_workers 20 \                
  --output grn_output.tsv \         
  --method grnboost2 \              
  counts.loom \                      
  allTFs_hg38.txt > grn.log 2>&1 &  

# 第二步 ctx
nohup pyscenic ctx \
  grn_output.tsv \                                                        
  hg19-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather \  
  --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \        
  --expression_mtx_fname counts.loom \                                     
  --mode "dask_multiprocessing" \                                         
  --output ctx.csv \                                                    
  --num_workers 20 \                                                     
  --mask_dropouts > ctx.log 2>&1 &                                    

# 第三步 AUCell
nohup pyscenic aucell \           
  counts.loom \                              
  ctx.csv \                                 
  --output aucell.loom \                  
  --num_workers 15 > aucell.log 2>&1 &  

# 结果可视化
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(pheatmap)
library(circlize)
library(gridExtra)

# 加载原始Seurat对象
load("BPC细胞鉴定.Rdata")

# 加载结果文件
loom <- open_loom("./aucell.loom")


# get_regulons_AUC返回计算好的每个细胞中转录因子的活性分数矩阵 
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
class(regulonAUC)
regulonAUC 

# 确保活性分数矩阵的细胞列名顺序和原始Seurat对象一致
sub_regulonAUC <- regulonAUC[,match(colnames(BPC),colnames(regulonAUC))]
identical(colnames(sub_regulonAUC), colnames(BPC))

# 随后将每个细胞的aucell值合并添加到seurat对象的meta.data中
BPC@meta.data = cbind(BPC@meta.data,t(sub_regulonAUC@assays@data$AUC))

scplotter::FeatureStatPlot(BPC, 
                           pt_size = 1.5,
                           plot_type = "dim", 
                           features = "IKZF1(+)", 
                           reduction = "umap") +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  

scplotter::FeatureStatPlot(BPC, 
                           pt_size = 1.5,
                           plot_type = "dim", 
                           features = "REL(+)", 
                           reduction = "umap") + 
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  

scplotter::FeatureStatPlot(BPC, 
                           pt_size = 1.5,
                           plot_type = "dim", 
                           features = "IRF8(+)", 
                           reduction = "umap") 
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  

scplotter::FeatureStatPlot(BPC, 
                           pt_size = 1.5,
                           plot_type = "dim", 
                           features = "CREB3L2(+)", 
                           reduction = "umap") + 
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  

scplotter::FeatureStatPlot(BPC, 
                           pt_size = 1.5,
                           plot_type = "dim", 
                           features = "XBP1(+)", 
                           reduction = "umap") + 
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  

scplotter::FeatureStatPlot(BPC, 
                           pt_size = 1.5,
                           plot_type = "dim", 
                           features = "ATF4(+)", 
                           reduction = "umap") + 
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  

# 热图展示每个细胞类型中top转录因子的aucell均值
auc <- getAUC(sub_regulonAUC)
auc_avg <- data.frame(t(auc))
auc_avg[,"group"] <- BPC@meta.data$celltype
auc_avg <- auc_avg %>%
  group_by(group) %>%
  summarise_all(mean)
auc_avg <- as.data.frame(auc_avg)
rownames(auc_avg ) <- auc_avg$group
colnames(auc_avg) <- c("group",rownames(auc)) 
top3_reg <- colnames(auc_avg)[2:60]
auc_avg_sub <- auc_avg[,top3_reg]
rownames(auc_avg) <- factor(rownames(auc_avg) , levels = c("BACH2","TSPAN2","CD99", "COPZ2", "ERN1", "IGCGAMMA"))

# 热图展示
pheatmap::pheatmap(t(auc_avg_sub), 
                   show_colnames = T, 
                   cluster_rows = T, 
                   cluster_cols = F, 
                   scale = "row" 
                   )

# RSS
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = BPC@meta.data$celltype)
rssPlot <- plotRSS(rss)
rssPlot$plot

BPC$celltype1 <- ifelse(
  BPC$celltype %in% c("BACH2", "TSPAN2", "CD99"),
  "B",
  ifelse(
    BPC$celltype %in% c("COPZ2", "ERN1", "IGCGAMMA"),
    "PC",
    NA
  )
)  

Idents(BPC) <- "celltype1"
Idents(BPC) <- "celltype"

auc <- getAUC(sub_regulonAUC)
auc_avg <- data.frame(t(auc))
auc_avg[,"group"] <- BPC@meta.data$celltype
auc_avg <- auc_avg %>%
  group_by(group) %>%
  summarise_all(mean)
auc_avg <- as.data.frame(auc_avg)
rownames(auc_avg ) <- auc_avg$group
colnames(auc_avg) <- c("group",rownames(auc)) 
top3_reg <- c("IKZF1(+)", "REL(+)", "IRF8(+)",
              "CREB3L2(+)", "XBP1(+)", "ATF4(+)")
auc_avg_sub <- auc_avg[,top3_reg]
rownames(auc_avg) <- factor(rownames(auc_avg), levels = c("BACH2","TSPAN2","CD99", "COPZ2", "ERN1", "IGCGAMMA"))

library(pheatmap)
pheatmap(t(auc_avg_sub), 
         show_colnames = T, 
         cluster_rows = F, 
         cluster_cols = F, 
         colorRampPalette(rev(brewer.pal(11, "RdBu")))(11),
         scale = "row" 
         )

# RSS
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = BPC@meta.data$celltype1)
rssPlot <- plotRSS(rss)
rssPlot$plot

# 分细胞类型展示RSS特异性分数
plotRSS_oneSet(rss,'B') 
plotRSS_oneSet(rss,'PC')

# 可视化TF调控的靶基因
sce_regulons <- read.csv("ctx.csv")
sce_regulons <- sce_regulons[-2, ] 
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1, ] 
colnames(sce_regulons) <- c("TF","ID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                            "Annotation","Context","TargetGenes","RankAtMax")

# 我这里关注IKZF1/IRF8/REL这3个TF
# IKZF1
IKZF1 <- subset(sce_regulons, TF=='IKZF1')
IKZF1 <- IKZF1[which(IKZF1$AUC>0.1),]
IKZF1 <- IKZF1[, c("TF","TargetGenes")]
IKZF1$TargetGenes <-gsub("\\[","",IKZF1$TargetGenes)
IKZF1$TargetGenes <-gsub("\\]","",IKZF1$TargetGenes)
IKZF1$TargetGenes <-gsub("\\(","",IKZF1$TargetGenes)
IKZF1$TargetGenes <-gsub("\\)","",IKZF1$TargetGenes)
IKZF1$TargetGenes <-gsub("\\'","",IKZF1$TargetGenes)
library(stringr)
split_IKZF1<-str_split(IKZF1$TargetGenes,",")
IKZF11 <- as.data.frame(split_IKZF1[[1]])
IKZF12<- as.data.frame(split_IKZF1[[2]])
IKZF13<-as.data.frame(split_IKZF1[[3]])

names(IKZF11) <- 'TF'
names(IKZF12) <- 'TF'
names(IKZF13) <- 'TF'

IKZF1 <- rbind(IKZF11,IKZF12,IKZF13)

IKZF1_target <- IKZF1[seq(1,nrow(IKZF1),2), ]
IKZF1_score <- IKZF1[seq(0,nrow(IKZF1),2), ]

IKZF1_gene <- data.frame(IKZF1_target,IKZF1_score)
IKZF1_gene <- IKZF1_gene[!duplicated(IKZF1_gene$IKZF1_target), ]
IKZF1_gene$gene <- 'IKZF1'
colnames(IKZF1_gene) <- c("target","score",'tf')

# 同理得到REL及其靶基因
REL <- subset(sce_regulons, TF=='REL')
REL <- REL[which(REL$AUC>0.1),]
REL <- REL[, c("TF","TargetGenes")]
REL$TargetGenes <-gsub("\\[","",REL$TargetGenes)
REL$TargetGenes <-gsub("\\]","",REL$TargetGenes)
REL$TargetGenes <-gsub("\\(","",REL$TargetGenes)
REL$TargetGenes <-gsub("\\)","",REL$TargetGenes)
REL$TargetGenes <-gsub("\\'","",REL$TargetGenes)
library(stringr)
split_REL<-str_split(REL$TargetGenes,",")
REL1 <- as.data.frame(split_REL[[1]])
names(REL1) <- 'TF'
REL <- rbind(REL1)
REL_target <- REL[seq(1,nrow(REL),2), ]
REL_score <- REL[seq(0,nrow(REL),2), ]
REL_gene <- data.frame(REL_target,REL_score)
REL_gene <- REL_gene[!duplicated(REL_gene$REL_target), ]
REL_gene$gene <- 'REL'
colnames(REL_gene) <- c("target","score",'tf')

# 同理得到IRF8及其靶基因
IRF8 <- subset(sce_regulons, TF=='IRF8')
IRF8 <- IRF8[which(IRF8$AUC>0.065),]
IRF8 <- IRF8[, c("TF","TargetGenes")]
IRF8$TargetGenes <-gsub("\\[","",IRF8$TargetGenes)
IRF8$TargetGenes <-gsub("\\]","",IRF8$TargetGenes)
IRF8$TargetGenes <-gsub("\\(","",IRF8$TargetGenes)
IRF8$TargetGenes <-gsub("\\)","",IRF8$TargetGenes)
IRF8$TargetGenes <-gsub("\\'","",IRF8$TargetGenes)
library(stringr)
split_IRF8<-str_split(IRF8$TargetGenes,",")
IRF81 <- as.data.frame(split_IRF8[[1]])
names(IRF81) <- 'TF'
IRF8 <- rbind(IRF81)
IRF8_target <- IRF8[seq(1,nrow(IRF8),2), ]
IRF8_score <- IRF8[seq(0,nrow(IRF8),2), ]
IRF8_gene <- data.frame(IRF8_target,IRF8_score)
IRF8_gene <- IRF8_gene[!duplicated(IRF8_gene$IRF8_target), ]
IRF8_gene$gene <- 'IRF8'
colnames(IRF8_gene) <- c("target","score",'tf')

# two thousand years later
# 合并3个转录因子
TF_target <- rbind(IKZF1_gene, IRF8_gene, REL_gene)
TF_target$score <- as.numeric(TF_target$score)

# 整理 TF-target 网络数据
library(dplyr)
TF_target <- TF_target %>%
  filter(!is.na(score)) %>%
  mutate(
    tf = trimws(tf),
    target = trimws(target)
  )

TF_target_top <- TF_target %>%
  group_by(tf) %>%
  group_modify(~{
    if (.y$tf == "IRF8") {
      slice_max(.x, order_by = score, n = 20, with_ties = FALSE)
    } else {
      .x 
    }
  }) %>%
  ungroup()
table(TF_target_top$tf)

# 构建nodes和edges
edges <- TF_target_top %>%
  select(from = tf, to = target, score)
nodes <- data.frame(
  name = unique(c(edges$from, edges$to))
)
nodes$type <- ifelse(nodes$name %in% c("IKZF1", "IRF8", "REL"), "TF", "Target")

# 画网络图
library(igraph)
library(ggraph)
library(ggplot2)

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

set.seed(128)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = score), alpha = 0.3) +
  geom_node_point(aes(size = type, color = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_size_manual(values = c("TF" = 8, "Target" = 4)) +
  scale_color_manual(values = c("TF" = "#E64B35", "Target" = "#4DBBD5")) +  # 👈 关键
  theme_void()

# 我们还可以可视化感兴趣的TF调控的靶基因
sce_regulons <- read.csv("ctx.csv")
sce_regulons <- sce_regulons[-2, ] 
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1, ] 
colnames(sce_regulons) <- c("TF","ID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                            "Annotation","Context","TargetGenes","RankAtMax")

# 我这里关注XBP1/CREB3L2/ATF4这3个TF
# XBP1
XBP1 <- subset(sce_regulons, TF=='XBP1')
XBP1 <- XBP1[which(XBP1$AUC>0.07),]
XBP1 <- XBP1[, c("TF","TargetGenes")]
XBP1$TargetGenes <-gsub("\\[","",XBP1$TargetGenes)
XBP1$TargetGenes <-gsub("\\]","",XBP1$TargetGenes)
XBP1$TargetGenes <-gsub("\\(","",XBP1$TargetGenes)
XBP1$TargetGenes <-gsub("\\)","",XBP1$TargetGenes)
XBP1$TargetGenes <-gsub("\\'","",XBP1$TargetGenes)
library(stringr)
split_XBP1<-str_split(XBP1$TargetGenes,",")
XBP11 <- as.data.frame(split_XBP1[[1]])
XBP12 <- as.data.frame(split_XBP1[[2]])
XBP13 <- as.data.frame(split_XBP1[[3]])
XBP14 <- as.data.frame(split_XBP1[[4]])
XBP15 <- as.data.frame(split_XBP1[[5]])
XBP16 <- as.data.frame(split_XBP1[[6]])
XBP17 <- as.data.frame(split_XBP1[[7]])
names(XBP11) <- 'TF'
names(XBP12) <- 'TF'
names(XBP13) <- 'TF'
names(XBP14) <- 'TF'
names(XBP15) <- 'TF'
names(XBP16) <- 'TF'
names(XBP17) <- 'TF'
XBP1 <- rbind(XBP11, XBP12, XBP13)

XBP1_target <- XBP1[seq(1,nrow(XBP1),2), ]
XBP1_score <- XBP1[seq(0,nrow(XBP1),2), ]

XBP1_gene <- data.frame(XBP1_target,XBP1_score)
XBP1_gene <- XBP1_gene[!duplicated(XBP1_gene$XBP1_target), ]
XBP1_gene$gene <- 'XBP1'
colnames(XBP1_gene) <- c("target","score",'tf')

# 同理得到REL及其靶基因
CREB3L2 <- subset(sce_regulons, TF=='CREB3L2')
CREB3L2 <- CREB3L2[which(CREB3L2$AUC>0.1),]
CREB3L2 <- CREB3L2[, c("TF","TargetGenes")]
CREB3L2$TargetGenes <-gsub("\\[","",CREB3L2$TargetGenes)
CREB3L2$TargetGenes <-gsub("\\]","",CREB3L2$TargetGenes)
CREB3L2$TargetGenes <-gsub("\\(","",CREB3L2$TargetGenes)
CREB3L2$TargetGenes <-gsub("\\)","",CREB3L2$TargetGenes)
CREB3L2$TargetGenes <-gsub("\\'","",CREB3L2$TargetGenes)
library(stringr)
split_CREB3L2<-str_split(CREB3L2$TargetGenes,",")
CREB3L21 <- as.data.frame(split_CREB3L2[[1]])
CREB3L22 <- as.data.frame(split_CREB3L2[[2]])
CREB3L23 <- as.data.frame(split_CREB3L2[[3]])
CREB3L24 <- as.data.frame(split_CREB3L2[[4]])
CREB3L25 <- as.data.frame(split_CREB3L2[[5]])
CREB3L26 <- as.data.frame(split_CREB3L2[[6]])
names(CREB3L21) <- 'TF'
names(CREB3L22) <- 'TF'
names(CREB3L23) <- 'TF'
names(CREB3L24) <- 'TF'
names(CREB3L25) <- 'TF'
names(CREB3L26) <- 'TF'
CREB3L2 <- rbind(CREB3L21, CREB3L22, CREB3L23,CREB3L24,CREB3L25,CREB3L26)
CREB3L2_target <- CREB3L2[seq(1,nrow(CREB3L2),2), ]
CREB3L2_score <- CREB3L2[seq(0,nrow(CREB3L2),2), ]
CREB3L2_gene <- data.frame(CREB3L2_target,CREB3L2_score)
CREB3L2_gene <- CREB3L2_gene[!duplicated(CREB3L2_gene$CREB3L2_target), ]
CREB3L2_gene$gene <- 'CREB3L2'
colnames(CREB3L2_gene) <- c("target","score",'tf')

# 同理得到IRF8及其靶基因
ATF4 <- subset(sce_regulons, TF=='ATF4')
ATF4 <- ATF4[which(ATF4$AUC>0.05),]
ATF4 <- ATF4[, c("TF","TargetGenes")]
ATF4$TargetGenes <-gsub("\\[","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\]","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\(","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\)","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\'","",ATF4$TargetGenes)
library(stringr)
split_ATF4<-str_split(ATF4$TargetGenes,",")
ATF41 <- as.data.frame(split_ATF4[[1]])
names(ATF41) <- 'TF'
ATF4 <- rbind(ATF41)
ATF4_target <- ATF4[seq(1,nrow(ATF4),2), ]
ATF4_score <- ATF4[seq(0,nrow(ATF4),2), ]
ATF4_gene <- data.frame(ATF4_target,ATF4_score)
ATF4_gene <- ATF4_gene[!duplicated(ATF4_gene$ATF4_target), ]
ATF4_gene$gene <- 'ATF4'
colnames(ATF4_gene) <- c("target","score",'tf')

# two thousand years later
# 合并3个转录因子
TF_target <- rbind(XBP1_gene, CREB3L2_gene, ATF4_gene)
TF_target$score <- as.numeric(TF_target$score)

# 整理 TF-target 网络数据
library(dplyr)
TF_target <- TF_target %>%
  filter(!is.na(score)) %>%
  mutate(
    tf = trimws(tf),
    target = trimws(target)
  )

TF_target_top <- TF_target %>%
  group_by(tf) %>%
  slice_max(order_by = score, n = 20, with_ties = FALSE) %>%
  ungroup()
table(TF_target_top$tf)

# 构建nodes和edges
edges <- TF_target_top %>%
  select(from = tf, to = target, score)
nodes <- data.frame(
  name = unique(c(edges$from, edges$to))
)
nodes$type <- ifelse(nodes$name %in% c("XBP1", "CREB3L2", "ATF4"), "TF", "Target")

# 画网络图
library(igraph)
library(ggraph)
library(ggplot2)

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

set.seed(128)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = score), alpha = 0.3) +
  geom_node_point(aes(size = type, color = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_size_manual(values = c("TF" = 8, "Target" = 4)) +
  scale_color_manual(values = c("TF" = "#e64b35", "Target" = "#4dbbd5")) + 
  scale_edge_width(range = c(0.5, 3)) +
  theme_void()


# =============================================
# 🎨 Fig. 6g
# =============================================

# 分组的细胞比例图
load("BPC细胞鉴定.Rdata")  
library(dplyr)
library(ggplot2)
library(gtools)
library(ggalluvial)

# 准备细胞比例输入数据
prop_df <- BPC@meta.data %>%
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

# 画图
ggplot(as.data.frame(prop_df), aes(x = Sample, y = Proportion, fill = Celltype, stratum = Celltype, alluvium = Celltype)) +
  geom_flow(width = 0.6, alpha = 0.3, knot.pos = 0.1) +  
  geom_col(width = 0.6) +  
  scale_y_continuous(expand = c(0, 0)) +   
  scale_fill_manual(values = c(
                  "#e64b35", 
                  "#4dbbd5", 
                  "#00a087",  
                  "#bc80bd",
                  "#f39b7f",
                  "#7384aa"
                  )) +  
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
# 🎨 Fig. 6i
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
load("BPC细胞鉴定.Rdata")

# 将Seurat对象转为SingleCellExperiment对象
sce <- as.SingleCellExperiment(BPC, assay = "RNA") 

# 运行主函数
sce_slingshot1 <- slingshot(data = sce,                    
                            reducedDim = 'UMAP',           
                            clusterLabels = sce$celltype,   
                            start.clus = NULL,              
                            end.clus = NULL,             
                            approx_points = 150             
                            )

# 保存
save(sce_slingshot1, file = "BPC拟时之Slingshot分析.Rdata")  
load("BPC拟时之Slingshot分析.Rdata")

# SlingshotDataSet函数查看轨迹信息
SlingshotDataSet(sce_slingshot1)
dim(slingPseudotime(sce_slingshot1)) 
slingPseudotime(sce_slingshot1)[1:3,1:2] 

# 以细胞类型为基础进行可视化
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$celltype <- sce_slingshot1$celltype
umap_df$group <- sce_slingshot1$group
cell_colors <- c("BACH2" = "#e64b35",  
                 "TSPAN2" = "#4dbbd5",  
                 "CD99" = "#00a087",
                 "COPZ2" = "#bc80bd",
                 "ERN1" = "#f39b7f",
                 "IGCGAMMA" = "#7384aa")

# 提取所有 Slingshot 拟合曲线
curves <- slingCurves(sce_slingshot1)

# 合并所有轨迹
curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  crv <- as.data.frame(curves[[i]]$s)
  crv$lineage <- paste0("Lineage", i)
  colnames(crv)[1:2] <- c("UMAP_1", "UMAP_2")
  return(crv)
}))

# 画总的UMAP轨迹图
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.3, alpha = 0.8, show.legend = FALSE) + 
  scale_color_manual(values = cell_colors) +
  geom_path(data = curve_df, 
            aes(x = UMAP_1, y = UMAP_2, group = lineage), 
            inherit.aes = FALSE,
            color = "black", size = 1) +
  theme_classic2() + 
  labs(x = "UMAP_1", y = "UMAP_2", color = "Cell Type") +
  theme(legend.position = "right")
dev.off()  

# 画第1条轨迹
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$pseudotime <- sce_slingshot1$slingPseudotime_1
curve_obj <- slingCurves(sce_slingshot1)[[1]]  
curve_df <- as.data.frame(curve_obj$s)
colnames(curve_df) <- c("UMAP_1", "UMAP_2")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1, alpha = 0.8, show.legend = FALSE) +
  scale_color_viridis_c(option = "D", na.value = "grey80") +
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
            inherit.aes = FALSE, color = "black", size = 1) +
  theme_void() + 
  labs(x = "UMAP_1", y = "UMAP_2", color = "Pseudotime")
dev.off()

# 画第2条轨迹
umap_df <- as.data.frame(reducedDims(sce_slingshot1)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$pseudotime <- sce_slingshot1$slingPseudotime_2
curve_obj <- slingCurves(sce_slingshot1)[[2]]  
curve_df <- as.data.frame(curve_obj$s)
colnames(curve_df) <- c("UMAP_1", "UMAP_2")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1, alpha = 0.8, show.legend = FALSE) +
  scale_color_viridis_c(option = "D", na.value = "grey80") +
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
            inherit.aes = FALSE, color = "black", size = 1) +
  theme_void() + 
  labs(x = "UMAP_1", y = "UMAP_2", color = "Pseudotime")
dev.off()

# 用scplotter包绘制试试
identical(colnames(sce_slingshot1),colnames(BPC))

# 直接提取拟时间值加到元数据列即可
slingPseudotime = as.data.frame(slingPseudotime(sce_slingshot1))
identical(rownames(slingPseudotime),colnames(BPC))
BPC$Lineage1 = slingPseudotime$Lineage1
BPC$Lineage2 = slingPseudotime$Lineage2

# 保存
save(BPC, file = "BPC细胞鉴定.Rdata")

# 开始用SCP绘制
scplotter::CellDimPlot(BPC, 
                       pt_size = 1,
                       pt_alpha = 0.5,
                       group_by = "celltype", 
                       reduction = "umap",
                       lineages = paste0("Lineage", 1:2),
                       lineages_trim = c(0, 1),
                       lineages_span = 0.75,
                       lineages_palcolor = c("red", "blue"),
                       lineages_arrow = arrow(length = unit(0.1, "inches")),
                       lineages_linewidth = 1.2,
                       lineages_line_bg = "white",
                       lineages_line_bg_stroke = 0.5,
                       lineages_whiskers = T,
                       lineages_whiskers_linewidth = 0.1,
                       lineages_whiskers_alpha = 0.1,
                       palcolor = c("#e64b35", 
                                    "#4dbbd5", 
                                    "#00a087",  
                                    "#bc80bd",
                                    "#f39b7f",
                                    "#7384aa") 
                       )
  

# =============================================
# 🎨 Fig. 6j
# =============================================

# 画拟时间细胞类型密度图
load("BPC细胞鉴定.Rdata")
load("BPC拟时之Slingshot分析.Rdata")
slingPseudotime = as.data.frame(slingPseudotime(sce_slingshot1))
identical(rownames(slingPseudotime),colnames(BPC))
BPC$Lineage1 = slingPseudotime$Lineage1
BPC$Lineage2 = slingPseudotime$Lineage2 
merge <- BPC@meta.data[, c("Lineage1", "celltype")]
merge$cell <- rownames(merge)
merge <- merge[!is.na(merge$Lineage1), ]
merge <- merge[merge$celltype != "ERN1", ]
merge <- merge[order(merge$Lineage1), ]

# 可视化密度图
ggplot(merge, aes(x=Lineage1,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#4dbbd5", 
                               "#e64b35",
                               "#00a087",  
                               "#bc80bd",
                               "#7384aa")) 

merge <- BPC@meta.data[, c("Lineage2", "celltype")]
merge$cell <- rownames(merge)
merge <- merge[!is.na(merge$Lineage2), ]
table(merge$celltype)
merge <- merge[merge$celltype != "IGCGAMMA", ]
merge <- merge[order(merge$Lineage2), ]
ggplot(merge, aes(x=Lineage2,y=celltype,fill=celltype)) +
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right') +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8)) +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = c("#4dbbd5", 
                               "#e64b35",
                               "#00a087",  
                               "#bc80bd",
                               "#f39b7f")) 


# =============================================
# 🎨 Fig. 6h
# =============================================
  
# 加载数据
load("BPC细胞鉴定.Rdata")
BPC$celltype1 <- as.character(BPC$celltype)
BPC$celltype1[BPC$celltype1 %in% c("COPZ2", "ERN1", "IGCGAMMA")] <- "PC"
Idents(BPC) <- BPC$celltype1
  
library(CytoTRACE) 
  
# 提取count矩阵
exp1 <- as.matrix(BPC@assays$RNA$counts)
dim(exp1)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
dim(exp1) 
  
# 运行CytoTRACE函数 
results <- CytoTRACE(exp1,                 
                     batch = NULL,          
                     enableFast = TRUE,     
                     ncores = 1,            
                     subsamplesize = 1000  
                     ) 

# 保存
save(results, file = "BPC拟时之CytoTRACE分析.Rdata")   
load("BPC拟时之CytoTRACE分析.Rdata")
  
# 提取细胞类型信息 并加上细胞名字 用于可视化
phenot <- BPC$celltype1
phenot <- as.character(phenot)
names(phenot) <- rownames(BPC@meta.data)
phenot[1:5]

# 提取UMAP坐标信息 用于可视化
emb <- BPC@reductions[["umap"]]@cell.embeddings
head(emb) 
  
# 自己画图
class(results)
CytoTRACE = as.data.frame(results[[1]]) 
identical(rownames(CytoTRACE),rownames(BPC@meta.data)) 
CytoTRACE$celltype <- BPC@meta.data$celltype1
colnames(CytoTRACE)[1] = "value"
colnames(CytoTRACE)
unique(CytoTRACE$celltype)
CytoTRACE$celltype <- factor(CytoTRACE$celltype, 
                             levels = c("TSPAN2", "BACH2", "CD99", "PC"))

# 箱线图
ggplot(CytoTRACE, aes(x = celltype,y = value)) + 
  geom_jitter(col="#00000020", 
              pch=19,           
              cex=0.8,           
              position = position_jitter(0.2)) + 
  geom_boxplot(position=position_dodge(0),  
               aes(color = factor(celltype),fill = factor(celltype)),
               outlier.shape = NA, 
               alpha = 0.75) + 
  scale_fill_manual(values = c("#4dbbd5",
                               "#e64b35", 
                               "#00a087",  
                               "#fa9fb5"))+
  scale_color_manual(values = c("#9a9a9a", "#9a9a9a", "#9a9a9a", "#9a9a9a"), guide = "none") +
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

# 差异分析
kruskal.test(value ~ celltype, data = CytoTRACE)$p.value 

