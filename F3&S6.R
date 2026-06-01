
# F3/S6: Clinical ketosis induces broad transcriptional and functional remodeling of hepatic immune cells.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 Fig. 3a
# =============================================

# 加载数据
load("CD45去双细胞后.Rdata")    
  
# 为细胞添加新的标签同时区分样本和细胞亚群
seurat_object[["sample.celltype"]]<-paste(seurat_object$orig.ident,
                                          Idents(seurat_object),
                                          sep = '_')  

# 使用AggregateExpression生成伪bulk表达矩阵
pseudo_bulk <- AggregateExpression(
  seurat_object,
  group.by = "orig.ident", 
  assays = "RNA",         
  slot = "counts"
  )  

# 提取UMI表达矩阵 
count_data <- pseudo_bulk$RNA 
count_data = as.data.frame(count_data)

# 数据过滤 主要针对低表达基因
keep_feature <- rowSums (count_data > 0) > 0.50*ncol(count_data) # 在至少50%的样本中表达的基因才会被保留
table(keep_feature) 
keep_sample <- colSums(count_data) / 1000000 >10 # 只保留总reads数超过10M的样本
table(keep_sample)
count_data <- count_data[keep_feature, keep_sample]
dim(count_data)
count_data[1:4, ]

# 转为CPM
dat = edgeR::cpm(count_data)
colSums(dat) 

# 添加分组信息
Group <- c(rep("Health", 3), rep("Ketosis", 5))
Group <- factor(Group, levels = c("Health", "Ketosis"))
Group

# 存储整理好的read counts表达矩阵 用于差异基因分析
save(count_data, Group, file = "Bulk原始矩阵.rdata")

# 存储整理好的TPM标准化表达矩阵 用于PCA/热图/基因比较
save(dat, Group, file = 'Bulk标准化矩阵.Rdata') 

rm(list = ls())
load("Bulk原始矩阵.rdata")      
load("Bulk标准化矩阵.Rdata")   

# PCA
library(FactoMineR) 
library(factoextra) 
boxplot(dat,las = 2)
dat = log2(dat + 1) 
boxplot(dat,las = 2)
dat = as.data.frame(t(dat))      
dat.pca = PCA(dat, graph = FALSE)
pca = prcomp(dat, scale. = TRUE)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
fviz_pca_ind(dat.pca,
             geom.ind = "point",                        
             pointsize = 3,                            
             pointshape = 21,                           
             fill = Group,                             
             col.ind = "gray", 
             addEllipses = TRUE,                      
             ellipse.level = 0.8,                     
             legend.title = "Group"                
             ) + 
  scale_fill_manual(values = c("Health" = "#4cdafe", 
                               "Ketosis" = "#ff6362")) +
  scale_color_manual(values = c("Health" = "#4cdafe", 
                                "Ketosis" = "#ff6362")) +  
  labs(title = "", x = "PC1 (32.10%)", y = "PC2 (23.96%)", color = "Group") +  
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        panel.grid.major = element_blank(),    
        panel.grid.minor = element_blank(),   
        axis.line = element_blank(),          
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.position = "none"
        )


# =============================================
# 🎨 Fig. 3b
# =============================================

# 进行差异基因分析
load("./Bulk原始矩阵.rdata")  

# 构建colData
levels(Group)
colData <- data.frame(
  row.names = colnames(count_data), 
  condition = Group 
  )
colData

# 构建DESeq2的DESeqDataSet对象
dds =  DESeqDataSetFromMatrix(
  countData = count_data,  
  colData = colData,  
  design = ~ condition) 
class(dds)

# 进行差异表达分析
dds2 <- DESeq(dds)
resultsNames(dds2)
tmp <- results(dds2)
head(tmp)
dim(tmp)  

# 按照FC值排序基因
DEG_deseq2 <- as.data.frame(tmp[order(tmp$log2FoldChange, decreasing = TRUE), ])
dim(DEG_deseq2)
table(is.na(DEG_deseq2))
DEG_deseq2 <- DEG_deseq2 %>% drop_na()
dim(DEG_deseq2)
head(DEG_deseq2)  

# 添加上调和下调基因列
logFC_cutoff = log2(1.5)  
pvalue = 0.05       
k1 = (DEG_deseq2$pvalue < pvalue) & (DEG_deseq2$log2FoldChange < -logFC_cutoff)
k2 = (DEG_deseq2$pvalue < pvalue) & (DEG_deseq2$log2FoldChange > logFC_cutoff)
DEG_deseq2$change = ifelse(k1, 'Down', ifelse(k2, 'Up', 'Stable'))
table(DEG_deseq2$change)
table(is.na(DEG_deseq2$change)) # 上调1523 下调1324

# 保存数据
save(DEG_deseq2, file = 'DEG_deseq2.Rdata') 
rm(list=ls()) 
load("DEG_deseq2.Rdata") 

# 画火山图
logFC_t = log2(1.5) 
p_t = 0.05 
head(DEG_deseq2)
df = DEG_deseq2
df$'-log10(P-value)' <- -log10(df$pvalue) 
df$change <- factor(df$change, levels = c("Up","Down","Stable"))

# 要标记的TOP基因
genes <- c("CSF3","SAA3","SELE","CCL17",
           "FCGR2A","MAPK12","CCL16","CCL2",
           "MS4A7","FCGR1A","MX1",
           "APOE","PLA2G2D4","TXNDC5","MARCO","VCAM1",
           "CD68","LYZ2","CTSZ","NUPR1","PLVAP",
           "CSF2RB","C3AR1","IL6","HMOX1","ELANE",
           "AIF1","IFI6","LPL","TREM2","XCR1","IL18","CSF1R","FABP5", # 上调标记
           
           "CTLA4", "FOXP3", "TNFSF4", "NCR1",  
           "RORA", "KLRD1", "S1PR5", "FCAR", "CAMK4",
           "TCF7", "GNLY", "THEMIS", "CD200", "LEF1",
           "EOMES", "BCL11B", "SH2D1A", "TXK", "RASGRP1", "IL7R",
           "IDO1", "TOX", "TRAT1", "IL17RB","IL2RA"# 下调标记
           )
gene_plot <- df[genes, ]

# 把超过阈值的点就设置为阈值 让其可以在图的框线上显示
xlim <- 4     
ylim <- 15    
df$log2FoldChange_capped <- pmax(pmin(df$log2FoldChange, xlim), -xlim)
df$`-log10P_capped` <- pmin(df$`-log10(P-value)`, ylim)

# 画图
ggplot(data = df, aes(x = log2FoldChange_capped, y = `-log10P_capped`, color = change)) + 
  geom_point(size = 2.2) +
  scale_x_continuous(expand = expansion(add = c(0, 0)), limits = c(-4, 4), breaks = seq(-4, 4, by = 2))+ 
  scale_y_continuous(expand = expansion(add = c(0, 0)), limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  scale_colour_manual(name = "", values = alpha(c("#ef4040","#608fe4","#d8d8d8"), 0.4)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_rect(color = "black", size = 1),  
        ) + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 14),
        legend.position = "none"
        ) + 
  geom_hline(yintercept = c(-log10(0.05)),size = 0.8,color = "gray",lty = "dashed") +          
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),size = 0.8,color = "gray",lty = "dashed") + 
  geom_point(data = df[rownames(df) %in% rownames(gene_plot),], stroke = 0.5, size = 3, shape = 16, color = "#4dec65") + 
  geom_text_repel(data = gene_plot, 
                  aes(label = rownames(gene_plot)), 
                  color="black", 
                  size=3.5, 
                  fontface="italic", 
                  arrow = arrow(ends="first", length = unit(0.01, "npc")), 
                  box.padding = 0.2,
                  point.padding = 0.3, 
                  segment.color = 'black', 
                  segment.size = 0.3, 
                  force = 1, 
                  max.iter = 3e3,
                  max.overlaps = Inf)

# 展示火山图的上调下调 上调1523 下调1324
data <- data.frame(
  category = c("A", "B"),
  value = c(1523, 1324)
  )
data$label <- data$value

ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = c("#f8b2b2", "#bfd2f4")) + 
  theme_void() + 
  theme(legend.position = "none")  


# =============================================
# 🎨 Fig. 3c
# =============================================

# 接下来进行功能分析
load("DEG_deseq2.Rdata")

# 做GSEA富集分析
CD45GSEAKEGG <- SCP::RunGSEA(
  geneID = rownames(DEG_deseq2),        
  geneScore = DEG_deseq2$log2FoldChange, 
  geneID_groups = rep("a",times=16330),  
                                        
  db = c("KEGG"),
  species = "Bos_taurus"
  )

CD45GSEAGO <- SCP::RunGSEA(
  geneID = rownames(DEG_deseq2),         
  geneScore = DEG_deseq2$log2FoldChange, 
  geneID_groups = rep("a",times=16330),  
                                       
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
  
# 保存数据
CD45GSEAKEGG = CD45GSEAKEGG[["enrichment"]]
save(CD45GSEAKEGG, file = 'CD45GSEAKEGG.Rdata') 
CD45GSEAGO = CD45GSEAGO[["enrichment"]]
save(CD45GSEAGO, file = 'CD45GSEAGO.Rdata')  

# 加载富集分析结果
load("~/奶牛肝脏解离单细胞3版/CD45GSEAGO.Rdata")
load("~/奶牛肝脏解离单细胞3版/CD45GSEAKEGG.Rdata")
CD45GSEAGO <- CD45GSEAGO[order(CD45GSEAGO$NES, decreasing = TRUE), ]
CD45GSEAKEGG <- CD45GSEAKEGG[order(CD45GSEAKEGG$NES, decreasing = TRUE), ]

# GO上调的TOP通路
CD45GSEAGO$Description[1:100] # antigen processing and presentation of peptide antigen via MHC class II
                              # positive regulation of leukocyte activation
                              # chemokine-mediated signaling pathway
                              # leukocyte chemotaxis
                              # leukocyte migration
                              # inflammatory response
                              # response to tumor necrosis factor
                              # response to type II interferon
                              # leukocyte cell-cell adhesion
                              # phagocytosis
                              # endocytosis
                              # B cell activation
                              # immunoglobulin production
                              # lipid catabolic process
                              # unsaturated fatty acid metabolic process
                              # icosanoid metabolic process

CD45GSEAGO = CD45GSEAGO[CD45GSEAGO$Description %in% 
                                    c("antigen processing and presentation of peptide antigen via MHC class II", 
                                      "positive regulation of leukocyte activation", 
                                      "chemokine-mediated signaling pathway",
                                      "leukocyte chemotaxis",
                                      "leukocyte migration",
                                      "inflammatory response", 
                                      "response to tumor necrosis factor", 
                                      "response to type II interferon", 
                                      "leukocyte cell-cell adhesion", 
                                      "phagocytosis", 
                                      "endocytosis", 
                                      "B cell activation",
                                      "immunoglobulin production",
                                      "lipid catabolic process", 
                                      "unsaturated fatty acid metabolic process",
                                      "icosanoid metabolic process"),] 
CD45GSEAGO = CD45GSEAGO[,c(2,5,6)]                          
CD45GSEAGO$'−log10(P-value)' <- -log10(CD45GSEAGO$pvalue)   
CD45GSEAGO <- CD45GSEAGO %>% arrange(NES)

# 手动改一下名字比较长的通路
CD45GSEAGO$Description[16] = "antigen processing and presentation via MHCII"
CD45GSEAGO$Description[16]

# 将通路设置为因子 保证绘图的排序和数据框一致
Description = CD45GSEAGO$Description
CD45GSEAGO <- 
  mutate(CD45GSEAGO, Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

# 画棒棒糖图
ggplot(CD45GSEAGO, aes(x = NES, y = Description)) + 
  geom_col(aes(fill = `−log10(P-value)`), width = 0.12, fill = "gray") +
  geom_point(aes(size = NES, fill = `−log10(P-value)`),
             shape = 21,
             stroke = 1.4,
             color = "gray11") +
  scale_size_continuous(range = c(2, 6), breaks = c(1.7, 1.75, 1.8, 1.85)) +
  scale_fill_distiller(palette = "RdPu", direction = 1,
                       limits = c(1, 4),
                       oob = scales::squish) +
  coord_cartesian(xlim = c(1.5, NA)) +
  theme_classic() +
  ylab('') +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 11.8, color = "black"),
    axis.text.y = element_text(size = 11.8, color = "black")
  )


# =============================================
# 🎨 Fig. 3d
# =============================================

# KEGG的GSEA条形码图展示
load("DEG_deseq2.Rdata")
CD45GSEAKEGG <- SCP::RunGSEA(
  geneID = rownames(DEG_deseq2),         
  geneScore = DEG_deseq2$log2FoldChange, 
  geneID_groups = rep("a",times=16330),  
                                        
  db = c("KEGG"),
  species = "Bos_taurus"
  )

# 可视化TOP通路
gseaplot2(CD45GSEAKEGG$results[["a-KEGG"]], 
          geneSetID = c("bta04613", "bta04932", "bta04148", "bta04514", "bta04657", "bta04621", "bta04620"),
          subplots = 1:3,
          color = c("#D73027", "#FC8D59", "#fec788", "#91BFDB", "#4575B4", "#66C2A5", "#1B9E77")
          )


# =============================================
# 🎨 Fig. 3e
# =============================================

# 8类细胞类型的两组间差异 已经做了差异基因分析 我们加载一下
# 以1倍以及0.05为阈值的话 差异基因如下
load("~/奶牛肝脏解离单细胞3版/group.dif.T.Rdata")  
table(group.dif.T$p_val < 0.05 & group.dif.T$avg_log2FC > log2(1)) 
table(group.dif.T$p_val < 0.05 & abs(group.dif.T$avg_log2FC) > log2(1)) # 2597 上调1089 下调1508
load("~/奶牛肝脏解离单细胞3版/group.dif.NK.Rdata") 
table(group.dif.NK$p_val < 0.05 & group.dif.NK$avg_log2FC > log2(1)) 
table(group.dif.NK$p_val < 0.05 & abs(group.dif.NK$avg_log2FC) > log2(1)) # 1592 上调604 下调988
load("~/奶牛肝脏解离单细胞3版/group.dif.B.Rdata")
table(group.dif.B$p_val < 0.05 & group.dif.B$avg_log2FC > log2(1)) 
table(group.dif.B$p_val < 0.05 & abs(group.dif.B$avg_log2FC) > log2(1)) # 1310 上调863 下调447
load("~/奶牛肝脏解离单细胞3版/group.dif.Plasma.Rdata")
table(group.dif.Plasma$p_val < 0.05 & group.dif.Plasma$avg_log2FC > log2(1)) 
table(group.dif.Plasma$p_val < 0.05 & abs(group.dif.Plasma$avg_log2FC) > log2(1)) # 1373 上调174 下调1199
load("~/奶牛肝脏解离单细胞3版/group.dif.Mono.Rdata")
table(group.dif.Mono$p_val < 0.05 & group.dif.Mono$avg_log2FC > log2(1)) 
table(group.dif.Mono$p_val < 0.05 & abs(group.dif.Mono$avg_log2FC) > log2(1)) # 3487 上调805 下调2682
load("~/奶牛肝脏解离单细胞3版/group.dif.MAC.Rdata")
table(group.dif.MAC$p_val < 0.05 & group.dif.MAC$avg_log2FC > log2(1)) 
table(group.dif.MAC$p_val < 0.05 & abs(group.dif.MAC$avg_log2FC) > log2(1)) # 1570 上调868 下调702
load("~/奶牛肝脏解离单细胞3版/group.dif.DC.Rdata")
table(group.dif.DC$p_val < 0.05 & group.dif.DC$avg_log2FC > log2(1)) 
table(group.dif.DC$p_val < 0.05 & abs(group.dif.DC$avg_log2FC) > log2(1)) # 2763 上调366 下调2397
load("~/奶牛肝脏解离单细胞3版/group.dif.Neu.Rdata")
table(group.dif.Neu$p_val < 0.05 & group.dif.Neu$avg_log2FC > log2(1)) 
table(group.dif.Neu$p_val < 0.05 & abs(group.dif.Neu$avg_log2FC) > log2(1)) # 697 上调135 下调562

# T
head(group.dif.T)
group.dif.T$gene <- rownames(group.dif.T)
group.dif.T$Group <- "T cell"
group.dif.T$Regulated <- ifelse(group.dif.T$p_val < 0.05 & group.dif.T$avg_log2FC > 0, "Up",
                                  ifelse(group.dif.T$p_val < 0.05 & group.dif.T$avg_log2FC < 0, "Down", "Stable"))
group.dif.T$log10fdr <- -log10(group.dif.T$p_val)
head(group.dif.T)
# NK
head(group.dif.NK)
group.dif.NK$gene <- rownames(group.dif.NK)
group.dif.NK$Group <- "NK"
group.dif.NK$Regulated <- ifelse(group.dif.NK$p_val < 0.05 & group.dif.NK$avg_log2FC > 0, "Up",
                                ifelse(group.dif.NK$p_val < 0.05 & group.dif.NK$avg_log2FC < 0, "Down", "Stable"))
group.dif.NK$log10fdr <- -log10(group.dif.NK$p_val)
head(group.dif.NK)
# B
head(group.dif.B)
group.dif.B$gene <- rownames(group.dif.B)
group.dif.B$Group <- "B cell"
group.dif.B$Regulated <- ifelse(group.dif.B$p_val < 0.05 & group.dif.B$avg_log2FC > 0, "Up",
                                ifelse(group.dif.B$p_val < 0.05 & group.dif.B$avg_log2FC < 0, "Down", "Stable"))
group.dif.B$log10fdr <- -log10(group.dif.B$p_val)
head(group.dif.B)
# PC
head(group.dif.Plasma)
group.dif.Plasma$gene <- rownames(group.dif.Plasma)
group.dif.Plasma$Group <- "PC"
group.dif.Plasma$Regulated <- ifelse(group.dif.Plasma$p_val < 0.05 & group.dif.Plasma$avg_log2FC > 0, "Up",
                                ifelse(group.dif.Plasma$p_val < 0.05 & group.dif.Plasma$avg_log2FC < 0, "Down", "Stable"))
group.dif.Plasma$log10fdr <- -log10(group.dif.Plasma$p_val)
head(group.dif.Plasma)
# Mono
head(group.dif.Mono)
group.dif.Mono$gene <- rownames(group.dif.Mono)
group.dif.Mono$Group <- "Mono"
group.dif.Mono$Regulated <- ifelse(group.dif.Mono$p_val < 0.05 & group.dif.Mono$avg_log2FC > 0, "Up",
                                     ifelse(group.dif.Mono$p_val < 0.05 & group.dif.Mono$avg_log2FC < 0, "Down", "Stable"))
group.dif.Mono$log10fdr <- -log10(group.dif.Mono$p_val)
head(group.dif.Mono)
# Mac
head(group.dif.MAC)
group.dif.MAC$gene <- rownames(group.dif.MAC)
group.dif.MAC$Group <- "Mac"
group.dif.MAC$Regulated <- ifelse(group.dif.MAC$p_val < 0.05 & group.dif.MAC$avg_log2FC > 0, "Up",
                                     ifelse(group.dif.MAC$p_val < 0.05 & group.dif.MAC$avg_log2FC < 0, "Down", "Stable"))
group.dif.MAC$log10fdr <- -log10(group.dif.MAC$p_val)
head(group.dif.MAC)
# DC
head(group.dif.DC)
group.dif.DC$gene <- rownames(group.dif.DC)
group.dif.DC$Group <- "DC"
group.dif.DC$Regulated <- ifelse(group.dif.DC$p_val < 0.05 & group.dif.DC$avg_log2FC > 0, "Up",
                                  ifelse(group.dif.DC$p_val < 0.05 & group.dif.DC$avg_log2FC < 0, "Down", "Stable"))
group.dif.DC$log10fdr <- -log10(group.dif.DC$p_val)
head(group.dif.DC)
# Neu
head(group.dif.Neu)
group.dif.Neu$gene <- rownames(group.dif.Neu)
group.dif.Neu$Group <- "Neu"
group.dif.Neu$Regulated <- ifelse(group.dif.Neu$p_val < 0.05 & group.dif.Neu$avg_log2FC > 0, "Up",
                                 ifelse(group.dif.Neu$p_val < 0.05 & group.dif.Neu$avg_log2FC < 0, "Down", "Stable"))
group.dif.Neu$log10fdr <- -log10(group.dif.Neu$p_val)
head(group.dif.Neu)

# 上下合并
diff_sc <- rbind(group.dif.T,
                 group.dif.NK,
                 group.dif.B,
                 group.dif.Plasma,
                 group.dif.Mono,
                 group.dif.MAC,
                 group.dif.DC,
                 group.dif.Neu)

# 设置diff_sc的水平信息
levels(diff_sc$Group)
diff_sc$Group <- factor(diff_sc$Group, levels = c("T cell", "NK", "B cell", "PC", "Mono", "Mac", "DC", "Neu"))
levels(diff_sc$Group)

# 展示上调基因
up_genes <- data.frame(Group = c(rep("T cell", 5), rep("NK", 5), rep("B cell", 5), rep("PC", 5), 
                                 rep("Mono", 5), rep("Mac", 5), rep("DC", 5), rep("Neu", 5)),
                       gene = c("CD69", "CXCR4", "XCL2", "ISG15", "IFI6",  # T cell
                                "XCL2", "CXCR4", "ISG15", "CSF2", "NUPR1", # NK
                                "FCRL3", "JCHAIN", "IGHM", "NFKBID", "ISG15", # B cell
                                "MKI67", "UBE2C", "CCL4", "PAX5", "CD83", # PC
                                "IFI6", "ISG15", "MX1", "ZBP1", "CCL16", # Mono
                                "LGALS3", "APOE", "FCGR2A", "ISG15", "ATF3", # Mac
                                "CLEC9A", "IL1A", "HMOX1", "PLAC8A", "ANPEP", # DC
                                "CXCL8", "CXCL2", "OSM", "GPR84", "HMOX1" # Neu
                                )
                       )

# 展示下调基因 
down_genes <- data.frame(Group = c(rep("T cell", 5), rep("NK", 5), rep("B cell", 5), rep("PC", 5),
                                   rep("Mono", 5), rep("Mac", 5), rep("DC", 5), rep("Neu", 5)),
                         gene = c("GNLY", "KLRG1", "FOXP1", "RORA", "PAG1",  # T cell
                                  "IL32", "DTX1", "SYTL2", "RICTOR", "PAG1", # NK
                                  "BACH2", "CD38", "DENND1B", "ZBTB20", "BOLA-DQB", # B cell
                                  "BTLA", "IL15RA", "NLRC5", "BOLA-DQA2", "IKZF2", # PC
                                  "CD36", "CLEC6A", "PELI2", "THBS1", "BOLA-DQB", # Mono
                                  "CD36", "CLEC4A", "CLEC6A", "BOLA-DQB", "DAB2", # Mac
                                  "BOLA-DQA5", "BOLA-DQB", "FSCN1", "CLEC7A", "CLEC4A", # DC
                                  "ABCA1", "FCGR3A", "FCAR", "PDK4", "IGSF6" # Neu
                                  )
                         )

# 从diff_sc中提取上下调挑选基因 
selected_up <- diff_sc %>%
  inner_join(up_genes, by = c("Group", "gene")) %>%
  arrange(Group, desc(avg_log2FC)) 
selected_down <- diff_sc %>%
  inner_join(down_genes, by = c("Group", "gene")) %>%
  arrange(Group, avg_log2FC)     

# 细胞类型的颜色
mycol <- c("#ae7eb5",
           "#f98177",
           "#aad393",
           "#568389",
           "#eeb066",
           "#f082a5",
           "#54c6a8",
           "#5a6fb5")

# 明确细胞顺序
cell_order <- c("T cell", "NK", "B cell", "PC", "Mono", "Mac", "DC", "Neu")
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
  ylab("-log10(padj)") +
  scale_x_continuous(limits = c(-5, 5)) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 0.8),
        strip.text.x = element_text(size = 10, face = 'bold'))

# 添加TOP基因
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
# 🎨 Fig. 3f
# =============================================

# 展示上下调基因个数
data <- data.frame(
  celltype = c("T cell", "NK", "B cell", "PC", 
               "Mono", "Mac", "DC", "Neu"),
  upregulated = c(1089, 604, 863, 174, 805, 868, 366, 135),
  downregulated = c(1508, 988, 447, 1199, 2682, 702, 2397, 562)
)
data

# 将数据转换为长格式
data_long <- reshape(data, 
                     varying = c("upregulated", "downregulated"), 
                     v.names = "gene_count", 
                     timevar = "gene_type", 
                     times = c("Upregulated", "Downregulated"),
                     direction = "long")
data_long$celltype <- factor(
  data_long$celltype,
  levels = c("T cell", "NK", "B cell", "PC", "Mono", "Mac", "DC", "Neu"))

# 绘制柱状图
ggplot(data_long, aes(x = celltype, y = gene_count, fill = factor(gene_type, levels = c("Upregulated", "Downregulated")))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 1) + 
  labs(x= NULL, y = "DEG count", fill = "Gene Regulation") +
  scale_fill_manual(values = c("Upregulated" = "#456baf", "Downregulated" = "#ef5688")) + 
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid = element_blank(),    
        axis.line = element_blank(),    
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13))


# =============================================
# 🎨 Fig. 3g-n
# =============================================

# 做GSEA分析  
# T
TGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.T),      
  geneScore = group.dif.T$avg_log2FC, 
  geneID_groups = rep("a",times=3364),     
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
TGSEAGO = TGSEAGO[["enrichment"]]
save(TGSEAGO, file = 'TGSEAGO.Rdata')
  
# NK
NKGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.NK),     
  geneScore = group.dif.NK$avg_log2FC,  
  geneID_groups = rep("a",times=3964),      
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
NKGSEAGO = NKGSEAGO[["enrichment"]]
save(NKGSEAGO, file = 'NKGSEAGO.Rdata')

# B
BGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.B),      
  geneScore = group.dif.B$avg_log2FC, 
  geneID_groups = rep("a",times=3896),
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
BGSEAGO = BGSEAGO[["enrichment"]]
save(BGSEAGO, file = 'BGSEAGO.Rdata')
  
# PC
PCGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.Plasma),      
  geneScore = group.dif.Plasma$avg_log2FC,  
  geneID_groups = rep("a",times=5407),       
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
PCGSEAGO = PCGSEAGO[["enrichment"]]
save(PCGSEAGO, file = 'PCGSEAGO.Rdata')

# Mono
MonoGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.Mono),      
  geneScore = group.dif.Mono$avg_log2FC,  
  geneID_groups = rep("a",times=6866),
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
    
# Mac
MacGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.MAC),     
  geneScore = group.dif.MAC$avg_log2FC,  
  geneID_groups = rep("a",times=2791),  
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
MacGSEAGO = MacGSEAGO[["enrichment"]]
save(MacGSEAGO, file = 'MacGSEAGO.Rdata')
  
# DC
DCGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.DC),      
  geneScore = group.dif.DC$avg_log2FC,  
  geneID_groups = rep("a",times=7406),      
  db = c("GO_BP"),
  species = "Bos_taurus"
  )  
DCGSEAGO = DCGSEAGO[["enrichment"]]
save(DCGSEAGO, file = 'DCGSEAGO.Rdata')
  
# Neu
NeuGSEAGO <- SCP::RunGSEA(
  geneID = rownames(group.dif.Neu),     
  geneScore = group.dif.Neu$avg_log2FC, 
  geneID_groups = rep("a",times=1434),  
  db = c("KEGG"),
  species = "Bos_taurus"
  )
NeuGSEAGO = NeuGSEAGO[["enrichment"]]
save(NeuGSEAGO, file = 'NeuGSEAGO.Rdata')


# T 上下调TOP通路各10条
load("TGSEAGO.Rdata")
TGSEAGO <- TGSEAGO[order(-TGSEAGO$NES), ]      
TGSEAGO <- TGSEAGO[TGSEAGO$pvalue <= 0.05, ]  
TGSEAGO = TGSEAGO[TGSEAGO$Description %in% c("positive regulation of T cell activation", 
                                      "adaptive immune response", 
                                      "immune effector process",
                                      "lymphocyte mediated immunity",
                                      "antigen processing and presentation",
                                      "antigen processing and presentation of peptide antigen", 
                                      "positive regulation of leukocyte cell-cell adhesion", 
                                      "oxidative phosphorylation", 
                                      "electron transport chain", 
                                      "generation of precursor metabolites and energy",
                                      "antigen receptor-mediated signaling pathway", 
                                      "signal transduction",
                                      "intracellular signal transduction",
                                      "cell communication",
                                      "small GTPase-mediated signal transduction",
                                      "regulation of GTPase activity",
                                      "Wnt signaling pathway",
                                      "regulation of transcription by RNA polymerase II",
                                      "positive regulation of RNA biosynthetic process",
                                      "protein deubiquitination"
                                      ),] 

TGSEAGO$'-log10pvalue' <- ifelse(TGSEAGO$NES < 0,
                                      -(-log10(TGSEAGO$pvalue)),
                                      -log10(TGSEAGO$pvalue)
                                      )

level_up <- TGSEAGO$Description[1:10]
level_down <- TGSEAGO$Description[11:20]
level <- c(level_up, level_down)
TGSEAGO$Description <- factor(TGSEAGO$Description, levels = rev(level))

ggplot(TGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") + 
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21,   
             stroke = 0.75,  
             color = "#939393",
             fill = ifelse(TGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(2, 6)) +
  ylab('') +
  scale_color_identity() +  
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)     

# NK 上下调TOP通路各10条
load("NKGSEAGO.Rdata")
NKGSEAGO <- NKGSEAGO[order(-NKGSEAGO$NES), ]      
NKGSEAGO <- NKGSEAGO[NKGSEAGO$pvalue <= 0.05, ]  
NKGSEAGO = NKGSEAGO[NKGSEAGO$Description %in% c("oxidative phosphorylation", 
                                      "cellular respiration", 
                                      "tricarboxylic acid cycle",
                                      "generation of precursor metabolites and energy",
                                      "immune response",
                                      "immune system process", 
                                      "positive regulation of leukocyte cell-cell adhesion", 
                                      "regulation of leukocyte cell-cell adhesion", 
                                      "antigen processing and presentation", 
                                      "antigen processing and presentation of peptide antigen",
                                      "intracellular signal transduction", 
                                      "immune response-activating cell surface receptor signaling pathway",
                                      "immune response-regulating signaling pathway",
                                      "secretion by cell",
                                      "cell surface receptor signaling pathway",
                                      "small GTPase-mediated signal transduction",
                                      "regulation of GTPase activity",
                                      "actin cytoskeleton organization",
                                      "cell motility",
                                      "exocytosis"
                                      ),] 

NKGSEAGO$'-log10pvalue' <- ifelse(NKGSEAGO$NES < 0,
                                      -(-log10(NKGSEAGO$pvalue)),
                                      -log10(NKGSEAGO$pvalue)
                                      )

level_up <- NKGSEAGO$Description[1:10]
level_down <- NKGSEAGO$Description[11:20]
level <- c(level_up, level_down)
NKGSEAGO$Description <- factor(NKGSEAGO$Description, levels = rev(level))

ggplot(NKGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") +
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21,  
             stroke = 0.75,   
             color = "#939393",
             fill = ifelse(NKGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(2, 6)) +
  ylab('') +
  scale_color_identity() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)      

# B 上下调通路各10
load("~/奶牛肝脏解离单细胞3版/BGSEAGO.Rdata")
BGSEAGO <- BGSEAGO[order(-BGSEAGO$NES), ]    
BGSEAGO <- BGSEAGO[BGSEAGO$pvalue <= 0.05, ]  
BGSEAGO = BGSEAGO[BGSEAGO$Description %in% c("oxidative phosphorylation", 
                                      "cellular respiration", 
                                      "tricarboxylic acid cycle",
                                      "generation of precursor metabolites and energy",
                                      "nucleotide biosynthetic process",
                                      "purine nucleotide biosynthetic process", 
                                      "reactive oxygen species metabolic process", 
                                      "Arp2/3 complex-mediated actin nucleation", 
                                      "response to external stimulus", 
                                      "apoptotic signaling pathway",
                                      "B cell mediated immunity", 
                                      "immunoglobulin production",
                                      "immunoglobulin mediated immune response",
                                      "antigen processing and presentation",
                                      "adaptive immune response",
                                      "regulation of lymphocyte activation",
                                      "positive regulation of lymphocyte activation",
                                      "immune effector process",
                                      "MAPK cascade",
                                      "small GTPase-mediated signal transduction"
                                      ),] 

BGSEAGO$'-log10pvalue' <- ifelse(BGSEAGO$NES < 0,
                                      -(-log10(BGSEAGO$pvalue)),
                                      -log10(BGSEAGO$pvalue)
                                      )

level_up <- BGSEAGO$Description[1:10]
level_down <- BGSEAGO$Description[11:20]
level <- c(level_up, level_down)
BGSEAGO$Description <- factor(BGSEAGO$Description, levels = rev(level))

ggplot(BGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") + 
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21,  
             stroke = 0.75,   
             color = "#939393",
             fill = ifelse(BGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(2, 5)) +
  ylab('') +
  scale_color_identity() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)    

# PC 上下调通路各10
load("PCGSEAGO.Rdata")
PCGSEAGO <- PCGSEAGO[order(-PCGSEAGO$NES), ]    
PCGSEAGO <- PCGSEAGO[PCGSEAGO$pvalue <= 0.05, ]  
PCGSEAGO = PCGSEAGO[PCGSEAGO$Description %in% c("immunoglobulin mediated immune response", 
                                      "B cell mediated immunity", 
                                      "adaptive immune response",
                                      "immune effector process",
                                      "inflammatory response",
                                      "defense response", 
                                      "response to cytokine", 
                                      "wound healing", 
                                      "B cell activation", 
                                      "antigen processing and presentation",
                                      "protein deubiquitination", 
                                      "protein modification by small protein conjugation or removal",
                                      "protein modification by small protein removal",
                                      "post-translational protein modification",
                                      "Notch signaling pathway",
                                      "regulation of transcription by RNA polymerase II",
                                      "transcription by RNA polymerase II",
                                      "cell-substrate adhesion",
                                      "DNA-templated transcription elongation",
                                      "regulation of cell projection organization"
                                      ),] 

PCGSEAGO$'-log10pvalue' <- ifelse(PCGSEAGO$NES < 0,
                                      -(-log10(PCGSEAGO$pvalue)),
                                      -log10(PCGSEAGO$pvalue)
                                      )

level_up <- PCGSEAGO$Description[1:10]
level_down <- PCGSEAGO$Description[11:20]
level <- c(level_up, level_down)
PCGSEAGO$Description <- factor(PCGSEAGO$Description, levels = rev(level))

ggplot(PCGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") +
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21,   
             stroke = 0.75,   
             color = "#939393",
             fill = ifelse(PCGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(2, 5)) +
  ylab('') +
  scale_color_identity() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)      

# Mono 上下调通路各10
load("~/奶牛肝脏解离单细胞3版/MonoGSEAGO.Rdata")
MonoGSEAGO <- MonoGSEAGO[order(-MonoGSEAGO$NES), ]   
MonoGSEAGO <- MonoGSEAGO[MonoGSEAGO$pvalue <= 0.05, ] 
MonoGSEAGO = MonoGSEAGO[MonoGSEAGO$Description %in% c("oxidative phosphorylation", 
                                      "cellular respiration", 
                                      "generation of precursor metabolites and energy",
                                      "translation",
                                      "ribosome biogenesis",
                                      "nucleoside triphosphate biosynthetic process", 
                                      "response to virus", 
                                      "defense response to virus", 
                                      "protein targeting to mitochondrion", 
                                      "multicellular organismal-level homeostasis",
                                      "antigen processing and presentation", 
                                      "antigen processing and presentation of peptide antigen via MHC class II",
                                      "antigen processing and presentation of exogenous antigen",
                                      "transforming growth factor beta receptor signaling pathway",
                                      "negative regulation of MAPK cascade",
                                      "receptor-mediated endocytosis",
                                      "regulated exocytosis",
                                      "small GTPase-mediated signal transduction",
                                      "regulation of cytoskeleton organization",
                                      "leukocyte cell-cell adhesion"
                                      ),] 

MonoGSEAGO$'-log10pvalue' <- ifelse(MonoGSEAGO$NES < 0,
                                      -(-log10(MonoGSEAGO$pvalue)),
                                      -log10(MonoGSEAGO$pvalue)
                                      )

level_up <- MonoGSEAGO$Description[1:10]
level_down <- MonoGSEAGO$Description[11:20]
level <- c(level_up, level_down)
MonoGSEAGO$Description <- factor(MonoGSEAGO$Description, levels = rev(level))

ggplot(MonoGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") + 
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21,   
             stroke = 0.75,   
             color = "#939393",
             fill = ifelse(MonoGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(3, 6)) +
  ylab('') +
  scale_color_identity() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)

# Mac 上下调通路各10
load("MacGSEAGO.Rdata")
MacGSEAGO <- MacGSEAGO[order(-MacGSEAGO$NES), ]      
MacGSEAGO <- MacGSEAGO[MacGSEAGO$pvalue <= 0.05, ]  
MacGSEAGO = MacGSEAGO[MacGSEAGO$Description %in% c("negative regulation of leukocyte activation", 
                                      "negative regulation of lymphocyte activation", 
                                      "negative regulation of T cell activation",
                                      "negative regulation of cell activation",
                                      "negative regulation of leukocyte cell-cell adhesion",
                                      "negative regulation of cell-cell adhesion", 
                                      "negative regulation of cell population proliferation", 
                                      "regulation of T cell proliferation", 
                                      "tissue development", 
                                      "animal organ development",
                                      "antigen processing and presentation", 
                                      "antigen processing and presentation of peptide antigen via MHC class II",
                                      "positive regulation of T cell activation",
                                      "positive regulation of lymphocyte activation",
                                      "immune effector process",
                                      "adaptive immune response",
                                      "lymphocyte mediated immunity",
                                      "positive regulation of leukocyte cell-cell adhesion",
                                      "receptor-mediated endocytosis",
                                      "regulation of GTPase activity"
                                      ),] 

MacGSEAGO$'-log10pvalue' <- ifelse(MacGSEAGO$NES < 0,
                                      -(-log10(MacGSEAGO$pvalue)),
                                      -log10(MacGSEAGO$pvalue)
                                      )

level_up <- MacGSEAGO$Description[1:10]
level_down <- MacGSEAGO$Description[11:20]
level <- c(level_up, level_down)
MacGSEAGO$Description <- factor(MacGSEAGO$Description, levels = rev(level))

ggplot(MacGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") + 
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21,  
             stroke = 0.75,  
             color = "#939393",
             fill = ifelse(MacGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(3, 4.5)) +
  ylab('') +
  scale_color_identity() +  
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)  # 只在x=0处添加虚线        

# DC 上下调通路各10
load("~/奶牛肝脏解离单细胞3版/DCGSEAGO.Rdata")
DCGSEAGO <- DCGSEAGO[order(-DCGSEAGO$NES), ]     
DCGSEAGO <- DCGSEAGO[DCGSEAGO$pvalue <= 0.05, ]
DCGSEAGO = DCGSEAGO[DCGSEAGO$Description %in% c("protein refolding", 
                                      "chaperone-mediated protein folding", 
                                      "'de novo' post-translational protein folding",
                                      "electron transport chain",
                                      "mitochondrial respiratory chain complex assembly",
                                      "nucleoside triphosphate biosynthetic process", 
                                      "pyrimidine nucleotide biosynthetic process", 
                                      "positive regulation of ERK1 and ERK2 cascade", 
                                      "regulation of Rho protein signal transduction", 
                                      "regulation of T cell proliferation",
                                      "antigen processing and presentation", 
                                      "antigen processing and presentation of peptide antigen via MHC class II",
                                      "antigen processing and presentation of exogenous antigen",
                                      "positive regulation of T cell activation",
                                      "positive regulation of lymphocyte activation",
                                      "adaptive immune response",
                                      "immune effector process",
                                      "pattern recognition receptor signaling pathway",
                                      "canonical NF-kappaB signal transduction",
                                      "actin cytoskeleton organization"
                                      ),] 

DCGSEAGO$'-log10pvalue' <- ifelse(DCGSEAGO$NES < 0,
                                      -(-log10(DCGSEAGO$pvalue)),
                                      -log10(DCGSEAGO$pvalue)
                                      )

level_up <- DCGSEAGO$Description[1:10]
level_down <- DCGSEAGO$Description[11:20]
level <- c(level_up, level_down)
DCGSEAGO$Description <- factor(DCGSEAGO$Description, levels = rev(level))

ggplot(DCGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") + 
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21, 
             stroke = 0.75, 
             color = "#939393",
             fill = ifelse(DCGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(3, 5.5)) +
  ylab('') +
  scale_color_identity() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)   

# Neu 上下调通路各10条
load("~/奶牛肝脏解离单细胞3版/NeuGSEAGO.Rdata")
NeuGSEAGO <- NeuGSEAGO[order(-NeuGSEAGO$NES), ]     
NeuGSEAGO <- NeuGSEAGO[NeuGSEAGO$pvalue <= 0.05, ]   
NeuGSEAGO = NeuGSEAGO[NeuGSEAGO$Description %in% c("IL-17 signaling pathway", 
                                      "TNF signaling pathway", 
                                      "NF-kappa B signaling pathway",
                                      "Cytokine-cytokine receptor interaction",
                                      "Toll-like receptor signaling pathway",
                                      "NOD-like receptor signaling pathway", 
                                      "Cytosolic DNA-sensing pathway", 
                                      "Complement and coagulation cascades", 
                                      "Efferocytosis", 
                                      "Lipid and atherosclerosis",
                                      "Phosphatidylinositol signaling system", 
                                      "Ras signaling pathway",
                                      "Notch signaling pathway",
                                      "Ubiquitin mediated proteolysis",
                                      "Mitophagy - animal",
                                      "Regulation of actin cytoskeleton",
                                      "Focal adhesion",
                                      "Adherens junction",
                                      "Gap junction",
                                      "FoxO signaling pathway"
                                      ),] 

NeuGSEAGO$'-log10pvalue' <- ifelse(NeuGSEAGO$NES < 0,
                                      -(-log10(NeuGSEAGO$pvalue)),
                                      -log10(NeuGSEAGO$pvalue)
                                      )

level_up <- NeuGSEAGO$Description[1:10]
level_down <- NeuGSEAGO$Description[11:20]
level <- c(level_up, level_down)
NeuGSEAGO$Description <- factor(NeuGSEAGO$Description, levels = rev(level))

ggplot(NeuGSEAGO,
             aes(x = NES, y = Description)) + 
  geom_col(width = 0.25, fill = "#d3d3d3") + 
  geom_point(aes(size = abs(`-log10pvalue`)),
             shape = 21,  
             stroke = 0.75,  
             color = "#939393",
             fill = ifelse(NeuGSEAGO$NES > 0, "#f38297", "#8ad0fa")) + 
  geom_point(aes(x = NES, y = Description), 
             color = "black", size = 0.1) + 
  scale_size_continuous(range = c(3, 5)) +
  ylab('') +
  scale_color_identity() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d3d3d3", size = 0.4)        

  
# =============================================
# 🎨 Supplementary Figure 6
# =============================================
  
library(aPEAR)
library(clusterProfiler)
library(org.Bt.eg.db)
library(DOSE) 
library(ggplot2)
library(grDevices) 

load("CD45GSEAGO.Rdata") 
set.seed(127) 
enrichmentNetwork(CD45GSEAGO[CD45GSEAGO$pvalue < 0.05,], 
                  drawEllipses = TRUE, 
                  repelLabels = T,
                  fontSize = 3
                  ) +
  scale_color_distiller(palette = "PiYG", direction = -1) 

load("CD45GSEAKEGG.Rdata") 
set.seed(125)  
enrichmentNetwork(CD45GSEAKEGG[CD45GSEAKEGG$pvalue < 0.05,], 
                  drawEllipses = TRUE, 
                  repelLabels = T,
                  fontSize = 3
                  ) +
  scale_color_distiller(palette = "RdPu", direction = 1) 


