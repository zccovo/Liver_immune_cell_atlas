
# F1/S1/S2/S3: Single-cell RNA sequencing maps the hepatic immune landscape.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 单细胞数据处理
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

# 数据导入
setwd("~/奶牛肝脏解离单细胞3版/NCBI三元文件")
data_dir <- "./H1"                      
list.files(data_dir)                    
H1 <- Read10X(data.dir = data_dir)      
dim(H1)                                 
class(H1)
str(H1)
H1[1:10,1:4]

data_dir <- "./H2"                      
list.files(data_dir)                    
H2 <- Read10X(data.dir = data_dir)      
dim(H2)                                 
class(H2)
str(H2)
H2[1:10,1:4]

data_dir <- "./H3"                      
list.files(data_dir)                    
H3 <- Read10X(data.dir = data_dir)      
dim(H3)                                
class(H3)
str(H3)
H3[1:10,1:4]

data_dir <- "./K1"                     
list.files(data_dir)                    
K1 <- Read10X(data.dir = data_dir)      
dim(K1)                                 
class(K1)
str(K1)
K1[1:10,1:4]

data_dir <- "./K2"                      
list.files(data_dir)                    
K2 <- Read10X(data.dir = data_dir)     
dim(K2)                                 
class(K2)
str(K2)
K2[1:10,1:4]

data_dir <- "./K3"                     
list.files(data_dir)                    
K3 <- Read10X(data.dir = data_dir)     
dim(K3)                                 
class(K3)
str(K3)
K3[1:10,1:4]

data_dir <- "./K4"                     
list.files(data_dir)                    
K4 <- Read10X(data.dir = data_dir)     
dim(K4)                                
class(K4)
str(K4)
K4[1:10,1:4]

data_dir <- "./K5"                      
list.files(data_dir)                    
K5 <- Read10X(data.dir = data_dir)     
dim(K5)                               
class(K5)
str(K5)
K5[1:10,1:4]

# 分别创建Seurat对象 再进行合并 创建时对细胞进行初步过滤
H1 <- CreateSeuratObject(counts = H1, 
                         project = "H1",      
                         min.cells = 3,        
                         min.features = 200   
                         ) # [1] 18185  9197
H2 <- CreateSeuratObject(counts = H2, 
                         project = "H2",       
                         min.cells = 3,        
                         min.features = 200 
                         ) # [1] 19007 13860
H3 <- CreateSeuratObject(counts = H3, 
                         project = "H3",       
                         min.cells = 3,        
                         min.features = 200 
                         ) # [1] 19182 15980
K1 <- CreateSeuratObject(counts = K1, 
                         project = "K1",       
                         min.cells = 3,        
                         min.features = 200 
                         ) # [1] 18709 15082
K2 <- CreateSeuratObject(counts = K2, 
                         project = "K2",       
                         min.cells = 3,        
                         min.features = 200 
                         ) # [1] 18913 11110
K3 <- CreateSeuratObject(counts = K3, 
                         project = "K3",       
                         min.cells = 3,        
                         min.features = 200 
                         ) # [1] 18630 13097
K4 <- CreateSeuratObject(counts = K4, 
                         project = "K4",       
                         min.cells = 3,        
                         min.features = 200 
                         ) # [1] 18913 12927
K5 <- CreateSeuratObject(counts = K5, 
                         project = "K5",       
                         min.cells = 3,        
                         min.features = 200 
                         ) # [1] 18630 10641

# 合并多个样本 给每个细胞Barcode信息加上前缀区分每个细胞
seurat_object <- merge(H1, 
                       y=c(H2,H3,K1,K2,K3,K4,K5), 
                       add.cell.ids = c("H1","H2","H3","K1","K2","K3","K4","K5"), 
                       project = "CD45")
seurat_object@assays$RNA$count[1:3,1:3]
dim(seurat_object) # [1] 20937 101894
table(seurat_object@meta.data[["orig.ident"]])

# 保存单细胞样本原始数据
save(seurat_object, file = "NCBI原始数据.Rdata") 
rm(list = ls())
load("NCBI原始数据.Rdata") 

#-------------------------------------------------------------------------------质控

# 计算每个细胞的线粒体基因转录本数的百分比(%) 牛的NCBI的线粒体基因是KEH36开头 如下
c("KEH36-t01", "KEH36-r02", "KEH36-t02", "KEH36-r01",
  "KEH36-t03", "KEH36-p13", "KEH36-t04", "KEH36-t05",
  "KEH36-t06", "KEH36-p12", "KEH36-t07", "KEH36-t08", 
  "KEH36-t09", "KEH36-t10", "KEH36-t11", "KEH36-p11", 
  "KEH36-t12", "KEH36-t13", "KEH36-p10", "KEH36-t14",
  "KEH36-p09", "KEH36-p08", "KEH36-p07", "KEH36-t15", 
  "KEH36-p06", "KEH36-t16", "KEH36-p05", "KEH36-p04", 
  "KEH36-t17", "KEH36-t18", "KEH36-t19", "KEH36-p03",
  "KEH36-p02", "KEH36-t20", "KEH36-p01", "KEH36-t21",
  "KEH36-t22") %in% rownames(seurat_object)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, 
                                                      pattern = "^KEH36" # 线粒体基因 提示受损/凋亡
                                                      )  

seurat_object[["percent.ribo"]] <- PercentageFeatureSet(seurat_object, 
                                                        pattern = "^RP[SL]" # 核糖体基因 一般不过滤 可能反应细胞活跃程度
                                                        )

hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ") # 血红蛋白基因 判断红系污染/红细胞成分
seurat_object[["percent.HB"]] <- PercentageFeatureSet(seurat_object, 
                                                      features = hb_genes
                                                      )

seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA) # log10GenesPerUMI 文库复杂度
  
# nFeature_RNA代表每个细胞测到的基因数目
# nCount_RNA代表每个细胞测到所有基因的UMI表达量之和
# percent.mt代表测到的线粒体基因的比例
VlnPlot(seurat_object, 
        features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI", "percent.mt", "percent.ribo", "percent.HB"), 
        ncol = 6
        ) 
  
# 过滤细胞 本数据为解离 查文献的线粒体阈值有的为 20/10/35/15
# 参数不是固定的 多少的都有 材料方法里表明即可  
dim(seurat_object) # 过滤前 [1] 20937 101894
seurat_object <- subset(seurat_object, 
                        subset = nFeature_RNA > 250 & 
                          nFeature_RNA < 4000 & 
                          nCount_RNA > 400 & 
                          nCount_RNA < 20000 & 
                          percent.mt < 25 & 
                          log10GenesPerUMI > 0.7 &
                          percent.HB < 5)   
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(seurat_object) # 过滤后 [1] 20785 94999  
  
# 对过滤后的QC metrics进行散点图可视化 高于0.9就没问题
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  
  
table(seurat_object@meta.data[["orig.ident"]])

#-------------------------------------------------------------------------------降维聚类

# 表达量数据标准化
seurat_object <- NormalizeData(seurat_object, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000 
                               )  

# 鉴定细胞间表达量高变的2000个基因 用于下游PCA分析
seurat_object <- FindVariableFeatures(seurat_object, 
                                      selection.method = "vst", 
                                      nfeatures = 2000
                                      ) 
head(VariableFeatures(seurat_object), 10) # 表达量变变化最高的10个基因

# 若不是研究关于细胞周期的课题 那么需要排除一下周期的影响 以免分出来的群是以周期分的
# Seurat R包自带的物种为人 其余物种的需要自己找文献或者按照人的转换过来
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# 同源基因转换 人转为牛
library(homologene) 
X = homologene(s.genes,        # 要转化的基因列表
               inTax = 9606,   # NCBI物种名 人
               outTax = 9913)  # NCBI物种名 牛
Y = homologene(g2m.genes,      # 要转化的基因列表
               inTax = 9606,   # NCBI物种名 人
               outTax = 9913)  # NCBI物种名 牛
# 提取结果
s.genes = X$"9913" ; g2m.genes = Y$"9913" 

# 细胞周期基因集打分(先合并S5分割矩阵)
seurat_object <- JoinLayers(seurat_object,assay = "RNA") 
seurat_object <- CellCycleScoring(seurat_object,
                                  s.features = s.genes,
                                  g2m.features = g2m.genes,
                                  set.ident = FALSE # 不将细胞周期状态设置为细胞身份active.ident
                                  )  
head(seurat_object[[]],20)

# 去除细胞周期对分群和降维的影响 同时对Seurat对象中的高可变基因进行归一化缩放处理
seurat_object <- ScaleData(seurat_object,                               # 对基因表达数据进行缩放和中心化处理
                           vars.to.regress = c("S.Score", "G2M.Score"), # 在缩放过程中 这两个评分的影响会被去除 以便减少细胞周期对数据分析的干扰
                           features = VariableFeatures(seurat_object)   # 指定要进行缩放的特征基因 这里使用的是在之前步骤中识别出的高变基因 这确保了仅对那些在细胞间具有显著变化的基因进行缩放处理
                           )

# 拿到2000个归一化的高变基因 但是不可能只拿这2000个基因去分群 相当于2000个维度 所以说要进行PCA线性降维 比如拿到PCA的前15个主成分 就可以代表这2000个高变基因95%权重 再进行后续分群即可
# 线性降维(PCA) 默认用高变基因集 但也可通过features参数自己指定 
seurat_object <- RunPCA(seurat_object,
                        features = VariableFeatures(object = seurat_object)) 
  
# 肘部图 主要目的是可视化每个主成分对数据总方差的贡献 
ElbowPlot(seurat_object, 50) 

# 基于PCA空间中的欧氏距离计算细胞之间的邻接关系 得到KNN和SNN 这是后面俩函数必不可少的输入信息
seurat_object <- FindNeighbors(seurat_object,
                               dims = 1:30 # 输入上一步得到的最优PC维数
                               )

# resolution参数决定下游聚类分析得到的分群数 相当于分辨率 推荐批量设置分辨率
seurat_object <- FindClusters(object = seurat_object,
                              resolution = c(seq(0.1, 1, 0.1)) # 这里的resolution控制分群数量
                              )
library(clustree)
clustree(seurat_object, prefix = "RNA_snn_res.") # 查看树状分群图
Idents(seurat_object) <- "RNA_snn_res.1"
seurat_object$seurat_clusters <- seurat_object@active.ident

# 查看前10个细胞的分群归属
head(Idents(seurat_object), 10) 

# UMAP非线性降维
seurat_object <- RunUMAP(seurat_object, dims = 1:25) # 这里的dims控制UMAP形状

# UMAP图
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(seurat_object, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)

# 保存数据
save(seurat_object, file = "NCBI降维聚类.Rdata") 
rm(list = ls())
load("NCBI降维聚类.Rdata")

# 大致看一下Marker
# CD45
Seurat::FeaturePlot(seurat_object, features = "PTPRC", reduction = "umap")
# T
Seurat::FeaturePlot(seurat_object, features = "CD3D", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD3E", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD3G", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CD3D","CD3E","CD3G"), group.by = "seurat_clusters")
# CD4T
Seurat::FeaturePlot(seurat_object, features = "CD4", reduction = "umap")
# CD8T
Seurat::FeaturePlot(seurat_object, features = "CD8A", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD8B", reduction = "umap")
# Pro_T
Seurat::FeaturePlot(seurat_object, features = "MKI67", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "TOP2A", reduction = "umap")
# B
Seurat::FeaturePlot(seurat_object, features = "MS4A1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD19", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD79B", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MS4A1","CD19","CD79B"), group.by = "seurat_clusters")
# Plasma_cell
Seurat::FeaturePlot(seurat_object, features = "MZB1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "JCHAIN", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MZB1","JCHAIN"), group.by = "seurat_clusters")
# NK
Seurat::FeaturePlot(seurat_object, features = "KLRB1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "KLRD1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "NCR1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "NKG7", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("KLRB1","KLRD1","NCR1"), group.by = "seurat_clusters")
# DC
Seurat::FeaturePlot(seurat_object, features = "CST3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FLT3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CLEC9A", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CST3","FLT3","CLEC9A"), group.by = "seurat_clusters")
# NEU
Seurat::FeaturePlot(seurat_object, features = "TGM3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CSF3R", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("TGM3","CSF3R"), group.by = "seurat_clusters")
# Mono
Seurat::FeaturePlot(seurat_object, features = "VCAN", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FCN1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD14", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CCR2", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("VCAN","FCN1"), group.by = "seurat_clusters")
# Macro
Seurat::FeaturePlot(seurat_object, features = "CD68", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD86", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD163", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CSF1R", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "C1QB", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "C1QC", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CD68","C1QB","C1QC"), group.by = "seurat_clusters")
# Hepatocyte
Seurat::FeaturePlot(seurat_object, features = "ALB", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "AFP", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CYP2E1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "HNF4A", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "ASGR1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "APOC3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FABP1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "APOA1", reduction = "umap")
# Cholangiocyte
Seurat::FeaturePlot(seurat_object, features = "KRT19", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "KRT7", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "EPCAM", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "SOX9", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("KRT19","KRT7","EPCAM","SOX9"), group.by = "seurat_clusters")
# LSEC
FeaturePlot(seurat_object, features = "LYVE1", reduction = "umap")
FeaturePlot(seurat_object, features = "STAB1", reduction = "umap")
FeaturePlot(seurat_object, features = "STAB2", reduction = "umap")
FeaturePlot(seurat_object, features = "FCGR2B", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("LYVE1","STAB1","STAB2","FCGR2B"), group.by = "seurat_clusters")
# Vascular EC
FeaturePlot(seurat_object, features = "VWF", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("VWF"), group.by = "seurat_clusters")
# VSMC
FeaturePlot(seurat_object, features = "MYH11", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MYH11"), group.by = "seurat_clusters")
# PF
FeaturePlot(seurat_object, features = "COL15A1", reduction = "umap")
FeaturePlot(seurat_object, features = "THY1", reduction = "umap")
FeaturePlot(seurat_object, features = "MSLN", reduction = "umap")
FeaturePlot(seurat_object, features = "ELN", reduction = "umap")
# HSC
FeaturePlot(seurat_object, features = "RGS5", reduction = "umap")
FeaturePlot(seurat_object, features = "COLEC11", reduction = "umap")
FeaturePlot(seurat_object, features = "HGF", reduction = "umap")

#-------------------------------------------------------------------------------移除非免疫细胞
  
# 现在通过细胞的Marker 可以知道不是免疫细胞的细胞
# 4/15/24/27 LSEC
# 31 VEC
# 23 Chol
# 30 VSMC

# 由于咱们是CD45分选的免疫细胞 所以把这些不是免疫细胞的细胞去掉
seurat_object <- subset(seurat_object, seurat_clusters %in% c("4","15","24","27","31","23","30"), invert = TRUE) # 反取子集  
  
# UMAP图
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(seurat_object, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)

# 保存数据
save(seurat_object, file = "CD45原始数据.Rdata") 
rm(list = ls())
load("CD45原始数据.Rdata")
  
#-------------------------------------------------------------------------------分群
  
# 初步尝试降维聚类分群已完成
# 样本分布极其相似 我决定不再用harmony进行去批次
# 重新跑标准流程 让UMAP图好看点
seurat_object <- FindNeighbors(seurat_object,
                               dims = 1:30 # 输入上一步得到的最优PC维数
                               )  

# resolution参数决定下游聚类分析得到的分群数 相当于分辨率 推荐批量设置分辨率
seurat_object <- FindClusters(object = seurat_object,
                              resolution = c(seq(0.1, 1, 0.1)) # 这里的resolution控制分群数量
                              )
library(clustree)
clustree(seurat_object, prefix = "RNA_snn_res.") # 查看树状分群图
Idents(seurat_object) <- "RNA_snn_res.0.8"
seurat_object$seurat_clusters <- seurat_object@active.ident

# 查看前10个细胞的分群归属
head(Idents(seurat_object), 10) 

# UMAP非线性降维
seurat_object <- RunUMAP(seurat_object, dims = 1:25) # 这里的dims控制UMAP形状

# UMAP图
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(seurat_object, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)

# 保存数据
save(seurat_object, file = "CD45降维聚类.Rdata") 
rm(list = ls())
load("CD45降维聚类.Rdata")   

# 大致看一下Marker
# CD45
Seurat::FeaturePlot(seurat_object, features = "PTPRC", reduction = "umap")
# T
Seurat::FeaturePlot(seurat_object, features = "CD3D", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD3E", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD3G", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CD3D","CD3E","CD3G"), group.by = "seurat_clusters")
# CD4T
Seurat::FeaturePlot(seurat_object, features = "CD4", reduction = "umap")
# CD8T
Seurat::FeaturePlot(seurat_object, features = "CD8A", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD8B", reduction = "umap")
# Pro_T
Seurat::FeaturePlot(seurat_object, features = "MKI67", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "TOP2A", reduction = "umap")
# B
Seurat::FeaturePlot(seurat_object, features = "MS4A1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD19", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD79B", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MS4A1","CD19","CD79B"), group.by = "seurat_clusters")
# Plasma_cell
Seurat::FeaturePlot(seurat_object, features = "MZB1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "JCHAIN", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MZB1","JCHAIN"), group.by = "seurat_clusters")
# NK
Seurat::FeaturePlot(seurat_object, features = "KLRB1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "KLRD1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "NCR1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "NKG7", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("KLRB1","KLRD1","NCR1"), group.by = "seurat_clusters")
# DC
Seurat::FeaturePlot(seurat_object, features = "CST3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FLT3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CLEC9A", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CST3","FLT3","CLEC9A"), group.by = "seurat_clusters")
# NEU
Seurat::FeaturePlot(seurat_object, features = "TGM3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CSF3R", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("TGM3","CSF3R"), group.by = "seurat_clusters")
# Mono
Seurat::FeaturePlot(seurat_object, features = "VCAN", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FCN1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD14", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CCR2", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("VCAN","FCN1"), group.by = "seurat_clusters")
# Macro
Seurat::FeaturePlot(seurat_object, features = "CD68", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD86", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD163", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CSF1R", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "C1QB", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "C1QC", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CD68","C1QB","C1QC"), group.by = "seurat_clusters")  
# Hepatocyte
Seurat::FeaturePlot(seurat_object, features = "ALB", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "AFP", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CYP2E1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "HNF4A", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "ASGR1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "APOC3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FABP1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "APOA1", reduction = "umap")
# Cholangiocyte
Seurat::FeaturePlot(seurat_object, features = "KRT19", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "KRT7", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "EPCAM", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "SOX9", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("KRT19","KRT7","EPCAM","SOX9"), group.by = "seurat_clusters")
# LSEC
FeaturePlot(seurat_object, features = "LYVE1", reduction = "umap")
FeaturePlot(seurat_object, features = "STAB1", reduction = "umap")
FeaturePlot(seurat_object, features = "STAB2", reduction = "umap")
FeaturePlot(seurat_object, features = "FCGR2B", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("LYVE1","STAB1","STAB2","FCGR2B"), group.by = "seurat_clusters")
# Vascular EC
FeaturePlot(seurat_object, features = "VWF", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("VWF"), group.by = "seurat_clusters")
# VSMC
FeaturePlot(seurat_object, features = "MYH11", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MYH11"), group.by = "seurat_clusters")
# PF
FeaturePlot(seurat_object, features = "COL15A1", reduction = "umap")
FeaturePlot(seurat_object, features = "THY1", reduction = "umap")
FeaturePlot(seurat_object, features = "MSLN", reduction = "umap")
FeaturePlot(seurat_object, features = "ELN", reduction = "umap")
# HSC
FeaturePlot(seurat_object, features = "RGS5", reduction = "umap")
FeaturePlot(seurat_object, features = "COLEC11", reduction = "umap")
FeaturePlot(seurat_object, features = "HGF", reduction = "umap")
  
#-------------------------------------------------------------------------------细胞鉴定
  
# 对未注释前细胞cluster进行差异基因分析 用来鉴定细胞
# FindAllMarkers函数旨在识别在不同细胞群体之间显著差异表达的基因 默认wilcox方法 将所有满足FC条件的基因筛选出来 没有设定P值范围
dif<-FindAllMarkers(seurat_object, 
                    group.by = seurat_object@meta.data$seurat_clusters, 
                    logfc.threshold = log2(1.2),                        
                    min.pct = 0.2,                                      
                    only.pos = T                                        
                    )                                                   
dif$pct_diff <- dif$pct.1 - dif$pct.2 
table(dif$cluster)                    
# 通过dplyr包对每组基因按差异倍数从大到小排序
dif<-dif %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)
# 保存差异分析结果
save(dif, file = "CD45整体cluster差异基因.Rdata") 
  
# 用经典Marker分大群 再加载细胞差异分析结果 找小Marker画图
load("~/奶牛肝脏解离单细胞3版/CD45整体cluster差异基因.Rdata")
  
# T 0/2/6/7/8/9/17/21
Seurat::FeaturePlot(seurat_object, features = "CD3D", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD3E", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD3G", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "ZAP70", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "BCL11B", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CD3D","CD3E","CD3G"), group.by = "seurat_clusters")
# NK 13
Seurat::FeaturePlot(seurat_object, features = "KLRB1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "KLRD1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "NCR1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "NKG7", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("KLRB1","KLRD1","NCR1","NKG7"), group.by = "seurat_clusters")
# B 12
Seurat::FeaturePlot(seurat_object, features = "MS4A1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD19", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD79B", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MS4A1","CD19","CD79B"), group.by = "seurat_clusters")
# Plasma_cell 18
Seurat::FeaturePlot(seurat_object, features = "MZB1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "JCHAIN", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("MZB1","JCHAIN"), group.by = "seurat_clusters")
# DC 15/22
Seurat::FeaturePlot(seurat_object, features = "CST3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FLT3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CLEC9A", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CST3","FLT3","CLEC9A"), group.by = "seurat_clusters")
# NEU 14/23
Seurat::FeaturePlot(seurat_object, features = "TGM3", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CSF3R", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("TGM3","CSF3R"), group.by = "seurat_clusters")
# Mono 5
Seurat::FeaturePlot(seurat_object, features = "VCAN", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "FCN1", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD14", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CCR2", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("VCAN","FCN1"), group.by = "seurat_clusters")
# Macro 1/3/4/10/11/16/17/19
Seurat::FeaturePlot(seurat_object, features = "CD68", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD86", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CD163", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "CSF1R", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "C1QB", reduction = "umap")
Seurat::FeaturePlot(seurat_object, features = "C1QC", reduction = "umap")
Seurat::DotPlot(seurat_object, features = c("CD68","C1QB","C1QC","CD86","CD163"), group.by = "seurat_clusters")    
  
# 大致知道了细胞类型组成
# 手动去掉已知的双细胞群体
# 17
seurat_subset <- subset(seurat_object, 
                        seurat_clusters %in% c("0","1","2","3","4","5","6","7","8","9","10","11","12","13",
                                               "14","15","16","18","19","20","21","22","23"))
# 查看提取后的UMAP图
Seurat::DimPlot(seurat_subset, reduction = "umap", pt.size = 0.1, label = T)
table(seurat_subset@active.ident)

# 明面上的双细胞群体已经删除了 现在正式进行细胞鉴定
Seurat::DimPlot(seurat_subset, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)    
  
# 假设我已经对每个Cluster鉴定完成了 接下来即可将数字替换为细胞名
# 为分群重新指定细胞类型 
new.cluster.ids <- c("T",
                     "Mac",
                     "T",
                     "Mac",
                     "Mac",
                     "Mono",
                     "T",
                     "T",
                     "T",
                     "T",
                     "Mac",
                     "Mac",
                     "B",
                     "NK",
                     "Neu",
                     "DC",
                     "Mac",
                     "PC",
                     "Mac",
                     "Mac",
                     "T",
                     "DC",
                     "Neu") # 自定义名称 按数字从小到大对应  
new.cluster.ids
names(new.cluster.ids) 
levels(seurat_subset)
names(new.cluster.ids) <- levels(seurat_subset) # 将seurat_object的水平属性赋值给new.cluster.ids的names属性 
names(new.cluster.ids)
new.cluster.ids
seurat_subset <- RenameIdents(seurat_subset, new.cluster.ids) # 给cluster亚群进行细胞命名  
  
# 为metadata添加一些有用信息
# 在meta.data中加上celltype卡槽
seurat_subset[["celltype"]]<-Idents(seurat_subset) 
# 在meta.data中加上group卡槽
# H1/H2/H3为对照组 7068+11807+13490 = 32365
# K1/K2/K3/K4/K5为模型组 12530+8884+10556+10337+8361 = 50668
unique(seurat_subset@meta.data[["orig.ident"]])
seurat_subset[["group"]]<- c(rep("Health", 32365), rep("Ketosis", 50668))
# 在meta.data中加上group.cluster卡槽 以同时区分分组名称和细胞名称
seurat_subset[["group.celltype"]]<-paste(seurat_subset$group, Idents(seurat_subset), sep = '_')   

# 查看一下统计结果
table(seurat_subset@meta.data[["orig.ident"]])
table(seurat_subset@meta.data[["group"]])
table(seurat_subset@meta.data[["celltype"]])
table(seurat_subset@meta.data[["group.celltype"]])  
  
# 绘制总umap图
Seurat::DimPlot(seurat_subset, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_subset, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)
seurat_object = seurat_subset
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

# 保存工作空间
save(seurat_object,file = "CD45细胞鉴定.Rdata") # 标注好 是细胞鉴定后的结果
rm(list = ls())
load("CD45细胞鉴定.Rdata")  

#------------------------------------------------------------------------------- 去双细胞

# 加载R包
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)

# 加载标准流程降维聚类后的数据
load("CD45细胞鉴定.Rdata") 

# 将Seurat对象转换为SingleCellExperiment对象
sce <- as.SingleCellExperiment(seurat_object)

# 如果是多个样本 需要为每个样本分别进行双细胞识别 则以下运行
sce <- scDblFinder(sce,
                   samples="orig.ident", # 按orig.ident列的样本分组分别处理双细胞
                                         # 确保人工生成的双细胞不会跨样本组合 现实中双细胞只可能来自同一样本的细胞 
                                         # 不同样本的双细胞比例可能不同 需分别估计
                   BPPARAM=MulticoreParam(10) # 使用10个CPU核心并行处理任务
                   ) # 结果会在sce@colData@listData中添加某些以scDblFinder为前缀的列 最重要的是以下两个 双细胞数和双细胞评分

# 查看双细胞数
table(sce@colData@listData[["scDblFinder.class"]]) # 判定分类 双细胞（"doublet"）与真实细胞（"singlet"）

# 查看双细胞评分
head(sce@colData@listData[["scDblFinder.score"]]) # 双细胞概率（0~1, 越接近1越可能是双细胞）

# 将分析结果添加至原始Seurat中而不是转成Seurat对象
seurat_object <- AddMetaData(seurat_object,
                             metadata = sce@colData@listData[["scDblFinder.class"]],
                             col.name = "scDblFinder.class")
seurat_object <- AddMetaData(seurat_object,
                             metadata = sce@colData@listData[["scDblFinder.score"]],
                             col.name = "scDblFinder.score")

# 保存双细胞分析结果
save(seurat_object,file = "CD45去双细胞前.Rdata") 
rm(list = ls())
load("CD45去双细胞前.Rdata")

# 在UMAP上标记双细胞
Seurat::DimPlot(seurat_object,
                reduction = "umap",                              # 使用UMAP降维
                group.by = "scDblFinder.class",                  # 指定双细胞分类列
                cols = c("singlet" = "gray", "doublet" = "red"), # 颜色自定义
                pt.size = 0.1) +                                 # 点大小
  ggtitle("UMAP: Singlets vs Doublets")                          # 名称

# 过滤双细胞 方便进行后续分析
seurat_scDblFinder <- subset(seurat_object, scDblFinder.class == "singlet")
Seurat::DimPlot(seurat_scDblFinder, reduction = "umap", group.by = "celltype", label = T) +
  Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "celltype", label = T)

# UMAP
seurat_object = seurat_scDblFinder
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "group", pt.size = 0.1, label = F)
Seurat::DimPlot(seurat_object, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

--------------------------------------------------------------------------------# 重新FindClusters

# resolution参数决定下游聚类分析得到的分群数 相当于分辨率 推荐批量设置分辨率
seurat_object <- FindClusters(object = seurat_object,
                              resolution = c(seq(0.1, 2, 0.1))
                              )
library(clustree)
clustree(seurat_object, prefix = "RNA_snn_res.") # 查看树状分群图
seurat_object$seurat_clusters <- seurat_object$RNA_snn_res.0.8
Idents(seurat_object) = seurat_object$celltype

# 查看前10个细胞的分群归属
head(Idents(seurat_object), 10) 

# 重新进行UMAP非线性降维 
seurat_object <- RunUMAP(seurat_object, dims = 1:30) # 这里的dims控制UMAP形状

# UMAP图
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", pt.size = 0.1, label = F)
Seurat::DimPlot(seurat_object, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, label = F)

# 删掉零零碎碎没有意义的小亚群
seurat_object <- subset(seurat_object, seurat_clusters %in% c("21", "23"), invert = TRUE) # 反取子集  

# UMAP图
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 把分辨率大于1的删除 占内存
seurat_object@meta.data[["RNA_snn_res.1.1"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.2"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.3"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.4"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.5"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.6"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.7"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.8"]] = NULL
seurat_object@meta.data[["RNA_snn_res.1.9"]] = NULL
seurat_object@meta.data[["RNA_snn_res.2"]] = NULL

# 重新FindClusters 设置最大值为1即可 
seurat_object <- FindClusters(object = seurat_object,
                              resolution = c(seq(0.1, 1, 0.1)) # 这里的resolution控制分群数量
                              )
library(clustree)
clustree(seurat_object, prefix = "RNA_snn_res.") # 查看树状分群图
seurat_object$seurat_clusters <- seurat_object$RNA_snn_res.0.9 # 23群
Idents(seurat_object) = seurat_object$celltype

# UMAP图
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 查看一下统计结果
table(seurat_object@meta.data[["orig.ident"]])
table(seurat_object@meta.data[["group"]])
table(seurat_object@meta.data[["celltype"]])
table(seurat_object@meta.data[["group.celltype"]])
table(seurat_object@meta.data[["seurat_clusters"]])

# 更改细胞名称的水平信息
Idents(seurat_object) <- factor(
  Idents(seurat_object),
  levels = c("T", "NK",
             "B", "PC",
             "Mono", "Mac", 
             "DC", "Neu"))
seurat_object[["celltype"]] <- Idents(seurat_object) 

# 保存数据
save(seurat_object,file = "CD45去双细胞中.Rdata") # 标注好 是scDblFinder方法去双细胞后的结果
rm(list = ls())
load("CD45去双细胞中.Rdata") 

#-------------------------------------------------------------------------------重新赋予细胞名称  

Seurat::DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = T)

# 需要重新改一下细胞名称 这样才能和Cluster对应
lvls <- levels(seurat_object$seurat_clusters)
lvls_new <- as.character(sort(as.numeric(lvls)))
seurat_object$seurat_clusters <- factor(
  seurat_object$seurat_clusters,
  levels = lvls_new
  )
levels(seurat_object$seurat_clusters)
table(seurat_object@meta.data[["seurat_clusters"]])
Idents(seurat_object) = seurat_object$seurat_clusters
# 假设我已经对每个Cluster鉴定完成了 接下来即可将数字替换为细胞名
# 为分群重新指定细胞类型 
new.cluster.ids <- c("Mac",
                     "Mac",
                     "T",
                     "T",
                     "T",
                     "Mono",
                     "Mac",
                     "T",
                     "T",
                     "T",
                     "T",
                     "T",
                     "B",
                     "Mac",
                     "NK",
                     "Neu",
                     "DC",
                     "Mac",
                     "Mac",
                     "PC",
                     "Mac",
                     "Mac",
                     "DC") 
new.cluster.ids
names(new.cluster.ids) 
levels(seurat_object)
names(new.cluster.ids) <- levels(seurat_object) 
names(new.cluster.ids)
new.cluster.ids
seurat_object <- RenameIdents(seurat_object, new.cluster.ids) 
  
# 为metadata添加一些有用信息
# 在meta.data中加上celltype卡槽
seurat_object[["celltype"]]<-Idents(seurat_object) 
# 在meta.data中加上group卡槽
# H1/H2/H3为对照组 6478+10360+11676 = 28514
# K1/K2/K3/K4/K5为模型组 10960+8113+9263+9334+7538 = 45208
unique(seurat_object@meta.data[["orig.ident"]])
table(seurat_object@meta.data[["orig.ident"]])
seurat_object[["group"]]<- c(rep("Health", 28514), rep("Ketosis", 45208))
# 在meta.data中加上group.cluster卡槽 以同时区分分组名称和细胞名称
seurat_object[["group.celltype"]]<-paste(seurat_object$group, Idents(seurat_object), sep = '_')   

# 查看一下统计结果
table(seurat_object@meta.data[["orig.ident"]])
table(seurat_object@meta.data[["group"]])
table(seurat_object@meta.data[["celltype"]])
table(seurat_object@meta.data[["group.celltype"]])  

# 细胞水平信息
Idents(seurat_object) <- factor(Idents(seurat_object),
                                levels = c("T", "NK", "B", "PC", "Mono", "Mac", "DC", "Neu"))
seurat_object[["celltype"]] <- Idents(seurat_object) 
  
# 绘制总umap图
Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T)
Seurat::DimPlot(seurat_object, reduction = "umap", split.by = "group", pt.size = 0.1, label = F)

# 重新保存工作空间
save(seurat_object,file = "CD45去双细胞后.Rdata") 
rm(list = ls())
load("CD45去双细胞后.Rdata")  

#-------------------------------------------------------------------------------整体细胞差异分析  

# 对注释后细胞亚群 进行差异基因分析
diff <- FindAllMarkers(seurat_object, 
                      group.by = seurat_object@active.ident, 
                      logfc.threshold = log2(1.2), 
                      min.pct = 0.2, 
                      only.pos = T 
                      )
diff$pct_diff <- diff$pct.1 - diff$pct.2 
table(diff$cluster) 
head(diff[diff$cluster == unique(diff$cluster)[1],], 50) $ gene 
# 通过dplyr包对每组基因按差异倍数从大到小排序
diff<-diff %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)
# 保存差异分析结果
save(diff, file = "CD45整体细胞类型差异基因.Rdata")  

#-------------------------------------------------------------------------------组间差异基因分析
  
# 对所有的组间细胞亚群进行差异基因分析
load("CD45去双细胞后.Rdata")
table(seurat_object@active.ident)
table(seurat_object@meta.data[["group.celltype"]])

# T cells
group.dif.T <- FindMarkers(seurat_object, 
                           group.by = "group.celltype", 
                           ident.1 = "Ketosis_T",     
                           ident.2 = "Health_T",     
                           logfc.threshold = 0,     
                           min.pct = 0.1,           
                           only.pos = FALSE         
                           ) 
group.dif.T$pct_diff <- group.dif.T$pct.1 - group.dif.T$pct.2 
group.dif.T<-group.dif.T %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.T, file = 'group.dif.T.Rdata')
write.table(group.dif.T,"group.dif.T.xls",row.names = T,col.names = NA,quote = F,sep = "\t") 
table(group.dif.T$p_val < 0.05 & group.dif.T$avg_log2FC > log2(1.2))
table(group.dif.T$p_val < 0.05 & abs(group.dif.T$avg_log2FC) > log2(1.2)) # 1145差异基因 256上调 889下调
table(group.dif.T$p_val < 0.05 & group.dif.T$avg_log2FC > 0)
table(group.dif.T$p_val < 0.05 & abs(group.dif.T$avg_log2FC) > 0) # 2597差异基因 1089上调 1508下调

# NK cells
group.dif.NK <- FindMarkers(seurat_object, 
                           group.by = "group.celltype", 
                           ident.1 = "Ketosis_NK",      
                           ident.2 = "Health_NK",      
                           logfc.threshold = 0,
                           min.pct = 0.1,
                           only.pos = FALSE 
                           ) 
group.dif.NK$pct_diff <- group.dif.NK$pct.1 - group.dif.NK$pct.2 
group.dif.NK<-group.dif.NK %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.NK, file = 'group.dif.NK.Rdata')
write.table(group.dif.NK,"group.dif.NK.xls",row.names = T,col.names = NA,quote = F,sep = "\t") 
table(group.dif.NK$p_val < 0.05 & group.dif.NK$avg_log2FC > log2(1.2))
table(group.dif.NK$p_val < 0.05 & abs(group.dif.NK$avg_log2FC) > log2(1.2)) # 1355差异基因 497上调 858下调
table(group.dif.NK$p_val < 0.05 & group.dif.NK$avg_log2FC > 0)
table(group.dif.NK$p_val < 0.05 & abs(group.dif.NK$avg_log2FC) > 0) # 1592差异基因 604上调 988下调

# B cells
group.dif.B <- FindMarkers(seurat_object, 
                           group.by = "group.celltype", 
                           ident.1 = "Ketosis_B",     
                           ident.2 = "Health_B",       
                           logfc.threshold = 0,
                           min.pct = 0.1,
                           only.pos = FALSE 
                           ) 
group.dif.B$pct_diff <- group.dif.B$pct.1 - group.dif.B$pct.2 
group.dif.B<-group.dif.B %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.B, file = 'group.dif.B.Rdata')
write.table(group.dif.B,"group.dif.B.xls",row.names = T,col.names = NA,quote = F,sep = "\t")
table(group.dif.B$p_val < 0.05 & group.dif.B$avg_log2FC > log2(1.2)) 
table(group.dif.B$p_val < 0.05 & abs(group.dif.B$avg_log2FC) > log2(1.2)) # 差异基因876 484上调 392下调
table(group.dif.B$p_val < 0.05 & group.dif.B$avg_log2FC > 0) 
table(group.dif.B$p_val < 0.05 & abs(group.dif.B$avg_log2FC) > 0) # 差异基因1310 863上调 447下调

# Plasma cells
group.dif.Plasma <- FindMarkers(seurat_object, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_PC",     
                            ident.2 = "Health_PC",       
                            logfc.threshold = 0,
                            min.pct = 0.1,
                            only.pos = FALSE 
                            ) 
group.dif.Plasma$pct_diff <- group.dif.Plasma$pct.1 - group.dif.Plasma$pct.2 
group.dif.Plasma<-group.dif.Plasma %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.Plasma, file = 'group.dif.Plasma.Rdata')
write.table(group.dif.Plasma,"group.dif.Plasma.xls",row.names = T,col.names = NA,quote = F,sep = "\t")
table(group.dif.Plasma$p_val < 0.05 & group.dif.Plasma$avg_log2FC > log2(1.2)) 
table(group.dif.Plasma$p_val < 0.05 & abs(group.dif.Plasma$avg_log2FC) > log2(1.2)) # 差异基因1173 153上调 1020下调
table(group.dif.Plasma$p_val < 0.05 & group.dif.Plasma$avg_log2FC > 0) 
table(group.dif.Plasma$p_val < 0.05 & abs(group.dif.Plasma$avg_log2FC) > 0) # 差异基因1373 174上调 1199下调

# Monocytes
group.dif.Mono <- FindMarkers(seurat_object, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_Mono",     
                            ident.2 = "Health_Mono",       
                            logfc.threshold = 0,
                            min.pct = 0.1,
                            only.pos = FALSE 
                            ) 
group.dif.Mono$pct_diff <- group.dif.Mono$pct.1 - group.dif.Mono$pct.2 
group.dif.Mono<-group.dif.Mono %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.Mono, file = 'group.dif.Mono.Rdata')
write.table(group.dif.Mono,"group.dif.Mono.xls",row.names = T,col.names = NA,quote = F,sep = "\t")
table(group.dif.Mono$p_val < 0.05 & group.dif.Mono$avg_log2FC > log2(1.2)) 
table(group.dif.Mono$p_val < 0.05 & abs(group.dif.Mono$avg_log2FC) > log2(1.2)) # 差异基因1684 490上调 1194下调
table(group.dif.Mono$p_val < 0.05 & group.dif.Mono$avg_log2FC > 0) 
table(group.dif.Mono$p_val < 0.05 & abs(group.dif.Mono$avg_log2FC) > 0) # 差异基因3487 805上调 2682下调

# Macrophages
group.dif.MAC <- FindMarkers(seurat_object, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_Mac",     
                            ident.2 = "Health_Mac",       
                            logfc.threshold = 0,
                            min.pct = 0.1,
                            only.pos = FALSE 
                            ) 
group.dif.MAC$pct_diff <- group.dif.MAC$pct.1 - group.dif.MAC$pct.2 
group.dif.MAC<-group.dif.MAC %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.MAC, file = 'group.dif.MAC.Rdata')
write.table(group.dif.MAC,"group.dif.MAC.xls",row.names = T,col.names = NA,quote = F,sep = "\t")
table(group.dif.MAC$p_val < 0.05 & group.dif.MAC$avg_log2FC > log2(1.2)) 
table(group.dif.MAC$p_val < 0.05 & abs(group.dif.MAC$avg_log2FC) > log2(1.2)) # 差异基因552 195上调 357下调
table(group.dif.MAC$p_val < 0.05 & group.dif.MAC$avg_log2FC > 0) 
table(group.dif.MAC$p_val < 0.05 & abs(group.dif.MAC$avg_log2FC) > 0) # 差异基因1570 868上调 702下调

# Dendritic cells
group.dif.DC <- FindMarkers(seurat_object, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_DC",     
                            ident.2 = "Health_DC",       
                            logfc.threshold = 0,
                            min.pct = 0.1,
                            only.pos = FALSE 
                            ) 
group.dif.DC$pct_diff <- group.dif.DC$pct.1 - group.dif.DC$pct.2 
group.dif.DC<-group.dif.DC %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.DC, file = 'group.dif.DC.Rdata')
write.table(group.dif.DC,"group.dif.DC.xls",row.names = T,col.names = NA,quote = F,sep = "\t")
table(group.dif.DC$p_val < 0.05 & group.dif.DC$avg_log2FC > log2(1.2)) 
table(group.dif.DC$p_val < 0.05 & abs(group.dif.DC$avg_log2FC) > log2(1.2)) # 差异基因1537 221上调 1316下调
table(group.dif.DC$p_val < 0.05 & group.dif.DC$avg_log2FC > 0) 
table(group.dif.DC$p_val < 0.05 & abs(group.dif.DC$avg_log2FC) > 0) # 差异基因2763 366上调 2397下调

# Neutrophil
group.dif.Neu <- FindMarkers(seurat_object, 
                            group.by = "group.celltype", 
                            ident.1 = "Ketosis_Neu",     
                            ident.2 = "Health_Neu",       
                            logfc.threshold = 0,
                            min.pct = 0.1,
                            only.pos = FALSE 
                            ) 
group.dif.Neu$pct_diff <- group.dif.Neu$pct.1 - group.dif.Neu$pct.2 
group.dif.Neu<-group.dif.Neu %>% arrange(desc(avg_log2FC), .by_group = TRUE) 
save(group.dif.Neu, file = 'group.dif.Neu.Rdata')
write.table(group.dif.Neu,"group.dif.Neu.xls",row.names = T,col.names = NA,quote = F,sep = "\t")
table(group.dif.Neu$p_val < 0.05 & group.dif.Neu$avg_log2FC > log2(1.2)) 
table(group.dif.Neu$p_val < 0.05 & abs(group.dif.Neu$avg_log2FC) > log2(1.2)) # 差异基因625 121上调 504下调
table(group.dif.Neu$p_val < 0.05 & group.dif.Neu$avg_log2FC > 0) 
table(group.dif.Neu$p_val < 0.05 & abs(group.dif.Neu$avg_log2FC) > 0) # 差异基因697 135上调 562下调

# =============================================
# 🎨 Fig. 1c
# =============================================

load("CD45去双细胞后.Rdata") 

# 加载R包
library(Seurat) 
library(circlize)
library(dplyr)
library(tibble)

# 暂时更改细胞名称的水平信息
Idents(seurat_object) <- factor(Idents(seurat_object),
                                levels = c("NK","T", "PC", "B","Neu", "Mono", "DC", "Mac"))
seurat_object[["celltype"]] <- Idents(seurat_object) 

# 我们的数据大群没有详细分亚群 所以这里直接将开始得seurat_clusers当作亚群
seurat_object$sub_type <- paste0("C",seurat_object$seurat_clusters)

df_data <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == "H1",] %>%
  distinct(celltype, sub_type, seurat_clusters) %>%
  dplyr::arrange(seurat_clusters) 

# 展示每个亚群的数量 当然直接展示数量柱状图太长 可以适用百分比
df_data$cellnumber <- prop.table(table(seurat_object$seurat_clusters)) * 100

# 按照plot顺序排序
df_data$celltype <- factor(df_data$celltype, levels = c("NK","T", "PC", "B","Neu", "Mono", "DC", "Mac"))
df_data <- df_data[order(df_data$celltype), ]

# 亚群名称 这里因为没有 所以自己结合了celltype和seurat_clusers
# 给细胞名称赋予全称
df_data$celltype1 <- c("NK",
                       "T","T","T","T","T","T","T","T",
                       "PC",
                       "B",
                       "Neu",
                       "Mono",
                       "DC","DC",
                       "Mac","Mac","Mac","Mac","Mac","Mac","Mac","Mac"
                       )

df_data$num <- seq(1:nrow(df_data))
df_data$cell_name <- paste0(df_data$sub_type,": ",df_data$celltype1 )

table(df_data$celltype)

group_colors <- c("NK" = "#f98177",
                  "T" = "#ae7eb5",
                  "PC"  = "#568389",
                  "B" = "#aad393",
                  "Neu"  = "#5a6fb5",
                  "Mono" = "#eeb066",
                  "DC" = "#54c6a8",
                  "Mac" = "#f082a5"
                  ) 

Seurat::DimPlot(seurat_object, reduction = "umap", pt.size = 0.1, label = T, cols = group_colors)

df_data$color <- group_colors[df_data$celltype]

DefaultAssay(seurat_object) <- "RNA"
markers <- c("KLRB1",     # NK
             "CD3E",      # T
             "MZB1",      # PC
             "MS4A1",     # B
             "CSF3R",     # Neu
             "FCN1",      # Mono
             "FLT3",      # DC
             "C1QB")      # Mac
             
Seurat::DotPlot(seurat_object, features = markers) + coord_flip()

circos.clear() # 清空当前作图 便于新的circle plot
group_size <- table(df_data$celltype) # 这个是每个细胞大群也就是分组的size 这里就是包含的亚群数目 需要注意这个涉及到后面扇形分区 所以顺序要对
group_size

# 设置布局
circos.par(start.degree = 270, cell.padding = c(0, 0, 0, 0), # 开口位置 扇区内行距为0
           gap.after = c(rep(2, length(group_size)-1),15),   # 设置每个扇区之间的gap 前面的扇区之间小一点 最后两个扇区也就是首尾的位置扇区开头大一点
           circle.margin = c(0.1, 0.1, 0.1, 0.1))            # 环形图距离画布的距离
# 初始化plot
circos.initialize(factors = df_data$celltype,  # 扇区scctor 这是已经排好序的数据
                  xlim = cbind(0, group_size)) # 每个扇区xlim 每个扇区元素不同 所以每个扇区的xlim是0到扇区元素长度

# 第一轨道
# 添加最外层每个细胞组的亚群
circos.track(
  ylim = c(0, 1), # y轴范围
  bg.border = NA, # 不要背景
  track.height = 0.01, # 高度
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index") # 获取当前扇区index
    group_size = group_size[sector_index]             # 获取当前扇区长度
    
    # 适用循环plot文字 因为是多个扇区
    for (i in 1:group_size) {
      circos.text(
        x = i - 0.5, # 位于中间
        y = 0.3,     # y轴位置
        labels = df_data$cell_name[df_data$celltype == sector_index][i], # 标注 索引到扇区对应的亚群
        col= df_data$color[df_data$celltype == sector_index][i],         # 颜色也是设置好的 索引即可
        
        font = 2, # 文字加粗
        facing = "reverse.clockwise", # 文字排布方式向外
        niceFacing = TRUE,
        adj = c(1, 0.5), 
        cex = 1) # 文字大小
    }
  }
)

# 第二轨道
# 因为我们是按照扇区来plot的 所以添加group注释就很简单了
circos.track(ylim = c(0, 1),
             bg.border = NA, 
             track.height = 0.06,
             bg.col=group_colors, # 分组注释背景颜色
             
             panel.fun=function(x, y) {

               xlim = get.cell.meta.data("xlim") # 获取当前扇区xy范围
               ylim = get.cell.meta.data("ylim")
               
               sector.index = get.cell.meta.data("sector.index") # 扇区索引
               circos.text(mean(xlim), # 取mean 文字位于中心
                           mean(ylim),
                           sector.index, # label就是扇区索引
                           col = "white", # 文字颜色
                           cex = 0.9, 
                           font = 2,
                           facing = 'bending.inside', 
                           niceFacing = TRUE)
               }
             )

# 第三轨道
# 确保celltype顺序和plot顺序一致
seurat_object$sub_type <- factor(seurat_object$sub_type, levels = df_data$sub_type)
seurat_object$group <- factor(seurat_object$group, levels = c("Health","Ketosis"))

# 将每组celltype数量添加到plot data
df_data$Health <- prop.table(table(seurat_object$group,seurat_object$sub_type),2)[1,] 
df_data$Ketosis <- prop.table(table(seurat_object$group,seurat_object$sub_type),2)[2,]

# 设置sample颜色，分组颜色
sample <- c("Health", "Ketosis")
sample_cols <- c("#82cffa", "#f27d92") %>% setNames(., sample)

# plot堆叠图 因为是面 使用circos.rect
# 和第一轨道plot一样 这里还是按照扇区添加
circos.track(
  
  # 基本设置
  ylim = c(0,1),
  bg.border = NA, 
  track.height = 0.06,

  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")         # 扇区索引
    group_data = df_data[df_data$celltype == sector_index, ]  # 获取当前扇区数据
    
    for (i in 1:nrow(group_data)) {
      circos.rect(
        xleft = i - 0.9,  # x轴左侧 -0.9是为了让柱状图宽度向右移动
        xright = i - 0.1, # x轴右侧 -0.1是为了让柱状图宽度向左移动 这样不同细胞得柱状图避免紧挨在一起
        ybottom = c(0, cumsum(as.vector(group_data[i, sample[1:(length(sample)-1)]]))), # ybottom是一个向量 因为是堆叠得rect 有三组数据 所以要设置三个底部 分别是开头的0 然后依次是第一个组得最高处是第二组得底部 第一个组+第二组是第三个组得底部
        ytop = cumsum(as.vector(group_data[i, sample])),#三个面的顶部分别是，三个组的比例依次叠加
        col = sample_cols, #对应得颜色设置
        border = NA
      )
    }
  }
)

# 第四轨道
# 确保celltype顺序和plot顺序一致
seurat_object$sub_type <- factor(seurat_object$sub_type, levels = df_data$sub_type)
seurat_object$Phase <- factor(seurat_object$Phase , levels = c("G1","G2M","S"))

# 将每组celltype 数量添加到plot data
df_data$G1 <- prop.table(table(seurat_object$Phase,seurat_object$sub_type),2)[1,] 
df_data$G2M <- prop.table(table(seurat_object$Phase,seurat_object$sub_type),2)[2,]
df_data$S <- prop.table(table(seurat_object$Phase,seurat_object$sub_type),2)[3,]

# 设置sample颜色 分组颜色
cyl <- c("G1","G2M","S")
cyl_cols <- c('#c7c1ff', '#8888ff', '#6b58ef') %>% setNames(., cyl)

# plot堆叠图 因为是面 使用circos.rect
# 和第一轨道plot一样 这里还是按照扇区添加
circos.track(
  
  # 基本设置
  ylim = c(0,1),
  bg.border = NA, 
  track.height = 0.06,
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")         # 扇区索引
    group_data = df_data[df_data$celltype == sector_index, ]  # 获取当前扇区数据
    
    for (i in 1:nrow(group_data)) {
      circos.rect(
        xleft = i - 0.9,  # x轴左侧，-0.9是为了让柱状图宽度向右移动
        xright = i - 0.1, # x轴右侧，-0.1是为了让柱状图宽度向左移动，这样不同细胞得柱状图避免紧挨在一起
        ybottom = c(0, cumsum(as.vector(group_data[i, cyl[1:(length(cyl)-1)]]))), # ybottom是一个向量 因为是堆叠得rect 有三组数据 所以要设置三个底部 分别是开头的0 然后依次是第一个组得最高处是第二组得底部 第一个组+第二组是第三个组得底部
        ytop = cumsum(as.vector(group_data[i, cyl])), # 三个面的顶部分别是 三个组的比例依次叠加
        col = cyl_cols, # 对应得颜色设置
        border = NA
      )
    }
  }
)

# 第五轨道 
# marker gene dot expression
DefaultAssay(seurat_object) <- "RNA"
markers <- c("KLRB1",     # NK
             "CD3E",       # T
             "MZB1",      # PC
             "MS4A1",     # B
             "CSF3R",     # Neu
             "FCN1",      # Mono
             "FLT3",      # DC
             "C1QB"      # Mac
             )
Idents(seurat_object) <- "sub_type"

# 首先使用seurat DotPlot自带函数计算marker表达 目的是为获取每个基因表达量和表达百分比
# 此外 需要呈现什么颜色就在此处设置好 便于后面提取颜色
p  = Seurat::DotPlot(seurat_object, features = markers, cols = c("lightgrey", "red"))

# 提取表达数据
dot_data <- p$data
colnames(dot_data)[4] <- 'sub_type'
dot_data$ypos <- rep(seq(1:length(markers)),23) # 23为亚群数目

# 合并表达量数据与分组数据
dot_data <- merge(dot_data, df_data, by = 'sub_type')
dot_data$sub_type <- factor(dot_data$sub_type, levels = unique(p$data$id)) # 设置数据的levels排序和plot点图一致 便于后面添加颜色
dot_data <- dot_data[order(dot_data$sub_type), ]

# 提取颜色
p_color <-  ggplot_build(p)$data[[1]]

# 将颜色添加到表达量数据
dot_data$dot_color <- p_color$colour

# 构建segment data
seg_data <- data.frame("cell_type" = c("NK",
                                       "T",
                                       "PC",
                                       "B",
                                       "Neu",
                                       "Mono",
                                       "DC",
                                       "Mac"
                                       ),
                       "marker_genes" = markers,
                       "seg_y" = 1:length(markers))

# 添加轨道并绘制点图
circos.track(
  ylim = c(0, length(markers)+1), # y轴范围，根据marker数量确定
  bg.border = "black",            # 扇区边界颜色
  track.height = 0.17,            # 轨道高度（在circos中是每一圈的宽度）
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")
    
    # 提取数据 for segment
    seg_sub <- seg_data[seg_data$cell_type == sector_index,]
    
    # 提取数据 for dotplot
    cell_subtypes = unique(dot_data$sub_type[dot_data$celltype == sector_index])
    
    # 标注segment线
    circos.segments(x0 = 0, 
                    y0 = seg_sub$seg_y, 
                    x1 = length(cell_subtypes), 
                    y1 = seg_sub$seg_y, 
                    lty=1, lwd=0.15)
    
    # 绘制点图
    for (i in 1:length(cell_subtypes)) {
      
      subtype_data <- dot_data[dot_data$sub_type == cell_subtypes[i], ]
      
      circos.points(
        x = i - 0.5, 
        y = subtype_data$ypos,
        pch = 16,
        col = subtype_data$dot_color,
        cex = subtype_data$pct.exp / 30  # 缩放除以30
      )}
    
    if(length(seg_sub$seg_y) >1){
      
      y1 <- c(seg_sub$seg_y[min(1:length(seg_sub$seg_y))]-1,seg_sub$seg_y[max(1:length(seg_sub$seg_y))]+1)
      
    }else{
      
      y1 <- seg_sub$seg_y
    }
    
    circos.text(x = length(cell_subtypes)/2,
                y = y1, 
                font = 3, 
                cex = 0.8, 
                labels = seg_sub$marker_genes, 
                col="black",
                facing="bending.inside", 
                niceFacing = TRUE) 
    }
  )

# 第六轨道 
circos.track(
  
  # 基本设置
  ylim = c(0, max(df_data$cellnumber)),
  bg.border = NA, 
  track.height = 0.1,
  
  # 标注segment线

  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")             # 扇区索引
    cellnumber_data = df_data[df_data$celltype == sector_index, ] # 获取当前扇区数据
    
    # 划分割celltype种类的虚线 由于我们的plot是按照扇区plot的 所以这个可有可无
    circos.segments(x0 = nrow(cellnumber_data), 
                    y0 = 0, 
                    x1 = nrow(cellnumber_data), 
                    y1 = 100, 
                    lty=2, lwd=0.35)
    
    
    for (i in 1:nrow(cellnumber_data)) {
      circos.rect(
        xleft = i - 0.8,  # x轴左侧 -0.9是为了让柱状图宽度向右移动
        xright = i - 0.2, # x轴右侧 -0.1是为了让柱状图宽度向左移动 这样不同细胞得柱状图避免紧挨在一起
        ybottom =0,
        ytop = cellnumber_data$cellnumber[i], # 三个面的顶部分别是 三个组的比例依次叠加
        col = "#e3782b", # 对应的颜色设置
        border = NA
      )
    }
  }
)

# 画中心的UMAP图
group_colors <- c("T" = "#ae7eb5",
                  "NK" = "#f98177",
                  "B" = "#aad393",
                  "PC"  = "#568389",
                  "Mono" = "#eeb066",
                  "Mac" = "#f082a5",
                  "DC" = "#54c6a8",
                  "Neu"  = "#5a6fb5"
                  ) 

Idents(seurat_object) <- 'celltype'
Seurat::DimPlot(seurat_object, 
        label = T, 
        label.size = 7, 
        cols = group_colors, 
        pt.size = 0.5)  + 
  theme_void() + NoLegend()

# 四个角为分类型的UMAP 目的cluster为彩色 其余细胞为灰色即可
# T/NK lineage
cluster_colors <- c("0"  = "grey75","1"  = "grey75","2"  = "#bef0b0","3"  = "#ff6f94","4"  = "grey75",
                    "5"  = "grey75","6"  = "#ffc761","7"  = "grey75","8"  = "grey75","9"  = "grey75",
                    "10"  = "grey75","11"  = "grey75","12"  = "#2990c0",  "13"  = "grey75","14" = "grey75",
                    "15"  = "grey75","16"  = "#ff9572","17"  = "#d65db1","18"  = "grey75",  "19"  = "grey75",
                    "20" = "#7ecdbb",  "21" = "#6f99ad", "22" = "#cab2d6")
Seurat::DimPlot(seurat_object,group.by = "seurat_clusters",cols = cluster_colors,pt.size = 1.5,raster = TRUE) + 
  theme_void() + NoLegend()+ labs(title = NULL)

# B-lineage 
cluster_colors <- c("0"  = "grey75","1"  = "grey75","2"  = "grey75","3"  = "grey75","4"  = "#00bbb1",
                    "5"  = "grey75","6"  = "grey75","7"  = "grey75", "8"  = "grey75", "9"  = "grey75",
                    "10"  = "grey75","11"  = "#f9b64b","12"  = "grey75",  "13"  = "grey75","14" = "grey75",
                    "15"  = "grey75","16"  = "grey75","17"  = "grey75","18"  = "grey75",  "19"  = "grey75",
                    "20" = "grey75", "21" = "grey75", "22" = "grey75")
Seurat::DimPlot(seurat_object,group.by = "seurat_clusters",cols = cluster_colors,pt.size = 1.5,raster = TRUE) + 
  theme_void() + NoLegend()+ labs(title = NULL)

# MPS 
cluster_colors <- c("0"  = "#e4c595","1"  = "#e88ac3","2"  = "grey75","3"  = "grey75","4"  = "grey75",
                    "5"  = "#abdda4","6"  = "grey75","7"  = "grey75","8"  = "#a1d955","9"  = "#8f9fcb",
                    "10"  = "#bebada","11"  = "grey75","12"  = "grey75",  "13"  = "#52b2d3","14" = "#fb8f64",
                    "15"  = "#c93103","16"  = "grey75","17"  = "grey75","18"  = "#64c2a5",  "19"  = "#fcc652",
                    "20" = "grey75", "21" = "grey75", "22" = "grey75")
Seurat::DimPlot(seurat_object,group.by = "seurat_clusters",cols = cluster_colors,pt.size = 1.5,raster = TRUE) + 
  theme_void() + NoLegend()+ labs(title = NULL)

# Granulocytes
cluster_colors <- c("0"  = "grey75","1"  = "grey75","2"  = "grey75","3"  = "grey75","4"  = "grey75",
                    "5"  = "grey75","6"  = "grey75","7"  = "#6471b8","8"  = "grey75","9"  = "grey75",
                    "10"  = "grey75","11"  = "grey75","12"  = "grey75",  "13"  = "grey75","14" = "grey75",
                    "15"  = "grey75","16"  = "grey75","17"  = "grey75","18"  = "grey75", "19"  = "grey75",
                    "20" = "grey75", "21" = "grey75", "22" = "grey75")
Seurat::DimPlot(seurat_object,group.by = "seurat_clusters",cols = cluster_colors,pt.size = 1.5,raster = TRUE) + 
  theme_void() + NoLegend()+ labs(title = NULL)

# =============================================
# 🎨 Fig. 1d
# =============================================

# Marker基因的气泡图
# 对注释后细胞亚群 进行差异基因分析 前面已完成 1.2倍 20%占比
# 查看合适的Marker画图
load("CD45去双细胞后.Rdata")
load("CD45整体细胞类型差异基因.Rdata")
head(diff[diff$cluster == unique(diff$cluster)[8],], 50) $ gene 
library(scplotter)
library(plotthis)
plotthis::show_palettes(type = "continuous", index = 1:30)
Seurat::DotPlot(seurat_object, 
                features = c("CD3D","CD3E","CD3G", "BCL11B",   # T
                             "NKG7","KLRD1","NCR1","KLRF1",    # NK
                             "MS4A1","CD79B","CD79A","PAX5",   # B
                             "JCHAIN","MZB1","TXNDC5","DERL3", # PC
                             "FCN1","VCAN","S100A8","S100A9",  # Mono
                             "C1QA","C1QB","C1QC","CD68",      # Mφ
                             "CST3","PLD4","FLT3","P2RY6",     # DC
                             "CSF3R","IL1B", "CXCL8","CXCR2"   # Neu
                             ),           
                scale = T, # 还是得加上归一化才好看 
                group.by = "celltype") +
  # 旋转细胞名
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10)) +
  theme(axis.text.y = element_text(face = "italic", size = 10)) +
  # 去除横纵坐标标签
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  # 设置颜色渐变
  scale_color_distiller(palette = "RdBu")  +
  # 添加四周框线 并去掉XY轴线
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), # 四周框线
        panel.grid = element_blank(),    # 去掉背景网格线
        axis.line = element_blank()) +   # 去掉重复的XY轴线
  # 转置坐标轴
  coord_flip()
dev.off()


# =============================================
# 🎨 Fig. 1e
# =============================================

# Marker基因的UMAP图
load("CD45去双细胞后.Rdata")
library(scplotter)
# T
FeatureStatPlot(seurat_object, plot_type = "dim", features = "CD3E", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())   
# NK
FeatureStatPlot(seurat_object, plot_type = "dim", features = "KLRF1", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   
# B
FeatureStatPlot(seurat_object, plot_type = "dim", features = "MS4A1", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   
# PC
FeatureStatPlot(seurat_object, plot_type = "dim", features = "MZB1", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   
# Mono
FeatureStatPlot(seurat_object, plot_type = "dim", features = "VCAN", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   
# MAC
FeatureStatPlot(seurat_object, plot_type = "dim", features = "C1QB", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   
# DC
FeatureStatPlot(seurat_object, plot_type = "dim", features = "FLT3", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   
# Neu
FeatureStatPlot(seurat_object, plot_type = "dim", features = "CXCR2", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   
# 图注
FeatureStatPlot(seurat_object, plot_type = "dim", features = "CXCR2", reduction = "umap", lower_cutoff = 0.5, upper_cutoff = 2) +
  theme(legend.position = "top",axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())   


# =============================================
# 🎨 Fig. 1g
# =============================================
# Marker基因棒棒糖图
load("CD45整体细胞类型差异基因.Rdata")
table(diff$cluster) 
diff = diff[diff$p_val < 0.05,]
table(diff$cluster)
deg_df <- data.frame(
  group = factor(c("T", "NK", "B", "PC", "Mono", "Mac", "DC", "Neu"),
                 levels = c("Neu", "DC", "Mac", "Mono", "PC", "B", "NK", "T")),
  n_deg = c(877, 878, 984, 885, 1869, 700, 1974, 533)
  )

group_colors <- c(
  "T"    = "#ae7eb5",
  "NK"   = "#f98177",
  "B"    = "#aad393",
  "PC"   = "#568389",
  "Mono" = "#eeb066",
  "Mac"  = "#f082a5",
  "DC"   = "#54c6a8",
  "Neu"  = "#5a6fb5"
  )

ggplot(deg_df, aes(x = n_deg, y = group)) +
  geom_segment(aes(x = 0, xend = n_deg, y = group, yend = group),
               color = "grey80", linewidth = 1) +
  geom_point(aes(fill = group), shape = 21, size = 5, color = "black", stroke = 1) +
  scale_fill_manual(values = group_colors) +
  labs(x = "Number of DEGs", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(size = 12)
    ) +
  scale_x_continuous(
    limits = c(0, 2250),
    breaks = c(0, 1000, 2000),
    expand = c(0, 0)
    )

# Marker基因的GO通路 只需要一个最典型的即可
load("CD45整体细胞类型差异基因.Rdata")
table(diff$cluster)
# 将diff数据框按cluster进行拆分为单个数据框
dif_list <- split(diff, diff$cluster)
names(dif_list)
T = dif_list[[1]]
NK = dif_list[[2]]
B = dif_list[[3]]
PC = dif_list[[4]]
Mono = dif_list[[5]]
Mac = dif_list[[6]]
DC = dif_list[[7]]
Neu = dif_list[[8]]

# 分别进行功能分析
GSEA_T <- SCP::RunGSEA(
  geneID = T$gene,                    
  geneScore = T$avg_log2FC,           
  geneID_groups = rep("a",times=877), 
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_T = GSEA_T[["enrichment"]]

GSEA_NK <- SCP::RunGSEA(
  geneID = NK$gene,                    
  geneScore = NK$avg_log2FC,          
  geneID_groups = rep("a",times=878), 
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_NK = GSEA_NK[["enrichment"]]

GSEA_B <- SCP::RunGSEA(
  geneID = B$gene,                    
  geneScore = B$avg_log2FC,           
  geneID_groups = rep("a",times=984),
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_B = GSEA_B[["enrichment"]]

GSEA_PC <- SCP::RunGSEA(
  geneID = PC$gene,             
  geneScore = PC$avg_log2FC,         
  geneID_groups = rep("a",times=885), 
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_PC = GSEA_PC[["enrichment"]]

GSEA_Mono <- SCP::RunGSEA(
  geneID = Mono$gene,            
  geneScore = Mono$avg_log2FC,          
  geneID_groups = rep("a",times=1869), 
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_Mono = GSEA_Mono[["enrichment"]]

GSEA_Mac <- SCP::RunGSEA(
  geneID = Mac$gene,                   
  geneScore = Mac$avg_log2FC,           
  geneID_groups = rep("a",times=700), 
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_Mac = GSEA_Mac[["enrichment"]]

GSEA_DC <- SCP::RunGSEA(
  geneID = DC$gene,            
  geneScore = DC$avg_log2FC,           
  geneID_groups = rep("a",times=1974), 
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_DC = GSEA_DC[["enrichment"]]

GSEA_Neu <- SCP::RunGSEA(
  geneID = Neu$gene,          
  geneScore = Neu$avg_log2FC,       
  geneID_groups = rep("a",times=533), 
  db = c("GO_BP"),
  species = "Bos_taurus"
  )
GSEA_Neu = GSEA_Neu[["enrichment"]]

# 保存数据
save(GSEA_T, file = 'GSEA_T.Rdata') 
save(GSEA_NK, file = 'GSEA_NK.Rdata')
save(GSEA_B, file = 'GSEA_B.Rdata')
save(GSEA_PC, file = 'GSEA_PC.Rdata')
save(GSEA_Mono, file = 'GSEA_Mono.Rdata')
save(GSEA_Mac, file = 'GSEA_Mac.Rdata')
save(GSEA_DC, file = 'GSEA_DC.Rdata')
save(GSEA_Neu, file = 'GSEA_Neu.Rdata')

# 选择最符合细胞特征的TOP通路
GSEA_T$Description[1:20] # 8 T cell activation
GSEA_NK$Description[1:20] # 1 immune effector process
GSEA_B$Description[1:20] # 1 cell surface receptor signaling pathway  
GSEA_PC$Description[1:20] # 6 endoplasmic reticulum unfolded protein response
GSEA_Mono$Description[1:20] # 16 response to lipopolysaccharide
GSEA_Mac$Description[1:20] # 2 lipid metabolic process
GSEA_DC$Description[1:20] # 2 antigen processing and presentation
GSEA_Neu$Description[1:20] # 2 chemotaxis

GSEA_T[8,]
GSEA_NK[1,]
GSEA_B[1,]
GSEA_PC[6,]
GSEA_Mono[16,]
GSEA_Mac[2,]
GSEA_DC[2,]
GSEA_Neu[2,]

gsea_plot_df <- dplyr::bind_rows(
  dplyr::mutate(as.data.frame(GSEA_T)[8, ], group = "T"),
  dplyr::mutate(as.data.frame(GSEA_NK)[1, ], group = "NK"),
  dplyr::mutate(as.data.frame(GSEA_B)[1, ], group = "B"),
  dplyr::mutate(as.data.frame(GSEA_PC)[6, ], group = "PC"),
  dplyr::mutate(as.data.frame(GSEA_Mono)[16, ], group = "Mono"),
  dplyr::mutate(as.data.frame(GSEA_Mac)[2, ], group = "Mac"),
  dplyr::mutate(as.data.frame(GSEA_DC)[2, ], group = "DC"),
  dplyr::mutate(as.data.frame(GSEA_Neu)[2, ], group = "Neu")
  )

colnames(gsea_plot_df)

gsea_plot_df$group <- factor(
  gsea_plot_df$group,
  levels = c("Neu", "DC", "Mac", "Mono", "PC", "B", "NK", "T")
  )

gsea_plot_df <- gsea_plot_df[order(gsea_plot_df$group), ]
gsea_plot_df$Description <- factor(gsea_plot_df$Description, levels = gsea_plot_df$Description)

group_colors <- c(
  "T"    = "#ae7eb5",
  "NK"   = "#f98177",
  "B"    = "#aad393",
  "PC"   = "#568389",
  "Mono" = "#eeb066",
  "Mac"  = "#f082a5",
  "DC"   = "#54c6a8",
  "Neu"  = "#5a6fb5")

ggplot(gsea_plot_df, aes(x = NES, y = Description, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(x = "NES", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(size = 12)
    )+
  scale_x_continuous(breaks = c(1, 1.3, 1.6, 1.9, 2.2)) +
  coord_cartesian(xlim = c(1, 2.3)) 

# =============================================
# 🎨 Supplementary Figure 3c
# =============================================

# 对每个亚群的1.2倍显著上调基因进行GO的GSEA富集分析 对TOP通路画气泡图
# 前面已经做了每个亚群的1.2倍显著上调基因进行GO的GSEA富集分析 这里直接加载用
load("~/奶牛肝脏解离单细胞3版/GSEA_T.Rdata")
load("~/奶牛肝脏解离单细胞3版/GSEA_PC.Rdata")
load("~/奶牛肝脏解离单细胞3版/GSEA_NK.Rdata")
load("~/奶牛肝脏解离单细胞3版/GSEA_Neu.Rdata")
load("~/奶牛肝脏解离单细胞3版/GSEA_Mono.Rdata")
load("~/奶牛肝脏解离单细胞3版/GSEA_Mac.Rdata")
load("~/奶牛肝脏解离单细胞3版/GSEA_DC.Rdata")
load("~/奶牛肝脏解离单细胞3版/GSEA_B.Rdata")
# 重新命名
T = GSEA_T
NK = GSEA_NK
B = GSEA_B
PC = GSEA_PC
Mono = GSEA_Mono
Mac = GSEA_Mac
DC = GSEA_DC
Neu = GSEA_Neu
# 选择每个亚群的TOP通路
T <- T[1:20,]
NK <- NK[1:20,]
B <- B[1:20,]
PC <- PC[1:20,]
Mono <- Mono[1:20,]
Mac <- Mac[1:20,]
DC <- DC[1:20,]
Neu <- Neu[1:20,]
# 为每个细胞群体添加标签
T$group <- "T"
NK$group <- "NK"
B$group <- "B"
PC$group <- "PC"
Mono$group <- "Mono"
Mac$group <- "Mac"
DC$group <- "DC"
Neu$group <- "Neu"
# 选择TOP通路
# T
select_T = c("T cell activation",
             "regulation of response to stimulus",
             "positive regulation of response to stimulus",
             "immune system process",
             "lymphocyte activation",
             "immune response",
             "regulation of immune system process",
             "cell surface receptor signaling pathway")
# NK
select_NK = c("immune effector process",
              "innate immune response",
              "response to biotic stimulus",
              "defense response to other organism",
              "response to other organism",
              "defense response to symbiont")
# B
select_B = c("immune response-regulating signaling pathway",
              "regulation of immune response",
              "activation of immune response",
             "anatomical structure development",
             "developmental process",
             "positive regulation of immune response",
             "regulation of immune system process")
# PC
select_PC = c("endoplasmic reticulum unfolded protein response",
               "cellular response to unfolded protein",
               "response to unfolded protein",
               "response to endoplasmic reticulum stress",
              "protein localization to endoplasmic reticulum")
# Mono
select_Mono = c("response to lipopolysaccharide",
                "defense response to bacterium",
                "inflammatory response",
                "response to bacterium",
                "response to stress",
                "response to molecule of bacterial origin")
# Mac
select_Mac = c("lipid metabolic process",
               "cell population proliferation",
               "cellular lipid metabolic process",
               "regulation of leukocyte cell-cell adhesion",
               "mononuclear cell proliferation",
               "leukocyte proliferation")
# DC
select_DC = c("antigen processing and presentation", 
              "production of molecular mediator of immune response",
              "antigen processing and presentation of peptide antigen via MHC class II",
              "antigen processing and presentation of peptide or polysaccharide antigen via MHC class II", 
              "antigen processing and presentation of peptide antigen")
# Neu
select_Neu = c("chemotaxis",
               "leukocyte migration", 
               "G protein-coupled receptor signaling pathway",
               "cellular response to chemical stimulus",
               "response to lipid",
               "taxis",
               "cell migration")

select = c(select_T, select_NK, select_B, select_PC,
           select_Mono, select_Mac, select_DC, select_Neu)
select

# 选择每个亚群的TOP通路
T <- T[T$Description %in% select,]
NK <- NK[NK$Description %in% select,]
B <- B[B$Description %in% select,]
PC <- PC[PC$Description %in% select,]
Mono <- Mono[Mono$Description %in% select,]
Mac <- Mac[Mac$Description %in% select,]
DC <- DC[DC$Description %in% select,]
Neu <- Neu[Neu$Description %in% select,]
# 生成新的P值列
T$`-log10pvalue` <- -log10(T$pvalue)
NK$`-log10pvalue` <- -log10(NK$pvalue)
B$`-log10pvalue` <- -log10(B$pvalue)
PC$`-log10pvalue` <- -log10(PC$pvalue)
Mono$`-log10pvalue` <- -log10(Mono$pvalue)
Mac$`-log10pvalue` <- -log10(Mac$pvalue)
DC$`-log10pvalue` <- -log10(DC$pvalue)
Neu$`-log10pvalue` <- -log10(Neu$pvalue)
# 合并所有数据
all <- rbind(T, NK, B, PC, Mono, Mac, DC, Neu)
all$Description <- gsub("-", " ", all$Description)
# 为Description列设置因子顺序
library(forcats)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)
# 指定分组顺序
My_levels <- c("T", "NK", "B", "PC", "Mono", "Mac", "DC", "Neu")
all$group <- factor(all$group, levels= My_levels)
# 绘制GO气泡图
ggplot(all, aes(group, Description)) +
  theme_bw() +
  geom_point(aes(fill = `-log10pvalue`, size = NES), shape = 21, colour = "black", alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_line(color = "#ececec", size = 0.5),
        panel.grid.minor = element_line(color = "#ececec", size = 0.5),
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 1)) +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_size_continuous(range = c(1, 5))


# =============================================
# 🎨 Fig. 1f
# =============================================

# TOP10 Marker基因热图
library(Seurat)
library(pheatmap)
library(ggplot2)
library(dplyr)

# 加载数据
load("CD45去双细胞后.Rdata")
load("CD45整体细胞类型差异基因.Rdata")
table(seurat_object@meta.data[["celltype"]])

# 随机抽样
table(seurat_object@meta.data[["celltype"]])
cell_types <- unique(seurat_object$celltype)

# 创建一个空的 list 来存储每个细胞类型的子集
subset_cells <- list()

# 对每个细胞类型进行处理
for (cell_type in cell_types) {
  cells_of_type <- WhichCells(seurat_object, expression = celltype == cell_type) 
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
new_seurat_object <- subset(seurat_object, cells = subset_all_cells)

# 查看新 Seurat对象
new_seurat_object

# TOP10基因
dif<-FindAllMarkers(seurat_object,
                    group.by = seurat_object@active.ident, 
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
# 🎨 Supplementary Figure 2c, d
# =============================================
  
# 质控前后信息对比
# 加载数据 质控后
load("CD45去双细胞后.Rdata")
VlnPlot(seurat_object, 
        features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI", "percent.mt", "percent.HB"), 
        ncol = 6,
        group.by = "orig.ident"
        ) 
# 加载数据 质控前
load("NCBI原始数据.Rdata")
seurat_object2 = seurat_object
VlnPlot(seurat_object, 
        features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI", "percent.mt", "percent.HB"), 
        ncol = 6,
        group.by = "orig.ident"
        ) 
  
seurat_object$orig.ident <- factor(
  seurat_object$orig.ident,
  levels = c("H1", "H2", "H3", "K1", "K2", "K3", "K4", "K5"),
  labels = c("Health1", "Health2", "Health3", "Ketosis1", "Ketosis2", "Ketosis3", "Ketosis4", "Ketosis5")
  )
seurat_object2$orig.ident <- factor(
  seurat_object2$orig.ident,
  levels = c("H1", "H2", "H3", "K1", "K2", "K3", "K4", "K5"),
  labels = c("Health1", "Health2", "Health3", "Ketosis1", "Ketosis2", "Ketosis3", "Ketosis4", "Ketosis5")
  )
table(seurat_object$orig.ident)
table(seurat_object2$orig.ident)

VlnPlot(
  seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI", "percent.mt", "percent.HB"),
  ncol = 5,
  group.by = "orig.ident",
  cols = c("#8ad0fa", "#8ad0fa", "#8ad0fa", "#f38297",
           "#f38297", "#f38297", "#f38297", "#f38297"),
  alpha = 0.2) &
  theme(
    axis.title.x = element_blank(),   # 去掉横坐标标题 Identity
    plot.title = element_blank(), # 去掉每个小图标题
    axis.text.x = element_text(size = 14.5),
    axis.text.y = element_text(size = 14.5)
  )

VlnPlot(
  seurat_object2,
  raster = FALSE,
  features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI", "percent.mt", "percent.HB"),
  ncol = 5,
  group.by = "orig.ident",
  cols = c("#8ad0fa", "#8ad0fa", "#8ad0fa", "#f38297",
           "#f38297", "#f38297", "#f38297", "#f38297"),
  alpha = 0.2) &
  theme(
    axis.title.x = element_blank(),   # 去掉横坐标标题 Identity
    plot.title = element_blank(), # 去掉每个小图标题
    axis.text.x = element_text(size = 14.5),
    axis.text.y = element_text(size = 14.5)
  )

# =============================================
# 🎨 Supplementary Figure 3b
# =============================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(corrplot)
library(RColorBrewer)

load("CD45去双细胞后.Rdata")
seurat_object$subtype <- paste0("Cluster_",seurat_object$seurat_clusters)
Idents(seurat_object) <- "subtype"
Markers <- FindAllMarkers(seurat_object, 
                          only.pos= T, 
                          min.pct = 0.2, 
                          logfc.threshold = log2(1.2), 
                          verbose = T, 
                          max.cells.per.ident = 25)

# 将一些有影响的基因去除 比如核糖体基因
Markers <- subset(Markers[grep("^RP[L|S]", Markers$gene, ignore.case = FALSE, invert=TRUE),], subset = p_val_adj < 0.05)
# 计算基因平均表达量
Markers_av <- AverageExpression(seurat_object,group.by = "subtype",features = unique(Markers$gene),assays = "RNA") 
Markers_av <- Markers_av$RNA
gene_cell_exp <- t(scale(t(Markers_av),scale = T,center = T))
# 计算相关性
cell_cor <- cor(as.matrix(gene_cell_exp), method = 'spearman')
# 显著性检验
cell_cor_p <- cor.mtest(cell_cor, conf.level = 0.95)$p

# 获取细胞亚群与大群对应关系
df <- seurat_object@meta.data[,c("celltype","subtype")] %>%
  distinct(celltype, subtype)
rownames(df) <- df$subtype

# 设置分组颜色
group_color <- c(
  "T"    = "#ae7eb5",
  "NK"   = "#f98177",
  "B"    = "#aad393",
  "PC"   = "#568389",
  "Mono" = "#eeb066",
  "Mac"  = "#f082a5",
  "DC"   = "#54c6a8",
  "Neu"  = "#5a6fb5")
df$color <- group_color[df$celltype]

lab_order <- colnames(p$corr)
lab_order 
lab_order <- gsub("-", "_", lab_order)
lab_order 
df <- df[lab_order,] 

corrplot(cell_cor, 
         method = 'circle',
         order = "hclust", 
         hclust.method='ward.D',
         type = "full", 
         addrect = 5,
         col = rev(brewer.pal(n=8, name="RdYlBu")),
         tl.cex=0.8,
         tl.col = df$color,
         cl.cex = 0.5,
         cl.pos = "b",
         cl.length = 5,
         cl.ratio=0.1,
         sig.level = 0.05,
         insig='blank',
         mar=c(3, 0, 2, 0)) 


df$subtype <- factor(df$subtype, levels = df$subtype) # 固定顺序

ggplot(df, aes(x=subtype,y=1,fill=celltype))+ 
  geom_tile() + 
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom")+
  scale_fill_manual(values = group_color)


# =============================================
# 🎨 Supplementary Figure 3a
# =============================================

# 使用基因表达数据绘制细胞亚群之间的相关性聚类树
library(ggdendro)
Markers <- FindAllMarkers(seurat_object,
                          only.pos = T, 
                          min.pct = 0.2, 
                          logfc.threshold = log2(1.2), 
                          verbose = T,
                          max.cells.per.ident = 2000
                          )
# 筛选差异表达基因: p值小于0.05且不以"RP"开头 "RP"为核糖体基因 常常是高表达但无生物学意义的基因
Markers <- subset(Markers[grep("^RP[L|S]", 
                               Markers$gene, 
                               ignore.case = FALSE, 
                               invert=TRUE),],
                  subset = p_val_adj < 0.05
                  )
# 计算每个基因在不同细胞类型中的平均表达
Markers_av <- AverageExpression(seurat_object,
                                group.by = "celltype",
                                features = unique(Markers$gene),
                                assays = "RNA"
                                )
# 提取RNA数据
Markers_av <- Markers_av$RNA
# 标准化表达数据
gene_cell_exp <- t(scale(t(Markers_av),scale = T,center = T))
# 层次聚类
hc1 = hclust(dist(t(as.matrix((gene_cell_exp)))),method = 'ward.D')
# 按细胞类型着色
hc1 = as.dendrogram(hc1)
labels_colors(hc1) <- c("T" = "#ae7eb5",
                        "NK" = "#f98177",
                        "B" = "#aad393",
                        "PC" = "#568389",
                        "Mono" = "#eeb066",
                        "Mac" = "#f082a5",
                        "DC" = "#54c6a8",
                        "Neu" = "#5a6fb5")[labels(hc1)]
# 根据层次聚类的结果绘制一个树状图 并在树状图的叶节点上添加标签
plot(dend, horiz = TRUE, main = "", xlab = "", col = "#8c6bb1")


# =============================================
# 🎨 Supplementary Figure 1b-m
# =============================================

# TG(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(0.82, 1.14, 1.59, 3.56, 4.14, 7.31, 6.76, 5.16) 
df <- data.frame(Group = group, Value = value)
set.seed(123)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 8.7) + 
  ylab("TG (% wet weight)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),  # 👈这里
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
t.test(Value ~ Group, data = df) # 0.00318 独立两组-正态-方差非齐性

# BW(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(638.82, 652.42, 655.18, 672.81, 649.18, 657.83, 641.37, 642.93)
df <- data.frame(Group = group, Value = value)
set.seed(123)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 679) +  
  ylab("BW (kg)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),  # 👈这里
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
t.test(Value ~ Group, data = df, var.equal = TRUE) # 0.6546 独立两组-正态-方差齐性

# Milk yield(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(31.72, 32.86, 34.09, 27.11, 28.37, 29.74, 26.21, 29.04)
df <- data.frame(Group = group, Value = value)
dev.off()
pdf("P43.pdf", width = 2.8, height = 3.3) # 新建画布
set.seed(123)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 35.8) +  # 设置 y 轴范围
  ylab("Milk yield (kg/d)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),  # 👈这里
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
t.test(Value ~ Group, data = df, var.equal = TRUE) # 0.002856 独立两组-正态-方差齐性

# DMI(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(21.37, 21.62, 22.47, 19.17, 19.56, 20.12, 18.82, 19.89)
df <- data.frame(Group = group, Value = value)
dev.off()
pdf("P44.pdf", width = 2.8, height = 3.3) # 新建画布
set.seed(123)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5, hjust = 0.5, vjust = 0.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 23.2) +  
  ylab("DMI (kg/d)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90),  
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
t.test(Value ~ Group, data = df, var.equal = TRUE) # 0.001142 独立两组-正态-方差齐性

# BCS(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(3.16, 3.29, 3.37, 3.39, 3.51, 3.58, 3.31, 3.63)
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 3.74) +  
  ylab("BCS") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),
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
dev.off()  
shapiro.test(df$Value[df$Group == "Health"])
shapiro.test(df$Value[df$Group == "Ketosis"])
leveneTest(Value ~ Group, data = df, center = "mean")
t.test(Value ~ Group, data = df, var.equal = TRUE) # 0.05947 独立两组-正态-方差齐性

# NEFA(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(0.26, 0.29, 0.37, 1.09, 1.13, 1.35, 1.01, 1.27)
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 1.59) +  
  ylab("Serum NEFA (mM)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),
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
dev.off()  
shapiro.test(df$Value[df$Group == "Health"])
shapiro.test(df$Value[df$Group == "Ketosis"]) 
leveneTest(Value ~ Group, data = df, center = "mean")
t.test(Value ~ Group, data = df, var.equal = TRUE) # 5.523e-05 独立两组-正态-方差齐性

# BHBA(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(0.58, 0.64, 0.73, 3.23, 3.62, 3.73, 3.06, 3.47)
df <- data.frame(Group = group, Value = value)
dev.off()
pdf("P47.pdf", width = 2.8, height = 3.3) 
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5, hjust = 0.5, vjust = 0.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 4.4) +
  ylab("Serum BHBA (mM)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90),
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
dev.off()  
shapiro.test(df$Value[df$Group == "Health"])
shapiro.test(df$Value[df$Group == "Ketosis"])
leveneTest(Value ~ Group, data = df, center = "mean") 
t.test(Value ~ Group, data = df, var.equal = TRUE) # 3.105e-06 独立两组-正态-方差齐性

# GLU(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(3.59, 3.74, 3.89, 2.24, 2.37, 2.51, 2.73, 2.54)
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 4.3) + 
  ylab("Serum Glucose (mM)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),
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
t.test(Value ~ Group, data = df, var.equal = TRUE) # 6.026e-05 独立两组-正态-方差齐性

# AST(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(63.96, 70.12, 72.38, 141.12, 151.87, 157.11, 138.18, 145.23)
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 175) +
  ylab("Serum AST (U/L)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
wilcox.test(Value ~ Group, data = df)
t.test(Value ~ Group, data = df, var.equal = TRUE) # 4.328e-06 独立两组-正态-方差齐性

# ALT(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(16.33, 18.01, 19.33, 68.88, 74.56, 89.07, 65.88, 78.81)
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 105) + 
  ylab("Serum ALT (U/L)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
dev.off()  
shapiro.test(df$Value[df$Group == "Health"])
shapiro.test(df$Value[df$Group == "Ketosis"]) 
leveneTest(Value ~ Group, data = df, center = "mean") 
t.test(Value ~ Group, data = df, var.equal = TRUE) # 4.342e-05 独立两组-正态-方差齐性

# 肝脏IKKβ活性(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(0.757587111,
           1.002622705,
           1.239790184,
           2.738104159,
           2.807793181,
           3.056200824,
           3.322592731,
           3.728362683
           )
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 4.35) +  
  ylab("Hepatic IKKβ activity (fold)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
t.test(Value ~ Group, data = df, var.equal = TRUE) # 0.0001871 独立两组-正态-方差齐性

# 肝脏Caspase-1活性(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(0.790984639,
           0.929992445,
           1.279022916,
           1.866784185,
           1.952908587,
           2.138000504,
           2.261143289,
           2.437169479
           )
df <- data.frame(Group = group, Value = value)
dev.off()
pdf("P50-2.pdf", width = 2.8, height = 3.3) 
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 2.88) + 
  ylab("Hepatic Caspase-1 activity (fold)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),
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
t.test(Value ~ Group, data = df, var.equal = TRUE) # 0.0006212 独立两组-正态-方差齐性


# =============================================
# 🎨 Supplementary Figure 1n
# =============================================

# 油红O(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(0.017, 0.008, 0.003, 0.201, 0.298, 0.377, 0.412, 0.445)
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 0.52) + 
  ylab("Oil red area/total area") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
t.test(Value ~ Group, data = df) # 0.001456 独立两组-正态-方差非齐性


# =============================================
# 🎨 Supplementary Figure 1o
# =============================================

# 脂肪细胞大小(H1/2/3 K1/2/3/4/5)
group <- c(rep("Health", 3), rep("Ketosis", 5))
value <- c(185, 197, 231, 161, 150, 127, 119, 101)
df <- data.frame(Group = group, Value = value)
set.seed(125)
ggplot(data = df, aes(x = Group, y = Value)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 6.5) +
  scale_fill_manual(values=c("#8ad0fa", "#f38297")) + 
  ylim(NA, 252) +  
  ylab("Adipocyte diameter (μm)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 19, color = "black"),  
        axis.text.y = element_text(size = 18, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
t.test(Value ~ Group, data = df, var.equal = TRUE) # 0.00602 独立两组-正态-方差齐性
