
# F2/S4: Clinical ketosis reshapes the hepatic immune cell composition in dairy cows.
# Author: Chenchen Zhao
# Date: 2026-06-01
# Contact: jluzhaocc@126.com


# =============================================
# 🎨 Fig. 2a
# =============================================

# 分组的细胞比例图 H1/2/3 K1/2/3/4/5
library(dplyr)
library(ggplot2)
library(gtools)
library(ggalluvial)
# 准备细胞比例输入数据
prop_df <- seurat_object@meta.data %>%
  dplyr::select(Sample = orig.ident, Celltype = celltype) %>%
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
  geom_flow(width = 0.6, alpha = 0.3, knot.pos = 0.35) +  
  geom_col(width = 0.6,) +  
  scale_y_continuous(expand = c(0, 0)) +   
  scale_fill_manual(values = c("#ae7eb5",
                               "#f98177",
                               "#aad393",
                               "#568389",
                               "#eeb066",
                               "#f082a5",
                               "#54c6a8",
                               "#5a6fb5")) +  
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
# 🎨 Fig. 2c
# =============================================

# 每个细胞亚群在两组间的比例点图
load("CD45去双细胞后.Rdata")
# 查看数据集中的原始标识（不同实验条件的分组信息）
table(seurat_object$orig.ident)
# 计算每个细胞群体在总细胞中的比例
prop.table(table(Idents(seurat_object)))
# 按细胞群体和原始标识查看交叉表
table(Idents(seurat_object), seurat_object$orig.ident)
# 计算每个细胞群体在健康和酮病奶牛中的比例
Cellratio <- prop.table(table(Idents(seurat_object), seurat_object$orig.ident), margin=2)
# 设置R输出科学计数法的格式为更易读的数字形式
options(scipen=200)
# 查看比例数据
Cellratio
# 将比例数据转化为数据框
Cellratio <- as.data.frame(Cellratio)
Cellratio
# 将结果写入文本文件
write.table(Cellratio, file="CD45_Cellratio.txt", sep="\t", row.names=FALSE, quote=FALSE)
# 读取文件中的数据
CellRa <- read.table("CD45_Cellratio.txt", sep = "\t", header = T, na.strings = "NA", stringsAsFactors = FALSE)
# 创建新的分组变量
CellRa$group <- CellRa$Var2
# 移除group变量中的最后1个字符（例如：可能是样本的标识符后缀）
CellRa$group <- gsub('.{1}$', '', CellRa$group)
# 设置分组变量的因子顺序
CellRa$group <- factor(CellRa$group, levels=c("H", "K"))

# T
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="T",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 80) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
 
shapiro.test(CellRa[CellRa$Var1=="T",]$Freq[CellRa[CellRa$Var1=="T",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="T",]$Freq[CellRa[CellRa$Var1=="T",]$group == "K"]) 
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="T",], center = "mean") 
t.test(Freq ~ group, data = CellRa[CellRa$Var1=="T",]) # 0.003794 独立两组-正态-方差非齐

# NK
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="NK",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 8.2) +
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),
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
shapiro.test(CellRa[CellRa$Var1=="NK",]$Freq[CellRa[CellRa$Var1=="NK",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="NK",]$Freq[CellRa[CellRa$Var1=="NK",]$group == "K"]) 
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="NK",], center = "mean") 
t.test(Freq ~ group, data = CellRa[CellRa$Var1=="NK",], var.equal = TRUE) # 0.0007442 独立两组-正态-方差齐性

# B
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="B",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 8) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
  
shapiro.test(CellRa[CellRa$Var1=="B",]$Freq[CellRa[CellRa$Var1=="B",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="B",]$Freq[CellRa[CellRa$Var1=="B",]$group == "K"]) 
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="B",], center = "mean") 
t.test(Freq ~ group, data = CellRa[CellRa$Var1=="B",]) # 0.2085 独立两组-正态-方差非齐

# PC
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="PC",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 6) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
shapiro.test(CellRa[CellRa$Var1=="PC",]$Freq[CellRa[CellRa$Var1=="PC",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="PC",]$Freq[CellRa[CellRa$Var1=="PC",]$group == "K"]) 
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="PC",], center = "mean") 
t.test(Freq ~ group, data = CellRa[CellRa$Var1=="PC",], var.equal = TRUE) # 0.1864 独立两组-正态-方差齐性

# Mono
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="Mono",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 10.8) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
shapiro.test(CellRa[CellRa$Var1=="Mono",]$Freq[CellRa[CellRa$Var1=="Mono",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="Mono",]$Freq[CellRa[CellRa$Var1=="Mono",]$group == "K"]) 
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="Mono",], center = "mean") 
t.test(Freq ~ group, data = CellRa[CellRa$Var1=="Mono",], var.equal = TRUE) # 0.8908 独立两组-正态-方差齐性

# Mac
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="Mac",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 73) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),  
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
shapiro.test(CellRa[CellRa$Var1=="Mac",]$Freq[CellRa[CellRa$Var1=="Mac",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="Mac",]$Freq[CellRa[CellRa$Var1=="Mac",]$group == "K"]) 
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="Mac",], center = "mean") 
t.test(Freq ~ group, data = CellRa[CellRa$Var1=="Mac",], var.equal = TRUE) # 0.0008379 独立两组-正态-方差齐性

# DC
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="DC",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 7.4) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
shapiro.test(CellRa[CellRa$Var1=="DC",]$Freq[CellRa[CellRa$Var1=="DC",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="DC",]$Freq[CellRa[CellRa$Var1=="DC",]$group == "K"]) 
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="DC",], center = "mean") 
wilcox.test(Freq ~ group, data = CellRa[CellRa$Var1=="DC",]) # 0.3929 独立两组-非正态

# Neu
set.seed(123)
ggplot(data=CellRa[CellRa$Var1=="Neu",], aes(group, Freq * 100)) +
  geom_jitter(aes(fill = group), position = position_jitter(0), shape = 21, size = 3.8) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 5.7) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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

shapiro.test(CellRa[CellRa$Var1=="Neu",]$Freq[CellRa[CellRa$Var1=="Neu",]$group == "H"])
shapiro.test(CellRa[CellRa$Var1=="Neu",]$Freq[CellRa[CellRa$Var1=="Neu",]$group == "K"])
leveneTest(Freq ~ group, data = CellRa[CellRa$Var1=="Neu",], center = "mean") 
t.test(Freq ~ group, data = CellRa[CellRa$Var1=="Neu",], var.equal = TRUE) # 0.1638 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 2a
# =============================================

# 分组密度图
load("CD45去双细胞后.Rdata")

# 提取UMAP坐标
umap_df <- Embeddings(seurat_object, reduction = "umap") %>%
  as.data.frame()

# 确保列名统一为UMAP_1和UMAP_2
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# 添加meta信息 比如分组和细胞类型
umap_df$group <- seurat_object$group
umap_df$cell_class <- seurat_object$celltype

#加载需要的R包
library(ggplot2)
library(ggrepel)
library(ggforce)
library(viridis)
library(ggpubr)

# 计算每组每类细胞的中心位置
label_centroids <- umap_df %>%
  group_by(group, cell_class) %>%
  summarize(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = 'drop'
    )

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE, n = 200) +
  geom_point(color = "#FFFFFFAA", size = 0.001, alpha = 0.2) +
  facet_wrap(~group, ncol = 2) +
  scale_fill_viridis(option = "magma", direction = 1) +
  coord_cartesian(expand = FALSE) +
  theme_void() +
  theme(strip.text = element_text(face = "bold", size = 14),
        legend.position = "none") 


# =============================================
# 🎨 Fig. 2d
# =============================================

# 流式数据绘制箱线图
library(tidyverse)
flow_data <- data.frame(
  SampleID = c("H1", "H2", "H3", "H4", "K1", "K2", "K3", "K4"),
  Group = c("Health", "Health", "Health", "Health",
            "Ketosis", "Ketosis", "Ketosis", "Ketosis"),
  T = c(31.42, 42.17, 29.89, 36.78, 21.84, 17.63, 13.92, 20.47), 
  NK = c(10.36, 8.41, 11.72, 6.64, 6.02, 5.38, 3.75, 2.29),      
  B = c(8.73, 6.92, 5.58, 8.84, 8.27, 9.63, 13.52, 11.18), 
  PC = c(1.18, 1.66, 2.67, 0.99, 3.36, 2.63, 3.11, 3.57),        
  Mono = c(8.44, 5.38, 7.02, 6.21, 7.66, 10.14, 10.88, 8.03),   
  Mac = c(27.83, 24.54, 22.67, 26.34, 39.53, 33.87, 44.65, 49.71),
  DC = c(2.86, 2.02, 1.37, 2.91, 2.72, 3.11, 3.58, 4.67),        
  Neu = c(1.58, 1.50, 3.13, 2.28, 4.05, 2.98, 3.39, 4.03)        
  )

# 差异分析
t.test(T ~ Group, data = flow_data, var.equal = TRUE) # 0.00236
t.test(NK ~ Group, data = flow_data, var.equal = TRUE) # 0.01233
t.test(B ~ Group, data = flow_data, var.equal = TRUE) # 0.06237
t.test(PC ~ Group, data = flow_data, var.equal = TRUE) # 0.01117
t.test(Mono ~ Group, data = flow_data, var.equal = TRUE) # 0.05611
t.test(Mac ~ Group, data = flow_data, var.equal = TRUE) # 0.003547
t.test(DC ~ Group, data = flow_data, var.equal = TRUE) # 0.05611
t.test(Neu ~ Group, data = flow_data, var.equal = TRUE) # 0.01767

# 数据转换
flow_long <- flow_data %>%
  pivot_longer(
    cols = c(T, NK, B, PC, Mono, Mac, DC, Neu),
    names_to = "CellType",
    values_to = "Value"
  ) %>%
  mutate(CellType = factor(CellType, levels = c("T", "NK", "B", "PC", "Mono", "Mac", "DC", "Neu")))

# 保存
write.csv(flow_long, "flow_data_long.csv", row.names = FALSE)

# 作图
set.seed(127)
ggplot(flow_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.3) +
  geom_jitter(
    width = 0.25,
    size = 2.2,
    color = "black",
    alpha = 0.6) +
  facet_wrap(~ CellType, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Health" = "#4cdafe", "Ketosis" = "#ff6362")) +
  theme_bw() +
  theme(legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12,color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
    strip.text = element_text(size = 12),            
    panel.grid = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.38))) 


# =============================================
# 🎨 Fig. 2e
# =============================================

# F2的CD3E的荧光定量
health <- c(0.84096024, 
            0.936234059, 
            1.092273068, 
            1.130532633)
ketosis <- c(0.600900225, 
             0.455363841, 
             0.438109527, 
             0.376594149)
# 转换为长格式
flow_long <- data.frame(
  Group = rep(c("Health", "Ketosis"), each = 4),
  Value = c(health, ketosis)
  )
flow_long
flow_long$CellType <- "CD3E"
# 开始画图
set.seed(123)
ggplot(flow_long, aes(x = Group, y = Value, fill = Group)) +
  geom_jitter(aes(fill = Group), position = position_jitter(0), shape = 21, size = 3.5) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0, 1.4) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),
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
t.test(Value ~ Group, data = flow_long, var.equal = TRUE) # 0.0006636 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 2g
# =============================================

# CD68的荧光定量
health <- c(0.738461538,
            0.912820513,
            1.107692308,
            1.241025641)
ketosis <- c(1.764102564,
             2.082051282,
             2.307692308,
             2.471794872)
# 转换为长格式
flow_long <- data.frame(
  Group = rep(c("Health", "Ketosis"), each = 4),
  Value = c(health, ketosis)
  )
flow_long
flow_long$CellType <- "CD68"
# 开始画图
set.seed(123)
ggplot(flow_long, aes(x = Group, y = Value, fill = Group)) +
  geom_jitter(aes(fill = Group), position = position_jitter(0), shape = 21, size = 3.5) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  scale_x_discrete(labels = c("H" = "Health", "K" = "Ketosis")) +
  ylim(0.5, 3) + 
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5), 
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
t.test(Value ~ Group, data = flow_long, var.equal = TRUE) # 0.0008639 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 2f
# =============================================

# NCR1的荧光定量
health <- c(1.294429708,
            1.145888594,
            0.806366048,
            0.75331565)
ketosis <- c(0.583554377,
             0.339522546,
             0.297082228,
             0.201591512)
# 转换为长格式
flow_long <- data.frame(
  Group = rep(c("Health", "Ketosis"), each = 4),
  Value = c(health, ketosis)
  )
flow_long
flow_long$CellType <- "NCR1"
# 开始画图
set.seed(123)
ggplot(flow_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(aes(fill = Group), width = 0.5, outlier.shape = NA, color = "black", linewidth = 0.3) +
  geom_jitter(width = 0.25, size = 2.2, color = "black", fill = "black", shape = 21, alpha = 0.6) +
  scale_fill_manual(values=c("#4cdafe", "#ff6362")) + 
  ylim(0, 1.6) +  
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  
        axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust = 0.5, vjust = 0.5),  
        legend.position = "none")
dev.off()  
t.test(Value ~ Group, data = flow_long, var.equal = TRUE) # 0.00583 独立两组-正态-方差齐性


# =============================================
# 🎨 Fig. 2h
# =============================================

# 公共主数据库33189277 用Marker基因表达来佐证细胞数量的变化情况
rm(list = ls())
load("标准化后的表达矩阵.Rdata") 
load("Deseq2差异分析结果.Rdata") 

# 这是产前产后配对数据 产前19头 产后死亡两头 剩余17头
# 首先我把死亡两头对应的产前数据删掉 列名为T1_5B和T4_5B
colnames(标准化后的表达矩阵)
mat = 标准化后的表达矩阵[,-c(22,36)]
colnames(mat)

# 重新设置分组
Group <- c(rep("After", 17), rep("Before", 17))
Group <- factor(Group, levels = c("Before", "After"))
Group

# 存储整理好的数据
save(mat, Group, file = "产前产后配对标准化表达矩阵.Rdata") 
load("产前产后配对标准化表达矩阵.Rdata")

boxplot(mat, las = 2)  
dat = log2(mat + 1) 
boxplot(dat, las = 2)  
dat[1:4, 1:4]

dat = t(dat) %>%              
  as.data.frame(.) %>%        
  rownames_to_column(.) %>%   
  mutate(.,Group = Group) 

# CD3E
cd3e_df <- dat %>%
  dplyr::select(Sample = rowname, CD3E) %>%
  mutate(Group = Group)  
cd3e_df
cd3e_df <- cd3e_df %>%
  mutate(ID = rep(1:17, 2))  

set.seed(123)
ggplot(cd3e_df, aes(x = Group, y = CD3E, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
    )

shapiro.test(cd3e_df$CD3E[cd3e_df$Group == "Before"])
shapiro.test(cd3e_df$CD3E[cd3e_df$Group == "After"])
wilcox.test(cd3e_df$CD3E[cd3e_df$Group == "Before"], cd3e_df$CD3E[cd3e_df$Group == "After"], paired = TRUE) # 0.009338

# BCL11B
df <- dat %>%
  dplyr::select(Sample = rowname, BCL11B) %>%
  mutate(Group = Group)  
df <- df %>%
  mutate(ID = rep(1:17, 2))  

set.seed(123)
ggplot(df, aes(x = Group, y = BCL11B, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$BCL11B[df$Group == "Before"])
shapiro.test(df$BCL11B[df$Group == "After"])
t.test(df$BCL11B[df$Group == "Before"], df$BCL11B[df$Group == "After"], paired = TRUE) # 0.003524

# NCR1
df <- dat %>%
  dplyr::select(Sample = rowname, NCR1) %>%
  mutate(Group = Group) 
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = NCR1, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
    )

shapiro.test(df$NCR1[df$Group == "Before"])
shapiro.test(df$NCR1[df$Group == "After"])
t.test(df$NCR1[df$Group == "Before"], df$NCR1[df$Group == "After"], paired = TRUE) # 0.03122

# KLRD1
df <- dat %>%
  dplyr::select(Sample = rowname, KLRD1) %>%
  mutate(Group = Group) 
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = KLRD1, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$KLRD1[df$Group == "Before"])
shapiro.test(df$KLRD1[df$Group == "After"])
t.test(df$KLRD1[df$Group == "Before"], df$KLRD1[df$Group == "After"], paired = TRUE) # 8.552e-05


# KLRF1
df <- dat %>%
  dplyr::select(Sample = rowname, KLRF1) %>%
  mutate(Group = Group) 
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = KLRF1, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$KLRF1[df$Group == "Before"])
shapiro.test(df$KLRF1[df$Group == "After"])
wilcox.test(df$KLRF1[df$Group == "Before"], df$KLRF1[df$Group == "After"], paired = TRUE) # 0.003845

# CD40
df <- dat %>%
  dplyr::select(Sample = rowname, CD40) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))  

set.seed(123)
ggplot(df, aes(x = Group, y = CD40, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$CD40[df$Group == "Before"])
shapiro.test(df$CD40[df$Group == "After"])
t.test(df$CD40[df$Group == "Before"], df$CD40[df$Group == "After"], paired = TRUE) # 0.03185

# MZB1
df <- dat %>%
  dplyr::select(Sample = rowname, MZB1) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))

set.seed(123)
ggplot(df, aes(x = Group, y = MZB1, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$MZB1[df$Group == "Before"])
shapiro.test(df$MZB1[df$Group == "After"])
t.test(df$MZB1[df$Group == "Before"], df$MZB1[df$Group == "After"], paired = TRUE) # 0.1171

# DERL3
df <- dat %>%
  dplyr::select(Sample = rowname, DERL3) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = DERL3, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )
dev.off()

shapiro.test(df$DERL3[df$Group == "Before"])
shapiro.test(df$DERL3[df$Group == "After"])
wilcox.test(df$DERL3[df$Group == "Before"], df$DERL3[df$Group == "After"], paired = TRUE) # 0.0008392

# SIRPA
df <- dat %>%
  dplyr::select(Sample = rowname, SIRPA) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))

set.seed(123)
ggplot(df, aes(x = Group, y = SIRPA, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$SIRPA[df$Group == "Before"])
shapiro.test(df$SIRPA[df$Group == "After"])
t.test(df$SIRPA[df$Group == "Before"], df$SIRPA[df$Group == "After"], paired = TRUE) # 0.0003707

# FCN1
df <- dat %>%
  dplyr::select(Sample = rowname, FCN1) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = FCN1, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$FCN1[df$Group == "Before"])
shapiro.test(df$FCN1[df$Group == "After"])
t.test(df$FCN1[df$Group == "Before"], df$FCN1[df$Group == "After"], paired = TRUE) # 0.0002916

# VCAN
df <- dat %>%
  dplyr::select(Sample = rowname, VCAN) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = VCAN, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$VCAN[df$Group == "Before"])
shapiro.test(df$VCAN[df$Group == "After"])
t.test(df$VCAN[df$Group == "Before"], df$VCAN[df$Group == "After"], paired = TRUE) # 0.002992

# CD80
df <- dat %>%
  dplyr::select(Sample = rowname, CD80) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))  

set.seed(123)
ggplot(df, aes(x = Group, y = CD80, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$CD80[df$Group == "Before"])
shapiro.test(df$CD80[df$Group == "After"])
t.test(df$CD80[df$Group == "Before"], df$CD80[df$Group == "After"], paired = TRUE) # 6.445e-07


# CD86
df <- dat %>%
  dplyr::select(Sample = rowname, CD86) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))

set.seed(123)
ggplot(df, aes(x = Group, y = CD86, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$CD86[df$Group == "Before"])
shapiro.test(df$CD86[df$Group == "After"])
t.test(df$CD86[df$Group == "Before"], df$CD86[df$Group == "After"], paired = TRUE) # 0.02936

# CD68
df <- dat %>%
  dplyr::select(Sample = rowname, CD68) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = CD68, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )
dev.off()

shapiro.test(df$CD68[df$Group == "Before"])
shapiro.test(df$CD68[df$Group == "After"])
t.test(df$CD68[df$Group == "Before"], df$CD68[df$Group == "After"], paired = TRUE) # 0.0002783

# FLT3
df <- dat %>%
  dplyr::select(Sample = rowname, FLT3) %>%
  mutate(Group = Group) 
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = FLT3, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )
dev.off()

shapiro.test(df$FLT3[df$Group == "Before"])
shapiro.test(df$FLT3[df$Group == "After"])
t.test(df$FLT3[df$Group == "Before"], df$FLT3[df$Group == "After"], paired = TRUE) # 0.05603

# PLD4
df <- dat %>%
  dplyr::select(Sample = rowname, PLD4) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2))  

set.seed(123)
ggplot(df, aes(x = Group, y = PLD4, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$PLD4[df$Group == "Before"])
shapiro.test(df$PLD4[df$Group == "After"])
t.test(df$PLD4[df$Group == "Before"], df$PLD4[df$Group == "After"], paired = TRUE) # 0.02476

# S100A12
df <- dat %>%
  dplyr::select(Sample = rowname, S100A12) %>%
  mutate(Group = Group)  
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = S100A12, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$S100A12[df$Group == "Before"])
shapiro.test(df$S100A12[df$Group == "After"])
wilcox.test(df$S100A12[df$Group == "Before"], df$S100A12[df$Group == "After"], paired = TRUE) # 0.007904

# CXCL8
df <- dat %>%
  dplyr::select(Sample = rowname, CXCL8) %>%
  mutate(Group = Group) 
df
df <- df %>%
  mutate(ID = rep(1:17, 2)) 

set.seed(123)
ggplot(df, aes(x = Group, y = CXCL8, fill = Group)) +
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
    axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust = 0.5, vjust = 0.5)
  )

shapiro.test(df$CXCL8[df$Group == "Before"])
shapiro.test(df$CXCL8[df$Group == "After"])
t.test(df$CXCL8[df$Group == "Before"], df$CXCL8[df$Group == "After"], paired = TRUE) # 0.003843


# =============================================
# 🎨 Supplementary Figure 4
# =============================================

# 加载数据
load("CD45去双细胞后.Rdata") 

# 这里我们查找免疫相关的GO基因集 加载函数 
source("自定义函数之GO获取.R")  
  
# 获取牛所有的GO基因集信息
Bovine_GO_DATA <- get_GO_data("org.Bt.eg.db",
                              "ALL",          
                              "SYMBOL"
                              )
# 保存牛的所有GO基因集信息
save(Bovine_GO_DATA,
     file = "Bovine_GO_DATA.Rdata"
     )
  
# 根据关键词查找相关通路的ID
load("Bovine_GO_DATA.Rdata")  
GO_DATA = Bovine_GO_DATA
findGO("inflammatory response", method = "key") 
findGO("immune response", method = "key")
findGO("hepatic immune response", method = "key")
  
# 根据GO通路编号 提取参与该GO通路的所有基因
inflammatory_response <- getGO("GO:0006954")  
immune_response <- getGO("GO:0006955")      
regulation_of_cytokine_production <- getGO("GO:0001817")
response_to_cytokine <- getGO("GO:0034097")
activation_of_immune_response <- getGO("GO:0002250")
cell_activation_involved_in_immune_response <- getGO("GO:0002263")
  
# 进行基因集打分
seurat_object <- AddModuleScore(seurat_object, features = inflammatory_response, name = "inflammatory_response")
seurat_object <- AddModuleScore(seurat_object, features = immune_response, name = "immune_response")
seurat_object <- AddModuleScore(seurat_object, features = regulation_of_cytokine_production, name = "regulation_of_cytokine_production")
seurat_object <- AddModuleScore(seurat_object, features = response_to_cytokine, name = "response_to_cytokine")
seurat_object <- AddModuleScore(seurat_object, features = activation_of_immune_response, name = "activation_of_immune_response")
seurat_object <- AddModuleScore(seurat_object, features = cell_activation_involved_in_immune_response, name = "cell_activation_involved_in_immune_response")

# 基因集打分存储在metadata信息中
colnames(seurat_object@meta.data) 

# 提取需要修改的列名
col_names <- colnames(seurat_object@meta.data)[-6:-1]
for (name in col_names) {
  new_name <- gsub("1", "", name) 
  colnames(seurat_object@meta.data) <- gsub(name, new_name, colnames(seurat_object@meta.data))
  }

colnames(seurat_object@meta.data)

# 可视化
mydata1 <- FetchData(seurat_object, vars = c("celltype","inflammatory_response")) # 可以
mydata2 <- FetchData(seurat_object, vars = c("celltype","immune_response")) # 可以
mydata4 <- FetchData(seurat_object, vars = c("celltype","regulation_of_cytokine_production")) # 可以
mydata5 <- FetchData(seurat_object, vars = c("celltype","response_to_cytokine")) # 可以
mydata10 <- FetchData(seurat_object, vars = c("celltype","activation_of_immune_response")) # 可以
mydata11 <- FetchData(seurat_object, vars = c("celltype","cell_activation_involved_in_immune_response")) # 可以

ggplot(mydata1, aes(x = celltype, y = inflammatory_response)) + 
  geom_boxplot(position = position_dodge(0),
               aes(color = factor(celltype), fill = factor(celltype)),
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = c("#ae7eb5", "#f98177", "#aad393", "#568389",
                               "#eeb066", "#f082a5", "#54c6a8", "#5a6fb5"),
                    guide = "none") + 
  scale_color_manual(values = rep("gray20", 8),
                     guide = "none") +
  labs(x = NULL,  
       y = NULL) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"),
        legend.position = "none") 

ggplot(mydata2, aes(x = celltype, y = immune_response)) + 
  geom_boxplot(position = position_dodge(0),
               aes(color = factor(celltype), fill = factor(celltype)),
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = c("#ae7eb5", "#f98177", "#aad393", "#568389",
                               "#eeb066", "#f082a5", "#54c6a8", "#5a6fb5"),
                    guide = "none") + 
  scale_color_manual(values = rep("gray20", 8),
                     guide = "none") +
  labs(x = NULL,  
       y = NULL) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
        legend.position = "none") 

ggplot(mydata4, aes(x = celltype, y = regulation_of_cytokine_production)) + 
  geom_boxplot(position = position_dodge(0),
               aes(color = factor(celltype), fill = factor(celltype)),
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = c("#ae7eb5", "#f98177", "#aad393", "#568389",
                               "#eeb066", "#f082a5", "#54c6a8", "#5a6fb5"),
                    guide = "none") + 
  scale_color_manual(values = rep("gray20", 8),
                     guide = "none") + 
  labs(x = NULL,  
       y = NULL) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
        legend.position = "none") 

ggplot(mydata5, aes(x = celltype, y = response_to_cytokine)) + 
  geom_boxplot(position = position_dodge(0),
               aes(color = factor(celltype), fill = factor(celltype)),
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = c("#ae7eb5", "#f98177", "#aad393", "#568389",
                               "#eeb066", "#f082a5", "#54c6a8", "#5a6fb5"),
                    guide = "none") +  
  scale_color_manual(values = rep("gray20", 8),
                     guide = "none") +  
  labs(x = NULL,  
       y = NULL) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
        legend.position = "none") 
dev.off()  

ggplot(mydata10, aes(x = celltype, y = activation_of_immune_response)) + 
  geom_boxplot(position = position_dodge(0),
               aes(color = factor(celltype), fill = factor(celltype)),
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = c("#ae7eb5", "#f98177", "#aad393", "#568389",
                               "#eeb066", "#f082a5", "#54c6a8", "#5a6fb5"),
                    guide = "none") +  
  scale_color_manual(values = rep("gray20", 8),
                     guide = "none") +  
  labs(x = NULL,  
       y = NULL) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
        legend.position = "none")  

ggplot(mydata11, aes(x = celltype, y = cell_activation_involved_in_immune_response)) + 
  geom_boxplot(position = position_dodge(0),
               aes(color = factor(celltype), fill = factor(celltype)),
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = c("#ae7eb5", "#f98177", "#aad393", "#568389",
                               "#eeb066", "#f082a5", "#54c6a8", "#5a6fb5"),
                    guide = "none") +  
  scale_color_manual(values = rep("gray20", 8),
                     guide = "none") +  
  labs(x = NULL,  
       y = NULL) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"),  
        legend.position = "none")  

# 接下来画打分的UMAP图
colnames(seurat_object@meta.data) 
levels(seurat_object@meta.data[["celltype"]])

scplotter::FeatureStatPlot(seurat_object, 
                           plot_type = "dim", 
                           features = "inflammatory_response", 
                           reduction = "umap",
                           bg_cutoff = -Inf,
                           theme = "theme_void",
                           palette = "RdYlBu",
                           pt_size = 0.2,
                           raster = T
                           )

scplotter::FeatureStatPlot(seurat_object, 
                           plot_type = "dim", 
                           features = "immune_response", 
                           reduction = "umap",
                           bg_cutoff = -Inf,
                           theme = "theme_void",
                           palette = "RdYlBu",
                           pt_size = 0.2,
                           raster = T
                           )

scplotter::FeatureStatPlot(seurat_object, 
                           plot_type = "dim", 
                           features = "regulation_of_cytokine_production", 
                           reduction = "umap",
                           bg_cutoff = -Inf,
                           theme = "theme_void",
                           palette = "RdYlBu",
                           raster = T,
                           pt_size = 0.2
                           )

scplotter::FeatureStatPlot(seurat_object, 
                           plot_type = "dim", 
                           features = "response_to_cytokine", 
                           reduction = "umap",
                           bg_cutoff = -Inf,
                           theme = "theme_void",
                           palette = "RdYlBu",
                           raster = T,
                           pt_size = 0.2
                           )

scplotter::FeatureStatPlot(seurat_object, 
                           plot_type = "dim", 
                           features = "activation_of_immune_response", 
                           reduction = "umap",
                           bg_cutoff = -Inf,
                           theme = "theme_void",
                           palette = "RdYlBu",
                           raster = T,
                           pt_size = 0.2
                           )

scplotter::FeatureStatPlot(seurat_object, 
                           plot_type = "dim", 
                           features = "cell_activation_involved_in_immune_response", 
                           reduction = "umap",
                           bg_cutoff = -Inf,
                           theme = "theme_void",
                           palette = "RdYlBu",
                           raster = T,
                           pt_size = 0.2
                           )

# 保存一下打分后的结果
save(seurat_object,file = "CD45去双细胞后打分版") 



