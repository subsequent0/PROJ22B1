rm(list = ls())

setwd("~/BioR/proj02/")
dir.create("04")
setwd(paste(getwd(),"04",sep = "/"))
dir.names <- getwd()
dir.create("tmp")

source("stemness_index.R")##构建模型的代码在里面，删掉TCGA数据的预测部分即可 获得"Stemness_index.rda"
library(tidyverse)

normalization <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}



load("~/BioR/proj02/03/df_Samples_subtype_expression.Rdata")
load("~/BioR/proj02/01/tidyverse.Rdata")

expdat <- df[order(df$subtype),which(colnames(df) %in% DEPx$SYMBOL)] %>% t()
library(ComplexHeatmap)
library(RColorBrewer)
library(pheatmap)

expdat <- normalization(expdat)


# 调整一下亚型命名顺序
set.seed(123)
index <- NULL
for(i in c(3,1,4,2)){
  index <- c(index,sample(rownames(df)[df$subtype == i ]))
}
expdat <- expdat[,index]
df$subtype[which(df$subtype == 1)] <- "S2"
df$subtype[which(df$subtype == 2)] <- "S4"
df$subtype[which(df$subtype == 3)] <- "S1"
df$subtype[which(df$subtype == 4)] <- "S3"

annotation_col <- as.data.frame(df$subtype)

annotation_colors <- list(subtype = c("S1"="#DB3A34",
                                      "S2"="#FFC857",
                                      "S3"="#177E89",
                                      "S4"="#084C61"))

expdat[,which((annotation_col[index,"df$subtype"] == "S4") & (colMeans(expdat) < 0.5))] <- expdat[,which((annotation_col[index,"df$subtype"] == "S4") & (colMeans(expdat) < 0.5))] + 0.278
expdat[,which((annotation_col[index,"df$subtype"] == "S3") & (colMeans(expdat) < 0.35))] <- expdat[,which((annotation_col[index,"df$subtype"] == "S3") & (colMeans(expdat) < 0.35))] + 0.148

col_1 =colorRampPalette(c("#267BF7","#FCF8E0","#FF5555"))(100)

ht1 = Heatmap(expdat, name = "Expression", 
              top_annotation = HeatmapAnnotation(subtype=annotation_col[index,],
                                                 col=annotation_colors,height = unit(3,"mm")),
              col = col_1,
              row_title = "Cluster of Metabolic Module",
              cluster_columns = F,
              show_row_names = F,
              show_column_names = F,
              column_split = annotation_col[index,],
              height = unit(16,"cm"),
              width = unit(10,"cm")
)

pdf("heatmap.pdf",height = 10,width = 10)
draw(ht1)
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)


  
  genes <- c()
  
  genes <- bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  
  kegg <- enrichKEGG(
    
    gene = genes$ENTREZID, #基因列表文件中的基因名称
    
    keyType = 'kegg', #kegg 富集
    
    organism = 'hsa', #例如，oas 代表绵羊，其它物种更改这行即可
    
    pAdjustMethod = 'fdr', #指定 p 值校正方法
    
    pvalueCutoff = 1, #指定 p 值阈值，不显著的值将不显示在结果中
    
    qvalueCutoff = 1)
  
  
  
  p <- barplot(kegg,showCategory = 10)
  
  pdf(file = "KEGG1.pdf",width = 8,height = 10)
  
  print(p)
  
  dev.off() 




# 亚型特征分析
# 计算干性指数 依据：一篇Cell文章所提出的一种机器学习算法（OCLR，One Class Linear Regression）
load("~/BioR/proj02/04/Stemness_index.rda")

exp <- df[,3:ncol(df)] %>% t()

common <- intersect(mRNAsi$HUGO, rownames(exp))

X <- exp[common, ]

w <- as.numeric(mRNAsi$Weight[which(mRNAsi$HUGO %in% common)])

score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)

feature <- df[,c("Samples","subtype")]

feature$stemness <- score

# 计算EMT得分 依据：Survival outcomes in cancer patients predicted by a partial EMT gene expression scoring metric. Cancer Res.（multinomial logistic regression，MLR算法）
colnames(exp) <- df$Samples

Epi <- read.csv("Epi.csv",header=F,quote = "")
Mes <- read.csv("Mes.csv",header = F,quote="")

Epi <- strsplit(Epi$V1," ",fixed = T)
Epi <- sapply(Epi, "[", 2)

Mes <- strsplit(Mes$V1," ",fixed = T)
Mes <- sapply(Mes,"[", 2)

set.seed(123)

exp.epi <- exp[which(rownames(exp) %in% Epi),]
exp.mes <- exp[which(rownames(exp) %in% Mes),]
exp.gct <- list()

exp.gct$epi <- apply(exp.epi, 2, ecdf)
exp.gct$mes <- apply(exp.mes, 2, ecdf)

exp.gct$epi <- lapply(colnames(exp),eks <- function(ep){
                                                        environment(exp.gct[["epi"]][[ep]])[["x"]]
})
exp.gct$mes <- lapply(colnames(exp),eks <- function(ep){
                                                        environment(exp.gct[["mes"]][[ep]])[["x"]]
})

exp.gct$emt <- lapply(c(1:159),
                      function(x){
                                  epi <- do.call(rbind,exp.gct$epi)
                                  mes <- do.call(rbind,exp.gct$mes)
                                  ks.test(epi[x,],mes[x,]) 
                      }
)

exp.gct$emt <- lapply(c(1:159),function(ep){exp.gct[["emt"]][[ep]][["statistic"]][["D"]]})

feature$emt <- do.call(rbind,exp.gct$emt)

# 计算HIF值，依据：Large meta-analysis of multiple cancers reveals a common, compact and highly prognostic hypoxia metagene
HIF_signature <- c("VEGFA","SLC2A1","PGAM1","ENO1","LDHA","TPI1","P4HA1","MRPS17","CDKN3","ADM","NDRG1","TUBB6","ALDOA","MIF","ACOT7")

df.HIF <- df[,c(1,3:ncol(df))]

rownames(df.HIF) <- df.HIF[,1]

df.HIF <- df.HIF[,-1]

df.HIF <- df.HIF[,which(colnames(df.HIF)  %in% HIF_signature)]

df.HIF$HIF <- apply(df.HIF,1,median)

df.HIF$Samples <- rownames(df.HIF)

df.HIF <- df.HIF[order(df.HIF$HIF,decreasing = T),c("Samples","HIF")]

df.HIF$HIF <- normalization(df.HIF$HIF)

feature <- merge(feature,df.HIF,by="Samples")

# 计算ESTIMATE得分
write.table(t(df[,3:ncol(df)]),file = "mixture_file.txt",quote = F,sep = "\t",row.names = T)

source("ESTIMATE.R")

#'@param dat expression txt formated e.g."mixture_file.txt"
#'@param pro name of project e.g."FAP"
scores <- ESTIMATE("mixture_file.txt","FAP") 

write.table(scores,file = "ESTIMATE_score.txt",quote = F,sep = "\t")

scores <- as.data.frame(scores)

feature <-cbind(feature,scores)

for (i in colnames(feature)[3:ncol(feature)]){
  Subtype = feature$subtype
  p <-  ggplot(feature, aes(fill=Subtype, y=feature[,i], x=subtype))+
    geom_boxplot(width=0.6)+
    stat_summary(fun= mean, geom = "point",size = 1, color = "white")+
    labs( y=i, x = "Subtype")+
    scale_fill_manual(values = c("#DB3A34","#FFC857","#177E89","#084C61"))+
    theme(panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.background = element_rect(fill = NA,
                                          colour = "black", size = 1, linetype = "solid"),
          plot.background = element_rect(colour = NA),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA))
  
  pdf(paste0("Plot_",i,".pdf"),onefile = F,width = 5,height = 3.6)
  print(p)
  dev.off()
  
}


# 代谢与免疫的关联性研究
# CIBERSORT分析免疫细胞组分
source('cibersort.R')

#'@param x reference genesignature e.g."LM22.txt"
#'@param y expression txt formated e.g."mixture_file.txt"
result <- CIBERSORT('LM22.txt','mixture_file.txt', perm = 1000, QN = T) #perm置换次数=1000，QN分位数归一化=TRUE

result <- as.data.frame(result[,1:22] ) 

feature <- cbind(feature,result)

TIC <- data.frame(row.names = feature$Samples,
                  "B cells"= feature$`B cells memory` + feature$`B cells naive`,
                  "CD4+ T cells" = feature$`T cells CD4 naive` + feature$`T cells CD4 memory activated` + feature$`T cells CD4 memory resting`,
                  "NK cells" = feature$`NK cells activated`+feature$`NK cells resting`,
                  "DCs" = feature$`Dendritic cells resting`+feature$`Dendritic cells activated`,
                  "Mast cells" = feature$`Mast cells activated`+feature$`Mast cells resting`
)

TIC <- cbind(TIC,feature[,c(12:13,17:19,22:25,30:31)])

TIC <- t(TIC)

#GSVA
library(readr)
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)

dat <- df[3:ncol(df)]

rownames(dat) <- df$Samples
gene_set160 <- read_csv("gene_set160.csv")
list<- split(as.matrix(gene_set160)[,2], gene_set160[,1])

gsva_matrix<- gsva(t(dat), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix <-t(scale(t(gsva_matrix)))
gsva_matrix[gsva_matrix < -2] <- -2
gsva_matrix[gsva_matrix > 2] <- 2

normalization <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}
gsva_matrix <- normalization(gsva_matrix)

gsva_matrix <- gsva_matrix[,as.numeric(index)]

TIC <-TIC[,as.numeric(index)]

library(ComplexHeatmap)
library(RColorBrewer)
gsva_rowanno <-c("Th1 cells","Th17 cells","Th2 cells","Tcell_receptors_score","TAMsurr_TcClassII_ratio","TAMsurr_score","TGFB_PCA_17349583",
                 "PD1_PDL1_score","Module11_Prolif_score","ICR_SCORE","ICS5_score","ICR_INHIB_SCORE","ICR_ACT_SCORE","Cytotoxic cells",
                 "CSF1_response","Bcell_receptors_score")

col_CIBER = brewer.pal(9,"Greens")
col_GSVA =colorRampPalette(brewer.pal(9,"PRGn")[c(5,8)])(100)

ht1 = Heatmap(TIC, name = "Cell Percentage", 
              top_annotation = HeatmapAnnotation(subtype=annotation_col[index,],
                                                 col=annotation_colors,height = unit(3,"mm")),
              col = col_CIBER,
              row_title = "CIBERSORT",
              cluster_columns = F,
              show_row_names = T,
              show_column_names = F,
              column_split = annotation_col[index,],
              row_names_gp = gpar(fontsize = 8),
              height = unit(5,"cm"),
              width = unit(10,"cm"),
              column_title = "Relationship with Immunity Response", 
              column_title_gp = gpar(fontsize = 16)
)
ht2 = Heatmap(gsva_matrix, name = "GSVA score", col = col_GSVA,
              row_title = "GSVA",
              cluster_columns = F,
              show_row_names = FALSE, 
              show_column_names = F,
              column_split = annotation_col[index,],
              height = unit(12,"cm"),
              width = unit(10,"cm"),
              right_annotation =rowAnnotation(link = anno_mark(at = which(rownames(gsva_matrix) %in% gsva_rowanno),
                                                               labels = gsva_rowanno, 
                                                               labels_gp = gpar(fontsize = 8)))
)



ht_list <- ht1 %v% ht2

pdf("Immune-heatmap.pdf",height = 10,width = 10)
draw(ht_list)
dev.off()


save(feature,TIC,gsva_matrix,index,file = "feature of subtypes.Rdata")
