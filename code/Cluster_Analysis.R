rm(list = ls())

setwd("~/BioR/proj02/")
dir.create("03")
setwd(paste(getwd(),"03",sep = "/"))
dir.names <- getwd()
dir.create("tmp")

library(factoextra)
library(NbClust)
library(ggplot2)
library(ggsci)
library(vegan)
library(fpc)
library(cluster)
library(ConsensusClusterPlus)
dir.names='ConsensusClusterPlus'


load("~/BioR/proj02/01/tidyverse.Rdata")
load("~/BioR/proj02/02/df.Rdata")

df.fap <- df[df$type ==1,which(colnames(df) %in% c("Samples",DEP))]

clustdat <- df.fap[,2:ncol(df.fap)]

## Calinsky criterion
calinski_clust <- cascadeKM(clustdat, 1, 10, iter = 1000)

calinski_clust$results

calinski.best <- as.numeric(which.max(calinski_clust$results[2,]))

calinski.best


## PAM Kmeans 
pamk.best <- pamk(clustdat)

pamk.best$nc

clusplot(pam(clustdat, pamk.best$nc),main = "Best Cluster",sub=NULL,col.clus="black")

## nbClust 
nc <- NbClust(clustdat, distance = "euclidean", min.nc = 2,max.nc = 15,method = "ward.D2")


## 综合上述结果确定最佳的聚类值，并修改下面的xintercept和centers
## WSS聚类分析 
fviz_nbclust(clustdat, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)

km.res <- kmeans(clustdat, centers =  4)

fviz_cluster(km.res, data = clustdat,geom = "point")+ 
  scale_color_npg()+
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



## ConsensusCluster 一致性聚类
clustdat <- t(clustdat)

results <- ConsensusClusterPlus(clustdat,
                                maxK=6,
                                reps=500,
                                pItem=0.8,
                                pFeature=1,
                                tmyPal = "#484452",
                                title=dir.names,
                                clusterAlg='pam',
                                distance='pearson',
                                seed=123,
                                innerLinkage = "complete",
                                plot='pdf')

icl = calcICL(results,
              title=dir.names,
              plot="pdf")

table(results[[4]]$consensusClass)

df.fap$subtype <- results[[4]]$consensusClass

df.fap <- subset(df.fap,select=c("Samples","subtype"))

load("~/BioR/proj02/FAP_expression.Rdata")

df <- t(FAP.all) %>% as.data.frame()

df$Samples <- rownames(df)

df <- df[which(df$Samples %in% df.fap$Samples),]

df <- merge(df.fap, df, by="Samples")

save(df,file="df_Samples_subtype_expression.Rdata")