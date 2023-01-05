
########FAP RNA-seq Data process

file_names <- paste("FAP_rawdata/GSE88945_RAW", dir("FAP_rawdata/GSE88945_RAW", pattern = glob2rx("GSM*_H_G*.txt")),sep = "/" )
data_list <- lapply(file_names, read.table, header = TRUE, quote="", sep= "\t") 

for(i in 2:length(data_list)){
  data_list[[i]] <- data_list[[i]][match(data_list[[1]][,1],data_list[[i]][,1]),]
}

GSE88945 <- do.call(cbind,data_list)

GSE88945[1:5,]

GSE88945 <- GSE88945[,c(1,5,12,19)]

sample <- strsplit(dir("FAP_rawdata/GSE88945_RAW", pattern = glob2rx("GSM*_H_G*.txt")),"_",fixed = T)

for (i in c(1:length(sample))){colnames(GSE88945)[i+1] <- paste0("H_",sample[[i]][3])}

colnames(GSE88945)[1] <- "Gene"

file_names <- paste("FAP_rawdata/GSE106500_RAW", dir("FAP_rawdata/GSE106500_RAW", pattern = glob2rx("GSM*_G*.txt")),sep = "/" )
data_list <- lapply(file_names, read.table, header = TRUE, quote="", sep= "\t") 

for(i in 2:length(data_list)){
  data_list[[i]] <- data_list[[i]][match(data_list[[1]][,1],data_list[[i]][,1]),]
}


GSE106500 <- do.call(cbind,data_list)

GSE106500[1:5,]

GSE106500 <- GSE106500[,c(1,which(colnames(GSE106500) == "expected_count"))]

sample <- strsplit(dir("FAP_rawdata/GSE106500_RAW", pattern = glob2rx("GSM*_G*.txt")),"_",fixed = T)

for (i in c(1:length(sample))){colnames(GSE106500)[i+1] <- paste0("H_",sample[[i]][3])}

colnames(GSE106500)[1] <- "Gene"

FAP <- merge(GSE106500,GSE88945,by="Gene")

GSE153385 <- read.table("FAP_rawdata/GSE153385_LS.FAP.patient.raw.read.counts.txt",header = TRUE, quote="", sep= "\t")

GSE153385[1:5,]

GSE153385 <- GSE153385[,c(1,which(!colnames(GSE153385) %in% colnames(FAP)))]

FAP <- merge(FAP,GSE153385,by="Gene")

GSE156172 <- read.table("FAP_rawdata/GSE156172_FAP_SulBex.raw.read.counts.txt",header = TRUE, quote="", sep= "\t")

GSE156172[1:5,]

colnames(GSE156172)[1] <- "Gene"

FAP <- merge(FAP,GSE156172,by="Gene")

library(GenomicFeatures)#载入包
library(RMariaDB)
txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",#设定种族
                            release=101,#设定版本
                            server="ensembldb.ensembl.org")#审定服务器地址

saveDb(txdb, file='normalizeGeneCounts/txdbensemble101.sqlite')#保存生存的TXDB文件

TxDb <- loadDb(file='normalizeGeneCounts/txdbensemble101.sqlite')#载入参考数据集

source("normalizeGeneCounts/normalizeGeneCounts.R")#载入本地脚本

FAP_RPKM <-FAP

rownames(FAP_RPKM) <- FAP_RPKM$Gene

FAP_RPKM <- FAP_RPKM[,-1]

FAP_RPKM <- normalizeGeneCounts(FAP_RPKM, TxDb, method = "RPKM")

FAP[,2:ncol(FAP)] <- FAP_RPKM

GSE94919 <- read.table("FAP_rawdata/GSE94919_43FAP_RPKM.txt",header = TRUE, quote="", sep= "\t")

GSE94919[1:5,]

FAP <- merge(FAP,GSE94919,by="Gene")

rownames(FAP) <- FAP$Gene

FAP <- FAP[,-1]

library(sva)
library(bladderbatch)
library("FactoMineR")
library("factoextra")

data_list <- list("GSE106500"=colnames(GSE106500)[2:ncol(GSE106500)],
                  "GSE153385"=colnames(GSE153385)[2:ncol(GSE153385)],
                  "GSE156172"=colnames(GSE156172)[2:ncol(GSE156172)],
                  "GSE88945"=colnames(GSE88945)[2:ncol(GSE88945)],
                  "GSE94919"=colnames(GSE94919)[2:ncol(GSE94919)])
colmeta <- reshape2::melt(data_list,variable.name="dataset",value.name = c("samples"))[,2:1]

colmeta <- colmeta[match(colnames(FAP),colmeta$samples),]

pca.plot = function(dat,col){
  
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups"
  )
}


save(FAP,file = "tmp/FAP_T+N_expression_df.Rdata")
save(colmeta,file = "tmp/FAP_T+N_colnames_meta.Rdata")

###### Chip Data Process

if(!require("GEOquery")) BiocManager::install("GEOquery")
library(GEOquery)
library(dplyr)

gse <- getGEO("GSE9689", GSEMatrix = TRUE)

show(gse)

metdata <-pData(gse[[1]])

expdata <- exprs(gse[[1]])

expdata <- as.data.frame(expdata)

GPL="GPL3408"##下载平台注释

gpl <- getGEO(GPL,destdir = './tmp')  

gpl <- Table(gpl) 

head(gpl)

colnames(gpl)

expdata <- expdata[which(gpl$Gene_Symbol != ""),]

expdata$Gene  <- gpl$Gene_Symbol[which(gpl$Gene_Symbol != "")]

expdata$Gene[grep("***",expdata$Gene,fixed = T)] <- substring(expdata$Gene[grep("***",expdata$Gene,fixed = T)],4)

expdata <- aggregate(x = expdata,by=list(expdata$Gene),FUN=max)

expdata <- expdata[,-1]

expGSE9689 <- expdata

metaGSE9689 <- metdata[,36:37]

metaGSE9689 <- metaGSE9689[which(metaGSE9689$`MUTATION:ch1`!="MYH"),]

colnames(metaGSE9689) <- c("type","source")

metaGSE9689$source <- "GSE9689"

expGSE9689 <- expGSE9689[,c(rownames(metaGSE9689),"Gene")]

expGSE79460 <- read.table("GSE79460_expression.txt",header = T,quote = "",sep = "\t")

colnames(expGSE79460)[1] <- "Gene"

metaGSE79640 <- read.table("FAP_rawdata/GSE79460_clinical.txt",header = T,quote = "",sep = "\t")

metaGSE79640 <- metaGSE79640[8:16,]

rownames(metaGSE79640) <- metaGSE79640$accession

metaGSE79640 <- metaGSE79640[,c(2,4)]

colnames(metaGSE79640)[1] <- "type"

metaGSE79640$type <- "Tumor"

metaGSE79640$source <- "GSE79640"

expGSE109812 <- read.table("FAP_rawdata/GSE109812_expression.txt",header = T,quote = "",sep = "\t")

colnames(expGSE109812)[1] <- "Gene"

metaGSE109812 <- data.frame(type=rep("Tumor",ncol(expGSE109812)-1),
                            source=rep("GSE109812",ncol(expGSE109812)-1),
                            row.names = colnames(expGSE109812)[2:ncol(expGSE109812)])

FAP_chip <- merge(expGSE109812,expGSE79460,by="Gene")

FAP_chip <- merge(FAP_chip,expGSE9689,by="Gene")

metaFAP_chip <-rbind(metaGSE109812,metaGSE79640)

metaFAP_chip <-rbind(metaFAP_chip ,metaGSE9689)

save(FAP_chip,file = "tmp/FAP_T+N_chip.Rdata")

save(metaFAP_chip,file = "tmp/FAP_T+N_chip_colDATA.Rdata")

rownames(FAP_chip) <- FAP_chip$Gene

FAP_chip <-FAP_chip[,-1]


######Intergrate and remove batch effect

boxplot(FAP, las=2)
FAP <- log(FAP+1)

boxplot(FAP, las=2)
boxplot(FAP_chip)

FAP_chip[,1:40] <- log(FAP_chip[,1:40])

boxplot(FAP_chip)

FAP$Gene <- rownames(FAP)
FAP_chip$Gene <- rownames(FAP_chip)
FAP.all <- merge(FAP,FAP_chip,by="Gene")
rownames(FAP.all) <- FAP.all$Gene
FAP.all <- FAP.all[,-1]
boxplot(FAP.all)

batch <- c(colmeta$L1,metaFAP_chip$source)

pca.plot(FAP.all,factor(batch))

combat_FAP.all<- ComBat(dat = as.matrix(FAP.all), batch = batch,
                    par.prior = F)



combat_FAP.all[which(is.na(combat_FAP.all)) ] <- 8

boxplot(combat_FAP.all)

pca.plot(combat_FAP.all,factor(batch))

FAP.all <- combat_FAP.all

rownames(colmeta) <- colmeta$samples

colmeta$samples <-"Tumor"

Normal <- c("H_G114" ,"H_G96" ,"H_G110","H_G84","H_G86","H_G93","H_G106","H_G95","H_G121","H_G1","H_G54","H_G4","H_G16","H_G26","H_G31","H_G40","H_G49","H_G60","H_G62")

colmeta$samples[which(rownames(colmeta) %in% Normal)] <- "Normal"

colmeta <- colmeta[,c(2,1)]

colnames(colmeta) <- c("type","source")

FAP.meta <- rbind(colmeta,metaFAP_chip)

save(FAP.all,file = "FAP_expression.Rdata")

save(FAP.meta,file = "FAP_metadata.Rdata")
