rm(list = ls())


library(dplyr)


df <- df[,c(1,ncol(df))]
expdata <- merge(df,expdata,by="Samples")

pdf("RS-PDL1.pdf",width = 4.5,height = 4)
p <- ggplot(expdata,aes(x=signature, y=CD274 ))+
  geom_point(size=2 , alpha=0.5)+
  labs( y="PDL1(CD274) ", x = "RiskScore")+
  geom_vline(xintercept = median(expdata$signature),linetype=2,cex=1)+
  stat_smooth(method = "gam",color="#00497D")+ 
  coord_cartesian(ylim = c(0, 1200))+
  geom_text(aes(x = -2.5, y =1200, 
                label = paste0("P = ",as.character(
                round(cor.test(expdata$signature,expdata$CD274)[["p.value"]],3)
                  )),
                parse = TRUE))+ 
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
p
dev.off()

pdf("RS-CTLA4.pdf",width = 4.5,height = 4)
p <- ggplot(expdata,aes(x=signature, y=CTLA4 ))+
  geom_point(size=2 , alpha=0.5)+
  labs( y="CTLA4 ", x = "RiskScore")+
  geom_vline(xintercept = median(expdata$signature),linetype=2,cex=1)+
  stat_smooth(method = "gam",color="#00497D")+ 
  coord_cartesian(ylim = c(0, 750))+
  geom_text(aes(x = -2.5, y =750, 
                label = paste0("P = ",as.character(
                  round(cor.test(expdata$signature,expdata$CTLA4)[["p.value"]],3)
                )),
                parse = TRUE))+ 
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
p
dev.off()


testMatrix <- expdata[,4:ncol(expdata)]
rownames(testMatrix) <- expdata$Samples
testMatrix <- t(testMatrix)
save(testMatrix,file = "testMatrix.Rdata")

library(parallel)
library(pRRophetic)
library(ggplot2)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(cgp2016ExprRma)
# possibleDrugs2016 <- unique( drugData2016$Drug.name)
# possibleDrugs2016 <- possibleDrugs2016[-c(25,41,48,52,54,58,89,107,180,183,197,199,213,217)]##不报错的药物
# save(possibleDrugs2016,file = "Worked_possibleDrug2016.Rdata")

load("Worked_possibleDrug2016.Rdata")

# 用system.time来返回计算所需时间
head(possibleDrugs2016)
system.time({ 
  cl <- makeCluster(8)  
  results <- parLapply(cl,possibleDrugs2016,
                       function(x){
                         library(pRRophetic) 
                         load("testMatrix.Rdata")
                         predictedPtype = pRRopheticPredict(
                           testMatrix,
                           drug = x,
                           tissueType = "urogenital_system", 
                           batchCorrect = "eb",
                           selection = 1,
                           dataset = "cgp2016")
                         return(predictedPtype)
                       }) # lapply的并行版本
  res.df <- do.call('rbind',results) # 整合结果
  stopCluster(cl) # 关闭集群
})

rownames(res.df) <- possibleDrugs2016
save(res.df,file = "CGP2016_drug_predict.Rdata")

df <- cbind(df,t(res.df))

cort <- function(x){R = cor.test(df$signature,df[,x])[["estimate"]][["cor"]]
            P = cor.test(df$signature,df[,x])[["p.value"]]
            COR <-c(R,P)
            return(COR)}
result <- lapply(colnames(df)[3:ncol(df)], cort)
result.cor1 <- do.call(rbind,result)
colnames(result.cor1) <- c("R","P")
rownames(result.cor1) <- colnames(df)[3:ncol(df)]
result.cor1 <- as.data.frame(result.cor1)



##oncoPredict on GDSCv2 (R>4.1.2required)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

dir='/DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=paste0(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = paste0(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 
load("testMatrix.Rdata")


calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testMatrix,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

DrugPredictions <- read.csv("calcPhenotype_Output/DrugPredictions.csv",header = T,quote = "",row.names = 1)

rownames(DrugPredictions) <- substr(rownames(DrugPredictions),2,17)

colnames(DrugPredictions) <- substring(colnames(DrugPredictions),3)

colnames(DrugPredictions) <-  do.call(rbind,strsplit(colnames(DrugPredictions),"_",fixed = T))[,1]

DrugPredictions <-t(DrugPredictions)

df <- df[,1:2] 
df <- cbind(df,t(DrugPredictions))
result <- lapply(colnames(df)[3:ncol(df)], cort)
result.cor2 <- do.call(rbind,result)
colnames(result.cor2) <- c("R","P")
rownames(result.cor2) <- colnames(df)[3:ncol(df)]
result.cor2 <- as.data.frame(result.cor2)


library(ComplexHeatmap)
library(RColorBrewer)
CGP_rowanno <- rownames(result.cor1[order(abs(result.cor1$R)),])[which(result.cor1$P < 0.001)][1:25]
GDSC_rowanno <-rownames(result.cor2[order(abs(result.cor2$R)),])[which(result.cor2$P < 0.001)][1:25]

col_CGP = colorRampPalette(brewer.pal(9,"PRGn")[c(2,5,8)])(100)
col_GDSC =colorRampPalette(brewer.pal(9,"PRGn")[c(2,5,8)])(100)

load("~/BioR/proj01/06/RiskScore+ImmuneScore.Rdata")
annotation_col <- data.frame(row.names = ImmuInfo$Samples,
                             RiskScore = ImmuInfo$signature_by2)
annotation_colors <- list(RiskScore=c("High-Score"='#5494cc',"Low-Score"='#f9cc52'))

res.df <- t(scale(t(res.df))) 

res.df[res.df < -2] <- -2
res.df[res.df > 2] <- 2

DrugPredictions <- t(scale(t(DrugPredictions)))

DrugPredictions[DrugPredictions < -2] <- -2
DrugPredictions[DrugPredictions > 2] <- 2

normalization <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}
res.df <-normalization(res.df)
DrugPredictions <- normalization(DrugPredictions)


ht1 = Heatmap(res.df, name = "CGP2016", 
              top_annotation = HeatmapAnnotation(df=annotation_col[order(annotation_col$RiskScore,decreasing = T),],
                                                 col=annotation_colors,height = unit(3,"mm")),
              col = col_CGP,
              row_title = "CGP2016",
              cluster_columns = F,
              show_row_names = F,
              show_column_names = F,
              height = unit(8.16,"cm"),
              column_title = "Relationship with Drug Response", 
              column_title_gp = gpar(fontsize = 16),
              right_annotation =rowAnnotation(link = anno_mark(at = which(rownames(res.df) %in% CGP_rowanno),
                                                               labels = CGP_rowanno, 
                                                               labels_gp = gpar(fontsize = 8)))
)

atix <-which(rownames(DrugPredictions) %in% GDSC_rowanno)
atix <- atix[!atix %in% which(duplicated(rownames(DrugPredictions)))]
ht2 = Heatmap(DrugPredictions, name = "GDSCv2", col = col_GDSC,
              row_title = "GDSC v2",
              cluster_columns = F,
              show_row_names = FALSE, 
              show_column_names = F,
              height = unit(8,"cm"),
              right_annotation =rowAnnotation(link = anno_mark(at = atix,
                                                               labels = GDSC_rowanno, 
                                                               labels_gp = gpar(fontsize = 8)))
)



ht_list <- ht1 %v% ht2

pdf("heatmap.pdf",height = 12,width = 10)
draw(ht_list)
dev.off()


library(corrplot)
load("~/BioR/proj01/02/Coefficiency.RData")
load("~/BioR/proj01/02/riskscore.Rdata")
corr_matrix1 <- df[,Final_GetFactors]

corr <- cor(corr_matrix1,method = "pearson")

sig <- cor.mtest(corr_matrix1)

pdf("corr-sig.pdf",width = 8, height = 6.4)

corrplot(corr, order="AOE", method = "circle", type = "lower", diag = FALSE, cl.lim = c(-1,1), cl.pos = "r",shade.col = "white", col=colorRampPalette(brewer.pal(10,"RdYlBu"))(10),tl.cex = 1, tl.col = "black", tl.pos = "l", mar = c(0,0,3,8))

#corrplot(corr, order="AOE",add = TRUE, method = "color", p.mat = sig$p, number.cex = 0.8, type = "upper", diag = FALSE, cl.lim = c(-1,1), cl.pos = "r", shade.col = "white", insig = "p-value",sig.level = -1, col="black",tl.cex = 0.4, tl.col = "black", tl.pos = "n")

dev.off()







