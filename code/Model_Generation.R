rm(list = ls())

setwd("~/BioR/proj02/")
dir.create("02")
setwd(paste(getwd(),"02",sep = "/"))
dir.names <- getwd()
dir.create("tmp")

library(tidyverse) ## 加载R包
library(rms)
library(Hmisc)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(foreign)
library(survMisc)
library(plyr)
library(DT)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(timeROC)
library(survminer)
library(ResourceSelection)
library(survivalsvm) 
library(randomForestSRC)
library(randomForest)
library(varSelRF)
library(pROC)

load("~/BioR/proj02/01/tidyverse.Rdata")
load("~/BioR/proj02/00/FAP_expression.Rdata")
load("~/BioR/proj02/00/FAP_metadata.Rdata")

df <- t(FAP.all) %>% as.data.frame()

df$Samples <- rownames(df)

FAP.meta$Samples <- rownames(FAP.meta)

df <- merge(df,FAP.meta,by="Samples")

df <- df %>% dplyr::select(Samples,type,source, everything()) ## 调整一下数据，对应后续的操作

df <- df[,-3]

df$type <- ifelse(df$type == "Tumor",1,0)

df <- df[,which(colnames(df) %in% c("Samples","type",DEPx$SYMBOL))]

dft <- df[-grep("GSM",df$Samples),]
dfv <- df[grep("GSM",df$Samples),]



##### 先用随机森林进行变量筛选
dft.rsf <- dft[,-1]
dfv.rsf <- dfv[,-1]
# dft.rsf[,1] <- as.factor(dft.rsf[,1])
# dfv.rsf[,1] <- as.factor(dfv.rsf[,1])

set.seed(1)

out.rsf <- randomForest(type ~ . ,
                        data = dft.rsf,
                        ntree = 800,
                        nsplit = 1,
                        important=TRUE,
                        proximity=TRUE)

plot(out.rsf)

# 交叉验证
##先清除NA的数据样本，交叉验证不允许有NA
# na.omit

# 交叉验证中添加随机数的训练集、分组、交叉验证的次数
result=rfcv(dft.rsf[,2:ncol(dft.rsf)],dft.rsf$type,cv.fold = 5)
# 绘制错误率曲线，观察错误率与使用Markers数量的变化
with(result,plot(n.var,error.cv,log="x",type="o",lwd=2))

##使用replicate多次交叉验证
result=replicate(5,rfcv(dft.rsf[,2:ncol(dft.rsf)],dft.rsf$type),simplify = FALSE)

error.cv=sapply(result,"[[","error.cv")

matplot(result[[1]]$n.var,cbind(rowMeans(error.cv),error.cv),type = "l",lwd = c(2,rep(1,ncol(error.cv))),col = 1,lty = 1,log="x",xlab = "Number of variables",ylab="CV Error")

varImpPlot(out.rsf, main = "variable importance")

(randomForestSRC_geneids <- colnames(dft.rsf[,2:ncol(dft.rsf)])[order(out.rsf[["importance"]],decreasing = T)][1:30]) 

pred <- predict(out.rsf,newdata = dfv.rsf)

outcome.rsf <- roc(response=dfv.rsf$type,predictor = as.numeric(as.character(pred))) 

plot.roc(outcome.rsf,col="red",print.auc = TRUE)

save(out.rsf, file="Randomforest_model.Rdata")


######LASSO回归分析
library("glmnet") ## 加载R包
library(tidyverse)


## 根据R包的要求，将数据需要筛选的部分提取转换为矩阵,lasso不需要应变量
dft.lasso <- dft[,-c(1:2)] %>% as.matrix()
dfv.lasso <- dfv[,-c(1:2)] %>% as.matrix()


#df.lasso <- df.lasso[,colnames(df.lasso) %in% randomForestSRC_geneids]
dft$type <-as.factor(dft$type)
dfv$type <-as.factor(dfv$type)

cvfit = cv.glmnet(dft.lasso,
                  dft$type, 
                  nfold=10,#10倍交叉验证，非必须限定条件，其他文献大多没提
                  family = "binomial"
) 

plot(cvfit) ## 画图

fit <- glmnet(dft.lasso, dft$type, 
              family = "binomial") 

plot(fit)

coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个

active.min = which(coef.min != 0 ) ## 找出那些回归系数没有被惩罚为0的

(lasso_geneids <- colnames(dft.lasso)[active.min]) ## 提取基因名称

save(cvfit,file = 'lasso_model.RData')

pred <- predict(cvfit,newx = dfv.lasso,type = 'response')

outcome.lasso <- roc(response=dfv$type,predictor = as.numeric(as.character(pred))) 

plot.roc(outcome.lasso,col="red",print.auc = TRUE)


# svm算法进行分类模型构建。####
library(e1071)
dft.svm <- dft[,-c(1:2)] %>% as.matrix()
dfv.svm <- dfv[,-c(1:2)] %>% as.matrix()

dft$type <-as.factor(dft$type)
dft$type <-as.numeric(as.character(dft$type))
dfv$type <-as.factor(dfv$type)
dfv$type <-as.numeric(as.character(dfv$type))

obj <- tune(svm, dft.svm, dft$type,
            ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
            tunecontrol = tune.control(sampling = "fix")
)

out.svm <-svm(dft.svm,dft$type, scale=FALSE, kernel='radial', gamma=obj$best.parameters$gamma,
                cost= obj$best.parameters$cost)

svm_geneids <- ogj$best.parameters##变量选择代码还需研究

pred <- predict(out.svm, dfv.svm)

outcome.svm <- roc(response=dfv$type,predictor = as.numeric(as.character(pred))) 

plot.roc(outcome.svm,col="red",print.auc = TRUE)

save(out.svm,file = "svm_model.Rdata")

target <- intersect(lasso_geneids,randomForestSRC_geneids, svm_geneids)

lasso_res <- data.frame(Genes = lasso_geneids,
                        coef =coef.min[active.min])

save(lasso_geneids,lasso_res,target,file = "Diagnosis_model.Rdata")

save(df,file = "df.Rdata")

# Model Evaluation模型评估####
library(corrplot)
library(RColorBrewer)

corr_matrix1 <- df[,lasso_geneids]

corr <- cor(corr_matrix1,method = "pearson")

pdf("corr-sig.pdf",width = 8, height = 6.4)

corrplot(corr, order="AOE", method = "circle", type = "full", diag = FALSE, col.lim = c(-1,1), cl.pos = "r",shade.col = "white", col=colorRampPalette(brewer.pal(10,"PiYG"))(40),tl.cex = 1, tl.col = "black", tl.pos = "lt", mar = c(0,0,3,8))

dev.off()


rownames(df) <- df$Samples
df <- df[,colnames(df) %in% lasso_geneids]

df$signature <- as.matrix(df) %*% lasso_res$coef

df$score <- ifelse(df$signature < median(df$signature),"High-Score","Low-Score")


load("~/BioR/proj02/01/tidyverse.Rdata")
load("~/BioR/proj02/FAP_expression.Rdata")


testMatrix <- as.matrix(FAP.all)

save(testMatrix,file = "testMatrix.Rdata")

library(parallel)
library(pRRophetic)
library(ggplot2)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(cgp2016ExprRma)
possibleDrugs2016 <- unique( drugData2016$Drug.name)
possibleDrugs2016 <- possibleDrugs2016[-c(25,41,48,52,54,58,89,107,180,183,197,199,213,217)]##不报错的药物
save(possibleDrugs2016,file = "Worked_possibleDrug2016.Rdata")

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

df <- cbind(df[,c("signature","score")],t(res.df))

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

dir='DataFiles/Training Data/'
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

annotation_col <- data.frame(row.names = rownames(df),
                             RiskScore = df$score)

annotation_col <- annotation_col[order(df$signature),]

annotation_colors <- list(RiskScore=c("High-Score"='#5494cc',"Low-Score"='#f9cc52'))

res.df <- t(scale(t(res.df))) 

res.df[res.df < -2] <- -2
res.df[res.df > 2] <- 2

res.df <- res.df[,order(df$signature)]

DrugPredictions <- t(scale(t(DrugPredictions)))

DrugPredictions[DrugPredictions < -2] <- -2
DrugPredictions[DrugPredictions > 2] <- 2

normalization <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}
res.df <-normalization(res.df)
DrugPredictions <- normalization(DrugPredictions)

DrugPredictions <- DrugPredictions[,order(df$signature)]

ht1 = Heatmap(res.df, name = "CGP2016", 
              top_annotation = HeatmapAnnotation(df = annotation_col,
                                                 col = annotation_colors,height = unit(3,"mm")),
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

