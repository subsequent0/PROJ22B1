rm(list = ls())


df <- read.csv("../00/TCGA_rawdata/dataFilt.COAD.final.patient.csv",header = T,quote = "",check.names = F,row.names = 1)

df <- as.data.frame(t(df))


df <- df[,colnames(df) %in% lasso_geneids]

df$signature <- as.matrix(df) %*% lasso_res$coef

df$score <- ifelse(df$signature > median(df$signature),"High-Score","Low-Score")

download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.GDC_phenotype.tsv.gz","TCGA-COAD.phenotype.tsv.gz")

clinicdata <- read.table(gzfile('TCGA-COAD.phenotype.tsv.gz',open="rt"),header=T,sep = '\t',quote = '')

download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.survival.tsv","TCGA-COAD.survival.tsv")

surviv.info <- read.table("TCGA-COAD.survival.tsv",header = T,sep = '\t',quote = '')

#提取生存状况
surviv.info <- surviv.info[,colnames(surviv.info) %in% c("sample", 
                                                         "OS", 
                                                         "OS.time")]

clinicdata <- clinicdata[,colnames(clinicdata) %in% c("submitter_id.samples",
                                                      "gender.demographic",
                                                      "age_at_initial_pathologic_diagnosis",
                                                      "pathologic_M",
                                                      "pathologic_N",
                                                      "pathologic_T",
                                                      "tumor_stage.diagnoses")]

clinicdata <- clinicdata[substring(clinicdata$submitter_id.samples,14,14) == 0 ,]

surviv.info <- surviv.info[substring(surviv.info$sample,14,14) == 0 ,]

names(clinicdata)[names(clinicdata) == 'submitter_id.samples'] <- 'sample'
names(surviv.info)[names(surviv.info) == 'OS'] <- 'fustatus'
names(surviv.info)[names(surviv.info) == 'OS.time'] <- 'futime'

df_omit_clear <- merge(clinicdata, surviv.info, by = "sample")
if (length(which(df_omit_clear[nrow(df_omit_clear),] == '')) == 0 ){
  x <- unique(which(df_omit_clear[,2:ncol(df_omit_clear)] == '') %% nrow(df_omit_clear))
}else{
  x <-c(unique(which(df_omit_clear[,2:ncol(df_omit_clear)] == '') %% nrow(df_omit_clear)),nrow(df_omit_clear))
}
df_omit_clear <- df_omit_clear[!( as.numeric(rownames(df_omit_clear) %in% x )),]
rownames(df_omit_clear) <- df_omit_clear$sample %>% as.character()#行名
df_omit_clear <- df_omit_clear[,-1]
df_omit_clear$fustatus <- as.factor(df_omit_clear$fustatus)

df_omit_clear <- df_omit_clear[match(rownames(df) ,rownames(df_omit_clear)),]
df_omit_clear$signature <- df$signature
df_omit_clear$score <- df$score

df <- cbind(df_omit_clear,df[,1:15])

df <- df[!is.na(df$futime),]

df$fustatus <- as.numeric(as.character(df$fustatus))

df_omit_clear <- df

df_omit_clear$GPS <- ifelse(df$GPS > mean(df$GPS),"High","Low")

fit1 <- survfit(Surv(futime, fustatus) ~ GPS,data = df_omit_clear)

p1 <- ggsurvplot(fit1, conf.int=F, pval=T, risk.table=T, legend.labs = c('High','Low'), 
                legend.title='GPR182+ PSC score', palette = c("dodgerblue2","#FFCC22"),
                risk.table.height = 0.3)

pdf('Plot_KMplot.pdf',onefile = F)
print(p1)
dev.off()


