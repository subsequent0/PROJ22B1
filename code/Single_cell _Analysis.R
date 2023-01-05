rm(list = ls())

library(Seurat)

raw.data <- read.table("GSE109308_Colon_FAP_TPM_new2.txt",header = T,sep = "\t",quote = "",row.names = 1)

sce.all <- CreateSeuratObject(counts = raw.data)

rm(raw.data)


sce.all=PercentageFeatureSet(sce.all, "^mt-", col.name = "percent_mito")

sce.all=PercentageFeatureSet(sce.all, "^Rp[sl]", col.name = "percent_ribo")

sce.all=PercentageFeatureSet(sce.all, "^Hb[^(p)]", col.name = "percent_hb")

selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 300)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]
sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)

selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 20)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo >= 0)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 0.1)
sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
save(sce.all.filt,file = "sce.all.filt.Rdata")

library(devtools)
library(clustree)
library(tidyverse)
library(gridExtra)
library(ggridges)
library(ggplot2)
library(ggExtra)
library(DoubletFinder)
library(cowplot)
library(dplyr)
library(data.table)
library(ggsci)
library(RColorBrewer)

sce.all=sce.all.filt
sce.all = NormalizeData(sce.all)
sce.all <- FindVariableFeatures(sce.all)

sce.all = ScaleData(sce.all 
                    #,vars.to.regress = c("nFeature_RNA", "percent_mito")
)  

sce.all <- RunPCA(sce.all, npcs = 10)
sce.all <- RunTSNE(sce.all, dims = 1:10)
sce.all <- RunUMAP(sce.all, dims = 1:10)


sce.all <- FindNeighbors(sce.all, dims = 1:10)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all <- FindClusters(sce.all, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
clustree(sce.all@meta.data, prefix = "RNA_snn_res.") 

plot_grid(ncol = 3,
          DimPlot(sce.all, reduction = "pca")+NoAxes()+ggtitle("PCA raw_data"),
          DimPlot(sce.all, reduction = "tsne")+NoAxes()+ggtitle("tSNE raw_data"),
          DimPlot(sce.all, reduction = "umap")+NoAxes()+ggtitle("UMAP raw_data")
          
)

DimPlot(sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.05", label = T,label.box = T)


def_color <- c("#D86D30","#F57A01","#25519B","#023459","#0A448D","#2D76A8","#007BBF","#075FF3","#9DC513")

p_sample <- DimPlot(sce.all, reduction = "tsne", group.by = "orig.ident")+scale_color_manual(values = def_color)

p_sample

library(SingleR)
library(celldex)

HGA <-  celldex::HumanPrimaryCellAtlasData()
#save(HGA,file = "HumanPrimaryCellAtlasData.Rdata")
#load("HumanPrimaryCellAtlasData.Rdata")

sce.singleR <- GetAssayData(sce.all,slot = "data")

clusters <- sce.all@meta.data$RNA_snn_res.0.1

pred.sce <- SingleR(test = sce.singleR, ref = HGA,
                    labels = HGA$label.main,
                    clusters = clusters)

celltype <- data.frame(clusterID = levels(sce.all@meta.data$RNA_snn_res.0.1),
                       celltype = pred.sce$labels)

sce.all@meta.data$singleR <-  celltype[match(clusters, celltype$clusterID),"celltype"]

p_cell <- DimPlot(sce.all, reduction = "tsne", group.by = "singleR")+scale_color_npg()

p_cell


cell.use <- row.names(sce.all@meta.data)[grep("FAP",sce.all@meta.data$orig.ident,fixed = T)]

sce_second <- subset(sce.all, cells=cell.use)

DefaultAssay(sce_second) <- "RNA"

sce_second <- NormalizeData(sce_second, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce_second <- FindVariableFeatures(sce_second, selection.method = 'vst', nfeatures = 2000)
sce_second <- ScaleData(sce_second)
sce_second <- RunPCA(sce_second, features = VariableFeatures(object = sce_second))
sce_second <- RunPCA(object = sce_second, pc.genes = VariableFeatures(sce_second)) 
sce_second <- RunTSNE(sce_second, dims = 1:30)
sce_second <- RunUMAP(sce_second, dims = 1:30)
sce_second <- FindNeighbors(sce_second, dims = 1:15)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce_second <- FindClusters(sce_second, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
clustree(sce_second@meta.data, prefix = "RNA_snn_res.") 

# DimPlot(sce_second, reduction = "tsne", group.by = "RNA_snn_res.0.05", label = T, label.box = T,pt.size = 1)

p_FAP <- DimPlot(sce_second, reduction = "tsne", group.by = "singleR")+scale_color_npg()

p_FAP


cell.use <- row.names(sce_second@meta.data)[which(sce_second@meta.data$singleR == "Epithelial_cells")]

sce_epi <- subset(sce_second, cells=cell.use)

DefaultAssay(sce_epi) <- "RNA"

sce_epi <- NormalizeData(sce_epi, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce_epi <- FindVariableFeatures(sce_epi, selection.method = 'vst', nfeatures = 2000)
sce_epi <- ScaleData(sce_epi)
sce_epi <- RunPCA(sce_epi, features = VariableFeatures(object = sce_epi))
sce_epi <- RunPCA(object = sce_epi, pc.genes = VariableFeatures(sce_epi)) 
sce_epi <- RunTSNE(sce_epi, dims = 1:30)
sce_epi <- RunUMAP(sce_epi, dims = 1:30)
sce_epi <- FindNeighbors(sce_epi, dims = 1:15)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce_epi <- FindClusters(sce_epi, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
clustree(sce_epi@meta.data, prefix = "RNA_snn_res.") 

def_color2 <- c("#FD5149","#80C03A","#023459","#DC9D10","#A27C77","#F23C68","#254840","#7A97AF","#6C69AA","#F0DA87","#C2DDA6","#B6312F","#487388")

p_epi <- DimPlot(sce_epi, reduction = "tsne", group.by = "RNA_snn_res.0.5")+scale_color_manual(values = def_color2)

p_epi

pdf("tsne_CRC_FAP_MAP.pdf",width = 5,height = 4)
p_sample
dev.off()

pdf("tsne_SingleR_all.pdf",width = 6,height = 4)
p_cell
dev.off()

pdf("tsne_SingleR_FAP.pdf",width = 6,height = 4)
p_FAP
dev.off()

pdf("tsne_Epi.pdf",width = 5,height = 4)
p_epi
dev.off()

###################
sce.stat <- table(sce.all@meta.data$singleR) %>% as.data.frame()
sce.stat$Perc <- sce.stat$Freq/sum(sce.stat$Freq)*100

sce.stat <- data.frame()


for (i in 1:6) {
  sce.stat <- rbind(sce.stat,table(sce_epi@meta.data$RNA_snn_res.0.5[which(sce_epi@meta.data$orig.ident == paste0("FAP",i))]) %>% as.data.frame())
}
sce.stat$ident <- rep(paste("FAP",seq(1:6),sep = ""),13)[order(rep(paste("FAP",seq(1:6),sep = ""),13)) ]


p_perc <- ggplot2::ggplot(data = sce.stat,mapping = aes(x=ident,y=Freq,fill=Var1))+
          geom_bar(stat = "identity",width = 0.6,position = position_fill(reverse = T))+ 
          theme(axis.line = element_line(linetype = "solid"),
          panel.grid.major = element_line(linetype = "blank"),
          panel.grid.minor = element_line(linetype = "blank"),
          axis.title = element_text(family = "Helvetica"),
          axis.text = element_text(family = "Helvetica", colour = "black"), panel.background = element_rect(fill = NA),
                      legend.position = "none") +labs(x = NULL, y = NULL)+
          scale_fill_manual(values = def_color2)
############


for (i in c("EPCAM","MUC5AC","S100P","MKI67","SOX4")) {
  p_vln <- VlnPlot(sce_epi, features = i, pt.size = 0 ,group.by = "RNA_snn_res.0.5")+scale_fill_manual(values = def_color2)
  print(p_vln)
} 


############
Idents(sce_epi) <- sce_epi@meta.data$RNA_snn_res.0.5

markers <- FindAllMarkers(sce_epi, only.pos = T, min.pct = 0.25, logfc.threshold = 0.1)

top_TC_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

DoHeatmap(sce_epi, features = top_TC_markers$gene, group.by = "RNA_snn_res.0.5", group.colors = def_color2) +
  scale_fill_viridis()

for (i in top_TC_markers$gene[which(top_TC_markers$cluster == 5)]) {
  
p <-FeaturePlot(sce_epi,features = i,reduction = "tsne",pt.size = 0.4)+
    scale_x_continuous("")+scale_y_continuous("")+
    theme_bw()+ 
    theme( 
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
      axis.ticks = element_blank(),axis.text = element_blank(), 
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5,size=14) 
    )
pdf(file=paste0(i,".pdf"),width = 4,height = 4.2)  
print(p)
dev.off()
}





############
library(monocle)
expr_matrix <- as(as.matrix(sce_epi@assays$RNA@counts),"sparseMatrix")
p_data <- sce_epi@meta.data
f_data <- data.frame(gene_short_name=row.names(sce_epi),row.names = row.names(sce_epi))
pd <- new("AnnotatedDataFrame",data = p_data)
fd <- new("AnnotatedDataFrame",data = f_data)
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())

#cds <- importCDS(sce_epi,import_all=T)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~RNA_snn_res.0.5",cores = 1)
deg <- subset(diff,qval < 0.01)
deg <- deg[order(deg$qval,decreasing = F),]
ordergene <- rownames(deg)[1:800]
cds <- setOrderingFilter(cds,ordergene)
plot_ordering_genes(cds)

cds <-  reduceDimension(cds,max_components = 2,reduction_method = "DDRTree", norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed")
cds <- orderCells(cds)
plot_cell_trajectory(cds,color_by = "Pseudotime")
plot_cell_trajectory(cds,color_by = "RNA_snn_res.0.5")+scale_color_manual(values = def_color2)
                                                                             
#cell_cycling score
pData(cds)$CD44 <- log2(exprs(cds)["CD44",]+1)
plot_cell_trajectory(cds,color_by = "CD44")+scale_colour_gsea()

#Stemness index
pData(cds)$CEACAM7 <- log2(exprs(cds)["CEACAM7",]+1)
plot_cell_trajectory(cds,color_by = "CEACAM7")+scale_colour_gsea()

s_genes <- c("LGR5","SOX9","CD44")
plot_genes_in_pseudotime(cds[s_genes,],color_by = "RNA_snn_res.0.5")+scale_color_manual(values = def_color2)

s_genes <- c("ERBB3","KLF4","PROX1")
plot_genes_in_pseudotime(cds[s_genes,],color_by = "RNA_snn_res.0.5")+scale_color_manual(values = def_color2)

s_genes <- c("PCNA","CCND1","CCND2")
plot_genes_in_pseudotime(cds[s_genes,],color_by = "RNA_snn_res.0.5")+scale_color_manual(values = def_color2)

pdf("pseudotime_heatmap.pdf",width = 6,height = 8)
plot_pseudotime_heatmap(cds[deg$gene_short_name[1:100],],show_rownames = T, hmcols = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(40))
dev.off()

######cellchat
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(SingleR)

cell.use.t <- row.names(sce_second@meta.data)[which(sce_second@meta.data$singleR == "T_cells")]
cell.use.B <- row.names(sce_second@meta.data)[which(sce_second@meta.data$singleR == "B_cell")]
cell.use.M <- row.names(sce_second@meta.data)[which(sce_second@meta.data$singleR == "Macrophage")]

sce_t <- subset(sce_second, cells=cell.use.t)
sce_B <- subset(sce_second, cells=cell.use.B)
sce_M <- subset(sce_second, cells=cell.use.M)

sce.singleR.t <- GetAssayData(sce_t,slot = "data")

clusters.t <- as.factor(as.character(sce_t@meta.data$RNA_snn_res.0.1))

pred.sce.t <- SingleR(test = sce.singleR.t, ref = HGA,
                    labels = HGA$label.fine,
                    clusters = clusters.t)

celltype.t <- data.frame(clusterID = levels(clusters.t),
                       celltype = pred.sce.t$labels)

sce_t@meta.data$singleR <-  celltype.t[match(clusters.t, celltype.t$clusterID),"celltype"]
#sce_t@meta.data$singleR[which(sce_t@meta.data$singleR == "B_cell:Memory")] <- "Treg"

p_t_cell <- DimPlot(sce_t, reduction = "tsne", group.by = "singleR")+scale_color_npg()
p_t_cell

sce.singleR.B <- GetAssayData(sce_B,slot = "data")

clusters.B <- as.factor(as.character(sce_B@meta.data$RNA_snn_res.0.1))

pred.sce.B <- SingleR(test = sce.singleR.B, ref = HGA,
                      labels = HGA$label.fine,
                      clusters = clusters.B)

celltype.B <- data.frame(clusterID = levels(clusters.B),
                         celltype = pred.sce.B$labels)

sce_B@meta.data$singleR <-  celltype.B[match(clusters.B, celltype.B$clusterID),"celltype"]

p_B_cell <- DimPlot(sce_B, reduction = "tsne", group.by = "singleR")+scale_color_npg()
p_B_cell

sce.singleR.M <- GetAssayData(sce_M,slot = "data")

clusters.M <- as.factor(as.character(sce_M@meta.data$RNA_snn_res.0.1))

pred.sce.M <- SingleR(test = sce.singleR.M, ref = HGA,
                      labels = HGA$label.fine,
                      clusters = clusters.M)

celltype.M <- data.frame(clusterID = levels(clusters.M),
                         celltype = pred.sce.M$labels)

sce_M@meta.data$singleR <-  celltype.M[match(clusters.M, celltype.M$clusterID),"celltype"]

p_M_cell <- DimPlot(sce_M, reduction = "tsne", group.by = "singleR")+scale_color_npg()
p_M_cell

#merge
cell.use.isc <- row.names(sce_epi@meta.data)[which(sce_epi@active.ident == "5")]
sce.tme <- subset(sce_epi, cells=cell.use.isc)
sce.tme <- merge(sce.tme,sce_t)
sce.tme <- merge(sce.tme,sce_B)
sce.tme <- merge(sce.tme,sce_M)


#cellchat
cellchat <- createCellChat(sce.tme@assays$RNA@data)
meta <- data.frame(cellType = sce.tme$singleR, row.names =  Cells(sce.tme))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
cellchat <- setIdent(cellchat, ident.use = "cellType")
cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

cellchat <- createCellChat(sce_second@assays$RNA@data)
meta <- data.frame(cellType = sce_second$singleR, row.names =  Cells(sce_second))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
cellchat <- setIdent(cellchat, ident.use = "cellType")
cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")



for (pathways.show in c("CXCL","TGFb")) {
  pdf(paste0(pathways.show,"_signalingrole.pdf"),width = 12,height = 8)
  Role <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  print(Role)
  dev.off()
  
}



ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("outin.cell.pdf",width = 12,height = 8)
ht1 + ht2
dev.off()

####data for CSOmap and cellphoneDB
cso <- as.data.frame(sce.tme@assays$RNA@data)
write.table(cso,"TPM.txt",row.names = T,quote = F,sep = "\t")
cso_label <-data.frame(cellType = sce.tme$singleR, row.names =  Cells(sce.tme))
write.table(cso_label,"label.txt",row.names = T,quote = F,sep = "\t")                       
