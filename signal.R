system("ls")
system("mkdir -p Signaling/{Data,Images,Output}")
system("pwd")
system("./Cell_Signaling/Cell_Signaling/fetch_data.sh")

#################### Installations ################
install.packages("anndata")
install.packages("Seurat")
install.packages("BiocManager")
install.packages("scran")
install.packages("SingleCellExperiment")
install.packages("Matrix")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("SingleCellSignalR")
devtools::install_github('satijalab/seurat-data')
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("sqjin/CellChat")
BiocManager::install("topGo")
BiocManager::install("biomaRt")
install.packages("tibble")


################# Load Libraries #####################


library(anndata)
library(Seurat)
library(SingleCellSignalR)
library(SingleCellExperiment)
library(Matrix)
library(methods)
library(CellChat)
library(patchwork)
library(dplyr)
library(topGO)
library(tibble)
library(biomaRt)
library(tidyr)

################ Load Data ##################

kidney_ann <- anndata::read_h5ad("Signaling/Data/local.h5ad")
data_rds <- readRDS("Signaling/Data/local.rds")

### View anndata data  and rds data ######

kidney_ann_data <- kidney_ann$X
kidney_count_data <- data_rds@assays$RNA@counts

head(kidney_ann_data)

str(data_rds)
summary(data_rds)
head(data_rds)
names(data_rds)

colnames(x = data_rds)
rownames(x = data_rds)
Cells(x = data_rds)
Idents(object = data_rds)
levels(x = data_rds)

#### Seurat Object #################

seurat_obj <- CreateSeuratObject(counts = kidney_count_data)

####################################### Pre-Processing PROCESS #########################

data_rds_copy <- data_rds

data_rds_copy[["old.ident"]] <- Idents(object = data_rds_copy)
data_rds_copy <- StashIdent(object = data_rds_copy, save.name = "old.ident")
 
### Set identity classes to cell type
Idents(object = data_rds_copy) <- "cell_type"


feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
pdf(file = "Signaling/Images/Feature_Cell_type.pdf", width = 15, height = 10)
VlnPlot(data_rds_copy, group.by = "cell_type", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()
dev.off() 

plot1 <- FeatureScatter(data_rds_copy, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_rds_copy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file = "Signaling/Images/Features.pdf",width = 15, height = 10)
plot1 + plot2
dev.off()

pdf(file = "Signaling/Images/FeaturesScatter.pdf",width = 15, height = 10)
FeatureScatter(data_rds_copy, "nCount_RNA", "nFeature_RNA", group.by = "cell_type", pt.size = 1)
dev.off()

data_rds_copy <- NormalizeData(data_rds_copy)
data_rds_copy <- FindVariableFeatures(data_rds_copy, selection.method = "vst")


top10 <- head(VariableFeatures(data_rds_copy), 10)
plot1 <- VariableFeaturePlot(data_rds_copy)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "Signaling/Images/variableFeatures.pdf",width = 15, height = 10)
plot1 + plot2
dev.off()

all.genes <- rownames(data_rds_copy)
data_rds_copy <- ScaleData(data_rds_copy)


data_rds_copy <- RunPCA(data_rds_copy, features = VariableFeatures(object = data_rds_copy))
print(data_rds_copy[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file = "Signaling/Images/PCA_1.pdf",width = 15, height = 10)
VizDimLoadings(data_rds_copy, dims = 1:5, reduction = "pca")
dev.off()
pdf(file = "Signaling/Images/PCA_2.pdf",width = 15, height = 10)
DimPlot(data_rds_copy, reduction = "pca")
dev.off()
pdf(file = "Signaling/Images/PCA_HM.pdf",width = 10, height = 8)
DimHeatmap(data_rds_copy, dims = 1:5, cells = 5, balanced = TRUE)
dev.off()

data_rds_copy <- FindNeighbors(data_rds_copy)
data_rds_copy <- FindClusters(data_rds_copy, resolution = .5)


data_rds_copy.markers <- FindAllMarkers(data_rds_copy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(data_rds_copy.markers, n = 10)
data_rds_copy.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#get the gene_name from ensemble ID
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=top10,mart= mart)

pdf(file = "Signaling/Images/Marker_10.pdf",width = 10, height = 8)
FeaturePlot(data_rds_copy, features = top10)
dev.off()

rownames(data_rds_copy@assays$RNA@counts) <- data_rds_copy@assays$RNA@meta.features$feature_name
rownames(data_rds_copy@assays$RNA@data) <- data_rds_copy@assays$RNA@meta.features$feature_name

data_rds_copy <- RunUMAP(data_rds_copy, dims = 1:10, nfeatures = 10)
pdf(file = "Signaling/Images/UMAP.pdf",width = 10, height = 8)
DimPlot(data_rds_copy, reduction = "umap",group.by = 'cell_type')
dev.off()
pdf(file = "Signaling/Images/Gene_List.pdf",width = 10, height = 8)
FeaturePlot(data_rds_copy, features = c(G_list$hgnc_symbol))
dev.off()

plot <- FeaturePlot(data_rds_copy, features = 'APOE')
Hplot <- HoverLocator(plot=plot, infromation=FetchData(data_rds_copy, vars=data_rds_copy@meta.data$cell_type))
htmlwidgets::saveWidget(widget = Hplot, file = 'Signaling/Images/APOE.html', selfcontained = TRUE)
#######################################################################

######################### CellChat ##########################################

genes <- c(rownames(data_rds_copy@assays$RNA@data))
data.input <- data_rds_copy@assays$RNA@data # normalized data matrix
Idents(data_rds_copy) <- 'cell_type'
labels <- Idents(data_rds_copy)
labels
meta <- data.frame(group = labels, row.names = names(labels))
head(meta)# create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use


cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)

##In case of error##
#length(cellchat@idents)
#dim(cellchat@data.signaling)
#unique(cellchat@idents)
#cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))


cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,weight.scale = T, label.edge = F, title.name = "Interactions")

setwd("Signaling/Images")

pathways.show.all <- cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  vertex.receiver = seq(1,5)
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_info.pdf"), plot=gg, width = 5, height = 8, units = 'in', dpi = 300)
}

protein = c("EGF","KIT")
for (i in 1:length(protein)){
  pathways.show <- protein[i] 
  vertex.receiver = seq(1,4) 
  par(mfrow=c(1,1))
  png(file = paste0(protein[i] ,".png"), width = 1500, height =  1000)
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = c("hierarchy"))
  dev.off()
  png(file = paste0(protein[i] ,"_circle.png"), width = 800, height =  1000)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  # Heatmap for analysing interactions
  png(file = paste0(protein[i] ,"_HM.png"),width = 1200, height = 800)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
}

png(file = "KIT_EGF.png")
netVisual_bubble(cellchat, signaling = c("KIT","EGF"))
dev.off()


## GET all the pairs ###
ligands <- unique(cellchat@LR$LRsig$ligand)
typeof(ligands)
write.csv(ligands, file = 'ligands.txt',row.names = FALSE)
receptors <- unique(cellchat@LR$LRsig$receptor)
write.csv(receptors, file = 'receptors.txt',row.names = FALSE)
LR_pair.all <- extractEnrichedLR(cellchat, signaling = pathways.show.all, geneLR.return = FALSE)
interactions <- LR_pair.all$interaction_name
write.csv(interactions, file = 'LR_interactions.txt',row.names = FALSE)

####################### SingleCellSignalR code ###################

data_scsr <- subset(seurat_obj, subset = nFeature_RNA > 1000)
data_scsr <- NormalizeData(data_scsr)
data_scsr <- FindVariableFeatures(data_scsr, selection.method = "vst", nfeatures = 1000)


genes <- rownames(data_scsr)
data_scsr <- ScaleData(data_scsr)
data_scsr <- RunPCA(data_scsr, features = VariableFeatures(object = data_scsr))
data_scsr <- FindNeighbors(data_scsr)
data_scsr <- FindClusters(data_scsr, resolution = 0.5)


cluster = as.numeric(Idents(data_scsr))
data = data.frame(data_scsr[["RNA"]]@data)
gc()
signal= cell_signaling(data=data,genes=genes,cluster=cluster)
gc()
visualize_interactions(signal)
intra = intra_network(data = data_chunck_10, signal = signal, genes = genes, cluster = cluster)
gc()

###No results---


#########################(Optional) kidney data set for SingelCellSignalR code (Alternative)################### 

cells_list <- list(rownames(kidney_ann_data))
genes <- colnames(kidney_ann_data)

cell_id_list <- list(kidney$obs["cell_type_ontology_term_id"][[1]])
cell_id_list[[1]]
data_mat <- data.frame(matrix(nrow = nrow(kidney_ann_data), ncol = ncol(kidney_ann_data)+1))

data_mat['cell_ids'] <- c(cell_id_list[[1]])
colnames(data_mat) <- c('cell_ids', genes)
length(colnames(data_mat))
data_mat[1:5,1:5]

data_mat[1,1]

## access the dgRMatrix format of anndata ???

kidney_data_new <- as(kidney_ann_data, "CsparseMatrix")

############## HubMAP Data ##########
hubmap_gene_data_list <- character()

hubmap_gene_list <- c('Foxd1','NPHS2','NPHS2','AQP2','SLC12A1','CRYAB','SERPINE2','AQP2','NPHS2','SERPINE2','VCAM1','NPHS2','EHD3','POSTN','CUBN','SLC5A12','SLC5A11','SLC22A7','CRYAB','CRYAB','CRYAB','CRYAB','CRYAB','SLC12A1','SLC12A1','SLC12A1','SLC12A1','SLC12A3','SLC12A3','SLC12A3','SLC8A1','AQP2','CA2','SLC4A1','SLC26A4','GATA3','AQP2','AQP2','AQP2','AQP2','FXYD4','CA2','SLC4A1','SLC26A7','SLC26A7','SLC4A9','SLC4A9','EMCN','SERPINE2','SERPINE2','RAMP3','MCTP1','NR2F2','MMRN1','TAGLN','MYH11','REN','DCN','DCN','TAGLN','CD163','ITGAL','NKG7','NKG7','APOE','ITGAX','IL3RA','C1QA','MSR1','CD96','GZMA','JCHAIN','MZB1','S100A9','KIT','KRT5')
hubmap_gene_list <- unique(hubmap_gene_list)
hubmap_cell_list <- c('CL:0000653','CL:0000653','CL:1000714','CL:1001108','CL:1001285','CL:1000718','CL:0000653','CL:1001285','CL:1000452','CL:0000653','CL:1001005','CL:1000742','CL:0002306','CL:4030009','CL:4030010 ','CL:4030011','CL:1001111','CL:4030012','CL:4030013','CL:4030014','CL:1001107','CL:1001106','CL:1001109','CL:1001108','CL:1000850','CL:1000849','CL:4030016','CL:4030017','CL:1000768','CL:4030018','CL:4030019','CL:4030020','CL:4030021','CL:1000454','CL:1001431','CL:1000714','CL:1000716','CL:1000718','CL:1001432','CL:4030015','CL:1000715','CL:1000717','CL:0002201','CL:1000715','CL:0000115','CL:1001096','CL:1001099','CL:1001033','CL:1001285','CL:1001131','CL:0002138','CL:1001318','CL:0000648','CL:1000692','CL:4030022','CL:1000691','CL:1000698','CL:0000875','CL:0000623','CL:0000814','CL:0001056','CL:0000990','CL:0001058','CL:0000576','CL:0000084','CL:0000910','CL:0000236','CL:0000786','CL:0000775','CL:0000097','CL:1000597')

X_matrix <- as.matrix(kidney$var["feature_name"])
X_matrix <- data.matrix(X_matrix)

X_matrix[1]

for (i in 1:length(hubmap_gene_list)){
  for (j in 1:length(X_matrix)){
    print(hubmap_gene_list[i])
    print(X_matrix[j])
    if (hubmap_gene_list[i] != X_matrix[j]) {
      next
    } else {
      hubmap_gene_data_list <- c(hubmap_gene_data_list, rownames(kidney$var["feature_name"][1])[j])
    }
  }
}


length(hubmap_gene_data_list)


## fill the new view data set with using the kidney data:
kidney_data_new <- as.matrix(kidney_data_new)
row.names(kidney_data_new) <- rownames(kidney$X)
colnames(kidney_data_new) <- colnames(kidney$X)

kidney_data_new <- data.frame(kidney_data_new)

for (i in 1:nrow(kidney_data_new)){
  for (j in 1:ncol(kidney_data_new)){
    if (kidney_data_new[i,j] != 0 ){
      ele <- kidney_data_new[i,j]
      col_name <- colnames(kidney_data_new)[j]
      print(col_name)
      if ((col_name %in% hubmap_gene_data_list) == TRUE){
        data_mat[i,j+1] = kidney_data_new[i,j]
        print(data_mat[i,j+1])
      }
    }
  }
}


### free up space ####
gc()
rm()
