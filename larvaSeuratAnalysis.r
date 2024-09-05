library(cowplot)
library(Matrix)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(stringr)
library(plyr)
library(tidyverse)
library(tidyft)
library("ComplexHeatmap")
#library("clustify")

library("limma")
library(data.table)
library("RColorBrewer")
library(scico)

# Plotting
library(ggraph)
library("clustree")
library(viridisLite)
library("viridis")
library(circlize)

sessionInfo()

wd1 <- "yourDirectory"
wd2 <- "DirectoryWithExpressionMatrices"
wd3 <- "outs/filtered_feature_bc_matrix/"    ### use filtered matrices (empty droplets are removed) 

setwd(wd1)
object <- "yourSeuratobject.RData"

ids1 <- c("larva_rep1",  "larva_rep2", "larva_rep3")

data <- sapply(ids1, function(i){
  d10x <- Read10X(file.path(wd2,i,wd3))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

data <- do.call("cbind", data)

larva <- CreateSeuratObject(counts =  data, 
                            project = "larva", 
                            min.cells = 3, 
                            min.features = 200,   
                            names.field = 2,
                            names.delim = "\\-")


larva@meta.data[["orig.ident"]] <- factor(larva@meta.data[["orig.ident"]], levels = ids1)

samplename = larva@meta.data$orig.ident
table(samplename)


###larva_rep1 8014
###larva_rep2 6443
###larva_rep3 6599


# KH annotation (Cao et al, Nature, 2019)

KHclusterID <- read_delim("KHAnnotationFile.txt")
KHclusterID[c("orig.ID","cell.name")] <- str_split_fixed(KHclusterID$NAME, "_", 2)
KHclusterID[c("stage","rep")] <- str_split_fixed(KHclusterID$orig.ID, "\\.", 2)

lvKHclusterID <- KHclusterID[ which(KHclusterID$stage == "lv"), ]

###larva_rep1 is lv.3 #3680 cells
###larva_rep2 is lv.4 #4749 cells
###larva_rep3 is lv.1 #7362 cells

larva_rep1KHclusterID <- lvKHclusterID[ which(lvKHclusterID$rep == 3), ]
larva_rep1KHclusterID <- subset(larva_rep1KHclusterID, select = c("Tissue Type", "cell.name"))

larva_rep2KHclusterID <- lvKHclusterID[ which(lvKHclusterID$rep == 4), ]
larva_rep2KHclusterID <- subset(larva_rep2KHclusterID, select = c("Tissue Type", "cell.name"))

larva_rep3KHclusterID <- lvKHclusterID[ which(lvKHclusterID$rep == 1), ]
larva_rep3KHclusterID <- subset(larva_rep3KHclusterID, select = c("Tissue Type", "cell.name"))

pdf(file="gene-RNA-raw.content-2.pdf", 
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeatureScatter(larva, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


###Low-quality / dying cells often exhibit extensive mitochondrial contamination
###We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features


larva[["percent.mt"]] <- PercentageFeatureSet(larva, pattern = "^mt.")
VlnPlot(larva, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf(file="gene-RNA-raw.content-1.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
VlnPlot(larva, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


VlnPlot(larva, features = "percent.mt") +
  geom_hline(yintercept=15, color = "red", size=2)


###Filtering cells per sample
##larva1 8014 cells

larva1 <- subset(larva, subset = orig.ident == "larva_rep1")
VlnPlot(larva1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

#### nFeature_RNA and nCount_RNA is a bimodal distribution find minimum

## nFeature_RNA
dFeature.lv1 <- density(larva1$nFeature_RNA)

hist(larva1$nFeature_RNA,prob=TRUE)
lines(dFeature.lv1 )
optimize(approxfun(dFeature.lv1 $x,dFeature.lv1 $y),interval=c(200,1200))$minimum

hist(larva1$nFeature_RNA,prob=TRUE)
lines(dFeature.lv1, col="red", lty=2)
vFeature.lv1  <- optimize(approxfun(dFeature.lv1$x,dFeature.lv1$y),interval=c(200,1200))$minimum
abline(v=vFeature.lv1, col="blue")


VlnPlot(larva1, features = "nFeature_RNA") +
  geom_hline(yintercept=optimize(approxfun(dFeature.lv1$x,dFeature.lv1$y),interval=c(200,1200))$minimum, color = "red", size=2)

## nCount_RNA
VlnPlot(larva1, features = "nCount_RNA")

max.nCount.lv1 <- mean(larva1$nCount_RNA) + 5*sd(larva1$nCount_RNA)
VlnPlot(larva1, features = "nCount_RNA") +
  geom_hline(yintercept=max.nCount.lv1, color = "red", size=2)

###Threshold based on nFeature threshold. Then nCount (mean + 5.sd) .
##Then too high mt.precent combined to < minimum nFeature

lv1 <- larva1@meta.data


for (x in 1:nrow(lv1)) 
  {
 
  if(lv1[x, "nFeature_RNA"] < vFeature.lv1)
  {
    lv1[x, "outliner"] <- 1
  }
  else
    if(lv1[x, "nCount_RNA"] > max.nCount.lv1)
    {lv1[x, "outliner"] <- 2}
  else
  {  
    if(lv1[x, "percent.mt"] > 20 | lv1[x, "percent.mt"] < 0.08)
    {lv1[x, "outliner"] <- 3}
    else
      {lv1[x, "outliner"] <- 0}
  }
} 
  

lv1$outliner <- as.character(lv1$outliner)

ggplot(lv1, aes(x = nFeature_RNA, y = nCount_RNA, colour = outliner)) +
  geom_point()


lv1.2 <- lv1[which(lv1$outliner == "0"),]


larva1 <- AddMetaData(larva1, lv1$outliner, col.name = "outliner")

larva1 <- subset(larva1, subset = outliner == "0")
VlnPlot(larva1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 


###After qc subsetting, 3950 cells

#### Larva2
larva2 <- subset(larva, subset = orig.ident == "larva_rep2")
VlnPlot(larva2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

#### nFeature_RNA and nCount_RNA is a bimodal distribution find minimum

## nFeature_RNA
dFeature.lv2 <- density(larva2$nFeature_RNA)

hist(larva2$nFeature_RNA,prob=TRUE)
lines(dFeature.lv2, col="red", lty=2)

vFeature.lv2  <- optimize(approxfun(dFeature.lv2$x,dFeature.lv2$y),interval=c(200,1200))$minimum

VlnPlot(larva2, features = "nFeature_RNA") +
  geom_hline(yintercept=optimize(approxfun(dFeature.lv2$x,dFeature.lv2$y),interval=c(200,1200))$minimum,
             color = "red", size=2)

## nCount_RNA
VlnPlot(larva2, features = "nCount_RNA")

max.nCount.lv2 <- mean(larva2$nCount_RNA) + 5*sd(larva2$nCount_RNA)
VlnPlot(larva2, features = "nCount_RNA") +
  geom_hline(yintercept=max.nCount.lv2, color = "red", size=2)


###Threshold based on nFeature threshold. Then nCount (mean + 5.sd) .
##Then too high mt.precent combined to < minimum nFeature

lv2 <- larva2@meta.data

for (x in 1:nrow(lv2)) 
{
  if(lv2[x, "nFeature_RNA"] < vFeature.lv2)
  {
    lv2[x, "outliner"] <- 1
  }
  else
   { 
     if(lv2[x, "nCount_RNA"] > max.nCount.lv2)
    {
      lv2[x, "outliner"] <- 2
    }
  else
      {
        if(lv2[x, "percent.mt"] > 20 | lv2[x, "percent.mt"] < 0.08)
        {lv2[x, "outliner"] <- 3}
        else
        {lv2[x, "outliner"] <- 0}
        }
      }
     }


lv2$outliner <- as.character(lv2$outliner)
ggplot(lv2, aes(x = nFeature_RNA, y = nCount_RNA, colour = outliner)) +
  geom_point()
lv2.2 <- lv2[which(lv2$outliner == "0"),]
larva2 <- AddMetaData(larva2, lv2$outliner, col.name = "outliner")

larva2 <- subset(larva2, subset = outliner == "0")
VlnPlot(larva2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

### after qc subsetting, 4965 cells

### Larva rep3
larva3 <- subset(larva, subset = orig.ident == "larva_rep3")
VlnPlot(larva3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

#### nFeature_RNA and nCount_RNA is a bimodal distribution find minimum

## nFeature_RNA
dFeature.lv3 <- density(larva3$nFeature_RNA)

hist(larva3$nFeature_RNA,prob=TRUE)
lines(dFeature.lv3, col="red", lty=2)

vFeature.lv3  <- optimize(approxfun(dFeature.lv3$x,dFeature.lv3$y),interval=c(500,1200))$minimum

VlnPlot(larva3, features = "nFeature_RNA") +
  geom_hline(yintercept=optimize(approxfun(dFeature.lv3$x,dFeature.lv3$y),interval=c(500,1200))$minimum,
             color = "red", size=2)

## nCount_RNA
VlnPlot(larva3, features = "nCount_RNA")
max.nCount.lv3 <- mean(larva3$nCount_RNA) + 5*sd(larva3$nCount_RNA)
VlnPlot(larva3, features = "nCount_RNA") +
  geom_hline(yintercept=max.nCount.lv3, color = "red", size=2)


###Threshold based on nFeature threshold. Then nCount (mean + 5.sd) .
##Then too high mt.precent combined to < minimum nFeature

lv3 <- larva3@meta.data

for (x in 1:nrow(lv3)) 
{
  if(lv3[x, "nFeature_RNA"] < vFeature.lv3)
  {
    lv3[x, "outliner"] <- 1
  }
  else
  {
    if(lv3[x, "nCount_RNA"] > max.nCount.lv3)
    {
      lv3[x, "outliner"] <- 2
    }
  else
    
  {
        if(lv3[x, "percent.mt"] > 20 | lv3[x, "percent.mt"] < 0.08)
        {lv3[x, "outliner"] <- 3}
        else
          {lv3[x, "outliner"] <- 0}}
    }
}


lv3$outliner <- as.character(lv3$outliner)
ggplot(lv3, aes(x = nFeature_RNA, y = nCount_RNA, colour = outliner)) +
  geom_point()
lv3.2 <- lv3[which(lv3$outliner == "0"),]
larva3 <- AddMetaData(larva3, lv3$outliner, col.name = "outliner")

larva3 <- subset(larva3, subset = outliner == "0")
VlnPlot(larva3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

larva <- merge(x = larva1, y = list(larva2, larva3))

pdf(file="gene-RNA.content-1.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
VlnPlot(larva, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(file="gene-RNA.content-2.pdf", 
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeatureScatter(larva, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

samplename = larva@meta.data$orig.ident
table(samplename)

### larva_rep1 3950 cells
### larva_rep2 4965 cells
### larva_rep3 5580 cells



# split the dataset into a list of three seurat objects

larva.list <- SplitObject(larva, split.by = "orig.ident")


# normalize and identify variable features for each dataset independently

larva.list <- lapply(X = larva.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})

# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = larva.list, nfeatures = 4000)

save.image(object)

#Perform integration_using cca method

larva.anchors <- FindIntegrationAnchors(object.list = larva.list, anchor.features = features)

save.image(object)

# this command creates an 'integrated' data assay

larva.combined <- IntegrateData(anchorset = larva.anchors)
save.image(object)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay

DefaultAssay(larva.combined) <- "integrated"


# Run the standard workflow for visualization and clustering

larva.combined <- ScaleData(larva.combined, verbose = FALSE)
larva.combined <- RunPCA(larva.combined, npcs = 40, verbose = FALSE)

larva.combined$orig.ident <- factor(x= larva.combined$orig.ident, levels = ids1)

save.image(object)

#DimHeatmap(larva.combined, dims = 1:10, cells = 500, balanced = TRUE)

pdf(file="PCA.pdf", 
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(larva.combined, reduction = "pca")
dev.off()

#Choosing number of PCAs for clustering

pdf(file="Elbowplot.pdf", 
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ElbowPlot(larva.combined, ndims = 40)
dev.off()


save.image(object)

#number of PCA for the analysis 16  

n.dims <- 16

#Project cells into UMAP and tSNE space

larva.combined <- RunTSNE(larva.combined, reduction = "pca", dims = 1:n.dims)

larva.combined <- RunUMAP(larva.combined, reduction = "pca",
                          dims = 1:n.dims, reduction.name = "UMAP")

col3 <- c("#35608D","#66CC33","#E31A1C") 

pdf(file="sampleID-tSNE_.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(larva.combined, reduction = "tsne", 
        group.by = "orig.ident", cols = alpha(col3, 0.5))+
  coord_fixed()
dev.off()

DimPlot(larva.combined, reduction = "tsne", 
        split.by = "orig.ident")

DimPlot(larva.combined, reduction = "UMAP", 
        group.by = "orig.ident", cols = alpha(col3, 0.5), label.size = 2, pt.size = 0.2)+
  coord_fixed()
ggsave("sampleID-UMAP.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm") 

#Perform clustering

larva.combined <- FindNeighbors(larva.combined, reduction = "pca", dims = 1:n.dims)
resolutions <- seq(0.5, 1, 0.1)
larva.combined <- FindClusters(larva.combined,
                               resolution = resolutions,
                               algorithm = 1) #(1 = original Louvain algorithm; 
                                              #2 = Louvain algorithm with multilevel refinement; 
                                              #3 = SLM algorithm; 
                                              #4 = Leiden algorithm). Leiden requires the leidenalg python.


pdf(file="cluster-tree.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
clustree(larva.combined)
dev.off()

save.image(object)

# take resolution 0.7
head(Idents(larva.combined))

larva.combined$seurat_clusters <- larva.combined$integrated_snn_res.0.7
Idents(larva.combined) <- "seurat_clusters"

DimPlot(larva.combined, reduction = "tsne",
        group.by = "integrated_snn_res.0.7", label = TRUE) + NoLegend()
ggsave("tSNE-cluster_res0.7.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm") 

pdf(file="UMAP_seurat-cluster_res0.7.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(larva.combined, reduction = "UMAP", label = TRUE, group.by = "integrated_snn_res.0.8") + NoLegend()
dev.off()

###Find few cell type
DefaultAssay(larva.combined) <- "RNA"  
FeaturePlot(larva.combined, features = "KY21:KY21.Chr12.803",
            reduction = "tsne", max.cutoff = "q90") ##Pigment cells
FeaturePlot(larva.combined, features = "KY21:KY21.Chr1.1601",
            reduction = "tsne", max.cutoff = "q90") ##Muscle


# Identify cluster types (ATTENTION: as to be performed on the RNA assay slot)

DefaultAssay(larva.combined) <- "RNA"  
larva.combined.markers <- FindAllMarkers(larva.combined, only.pos = TRUE,
                                         min.pct = 0.6,
                                         logfc.threshold = 0.5,
                                         test.use = "roc",
                                         slot = "data")

save.image(object)

#annotate markers
human.homo <- read.csv("FileWithGeneAnnotation.csv", header = TRUE, row.names = 1)

human.homo <- as.data.table(human.homo)

human.homo <- human.homo %>% mutate(human.homolog=coalesce(human.homolog,gene))

larva.combined.markers <- tidyft::left_join(larva.combined.markers,human.homo,
                                            by = "gene")
write.table(larva.combined.markers, file="DEG_combined_res0.7.txt",
            quote=F, sep="\t", col.names=NA)

save.image(file = object)

#see which tissue type cells are with KH annotation (Cao et al, Nature, 2019)

metadata <- larva.combined@meta.data

###larva_rep1 is lv.3
###larva_rep2 is lv.4
###larva_rep3 is lv.1

md_larva_rep1 <- metadata[ which(metadata$orig.ident == "larva_rep1"), ]


md_larva_rep1 <- rownames_to_column(md_larva_rep1)
md_larva_rep1[c("cell.name", "orig.ID2")] <-  str_split_fixed(md_larva_rep1$rowname, "-", 2)

md_larva_rep1 <- subset(md_larva_rep1, select = - 15)
md_larva_rep1 <- left_join(md_larva_rep1, larva_rep1KHclusterID, by = "cell.name")

md_larva_rep1 <- column_to_rownames(md_larva_rep1)
md_larva_rep1 <- subset(md_larva_rep1, select = "Tissue Type")


##Rep2
md_larva_rep2 <- metadata[ which(metadata$orig.ident == "larva_rep2"), ]

md_larva_rep2 <- rownames_to_column(md_larva_rep2)
md_larva_rep2[c("cell.name", "orig.ID2")] <-  str_split_fixed(md_larva_rep2$rowname, "-", 2)

md_larva_rep2 <- subset(md_larva_rep2, select = - 15)
md_larva_rep2 <- left_join(md_larva_rep2, larva_rep2KHclusterID, by = "cell.name")

md_larva_rep2 <- column_to_rownames(md_larva_rep2)
md_larva_rep2 <- subset(md_larva_rep2, select = "Tissue Type")

##Rep3
md_larva_rep3 <- metadata[ which(metadata$orig.ident == "larva_rep3"), ]


larva_rep3KHclusterID[c("cell.name2", "orig.ID2")] <-  str_split_fixed(larva_rep3KHclusterID$cell.name, "-", 2)
larva_rep3KHclusterID <- subset(larva_rep3KHclusterID, select = c(1,3))
colnames(larva_rep3KHclusterID)[2]<- "cell.name"

larva_rep3KHclusterID$duplicate.cell <- duplicated(larva_rep3KHclusterID$cell.name)
larva_rep3KHclusterID <- larva_rep3KHclusterID[which(larva_rep3KHclusterID$duplicate.cell == FALSE),]
larva_rep3KHclusterID <- subset(larva_rep3KHclusterID, select = -3)


md_larva_rep3 <- rownames_to_column(md_larva_rep3)
md_larva_rep3[c("cell.name", "orig.ID2")] <-  str_split_fixed(md_larva_rep3$rowname, "-", 2)

md_larva_rep3 <- subset(md_larva_rep3, select = - 15)
md_larva_rep3 <- left_join(md_larva_rep3, larva_rep3KHclusterID, by = "cell.name")

md_larva_rep3 <- column_to_rownames(md_larva_rep3)
md_larva_rep3 <- subset(md_larva_rep3, select = "Tissue Type")

md <- rbind(md_larva_rep1, md_larva_rep2, md_larva_rep3)
colnames(md) <- "KH.Tissue.Type"

larva.combined <- AddMetaData(larva.combined, md)

DimPlot(larva.combined, reduction = "tsne", group.by = "KH.Tissue.Type")
ggsave("tSNE-KH.tissue.type.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm")

save.image(object)


#Assign NS cell types from Cao, et al 2019
#see which tissue type cells are with KH annotation (Cao et al, Nature, 2019)

KHclusterID <- read_delim("KH.CNS.lv.annotation.txt")
KHclusterID[c("orig.ID","cell.name", "CNS")] <- str_split_fixed(KHclusterID$NAME, "_", 3)
KHclusterID[c("stage","rep")] <- str_split_fixed(KHclusterID$orig.ID, "\\.", 2)

metadata <- larva.combined@meta.data

##Rep1
md_larva_rep1 <- metadata[ which(metadata$orig.ident == "larva_rep1"), ]
larva_rep1KHclusterID <- KHclusterID[ which(KHclusterID$rep == 3), ]
larva_rep1KHclusterID <- subset(larva_rep1KHclusterID, select = c("Tissue Type", "cell.name"))
md_larva_rep1 <- rownames_to_column(md_larva_rep1)
md_larva_rep1[c("cell.name", "orig.ID2")] <-  str_split_fixed(md_larva_rep1$rowname, "-", 2)

md_larva_rep1 <- left_join(md_larva_rep1, larva_rep1KHclusterID, by = "cell.name")

md_larva_rep1 <- column_to_rownames(md_larva_rep1)
md_larva_rep1 <- subset(md_larva_rep1, select = "Tissue Type")


##Rep2
md_larva_rep2 <- metadata[ which(metadata$orig.ident == "larva_rep2"), ]
larva_rep2KHclusterID <- KHclusterID[ which(KHclusterID$rep == 4), ]
larva_rep2KHclusterID <- subset(larva_rep2KHclusterID, select = c("Tissue Type", "cell.name"))
md_larva_rep2 <- rownames_to_column(md_larva_rep2)
md_larva_rep2[c("cell.name", "orig.ID2")] <-  str_split_fixed(md_larva_rep2$rowname, "-", 2)

md_larva_rep2 <- left_join(md_larva_rep2, larva_rep2KHclusterID, by = "cell.name")

md_larva_rep2 <- column_to_rownames(md_larva_rep2)
md_larva_rep2 <- subset(md_larva_rep2, select = "Tissue Type")

##Rep3
md_larva_rep3 <- metadata[ which(metadata$orig.ident == "larva_rep3"), ]
larva_rep3KHclusterID <- KHclusterID[ which(KHclusterID$rep == 1), ]
larva_rep3KHclusterID <- subset(larva_rep3KHclusterID, select = c("Tissue Type", "cell.name"))
larva_rep3KHclusterID[c("cell.name","cell.name2")] <- str_split_fixed(larva_rep3KHclusterID$cell.name, "-", 2)
larva_rep3KHclusterID$duplicate <- duplicated(larva_rep3KHclusterID$cell.name) 
larva_rep3KHclusterID <- subset(larva_rep3KHclusterID, select = c("Tissue Type", "cell.name"))
md_larva_rep3 <- rownames_to_column(md_larva_rep3)
md_larva_rep3[c("cell.name", "orig.ID2")] <-  str_split_fixed(md_larva_rep3$rowname, "-", 2)

md_larva_rep3 <- left_join(md_larva_rep3, larva_rep3KHclusterID, by = "cell.name")

md_larva_rep3 <- column_to_rownames(md_larva_rep3)
md_larva_rep3 <- subset(md_larva_rep3, select = "Tissue Type")

md <- rbind(md_larva_rep1, md_larva_rep2, md_larva_rep3)
colnames(md) <- "KH.NS.Type"

larva.combined <- AddMetaData(larva.combined, md)

DimPlot(larva.combined, reduction = "tsne",
        group.by = "KH.NS.Type")
ggsave("tSNE-KH.NS.type.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

save.image(object)

head(Idents(larva.combined))

#investigate gene expression/marker genes

FeaturePlot(larva.combined, features = "KY21:KY21.Chr2.862", 
            reduction = "tsne", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(larva.combined, features = c("KY21:KY21.Chr4.186", #"KH2012:KH.C4.697", #ATM (cluster19)
                                        "KY21:KY21.Chr8.285", # "KH2012:KH.C8.738", #ATM TVCs sensory vesicle
                                        "KY21:KY21.Chr8.449", # "KH2012:KH.C8.796", #Rudder cells
                                         "KY21:KY21.Chr7.170",#"KH2012:KH.L141.73", #muscle Rubber cells (cluster9)
                                         "KY21:KY21.Chr1.8"), #"KH2012:KH.C1.344"), #Stronger in b muscle
            reduction = "tsne", ncol=3,
            min.cutoff = "q10", max.cutoff = "q90")


#Germ cells

col1 = c("lightgrey", "#330066")

pdf(file="Germ cells-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(larva.combined, features = c( "KY21:KY21.Chr11.316", #    "KH2012:KH.C11.383",# Brsk
                                         "KY21:KY21.Chr1.618", #"KH2012:KH.C1.755", #Pem1
                                         "KY21:KY21.Chr8.476", #"KH2012:KH.C8.370"?,
                                        "KY21:KY21.Chr4.1200" , #"KH2012:KH.S852.2", #Pem3
                                        "KY21:KY21.Chr11.697" , #"KH2012:KH.L154.37",#Haus1 Ci-ZF114 // Zfl
                                        "KY21:KY21.Chr2.1118" ,#"KH2012:KH.C2.261"),
                                        "KY21:KY21.Chr1.230", #Pem4
                                        "KY21:KY21.Chr3.749", #Pem5
                                        "KY21:KY21.Chr3.118", #Pem6
                                        "KY21:KY21.Chr2.302", #Pen1
                                        "KY21:KY21.Chr14.518", #Pem12
                                        "KY21:KY21.Chr2.1111", #Pem13
                                        "KY21:KY21.Chr2.1153", #Epb4115
                                        "KY21:KY21.Chr8.1130", #Midn
                                        "KY21:KY21.Chr12.615", #Rev3l
                                        "KY21:KY21.Chr5.966", #Syde2
                                        "KY21:KY21.Chr4.491", # PTP-like
                                        "KY21:KY21.Chr9.1133", #Tdrd-r.a
                                        "KY21:KY21.Chr13.286", #Naa40
                                        "KY21:KY21.Chr4.118"),  #Pem2
            reduction = "tsne", ncol=5,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

FeaturePlot(larva.combined, features = c( "KY21:KY21.Chr11.316", #    "KH2012:KH.C11.383",# Brsk
                                          "KY21:KY21.Chr1.618", #"KH2012:KH.C1.755", #Pem1
                                          "KY21:KY21.Chr8.476", #"KH2012:KH.C8.370"?,
                                          "KY21:KY21.Chr4.1200" , #"KH2012:KH.S852.2", #Pem3
                                          "KY21:KY21.Chr11.697" , #"KH2012:KH.L154.37",#Haus1 Ci-ZF114 // Zfl
                                          "KY21:KY21.Chr2.1118") ,#"KH2012:KH.C2.261"),
            reduction = "tsne", ncol=3, cols = col1,
            min.cutoff = "q05", max.cutoff = "q80")
ggsave("tSNE-germcells.pdf", device= "pdf", width = 48, 
        height = 24, units = "cm")

save.image(object)

##Germ cells in cluster 16

germgenes <- c( "KY21:KY21.Chr11.316", #    "KH2012:KH.C11.383",# Brsk
                "KY21:KY21.Chr1.618", #"KH2012:KH.C1.755", #Pem1
                "KY21:KY21.Chr8.476", #"KH2012:KH.C8.370"?,
                "KY21:KY21.Chr4.1200" , #"KH2012:KH.S852.2", #Pem3
                "KY21:KY21.Chr11.697" , #"KH2012:KH.L154.37",#Haus1 Ci-ZF114 // Zfl
                "KY21:KY21.Chr2.1118") #"KH2012:KH.C2.261")
                
exp.mat <- larva.combined[["RNA"]]@data[germgenes,]
exp.mat <- as.data.frame(exp.mat)

#assign positive cells at 1

md <- rownames_to_column(md)

germ.cells <- as.data.frame( md$rowname)
colnames(germ.cells) <- "cell.name"
germ.cells$germ <- "FALSE"

for (x in 1:nrow(germ.cells))
{
  a <- as.character(germ.cells[x,"cell.name"])
  if (exp.mat["KY21:KY21.Chr11.316", a] > 1 & exp.mat["KY21:KY21.Chr1.618", a] > 1 &
      exp.mat["KY21:KY21.Chr8.476", a] > 1 &  exp.mat["KY21:KY21.Chr4.1200", a] > 1 &
      exp.mat["KY21:KY21.Chr11.697", a] > 1)
  {germ.cells[x,"germ"]<- TRUE}
}

germ.cells <- germ.cells[which(germ.cells$germ == TRUE),]

DimPlot(larva.combined, cells.highlight= list(germ.cells$cell.name),
        cols.highlight = c("darkblue"), reduction = "tsne")

save.image(object)
##rubber cells in cluster 16 

rubbergenes <- c("KY21:KY21.Chr1.1701", #Hand2
                 "KY21:KY21.Chr7.840",
                 "KY21:KY21.Chr14.457",
                 "KY21:KY21.Chr4.737")

FeaturePlot(larva.combined, reduction = "tsne", features = rubbergenes,
            min.cutoff = "q05", max.cutoff = "q80")

exp.mat <- larva.combined[["RNA"]]@data[rubbergenes,]
exp.mat <- as.data.frame(exp.mat)

rubber.cells <- as.data.frame( md$rowname)
colnames(rubber.cells) <- "cell.name"
rubber.cells$rubber <- "FALSE"

for (x in 1:nrow(rubber.cells))
{
  a <- as.character(rubber.cells[x,"cell.name"])
  if (exp.mat["KY21:KY21.Chr1.1701", a] > 1 & exp.mat["KY21:KY21.Chr7.840", a] > 1 &
      exp.mat["KY21:KY21.Chr14.457", a] > 1 &  exp.mat["KY21:KY21.Chr4.737", a] > 1)
  {rubber.cells[x,"rubber"]<- TRUE}
}

rubber.cells <- rubber.cells[which(rubber.cells$rubber == TRUE),]

DimPlot(larva.combined, cells.highlight= list(rubber.cells$cell.name),
        cols.highlight = c("darkblue"), reduction = "tsne")

save.image(object)


#assign tissue types and reassign germ cells, rubber cells

tissue.type1 <- read.csv("tissue.types.csv",header=TRUE)
tissue.type1 <- subset(tissue.type1, select = c("cluster","tissue.type.1","tissue.type.2","tissue.type.3"))
colnames(tissue.type1)<-c("cluster","tissue","cell.type.1", "cell.type.2")

metadata <- larva.combined@meta.data

md1 <- subset(metadata,select = seurat_clusters)
md1 <- rownames_to_column(md1)
head(md1)
class(md1$seurat_clusters)
class(tissue.type1$cluster)
tissue.type1[,"cluster"]<-factor(tissue.type1[,"cluster"])
class(tissue.type1)

md1 <- left_join(md1,tissue.type1, by = c("seurat_clusters" = "cluster"))

head(md1)

#Assign germ cells and rubber cells

for (x in 1:nrow(md1))
{
   if( any(germ.cells$cell.name == md1[x,"rowname"]))
  {
    md1[x,c("tissue", "cell.type.1", "cell.type.2")] <- "germ cells"
   }
  else
  {
    if(any(rubber.cells$cell.name == md1[x,"rowname"]))
       {
         md1[x,"tissue"] <- "nervous system"
         md1[x,"cell.type.1"] <- "rubber cells"
         md1[x,"cell.type.2"] <- "nervous system_8"
    }
  }
    
}

md1 <- column_to_rownames(md1)

#Add metadata
md1 <- md1[,-1]

larva.combined <- AddMetaData(larva.combined, md1)

save.image(object)

#Paired palette
col2epi <-"#62A3CB"       #Blue
col2NS <- c("#6FBE58") #Green
col2endo <- "#EB4445" #Red
col2mes <-"#FDB55F"#Orange
col2TVCmusc <- c("#CAB2D6", "#6A3D9A") #Violet
col2germ<- "#B15928" #Brown
col2noto <- "#FFFF99" #Yellow

barplot(1:2, col = c(col2germ, col2noto))

col2 <- c(col2epi, col2NS, col2endo, col2mes, col2TVCmusc, col2noto, col2germ)

larva.combined$tissue <- factor(x= larva.combined$tissue, levels = c("epidermis",
                                                                             "nervous system",
                                                                             "endoderm",
                                                                             "mesenchyme",
                                                                             "TVCs",
                                                                             "muscle",
                                                                             "notochord",
                                                                             "germ cells"))

DimPlot(larva.combined, reduction = "tsne", 
        group.by = "tissue", label.size = 2, pt.size = 0.2,
        label = FALSE, cols = col2)+
  coord_fixed()
ggsave("tsne.tissue.larva.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm") 
save.image(object)

#Subset Nervous system (based on tissue)
head(larva.combined@meta.data)

Idents(larva.combined) <- "tissue"
head(Idents(larva.combined))

NS <- subset(larva.combined, idents = "nervous system")  
saveRDS(NS, file = "NS.rds")
save.image(object)

DimPlot(NS, reduction = "tsne", 
        group.by = "cell.type.1", label.size = 2, pt.size = 0.2,
        label = FALSE)+
  coord_fixed()
ggsave("tsne.cell.type1NS.larva.pdf", device= "pdf", width = 30, 
       height = 20, units = "cm") 

save.image(object)

#Subset Nervous system (10, 14, 23, 24, 26, 27, 28) with some epidermal cluster (17)
Idents(larva.combined) <- "seurat_clusters"
epiNS <- subset(larva.combined, idents = c(10, 14, 23, 24, 26, 27, 28, 17))  
saveRDS(epiNS, file = "epiNS.rds")

save.image(object)

#Verify if cluster present in dataset
Idents(epiNS) <- "KH.NS.Type"
head(Idents(epiNS))
DimPlot(epiNS, reduction = "tsne")
DimPlot(epiNS, reduction = "tsne", 
        cells.highlight= WhichCells(epiNS, idents = "Six3/6+ pro-aSV"),
        cols.highlight = c("darkblue"))

md2 <- epiNS@meta.data

save.image(object)

### Get cell ID for each cluster and tissue type.

metadata <- larva.combined@meta.data

cellIDtypes <- metadata[,c("tissue", "seurat_clusters")]

write.table(cellIDtypes, file="Table1_tissueTypes-wholeLarvae.csv",
            quote=F, sep=",", col.names=NA)

save.image(object)
