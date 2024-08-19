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
library(cowplot)

wd1 <- "yourDirectory"
setwd(wd1)
object <- "yourFileName.RData"

NS <- readRDS( file = "epiNS.rds")          ### include epidermal cluster which might have some NS identity
col1 <- c("#35608D", "#66CC33", "#E31A1C" )

DimPlot(NS, group.by = "seurat_clusters")
DimPlot(NS, reduction = "tsne")

DimPlot(NS, group.by = "KH.NS.Type")

DimPlot(NS, reduction = "tsne", 
        group.by = "KH.NS.Type")
ggsave("tSNE-KH.NS.type.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

Idents(NS) <- "KH.NS.Type"
DimPlot(NS, reduction = "tsne", 
        cells.highlight= WhichCells(NS, idents = "Ci-hh2+ SV"),
        cols.highlight = c("darkblue"))

head(Idents(NS))

DimPlot(NS, reduction = "tsne", 
        cells.highlight= WhichCells(NS, idents = "pigment cells"),
        cols.highlight = c("darkblue"))

save.image(object)

#Clustering NS (start workflow again)
##Find variable features in the subset using the RNA assay
DefaultAssay(NS) <- "RNA" 
NS <- FindVariableFeatures(NS, nfeatures = 4000) 
NS.var <- NS@assays[["RNA"]]@var.features

##Run ScaleData on the integrate assay on the new set of variable features

DefaultAssay(NS) <-"integrated"
all.genes <- rownames(NS)
NS <- ScaleData(NS, features = NS.var)

##Run PCA on the integrated assay using the new set of variable features

NS <- RunPCA(NS, features = NS.var, return.only.var.genes=FALSE,
                reduction.key = "NSpca_",
             reduction.name = "NSPCA", npcs = 40)

pdf(file = "NS.PCAs.pdf",
    width = 8,
    height = 8)
DimPlot(NS, reduction = "NSPCA", group.by = "orig.ident")
dev.off()

pdf(file="NS.ElbowPlot.pdf", 
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ElbowPlot(NS, ndims = 40, reduction = "NSPCA")
dev.off()

NS@reductions[["NSPCA"]]@stdev

save.image(object)

##Run FindNeighbors and FindClusters using the new PC dimensions: pick 24

n.dims <- 24

NS <- FindNeighbors(NS, dims = 1:n.dims, 
                       graph.name = c("NS_nn","NS_snn"), reduction = "NSPCA")


resolutions <- seq(2, 10, 2)
NS <- FindClusters(NS, reduction.type = "NSPCA",
                   dims.use = 1:n.dims, 
                   resolution = resolutions,
                   graph.name = "NS_snn",
                   algorithm = 1)



pdf(file="NS.cluster.tree.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
clustree(NS,  prefix = "NS_snn_res.")
dev.off()

save.image(object)



#Pick resolution 10 (52 clusters)

NS@meta.data[["seurat_clusters"]] <- NS@meta.data[["NS_snn_res.10"]]
metadata <- NS@meta.data

NS <- RunUMAP(NS, reduction = "NSPCA",
              dims = 1:n.dims, reduction.name = "umap")

NS@meta.data[["Whole_larva"]] <- NS@meta.data[["integrated_snn_res.0.7"]]

head(Idents(NS))
Idents(NS) <- "seurat_clusters"

pdf(file="NS.umap_clusters.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(NS, reduction = "umap", 
        label = TRUE) + NoLegend()
dev.off()


DimPlot(NS, reduction = "umap", group.by = "KH.NS.Type")
ggsave("NS.umap-KH.NS.type.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

pdf(file="NS.UMAP_previous clustering.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(NS, reduction = "umap", group.by = "Whole_larva",
        label = TRUE) + NoLegend()
dev.off()

pdf(file="NS.UMAP_tissue.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(NS, reduction = "umap", group.by = "tissue", label = TRUE)
dev.off()

pdf(file="NS.UMAP_orig.ident.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(NS, reduction = "umap", group.by = "orig.ident")
dev.off()

save.image(object)

NS <- RunTSNE(NS, reduction = "NSPCA", dims = 1:n.dims,
              reduction.name = "NStsne")

DimPlot(NS, reduction = "NStsne", 
        label = TRUE) + NoLegend()
ggsave("NS.tsne.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

DimPlot(NS, reduction = "NStsne", group.by = "KH.NS.Type")
ggsave("NS.tsne.KH.NS.Type.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

DimPlot(NS, reduction = "NStsne", group.by = "Whole_larva")
ggsave("NS.tsne.previous clustering.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

#### Identify clusters

DefaultAssay(NS) <- "RNA"  

head(Idents(NS))
Idents(NS) <- "seurat_clusters"

# Identify cluster types (ATTENTION: as to be performed on the RNA assay slot)

NS.markers <- FindAllMarkers(NS, only.pos = TRUE,
                             min.pct = 0.6,
                             logfc.threshold = 0.4,
                             test.use = "roc",
                             slot = "data")

#annotate markers
human.homo <- read.delim("FileWithGeneAnnotation.csv", 
                         header = TRUE, row.names = 1, sep = ",")

human.homo <- as.data.table(human.homo)

larva.combined.markers <- tidyft::left_join(NS.markers,human.homo,
                                            by = "gene")

NS.markers <- tidyft::left_join(NS.markers,human.homo,
                                   by = "gene")
write.table(NS.markers, file="DEG_NS.res10.txt",
            quote=F, sep="\t", col.names=NA)

save.image(object)


### Identify Hh2+ cells
Idents(NS) <- "KH.NS.Type"
DimPlot(NS, reduction = "umap", 
        cells.highlight= WhichCells(NS, idents = "Ci-hh2+ SV"),
        cols.highlight = c("darkblue"))
DimPlot(NS, reduction = "umap", 
        cells.highlight= WhichCells(NS, idents = "pigment cells"),
        cols.highlight = c("darkblue"))

DimPlot(NS, reduction = "NStsne", 
        cells.highlight= WhichCells(NS, idents = "Ci-hh2+ SV"),
        cols.highlight = c("darkblue"))
DimPlot(NS, reduction = "NStsne", 
        cells.highlight= WhichCells(NS, idents = "pigment cells"),
        cols.highlight = c("darkblue"))

metadata <- NS@meta.data

###Marker

FeaturePlot(NS, features = c("KY21:KY21.Chr6.250", "KY21:KY21.Chr14.324", "KY21:KY21.Chr14.325"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### Lens cells

FeaturePlot(NS, features = "KY21:KY21.Chr6.22", 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### Glial cells -> cluster 22 and 33

FeaturePlot(NS, features = c("KY21:KY21.Chr5.389","KY21:KY21.Chr12.803", "KY21:KY21.Chr8.260"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### pigment cells

FeaturePlot(NS, features = c("KY21:KY21.Chr2.793", "KY21:KY21.Chr4.1089", "KY21:KY21.Chr7.1221"), 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90") ##Gabaergic

FeaturePlot(NS, features = "KY21:KY21.Chr3.1172", 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90") ### Glutamimergic

FeaturePlot(NS, features = "KY21:KY21.Chr1.783", 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90") ### Cholinergic

FeaturePlot(NS, features = "KY21:KY21.Chr7.480", 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90") ### FoxP

FeaturePlot(NS, features = "KY21:KY21.Chr6.58", 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90") ### Etr

FeaturePlot(NS, features = "KY21:KY21.Chr1.2012", 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90") ### Tff1

FeaturePlot(NS, features = c("KY21:KY21.Chr11.257", "KY21:KY21.Chr6.340"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### VP RNs

FeaturePlot(NS, features = "KY21:KY21.Chr11.1031", 
            reduction = "umap", min.cutoff = "q10", max.cutoff = "q90") ### Prop

FeaturePlot(NS, features = "KY21:KY21.Chr12.471", 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### VPR RNs

FeaturePlot(NS, features = c("KY21:KY21.Chr13.415", "KY21:KY21.Chr2.892", "KY21:KY21.Chr1.2268"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### PSCs vs PSC-related

FeaturePlot(NS, features = "KY21:KY21.Chr8.704", 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### arx nerve chord

FeaturePlot(NS, features = c("KY21:KY21.Chr3.1500", "KY21:KY21.Chr11.1194"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ### collocyte

FeaturePlot(NS, features = c("KY21:KY21.Chr5.501", "KY21:KY21.Chr7.541"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ###En Ependymal cells

FeaturePlot(NS, features = c("KY21:KY21.Chr3.541", "KY21:KY21.Chr8.1223"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ###Dll+ ANB Pitx ANB

FeaturePlot(NS, features = c("KY21:KY21.Chr9.1166", "KY21:KY21.Chr4.347"), 
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90") ###Pitx ANB

FeaturePlot(NS, features = c("KY21:KY21.Chr1.1326",
                            "KY21:KY21.Chr1.2230",
                             "KY21:KY21.Chr1.783"), ncol = 3,
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90")

FeaturePlot(NS, features = c("KY21:KY21.Chr11.890",
                             "KY21:KY21.Chr9.574",
                             "KY21:KY21.Chr14.837",
                             "KY21:KY21.Chr13.457",
                             "KY21:KY21.Chr10.378",
                             "KY21:KY21.Chr4.865",
                             "KY21:KY21.Chr6.231"), ncol = 3,
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90")



FeaturePlot(NS, features = c("KY21:KY21.Chr1.1766"), ncol = 1,
            reduction = "NStsne", min.cutoff = "q10", max.cutoff = "q90")

DotPlot(NS, features = c("KY21:KY21.Chr6.250", "KY21:KY21.Chr14.324", "KY21:KY21.Chr14.325")) ### Lens cells


### 
md1 <- metadata[,c("seurat_clusters","KH.NS.Type")]
write.table(md1, file="correspondance.txt",
            quote=F, sep="\t", col.names=NA)
save.image(object)

###
Idents(NS) <- "seurat_clusters"
DimPlot(NS, reduction = "NStsne", 
        cells.highlight= WhichCells(NS, idents = 10),
        cols.highlight = c("darkblue"))

DimPlot(NS, reduction = "NStsne", label = TRUE) + NoLegend()

###Assign identified tissue types

tissue.type1 <- read.csv("tissue.types.csv",header=TRUE)
tissue.type1 <- subset(tissue.type1, select = c("cluster","tissue.1"))
colnames(tissue.type1)<-c("cluster", "cell.type.3")

metadata <- NS@meta.data

md2 <- subset(metadata,select = seurat_clusters)
md2 <- rownames_to_column(md2)
head(md2)
class(md2$seurat_clusters)
class(tissue.type1$cluster)

md2 <- left_join(md2,tissue.type1, by = c("seurat_clusters" = "cluster"))

head(md2)

md2 <- column_to_rownames(md2)

md2 <- subset(md2, select = "cell.type.3")

NS <-  AddMetaData(NS, md2)

DimPlot(NS, reduction = "NStsne",
        group.by = "cell.type.3")
ggsave("NS.tsne.cell.type.3.pdf", device= "pdf", width = 40, 
      height = 20, units = "cm")

### subcluster
## 7-8-12-13-41
## Not 15, leave as CESNs

NS$subcluster <- NULL

Idents(NS) <- "seurat_clusters"

head(Idents(NS))

NS <- FindSubCluster(NS,
  7,
  graph.name = "NS_snn",
  subcluster.name = "sub.cluster7",
  resolution = 1,
  algorithm = 1)

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("7")),
        group.by = "sub.cluster7")

FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("7")),
            features = c("KY21:KY21.Chr11.1031","KY21:KY21.Chr11.257"))


Idents(NS) <- "sub.cluster7"

NS <- FindSubCluster(NS,
                     8,
                     graph.name = "NS_snn",
                     subcluster.name = "sub.cluster8",
                     resolution = 1,
                     algorithm = 1)

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("8")),
        group.by = "sub.cluster8")

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("8")),
        group.by = "KH.NS.Type")

FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("8")),
            features = c("KY21:KY21.Chr5.641","KY21:KY21.Chr5.389","KY21:KY21.Chr8.260", "KY21:KY21.Chr12.803"),
            min.cutoff = "q10", max.cutoff = "q90")

Idents(NS) <- "sub.cluster8"

NS <- FindSubCluster(NS,
                     "8_1",
                     graph.name = "NS_snn",
                     subcluster.name = "sub.cluster8.1",
                     resolution = 1,
                     algorithm = 1)

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("8_0", "8_1")),
        group.by = "sub.cluster8.1")

Idents(NS) <- "sub.cluster8.1"

NS <- FindSubCluster(NS,
                     12,
                     graph.name = "NS_snn",
                     subcluster.name = "sub.cluster12",
                     resolution = 1,
                     algorithm = 1)

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("12")),
        group.by = "sub.cluster12")
DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("12")),
        group.by = "KH.NS.Type")

Idents(NS) <- "sub.cluster12"

NS <- FindSubCluster(NS,
                     13,
                     graph.name = "NS_snn",
                     subcluster.name = "sub.cluster13",
                     resolution = 0.8,
                     algorithm = 1)

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13")),
        group.by = "sub.cluster13")

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13")),
        group.by = "KH.NS.Type")

FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13")),
            features = c("KY21:KY21.Chr14.324",
                         "KY21:KY21.Chr3.430",
                         "KY21:KY21.Chr14.325",
                         "KY21:KY21.Chr14.267",
                         "KY21:KY21.Chr12.803",
                         "KY21:KY21.Chr8.260",
                         "KY21:KY21.Chr5.389"),
            min.cutoff = "q10", max.cutoff = "q90")

Idents(NS) <- "sub.cluster13"

NS <- FindSubCluster(NS,
                     "13_1",
                     graph.name = "NS_snn",
                     subcluster.name = "sub.cluster13.1",
                     resolution = 1,
                     algorithm = 1)

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13_0","13_1")),
        group.by = "sub.cluster13.1")


Idents(NS) <- "sub.cluster13.1"

NS <- FindSubCluster(NS,
                     41,
                     graph.name = "NS_snn",
                     subcluster.name = "sub.cluster",
                     resolution = 1,
                     algorithm = 1)

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("41")),
        group.by = "sub.cluster")

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("41")),
        group.by = "KH.NS.Type")

save.image(object)


##Find marker genes of subcluster

DimPlot(NS, reduction = "NStsne", group.by = "sub.cluster")

ids <- c("7_0","7_1", "7_2", "8_0", "8_1_0", "8_1_1", "12_0", "12_1", "12_2",
         "13_0", "13_1_0", "13_1_1", "41_0", "41_1")

Idents(NS) <- "sub.cluster"
head(Idents(NS))

marker <- lapply(ids, function(i){
  subclusterMarker <- FindMarkers(NS, ident.1 = i, 
                                  only.pos = TRUE,
                                  min.pct = 0.5,
                                  logfc.threshold = 0.4,
                                  test.use = "roc",
                                  slot = "data")
  subclusterMarker <- rownames_to_column(subclusterMarker)
  colnames(subclusterMarker)[1] <- "gene"
  subclusterMarker[,"cluster"] <- i
  subclusterMarker
})

marker <- do.call("rbind", marker)

marker <- tidyft::left_join(marker,human.homo,
                                by = "gene")
write.table(marker, file="DEG_NS_subcluster.txt",
            quote=F, sep="\t", col.names=NA)

save.image(object)

#### correspondance between subcluster and cell types

metadata <- NS@meta.data
md3 <- subset(metadata, select = c("sub.cluster", "KH.NS.Type"))

write.table(md3, file="correspondance_KH_subcluster.txt",
            quote=F, sep="\t", col.names=NA)


### Assign identity

Idents(NS) <-"seurat_clusters"
DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = 12),
        group.by = "sub.cluster")

Idents(NS) <- "sub.cluster"
DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("12_0", "12_1", "12_2")),
        group.by = "KH.NS.Type")

FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("12_0", "12_1", "12_2")),
            features = c("KY21:KY21.Chr14.324",
                         "KY21:KY21.Chr3.430",
                         "KY21:KY21.Chr14.325",
                         "KY21:KY21.Chr14.267",
                         "KY21:KY21.Chr6.250",
                         "KY21:KY21.Chr5.641"), ncol = 3,
            min.cutoff = "q10", max.cutoff = "q90")


DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13_0", "13_1_0", "13_1_1")),
        group.by = "sub.cluster")

FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13_0", "13_1_0", "13_1_1")),
            features = c("KY21:KY21.Chr14.324",
                         "KY21:KY21.Chr3.430",
                         "KY21:KY21.Chr14.325",
                         "KY21:KY21.Chr14.267",
                         "KY21:KY21.Chr6.250",
              "KY21:KY21.Chr12.306", "KY21:KY21.Chr12.803",
                         "KY21:KY21.Chr5.389", "KY21:KY21.Chr8.260"), ncol = 3,
            min.cutoff = "q10", max.cutoff = "q90")



DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("41_0", "41_1")),
        group.by = "sub.cluster")
FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("41_0", "41_1")),
            features = c("KY21:KY21.Chr4.1089",
                         "KY21:KY21.Chr11.341"), ncol = 2,
            min.cutoff = "q10", max.cutoff = "q90")

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("7_0", "7_1", "7_2")),
        group.by = "sub.cluster")
FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("7_0", "7_1", "7_2")),
            features = c("KY21:KY21.Chr4.1089",
                         "KY21:KY21.Chr11.341",
                         "KY21:KY21.Chr14.647",
                         "KY21:KY21.Chr11.1031",
                         "KY21:KY21.Chr11.1275"), ncol = 2,
            min.cutoff = "q10", max.cutoff = "q90")

DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("8_0", "8_1_0", "8_1_1")),
        group.by = "sub.cluster")

FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("8_0", "8_1_0", "8_1_1")),
            features = c("KY21:KY21.Chr14.324",
                         "KY21:KY21.Chr3.430",
                         "KY21:KY21.Chr14.325",
                         "KY21:KY21.Chr14.267",
                         "KY21:KY21.Chr6.250",
                         "KY21:KY21.Chr5.641", "KY21:KY21.Chr12.803",
                         "KY21:KY21.Chr5.389", "KY21:KY21.Chr8.260"), ncol = 3,
            min.cutoff = "q10", max.cutoff = "q90")


FeaturePlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13_0", "13_1_0", "13_1_1", "8_0", "8_1_0", "8_1_1")),
            features = c("KY21:KY21.Chr1.2268",
                         "KY21:KY21.Chr2.1075"), ncol = 2,
            min.cutoff = "q10", max.cutoff = "q90")
DimPlot(NS, reduction = "NStsne", cells = WhichCells(NS, idents = c("13_0", "13_1_0", "13_1_1", "8_0", "8_1_0", "8_1_1")),
            group.by = "sub.cluster")

save.image(object)

#assign cell type subcluster

tissue.type2 <- read.csv("tissue.type.subcluster.csv",header=TRUE)
tissue.type2 <- subset(tissue.type2, select = c("cluster","cell.type.4"))

metadata <- NS@meta.data

md4 <- subset(metadata,select = sub.cluster)
md4 <- rownames_to_column(md4)
head(md4)
class(md4$sub.cluster)
class(tissue.type2$cluster)

md4 <- left_join(md4,tissue.type2, by = c("sub.cluster" = "cluster"))

head(md4)

md4 <- column_to_rownames(md4)

md4 <- subset(md4, select = "cell.type.4")

NS <-  AddMetaData(NS, md4)

DimPlot(NS, reduction = "NStsne",
        group.by = "cell.type.4")

metadata <- NS@meta.data

for (x in 1:nrow(metadata)) 
{
  if(any(ids == metadata[x,"sub.cluster"]))
  {
    metadata[x,"cell.type.5"] <- metadata[x,"cell.type.4"]
  }

else
{
  metadata[x,"cell.type.5"] <- metadata[x,"cell.type.3"]
}
}

## Modify CESNs + pATENs to CESNs

for (x in 1:nrow(metadata)) 
{
  a <- metadata[x,"cell.type.5"]
  if(a == "CESNs + pATENs")
  {metadata[x,"cell.type.5"] <- "CESNs"}
}

NS <- AddMetaData(NS, metadata$cell.type.5, col.name = "cell.type.5")
DimPlot(NS, reduction = "NStsne",
        group.by = "cell.type.5")
ggsave("NS.tsne.cell.type.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

DimPlot(NS, reduction = "umap",
        group.by = "cell.type.5")
ggsave("NS.umap.cell.type.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

save.image(object)

### Remove endoderm, glia cells and mesenchyme

Idents(NS) <- "cell.type.5"
head(Idents(NS))

NS.2 <- subset(NS, idents =  c("endoderm", "Glia cells", "mesenchyme"), invert = TRUE)
head(Idents(NS.2))


DimPlot(NS.2, reduction = "NStsne",
        group.by = "cell.type.5")
ggsave("NS.no.conta.tsne.cell.type.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

DimPlot(NS.2, reduction = "NStsne", cells.highlight= list(WhichCells(NS.2, idents = "pigment cells"),
                                                          WhichCells(NS.2, idents =  "Hh2+ cells")),
cols.highlight = c("darkblue", "darkgreen"))
ggsave("NS.no.conta.pigment.HH2.cell.type.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

cell.type.order <- NS.2@meta.data$cell.type.5

cell.type.order <- unique(cell.type.order)

cell.type.order <- sort(cell.type.order)

NS.2$cell.type.5 <- factor(x= NS.2$cell.type.5, levels = cell.type.order)

### Identify Hh2+ cells highly expressed genes 

Hh2.marker <- FindMarkers(NS.2, ident.1 = "Hh2+ cells",
                    min.pct = 0.25,
                    logfc.threshold = 0.5, only.pos = TRUE)

Hh2.marker <- rownames_to_column(Hh2.marker)
colnames(Hh2.marker)[1]<-"gene"

Hh2.marker <- rownames_to_column(Hh2.marker)
colnames(Hh2.marker)[1]<-"DEG.order"

Hh2.marker <- tidyft::left_join(Hh2.marker,human.homo,
                            by = ("gene"))
write.table(Hh2.marker, file="DEG_Hh2_wilcox.txt",
            quote=F, sep="\t", col.names=NA)

Hh2.marker$DEG.order <- as.integer(Hh2.marker$DEG.order)

top20 <- Hh2.marker[which(Hh2.marker$DEG.order <21),]
top20 <- top20[order(top20$DEG.order),]

top20gene <- top20$gene

DotPlot(NS.2, features = rev(top20gene), cols = c("lightgrey", "#330066"), cluster.idents = FALSE) +coord_flip()+
  #scale_x_discrete(position = "top")+
  scale_y_discrete(limits = cell.type.order)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=2/3)+
  labs( y = "" , x= "")
ggsave("dotplot.top20.Hh2cells.alphabetic.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

DotPlot(NS.2, features = rev(top20gene), cols = c("lightgrey", "#330066"), cluster.idents = TRUE) +coord_flip()+
  #scale_x_discrete(position = "top")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=2/3)+
  labs( y = "" , x= "")
ggsave("dotplot.top20.Hh2cells.cluster.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

top10 <- Hh2.marker[which(Hh2.marker$DEG.order <11),]
top10 <- top10[order(top10$DEG.order),]

top10gene <- top10$gene


DotPlot(NS.2, features = rev(top10gene), cols = c("lightgrey", "#330066"), cluster.idents = TRUE) +coord_flip()+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/3)+
  labs( y = "" , x= "")
ggsave("dotplot.top10.Hh2cells.cluster.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

DotPlot(NS.2, features = rev(top10gene), cols = c("lightgrey", "#330066"), cluster.idents = FALSE) +coord_flip()+
  scale_y_discrete(limits = cell.type.order)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/3)+
  labs( y = "" , x= "")
ggsave("dotplot.top10.Hh2cells.alphabetic.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

save.image(object)


##Pigment cells most differentially expressed genes

pigment.marker <- FindMarkers(NS.2, ident.1 = "pigment cells",
                          min.pct = 0.25,
                          logfc.threshold = 0.5, only.pos = TRUE)

pigment.marker <- rownames_to_column(pigment.marker)
colnames(pigment.marker)[1]<-"gene"

pigment.marker <- rownames_to_column(pigment.marker)
colnames(pigment.marker)[1]<-"DEG.order"

pigment.marker <- tidyft::left_join(pigment.marker,human.homo,
                                by = ("gene"))
write.table(pigment.marker, file="DEG_pigment_wilcox.txt",
            quote=F, sep="\t", col.names=NA)

pigment.marker$DEG.order <- as.integer(pigment.marker$DEG.order)

pigment.top20 <- pigment.marker[which(pigment.marker$DEG.order <21),]
pigment.top20 <- pigment.top20[order(pigment.top20$DEG.order),]

pigment.top20gene <- pigment.top20$gene

DotPlot(NS.2, features = rev(pigment.top20gene), cols = c("lightgrey", "#330066"), cluster.idents = TRUE) +coord_flip()+
  #scale_x_discrete(position = "top")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=2/3)+
  labs( y = "" , x= "")
ggsave("dotplot.top20.pigmentcells.cluster.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

save.image(object)

##Number of cells in sensory vesicle

head(Idents(NS.2))

SVtypes <- c("coronet cells","FoxP RNs", "Gabaergic interneurons", "Hh2+ cells" ,
             "lens cells",  "OACCs", "pigment cells",            "Rx+ PRCs",
             "Six3/6 pro-ASV", "Stum+ PRCs", "Eminens",
            "switch neurons" , "VP RNs"  ,  "VPR+ RNs"    )

SV <- subset(NS.2, idents = SVtypes) 

save.image(object)
