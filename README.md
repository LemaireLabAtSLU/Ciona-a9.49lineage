# Ciona-a9.49lineage
This repository contains the R scripts to re-analyzed the single cell RNA-seq data at the larva stage (GEO accession number GSE131155). The reads of the three larval replicates were remapped on the KY21 gene model using CellRanger 7.0.0. The obtained filtered gene expression matrices are available in the folder GeneExpressionMatrice. The analysis was performed in R version 4.2.3 with Seurat version 4.3.0. The version of the other packages used in the analysis as well as the operating system is indicated in sessionInfo().txt.
No non-standard hardware is required for the analysis.
R and Seurat pacakge installation can be found on the following website https://www.r-project.org/ and https://satijalab.org/seurat/articles/install respectively. The installation should take typically a few minutes.


Demo and intructions for use:
The filtered gene expression matrices should be first analyzed using larvaSeuratAnalysis.r to identify and isolate the nervous system. The file tissue.type.SeuratAnalysis.csv contains the correspondence between the clusters and the identified cell types. The neural cells are then subsetted as a new seurat object and saved as an RDS file. This file is then used to run the script NS.larva.SeuratAnalysis.r. The files tissue.types.NSlarva.csv and tissue.type.subcluster.NSlarva.csv contain the correspondence between the clusters and the identified neural cell types. The files Supp Table1_DEG_Hh2_wilcox.xlsx contains the differentially expressed genes of the Hh2+ ependymal cells. The scripts should run in a few hours.

