## Clustering analysis for E11.5-E15.5 cells with more than 6,000 expressed genes 
## and identification of marker genes of each cell cluster.

library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)

# Load the data
mat <- read.table("all.tube.sc.rpkm.txt",header = T,row.names = 1)
cell_anno <- read.table("mus_all_sc_anno.txt",header = T,row.names = 1,sep = "\t")

# Get Foxa2+ cells from E11.5-E15.5 embryo
cell_anno <- cell_anno[(cell_anno$Stage=="E11.5"|cell_anno$Stage=="E12.5"|cell_anno$Stage=="E13.5"|cell_anno$Stage=="E14.5"|cell_anno$Stage=="E15.5"|cell_anno$Stage=="E14.5_GFPn")&cell_anno$num_genes_expressed>6000,]
ref <- read.table("Mus_ref.txt",header=T,row.names=1,sep = "\t")
pos <- rownames(mat) %in% rownames(ref)
data <- mat[pos,rownames(cell_anno)]

# Create Seurat object
mus <- CreateSeuratObject(raw.data=data,min.cells = 5, min.genes=1000, is.expr =1,project="Mus_cluster")
mus <- NormalizeData(object = mus,normalization.method = "LogNormalize",scale.factor = 10000)
ref <- ref[!grepl("ERCC-",rownames(ref)), ]
ref <- ref[grepl("protein_coding|^IG_|^TR_",ref$Biotpye), ]
ref <- ref[!duplicated(ref$GeneName),]
pos <- rownames(mus@data) %in% rownames(ref)
mus@data <- mus@data[pos,]
rownames(mus@data) <-  ref[rownames(mus@data),5]

# Identification of highly variable genes
mus <- FindVariableGenes(object = mus, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.5)
length(mus@var.genes)

# Scaling the data
mus <- ScaleData(object = mus)

# Perform linear dimensional reduction
mus <- RunPCA(object = mus, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
mus <- ProjectPCA(object = mus, do.print = FALSE)
PCAPlot(object = mus, dim.1 = 1, dim.2 = 2,pt.size = 1)
PCHeatmap(mus, pc.use = 1:15, cells.use = 400, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCElbowPlot(mus)

# Cluster the cells
mus <- FindClusters(mus, reduction.type = "pca",dims.use = 1:8, save.SNN = T,resolution = 0.7,temp.file.location="/hwfssz1/ST_PRECISION/USER/zhongyu/1.Project/Mus_Liver/seurat/results/")
PrintFindClustersParams(object = mus)

# Run non-linear dimensional reduction (tSNE)
mus <- RunTSNE(mus, dims.use = 1:8, do.fast = T)
TSNEPlot(mus,pt.size = 1.5,colors.use = colors(8,"Dark2"))

# Find markers for every cluster compared to all remaining cells
cluster.markers <- FindAllMarkers(mus, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5)
cluster.markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20
DoHeatmap(mus, genes.use = top20$gene, slim.col.label = TRUE, remove.key = FALSE,cex.row=8)

save.image("E11_15_sc_seurat_analysis.RData")

# Lable cells with embryo stages
cluster <- factor(cell_anno$Stage,levels=c("E11.5","E12.5","E13.5","E14.5","E14.5_GFPn","E15.5"))
names(cluster)<-rownames(cell_anno)
mus@ident <- cluster
TSNEPlot(mus,pt.size = 1.5,colors.use =  c("#FFFF33","#A65628","#F781BF","#191970","#90EE90","#E7298A"))
