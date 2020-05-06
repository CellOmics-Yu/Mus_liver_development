library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)

# Load the data
mat <- read.table("all.tube.sc.rpkm.txt",header = T,row.names = 1)
cell_anno <- read.table("mus_all_sc_anno.txt",header = T,row.names = 1,sep = "\t")

# Get Foxa2+ cells from E9.5-E10.5 embryo
cell_anno <- cell_anno[(cell_anno$Stage=="E9.5"|cell_anno$Stage=="E10.0"|cell_anno$Stage=="E10.5")&cell_anno$num_genes_expressed>6000&cell_anno$Foxa2_EGFP=="Express",]
ref <- read.table("Mus_ref.txt",header=T,row.names=1,sep = "\t")
pos <- rownames(mat) %in% rownames(ref)
data <- mat[pos,rownames(cell_anno)]

# Create Seurat object
mus <- CreateSeuratObject(raw.data=data,min.cells = 5, min.genes=1000, is.expr =1,project="Mus_cluster")
mus <- NormalizeData(object = mus,normalization.method = "LogNormalize",scale.factor = 10000)
ref <- ref[!grepl("ERCC-",rownames(ref)), ]
ref <- ref[grepl("protein_coding",ref$Biotpye), ]
ref <- ref[!grepl("^mt-|^Hist",ref$GeneName), ]
ref <- ref[!duplicated(ref$GeneName),]
pos <- rownames(mus@data) %in% rownames(ref)
mus@data <- mus@data[pos,]
rownames(mus@data) <-  ref[rownames(mus@data),5]

# Identification of highly variable genes
mus <- FindVariableGenes(object = mus, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1)
length(mus@var.genes)

# Scaling the data
mus <- ScaleData(object = mus)

# Perform linear dimensional reduction
mus <- RunPCA(object = mus, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
mus <- ProjectPCA(object = mus, do.print = FALSE)
PCAPlot(object = mus, dim.1 = 1, dim.2 = 2,pt.size = 3)
PCHeatmap(mus, pc.use = 1:20, cells.use = 100, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCElbowPlot(mus)

# Cluster the cells
mus <- FindClusters(mus, reduction.type = "pca",dims.use = 1:5, save.SNN = T,resolution =0.8,temp.file.location="/hwfssz1/ST_PRECISION/USER/zhongyu/1.Project/Mus_Liver/seurat/results/",force.recalc=T,k.param=15,k.scale=5)
PrintFindClustersParams(object = mus)

# Run non-linear dimensional reduction (tSNE)
mus <- RunTSNE(mus, dims.use = 1:5, do.fast = T,perplexity=15,theta=0)
TSNEPlot(mus,pt.size = 1.5,colors.use = colors(8,"Dark2"))

# Find markers for every cluster compared to all remaining cells
cluster.markers <- FindAllMarkers(mus, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5)
cluster.markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20
DoHeatmap(mus, genes.use = top20$gene, slim.col.label = TRUE, remove.key = FALSE,cex.row=8)

save.image("E9_10_foxa2p_seurat_analysis.RData")

# Lable cells with embryo stages
cluster <- factor(cell_anno$Stage,levels=c("E9.5","E10.0","E10.5"))
names(cluster)<-rownames(cell_anno)
mus@ident <- cluster
TSNEPlot(mus,pt.size = 1.5,colors.use =  c("#4DAF4A","#984EA3","#FF7F00"))
