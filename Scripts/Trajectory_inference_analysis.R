## Build single-cell trajectories from liver primordium to liver maturation

library("monocle")
library("reshape2")

## load data
expr_matrix <- read.table("all.tube.sc.rpkm.txt",header=T,row.names=1)
sample_sheet <- read.table("mus_foxa2p_all_sc.txt",header=T,row.names=1,sep="\t")
sample_sheet <- sample_sheet[sample_sheet$Cell_Cluster=="Liver primordium"|sample_sheet$Cell_Cluster=="Liver bud"|sample_sheet$Cell_Cluster=="Hepatocyte",]
sample_sheet$Stage <-factor(sample_sheet$Stage,levels = c("E9.5", "E10.0", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5"))
sample_sheet$Cell_Cluster <- factor(sample_sheet$Cell_Cluster,levels = c("Liver primordium","Liver bud","Hepatocyte"))

ref <- read.table("Mus_ref.txt",header=T,row.names=1,sep="\t")
ref <- ref[!grepl("ERCC-",rownames(ref)), ]
ref <- ref[grepl("protein_coding",ref$Biotpye), ]
ref <- ref[!grepl("^mt-|^Hist",ref$gene_short_name), ]
ref <- ref[!duplicated(ref$gene_short_name),]
"load is done"

## creat dataset and retain genes expressed in at least 1% cells with RPKM > 1
pos = which(rownames(ref) %in% rownames(expr_matrix))
ref = ref[pos,]
expr_matrix = expr_matrix[rownames(ref),rownames(sample_sheet)]
dim(expr_matrix)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = ref)
MICE <- newCellDataSet(as.matrix(expr_matrix), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower = 0.1))
rpc_matrix <- relative2abs(MICE, method = "num_genes")
MICE <- newCellDataSet(as(as.matrix(rpc_matrix),"sparseMatrix"), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size(), lowerDetectionLimit = 0.5)

MICE <- estimateSizeFactors(MICE)
MICE <- estimateDispersions(MICE)
MICE <- detectGenes(MICE, min_expr = 1)
expressed_genes <- row.names(subset(fData(MICE), num_cells_expressed > 10))
paste("dataset is created at", date(),sep=" ")

## find differentially expressed genes
diff_test_res <- differentialGeneTest(MICE[expressed_genes,],fullModelFormulaStr = "~Stage",cores = 6)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.001))
ordering_genes <- intersect(ordering_genes, expressed_genes)
paste("finding differentially expressed genes is done at", date(),sep=" ")

## order cells by progress
MICE <- setOrderingFilter(MICE, ordering_genes)
MICE <- reduceDimension(MICE,max_components = 2) ## Reducing to independent components so Monocle will work better
#MICE <- orderCells(MICE,reverse=TRUE) ##  num_paths allows Monocle to assign cells to one of several alternative fates
MICE <- orderCells(MICE)
plot_cell_trajectory(MICE, show_cell_names = F, color_by = "Stage")
plot_cell_trajectory(MICE, show_cell_names = F, color_by = "Cell_Cluster")
plot_cell_trajectory(MICE, show_cell_names = F, color_by = "State")
plot_cell_trajectory(MICE, show_cell_names = F, color_by = "Pseudotime")

dev.off()
paste("ordering is done at", date(),sep=" ")

## Finding genes that change as a function of pseudotime
diff_test_res <- differentialGeneTest(MICE[expressed_genes,],fullModelFormulaStr ="~sm.ns(Pseudotime)",cores = 6)
paste("finding differentially expressed genes is done at", date(),sep=" ")

## Clustering genes by pseudotemporal expression pattern
diff_test_res <- diff_test_res[order(diff_test_res$qval),]
filt <- grepl(("^Gm|^Mir|Rik$"), diff_test_res$gene_short_name) # Filtering genes with unknown function
diff_test_res <- diff_test_res[!filt,]
#sig_gene_names <- row.names(diff_test_res[1:200,]) # Select top 200 gene used to cluster
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
plot_pseudotime_heatmap(MICE[sig_gene_names,],num_clusters = 2,cores = 6,show_rownames = F) # you can set different cluster number with num_clusters
paste("Clustering genes by pseudotemporal is done at", date(),sep=" ")
save.image("E9_15_liver_monocle.Rdata")
