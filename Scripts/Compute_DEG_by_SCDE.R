## Usage: Rscript <cell_anno> <group1> <group2> <out_prefix>
library(methods)
library(scde)

args = commandArgs(T)
prefix = args[4]

# load dataset
count <- read.table("reads_count_of_sc.txt", header = T,row.names = 1,stringsAsFactors = F)
anno <- read.table(args[1], header = T,row.names = 1,stringsAsFactors = F,sep="\t")
anno <- anno[anno$cluster==args[2]|anno$cluster==args[3],]
group1 <- anno[anno$cluster==args[2],]
group2 <- anno[anno$cluster==args[3],]
count <- count[,rownames(anno)]

# factor determining cell types
sg <- factor(anno$cluster,levels =c(args[2],args[3]))
# the group factor should be named accordingly
names(sg) <- rownames(anno)

ref <- read.table("Mus_ref.txt",header = T,row.names = 1,stringsAsFactors = F)
ref <- ref[,c(5,7)]

# clean up the dataset
cd <- clean.counts(count, min.lib.size=1000, min.reads = 1, min.detected = 2)

# calculate models
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 6, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)

# define two groups of cells
groups <- factor(anno[rownames(o.ifm),11],levels =c(args[2], args[3]))
names(groups) <- row.names(o.ifm)

# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  6, verbose  =  1)

mat1 <- count[rownames(ediff),rownames(group1)]
mat2 <- count[rownames(ediff),rownames(group2)]
fc <- log2(rowMeans(mat1+1)/rowMeans(mat2+1))

# write out a table with all the results, showing most significantly different genes (in both directions) on the top
ediff <- cbind(ediff,ref[rownames(ediff),],group1=rowMeans(mat1),group2=rowMeans(mat2),log2FC=fc)
colnames(ediff)[9:10]=c(args[3],args[4])
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = paste(prefix, "DE_results.txt",sep="_"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
