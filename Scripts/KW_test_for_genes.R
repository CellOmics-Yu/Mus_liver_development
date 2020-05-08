# Kruskal-Wallis test for selecting variable genes across cell type
# usage:  Rscript KW_test.R <FI: matrix> <FI: cluster info> <FO: output>
args <- commandArgs(T)

MAT <- read.table(args[1], header=TRUE)
MAT=log2(MAT + 1)
sample <- read.table(as.character(args[2]), header=TRUE, row.names=1)
MAT = MAT[, rownames(sample)]

# Filter genes exressed in less than 1% cells
num = trunc(ncol(MAT) * 0.01)
filt = apply(MAT, 1, function(x) length(x[x>1]) > num)
MAT = MAT[filt,]

MAT <- data.frame(MAT, p.value=NA, q.value=NA, KW.value=NA)
sample <- sample[match(colnames(MAT)[1:(ncol(MAT)-3)], rownames(sample)),]

# Kruskal-Wallis test for selecting variable genes across cell type
for (marker in c(1:nrow(MAT))){
  DAT <- data.frame(CN=as.vector(t(MAT[marker,1:(ncol(MAT)-3)])),Type=sample)
  res <- kruskal.test(CN ~ Type, data = DAT)
  MAT$KW.value[marker] <- res[1]$statistic
  MAT$p.value[marker] <- res[3]$p.value
}

MAT$q.value = p.adjust(MAT$p.value, method="fdr", n=length(MAT$p.value))
MATo <- data.frame(gene=rownames(MAT), p.value=MAT$p.value, q.value=MAT$q.value, KW.value=MAT$KW.value)

Gene = read.table("../Data/Gene_info.txt", header=T, row.names=1)
MATo$name=Gene[rownames(MAT),1]

write.table(MATo, file=as.character(args[3]),  quote=FALSE, sep="\t", row.names=FALSE)
