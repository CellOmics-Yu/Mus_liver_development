## load class definition and functions
source("./RaceID_class.R")

mat <- read.table("../Data/E9_11_GFPp.sc.rpkm.KW.e2.txt", header=T)
prdata <- log2(mat+1)

RaceID = function (prdata, downsample=F, cln=0, main="All_cells", test=FALSE){
  pdf(paste(main,"pdf",sep="."))
  sc <- SCseq(prdata)
  nums=2
  #14: quantile: 0.9999
  sc <- filterdata(sc, mintotal=3000, minexpr=1, minnumber=nums, maxexpr=14,
                   downsample=as.logical(downsample), dsn=1, rseed=17000)
  sc <- clustexp(sc, clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,
                 SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=as.numeric(cln),
                 rseed=17000)
  plotgap(sc)
  plotsilhouette(sc)
  plotjaccard(sc)
  abline(h=0.6)
  
  if(!test){
    clustheatmap(sc,final=FALSE,hmethod="single")
    sc <- comptsne(sc,rseed=15555)
    plottsne(sc,final=FALSE)
    stages = unlist(strsplit(names(sc@ndata), "_"))[seq(2,length(names(sc@ndata))*3, by=3)]
    plotsymbolstsne(sc,types=stages)
    mat_tmp = prdata[1:2, names(sc@ndata)]
    filt = apply(mat_tmp, 2, function(x) all(x>1))
    Foxa2_EGFP=c()
    Foxa2_EGFP[filt] = "Express"
    Foxa2_EGFP[!filt] = "No_express"
    plotsymbolstsne(sc,types=Foxa2_EGFP)
    plotlabelstsne(sc, labels=names(sc@ndata))
    
    x <- data.frame(CELLID=names(sc@kmeans$kpart),cluster=sc@kmeans$kpart)
    write.table(x[order(x$cluster,decreasing=FALSE),],paste(main,"cell_clust.txt", sep="_"),
                row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  }
  dev.off()
  sc
}

plot_genes = function(sc, genes){
  for(i in 1:nrow(genes)){
    g = genes[i,2]
    if(!rownames(sc@expdata[g,])=="NA"){
      plotexptsne(sc,g,n=genes[i,1],logsc=F)
    }else{
      message(paste(genes[i,1],"is filtered!",sep=" "))
    }
  }
}

prefix="E9_11_GFPp"
sc5=RaceID(prdata, main=paste(prefix, "5",sep="_"), cln=5)
num = 5

# computes differentially expressed genes in a cluster by comparing to all remaining cells outside of the cluster 
# based on a negative binomial model of gene expression
sample <- read.table("Sampe_cell_cluster.txt", header=TRUE, row.names=1)
sample <- sample[names(sc5@kmeans$kpart),]
sc5@kmeans$kpart <- sample
cdiff <- clustdiffgenes(sc5, pvalue=.05)

genes = read.table("../Data/Gene_info.txt", header=T, row.names=1, stringsAsFactors=F)
tfs=read.table("../Data/Mus_musculus_transcription_factors_gene_list.txt", header=T, row.names=1, stringsAsFactors=F)

info=cbind(cdiff$cl.1, as.character(genes[rownames(cdiff$cl.1), 1]), "cluster1")
pos = which(rownames(cdiff$cl.1) %in% rownames(tfs))
info[pos,8] = "tfs"
info[-pos,8] = "-"
colnames(info)[6:8] = c("Symbol", "Cluster", "TFs")
write.table(info, file=paste(prefix, num, "genes", sep="_"), quote=F, sep="\t")
for(i in 2:num){
  if(nrow(cdiff[[i]]) > 0){
    tmp = cbind(cdiff[[i]], as.character(genes[rownames(cdiff[[i]]), 1]), paste("cluster", i, sep=""))
    pos = which(rownames(cdiff[[i]]) %in% rownames(tfs))
    tmp[pos,8] = "tfs"
    tmp[-pos,8] = "-"
    write.table(tmp, file=paste(prefix, num, "genes", sep="_"), quote=F, sep="\t", col.names=F, append=T)
    colnames(tmp)[6:8] = c("Symbol", "Cluster", "TFs")
    info = rbind(info, tmp)
  }
}

filt=grep("Gm|Rik|Mir|mir|RP|Rp|Rn|mt", info[,6])
write.table(info[-filt,], file=paste(prefix, num, "genes_filt", sep="_"), quote=F, sep="\t")
markers = cbind(as.character(info[-filt,6]), as.character(rownames(info[-filt,])))
pdf(paste(prefix,num, "markers.pdf", sep="_"))
plot_genes(sc5, markers)
dev.off()

genes = read.table("../9.all_tube/marker.info.reorder.txt", header=T, stringsAsFactors=F)
pdf(paste(prefix, num, "add.pdf",sep="_"))
plot_genes(sc5, genes)
dev.off()
