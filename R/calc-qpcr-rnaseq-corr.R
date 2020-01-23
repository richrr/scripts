library(hash)

args = commandArgs(trailingOnly=TRUE)


source("/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/comp-correl-delta-paired-util-functs.R")


# cd /nfs3/PHARM/Morgun_Lab/richrr/Don_Jump/analysis/djump-oct-2018/relativised/transnet
# Rscript ~/Morgun_Lab/richrr/scripts/R/calc-qpcr-rnaseq-corr.R ~/Morgun_Lab/richrr/Don_Jump/Data/depner-qpcr-ensembl-gene-pair.txt ../depner_rnaseq_rel.mill._quantile.txt.pheno.csv

pairs = read.delim(args[1], header=F,sep='\t')
#head(pairs)

df = read.csv(args[2], header = TRUE,check.names=FALSE,na.strings=c("","na","NA", "Na", "NaN"), row.names=1)
#head(df)


dict = hash()
dict['categ'] = colnames(df)[-1]
#print(dict)

res = calculateCorrelation(pairs, df, 'categ', dict, 'spearman', 1)


CalcMedExpress = function(pair, df){
	c1 = numericizeVector(as.vector(df[pair[1], ]))
    c2 = numericizeVector(as.vector(df[pair[2], ]))
    med = cbind(median(c1, na.rm=T), median(c2, na.rm=T))
	med
}


out = apply(pairs, 1, CalcMedExpress, df)
out = t(out)
colnames(out) = c("Median expression QPCR", "Median expression RNA-Seq")
out

res = cbind(res, out)
head(res)

write.csv(res,"qpcr-rnaseq-corr-output.csv",row.names=FALSE)


#cdf = res[, c(5,2)] # qpcr
cdf = res[, c(6,2)] # rnaseq
cdf = apply(cdf, 2, function(x){as.numeric(as.vector(x))})
cdf <- cdf[order(cdf[,1]),] 
#head(cdf)

pdf("qpcr-rnaseq-corr-plot-output-log-rnaseq.pdf")
plot(cdf,log = "x")
dev.off()
