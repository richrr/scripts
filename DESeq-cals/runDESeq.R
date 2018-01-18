###
# Skeleton of this code was picked from Dr. Stephen Turner
# https://gist.github.com/stephenturner/f60c1934405c127f09a6
# this link has good explanation: https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
# this link has suggestions about how to do multifactorial deisgns: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Sep 28 2017
###

# cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/DESEQ
# [Linux@transkingdom DESEQ]$ Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/DESeq-cals/runDESeq.R ../summarize_per_sample/summed-values-abx2.csv ncd expt1_hfhs_ncd_
# [Linux@transkingdom DESEQ]$ Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/DESeq-cals/runDESeq.R ../summarize_per_sample/summed-values-abx3.csv ncd expt2_hfhs_ncd_



### things to consider:
# (i) do we normalize per compairson or normalize all samples once then do the comparisons? per comparison
# (ii) the way to handle outliers. currently for less than 7 replicates, even a single outlier sample completely removes the gene form diff analysis
#		This filtering can be turned off with results(dds, cooksCutoff=FALSE). May be ok since we do meta-analysis? May increase false positives.
#
# (iii) where is the independentfiltering option? Apparently done in the results method and by default it is true. 
# (iv) why is the heatmap not plotting pvalues


args = commandArgs(trailingOnly=TRUE)

# non-normalized rna-seq raw counts obtained from summing different lanes for the sample.


expressionDataFile = args[1]

expressionData = read.csv(expressionDataFile,header = TRUE,check.names=FALSE,na.strings=c("","na","NA", "Na", "NaN"), row.names=1)

controlstr = args[2] # e.g. ncd (case sensitive)
outstr = args[3] # e.g. hfhs_ncd

# for testing purposes only used the ncd and hfhs
countdata <- expressionData[ ,1:10]

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first 5 are controls, second 5 contain the expansion)
(condition <- factor(c(rep("ncd", 5), rep("hfhs", 5))))
print("Input")
condition


#### something like this HAS to be used, else the the lexicographic sort gives ncd/hfhs
#### that may also affect the following results, so be careful that this relevel works
#### the way you want it to and samples are not incorrectly assigned.
#####samples$condition <- relevel(samples$condition, controlstr)
condition = relevel(condition, controlstr)
print("Releveled. Note the order of the Levels: ")
condition

library(DESeq2)

##################
# there has to be a check to make sure the samples are ordered in the exact
# way in the sample and mapping files.
if(FALSE){
cat("In the three examples above, we applied the function factor to the column of interest in colData,
supplying a character vector of levels. It is important to supply levels (otherwise the levels are chosen
in alphabetical order) and to put the control or untreatedlevel as the first element, so that the log2
fold changes and results will be most easily interpretable. A helpful R function for easily changing the
base level is relevel. An example of setting the base level with relevel is:
colData(dds)$condition <- relevel(colData(dds)$condition, control)
The reason for the importance of the specifying the base level is that the function model.matrix
is used by the DESeq2 package to build model matrices, and these matrices will be used to compare
all other levels to the base level.")
}
###################


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
#coldata$condition <- relevel(coldata$condition, controlstr)

coldata$condition



dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds




# Run the DESeq pipeline
dds <- DESeq(dds)


# Plot dispersions
png(paste0(outstr, "qc-dispersions.png"), 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds) #blind=TRUE by default.  blind=TRUE should be used for comparing samples in an manner unbiased by prior information on
          #samples, for example to perform sample QA (quality assurance).

head(assay(rld)) # assay seems to be a function in deseq2 that allows tabular display of the summarized experiment object which contains the regularized log-transformed values in its assay slot.

pdf(paste0(outstr, "assay-hist.pdf"))
hist(assay(rld))
dev.off()




# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png(paste0(outstr, "qc-heatmap-samples.png"), w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()


# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")  ## shouldn't this have PC2?
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png(paste0(outstr, "qc-pca.png"), 1000, 1000, pointsize=20)
#rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))   
plotPCA(rld, intgroup="condition")
dev.off()


# Get differential expression results
res <- results(dds) # the default is independentFiltering = TRUE . 
					# The pvalue calc or number of rows returned does not change if turned OFF (independentFiltering =FALSE), 
					# only, the adj p value improves when the filtering is ON (independentFiltering =TRUE).
					# as per (http://seqanswers.com/forums/showthread.php?t=62871), if you want consistent filtering across different comparisons
					# you could turn independentFiltering=FALSE 
					#
					# Considering the fact that we perform meta-analysis and then calc our own FDR without using the FDR directly from DESEQ, it doesn't matter
					# if we keep independentFiltering ON or OFF. Infact keeping the default may allow us to directly look at the DESEq outputs "as is".
					#(https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/results)
					# The adjusted p-values for the genes which do not pass the filter threshold are set to NA
					# By default, results assigns a p-value of NA to genes containing count outliers, as identified using Cook's distance.
					
					# finally, although the genes can be filtered just before fdr calc (using filter=rv) (row variance), https://support.bioconductor.org/p/61620/
					# it does not help our case
					# since we do not care about the fdr, plus it would be hard to compare the expts. and other comparisons if some of the genes are removed in results

					
#res = results(dds, independentFiltering = FALSE)

# print out the genes less than 5% fdr
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]


## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file=paste0(outstr, "diffexpr-results.csv"))

pdf(paste0(outstr, "hist-pval.pdf"))
## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")
dev.off()


## Examine independent filtering
#attr(res, "filterThreshold")  # this is for the old version
#metadata(res)
#metadata(res)$filterThreshold  #https://www.biostars.org/p/174695/
pdf(paste0(outstr, "independentFilterPlot.pdf"))
plot(metadata(res)$filterNumRej, type="b", xlab="quantiles of baseMean", ylab="number of rejections")
dev.off()




## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png(paste0(outstr, "diffexpr-maplot.png"), 1500, 1000, pointsize=20)
#maplot(resdata, main="MA Plot")
plotMA(dds, ylim=c(-1,1), cex=1)
dev.off()


## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  require(calibrate)
  require(RColorBrewer)
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png(paste0(outstr, "diffexpr-volcanoplot.png"), 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()



print("Done")




