library(Rcmdr)

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)



getwd()
setwd("D:/MORGUN LAB/11_TEMP_ANALYSIS/R_wd")

dds <- makeExampleDESeqDataSet(m=12)
colData(dds)
colnames(dds)
# make data with two technical replicates for three samples
dds$sample <- factor(sample(paste0("sample",rep(1:9, c(2,1,1,2,1,1,2,1,1)))))
dds$run <- paste0("run",1:12)
colData(dds)
colnames(dds)
ddsColl <- collapseReplicates(dds, dds$sample, dds$run)
# examine the colData and column names of the collapsed data
colData(ddsColl)
colnames(ddsColl)

#?do not understand the last line:
# check that the sum of the counts for "sample1" is the same
# as the counts in the "sample1" column in ddsColl
matchFirstLevel <- dds$sample == levels(dds$sample)
stopifnot(all(rowSums(counts(dds[,matchFirstLevel])) == counts(ddsColl[,1])))

vignette("DESeq2")
