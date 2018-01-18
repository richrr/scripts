source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
library(gtools)


getwd()
setwd("D:/MORGUN LAB/11_TEMP_ANALYSIS/R_wd")
directory <- "D:/MORGUN LAB/11_TEMP_ANALYSIS/R_wd/CC_raw"
files.list <-  c(list.files(path=directory))
files.list <- mixedsort(files.list)
singl_condit <- c(rep("P. disiens", 4), rep("P. bivia", 4), rep("L. crispatus", 4), rep("No ctrl", 4))
sample_conditions <- c(rep(singl_condit, 5))
sampleTable <- data.frame(sampleName = files.list, fileName = files.list, condition = sample_conditions)

#I set the formula to condition, but actually I think we need to include experiment too
#to include it into design! the instructions say it is possible later but I did not figure it out yet
ddsHTSeq_1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)

sample_name <- rep(1:20, each=4)
ddsHTSeq_1$sample <- factor(paste("CC",rep(1:20, each=4)))
ddsHTSeq_1$experiment <- factor (rep(1:5, each=16))
ddsHTSeq_1$file <- rep (1:80)
ddsHTS1_coll <- collapseReplicates(ddsHTSeq_1, ddsHTSeq_1$sample, ddsHTSeq_1$file)

#releveling for analysis - pay attantion that I relevel not collapsed data! should be fixed if needed.
#colData(ddsHTSeq_1)$condition <- factor(colData(ddsHTSeq_1)$condition,
 #                                     levels=c("L. crispatus","No ctrl", "P. disiens", "P. bivia"))
#colData(ddsHTSeq_1)$condition <- relevel(colData(ddsHTSeq_1)$condition, "L. crispatus")

#how to access the counts. will use to access normalized counts too
colData(ddsHTSeq_1)
merged_counts <- counts(ddsHTSeq_1)
counts(ddsHTS1_coll)
merged_counts_collapsed <- counts(ddsHTS1_coll)

#example of export of counts
#write.csv(merged_counts, "merged_counts.csv") 

colData(ddsHTS1_coll)
colnames(ddsHTS1_coll)


#here is collapsed datased is made to have experimental design as experiment+condition+experiment:condition
ddsHTSeq_2_coll <- DESeqDataSet(ddsHTS1_coll, design = ~condition + experiment)
colData(ddsHTSeq_2_coll)

#next I need to normalize everything and get normalized data.
rlog_norm_blind <- rlogTransformation(ddsHTSeq_2_coll, blind=TRUE)
vsd_norm_blind <- varianceStabilizingTransformation(ddsHTSeq_2_coll, blind=TRUE)
rlog_norm_design <- rlogTransformation(ddsHTSeq_2_coll, blind=F)
vsd_norm_design <- varianceStabilizingTransformation(ddsHTSeq_2_coll, blind=F)

#writing in the files normalized data
#write.csv(assay(rlog_norm_blind), "test.csv") 
write.csv(assay(rlog_norm_blind), "rlog_norm_blind.csv") 
write.csv(assay(vsd_norm_blind), "vsd_norm_blind.csv") 
write.csv(assay(rlog_norm_design), "rlog_norm_design.csv") 
write.csv(assay(vsd_norm_design), "vsd_norm_design.csv") 

#analyzing in DESeq
colData(ddsHTSeq_2_coll)$condition <- factor(colData(ddsHTSeq_2_coll)$condition, levels=c("L. crispatus","No ctrl", "P. disiens", "P. bivia"))
colData(ddsHTSeq_2_coll)$condition <- relevel(colData(ddsHTSeq_2_coll)$condition, "No ctrl")


ddsHTSeq_2_coll_prosecced <- DESeq(ddsHTSeq_2_coll)
resulsts1 <- results(ddsHTSeq_2_coll_prosecced)
#exporting normalized by DESeq_function counts
write.csv(counts(ddsHTSeq_2_coll_prosecced, normalized = T), "DESeq2_normalization.csv")

#colData(ddsHTSeq_1)$condition <- factor(colData(ddsHTSeq_1)$condition,
#                                     levels=c("L. crispatus","No ctrl", "P. disiens", "P. bivia"))
#colData(ddsHTSeq_1)$condition <- relevel(colData(ddsHTSeq_1)$condition, "L. crispatus")


#I try to do DESeq on specified groups: bact1/2 vs no cont vs l crispatus. so far not successful

res_LvsBivia <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "L. crispatus", "P. bivia"))
write.csv(res_LvsBivia, "res_LvsBivia.csv")

res_LvsDisiens <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "L. crispatus", "P. disiens"))
write.csv(res_LvsDisiens, "res_LvsDisiens.csv")

res_NvsBivia <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "No ctrl", "P. bivia"))
write.csv(res_NvsBivia, "res_NvsBivia.csv")

res_NvsDisiens <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "No ctrl", "P. disiens"))
write.csv(res_NvsDisiens, "res_NvsDisiens.csv")

res_NvsL <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "No ctrl", "L. crispatus"))
write.csv(res_NvsL, "res_NvsL.csv")

res_LvsN <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "L. crispatus", "No ctrl"))
write.csv(res_LvsN, "res_LvsN.csv")

?results
#redoing the tables by making T/C ratio

res_BvsL <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "P. bivia", "L. crispatus"))
write.csv(res_BvsL, "res_BvsL.csv")

res_DvsL <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "P. disiens", "L. crispatus"))
write.csv(res_DvsL, "res_DvsL.csv")

res_BvsN <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "P. bivia", "No ctrl"))
write.csv(res_BvsN, "res_BvsN.csv")

res_DvsN <- results(ddsHTSeq_2_coll_prosecced, contrast=c("condition", "P. disiens", "No ctrl"))
write.csv(res_DvsN, "res_DvsN.csv")
