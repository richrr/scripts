#####################################################################################################################
# script written by Lina Thomas
# date: 12/11/2013

# run this script when output.option="one_table_each" in analysis_AE.R

######################################################################################################################

# command in cgrb
# 
# SGE_Batch -c 'R < scripts/analysis_AE/table_cor_between_differences.R "cgrb" "/capecchi/pharmacy/morgunlab/lina/ArrayExpress" "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/DAPs_Dif_first100" "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/DEGs_Dif" "DAPs_Dif_table_first100.txt" "DEGs_Dif_table.txt" "correlacao_diferencas_first100daps.txt"   --no-save' -r cor_differences_first100daps -f 4G -m 4G 
##############################################
#libraries

require(utils)
library(plyr)

#functions
########################################################################################################################################
source("/raid1/home/pharmacy/thomasli/scripts/analysis_AE/functions_daps_analysis.R")


# Define a function to calculate sample sizes for gene-gene correlations -- courtesy of Vivek Gopalan

# x           is the sequence from 1 to G (i.e. 1:nrow(aa)), where G = number of genes
# aa          is a G x S matrix of gene expression values, where G = number of genes and S = number of samples

# WARNING!!!  The input aa must be a matrix -- 'aa' cannot be a data frame.

min_size_cal         <- function(x,aa){rowSums(!is.na(t(t(aa)*aa[x,])))}

####################################################################################################

# Define a function to calculate sample sizes for gene-gene correlations from 2 different matrices

# x           is each row of B, when calling, write 1:row(B)
# A          is a G_A x S matrix of gene expression values, where G = number of genes and S = number of samples
# B          is a G_B x S matrix of gene expression values, where G = number of genes and S = number of samples

# WARNING!!!  A and B must have the same number of columns!!!!

# WARNING!!!  The inputs A and B must be matrices -- A and B cannot be a data frames.

min_size_cal_2m <- function(x,A,B){rowSums(!is.na(t(t(A)*B[x,])))}

###########################################################################################
args=commandArgs(trailingOnly = FALSE)

local=args[2]

local
path_final=args[3]
path_daps=args[4]
path_degs=args[5]
daps_difcor_table.file=args[6]
degs_difmed_table.file=args[7]
table.out = args[8]

print("did the arguments")
#################################

#correlation between delta cor in DAPS and delta median in DEGs

cor_met="pearson"

setwd(path_daps)
print(path_daps)
print(daps_difcor_table.file)
daps_orig=read.table(daps_difcor_table.file, header = TRUE, sep = "\t")

setwd(path_degs)
degs_orig=read.table(degs_difmed_table.file, header = TRUE)

daps.genesymbol = daps_orig[, 1:3]
daps.genesymbol[,"pair_probe.id"] = paste(daps.genesymbol[,"ID.x"], daps.genesymbol[,"ID.y"])

degs.genesymbol = degs_orig[, 1:2]

daps = t(daps_orig[,c(6:ncol(daps_orig))])
colnames(daps)=paste(daps_orig[,"ID.x"], daps_orig[,"ID.y"])

degs = t(degs_orig[,-c(1:2)])
colnames(degs)=degs_orig[,1]

# se o numero de samples for diferente nas bases de degs de daps.
cond1 = which(!(row.names(degs)%in% row.names(daps)))
if(length(cond1)>0)
  degs=degs[-cond1,]

cond2 = which(!(row.names(daps)%in% row.names(degs)))
if(length(cond2)>0)
  daps=daps[-cond2,]


#complete pairwise correlation between each daps and degs columns
cor=cor(daps,degs,use="pairwise.complete.obs",method=cor_met)

#sample size used to calculate complete pairwise correlation
min.mat=sapply(1:nrow(t(degs)),FUN=min_size_cal_2m,as.matrix(t(daps)),as.matrix(t(degs)))
colnames(min.mat)=colnames(degs)

# t statistics
t <- (sqrt(min.mat-2)*cor)/sqrt(1-cor^2)

#pvalue
pvalue <- 2*pt(-abs(t), min.mat-2)


#####################################################################################################
#stacking

cor=stack(as.data.frame(cor))

pvalue=stack(as.data.frame(pvalue))

fcor=data.frame(rep(colnames(daps),ncol(degs)),cor[,2],cor[,1],pvalue[,1])

fcor[,5] = p.adjust(as.numeric(fcor[,4]),method="fdr")

colnames(fcor)=c("DAPs","DEGs","cor","pvalue","fdr")

#############################################################################3
# getting gene symbol

fcor = merge(fcor, degs.genesymbol, by.x = "DEGs", by.y = 1, all.x = T)

fcor = merge(fcor, daps.genesymbol[,c("pair_probe.id", "pair")], by.x = "DAPs", by.y = "pair_probe.id", all.x = T)

fcor = fcor[ !is.na(fcor$fdr) , ]

fcor = fcor[ sort(fcor$fdr) , ]

setwd(path_final)
write.table(fcor[,c(6,7,1:5)],file=table.out, col.names = c("DAPs.Gene.symbol", "DEGs.Gene.symbol","DAPs.probe.id","DEGs.probe.id", "correlation between differences", "pvalue","fdr"))

