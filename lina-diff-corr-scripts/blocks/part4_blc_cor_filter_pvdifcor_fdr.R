######################################################################################################

# Lina Thomas
# date: 1/22/14

########################################################################################################################################
# 
####################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < scripts/blocks/part4_blc_cor_filter_pvdifcor_fdr.R --no-save "cgrb" 3 "Project_03mp" NULL 0.2 0.2 "Identifiers.csv" "gene_identificators_probeid_genesymbol.csv" "dados1" "dados2" "dados3" "table_same_dir_difcor.txt"'  -r part4_imsys -f 256G -m 256G -q samwise
# only two datasets 10/16/2015
# SGE_Batch -c 'R < scripts/blocks/part4_blc_cor_filter_pvdifcor_fdr.R --no-save "cgrb" 2 "Project_03mp_2datasets" NULL 0.2 0.1 "Identifiers.csv" "gene_identificators_probeid_genesymbol.csv" "dados1" "dados3" "table_same_dir_difcor.txt"'  -r part4_imsys_2datasets -f 256G -m 256G -q transkingdom
# 
# tests
# SGE_Batch -c 'R < scripts/blocks/part4_blc_cor_filter_pvdifcor_fdr.R --no-save "cgrb" 3 "Project_test" NULL 0.2 0.2 "Identifiers.csv" "gene_identificators_probeid_genesymbol.csv" "amostra1" "amostra2" "amostra3" "table_same_dir_difcor.txt"'  -r part4_imsys -f 1G -m 1G -q samwise
# at lab: R < Dropbox/Tese/scripts/DAPfinder/blocks/part4_blc_cor_filter_pvdifcor_fdr.R --no-save "lab" 3 "Project_teste" NULL 0.2 0.2 "Identifiers.csv" "gene_identificators_probeid_genesymbol.csv" "amostra1" "amostra2" "amostra3" "table_same_dir_difcor.txt" 
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

local = args[3]

print(local)

total_studies = as.numeric(args[4]) #number of datasets
newfoldername = args[5]
method = args[6]
pv_thrsd_cor = args[7]
pv_thrsd_difcor = args[8]
class.file = args[9]
gene.symbol.file = args[10]
# by.study = args[8]
 
study.names = vector()
for(counter1 in 1:total_studies)
  study.names[counter1] = args[counter1 + 10] #name of study file
counter1 = 10 + counter1

table.file = args[counter1 + 1]

##############################################################################################
# 
# if(!"psych" %in% installed.packages()) install.packages("psych")
# library("psych")

################################################################################

if(local=="home")
  path_dados="C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO" else
    if(local=="lab")
      path_dados="/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO" else
        if(local=="cgrb")
          path_dados="/capecchi/pharmacy/morgunlab/lina/imsys" 

if(local == "cgrb")
  path_final=paste(path_dados, newfoldername ,sep="/") else
    if(local == "lab")
      path_final=paste(path_dados, "DAPfinder", newfoldername ,sep="/")

#################################################################################################################

#functions

if(local=="home")
  source("C:/Users/Lina/Dropbox/Tese/scripts/functions.R") else
    if(local=="lab")
      source("/home/thomasli/Dropbox/Tese/scripts/functions.R") else
        if(local=="cgrb") 
          source("/raid1/home/pharmacy/thomasli/scripts/functions.R")

print("pegou as funcoes")

coherent.cor.directions.one.state=function(final_table,pv_thrsd_cor = 0.2, total_studies=1, study.names = "gr1" ,class_type=c(0,1), method = "pearson"){
  
  #checking if final_table already has cor directions assigned
  dir_assigned = grep("dircor",colnames(final_table))
  
  if(length(dir_assigned) == 0)
    final_table = cor.directions(final_table, pv_thrsd_cor, total_studies, study.names , class_type, method = method)
  
  # only keeping the pairs with the same correlation directions in all stydies for both classes
  
  dircor0_index = grep("dircor0", colnames(final_table))
  sum_dircor0 = apply(final_table[,dircor0_index],1,sum)
  
  dircor1_index = grep("dircor1", colnames(final_table))
  sum_dircor1 = apply(final_table[,dircor1_index],1,sum)
  
  cond = abs(sum_dircor0) == total_studies | abs(sum_dircor1) == total_studies 
  final_table = final_table[cond,]
  
  return(final_table)
}


#################################################################################################
setwd(path_dados)
class.table = read.csv(class.file,header=TRUE)
class_type = class.table[!duplicated(class.table[,"class"]),"class"] 

setwd(path_final)
final_table = read.table(table.file, header=T, stringsAsFactors=F)

# we only want pairs where correlations are significant in at least one state
final_table = coherent.cor.directions.one.state(final_table, pv_thrsd_cor, total_studies, study.names, class_type)

final_table = pv.dif.cor.table(final_table, study.names, class_type)

final_table = coherent.dircor.pvalue(final_table, pv_thrsd_difcor)

final_table[, "fisher_pv"] = pv.meta.analysis.v2( final_table[, grep("pv_difcor", colnames(final_table))] )

final_table[, "fisher_fdr"] = p.adjust(final_table[,"fisher_pv"], method = "fdr")

setwd(path_final)
write.table(final_table, paste("final_table_cor",gsub("[[:punct:]]","",as.character(pv_thrsd_cor)),"_","difcor_", gsub("[[:punct:]]","",as.character(pv_thrsd_difcor)),".txt",sep=""), row.names = FALSE, sep="\t")

#################################################################################################
# Gene.symbol
setwd(path_dados)
gene.symbol = read.table(gene.symbol.file, header = T, sep = "\t")

final_table[, "gene.x"] = substr(final_table[ , "gene.x"], 2, nchar(final_table[ , "gene.x"]))
final_table[, "gene.y"] = substr(final_table[ , "gene.y"], 2, nchar(final_table[ , "gene.y"]))

final_table = merge(final_table, gene.symbol, by.x = "gene.x", by.y = "probe.id", all.x = T)
final_table = merge(final_table, gene.symbol, by.x = "gene.y", by.y = "probe.id", all.x = T)

setwd(path_final)
write.table(final_table[,c(ncol(final_table) - 1, ncol(final_table), 2, 1, seq(3,ncol(final_table)-2))], paste("final_table_cor",gsub("[[:punct:]]","",as.character(pv_thrsd_cor)),"_","difcor_", gsub("[[:punct:]]","",as.character(pv_thrsd_difcor)),"_with_genesymbol.txt",sep=""), row.names = FALSE, sep="\t")
###############################################################################################
# 
# setwd(path_final)
# 
# pv_thrsd_cor = 0.2
# pv_thrsd_difcor = 0.2
# 
# final_table = read.table(paste("final_table_cor",gsub("[[:punct:]]","",as.character(pv_thrsd_cor)),"_","difcor_", gsub("[[:punct:]]","",as.character(pv_thrsd_difcor)),"_with_genesymbol.csv",sep=""), header = T, sep = ",", stringsAsFactors=F)

#########################################################################
pv_thrsd_difcor = 0.1

final_table = coherent.dircor.pvalue(final_table, pv_thrsd_difcor)

final_table[, "fisher_pv"] = pv.meta.analysis.v3( final_table[, grep("pv_difcor", colnames(final_table))] )

final_table[, "fisher_fdr"] = p.adjust(final_table[,"fisher_pv"], method = "fdr")

final_table = final_table[final_table[,"Gene.symbol.x"] != "" & final_table[,"Gene.symbol.y"] != "",]

write.table(final_table, paste("final_table_cor",gsub("[[:punct:]]","",as.character(pv_thrsd_cor)),"_","difcor_", gsub("[[:punct:]]","",as.character(pv_thrsd_difcor)),".txt",sep=""), row.names = FALSE, sep="\t")

#fdr less than 0.01
daps = final_table[final_table[,"fisher_fdr"] < 0.01,]

daps_table = table(c(daps[,"Gene.symbol.x"], daps[,"Gene.symbol.y"]))

write.table(daps_table,"daps_fdr001_freq_table", sep = "\t", row.names = F)

#fdr less than 0.001
daps = final_table[final_table[,"fisher_fdr"] < 0.001,]

daps_table = table(c(daps[,"Gene.symbol.x"], daps[,"Gene.symbol.y"]))

write.table(daps_table,"daps_fdr0001_freq_table", sep = "\t", row.names = F)


# ##########################################################################
# 
# pv_thrsd_cor = 0.1
# 
# final_table = final_table[, c(1:25)]
# 
# final_table = coherent.cor.directions.one.state(final_table, pv_thrsd_cor, total_studies, study.names, class_type)
# 
# final_table = pv.dif.cor.table(final_table, study.names, class_type)
# 
# final_table = coherent.dircor.pvalue(final_table, pv_thrsd_difcor)
# 
# final_table[, "fisher_pv"] = pv.meta.analysis.v2( final_table[, grep("pv_difcor", colnames(final_table))] )
# 
# final_table[, "fisher_fdr"] = p.adjust(final_table[,"fisher_pv"], method = "fdr")
# 
# setwd(path_final)
# write.table(final_table, paste("final_table_cor",gsub("[[:punct:]]","",as.character(pv_thrsd_cor)),"_","difcor_", gsub("[[:punct:]]","",as.character(pv_thrsd_difcor)),".txt",sep=""), row.names = FALSE, sep="\t")
