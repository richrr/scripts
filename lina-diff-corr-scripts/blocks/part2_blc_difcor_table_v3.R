######################################################################################################
#
# Lina Thomas
# date: 1/22/14
#
########################################################################################################################################
# 
####################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < scripts/blocks/part2_blc_difcor_table_v3.R --no-save FALSE NULL "cgrb" 3 TRUE "Project_03mp" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE "dados1.txt" "dados2.txt" "dados3.txt" "Identifiers.csv"' -r dapfinder_part2 -f 256G -m 256G -q samwise
# SGE_Batch -c 'R < scripts/blocks/part2_blc_difcor_table_v3.R --no-save FALSE NULL "cgrb" 2 TRUE "Project_mp03_2datasets" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE " "Gene_Expression_balbc_num.txt" "Gene_Expression_B10A_lit_num.txt" "Identifiers.csv" ' -r part2_imsys_2datasets -f 100G -m 100G -q transkingdomSGE_Batch -c 'R < scripts/blocks/part2_blc_difcor_table.R --no-save FALSE NULL "cgrb" 2 TRUE "Project_mp03_2datasets" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE " "Gene_Expression_balbc_num.txt" "Gene_Expression_B10A_lit_num.txt" "Identifiers.csv" ' -r part2_imsys_2datasets -f 100G -m 100G -q transkingdom
# in the lab
# SGE_Batch -c 'R < /media/thomasli/scripts/DAPfinder/blocks/part2_blc_difcor_table_v3.R --no-save FALSE NULL "cgrb" 2 TRUE "Project_mp03_2datasets" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE " "Gene_Expression_balbc_num.txt" "Gene_Expression_B10A_lit_num.txt" "Identifiers.csv" ' -r part2_imsys_2datasets -f 100G -m 100G -q transkingdomSGE_Batch -c 'R < scripts/blocks/part2_blc_difcor_table.R --no-save FALSE NULL "cgrb" 2 TRUE "Project_mp03_2datasets" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE " "Gene_Expression_balbc_num.txt" "Gene_Expression_B10A_lit_num.txt" "Identifiers.csv" ' -r part2_imsys_2datasets -f 100G -m 100G -q transkingdom

# tests
# SGE_Batch -c 'R < scripts/blocks/part2_blc_difcor_table_v3.R --no-save FALSE NULL "cgrb" 1 TRUE "Project_test" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE "amostra1.txt" "Identifiers.csv"' -r cor_difcor_imsys_test -f 4G -m 4G -q samwise
# at lab: R < Dropbox/Tese/scripts/DAPfinder/blocks/part2_blc_difcor_table_v3.R --no-save FALSE NULL "lab" 1 TRUE "Project_teste" "filter_teste.txt" TRUE "pearson" "pairwise.complete.obs" TRUE "amostra1.txt" "Identifiers.csv"
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

filtro_BRB = args[3]
mp = as.numeric(args[4])
local = args[5]
total_studies = as.numeric(args[6]) #number of datasets
class_same_file = args[7]
newfoldername = args[8]
filterfile = args[9]
intersect.filter = args[10]
method = args[11]
use = args[12]
by.study = args[13]

study.file = vector()
for(counter1 in 1:total_studies)
  study.file[counter1] = args[13+counter1] #name of study file
counter1 = 14 + counter1
counter1

study.names = vector()
for(counter2 in 1:total_studies)
  study.names[counter2] = args[counter1 + counter2] #name of study file
counter2 = counter1 + counter2

class.file = vector()
if(class_same_file == TRUE) class.file[1] = args[counter2+1] else
  for(counter3 in 1:total_studies)
    class.file[counter3] = args[counter3 + counter2] #name of study file
# counter3 = counter3+counter2

class.file

################################################################################

if(local=="home")
  path_dados="C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder" else
    if(local=="lab")
      path_dados="/media/thomasli/Lina-Andrey/raw-data-BcKO/DAPfinder" else
        if(local=="cgrb")
          path_dados="/capecchi/pharmacy/morgunlab/lina/imsys" 

path_final=paste(path_dados, newfoldername ,sep="/")

if(intersect.filter == FALSE) dir.create(path=path_final,showWarnings = TRUE)

#################################################################################################################

#functions
# 
# if(local=="home")
#   source("C:/Users/Lina/Dropbox/Tese/scripts/functions.R") else
#     if(local=="lab")
#       source("/media/thomasli/scripts/functions.R") else
#         if(local=="cgrb") 
#           source("/raid1/home/pharmacy/thomasli/scripts/functions.R")
# 
# print("pegou as funcoes")
# 
# table_cor_dif_blocks=function(cor,num_smp, pv_cor, class_type,method = "pearson",use = "pairwise.complete.obs",by.study = FALSE, study.names = "gr1"){
#   
#   if(!is.list(cor)) cor=list(list(cor))
#   if(length(cor) != length(study.names) & length(study.names) == 1) study.names = c(rep(study.names,length(study.names)))
#   genes=colnames(cor[[1]][[1]])
#   
#   comb=combn(rev(genes),2)
#   
#   final_table= data.frame(gene.x = rev(comb[2,]), gene.y = rev(comb[1,]))
#   
#   for(j in 1:length(cor))
#   {
#     for(i in 1:length(cor[[j]]))
#       final_table[,paste("cor",class_type[i],"_",study.names[j],sep="")] = cor[[j]][[i]][upper.tri(cor[[j]][[i]])]
#     if(by.study == TRUE) rm(cor) else cor[[j]][[i]] = NA
#     
#     if(length(pv_cor) != 0)
#     {
#       for(i in 1:length(cor[[j]]))
#         final_table[,paste("pv_cor",class_type[i],"_",study.names[j],sep="")] = pv_cor[[j]][[i]][upper.tri(pv_cor[[j]][[i]])]
#       if(by.study == TRUE) rm(pv_cor) else pv_cor[[j]][[i]] = NA
#     }
#     
#     for(i in 1:length(cor[[j]]))
#       final_table[,paste("num_smp_cor",class_type[i],"_",study.names[j],sep="")] = num_smp[[j]][[i]][upper.tri(num_smp[[j]][[i]])]
#     if(by.study == TRUE) rm(num_smp) else num_smp[[j]][[i]] = NA
#     
#     #difference of correlations in both conditions
#     final_table[,paste("difcor_",study.names[j],sep="")] = 
#       final_table[,paste("cor",class_type[1],"_",study.names[j],sep="")] - 
#       final_table[,paste("cor",class_type[2],"_",study.names[j],sep="")]
#     
#   }
#   if(by.study == FALSE) rm(cor, pv_cor, num_smp)
#   return(final_table)
# }
# 
# 
# #################################################################################################
# 
# setwd(path_dados)
# 
# data = list()
# 
# study.names = substr(study.file,1,nchar(study.file)-4) 
# print(study.names)
# 
# for(j in 1:total_studies)
#   data[[study.names[j]]]=read.table(study.file[j],row.names=1,header=TRUE)
# 
# print("leu os dados")
# 
# if(intersect.filter == FALSE)
#   data = lapply(data,filtering_missing,filtro_BRB=filtro_BRB,mp=mp) else
#   {
#     print("intersect.filter=TRUE")
#     print(path_final)
#     setwd(path_final)
#     genefilter=read.table(filterfile)
#     genefilter = data.frame(lapply(genefilter, as.character), stringsAsFactors=FALSE)
#     data=lapply(data, function(data) data[,colnames(data) %in% genefilter[,1]])
#   }
# print("filtrou os dados")
# 
# setwd(path_dados)
# 
# if(class_same_file==TRUE)
#   class.table = read.csv(class.file,header=TRUE) else
#   {
#     class.table=list()
#     for(j in 1:length(data))
#       class.table[[names(data)[j]]] = read.csv(class.file[j],header=TRUE)
#   }
# 
# if(!"class" %in% colnames(data[[1]]))
# {
#   data=lapply(data, function(data) merge(data,class.table[,c(1,3)],by.x=0,by.y=1,all.x=T))
#   for(j in 1: total_studies) row.names(data[[j]])=data[[j]][,"Row.names"]
#   
#   data=lapply(data, function(data) data=data[,-c(1)])
# }
# 
# #each item is a matrix with columns referring to a state classification
# if(class_same_file==TRUE) class_type = class.table[!duplicated(class.table[,3]),3] else
#   class_type = class.table[[1]][!duplicated(class.table[[1]][,3]),3]
# 
# data.list=separate.classes(data,study.file,class_type)
# 
# rm(data)
# 
# cor=lapply(data.list, function(data) lapply(data, cor, use=use, method=method))    
# 
# num_smp=lapply(data.list, function(data) lapply(data, num.pairwise.sample))
# 
# # pv_cor = cor.pv.pairwise3(cor,num_smp)
# 
# rm(data.list)
# 
# final_table = table_cor_dif_blocks(cor, num_smp, pv_cor=NULL, class_type, study.names = study.names)
# 
# setwd(path_final)
# 
# if(by.study == FALSE) 
#   write.table(final_table,"cor_dircor_table_byblocks.txt",row.names=F,sep="\t") else
#     write.table(final_table,paste("cor_dircor_table_byblocks_", study.names, ".txt", sep = ""),row.names=F,sep="\t")
# 
# 
# 
