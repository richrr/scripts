######################################################################################################

# Lina Thomas
# date: 1/22/14

########################################################################################################################################
# REQUIREMENTS: if there are two or more datasets with the same chip (gene.symbol) put them first in the command line!
# WARNING: Code only handles one group of genes from the same chip. For example, if you have datasets A and B from chip1
#          and datasets C and D from chip2 and want to compare the probe.id within those groups before comparing 
#          gene.symbols, gotta improve the code!
####################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < scripts/blocks/part3_blc_filter_dirdifcor_v2.R --no-save "cgrb" 5 3 3 "Project_mine_03mp" "pearson" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer/Project_mine_03mp /capecchi/pharmacy/morgunlab/lina/CervicalCancer Pyeon Scotto Zhai Biewenga Igen "cor_dircor_table_byblocks_mine_Pyeon.txt" "cor_dircor_table_byblocks_mine_Scotto.txt" "cor_dircor_table_byblocks_mine_Zhai.txt" "cor_dircor_table_byblocks_mine_Biewenga.txt" "cor_dircor_table_byblocks_mine_Igen.txt" Pyeon_Annotation.txt Scotto_Annotation.txt Zhai_genesymbol.txt Biewenga_genesymbol.txt Igen_genesymbol.txt "table_same_dir_difcor.txt"' -r filter_difdifcor_mine -f 256G -m 256G -q samwise
#
# tests
# SGE_Batch -c 'R < scripts/blocks/part3_blc_filter_dirdifcor_v2.R --no-save "cgrb" 5 3 3 "Project_teste_mine" "pearson" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer/Project_teste_mine /capecchi/pharmacy/morgunlab/lina/CervicalCancer Pyeon Scotto Zhai Biewenga Igen "cor_dircor_table_byblocks_mine_Pyeon.txt" "cor_dircor_table_byblocks_mine_Scotto.txt" "cor_dircor_table_byblocks_mine_Zhai.txt" "cor_dircor_table_byblocks_mine_Biewenga.txt" "cor_dircor_table_byblocks_mine_Igen.txt" Pyeon_Annotation.txt Scotto_Annotation.txt Zhai_genesymbol.txt Biewenga_genesymbol.txt Igen_genesymbol.txt "table_same_dir_difcor.txt"' -r filter_difdifcor_mine_test -f 4G -m 4G -q samwise
# at lab: R < Dropbox/Tese/scripts/DAPfinder/blocks/part3_blc_filter_dirdifcor_v2.R --no-save "lab" 5 3 3 "Project_teste_mine" "pearson" TRUE /home/thomasli/Dropbox/Tese/Undergraduate_Internship/Justin_Lina/Mine/tables/100_genes/Project_teste_mine /home/thomasli/Dropbox/Tese/Undergraduate_Internship/Justin_Lina/Mine/tables Pyeon Scotto Zhai Biewenga Igen "cor_dircor_table_byblocks_mine_Pyeon.txt" "cor_dircor_table_byblocks_mine_Scotto.txt" "cor_dircor_table_byblocks_mine_Zhai.txt" "cor_dircor_table_byblocks_mine_Biewenga.txt" "cor_dircor_table_byblocks_mine_Igen.txt" pyeon_annotation.txt Scotto_Annotation.txt Zhai_genesymbol.txt Biewenga_genesymbol.txt Igen_Annotation.txt "merged_table.txt"
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

local = args[3]
total_studies = as.numeric(args[4]) #number of datasets
min_studies = as.numeric(args[5])
same.gene.symbol = as.numeric(args[6]) #how many datasets have the same chip? If they are all from different chips write 1
newfoldername = args[7]
method = args[8]
by.study = args[9]

path_tables = args[10]
path_annotation = args[11]

study.names = vector()
for(counter1 in 1:total_studies)
  study.names[counter1] = args[counter1 + 11] #name of study file
counter1 = 11 + counter1

if(by.study == FALSE)
{
  table.file = args[counter1 + 1]
  counter2 = counter1 + 1 
} else
{
  table.file = vector()
  for(counter2 in 1:total_studies)
    table.file[counter2] = args[counter1+counter2] #name of study file
  counter2 = counter1 + counter2
} 
# table.file=c("cor_dircor_table_byblocks_mine_Scotto.txt", "cor_dircor_table_byblocks_mine_Zhai.txt", "cor_dircor_table_byblocks_mine_Biewenga.txt", "cor_dircor_table_byblocks_mine_Igen.txt")

gene.annotation.file = vector()
for(counter3 in 1 : total_studies)
  gene.annotation.file[counter3] = args[counter2 + counter3]
counter3 = counter2 + counter3

final_table.file = args[counter3 + 1]
gene.annotation.file

##############################################################################################
# 
# if(!"psych" %in% installed.packages()) install.packages("psych")
# library("psych")

#################################################################################################################

#functions

if(local=="home")
  source("C:/Users/Lina/Dropbox/Tese/scripts/functions.R") else
    if(local=="lab")
      source("/home/thomasli/Dropbox/Tese/scripts/functions.R") else
        if(local=="cgrb") 
          source("/raid1/home/pharmacy/thomasli/scripts/functions.R")

print("pegou as funcoes")


#################################################################################################

path_tables
table.file

final_table.list = list()
if(by.study == FALSE) 
  final_table.list[[1]] = read.table(table.file, header = T) else
  {
    gene.annotation = list()
    for(j in 1:total_studies)
    {
      print(j)
      setwd(path_tables)
      final_table.list[[j]] = read.table(table.file[j] , header = T, stringsAsFactors = F, sep = "\t")

      if(length(grep("gene.x", colnames(final_table.list[[j]]))) != 0)
      {
        final_table.list[[j]]$gene.x = substr(final_table.list[[j]]$gene.x, 2, nchar(final_table.list[[j]]$gene.x))
        final_table.list[[j]]$gene.y = substr(final_table.list[[j]]$gene.y, 2, nchar(final_table.list[[j]]$gene.y))
        colnames(final_table.list[[j]])[which(colnames(final_table.list[[j]]) == "gene.x")] = "probe.id.x"
        colnames(final_table.list[[j]])[which(colnames(final_table.list[[j]]) == "gene.y")] = "probe.id.y"
      }

      if(length(grep("Gene.symbol", colnames(final_table.list[[j]]))) == 0)
      {
        setwd(path_annotation)
        gene.annotation[[j]]=read.table(gene.annotation.file[j], header = TRUE, sep = "\t", stringsAsFactors = F)
      }
    }
  }
# 

final_table = final_table.list[[1]]

by.study == TRUE
same.gene.symbol < total_studies

grep("Gene.symbol", colnames(final_table.list[[j]]))



if(by.study == TRUE)
{
  if(same.gene.symbol < total_studies)
    for(j in 1:total_studies)
    if(length(grep("Gene.symbol", colnames(final_table.list[[j]]))) == 0)
    {
      final_table.list[[j]] = merge(final_table.list[[j]],gene.annotation[[j]][,c("probe.id", "Gene.symbol")], by.x="probe.id.x", by.y = "probe.id", all.x = T)
      final_table.list[[j]] = final_table.list[[j]][!is.na(final_table.list[[j]]$Gene.symbol),]
      final_table.list[[j]] = merge(final_table.list[[j]],gene.annotation[[j]][,c("probe.id", "Gene.symbol")], by.x="probe.id.y", by.y = "probe.id", all.x = T)
      final_table.list[[j]] = final_table.list[[j]][!is.na(final_table.list[[j]]$Gene.symbol.y),]
    }
  
  if(grep("cor_dircor", table.file) == total_studies)
    final_table.list = lapply(final_table.list, ordpair.x_v2, TRUE) else
      for(j in 1:total_studies)
	if(length(grep("cor_dircor", table.file[j]))>0)
	  final_table.list[[j]] = ordpair.x_v2(final_table.list[[j]], TRUE)
  
  final_table = final_table.list[[1]]
#   final_table = merge(final_table,final_table.list[[j]],by = intersect(c("probe.id.x","probe.id.y"),c("probe.id.x","probe.id.y")),all=T)
    
#   for(j in 2:4) 
  for(j in 2:total_studies) 
  { 
    if(same.gene.symbol >= 2 & j <= same.gene.symbol)
      final_table = merge(final_table,final_table.list[[j]],by = intersect(c("probe.id.x","probe.id.y","Gene.symbol.x","Gene.symbol.y"),c("probe.id.x","probe.id.y", "Gene.symbol.x","Gene.symbol.y")),all=T)
    
    if( same.gene.symbol >= 1 & j > same.gene.symbol )
    { 
      probeid = grep("probe.id", colnames(final_table))
      if(length(probeid) > 0)
        final_table = final_table[, - probeid]
        
      final_table = merge(final_table,final_table.list[[j]],by = intersect(c("Gene.symbol.x","Gene.symbol.y"),c("Gene.symbol.x","Gene.symbol.y")),all=T)
    }
    print("fez merge com tabela ",j)
    setwd(path_tables)
    write.table(final_table, final_table.file , row.names = FALSE ,sep = "\t") 
 
  }
}

#final_table = filter.dircor.v2(final_table,study.names[1:j])

#write.table(final_table, "table_same_dir_difcor.txt", row.names = FALSE ,sep = "\t")


#final_table = final_table[final_table$number.no.na >= min_studies, ]

#nrow(final_table)
# stop()  


#final_table = final_table[final_table$number.no.na >= min_studies, ]
#write.table(final_table, "table_same_dir_difcor_minstudies3.txt", row.names = FALSE ,sep = "\t")


