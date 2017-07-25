
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
# SGE_Batch -c 'R < scripts/blocks/part3_blc_filter_dirdifcor_v3.R --no-save "cgrb" 2 1 "Project_03mp_alias" "pearson" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer/Project_03mp_alias /capecchi/pharmacy/morgunlab/lina/CervicalCancer Biewenga Igen "cor_dircor_table_byblocks_mine_Biewenga.txt" "cor_dircor_table_byblocks_mine_Igen.txt"  Biewenga_genesymbol3.txt Igen_genesymbol3.txt Merged_difcor_table_Biewenga_Igen.txt' -r filter_difdifcor_mine_biewenga_igen -f 200G -m 200G -q samwise
# SGE_Batch -c 'R < scripts/blocks/part3_blc_filter_dirdifcor_v3.R --no-save "cgrb" 2 1 "Project_03mp_alias" "pearson" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer/Project_03mp_alias /capecchi/pharmacy/morgunlab/lina/CervicalCancer Pyeon Scotto "cor_dircor_table_byblocks_mine_Pyeon.txt" "cor_dircor_table_byblocks_mine_Scotto.txt"  Pyeon_genesymbol3.txt Scotto_genesymbol3.txt Merged_difcor_table_Pyeon_Scotto.txt' -r filter_difdifcor_mine_pyeon_scotto -f 200G -m 200G -q samwise

# tests
# SGE_Batch -c 'R < scripts/blocks/part3_blc_filter_dirdifcor_v3.R --no-save "cgrb" 5 3 "Project_03mp_teste" "pearson" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer/Project_teste_mine /capecchi/pharmacy/morgunlab/lina/CervicalCancer Pyeon Scotto Zhai Biewenga Igen "cor_dircor_table_byblocks_mine_Pyeon.txt" "cor_dircor_table_byblocks_mine_Scotto.txt" "cor_dircor_table_byblocks_mine_Zhai.txt" "cor_dircor_table_byblocks_mine_Biewenga.txt" "cor_dircor_table_byblocks_mine_Igen.txt" Pyeon_Annotation.txt Scotto_Annotation.txt Zhai_genesymbol.txt Biewenga_genesymbol.txt Igen_genesymbol.txt' -r filter_difdifcor_mine_test -f 64G -m 64G -q samwise
# at lab: R < Dropbox/Tese/scripts/DAPfinder/blocks/part3_blc_filter_dirdifcor_v3.R --no-save "lab" 5 3 3 "Project_03mp_teste" "pearson" TRUE /home/thomasli/Dropbox/Tese/Undergraduate_Internship/Justin_Lina/Mine/tables/Project_03mp_teste /home/thomasli/Dropbox/Tese/Undergraduate_Internship/Justin_Lina/Mine/tables Pyeon Scotto Zhai Biewenga Igen "cor_dircor_table_byblocks_mine_Pyeon.txt" "cor_dircor_table_byblocks_mine_Scotto.txt" "cor_dircor_table_byblocks_mine_Zhai.txt" "cor_dircor_table_byblocks_mine_Biewenga.txt" "cor_dircor_table_byblocks_mine_Igen.txt" pyeon_annotation.txt Scotto_Annotation.txt Zhai_genesymbol.txt Biewenga_genesymbol.txt Igen_Annotation.txt
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

local = args[3]
total_studies = as.numeric(args[4]) #number of datasets
min_studies = as.numeric(args[5])
newfoldername = args[6]
method = args[7]
by.study = args[8]

path_tables = args[9]
path_annotation = args[10]

study.names = vector()
for(counter1 in 1:total_studies)
  study.names[counter1] = args[counter1 + 10] #name of study file
counter1 = 10 + counter1

if(by.study == FALSE)
{
  table.file = args[counter1 + 1]
  counter3 = counter1 + 2 
} else
{
  table.file = vector()
  for(counter2 in 1:total_studies)
    table.file[counter2] = args[counter1+counter2] #name of study file
  counter2 = counter1 + counter2
  gene.annotation.file = vector()
  for(counter3 in 1 : total_studies)
    gene.annotation.file[counter3] = args[counter2 + counter3]
  counter3 = counter2 + counter3
}  
# table.file=c("cor_dircor_table_byblocks_mine_Scotto.txt", "cor_dircor_table_byblocks_mine_Zhai.txt", "cor_dircor_table_byblocks_mine_Biewenga.txt", "cor_dircor_table_byblocks_mine_Igen.txt")


out.file = args[counter3 + 1]
out.file
#stop()

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

path_annotation

final_table.list = list()

setwd(path_tables)
if(by.study == FALSE) 
  final_table.list[[1]] = read.table(table.file, header = T, stringsAsFactors = F, sep = "\t") else
  {
    for(j in 1:total_studies)
    {
      print(j)
      final_table.list[[j]] = read.table(table.file[j] , header = T, stringsAsFactors = F, sep = "\t")
    }
  }
# 

final_table = final_table.list[[1]]

if(by.study == TRUE)
{
  gene.annotation=list()
  for(j in 1:total_studies)
  {
#     print(colnames(final_table.list[[j]]))
    if(!"pair" %in% colnames(final_table.list[[j]]))
    {
      setwd(path_annotation)
      gene.annotation[[j]]=read.table(gene.annotation.file[j], header = TRUE, sep = "\t", stringsAsFactors = F)
#       print(colnames(gene.annotation[[j]]))
      colnames(final_table.list[[j]])[1:2]=c("probe.id.x", "probe.id.y")
#       print(colnames(final_table.list[[j]]))
      final_table.list[[j]] = merge(final_table.list[[j]],gene.annotation[[j]][,c("probe.id", "Gene.symbol.official")], by.x="probe.id.x", by.y = "probe.id", all.x = T)
      final_table.list[[j]] = final_table.list[[j]][!is.na(final_table.list[[j]]$Gene.symbol.official),]
      final_table.list[[j]] = merge(final_table.list[[j]],gene.annotation[[j]][,c("probe.id", "Gene.symbol.official")], by.x="probe.id.y", by.y = "probe.id", all.x = T)
      final_table.list[[j]] = final_table.list[[j]][!is.na(final_table.list[[j]]$Gene.symbol.official.y),]

      final_table.list[[j]]$pair = ifelse(final_table.list[[j]]$Gene.symbol.official.x < final_table.list[[j]]$Gene.symbol.official.y, paste(final_table.list[[j]]$Gene.symbol.official.x , final_table.list[[j]]$Gene.symbol.official.y), paste(final_table.list[[j]]$Gene.symbol.official.y , final_table.list[[j]]$Gene.symbol.official.x))

      separated_pair = strsplit(final_table.list[[j]]$pair, split = " ")

      separated_pair = matrix(unlist(separated_pair), ncol = 2,byrow = T)
      colnames(separated_pair) = c("Gene.symbol.x","Gene.symbol.y")
      final_table.list[[j]][, c("Gene.symbol.x", "Gene.symbol.y")] = separated_pair[, c("Gene.symbol.x", "Gene.symbol.y")]
      
      final_table.list[[j]] = final_table.list[[j]][ order( final_table.list[[j]]$Gene.symbol.x, final_table.list[[j]]$Gene.symbol.y, abs(final_table.list[[j]][ , paste("difcor_", study.names[j], sep="")]), decreasing = T),]
      final_table.list[[j]] = final_table.list[[j]][!duplicated(final_table.list[[j]]$pair),]
      


    }
  }

  print("merged tables")
      
  print(head(final_table.list[[1]]))


  final_table = final_table.list[[1]][,c(which(colnames(final_table.list[[1]]) == "pair"),grep("cor", colnames(final_table.list[[1]])))]
  print(head(final_table))

  for(j in 2:total_studies) 
    final_table = merge(final_table,final_table.list[[j]][,c(which(colnames(final_table.list[[j]]) == "pair"),grep("cor", colnames(final_table.list[[j]])))],by = "pair", all=T)

  setwd(path_tables)
  write.table(final_table, out.file, row.names = FALSE ,sep = "\t")

# still testing
  if(min_studies > 1)
  {
    print("min_studies > 1")
    final_table = filter.dircor.v2(final_table, c("Biewenga", "Igen", "Zhai", "Pyeon", "Scotto"))
    print("filtered same direction of difference of correlation")
    final_table = final_table[final_table$number.no.na >= min_studies, ]
    print(paste("filtered only pairs that appear in at least", min_studies, "datasets"))
    
    write.table(final_table, paste(substr(out.file, 1, nchar(out.file)-4), "_filtered.txt", sep=""), row.names = FALSE ,sep = "\t")
  }
  
    
        
} else
{
  final_table = filter.dircor.v2(final_table,study.names[1:total_studies])
  setwd(path_tables)
  write.table(final_table, "table_same_dir_difcor.txt", row.names = FALSE ,sep = "\t")

  final_table = final_table[final_table$number.no.na >= min_studies, ]
  write.table(final_table, "table_same_dir_difcor_minstudies3.txt", row.names = FALSE ,sep = "\t")

}
