######################################################################################################

# Lina Thomas
# date: 1/22/14

########################################################################################################################################
# 
####################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < scripts/blocks/part3_blc_filter_dirdifcor.R --no-save "cgrb" 3 "Project_03mp" "pearson" TRUE "dados1" "dados2" "dados3" "cor_dircor_table_byblocks_dados1.txt" "cor_dircor_table_byblocks_dados2.txt" "cor_dircor_table_byblocks_dados3.txt"' -r filter_difdifcor_imsys -f 128G -m 128G -q samwise
#
# tests
# SGE_Batch -c 'R < scripts/blocks/part3_blc_filter_dirdifcor.R --no-save "cgrb" 3 "Project_test" "pearson" TRUE "amostra1" "amostra2" "amostra3" "cor_dircor_table_byblocks_amostra1.txt" "cor_dircor_table_byblocks_amostra2.txt" "cor_dircor_table_byblocks_amostra3.txt"' -r filter_difdifcor_test -f 4G -m 4G -q samwise
# at lab: R < Dropbox/Tese/scripts/DAPfinder/blocks/part3_blc_filter_dirdifcor.R --no-save "lab" 3 "Project_teste" "pearson" TRUE "amostra1" "amostra2" "amostra3" "cor_dircor_table_byblocks_amostra1.txt" "cor_dircor_table_byblocks_amostra2.txt" "cor_dircor_table_byblocks_amostra3.txt"
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

local = args[3]
total_studies = as.numeric(args[4]) #number of datasets
newfoldername = args[5]
method = args[6]
by.study = args[7]

study.names = vector()
for(counter1 in 1:total_studies)
  study.names[counter1] = args[counter1 + 7] #name of study file
counter1 = 7 + counter1

if(by.study == FALSE)
  table.file = args[counter1 + 1] else
  {
    table.file = vector()
    for(counter2 in 1:total_studies)
      table.file[counter2] = args[counter1+counter2] #name of study file
    counter2 = counter1 + counter2
  } 
    

##############################################################################################
# 
# if(!"psych" %in% installed.packages()) install.packages("psych")
# library("psych")

################################################################################

if(local=="home")
  path_dados="C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder" else
    if(local=="lab")
      path_dados="/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder" else
        if(local=="cgrb")
          path_dados="/capecchi/pharmacy/morgunlab/lina/imsys" 

path_final=paste(path_dados, newfoldername ,sep="/")

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

setwd(path_final)

print(table.file)

final_table.list = list()
if(by.study == FALSE) 
  final_table.list[[1]] = read.table(table.file, header = T) else
    for(j in 1:total_studies)
      final_table.list[[j]] = read.table(table.file[j] , header = T)

final_table = final_table.list[[1]]

if(by.study == TRUE)
  for(j in 2:total_studies) 
    final_table = merge(final_table,final_table.list[[j]],by = intersect(c("gene.x","gene.y"),c("gene.x","gene.y")),all=T)


final_table = filter.dircor(final_table,study.names)

write.table(final_table, "table_same_dir_difcor.txt", row.names = FALSE ,sep = "\t")




