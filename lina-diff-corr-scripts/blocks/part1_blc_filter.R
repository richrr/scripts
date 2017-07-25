######################################################################################################

# Lina Thomas
# date: 1/13/14

########################################################################################################################################
# 
####################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < part1_blc1_filter.R --no-save FALSE 0.3 "cgrb" 3 "Project_mp03" "filter.txt" "dados1.txt" "dados2.txt" "dados3.txt"  --no-save' -r filter_imsys  -f 12G -m 12G -q samwise
#
# SGE_Batch -c 'R < part1_blc1_filter.R --no-save FALSE 0.3 "cgrb" 3 "Project_test" "filter.txt" amostra1.txt" "amostra2.txt" "amostra3.txt"' -r filter_imsys -f 4G -m 4G -q samwise
######################################################################################################
# arguments from terminal

args=commandArgs(trailingOnly = FALSE)

args

filtro_BRB=args[3]
mp=as.numeric(args[4])
local = args[5]
total_studies=as.numeric(args[6]) #number of datasets
newfoldername=args[7]
filterfile=args[8]

study.file=vector()
for(j in 1:total_studies)
  study.file[j]=args[8+j] #name of study file

print(study.file)
##############################################################################################

if(local=="home")
  path_dados="C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder" else
    if(local=="lab")
      path_dados="/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder" else
        if(local=="cgrb")
          path_dados="/capecchi/pharmacy/morgunlab/lina/imsys" 

path_final=paste(path_dados, newfoldername ,sep="/")
dir.create(path=path_final,showWarnings = TRUE)

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

setwd(path_dados)

data=list()

for(j in 1:total_studies)
{
  data[[j]] = t(read.table(study.file[j], row.names = 1, header = TRUE, sep = "\t", stringsAsFactors=F))
  if("class" %in% colnames(data[[j]])) data[[j]]=data[[j]][,-which(colnames(data[[j]])=="class")] 
}
print("leu os dados")

data=lapply (data,filtering_missing,filtro_BRB=filtro_BRB,mp=mp)
print("filtrou os dados")

if(total_studies>1)
{
  genes=intersect(colnames(data[[1]]),colnames(data[[2]]))
  
  if(total_studies>2)
    for(j in 3:total_studies)
      genes=intersect(genes,colnames(data[[j]]))
  
} else genes=colnames(data[[1]])

setwd(path_final)
write.table(genes,filterfile,row.names=F,col.names=F)





