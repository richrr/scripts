######################################################################################################

# Lina Thomas
# date: 1/13/14

########################################################################################################################################
# 
####################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < part1_blc1_filter_v2.R --no-save FALSE 0.3 "cgrb" 3 "Project_mp03" "filter.txt" "dados1.txt" "dados2.txt" "dados3.txt"  --no-save' -r filter_imsys  -f 12G -m 12G -q samwise
#
# SGE_Batch -c 'R < part1_blc1_filter_v2.R --no-save FALSE 0.3 "cgrb" 3 "Project_test" "filter.txt" "amostra1.txt" "amostra2.txt" "amostra3.txt"' -r filter_imsys -f 4G -m 4G -q samwise
# lab:    R < /home/thomasli/Dropbox/Tese/scripts/DAPfinder/blocks/part1_blc_filter_v2.R --no-save FALSE 0.3 "lab" 5 3 "Project_teste_mine" "filter_teste.txt" "number_studies_present_Mine.txt" FALSE NULL /home/thomasli/Dropbox/Tese/Undergraduate_Internship/Justin_Lina/Mine/tables/100_genes "/home/thomasli/Dropbox/Tese/Undergraduate_Internship/Justin_Lina/Mine/tables" Biewenga Igen Pyeon Scotto Zhai Biewenga_100Genes.csv Igen_100Genes.csv Pyeon_100Genes.csv Scotto_100Genes.csv Zhai_100Genes.csv Biewenga_genesymbol.txt Igen_genesymbol.txt pyeon_annotation.txt Scotto_Annotation.txt Zhai_genesymbol.txt

######################################################################################################
# arguments from terminal

args=commandArgs(trailingOnly = FALSE)

args

filtro_BRB = args[3]
mp = as.numeric(args[4])
local = args[5]
total_studies = as.numeric(args[6]) #number of datasets
min_studies = args[7]
newfoldername = args[8]
filter.file = args[9]
presence.file = args[10]
pv_filter = args[11]
gene_max = as.numeric(args[12])

path_dados = args[13]
path_annotation = args[14]

print(path_dados)
print(path_annotation)

study.names = vector()
for(j in 1 : total_studies)
  study.names[j] = args[14+j] #name of study file
counter1 = 14 + j


study.file = vector()
for(j in 1 : total_studies)
  study.file[j] = args[counter1+j] #name of study file
counter2 = counter1+j

print(study.file)

gene.annotation.file = vector()
for(j in 1 : total_studies)
  gene.annotation.file[j] = args[counter2+j] #name of study file

##############################################################################################

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
gene.annotation = list()

print(gene.annotation.file)

for(j in 1:total_studies)
{
  setwd(path_dados)
  if(substr(study.file[j], nchar(study.file[j])-3,nchar(study.file[j])) == ".csv")
    data[[j]]=read.table(study.file[j],row.names=1,header=TRUE, sep = ",", stringsAsFactors=F) else
      if(substr(study.file[j],nchar(study.file[j])-3,nchar(study.file[j])) == ".txt")
        data[[j]]=read.table(study.file[j], row.names = 1, header = TRUE, sep = "\t", stringsAsFactors=F)
  
  setwd(path_annotation)
  print(paste("j = ", j))
  if(substr(gene.annotation.file[j],nchar(gene.annotation.file[j])-3,nchar(gene.annotation.file[j])) == ".csv")
    gene.annotation[[j]]=read.table(gene.annotation.file[j], row.names = 1, header = TRUE, sep = ",", stringsAsFactors=F) else
      if(substr(gene.annotation.file[j], nchar(gene.annotation.file[j]) -3, nchar(gene.annotation.file[j])) == ".txt")
        gene.annotation[[j]]=read.table(gene.annotation.file[j],header=TRUE, sep="\t", stringsAsFactors=F)

  
}
print("leu os dados")

data = lapply (data,filtering_missing, filtro_BRB = filtro_BRB, mp = mp)
print("filtrou os dados")
data = lapply(data, function(x) t(x))
print("fez transposta")

for(j in 1:total_studies)
{
  gene.annotation[[j]]$probe.id = paste("X", gene.annotation[[j]]$probe.id, sep = "")
  
  data[[j]] = merge(data[[j]], gene.annotation[j], all.x = T, by.x = 0, by.y = "probe.id")
  colnames(data[[j]])[1] = "probe.id"
}

data = lapply(data, function(x) x[!is.na(x$Gene.symbol),])
print("filtrou gene symbol na")

# for(j in 1:total_studies)
#   print(data[[j]]$Gene.symbol)


if(total_studies>1)
{
  
   genes = c( data[[1]]$Gene.symbol, data[[2]]$Gene.symbol)
     
  if(total_studies>2)
    for(j in 3:total_studies)
    {
      genes = c(genes, data[[j]]$Gene.symbol)
#       print(genes)
    }
  
} else genes=colnames(data[[1]])

genes = genes[!duplicated(genes)]

genes = data.frame( Gene.symbol = genes, stringsAsFactors=F)

class(genes)

for(j in 1:total_studies)
  genes[,study.names[j]] = ifelse(genes$Gene.symbol %in% data[[j]]$Gene.symbol, 1, 0)


colnames(genes)
head(genes)

sum(genes[,c((ncol(genes)-4):ncol(genes))])

genes[,"number.of.studies.present"] = apply(genes[,c((ncol(genes)-4):ncol(genes))], 1, sum)

setwd(path_final)
write.table(genes, presence.file, row.names=F, col.names=F, sep = "\t")

genestofilter = genes$Gene.symbol[genes$number.of.studies.present >= min_studies]
genestofilter

  # 
  # if(pv_filter == TRUE)
  # {
  #   data=lapply(data, function(data) data[,colnames(data) %in% genes)
  #   setwd(path_dados)
  #   
  #   if(class_same_file==TRUE)
  #     class.table = read.csv(class.file,header=TRUE) else
  #     {
  #       class.table=list()
  #       for(j in 1:length(data))
  #         class.table[[names(data)[j]]] = read.csv(class.file[j],header=TRUE)
  #     }
  #   
  #   if(!"class" %in% colnames(data[[1]]))
  #   {
  #     data=lapply(data, function(data) merge(data,class.table[,c(1,3)],by.x=0,by.y=1,all.x=T))
  #     for(j in 1: total_studies) row.names(data[[j]])=data[[j]][,"Row.names"]
  #     
  #     data=lapply(data, function(data) data=data[,-c(1)])
  #   }
  #   
  #   #each item is a matrix with columns referring to a state classification
  #   if(class_same_file==TRUE) class_type = class.table[!duplicated(class.table[,3]),3] else
  #     class_type = class.table[[1]][!duplicated(class.table[[1]][,3]),3]
  #   
  #   data.list=separate.classes(data,study.file,class_type)
  #   
  #   data.list = lapply (data.list, filter.by.cv, gene_max)
  #     
  #   if(total_studies>1)
  #   {
  #     genes=intersect(colnames(data.list[[1]][[1]]),colnames(data.list[[2]][[1]]))
  #     
  #     if(total_studies>2)
  #       for(j in 3:total_studies)
  #         genes=intersect(genes,colnames(data[[j]][[1]]))
  #     
  #   } else genes=colnames(data.list[[1]][[1]])
  #   
  # }


setwd(path_final)
write.table(genestofilter,filter.file,row.names=F,col.names=F, sep = "\t")





