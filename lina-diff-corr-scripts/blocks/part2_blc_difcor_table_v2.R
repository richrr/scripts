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
# SGE_Batch -c 'R < ~/scripts/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "cgrb" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer  Biewenga_Expression5.txt Biewenga "Mine_Class_all.txt" Biewenga_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt"' -r cor_difcor_biewenga -f 100G -m 100G -q samwise
# SGE_Batch -c 'R < ~/scripts/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "cgrb" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer  Igen_Expression5.txt Igen "Mine_Class_all.txt" Igen_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt"' -r cor_difcor_igen -f 100G -m 100G -q samwise
# SGE_Batch -c 'R < ~/scripts/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "cgrb" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer  Pyeon_Expression5.txt Pyeon "Mine_Class_all.txt" Pyeon_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt"' -r cor_difcor_pyeon -f 120G -m 120G -q samwise
# SGE_Batch -c 'R < ~/scripts/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "cgrb" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer  Scotto_Expression5.txt Scotto "Mine_Class_all.txt" Scotto_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt"' -r cor_difcor_scotto -f 120G -m 120G -q samwise
# SGE_Batch -c 'R < ~/scripts/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "cgrb" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer  Zhai_Expression5.txt Zhai "Mine_Class_all.txt" Zhai_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt"' -r cor_difcor_zhai -f 120G -m 120G -q samwise
# tests
# SGE_Batch -c 'R < ~/scripts/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "cgrb" 5 TRUE "Project_teste_mine" "filter_mine.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer /capecchi/pharmacy/morgunlab/lina/CervicalCancer  Biewenga_100Genes.csv Igen_100Genes.csv Pyeon_100Genes.csv Scotto_100Genes.csv Zhai_100Genes.csv Biewenga Igen Pyeon Scotto Zhai "Mine_Class_all_v2.txt" Biewenga_genesymbol.txt Igen_genesymbol.txt Pyeon_Annotation.txt Scotto_Annotation.txt Zhai_genesymbol.txt "cor_dircor_table_byblocks_mine.txt"' -r cor_difcor_imsys_test -f 4G -m 4G -q samwise
# at lab: 
# R < Dropbox/Tese/scripts/DAPfinder/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "lab" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /home/thomasli/Dropbox/Tese/Dados/Mine/tables "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/media/Storage/CervicalCancer/Project_03mp_alias" Biewenga_Expression5.txt Biewenga "Mine_Class_all.txt" Biewenga_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt"
# R < Dropbox/Tese/scripts/DAPfinder/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "lab" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /home/thomasli/Dropbox/Tese/Dados/Mine/tables "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/media/Storage/CervicalCancer/Project_03mp_alias"  Igen_Expression5.txt Igen "Mine_Class_all.txt" Igen_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt"
# nohup R < Dropbox/Tese/scripts/DAPfinder/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "lab" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /home/thomasli/Dropbox/Tese/Dados/Mine/tables "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/media/Storage/CervicalCancer/Project_03mp_alias" Pyeon_Expression5.txt Pyeon "Mine_Class_all.txt" Pyeon_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt" &
# nohup R < Dropbox/Tese/scripts/DAPfinder/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "lab" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /home/thomasli/Dropbox/Tese/Dados/Mine/tables "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/media/Storage/CervicalCancer/Project_03mp_alias" Scotto_Expression5.txt Scotto "Mine_Class_all.txt" Scotto_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt" &
# nohup R < Dropbox/Tese/scripts/DAPfinder/blocks/part2_blc_difcor_table_v2.R --no-save FALSE NULL "lab" 1 TRUE "Project_03mp_alias" "filter.txt" TRUE "pearson" "pairwise.complete.obs" TRUE /home/thomasli/Dropbox/Tese/Dados/Mine/tables "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/home/thomasli/Dropbox/Tese/Dados/Mine/tables" "/media/Storage/CervicalCancer/Project_03mp_alias" Zhai_Expression5.txt Zhai "Mine_Class_all.txt" Zhai_genesymbol3.txt "cor_dircor_table_byblocks_mine.txt" &



######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

filtro_BRB = args[3]
mp = as.numeric(args[4]) # caso queira passar filtro missing (intersect.filter=FALSE)
local = args[5]
total_studies = as.numeric(args[6]) #number of datasets
class_same_file = args[7] # TRUE quando informacao da classe de todos os estudos estao em um arquivo
newfoldername = args[8]
filter.file = args[9]
intersect.filter = args[10] # = FALSE caso nao quer usar filtro da parte 1 e TRUE cc
method = args[11] # tem que por "pearson"
use = args[12] # por "pairwise.complete.obs"
by.study = args[13] # TRUE quando salva estudos separado. Sempre por TRUE

path_dados = args[14]
path_annotation = args[15]
path_class = args[16]

study.file = vector()
for(counter1 in 1:total_studies)
  study.file[counter1] = args[16+counter1] #name of study file
counter1 = 16 + counter1
counter1

study.names = vector()
for(counter2 in 1:total_studies)
  study.names[counter2] = args[counter1 + counter2] #name of study file
counter2 = counter1 + counter2

class.file = vector()
if(class_same_file == TRUE) 
{
  class.file[1] = args[counter2+1]
  counter3=counter2+1
} else 
{
  for(counter3 in 1:total_studies)
    class.file[counter3] = args[counter3 + counter2] #name of study file
  counter3 = counter2 + counter3
}
# counter3 = counter3+counter2

gene.annotation.file = vector()
for(counter4 in 1 : total_studies)
  gene.annotation.file[counter4] = args[counter3 + counter4]
counter4 = counter3 + counter4

gene.annotation.file

out.file = args[counter4+1] # nome da tabela de saida

out.file

################################################################################

path_final=paste(path_dados, newfoldername ,sep="/")

if(intersect.filter == FALSE) dir.create(path=path_final,showWarnings = TRUE)

#################################################################################################################

#functions

if(local=="home")
  source("C:/Users/Lina/Dropbox/Tese/scripts/functions.R") else
    if(local=="lab")
      source("/home/thomasli/Dropbox/Tese/scripts/functions.R") else
        if(local=="cgrb") 
          source("/raid1/home/pharmacy/thomasli/scripts/functions.R")

print("pegou as funcoes")

table_cor_dif_blocks_v2 = function(cor, num_smp, pv_cor, class_type,method = "pearson",use = "pairwise.complete.obs", by.study = FALSE, study.names = "gr1"){
  
  if(!is.list(cor)) cor=list(list(cor))
  if(!is.list(num_smp)) cor=list(list(num_smp))
  
  if(length(cor) != length(study.names) & length(study.names) == 1) study.names = c(rep(study.names,length(study.names)))
  
  if(by.study == FALSE)
  {
    genes=colnames(cor[[1]][[1]])
    
    comb=combn(rev(genes),2)
    
    final_table= data.frame(gene.x = rev(comb[2,]), gene.y = rev(comb[1,]))
    
    for(j in 1:length(cor))
    {
      for(i in 1:length(cor[[j]]))
      {
        final_table[,paste("cor",class_type[i],"_",study.names[j],sep="")] = cor[[j]][[i]][upper.tri(cor[[j]][[i]])]
        
        if(by.study == TRUE) rm(cor) else cor[[j]][[i]] = NA
      }
      if(length(pv_cor) != 0)
        for(i in 1:length(cor[[j]]))
        {          
          final_table[,paste("pv_cor",class_type[i],"_",study.names[j],sep="")] = pv_cor[[j]][[i]][upper.tri(pv_cor[[j]][[i]])]
          if(by.study == TRUE) rm(pv_cor) else pv_cor[[j]][[i]] = NA
        }
      
      for(i in 1:length(cor[[j]]))
      {
        final_table[,paste("num_smp_cor",class_type[i],"_",study.names[j],sep="")] = num_smp[[j]][[i]][upper.tri(num_smp[[j]][[i]])]
        if(by.study == TRUE) rm(num_smp) else num_smp[[j]][[i]] = NA
      }
      #difference of correlations in both conditions
      final_table[,paste("difcor_",study.names[j],sep="")] = 
        final_table[,paste("cor",class_type[1],"_",study.names[j],sep="")] - 
        final_table[,paste("cor",class_type[2],"_",study.names[j],sep="")]
      
    }
    
  }  else
  {
    genes=colnames(cor[[1]])
    
    comb=combn(rev(genes),2)
    
    final_table= data.frame(gene.x = rev(comb[2,]), gene.y = rev(comb[1,]))
    
    for(i in 1:length(cor))
    {
      final_table[,paste("cor",class_type[i],"_",study.names,sep="")] = cor[[i]][upper.tri(cor[[i]])]
          
      if(length(pv_cor) != 0)
        final_table[,paste("pv_cor",class_type[i],"_",study.names,sep="")] = pv_cor[[i]][upper.tri(pv_cor[[i]])]
    
      final_table[,paste("num_smp_cor",class_type[i],"_",study.names,sep="")] = num_smp[[i]][upper.tri(num_smp[[i]])]
      
#       print(paste("num_smp_cor",class_type[i],"_",study.names,sep=""))
#       print(num_smp[[i]][upper.tri(num_smp[[i]])])
      
    }
    #difference of correlations in both conditions
    final_table[,paste("difcor_",study.names,sep="")] = 
      final_table[,paste("cor",class_type[1],"_",study.names,sep="")] - 
      final_table[,paste("cor",class_type[2],"_",study.names,sep="")]
    
  }
    
  return(final_table)
}


#################################################################################################

setwd(path_dados)

data = list()
gene.annotation = list()


for(j in 1:total_studies)
{
  setwd(path_dados)
  if(substr(study.file[j], nchar(study.file[j])-3,nchar(study.file[j])) == ".csv")
    data[[j]]=t(read.table(study.file[j],row.names=1,header=TRUE, sep = ",", stringsAsFactors=F)) else
      if(substr(study.file[j],nchar(study.file[j])-3,nchar(study.file[j])) == ".txt")
        data[[j]]=t(read.table(study.file[j], row.names = 1, header = TRUE, sep = "\t", stringsAsFactors=F))
  

  setwd(path_annotation)
  print(paste("j = ", j))
  if(substr(gene.annotation.file[j], nchar(gene.annotation.file[j]) - 3, nchar(gene.annotation.file[j])) == ".csv")
    gene.annotation[[j]]=read.table(gene.annotation.file[j], row.names = 1, header = TRUE, sep = ",", stringsAsFactors = F) else
      if(substr(gene.annotation.file[j], nchar(gene.annotation.file[j]) -3, nchar(gene.annotation.file[j])) == ".txt")
        gene.annotation[[j]]=read.table(gene.annotation.file[j], header = TRUE, sep = "\t", stringsAsFactors = F)

}

print("leu os dados")

if(intersect.filter == FALSE)
  data = lapply(data,filtering_missing,filtro_BRB=filtro_BRB,mp=mp) else
  {
    print("intersect.filter=TRUE")
    print(path_final)
    setwd(path_final)
    
    if(substr(filter.file,nchar(filter.file)-3,nchar(filter.file)) == ".csv")
      genefilter=read.table(filter.file, sep = ",", stringsAsFactors = F) else
        if(substr(filter.file,nchar(filter.file)-3,nchar(filter.file)) == ".txt")
          genefilter=read.table(filter.file, stringsAsFactors = F)
    
    print(ncol(genefilter))
 
    probe.id.in.genefilter = lapply(gene.annotation, function(x) x$probe.id[ x$Gene.symbol.official %in% genefilter[,1]])
#     probe.id.in.genefilter = lapply(probe.id.in.genefilter, function(x) paste("X", x, sep=""))

    for(j in 1:total_studies) # escolha prob.id com gene symbol de genefilter
      data[[j]] = data[[j]][, colnames(data[[j]]) %in% probe.id.in.genefilter[[j]]]

  }

print("filtrou os dados")

# 

setwd(path_class)

if(class_same_file==TRUE)
{
  if(substr(class.file,nchar(class.file)-3,nchar(class.file)) == ".csv")
    class.table=read.table(class.file, sep = ",", stringsAsFactors = F, header = TRUE) else
      if(substr(class.file,nchar(class.file)-3,nchar(class.file)) == ".txt")
        class.table=read.table(class.file, header = TRUE, stringsAsFactors = F, sep = "\t")
} else
  {
    class.table=list()
    for(j in 1:length(data))
      class.table[[names(data)[j]]] = read.csv(class.file[j],header=TRUE)
  }
  
head(class.table)

if(!"class" %in% colnames(data[[1]]))
{
  data=lapply(data, function(data) merge(data,class.table[,c("Unique.id","class")],by.x=0,by.y="Unique.id",all.x=T))
  for(j in 1: total_studies) row.names(data[[j]])=data[[j]][,"Row.names"]
  
  data=lapply(data, function(data) data=data[,-c(1)])
}


#each item is a matrix with columns referring to a state classification
if(class_same_file==TRUE) class_type = class.table[!duplicated(class.table[,"class"]),"class"] else
  class_type = class.table[[1]][!duplicated(class.table[[1]][,"class"]),"class"]

class_type = class_type[order(class_type)]
class_type

data.list=separate.classes2(data,study.names,class_type)

rm(data)

cor=lapply(data.list, function(data) lapply(data, cor, use=use, method=method))    

num_smp=lapply(data.list, function(data) lapply(data, num.pairwise.sample))

pv_cor = cor.pv.pairwise3(cor,num_smp)

rm(data.list)

if(by.study == T)
{
  final_table = list()
  for(j in 1:total_studies)
  {
    final_table = table_cor_dif_blocks_v2(cor[[j]], num_smp[[j]], pv_cor = NULL, class_type, study.names = study.names[j], by.study = by.study)
    print(colnames(final_table))
    
    setwd(path_final)
    write.table(final_table, paste(substr(out.file, 1, nchar(out.file) - 4 ), "_",
                                  study.names[j], substr(out.file, nchar(out.file) - 3, nchar(out.file)), sep = ""), row.names = F, sep = "\t")
    
  } 
} else   
{
  final_table = lapply(cor, function(cor) table_cor_dif_blocks(cor, num_smp, pv_cor=NULL, class_type, study.names = study.names))
  setwd(path_final)
  write.table(final_table,out.file,row.names=F,sep="\t")
}

