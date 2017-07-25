######################################################################################################

# Lina Thomas
# date: 12/13/13

# REQUIREMENT: The probe.ids must be in the same order in each study expression table!!!!!!!!!!!!!!!!
########################################################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < scripts/degs_correlation_v4.R --no-save $SGE_TASK_ID 3000 "pearson" "cgrb" TRUE FALSE "all" "degs_gbacc.txt" 3 "dados1.txt" "dados2.txt" "dados3.txt" "dados1" "dados2" "dados3" ' -r degs_correlation_imsys -t 1 -f 8G -m 8G -q samwise
#
# R < /home/thomasli/Dropbox/Tese/scripts/Networks/correlation/degs_correlation_v4.R --no-save 1 3000 "pearson" "lab" TRUE FALSE "all" "degs_gbacc.txt" 3 "dados1.txt" "dados2.txt" "dados3.txt" "dados1" "dados2" "dados3"  
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

print(args)

pi = as.numeric(args[3]) #which part
p = as.numeric(args[4]) #how many parts
method = args[5]
local = args[6]
Parallel = args[7] #running parallel or alltogether
pv_filter = args[8]

state = args[9] # case, control or all
if(state %in% c("case", "control"))
{
  print("todo: option of case and control separate")
  print("WARNING: considering state = all")
  state = "all"
}

if(state == "all") state = c("case", "control")

degs_orig = args[10] # file containing degs information

total_studies = as.numeric(args[11]) #number of datasets

# by.study=args[12]

study.file = vector()
for(j in 1:total_studies)
 study.file[j] = args[11 + j] #name of study file
counter1 = 11 + j
  
study.names = vector()
for(j in 1:total_studies)
  study.names[j] = args[counter1 + j] #name of study file
counter2 = counter1 + j

####################################################################################################

if(local=="home")
  path_dados="C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder" else
    if(local=="lab")
      path_dados="/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder" else
        if(local=="cgrb")
          path_dados="/capecchi/pharmacy/morgunlab/lina/imsys" 

if(local=="home")
  path_degs="C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DEGs" else
    if(local=="lab")
      path_degs="/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DEGs" else
        if(local=="cgrb")
          path_degs="/capecchi/pharmacy/morgunlab/lina/imsys/degs" 

if(local=="home")
  path="C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/AE_data" else
    if(local=="lab")
      path="/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/AE_data" else
        if(local=="cgrb")
          path="/capecchi/pharmacy/morgunlab/lina/ArrayExpress" 

#################################################################################################################

#functions

# source("/home/thomasli/Dropbox/Tese/scripts/functions.R")

if(local=="home")
  source("C:/Users/Lina/Dropbox/Tese/scripts/functions.R") else
    if(local=="lab")
      source("/home/thomasli/Dropbox/Tese/scripts/functions.R") else
        if(local=="cgrb") 
          source("/raid1/home/pharmacy/thomasli/scripts/functions.R")

table.perpair.cor.pv = function(data, pairs, method = "pearson", use = "pairwise.complete.obs"){
  
  final_table = data.frame(pairs) 
  colnames(final_table) = c("ID.x","ID.y")
  
  for(pair_index in 1:nrow(pairs))
  {
    
    temp = cor.perpair(data, pair = pairs[pair_index ,], method = method, use = use)
    
    for(j in 1: total_studies)
      final_table[pair_index, "cor"] = temp[1, 1]
    for(j in 1:total_studies)
      final_table[pair_index, "pv"] = temp[2, 1]
    
    if(method == "pearson")
      for(j in 1: total_studies)
        final_table[pair_index, "pairwise_num_smp"] = temp[3,1]
  }

  for(j in 1: total_studies)
    final_table[, "dircor"] = ifelse( final_table[, "cor"] > 0, 1, -1)   
  
  return(final_table)
  
}


cor.pv.ftb=function(data,n=1,gene.symbol=NULL){
  
  cor=lapply(data, function(x) cor(x,use="pairwise.complete.obs"))
  
  min.mat= pairwise_num_obs(data, n)
  
  pvalue=cor.pv.pairwise(cor,min.mat,n)
  
  if(n==1) dados_finais=final_table(cor,pvalue,min.mat,gene.symbol=gene.symbol) else
    dados_finais=final_table(cor,pvalue,min.mat,n,gene.symbol)
  
  
  return(dados_finais)
}

#################################################################################################################
setwd(path_degs)

if(substr(degs_orig,nchar(degs_orig)-3,nchar(degs_orig)) == ".csv")
  degs=read.csv(degs_orig,header=TRUE) else
    if(substr(degs_orig,nchar(degs_orig)-3,nchar(degs_orig)) == ".txt")
      degs =read.table(degs_orig,header=TRUE) else
      {
        print("File type different from csv or txt. Please select either one.")
        stop()
      }

degs[,1]=paste("X",degs[,1],sep="")

print(head(degs))

setwd(path_dados)

dados=list()

for(j in 1:total_studies)
  dados[[j]]=read.table(study.file[j],row.names=1,header=TRUE)
names(dados) = study.names

class_type = dados[[j]][!duplicated(dados[[j]][,"class"]), "class"]
class_type

if(p != 1)
{
  genes = degs[,"probe.id"]
  pairs = t(combn(rev(genes),2))
  
  row.names(pairs)=paste(pairs[,1],pairs[,2])
  pair=subsetting(row.names(pairs),p,pi)
  
  pairs=pairs[which(row.names(pairs) %in% pair),c(1:2)]
  degs_selected = c(pairs[,1], pairs[,2])[!duplicated(c(pairs[,1], pairs[,2]))]
  
} else if(p == 1 & pi == 1)
  degs_selected = degs[, "probe.id"] else
  {
    print("p == 1, but pi greater than 1.")
    stop()
  }
  
print(degs_selected)

# gives a list with degs expression only where components are case and control
dados_degs = lapply(dados, function(x) x[,colnames(x) %in% degs_selected])

for(j in 1:total_studies)
  if("class" %in% colnames(dados[[1]]))
    dados_degs[[j]][, "class"] = dados[[j]][, "class"] else 
      print("todo: get class names from file")

dados_degs = separate.classes2(dados_degs, study.names, class_type=c(0,1))

names(dados_degs)
names(dados_degs[[1]])
head(dados_degs[[1]][[1]])

dif_median_degs = dif.median.classes(dados_degs)
head(dif_median_degs[[1]])

regulation = lapply(dif_median_degs, function(x) ifelse(x > 0 ,"U","D"))
head(regulation[[1]])

if(Parallel==TRUE)
{
  table_list = lapply(dados_degs, function(x) lapply(x, table.perpair.cor.pv, pairs = pairs, method = method ))
  
  print(head(table-list[[1]][[1]]))
  
  dados_finais = as.data.frame(table_list)
  colnames(dados_finais)[1:2] = c("ID.x", "ID.y")
  print(dados_finais)
  
  dados_finais = dados_finais[, - grep("ID", colnames(dados_finais))[-c(1:2)]]

} else print("todo: Parallel == FALSE")

for(j in 1 : total_studies)
{
  dados_finais = merge(dados_finais, t(dif_median_degs[[j]]), by.x = "ID.x", by.y = 0, all.x = T)
  colnames(dados_finais)[ncol(dados_finais)] = paste("dif_median.x", study.names[j], sep=".") 
  
  dados_finais = merge(dados_finais, t(dif_median_degs[[j]]), by.x = "ID.y", by.y = 0, all.x = T)
  colnames(dados_finais)[ncol(dados_finais)] = paste("dif_median.y", study.names[j], sep=".") 
  
  dados_finais = merge(dados_finais, t(regulation[[j]]), by.x = "ID.x", by.y = 0, all.x = T)
  colnames(dados_finais)[ncol(dados_finais)] = paste("regulation.x", study.names[j], sep=".") 
  
  dados_finais = merge(dados_finais, t(regulation[[j]]), by.x = "ID.y", by.y = 0, all.x = T)
  colnames(dados_finais)[ncol(dados_finais)] = paste("regulation.y", study.names[j], sep=".") 
  
  
  dados_finais[paste("interaction",study.names[j] , sep = ".")] = 
    paste(dados_finais[,paste("regulation.x", study.names[j], sep=".")], dados_finais[,paste("regulation.y", study.names[j], sep=".")], sep = "")
  
  colnames(dados_finais) = gsub("class1", paste("class",class_type[1],sep = "."), colnames(dados_finais))
  colnames(dados_finais) = gsub("class2", paste("class",class_type[2],sep = "."), colnames(dados_finais))
  
}

print(dados_finais)

setwd(path_degs)

if(Parallel == TRUE)
{
  if(length(study.names) == 1)
  {
    if(paste("degs_correlation_parallel_", study.names[1], method, ".txt", sep = "" ) %in% list.files(path_degs))
      write.table(dados_finais, paste("degs_correlation_parallel_", study.names[1], method, ".txt", sep = "" ), col.names = F , row.names = F, sep = "\t", append = T) else
        write.table(dados_finais, paste("degs_correlation_parallel_", study.names[1], method, ".txt", sep = "" ), col.names = T , row.names = F, sep = "\t", append = T)
  }else
  {
    study.names.pasted = study.names[1]
    
    for(j in 1 : total_studies)
      study.names.pasted = paste(study.names.pasted, study.names[j],sep = "")
    
    if(paste("degs_correlation_parallel_", study.names.pasted, "_", method, ".txt", sep = "") %in% list.files(path_degs))
      write.table(dados_finais, paste("degs_correlation_parallel_", study.names.pasted, "_", method, ".txt", sep = ""), col.names = F, row.names = F, sep = "\t" , append = T) else
        write.table(dados_finais, paste("degs_correlation_parallel_", study.names.pasted, "_", method, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t" , append = T)
  }
  
} else print("todo: write table Parallel == FALSE")


# if(pv_filter == TRUE)
# {
#   dados_pv_filter = lapply(dados_finais, pv.filter, n = total_studies, meta = TRUE, fdr = TRUE)
#   
#   for(j in 1:length(dados_pv_filter))
#     write.table(dados_pv_filter[[j]],paste("degs_correlation_pv_filter_02_",names(dados_pv_filter)[j],".txt" ,sep=""),col.names=T,row.names=F, sep="\t")
# }
# 
# rm(pv_filter)
 
