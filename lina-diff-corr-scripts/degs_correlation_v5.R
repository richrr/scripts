######################################################################################################

# Lina Thomas
# date: 12/13/13

# REQUIREMENT: The probe.ids must be in the same order in each study expression table!!!!!!!!!!!!!!!!
########################################################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < scripts/degs_correlation_v5.R --no-save 1 1 "pearson" TRUE FALSE "all"  ~/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/Network/Correlation "genes_for_correlation_network.txt" 3 /home/thomasli/Dropbox/Tese/scripts ~/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/ "dados1.txt" "dados2.txt" "dados3.txt" "dados1" "dados2" "dados3"  TRUE Identifiers.txt  ' -r degs_correlation_imsys -t 1 -f 8G -m 8G -q samwise
# test
# SGE_Batch -c 'R < scripts/degs_correlation_v5.R --no-save $SGE_TASK_ID 3 "pearson" TRUE FALSE "all"  /capecchi/pharmacy/morgunlab/lina/imsys/Network/Correlation "genes_for_correlation_network.txt" 1 ~/scripts /capecchi/pharmacy/morgunlab/lina/imsys /capecchi/pharmacy/morgunlab/lina/imsys/Network/Correlation "amostra1.txt" "dados1" TRUE Identifiers.txt  ' -r degs_daps_correlation_imsys -t 1-3 -f 1G -m 1G -q samwise
# R < /home/thomasli/Dropbox/Tese/scripts/Networks/correlation/degs_correlation_v5.R --no-save 1 1 "pearson" TRUE FALSE "all"  ~/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/Network/Correlation "genes_for_correlation_network.txt" 3 /home/thomasli/Dropbox/Tese/scripts  ~/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/ ~/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/Network/Correlation "amostra1.txt" "amostra2.txt" "amostra3.txt" "dados1" "dados2" "dados3"  TRUE Identifiers.txt 
# 
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

print(args)

pi = as.numeric(args[3]) #which part
p = as.numeric(args[4]) #how many parts
method = args[5]
# local = args[6]
Parallel = args[6] #running parallel or alltogether
pv_filter = args[7]

state = args[8] # case, control or all
if(state %in% c("case", "control"))
{
  print("todo: option of case and control separate")
  print("WARNING: considering state = all")
  state = "all"
}

if(state == "all") state = c("case", "control")

path_degs = args[9]
degs_orig = args[10] # file containing degs information

total_studies = as.numeric(args[11]) #number of datasets

# by.study=args[12]
path_functions = args[12]
if(nchar(path_functions) == "/") path_functions = substr(path_functions, 1, nchar(path_functions)-1)
path_data = args[13]
path_final = args[14]

study.file = vector()
for(j in 1:total_studies)
 study.file[j] = args[14 + j] #name of study file
counter1 = 14 + j
  
study.names = vector()
for(j in 1:total_studies)
  study.names[j] = args[counter1 + j] #name of study file
counter2 = counter1 + j

class_same_file = args[counter2 + 1] # TRUE or FALSE
counter2 = counter2 + 1

class.file = vector()
if(class_same_file == TRUE) class.file[1] = args[counter2+1] else
  for(j in 1:total_studies)
    class.file[j] = args[j + counter2] #name of study file
# counter3 = j + counter2

###############################################################################################################

#functions

source(paste(path_functions, "functions.R", sep = "/"))

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
  
  if(n==1) final_data=final_table(cor,pvalue,min.mat,gene.symbol=gene.symbol) else
    final_data=final_table(cor,pvalue,min.mat,n,gene.symbol)
  
  
  return(final_data)
}

#################################################################################################################
setwd(path_degs)
print(path_degs)
list.files(path_degs)
degs_orig

if(substr(degs_orig,nchar(degs_orig)-3,nchar(degs_orig)) == ".csv")
  degs=read.csv(degs_orig,header=TRUE) else
    if(substr(degs_orig,nchar(degs_orig)-3,nchar(degs_orig)) == ".txt")
      degs =read.table(degs_orig, header=TRUE) else
      {
        print("File type different from csv or txt. Please select either one.")
        stop()
      }

degs[,1]=paste("X",degs[,1],sep="")

print(head(degs))

setwd(path_data)

data=list()

for(j in 1:total_studies)
  data[[j]]=read.table(study.file[j],row.names=1,header=TRUE)
names(data) = study.names

class_type = data[[j]][!duplicated(data[[j]][,"class"]), "class"]
class_type

genes = degs[, "probe.id"]
pairs = t(combn(rev(genes),2))

if(p != 1)
{
  pair = subsetting(paste(pairs[, 1], pairs[, 2]), p, pi)
      
  pairs = pairs[which(paste(pairs[, 1], pairs[, 2]) %in% pair), c(1 : 2)]
  degs_selected = c(pairs[, 1], pairs[, 2])[!duplicated(c(pairs[, 1], pairs[, 2]))]
  print(pairs)
} else if(p == 1 & pi == 1)
  degs_selected = degs[, "probe.id"] else
  {
    
    print("p == 1, but pi greater than 1.")
    stop()
  }   
  
if(substr(degs_selected,1,1) != "X") 
  degs_selected = paste("X", degs_selected, sep = "")

# degs_selected
# gives a list with degs expression only where components are case and control
data_degs = lapply(data, function(x) x[,colnames(x) %in% degs_selected])

for(j in 1:total_studies)
  if("class" %in% colnames(data[[1]]))
    data_degs[[j]][, "class"] = data[[j]][, "class"] else 
    {
      if(class_same_file == TRUE)
      {
        if(substr(class.file, nchar(class.file) - 3, nchar(class.file)) == ".txt")
          class = read.table(class.file, header = T) else
            if(substr(class.file, nchar(class.file) - 3, nchar(class.file)) == ".csv")
              class = read.csv(class.file, header = T) else
              {
                print("ERROR: class file is not txt or csv")
                stop()
              }
        
        for(j in 1 : total_studies)
        {
          data_degs[[j]] = merge(data_degs[[j]], class[,c("probe.id","class")], by.x = 0, by.y = "probe.id", all.x = TRUE)
          row.names(data_degs[[j]]) = data_degs[[j]][,"Row.names"] 
        }
        
        data_degs = lapply(data_degs, function(x) x[, - which(colnames(x) == "Row.names")] )
        
      } else
      {
        class = vector()
        for(j in 1 : total_studies)
          if(substr(class.file[j], nchar(class.file[j]) - 3, nchar(class.file[j])) == ".txt")
            class[j] = read.table(class.file[j], header = T) else
              if(substr(class.file[j], nchar(class.file[j]) - 3, nchar(class.file[j])) == ".csv")
                class[j] = read.csv(class.file[j], header = T) else
                {
                  print("ERROR: class file is not txt or csv")
                  stop()
                }
        
        for(j in 1 : total_studies)
        {
          data_degs[[j]] = merge(data_degs[[j]], class[j][,c("probe.id","class")], by.x = 0, by.y = "probe.id", all.x = TRUE)
          row.names(data_degs[[j]]) = data_degs[[j]][,"Row.names"] 
        }
        
        data_degs = lapply(data_degs, function(x) x[, - which(colnames(x) == "Row.names")] )
        
        
      }
                  
    }

data_degs = separate.classes2(data_degs, study.names, class_type=c(0,1))

# print(data_degs[[1]][[2]])

dif_median_degs = dif.median.classes(data_degs)
dif_median_degs[[1]]

regulation = lapply(dif_median_degs, function(x) ifelse(x > 0 ,"U","D"))
head(regulation[[1]])

if(length(grep("amostra",study.file))!=0)
{
  genes = c(colnames(data_degs[[1]][[1]]))
  pairs = t(combn(rev(genes),2))
  if(p != 1)
  {
    pair = subsetting(paste(pairs[, 1], pairs[, 2]), p, pi)
    
    pairs = pairs[which(paste(pairs[, 1], pairs[, 2]) %in% pair), c(1 : 2)]
    degs_selected = c(pairs[, 1], pairs[, 2])[!duplicated(c(pairs[, 1], pairs[, 2]))]
    print(pairs)
  }
}

if(Parallel==TRUE)
{
  table_list = lapply(data_degs, function(x) lapply(x, table.perpair.cor.pv, pairs = pairs, method = method ))
  print(table_list)
  final_data = as.data.frame(table_list)
  colnames(final_data)[1:2] = c("ID.x", "ID.y")
  print(final_data)
  
  final_data = final_data[, - grep("ID", colnames(final_data))[-c(1:2)]]

} else 
  stop("to do")


for(j in 1 : total_studies)
{
  print(final_data)
  final_data = merge(final_data, t(dif_median_degs[[j]]), by.x = "ID.x", by.y = 0, all.x = T)
  colnames(final_data)[ncol(final_data)] = paste("dif_median.x", study.names[j], sep=".") 

  
  final_data = merge(final_data, t(dif_median_degs[[j]]), by.x = "ID.y", by.y = 0, all.x = T)
  colnames(final_data)[ncol(final_data)] = paste("dif_median.y", study.names[j], sep=".") 
  
  final_data = merge(final_data, t(regulation[[j]]), by.x = "ID.x", by.y = 0, all.x = T)
  colnames(final_data)[ncol(final_data)] = paste("regulation.x", study.names[j], sep=".") 
  
  final_data = merge(final_data, t(regulation[[j]]), by.x = "ID.y", by.y = 0, all.x = T)
  colnames(final_data)[ncol(final_data)] = paste("regulation.y", study.names[j], sep=".") 
  
  
  final_data[paste("interaction",study.names[j] , sep = ".")] = 
    paste(final_data[,paste("regulation.x", study.names[j], sep=".")], final_data[,paste("regulation.y", study.names[j], sep=".")], sep = "")
  
  colnames(final_data) = gsub("class1", paste("class",class_type[1],sep = "."), colnames(final_data))
  colnames(final_data) = gsub("class2", paste("class",class_type[2],sep = "."), colnames(final_data))
  
}

# print(final_data)

setwd(path_degs)

if(Parallel == TRUE)
{
  if(length(study.names) == 1)
  {
    if(paste("degs_daps_correlation_parallel_", study.names[1], method, ".txt", sep = "" ) %in% list.files(path_degs))
      write.table(final_data, paste("degs_daps_correlation_parallel_", study.names[1], method, ".txt", sep = "" ), col.names = F , row.names = F, sep = "\t", append = T) else
        write.table(final_data, paste("degs_daps_correlation_parallel_", study.names[1], method, ".txt", sep = "" ), col.names = T , row.names = F, sep = "\t", append = T)
  }else
  {
    study.names.pasted = study.names[1]
    
    for(j in 1 : total_studies)
      study.names.pasted = paste(study.names.pasted, study.names[j],sep = "")
    
    if(paste("degs_correlation_parallel_", study.names.pasted, "_", method, ".txt", sep = "") %in% list.files(path_degs))
      write.table(final_data, paste("degs_daps_correlation_parallel_", study.names.pasted, "_", method, ".txt", sep = ""), col.names = F, row.names = F, sep = "\t" , append = T) else
        write.table(final_data, paste("degs_daps_correlation_parallel_", study.names.pasted, "_", method, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t" , append = T)
  }
  
} else print("todo: write table Parallel == FALSE")


# if(pv_filter == TRUE)
# {
#   data_pv_filter = lapply(final_data, pv.filter, n = total_studies, meta = TRUE, fdr = TRUE)
#   
#   for(j in 1:length(data_pv_filter))
#     write.table(data_pv_filter[[j]],paste("degs_correlation_pv_filter_02_",names(data_pv_filter)[j],".txt" ,sep=""),col.names=T,row.names=F, sep="\t")
# }
# 
# rm(pv_filter)
 
