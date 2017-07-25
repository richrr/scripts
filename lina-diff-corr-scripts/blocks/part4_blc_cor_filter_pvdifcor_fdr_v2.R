######################################################################################################

# Lina Thomas
# date: 1/22/14

########################################################################################################################################
# 
####################################################################################################
# CGRB command:
# 
# SGE_Batch -c 'R < scripts/blocks/part4_blc_cor_filter_pvdifcor_fdr.R --no-save "cgrb" 3 "Project_03mp" NULL 0.2 0.2 "Identifiers.csv" "gene_identificators_probeid_genesymbol.csv" "dados1" "dados2" "dados3" "table_same_dir_difcor.txt"'  -r part4_imsys -f 256G -m 256G -q samwise
#
# tests
# SGE_Batch -c 'R < scripts/blocks/part4_blc_cor_filter_pvdifcor_fdr_v2.R --no-save 5 3 "Project_03mp_alias" NULL TRUE 0.2 0.2 "Mine_Class_all.txt" "Pyeon" "Scotto" "Zhai" "Biewenga" "Igen" Merged_difcor_table_all5_filtered_smp.txt /capecchi/pharmacy/morgunlab/lina/CervicalCancer /raid1/home/pharmacy/thomasli/scripts'  -r part4_mine_smp -f 2G -m 2G -q samwise
# at lab: R < Dropbox/Tese/scripts/DAPfinder/blocks/part4_blc_cor_filter_pvdifcor_fdr_v2.R 
######################################################################################################
# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

total_studies = as.numeric(args[3]) #number of datasets
min_studies = as.numeric(args[4])
newfoldername = args[5]
method = args[6]
difcor.filtered = args[7]
pv_thrsd_cor = args[8]
pv_thrsd_difcor = args[9]
class.file = args[10]
# gene.symbol.file = args[10]
 
study.names = vector()
for(counter1 in 1:total_studies)
  study.names[counter1] = args[counter1 + 10] #name of study file
counter1 = 10 + counter1

# gene.annotation.file = vector()
# for(counter2 in 1 : total_studies)
#   gene.annotation.file[counter2] = args[counter2 + counter1]
# counter2 = counter2 + counter1

table.file = args[counter1 + 1]

path_dados = args[counter1 + 2]
path_final = paste(path_dados, newfoldername ,sep="/")
path_function = args[counter1 + 3]

path_dados
path_function
# gene.annotation.file

##############################################################################################
# 
# if(!"psych" %in% installed.packages()) install.packages("psych")
# library("psych")

#################################################################################################################

#functions

source(paste(path_function, "functions.R", sep = "/"))
print("pegou as funcoes")


coherent.cor.directions.one.state=function(final_table,pv_thrsd_cor = 0.2, total_studies=1, study.names = "gr1" ,class_type=c(0,1), method = "pearson"){
  
  #checking if final_table already has cor directions assigned
  dir_assigned = grep("dircor",colnames(final_table))
  
  if(length(dir_assigned) == 0)
    final_table = cor.directions(final_table, pv_thrsd_cor, total_studies, study.names , class_type, method = method)
  
  # only keeping the pairs with the same correlation directions in all stydies for both classes

  dircor0_index = grep("dircor0", colnames(final_table))
  sum_dircor0 = apply(final_table[,dircor0_index],1,sum, na.rm = T)
  
  dircor1_index = grep("dircor1", colnames(final_table))
  
  sum_dircor1 = apply(final_table[,dircor1_index],1,sum, na.rm = T)
  
  sum.not.na = function(vector) {
    sum.not.na = sum(!is.na(vector))
    return(sum.not.na)
  }

  cond_sum0 = apply(final_table[,dircor0_index], 1, sum.not.na)
  cond_sum1 = apply(final_table[,dircor1_index], 1, sum.not.na)


  cond = abs(sum_dircor0) == cond_sum0 | abs(sum_dircor1) == cond_sum1
  final_table = final_table[cond,]
  
  return(final_table)
}


#################################################################################################
setwd(path_dados)
class.table = read.csv(class.file,header=TRUE, sep = "\t")

class_type = class.table[!duplicated(class.table[,"class"]),"class"] 

setwd(path_final)
final_table = read.table(table.file, header=T, stringsAsFactors=F, sep="\t")

if(difcor.filtered == FALSE)
{
  final_table = filter.dircor.v2(final_table,study.names[1:total_studies])
  final_table = final_table[final_table$number.no.na >= min_studies, ]
}

# head(final_table)

# we only want pairs where correlations are significant in at least one state
final_table = coherent.cor.directions.one.state(final_table, pv_thrsd_cor, total_studies, study.names, class_type)

final_table = pv.dif.cor.table(final_table, study.names, class_type)


### ATTENTION: PROBLEM SOLVED ONLY FOR IGEN. FIX FOR ALL STUDIES.
final_table$number.no.na = ifelse(is.na(final_table$difcor_Igen), final_table$number.no.na,
                                  ifelse(is.na(final_table$pv_difcor_Igen), final_table$number.no.na-1, final_table$number.no.na)                                  )

final_table = final_table[final_table$number.no.na >= 3,]

final_table = coherent.dircor.pvalue(final_table, pv_thrsd_difcor)

final_table[, "fisher_pv"] = pv.meta.analysis.v3( final_table[, grep("pv_difcor", colnames(final_table))] )

final_table[, "fisher_fdr"] = p.adjust(final_table[,"fisher_pv"], method = "fdr")

setwd(path_final)
write.table(final_table, paste("final_table_cor",gsub("[[:punct:]]","",as.character(pv_thrsd_cor)),"_","difcor_", gsub("[[:punct:]]","",as.character(pv_thrsd_difcor)),".txt",sep=""), row.names = FALSE, sep="\t")

#################################################################################################


