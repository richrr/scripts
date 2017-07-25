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
# SGE_Batch -c 'R < ~/scripts/blocks/part2_2_organizing_na.R --no-save ' -r cor_difcor_mine -f 128G -m 128G -q samwise
#
# tests
# SGE_Batch -c 'R < ~/scripts/blocks/part2_2_organizing_na.R --no-save 
# at lab: R < Dropbox/Tese/scripts/DAPfinder/blocks/part2_2_organizing_na.R --no-save "/home/thomasli/Dropbox/Tese/Undergraduate_Internship/Justin_Lina/Mine/tables/Project_03mp_teste table_same_dir_difcor_minstudies3.txt" 3 Pyeon Scotto Zhai table_same_dir_difcor_treate.txt
######################################################################################################

# functions

organizing.na = function(table2, study.names){
  for(j in 1:length(study.names))
  {
    row.na = which(is.na(table2[, paste( "cor0_", study.names[j], sep = "")]))
    num.na = length(row.na)
    if(num.na == 0 | num.na == nrow(table2)) next
    
    row.dif.cor = which(!duplicated(table2[, paste( "cor0_", study.names[j], sep = "")]) & !is.na(table2[, paste( "cor0_", study.names[j], sep = "")]))
    num.dif.cor = length (row.dif.cor)
    
    table2[row.na, grep(study.names[j], colnames(table2))] = table2[ row.dif.cor[1], grep(study.names[j], colnames(table2))]
    
    if(num.dif.cor > 1)
      for(i.dif in 2:num.dif.cor)
      {
        table2[ (nrow(table2)+1) : (nrow(table2) + num.na), ] = table2[row.na,]
        table2[ (nrow(table2)-num.na+1) : nrow(table2), grep(study.names[j], colnames(table2)) ] = table2[row.dif.cor[i.dif], grep(study.names[j], colnames(table2))]
      }
  }
  return(table2)
}

organizing.na.table = function(table, study.names)
{
  pairs = paste(table$Gene.symbol.x, table$Gene.symbol.y)
  pairs= pairs[!duplicated(pairs)]
  
  for(i in 1:length(pairs))
  {
    table2 = table[paste(table$Gene.symbol.x, table$Gene.symbol.y) == pairs[i],]
    num.row = nrow(table2)
    if(num.row == 1) next
    
    table2 = organizing.na(table2, study.names)
    
    table[row.names(table)[1:num.row],] = table2[1:num.row,]
    
    if(nrow(table2) > num.row)
    {
      table = list(table, table2[(num.row+1) : nrow(table2),])  
      table = do.call(rbind, table)
    }
    
  }
  return(table)
}


##########
## args from terminal

args = commandArgs(trailingOnly = FALSE)
path_table = args[3]
table.file = args[4]
probe.id.studies = as.numeric(args[5])

study.names = vector()
for(j in 1:probe.id.studies)
  study.names[j] = args[5+j]
counter = 5+j

out.file = args[counter+1]
############3
##main
setwd(path_table)
table = read.table(table.file, header = T)
study.names = c("Pyeon", "Scotto", "Zhai")

table = organizing.na.table(table, study.names)

write.table(table, out.file, row.names = F, sep = "\t")