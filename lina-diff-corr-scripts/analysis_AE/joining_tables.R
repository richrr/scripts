#####################################################################################################################
# script written by Lina Thomas
# date: 12/11/2013

# run this script when output.option="one_table_each" in analysis_AE.R

######################################################################################################################

# command in cgrb
# 
# SGE_Batch -c 'R < /scripts/analysis_AE/joining_tables.R "cgrb" "/capecchi/pharmacy/morgunlab/lina/ArrayExpress" "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/DAPs_Dif" "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/DEGs_Dif" "all" "DAPs_Dif_table_first100.txt" NULL --no-save' -r joining_tables_first100daps -f 1G -m 1G 


#######################################################################################################################
#functions

reading_tables=function(daps_degs, path_daps, path_degs){
  
  if(daps_degs=="daps")
  {
    ls_linux= paste("ls -F ",path_daps," | grep 'DAPs_Dif_single'",sep="")
    
    Dif_DAPs=system(ls_linux, intern =TRUE)
    Dif_DEGs=NA
    
    Dif=list(DAPs=Dif_DAPs,DEGs=Dif_DEGs)
    
    return(Dif)
    
  } else if(daps_degs=="degs")
  {
    ls_linux= paste("ls -F ",path_degs," | grep 'DEGs_Dif_single'",sep="")
    
    Dif_DAPs=NA
    
    Dif_DEGs=system(ls_linux, intern =TRUE)
    
    Dif=list(DAPs=Dif_DAPs,DEGs=Dif_DEGs)
    
    return(Dif)
    
  } else if(daps_degs=="all")
  {
    ls_linux=paste("ls -F ",path_daps," | grep 'DAPs_Dif_single'",sep="")
    Dif_DAPs=system(ls_linux, intern =TRUE)
    
    ls_linux= paste("ls -F ",path_degs," | grep 'DEGs_Dif_single'",sep="")
    Dif_DEGs=system(ls_linux, intern =TRUE)
    
    Dif=list(DAPs=Dif_DAPs,DEGs=Dif_DEGs)
    
    return(Dif)
  } else
  {
    print("incorrect option: choose either daps or degs")
    stop
  }
      
		  
}
		
		#######################################################

merging_singles=function(daps_degs, data1, data2, file.out){
  temp=read.table(data1,header=T)
  if(daps_degs=="daps") data2=merge(data2,temp,by=c(1,2,3,4,5),all=TRUE) else
    if(daps_degs=="degs") data2=merge(data2,temp,by=c(1,2),all=TRUE)
  write.table(data2,file.out,col.names=TRUE,row.names=FALSE,sep="\t")
  
  return(data2)
  
}

		#########################################################

  joining_table=function(daps_degs,data,file.out,path){
    
    setwd(path)
    
    ls_table= paste("ls -F ",path," | grep '",file.out,"'",sep="")
    fdif=system(ls_table, intern =TRUE)
    
    if(length(fdif)==1)  fdif=read.table(fdif,header=T) else
      if(length(fdif)>1)  return(print("something strange in listing the file.out"))
  
    for(i in 1:length(data))
    {
      
      if(i==1)
      {
        if(length(fdif)==0) fdif=read.table(data[i],header=T) else
          fdif=merging_singles(daps_degs,data[i],fdif,file.out)       
                
      } else  fdif=merging_singles(daps_degs,data[i],fdif,file.out)
      
    }
    
  }

###########################################################################################
args=commandArgs(trailingOnly = FALSE)

local=args[2]

local
path=args[3]
path_daps=args[4]
path_degs=args[5]
daps_degs=args[6] # "daps", "degs" or "all"
daps.out=args[7]
degs.out=args[8]

print("did the arguments")
###########################################################################################

if(local=="lab")
  path="/home/thomasli/Dropbox/Tese/Dados/ArrayExpress" else
    if(local=="home")
      path="C:/Users/Lina/Dropbox/Tese/Dados/ArrayExpress"

          
Dif_single = reading_tables(daps_degs, path_daps, path_degs) 
print("read the name of the tables")

setwd(path)

if(!(1 %in% which(is.na(Dif_single)))) 
{
  joining_table("daps",Dif_single$DAPs,daps.out,path_daps)
  print("done for DAPs")
} else
 print("not done for DAPs")

if(!(2 %in% which(is.na(Dif_single)))) 
{
  joining_table("degs",Dif_single$DEGs,degs.out,path_degs)
  print("done for DEGs")
} else
  print("not done for DEGs")


