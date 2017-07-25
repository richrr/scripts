#####################################################################################################################
# script written by Lina Thomas
# date: 12/10/2013

# IMPORTANT!!!!!!!!!!! HAS TO INSTALL CHIP PACKAGES ON CRICK BEFORE RUNNING ON CGRB!!!!!!!

# quando tem control em um transition state, nem sempre pega a classe equivalente a control, mas sim prioriza as 2
# classes com maior numero de samples

# when output.option="one_table_each", run joining_tables.R to get one general table after 
# all single tables are done!
######################################################################################################################

# command in cgrb
# 
# SGE_Batch -c 'R < scripts/analysis_AE/analysis_AE_v2.R $SGE_TASK_ID 1 1 1 "cgrb" 10 TRUE FALSE "one_table_each" "status.txt" "AEmusmusculus20smp_affy_chipinst.txt" TRUE "in" "done.txt" "alldegs_probe.id.txt" "final_table_daps_fdr001.csv" --no-save' -r daps_analysis_done -t 1 -f 8G -m 8G -q samwise


#######################################################################################################################
args=commandArgs(trailingOnly = FALSE)

pi=as.numeric(args[2]) #which part
p=as.numeric(args[3]) #how many parts we are going to devide studies
start_studies=as.numeric(args[4])
num_studies=args[5] # either a number or "total"
if(num_studies!="total") num_studies=as.numeric(num_studies)
if(is.na(num_studies))
{
  print("ERROR: num_studies must be either a number or total")
  stop
}

num_studies

local= args[6]      #three possible options: "cgrb", "lab", "home"
min_n=as.numeric(args[7]) #minimum number of samples in each state
do.degs= args[8]
do.daps= args[9]
output.option=args[10]	#two possible options: "one_table_each" and "single table"
status_data= args[11]   # name of the file to save the status for each study
run_data=args[12]     # name of the file with the studies to run
rerun=args[13]       #two possible options: "yes" or "no"
select.type=args[14]	# two possible options: "in" or "out" of the rerun_data file
rerun_data=args[15]	# name of the file with the studies to filter/rerun

degs.file = args[16]
daps.file = args[17]

# path_daps # path of the 1 column DAP txt file 


#########################################################################################
if(local=="home")
{
  path="C:/Users/Lina/Dropbox/Tese/Dados/ArrayExpress"
  path_daps = "C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder/Project_03mp"
  path_degs = "C:/Users/Lina/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DEGs"
  source("C:/Users/Lina/Dropbox/Tese/scripts/100_sets/functions_daps_analysis.R")
  
}else if(local=="lab")
{
  path="/home/thomasli/Dropbox/Tese/Dados/ArrayExpress/musmusculus"
  path_daps = "/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DAPfinder/Project_03mp"
  path_degs = "/home/thomasli/Dropbox/Tese/Lina-Andrey/raw-data-BcKO/DEGs"
  source("/home/thomasli/Dropbox/Tese/scripts/100_sets/functions_daps_analysis.R")
} else  if(local=="cgrb")
{
  path = "/capecchi/pharmacy/morgunlab/lina/ArrayExpress"
  path_studies="/capecchi/pharmacy/morgunlab/lina/ArrayExpress/RData/MusMusculus_done" 
  path_daps = "/capecchi/pharmacy/morgunlab/lina/imsys/Project_03mp"
  path_degs = "/capecchi/pharmacy/morgunlab/lina/imsys/degs"
    
  source("/raid1/home/pharmacy/thomasli/scripts/analysis_AE/functions_daps_analysis.R")
  
#   ordpair=function(A){
#     A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
#     genecol1=unique(c(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]))
#     genecol1=sort(genecol1)
#     
#     for(k in 1:length(genecol1))
#     {
#       
#       if(k==1)  where=which(A[,"Gene.symbol.y"]==genecol1[k]) else
#         where=which(A[,"Gene.symbol.y"]==genecol1[k] & !(A[,"Gene.symbol.x"]%in% genecol1[1:k]))
#       
#       who=A[where,"Gene.symbol.x"]
#       A[where,"Gene.symbol.x"]=c(rep(genecol1[k],length(who)))
#       A[where,"Gene.symbol.y"]=who
#       A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
#       
#     } 
#     
#     return (A)
#   }
#   
#   ##########################################################################################
#   # function that divides the vector in the pith part
#   # studies : vector
#   # p : how many parts to divide
#   # pi: the pith part
#   
#   subsetting=function(studies,p,pi){
#     
#     if(!is.vector(studies)) return("object is nor a vector!!")
#     
#     if(p>length(studies)) return("p is greater than length of studies!!")
#     
#     if(pi>p) return("impossible to get a pith part greater than the total number of parts!!")
#     
#     #divisao inteira
#     n=length(studies)%/%p
#     # #resto da divisao
#     f=(length(studies)%%p)/p
#     
#     if (pi <= f*p) studies=studies[((pi-1)*n+pi):(pi*(n+1))] else
#       studies=studies[((pi-1)*n+1+f*p):(pi*n+f*p)]
#     
#     return(studies)
#     
#   }
#   ###########################################################################################
#   
#   loading_AEdata=function(study){
#     
#     files=list.files(paste(path,study,sep="/"))
#     #fileslist=as.list(c(paste(path,studies[i],sep="/"),files))
#     
#     raw=grep("raw" ,files,value=T)
#     
#     if(length(raw)==0)
#     {
#       write.table(t(c(studies[i],1,"no raw data1")),status_data,append=T,row.names=F,col.names=F,sep="\t")
#       return("no raw data")
#     }
#     
#     raw=grep("zip",raw,value=T)
#     
#     if(length(raw)==0)
#     {
#       write.table(t(c(studies,1,"no raw data2")),status_data,append=T,row.names=F,col.names=F,sep="\t")
#       return("no rawdata")
#       
#     } else if(length(raw)==1)
#       rawfiles=data.frame(unzip(paste(path,studies,raw,sep="/"),list=T)[,1],stringsAsFactors = FALSE) else
#       {
#         print("tem mais de um estudo")
#         rawfiles=unzip(paste(path,studies,raw[1],sep="/"),list=T)
#         rawfiles <- data.frame(lapply(rawfiles, as.character), stringsAsFactors=FALSE) 
#         for(raw_id in 2:length(raw))
#         {
#           temp= unzip(paste(path,studies,raw[raw_id],sep="/"),list=T)
#           temp <- data.frame(lapply(temp, as.character), stringsAsFactors=FALSE) 
#           rawfiles=data.frame(c(rawfiles[,1],temp[,1]),stringsAsFactors = FALSE)
#         }
#       }
#  
#     
#     AEfiles=list(sdrf=grep("sdrf",files,value=T),adf=grep("adf",files,value=T),idf=grep("idf",files,value=T),
#                  path=paste(path,studies,sep="/"),rawFiles=rawfiles[,1])
#     
#     
#     AEset=ae2bioc(mageFiles = AEfiles)
#     
#     if(class(AEset)!="list") AEset=list(AEset)
#     
#     return(AEset)
#     
#   }
#   
#   ######################################################################################
#   loading.studies.chip3=function(run_data, p, pi, rerun, select.type, rerun_data){
#     
#     studies_chip = read.table(run_data, header = F)
#     
#     if(rerun == TRUE)
#     {
#       rerun2=read.table(rerun_data,header=F,stringsAsFactors=F)[,1]
#       if(select.type=="out") studies_chip=studies_chip[!(studies_chip[,1] %in% rerun2),] else
#         if(select.type=="in") studies_chip=studies_chip[studies_chip[,1] %in% rerun2,] else
#         {
#           print("choose in ou out for select.type")
#           stop 
#         }
#     }
#     
#     chips=unique(studies_chip[,3])
#     chips <- as.character(sapply(chips, as.character), stringsAsFactors=FALSE)
#     
#     for(i in 1:length(chips))
#       library(paste(chips[i],".db",sep=""),character.only=T)
#     
#     studies = unique(studies_chip[,1])
#     studies = as.character(sapply(studies, as.character), stringsAsFactors=FALSE)
#     
#     studies = subsetting(studies,p,pi)
#     
#     return(studies)
#     
#   }
  
}
          

setwd(path_studies)


####################################################################################################################
# loading functions from other file



######################################################################################################################

#the fist time you run this code start the following table. The idea is to append
if(local=="lab")
  if(length(system("ls -F /home/thomasli/Dropbox/Tese/Dados/ArrayExpress/musmusculus | grep 'status.txt'", intern =TRUE))==0)  
    write.table(t(c("Array.Express.Data","group","reason")),status_data, row.names=F, col.names=F,sep="\t") else
      print("Warning: status File already exists, appending")

if(local=="home")
  if(length(system("ls -F C:/Users/Lina/Dropbox/Tese/Dados/ArrayExpress | grep 'status.txt'", intern =TRUE))==0)  
    write.table(t(c("Array.Express.Data","group","reason")),status_data, row.names=F, col.names=F,sep="\t") else
      print("Warning: status File already exists, appending")

if(local=="cgrb")
  if(length(system("ls -F /capecchi/pharmacy/morgunlab/lina/ArrayExpress | grep 'status.txt'", intern =TRUE))==0)  
    write.table(t(c("Array.Express.Data","group","reason")),status_data, row.names=F, col.names=F,sep="\t") else
      print("Warning: status File already exists, appending")


##################################################################################################################

setwd(path)

#######################################################################
#install packages in crick!!!

if (! "annotate" %in% installed.packages()) 
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("annotate")
}

library("annotate")

if (! "ArrayExpress" %in% installed.packages()) 
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("ArrayExpress")
}
library("ArrayExpress")

if (! "affy" %in% installed.packages()) 
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("affy")
}
library("affy")

if(local=="cgrb")
{
  print(path)
  setwd(path)
  studies=loading.studies.chip3(run_data,p,pi,rerun,select.type,rerun_data)
  
}else
  if(local=="lab")
  {
    #reading all the folders from file
    studies=system("ls -F /home/thomasli/Dropbox/Tese/Dados/ArrayExpress | grep '/'", intern =TRUE)
    #getting rid of the character /
    studies=gsub("/","",studies)
    studies=grep("E-",studies,value=T)
  }else if(local=="home")
  {
    studies=grep("E-",list.dirs(path, full.names=FALSE,recursive=FALSE),value=T)
    studies=substr(studies,(nchar(path)+2),nchar(studies))
  }
    
#####################################################################################################################
fez=data.frame("Array Express Data"=NA,"group"=NA)
f=1

##########################################################################################

if(num_studies=="total") num_studies=length(studies)
num_studies

# max(sapply(AEfiles,length))
for(i in start_studies:num_studies)
{
  anchor=i

  print(paste("i = ",i))
  print(studies[i])
  
  RData = grep(paste(studies[i],"_",1,".RData",sep=""),list.files(path_studies),value=TRUE)
  print(RData)
  
  if(length(RData) == 1)
  {
    if(RData==paste(studies[i],"_",1,".RData",sep=""))
    {
      warning("ALREADY HAD RData") 
      setwd(path_studies)
      load(paste(studies[i],"_",1,".RData",sep=""))

      i=anchor
	
      args=commandArgs(trailingOnly = FALSE)
      
      args
      
      pi=as.numeric(args[2]) #which part
      p=as.numeric(args[3]) #how many parts we are going to devide studies
      start_studies=as.numeric(args[4])
      num_studies=args[5] # either a number or "total"
      if(num_studies!="total") num_studies=as.numeric(num_studies)
      if(is.na(num_studies))
      {
        print("ERROR: num_studies must be either a number or total")
        stop
      }
      local= args[6]      #three possible options: "cgrb", "lab", "home"
      min_n=as.numeric(args[7]) #minimum number of samples in each state
      do.degs= args[8]
      do.daps= args[9]
      output.option=args[10]  #two possible options: "one_table_each" and "single table"
      status_data= args[11]   # name of the file to save the status for each study
      run_data=args[12]     # name of the file with the studies to run
      rerun=args[13]       #two possible options: "yes" or "no"
      select.type=args[14]  # two possible options: "in" or "out" of the rerun_data file
      rerun_data=args[15]	# name of the file with the studies to filter/rerun
      degs.file = args[16]
      daps.file = args[17]
      
      if(local=="cgrb")
      {
        path = "/capecchi/pharmacy/morgunlab/lina/ArrayExpress"
        path_studies="/capecchi/pharmacy/morgunlab/lina/ArrayExpress/RData/MusMusculus_done" 
        path_daps = "/capecchi/pharmacy/morgunlab/lina/imsys/Project_03mp"
        path_degs = "/capecchi/pharmacy/morgunlab/lina/imsys/degs"
        ##function which orders the pair of genes. The matrix A must have two columns named "gene.symbol 1" and "gene.symbol 2"
        
        source("/raid1/home/pharmacy/thomasli/scripts/analysis_AE/functions_daps_analysis.R")
      }
      
      setwd(path)
    	studies=loading.studies.chip3(run_data,p,pi,rerun,select.type,rerun_data)
      
    	RData=grep(paste(studies[i],"_",1,".RData",sep=""),list.files(path_studies),value=TRUE)
      
    } else if(RData==paste(studies[i],"_",1,".RDataTmp",sep=""))
    {
      print(RData)
      file.remove(paste(studies[i],"_",1,".RDataTmp",sep=""))
      RData=grep(paste(studies[i],"_",1,".RData",sep=""),list.files(path_studies),value=TRUE)
    }
      
  }
  
  if(length(RData)==0)
  {
    print("there is no RData yet")
    
    setwd(path)
    AEset=loading_AEdata(studies[i])
    
    if(is.character(AEset))
    {
      print(AEset)
      next
    }
    
    for(j in 1:length(AEset))
    {
      print(paste("j =",j))
            
      if(class(AEset[[j]])!="AffyBatch")
      {
        write.table(t(c(studies[i],j,"it is not Affymatrix data")),status_data,append=T,row.names=F,col.names=F,sep="\t")
        print("it is not Affymatrix data")
        next
      } 
    
      #rma is only for affymatrix data
      AEsetnorm = rma(AEset[[j]])
 
      if(exists("anchor")) rm(anchor)
      
      setwd(path_studies)
      save.image(file=paste(studies[i],"_",j,".RData",sep=""))
      anchor=i
      
      RDataj=grep(paste(studies[i],"_",j,".RData",sep=""),list.files(path_studies),value=TRUE)
      if(length(RDataj)==1)
        if(RDataj==paste(studies[i],"_",j,".RDataTmp",sep=""))
        {
          file.remove(RDataj)
          print("problem with RDat")
          write.table(t(c(studies[i],j,"problem with RData")),status_data,append=T,row.names=F,col.names=F,sep="\t")
          next
        }
    }
    
    RData=grep(paste(studies[i],"_",1,".RData",sep=""),list.files(path),value=TRUE)
    
    load(paste(studies[i],"_",1,".RData",sep=""))
    
    i=anchor
    
    args=commandArgs(trailingOnly = FALSE)
    
    pi=as.numeric(args[2]) #which part
    p=as.numeric(args[3]) #how many parts we are going to devide studies
    start_studies=as.numeric(args[4])
    num_studies=args[5] # either a number or "total"
    args=commandArgs(trailingOnly = FALSE)
    
    args
    
    pi=as.numeric(args[2]) #which part
    p=as.numeric(args[3]) #how many parts we are going to devide studies
    start_studies=as.numeric(args[4])
    num_studies=args[5] # either a number or "total"
    if(num_studies!="total") num_studies=as.numeric(num_studies)
    if(is.na(num_studies))
    {
      print("ERROR: num_studies must be either a number or total")
      stop
    }
    
    local= args[6]      #three possible options: "cgrb", "lab", "home"
    min_n=as.numeric(args[7]) #minimum number of samples in each state
    do.degs= args[8]
    do.daps= args[9]
    output.option=args[10]  #two possible options: "one_table_each" and "single table"
    status_data= args[11]   # name of the file to save the status for each study
    run_data=args[12]     # name of the file with the studies to run
    rerun=args[13]       #two possible options: "yes" or "no"
    select.type=args[14]  # two possible options: "in" or "out" of the rerun_data file
    rerun_data=args[15]	# name of the file with the studies to filter/rerun
    degs.file = args[16]
    daps.file = args[17]
    
    if(local=="cgrb")
    {
      path = "/capecchi/pharmacy/morgunlab/lina/ArrayExpress"
      path_studies="/capecchi/pharmacy/morgunlab/lina/ArrayExpress/RData/MusMusculus_done" 
      path_daps = "/capecchi/pharmacy/morgunlab/lina/imsys/Project_03mp"
      path_degs = "/capecchi/pharmacy/morgunlab/lina/imsys/degs"
      ##function which orders the pair of genes. The matrix A must have two columns named "gene.symbol 1" and "gene.symbol 2"
      
      source("/raid1/home/pharmacy/thomasli/scripts/analysis_AE/functions_daps_analysis.R")
    }
    
    print(path)
    setwd(path)
    studies=loading.studies.chip3(run_data,p,pi,rerun,select.type,rerun_data)
    
    RData=grep(paste(studies[i],"_",1,".RData",sep=""),list.files(path),value=TRUE)
    
    
  } else  if(length(RData)>1) print("Warning: More than one RData files")
  
  
  for(j in 1:length(AEset))
  {  
    print(paste("j = ",j))
    
    if(j==1 & length(RData)==0)
    {
      print("RData 1 does not exist")
      write.table(t(c(studies[i],j,"RData 1 does not exist")),status_data,append=T,row.names=F,col.names=F,sep="\t")
      next
    }
      
    
    if(j>1)
    {
      RData=grep(paste(studies[i],"_",j,".RData",sep=""),list.files(path_studies),value=TRUE)
      
      if(length(RData)==0)
      {
        print(paste("RData", j, "does not exist"))
        write.table(t(c(studies[i],j,paste("RData", j, "does not exist"))),status_data,append=T,row.names=F,col.names=F,sep="\t")
        next
      } else if(length(RData)==1)
      {
        if(RData!=paste(studies[i],"_",j,".RData",sep="")) 
        {
          if(!exists("AEset")) 
          {
            AEset=loading_AEdata(studies[i])
            if(is.character(AEset))
            {
              print(AEset)
              next
            }
          }
          AEsetnorm = rma(AEset[[j]])
          
        }else 
        {
          setwd(path_studies)
          load(paste(studies[i],"_",j,".RData",sep=""))


	        i=anchor	
	
	        args=commandArgs(trailingOnly = FALSE)
                     
          pi=as.numeric(args[2]) #which part
          p=as.numeric(args[3]) #how many parts we are going to devide studies
          start_studies=as.numeric(args[4])
          num_studies=args[5] # either a number or "total"
          
          if(num_studies!="total") num_studies=as.numeric(num_studies)
          if(is.na(num_studies))
          {
            print("ERROR: num_studies must be either a number or total")
            stop
          }
        
          local= args[6]      #three possible options: "cgrb", "lab", "home"
          min_n=as.numeric(args[7]) #minimum number of samples in each state
          do.degs= args[8]
          do.daps= args[9]
          output.option=args[10]  #two possible options: "one_table_each" and "single table"
          status_data= args[11]   # name of the file to save the status for each study
          run_data=args[12]     # name of the file with the studies to run
          rerun=args[13]       #two possible options: "yes" or "no"
          select.type=args[14]  # two possible options: "in" or "out" of the rerun_data file
          rerun_data=args[15]	# name of the file with the studies to filter/rerun
          degs.file = args[16]
          daps.file = args[17]
          
          if(local=="cgrb")
          {
            path = "/capecchi/pharmacy/morgunlab/lina/ArrayExpress"
            path_studies="/capecchi/pharmacy/morgunlab/lina/ArrayExpress/RData/MusMusculus_done" 
            path_daps = "/capecchi/pharmacy/morgunlab/lina/imsys/Project_03mp"
            path_degs = "/capecchi/pharmacy/morgunlab/lina/imsys/degs"
            ##function which orders the pair of genes. The matrix A must have two columns named "gene.symbol 1" and "gene.symbol 2"
            
            source("/raid1/home/pharmacy/thomasli/scripts/analysis_AE/functions_daps_analysis.R")
          }
          
          print(path)
          setwd(path)
          studies=loading.studies.chip3(run_data,p,pi,rerun,select.type,rerun_data) 
          
        }
        
      }else
      {
        if(!exists("AEset")) 
        {
          AEset=loading_AEdata(studies[i])
          if(is.character(AEset))
          {
            print(AEset)
            next
          }
        }
        AEsetnorm = rma(AEset[[j]])
      }
      
    }
      
                  
    
    #getting the possible states
    fac = grep("Factor.Value",colnames(pData(AEsetnorm)), value=T)
    if(length(fac)==0) fac=grep("FactorValue",colnames(pData(AEsetnorm)), value=T)
    if(length(fac)==0) fac=grep("Description",colnames(pData(AEsetnorm)), value=T)
    if(length(fac)==0)
    {
      print("no classification of states described as factor value")
      write.table(t(c(studies[i],j,"There is no classification of states described as factor value")),status_data,append=T,row.names=F,col.names=F,sep="\t")
      next
    }
    
    #getting classification of all samples into the each possible state transition 
    facs = pData(AEsetnorm)[,fac]

    
    if(length(fac)==1) 
    {
      num_smp=length(facs)
      if(length(table(facs))==1)
      {	
        print("has only one state")
        write.table(t(c(studies[i],j,"There is only on state")),status_data,append=T,row.names=F,col.names=F,sep="\t")
        next

      } 
    } else num_smp=nrow(facs)
    num_smp

    if(num_smp < 2*min_n)
    {
      print("not enough samples1")
      
      setwd(path)
      write.table(t(c(studies[i],j,paste("not enough samples (",min_n," min) remove study" ,sep=""))),status_data,append=T,row.names=F,col.names=F,sep="\t") 
        
        
      next
    }


    
    if(length(fac)==1 & length(names(facs))==0)
      names(facs)=pData(AEsetnorm)[,"Array.Data.File"]
    
    #getting names of possible classification for each state transition
    if(length(fac) != 1) labels=apply(facs,FUN=unique,2) else
      labels=fac
    
    print("fez labels")
  
    
    if(!is.matrix(labels)) nspec=grep("not specified",labels) else
    {
      nspec=vector()
      for(index_nspec in 1:ncol(labels))
        nspec[index_nspec]=grep("not specified",labels[,index_nspec])
      
      nspec=which(!is.na(nspec))
    }
    
    nspec
    
    if (length(nspec)== length(labels))
    {
      print("states not specified")
      write.table(t(c(studies[i],j,"states not specified")),"status.txt",append=T,row.names=F,col.names=F,sep="\t")
      next
    }else if(length(nspec)!= 0) 
    {
      print("tem not specified")
      nspec=grep("not specified",labels)
      facs = facs[,-nspec]
      fac=fac[-nspec]
      
      #update labels
      if(length(fac) != 1) labels=apply(facs,FUN=unique,2) else
        labels=fac
    }
    
    
         
    if(length(facs) != 0)
    {
      if(exists("freq_states")) rm(freq_states)
      
      print("tem classificacao de estados")
      
          
      if(!is.matrix(labels)) control=c(grep(c("control"),labels),grep(c("Control"),labels)) else
	    {
	      control=data.frame()
      
	      for(control_index in 1:length(fac))
	      {
	        control_status=c(grep(c("control"),labels[,control_index]),grep(c("Control"),labels[,control_index]))
	        if(length(control_status)!=0) control[control_index,1]=control_status else
	          control[control_index,1]=NA	      
	      }
	        
	      control=which(!is.na(control[,1]))

	    }
      print(control)
      
      min_num_states=9999999
      
      print("fez o control")
    
      if(length(control)!=0)
      {
        print("e tem control!")
        
        print(class(facs))
        
        if(!is.vector(facs)) facs2=facs[,control] else
          facs2=facs[control]
        
        fac2=fac[control]
        
        print(facs2)

      	if(length(fac2)==1)
      	{
	        if(length(table(facs2))==1)
	        {	
            print("has only one state")
            write.table(t(c(studies[i],j,"There is only on state")),status_data,append=T,row.names=F,col.names=F,sep="\t")
            next
          }
      	} 	


        print("filtrou so o control")
        
        print(class(facs2))
        print(facs2)
        print(colnames(facs2))
        print(names(facs2))
             
        if(length(control)==1)
          min_num_states=min(table(facs2)) else
          {
            print("mais de um control")
            
            freq_states=list()
            
            for(ind_freq in 1:length(control))
              freq_states[[fac2[ind_freq] ]]=table(facs2[,ind_freq])
            
            print(freq_states)
            
            if(length(freq_states)==0) 
            {
              print("NA states")
              write.table(t(c(studies[i],j,"NA states")),status_data,append=T,row.names=F,col.names=F,sep="\t")
              next
            }
            
            freq_states=lapply(freq_states,sort,decreasing=T)
            
#             for(ind_freq in 1:length(freq_states))
#               freq_states[[ind_freq]]=freq_states[[ind_freq]][-(names(freq_states[[ind_freq]])=="  ")]
            
            for(ind_freq in 1:length(freq_states))
            {
              if(sum(names(freq_states[[ind_freq]])=="  ")!=0)
                freq_states[[ind_freq]]=freq_states[[ind_freq]][-(names(freq_states[[ind_freq]])=="  ")]
              
              if(sum(names(freq_states[[ind_freq]])=="n/a")!=0)
                freq_states[[ind_freq]]=freq_states[[ind_freq]][-(names(freq_states[[ind_freq]])=="n/a")]
              
              if(sum(names(freq_states[[ind_freq]])=="not specified")!=0)
                freq_states[[ind_freq]]=freq_states[[ind_freq]][-(names(freq_states[[ind_freq]])=="not specified")]
              
            }
            
            print(freq_states)
            
                      
            min_num_states=sapply(freq_states, function(x) x[2])
            
            print(min_num_states)
          }
        
      }
      
      if(length(control)==0 | max(min_num_states)<min_n)
      {
        if(length(control)==0) print("mas nao tem control")
        
        if(max(min_num_states)<min_n) print("tem control mas nao o numero suficiente de samples")
        
        facs2=facs
        fac2=fac
       
        
        #checking how many samples we have in each state per type of state transition
        if(exists("freq_states")) rm(freq_states)
        
        if(length(fac2)==1)
          min_num_states=min(table(facs2)) else
          {
            print("tem mais de um tipo de transicao")
            freq_states=apply(facs2, FUN=table, 2)
            
            if(length(freq_states)==0) 
            {
              print("NA states")
              write.table(t(c(studies[i],j,"NA states")),status_data,append=T,row.names=F,col.names=F,sep="\t")
              next
            }
              
            if(!(is.list(freq_states))) freq_states=split(freq_states, rep(1:ncol(freq_states), each = nrow(freq_states)))
            freq_states=lapply(freq_states,sort,decreasing=T)
            
            for(ind_freq in 1:length(freq_states))
            {
              if(sum(names(freq_states[[ind_freq]])=="  ")!=0)
                freq_states[[ind_freq]]=freq_states[[ind_freq]][-(names(freq_states[[ind_freq]])=="  ")]
              
              if(sum(names(freq_states[[ind_freq]])=="n/a")!=0)
                freq_states[[ind_freq]]=freq_states[[ind_freq]][-(names(freq_states[[ind_freq]])=="n/a")]
              
              if(sum(names(freq_states[[ind_freq]])=="not specified")!=0)
                freq_states[[ind_freq]]=freq_states[[ind_freq]][-(names(freq_states[[ind_freq]])=="not specified")]
              
            }
                            
            
            min_num_states=sapply(freq_states, function(x) x[2])
          } 
       }
      
      print("fez as frequencias")
      #if non state transition has at least 2 states with at least min_n samples, stop
      if(max(min_num_states,na.rm=T)<min_n) 
      {
        print("not enough samples2")
        
        setwd(path)
        if(length(grep("Factor.Value",colnames(pData(AEsetnorm))))==0 & length(grep("FactorValue",colnames(pData(AEsetnorm))))==0)
          write.table(t(c(studies[i],j,paste("not enough samples (",min_n," min) e states not classified as factor",sep=""))),status_data,append=T,row.names=F,col.names=F,sep="\t") else
            {write.table(t(c(studies[i],j,paste("not enough samples (",min_n," min)",sep=""))),status_data,append=T,row.names=F,col.names=F,sep="\t")}
        
        next
        
      }
      
      
      #getting state transition with more number of samples considering only two states per transition
      if(exists("freq_states"))
        max_smp_index=grep(max(sapply(freq_states, function(x) x[1]+x[2]),na.rm=T),sapply(freq_states, function(x) x[1]+x[2])) else
          max_smp_index=1
      
      
      if(length(max_smp_index)>1) 
        states_index=max_smp_index[grep(max(min_num_states[max_smp_index]),min_num_states[max_smp_index])] else
          states_index=max_smp_index
      
      print("fez o state index")
      
      if(!is.data.frame(facs2)) facs2=as.data.frame(facs2, stringsAsFactors=F) 
      if (length(states_index)>1) states_index=states_index[1]
      #facs[states_index]
      state1=unique(facs2[states_index])[1,1]
      state2=unique(facs2[states_index])[2,1]
      
      print("separou os estados")
      #getting gene.symbol
      dados=as.data.frame(exprs(AEsetnorm))
      chip=paste(AEsetnorm@annotation,".db",sep="")
      
      # IMPORTANT!!!!!!!!!!! HAS TO INSTALL PACKAGES ON CRICK BEFORE RUNNING ON CGRB!!!!!!!
      
#       if (!chip %in% installed.packages()) 
#       {
#         source("http://bioconductor.org/biocLite.R")
#         biocLite(chip)
#       }
#       
#       library(chip,character.only=T)
           
      if(chip %in% installed.packages()) library(chip,character.only=T) else
      {
        print("pacote nao instalado")
        write.table(t(c(studies[i],j,paste("pacote ",chip," not installed"))),status_data,append=T,row.names=F,col.names=F,sep="\t")
        next
      }
      
      
      dados$Gene.symbol=getSYMBOL(row.names(dados),chip)
      
      print("pegou gene symbol")
      
      if(is.character(facs2))
      {
        if(length(grep(".CEL",names(facs2)))==0) samples=names(facs2) else
          samples=substr(names(facs2),1,(nchar(names(facs2))-4))
      }else
        if(length(grep(".CEL",row.names(facs2)))==0) samples=row.names(facs2) else
          samples=substr(row.names(facs2),1,(nchar(row.names(facs2))-4))
 
      
      if(do.degs==TRUE)
      {
        # reading DEGs
        
        setwd(path_degs)
        
        degs = read.table(degs.file, header=TRUE, stringsAsFactors=F)
                
        dados_degs=dados[dados$Gene.symbol %in% degs[,"Gene.symbol"],]
        #ident=data.frame("ID"=row.names(dados_degs),"Gene.symbol"=dados_degs[,"Gene.symbol"],row.names=NULL)
        ident=data.frame("Gene.symbol"=dados_degs[,"Gene.symbol"],row.names=row.names(dados_degs))
        dados_degs=data.matrix(dados_degs[,-ncol(dados_degs)])
        
        dl=list()
        if(sum(colnames(dados_degs)==samples)==(ncol(dados_degs)))
        {
          dl[[1]]= t(dados_degs[,t(facs2[states_index]) %in% state1])
          dl[[2]]= t(dados_degs[,t(facs2[states_index]) %in% state2])
        } else
        {
          print("problem with sample names in degs")
          write.table(t(c(studies[i],j,"problem with sample names in degs")),status_data,append=T,row.names=F,col.names=F,sep="\t")
          next
        }
        print("separou as classes in degs")
        
        med=lapply(dl,function(x) apply(x,2,median,na.rm=TRUE))
        print("fez a mediana")
        
        dif_degs=as.data.frame(med[[1]]-med[[2]])
        dif_degs=merge(ident,dif_degs,by=0)
        
#         dif_degs <- dif_degs[order(dif_degs[,2],- abs(dif_degs[,3]) ), ] #sort by id and abs(value)
#         dif_degs= dif_degs[ !duplicated(dif_degs[,2]), ][,-1]   
        
        print("fez a diferenca das medianas")
        
        if(nrow(dif_degs)!=0) 
        {
          fez[f,c("Array.Express.Data","group")]=c(studies[i],j)
          write.table(t(c(studies[i],j,"DEGs done")),status_data,append=T,row.names=F,col.names=F,sep="\t")
        }else
        {
          write.table(t(c(studies[i],j,"no DEGs in this AE dataset")),status_data,append=T,row.names=F,col.names=F,sep="\t")
          next
        }
        
        if(output.option=="single table")
        {
          if(nrow(fez)==1)
          {
            fdif_deg=as.data.frame(dif_degs)
            colnames(fdif_deg)[3]=studies[i]
                       
          } else
          {  
            
            fdif_deg=merge(fdif_deg,dif_degs,by=1,all=TRUE)
            colnames(fdif_deg)[ncol(fdif_deg)]=studies[i]
            write.table(fdif_deg,file="dif_degs.txt",col.names=TRUE,row.names=F,sep="\t")
            
          }
        } else if(output.option=="one_table_each")
        {
	  print(getew())
          fdif_deg=data.frame(dif_degs)
          colnames(fdif_deg)[3]=paste(studies[i],j,sep=".")
          write.table(fdif_deg,paste("DEGs_Dif_single_",studies[i],"_",j,".txt",sep=""),row.names=F, col.names=T,sep="\t")
          print("escreveu DEGs single table")
        }
        
      }
      
      if(do.daps==TRUE)
      {
        # reading DAPs
        
        print(path_daps)
        setwd(path_daps)
        
        daps = read.table(file = daps.file, header = TRUE, sep = ",", stringsAsFactors = F)

        print(colnames(daps))
        
        compdap=stack(daps,select= which(colnames(daps) %in% c("Gene.symbol.x", "Gene.symbol.y")))
        compdap=compdap[!duplicated(compdap[,1]),]
        
        #ordering the pair of DAPs
        pares=daps[, which(colnames(daps) %in% c("Gene.symbol.x", "Gene.symbol.y"))]
        colnames(pares)=c("Gene.symbol.x","Gene.symbol.y")
        
        pares=ordpair(pares)
        
        pares=paste(pares[,1],pares[,2])
        
        
        #subsetting only daps genes
        dados_daps=dados[dados$Gene.symbol %in% compdap[,1],]
        ident=data.frame("ID"=row.names(dados_daps),"Gene.symbol"=dados_daps[,"Gene.symbol"],row.names=NULL)
        dados_daps=data.matrix(dados_daps[,-ncol(dados_daps)])
        
        #adding the states
        #if(sum(colnames(dados_daps[,-ncol(dados_daps)])==row.names(facs))==(ncol(dados_daps)-1))
        
        print("pegou so os genes q formam daps")
                
        
        dl=list()
        if(sum(colnames(dados_daps)==samples)==(ncol(dados_daps)))
        {
          dl[[1]]= t(dados_daps[,t(facs2[states_index]) %in% state1])
          dl[[2]]= t(dados_daps[,t(facs2[states_index]) %in% state2])
        } else
        {
          print("problem with sample names")
          write.table(t(c(studies[i],j,"problem with sample names")),status_data,append=T,row.names=F,col.names=F,sep="\t")
          next
        }
        print("separou as classes")
        
        
        cor=lapply(dl, function(x) cor(x[-nrow(x),],use="pairwise.complete.obs"))
        print("fez a correlacao")
        
        dif_daps=t(cor[[1]]-cor[[2]])
        
        genes=colnames(dif_daps)
        
        comb=combn(rev(genes),2)
        print("fez as combinacoes dos pares dos genes")
        
        temp=data.frame("ID.x"=rev(comb[2,]),"ID.y"=rev(comb[1,]))
        temp[,"dif.cor"]=dif_daps[upper.tri(dif_daps)]
        
        #adding gene symbol
        temp=merge(temp,ident,by.y="ID",by.x="ID.x" )
        temp=merge(temp,ident,by.y="ID",by.x="ID.y" )
        
        #removing pairs with same gene symbol
        temp=subset(temp,temp[,5]!=temp[,4])
        
        temp <- data.frame(lapply(temp, as.character), stringsAsFactors=FALSE)
        
        temp=ordpair(temp)
        
        temp[,"pair"]=paste(temp[,4],temp[,5])
        
        temp[,3]=as.numeric(temp[,3])
        temp = temp[order(temp[,6],- abs(temp[,3]) ), ] #sort by id and abs(value)
        #no need for affymatrix. they should have the same ID
        #temp = temp[ !duplicated(temp[,6]), ] 
        
        temp=subset(temp, temp[,6] %in% pares)
        
        temp=temp[,c("pair","ID.x","ID.y","Gene.symbol.x","Gene.symbol.y","dif.cor")]
        
        if(nrow(temp)!=0) 
        {
          fez[f,c("Array.Express.Data","group")]=c(studies[i],j)
          write.table(t(c(studies[i],j,"done")),status_data,append=T,row.names=F,col.names=F,sep="\t")
        }else
        {
          write.table(t(c(studies[i],j,"no DAPs in this AE dataset")),status_data,append=T,row.names=F,col.names=F,sep="\t")
          next
        }
        
        if(output.option=="single table")
        {
          if(nrow(fez)==1)
          {
            fdif_dap=data.frame(temp)
            colnames(fdif_dap)[nrow(fez)+5]=paste(studies[i],j,sep=".")
          } else
          {  
            
            fdif_dap=merge(fdif_dap,temp,by=c(1,2,3,4,5),all=TRUE)
            colnames(fdif_dap)[nrow(fez)+5]=paste(studies[i],j,sep=".")
            rm(temp)
          }
        } else if(output.option=="one_table_each")
        {
          fdif_dap=data.frame(temp)
          colnames(fdif_dap)[6]=paste(studies[i],j,sep=".")
          setwd(path)
          write.table(fdif_dap,paste("DAPs_Dif_single_",studies[i],"_",j,".txt",sep=""),row.names=F, col.names=T,sep="\t")
          print("escreveu DAPs single table")
        }
      }
      
      
    }
    rm(states_index,facs2,facs,fac2,fac,min_num_states,max_smp_index,labels,AEsetnorm,nspec,control,dl,cor,med,dados,dados_daps)
    if(exists("freq_states")) rm(freq_states)
    
  }
}
      
