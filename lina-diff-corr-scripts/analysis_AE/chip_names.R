#####################################################################################################################
# script written by Lina Thomas
# date: 11/15/2013
# gets the chip names for all data set we are going to work with
# 
# REQUIREMENTS:
# the folder containing the Array Express data folders MUST NOT have any other folders
# 
# IMPROVEMENTS:
# Still needs to keep going regardless of possible errors
######################################################################################################################

subsetting=function(studies,p,pi){
  
  if(!is.vector(studies)) return("object is nor a vector!!")
  
  if(p>length(studies)) return("p is greater than length of studies!!")
  
  if(pi>p) return("impossible to get a pith part greater than the total number of parts!!")
  
  #divisao inteira
  n=length(studies)%/%p
  # #resto da divisao
  f=(length(studies)%%p)/p
  
  if (pi <= f*p) studies=studies[((pi-1)*n+pi):(pi*(n+1))] else
    studies=studies[((pi-1)*n+1+f*p):(pi*n+f*p)]
  
  return(studies)
  
}

loading_AEdata=function(studies){
  
  files=list.files(paste(path,studies,sep="/"))
  #fileslist=as.list(c(paste(path,studies[i],sep="/"),files))
  
  raw=grep("raw" ,files,value=T)
  
  if(length(raw)==0)
  {
    write.table(t(c(studies[i],1,"no raw data1")),status_data,append=T,row.names=F,col.names=F,sep="\t")
    return("no raw data")
  }
  
  raw=grep("zip",raw,value=T)
  
  if(length(raw)==0)
  {
    write.table(t(c(studies,1,"no raw data2")),status_data,append=T,row.names=F,col.names=F,sep="\t")
    return("no rawdata")
    
  } else if(length(raw)==1)
    rawfiles=data.frame(unzip(paste(path,studies,raw,sep="/"),list=T)[,1],stringsAsFactors = FALSE) else
    {
      print("tem mais de um estudo")
      rawfiles=unzip(paste(path,studies,raw[1],sep="/"),list=T)
      rawfiles <- data.frame(lapply(rawfiles, as.character), stringsAsFactors=FALSE) 
      for(raw_id in 2:length(raw))
      {
        temp= unzip(paste(path,studies,raw[raw_id],sep="/"),list=T)
        temp <- data.frame(lapply(temp, as.character), stringsAsFactors=FALSE) 
        rawfiles=data.frame(c(rawfiles[,1],temp[,1]),stringsAsFactors = FALSE)
      }
    }
  
  
  AEfiles=list(sdrf=grep("sdrf",files,value=T),adf=grep("adf",files,value=T),idf=grep("idf",files,value=T),
               path=paste(path,studies,sep="/"),rawFiles=rawfiles[,1])
  
   
  AEset=ae2bioc(mageFiles = AEfiles)
  
  if(class(AEset)!="list") AEset=list(AEset)
  
  return(AEset)
  
}

###########################################################################################
args=commandArgs(trailingOnly = FALSE)
args
pi=as.numeric(args[2])
p=as.numeric(args[3])
start_studies=as.numeric(args[4])
num_studies=args[5] # either a number or "total"
if(num_studies!="total") num_studies=as.numeric(num_studies)
if(is.na(num_studies))
{
  print("ERROR: num_studies must be either a number or total")
  stop
}
###############################################################################################################################
# libraries

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
##################################################################################################################################

local= "cgrb"

if(local=="home")
  path="C:/Users/Lina/Dropbox/Tese/Dados/Array Express" else
    if(local=="lab")
      path="/home/thomasli/Dropbox/Tese/Dados/Array Express" else
        if(local=="cgrb")
          path="/capecchi/pharmacy/morgunlab/lina/ArrayExpress" 

setwd(path)

#reading all the folders from file
studies=system("ls -F /capecchi/pharmacy/morgunlab/lina/ArrayExpress | grep '/' | grep 'E-'", intern =TRUE)
#getting rid of the character /
studies=gsub("/","",studies)

studies=subsetting(studies,p,pi)

if(num_studies=="total") num_studies=length(studies)

for(i in start_studies:num_studies)
{
  print(paste("study",studies[i]))
  
  AEset=loading_AEdata(studies[i])
      
  if(class(AEset[[1]])=="AffyBatch")  print("fez AEset") else next
  
  for(j in 1: length(AEset))
  {
    print(AEset[[j]]@annotation)
    write.table(t(c(studies[i],j,AEset[[j]]@annotation)),"chip_information2.txt",append=T,row.names=F,col.names=F,sep="\t")
    print("fez o write.table")
  }
}

# # 2nd part
#
# chip_table=read.table("chip_information2.txt",header=T)
# 
# chip=unique(chip_table[,3])
# 
# chip <- data.frame(lapply(chip, as.character), stringsAsFactors=FALSE) 
# 
# chip=chip[1,]
# 
# chip=chip[-grep("Error",chip)]
# chip=chip[-grep("ArrayExpress:",chip)]
# chip=chip[-grep("No rawdata",chip)]
# chip=chip[-grep("package",chip)]
# chip=chip[!is.na(chip)]
# 
# 
# for(i in 1: length(chip))
# {
#   if (!paste(chip[i],".db",sep="") %in% installed.packages()) 
#   {
#     source("http://bioconductor.org/biocLite.R")
#     biocLite(paste(chip[i],".db",sep=""))
#   }
# }
# 
# chip_table2=chip_table[chip_table[,3] %in% chip[paste(chip,".db",sep="") %in% installed.packages()],]
# 
# write.table(chip_table2, "AEmusmusculus20smp_affy_chipinst.txt",row.names=F)
