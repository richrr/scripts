##function which orders the pair of genes. The matrix A must have two columns named "gene.symbol 1" and "gene.symbol 2"

ordpair=function(A){
  A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
  genecol1=unique(c(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]))
  genecol1=sort(genecol1)
  
  for(k in 1:length(genecol1))
  {
    
    if(k==1)  where=which(A[,"Gene.symbol.y"]==genecol1[k]) else
      where=which(A[,"Gene.symbol.y"]==genecol1[k] & !(A[,"Gene.symbol.x"]%in% genecol1[1:k]))
    
    who=A[where,"Gene.symbol.x"]
    A[where,"Gene.symbol.x"]=c(rep(genecol1[k],length(who)))
    A[where,"Gene.symbol.y"]=who
    A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
    
  } 
  
  return (A)
}

##########################################################################################
# function that divides the vector in the pith part
# studies : vector
# p : how many parts to divide
# pi: the pith part

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

###########################################################################################
loading.studies.chip=function(rerun,p,pi){
  
  studies_chip=read.table("AEmusmusculus20smp_affy_chipinst.txt", header=T)
  
  if(rerun=="yes")
  {
    rerun2=read.table("AE_data_refazer.txt",header=F,stringsAsFactors=F)[,1]
    studies_chip=studies_chip[studies_chip[,1] %in% rerun2,]
  }
  
  chips=unique(studies_chip[,3])
  chips <- as.character(sapply(chips, as.character), stringsAsFactors=FALSE)
  
  for(i in 1:length(chips))
    library(paste(chips[i],".db",sep=""),character.only=T)
  
  studies = unique(studies_chip[,1])
  studies = as.character(sapply(studies, as.character), stringsAsFactors=FALSE)
  
  studies = subsetting(studies,p,pi)
  
  return(studies)
  
}

##########################################################################################

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

########################################################################################

loading.studies.chip3=function(run_data, p, pi, rerun, select.type, rerun_data){
    
    studies_chip = read.table(run_data, header = T)
    
    if(rerun == TRUE)
    {
      rerun2=read.table(rerun_data,header=F,stringsAsFactors=F)[,1]
      if(select.type=="out") studies_chip=studies_chip[!(studies_chip[,1] %in% rerun2),] else
        if(select.type=="in") studies_chip=studies_chip[studies_chip[,1] %in% rerun2,] else
        {
          print("choose in ou out for select.type")
          stop 
        }
    }
    
    chips=unique(studies_chip[,3])
    chips <- as.character(sapply(chips, as.character), stringsAsFactors=FALSE)
    
    for(i in 1:length(chips))
      library(paste(chips[i],".db",sep=""),character.only=T)
    
    studies = unique(studies_chip[,1])
    studies = as.character(sapply(studies, as.character), stringsAsFactors=FALSE)
    
    studies = subsetting(studies,p,pi)
    
    return(studies)
    
  }

