####################################################################################################################
# script written by Lina Thomas
# date: 2/8/2014
# downloads data from Array Express website into the labs computer and transfers to cgrb: capecchi
# 
# it is developed to devide the file with all studies into p parts and run the pith part separately and set to run in 
# CentOS terminal. 1st arg to type is which part and 2nd arg to type is how many parts

## R < ~/Dropbox/Tese/scripts/studies/1_Prepare_for_analysis/combo.R --no-save 1 10 1273 total  "chip_information_homosapiens.txt"



###################################################################################################################

# deviding studies in p parts to run faster!!!

args = commandArgs(trailingOnly = F)

#which part
pi=as.numeric(args[3])

#how many parts we are going to devide studies
p = as.numeric(args[4])

#start1 = as.numeric(args[5])

#end = args[6]
#if(end != "total") end = as.numeric(end) else
 # end = nrow(studies)

start2 = args[5]
if(start2 != "lastone") start2 = as.numeric(start2)

local = args[6]
chip_info = args[7]
# chip_info = "chip_information_homosapiens.txt"


######################################################################################################################

library("ArrayExpress")
library("affy")

if(local=="home")
  path="C:/Users/Lina/Dropbox/Tese/Dados/ArrayExpress" else
    if(local=="lab")
      path="/home/thomasli/Dropbox/Tese/Dados/ArrayExpress/homosapiens" else
        if(local=="cgrb")
          path="/capecchi/pharmacy/morgunlab/lina/ArrayExpress" 

# labs folder path

if(local == "home")
  path_dir = "C:/Users/Lina/Dropbox/Tese/Dados/ArrayExpress" else
    if(local == "lab")
      path_dir = "/media/Storage/ArrayExpress" else
        if(local == "cgrb")
          path_dir = "/capecchi/pharmacy/morgunlab/lina/ArrayExpress" 

#############################################################################
# functions
if(local == "lab")
  source("/home/thomasli/Dropbox/Tese/scripts/functions.R") else
    if(local == "cgrb")
      source("/raid1/home/pharmacy/thomasli/scripts/functions.R")

AE.loading.from.folder = function(study, path, chip_info){
  
  files = list.files(paste(path,study,sep = "/"))
  #fileslist = as.list(c(paste(path,study,sep = "/"),files))
  
  raw = grep("raw" , files,value = T)
  if(length(raw) == 0)
  {
    write.table(t(c(study, 1, "no raw data")), chip_info,append = T, row.names = F,col.names = F)
    return("no raw data")
  }
  raw = grep("zip", raw,value = T)
  
  print(length(raw))
  if(length(raw) == 0)
  {
    write.table(t(c(study,1,"no raw data")), chip_info,append = T, row.names = F, col.names = F)
    return("no raw data")
  }
  
  if(length(raw) == 1)
    rawfiles = data.frame(unzip(paste(path, study, raw,sep = "/"), list = T)[, 1],stringsAsFactors  =  FALSE) else
    {
      rawfiles = unzip(paste(path, study, raw[1], sep = "/"), list = T)
      rawfiles = data.frame(lapply(rawfiles, as.character), stringsAsFactors = FALSE) 
      for(raw_id in 2:length(raw))
      {
        temp =  unzip(paste(path, study, raw[raw_id], sep = "/"),list = T)
        temp = data.frame(lapply(temp, as.character), stringsAsFactors = FALSE) 
        rawfiles = data.frame(c(rawfiles[, 1], temp[,1]), stringsAsFactors  =  FALSE)
      }
    }
  #print("erro aqui 1")
  
  AEfiles = list(sdrf = grep("sdrf",files,value = T), adf = grep("adf",files,value = T),idf = grep("idf",files,value = T),
                 path = paste(path, study, sep = "/"), rawFiles = rawfiles[,1])
  print(AEfiles)
  
  AEset = ae2bioc(mageFiles  =  AEfiles)
  
  return(AEset)
}


##########################################################################################

setwd(path)

studies=read.csv("AE_homosapiens20smp_affy.csv")

studies <- data.frame(lapply(studies, as.character), stringsAsFactors=FALSE)

# ##############################################################################################
# 
# start1 = 1273
# end = nrow(studies)
# p = 1
# pi = 1

folders.capecchi = list.dirs("/capecchi/pharmacy/morgunlab/lina/ArrayExpress", full.names = F)
folders.capecchi = substr(folders.capecchi, nchar(folders.capecchi[1]) + 2, nchar(folders.capecchi))
folders.capecchi = folders.capecchi[-grep("Florey", folders.capecchi)]

#print(folders.capecchi)

downloaded.studies = which(studies[,1] %in% folders.capecchi)
#print(downloaded.studies)

studies = studies[ downloaded.studies, 1]
print(studies) 

#studies = subsetting(studies[start1:end, 1], p, pi)
# start2 = 519
#########################################################################################################################

if(start2 == "lastone") start2 = length(studies)

#start2 = which(studies == "E-MEXP-3252")+1

setwd(path_dir)
for(i in start2:length(studies))
{
  print(paste("study", studies[i]))
  #if(studies[i] %in% c("E-GEOD-45565", "E-GEOD-46969", "E-GEOD-42731", "E-GEOD-45328", "E-GEOD-44170","E-GEOD-44874","E-GEOD-43079","E-GEOD-37138", "E-GEOD-39280", "E-GEOD-42986", "E-GEOD-47968", "E-GEOD-40967", "E-GEOD-40966", "E-GEOD-44994","E-GEOD-44739", "E-GEOD-42200", "E-GEOD-39368","E-GEOD-39367", "E-GEOD-44434", "E-MTAB-1215",  "E-GEOD-36059", "E-GEOD-41646", "E-GEOD-42869", "E-GEOD-38073", "E-GEOD-41117", "E-GEOD-42245","E-GEOD-41130","E-GEOD-36908", "E-GEOD-41475", "E-GEOD-33477", "E-GEOD-36458", "E-GEOD-37384","E-GEOD-36403","E-GEOD-36190", "E-GEOD-36139","E-GEOD-34996","E-GEOD-34992","E-GEOD-33377","E-GEOD-31620", "E-GEOD-30777", "E-GEOD-18876","E-GEOD-28015","E-GEOD-28013", "E-GEOD-27940","E-GEOD-32519","E-GEOD-29449","E-GEOD-31174","E-GEOD-24188","E-GEOD-25104","E-GEOD-25099","E-GEOD-24391","E-GEOD-24135","E-MEXP-2777","E-GEOD-31424","E-GEOD-26736","E-GEOD-29619","E-GEOD-30483","E-GEOD-30453","E-GEOD-30422","E-GEOD-18668","E-GEOD-28111","E-GEOD-22455","E-GEOD-24551","E-GEOD-24550","E-GEOD-24549","E-GEOD-29156","E-GEOD-26841","E-GEOD-28030","E-GEOD-28814","E-GEOD-22874","E-GEOD-22862","E-GEOD-26020","E-GEOD-16008","E-GEOD-27608","E-GEOD-27342","E-MEXP-2644","E-GEOD-27388","E-GEOD-27251","E-GEOD-27002","E-GEOD-26173","E-GEOD-26440","E-GEOD-24779","E-GEOD-21780","E-GEOD-21541","E-GEOD-21450","E-MTAB-336","E-GEOD-21107","E-MTAB-181","E-GEOD-21713","E-GEOD-24558","E-GEOD-24182","E-MEXP-2767","E-GEOD-16778","E-GEOD-21140","E-GEOD-23980","E-GEOD-16157","E-GEOD-7822","E-TABM-779","E-GEOD-17895","E-GEOD-19612","E-GEOD-21163","E-GEOD-3929","E-GEOD-12815","E-GEOD-16125","E-GEOD-14352","E-GEOD-14471","E-MEXP-2361","E-GEOD-19475","E-GEOD-14264","E-GEOD-19090","E-GEOD-18111","E-GEOD-17602","E-GEOD-17601","E-GEOD-13159", "E-GEOD-25219","E-GEOD-17778","E-GEOD-16226", "E-GEOD-14771","E-GEOD-10592","E-GEOD-13989","E-GEOD-10300","E-GEOD-6269","E-GEOD-14535","E-GEOD-14323","E-GEOD-9649","E-GEOD-13373","E-GEOD-13372","E-GEOD-12711","E-GEOD-12709","E-GEOD-11915", "E-GEOD-12627","E-MEXP-1328","E-GEOD-12236","E-GEOD-11909","E-GEOD-12760","E-GEOD-8126","E-GEOD-8798","E-GEOD-8721", "E-GEOD-8650","E-GEOD-8192","E-GEOD-7880","E-GEOD-7152","E-GEOD-7127","E-GEOD-6751","E-GEOD-6521", "E-GEOD-5949","E-GEOD-6120","E-GEOD-2350","E-GEOD-10325","E-GEOD-9385","E-TABM-314","E-GEOD-473","E-GEOD-5388","E-GEOD-7785","E-TABM-244","E-MEXP-44","E-EMBL-6","E-GEOD-2109", "E-GEOD-8921", "E-GEOD-7952")) next  

  if(local == "lab")
    chip_table = read.table("/media/Storage/ArrayExpress/chip_information_homosapiens.txt") else
      if(local == "cgrb")
	chip_table = read.table("/capecchi/pharmacy/morgunlab/lina/ArrayExpress/chip_information_homosapiens.txt")

  colnames( chip_table ) = c("study", "group", "chip")

  #print(head(chip_table))
  
  #if( studies[i] %in% chip_table[, "study"]) next
 
  if(local == "lab") 
    if(!studies[i] %in% list.dirs(path_dir, full.names = F) )
    {
      # streaming from website
      print(studies[i])
      dir.create(path = paste(path_dir,studies[i],sep = "/"),showWarnings = TRUE)
      AEfiles = getAE(studies[i],path = paste(path_dir,studies[i],sep = "/"))
    }
  
  # loanding and getting chip

  result = tryCatch({
    AEset = AE.loading.from.folder(studies[i], path_dir, chip_info)
  }, warning = function(war){
    print(paste(war))
  }, error = function(err){ 
    print(paste(err))
  }, finally = {
    if(exists("AEset")) print("fez AEset")
  } 
  )

  print(exists("AEset"))

  if(!exists("AEset")) 
  {
    print("aqui1")
    allocate = grep("Error in ae2bioc", result, value = T)
    print("aqui1.2")
    
    if( length(allocate) == 0)
    {
print("aqui 2")
      if(length(result) == 1)
        write.table( t(c(studies[i], "deleted", result)), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t") else
          if("message" %in% names(result))
            write.table( t(c(studies[i], "deleted", result$message)), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t") else
            {
print("aqui 3")
              paste_result = result[[1]]
              
              for(index_result in 2 : length(result))
                paste(paste_result, result[[index_result]])
              
              write.table( t(c(studies[i], "deleted", result)), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t")
            } 
      print("aqui 4")
	
      remove = paste("rm -rf ",path_dir,"/",studies[i], sep = "")
	print(remove)

      system(remove)
      
      rm(allocate)
    
    } else
    {
       print("aqui 5")
       write.table( t(c(studies[i], "hand check", result)), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t")
    }
  rm(AEset, result)
    
  } else if(is.character(AEset))
  {
   if(length(result) == 1)
        write.table( t(c(studies[i], "deleted", result)), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t") else
          if("message" %in% names(result))
            write.table( t(c(studies[i], "deleted", result$message)), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t") else
            {
print("aqui 3")
              paste_result = result[[1]]
              
              for(index_result in 2 : length(result))
                paste(paste_result, result[[index_result]])
              
              write.table( t(c(studies[i], "deleted", result)), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t")
            } 

  } else
  {
    print("aqui 6")
    print(AEset)
    
    if(class(AEset)!= "list") AEset = list(AEset)
    
    #getting chip infor
    
    for(j in 1: length(AEset))
    {
      check = 0
      print(paste("j =",j))
      AEset_annotation=AEset[[j]]@annotation
      print(AEset_annotation)
      if( length(AEset_annotation) == 0) AEset_annotation = NA
      
      if(!studies[i] %in% chip_table[,"study"])
        write.table(t(c(studies[i],j,AEset_annotation)), chip_info, row.names = F, col.names = F, append = T, sep = "\t")    
      print("done!")
      
      if (!paste(AEset_annotation,".db",sep="") %in% installed.packages()) 
      {
        source("http://bioconductor.org/biocLite.R")
        biocLite(paste(AEset_annotation,".db",sep=""))
      }
      
      if (paste(AEset_annotation,".db",sep="") %in% installed.packages()) 
      {
        AEsetnorm = rma(AEset[[j]])
        setwd(path_dir)
        save.image(file = paste(studies[i],"_",j,".RData",sep=""))

        write.table( t(c(studies[i], "RData in capecchi", "might be useful")), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t")
        
        copy = paste("mv ",path_dir,"/",studies[i],"_",j,".RData","/capecchi/pharmacy/morgunlab/lina/ArrayExpress/RData/HomoSapiens", sep="")
        system(copy)
        
        #remove = paste("rm ",path_dir,"/",studies[i],"_",j,".RData", sep = "")
        #system(remove)
    
	rm(AEsetnorm, copy, remove)
	
      }  else 
      {
	 check = check + 1
	 write.table( t(c(studies[i], "deleted", "unable to install chip")), "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/AE_experiments_location.txt", append = T, col.names = F, row.names = F, sep = "\t")

	 remove = paste("rm ",path_dir,"/",studies[i],"_",j,".RData", sep = "")
         system(remove)

	 remove = paste("rm -rf ",path_dir,"/",studies[i], sep = "")
         system(remove)

         next
      }
      rm(AEset_annotation)
    }
    
    remove = paste("rm -rf ",path_dir,"/",studies[i], sep = "")
    system(remove)
    
    rm(AEset,result)  
    
  }
}


