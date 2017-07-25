#####################################################################################################################
# script written by Lina Thomas
# date: 11/7/2013
# downloads data from Array Express website into the labs computer and transfers to cgrb: capecchi
# 
# it is developed to devide the file with all studies into p parts and run the pith part separately and set to run in 
# CentOS terminal. 1st arg to type is which part and 2nd arg to type is how many parts
#
# cgrb: SGE_Batch -c 'R < scripts/transfer_file_linux.R --no-save 1 10 cgrb 3 total file "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/" "AEhomosapiens20smp_affy_chipinst.txt"  "thomasli@florey.cgrb.oregonstate.edu" "~/ArrayExpress/data/HomoSapiens/"' -r transfer_file -f 8G -m 8G -q samwise
# lab: R < ~/Dropbox/Tese/scripts/100_sets/transfer_file_linux.R --no-save 1 10 lab 1 total folder "thomasli@files.cgrb.oregonstate.edu" "/capecchi/pharmacy/morgunlab/lina/ArrayExpress/"  "thomasli@florey.cgrb.oregonstate.edu" "~/ArrayExpress/data/HomoSapiens/"
######################################################################################################################

args <- commandArgs(trailingOnly = F)

args

#how many parts we are going to devide studies  
pi = as.numeric(args[3])
#which part
p = as.numeric(args[4])
local = args[5]

start = args[6]
end = args[7]
if(end != "total") end = as.numeric(end)

from.where = args[8]
from.folder = args[9]
from.file = args[10]
to.host = args[11]
to.folder = args[12]


###############################################################################################

#install packages in crick!!!

library("ArrayExpress")

if(local=="home")
  path="C:/Users/Lina/Dropbox/Tese/Dados/Array Express" else
    if(local=="lab")
    {
      path="/home/thomasli/ArrayExpress" 
      source("/home/thomasli/Dropbox/Tese/scripts/functions.R")
    } else if(local=="cgrb")
    {
      path="/capecchi/pharmacy/morgunlab/lina/ArrayExpress" 
      source("/raid1/home/pharmacy/thomasli/scripts/functions.R")
    }
      

setwd(path)

if(from.where == "file") 
{
  studies = read.table(paste(from.folder,from.file,sep=""),header=T)
  studies = data.frame(lapply(studies, as.character), stringsAsFactors=FALSE)
  studies = studies[!duplicated(studies[,1]),1]
  print(studies)
} else
  if(from.where == "folder")
    if(local == "cgrb") 
    {
      files_from_folder = paste("ls -Fx ",from.folder, " | grep '/' | grep 'E-' ", sep="")
      
      print(files_from_folder)
      
      studies = system(files_from_folder, intern =TRUE)
      
      studies = gsub("/","",studies)
      
    }else if(local == "lab")
    {
      studies = list.dirs(from.folder, full.names=FALSE)
      studies = grep("E-",studies,value=T)
    }

###################################################################################################################

# deviding studies in p parts to run faster!!!
length(studies)
if(end == "total") end = length(studies)

studies= subsetting(studies,p,pi)
print(studies)
#########################################################################################################################

# transfering

for(i in start:end)
{
  
  print(studies[i])
  
  copy=paste("scp -r ",from.folder , studies[i], " ", to.host ,":",to.folder, sep="")
  print(copy)
  
  #copying files to /capecchi using scp
  system(copy)
  
  print("copied!")
  #deleting files from lab computer
  remove=paste("rm -rf ",from.folder, studies[i],"/", sep="")
  print(remove)
  system(remove)
  print("removed!")
  
}


#########################################################################################################################






