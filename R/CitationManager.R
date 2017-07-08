#setwd("~/EndNote Project/")
args = commandArgs(trailingOnly=TRUE)


#library(devtools)
library(RefManageR)
library(RCurl) #dependency: bitops

#Read input text file for DOI or PubMed ID
#Change file1 to open your text file; assuming no header
file1 = args[1]
d = readLines(file1)
#print(d)

#To Get BibEntry with DOI 
#if (interactive() && url.exists("http://dx.doi.org/")){
#  GetBibEntryWithDOI(d,"trial2.txt", FALSE)
#}


res = GetPubMedByID(d, db = "pubmed")
WriteBib(res, "citations.bib")


if(FALSE){
for(i in 1:length(d)){
    id = d[i]
	print(id)
	res = GetPubMedByID(id, db = "pubmed")
    print(res)
	

}
}