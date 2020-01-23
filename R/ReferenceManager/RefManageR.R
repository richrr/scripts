#Takes Input File 
args = commandArgs(trailingOnly=TRUE)


#Installs Packages
if (!require (RefManageR)){
  system("mkdir ~/R/x86_64-redhat-linux-gnu-library/3.3")
  install.packages(c("RefManageR", "RCurl","stringr","XML","entrez","rentrez","reutils"), repos='http://cran.us.r-project.org', lib = "~/R/x86_64-redhat-linux-gnu-library/3.3")
}


#Obtains functions of packages
library(RefManageR)
library(RCurl)
library(stringr)
library(rentrez)
library(reutils)
library(XML)

###Read Manuscript 
file1 = args[1]
rf1 = readLines(file1)

###Obtain DOI, PMID, PMCID 
##DOI
Temp1 = unlist(str_extract_all(rf1, "\\{DOI:.*?\\}"))
AllDOIs = gsub("\\{DOI:\\s*(.*)}", "\\1", Temp1)
AllDOIs = unique(AllDOIs)
#print(AllDOIs)
write(AllDOIs, file = "DOIS.txt", ncolumns = 1, append = TRUE)

##doi
Temp1 = unlist(str_extract_all(rf1, "\\{doi:.*?\\}"))
AllDOIs = gsub("\\{doi:\\s*(.*)}", "\\1", Temp1)
#print(AllDOIs)
AllDOIs = unique(AllDOIs)
write(AllDOIs, file = "DOIS.txt", ncolumns = 1, append = TRUE)

##DOI Without Header
Temp2 = unlist(str_extract_all(rf1, "\\{10\\..*\\/.*\\}"))
AllDOIs2 = gsub("\\{\\s*(.*)}", "\\1", Temp2)
#print(AllDOIs2)
AllDOIs2 = unique(AllDOIs2)
write(AllDOIs2, file = "DOIS.txt", ncolumns = 1, append = TRUE)


##PMID
Temp3 = unlist(str_extract_all(rf1, "\\{PMID:.*?\\}"))
AllPMIDS = gsub("\\{PMID:\\s*(.*)}", "\\1", Temp3)
AllPMIDS = unique(AllPMIDS)
#print(AllPMIDS)
write(AllPMIDS, file = "PMIDS.txt", ncolumns = 1, append = TRUE)


#For PMID without header 
Temp4 = unlist(str_extract_all(rf1, "\\{[0-9]+.*?\\}"))
AllPMIDS2 = gsub("\\{\\s*(.*)}", "\\1", Temp4)
AllPMIDS2 = unique(AllPMIDS2)

# remove if they have '/' e.g. in doi:
AllPMIDS2 = AllPMIDS2[!grepl('/', AllPMIDS2)]

#print(AllPMIDS2)
write(AllPMIDS2, file = "PMIDS.txt", ncolumns = 1, append = TRUE)

##PMCID Without Header  
Temp6 = unlist(str_extract_all(rf1, "\\{PMC.*?\\}"))
AllPMCIDS2 = gsub("\\{PMCID:\\s*(.*)}", "\\1",Temp6)
AllPMCIDS3 = gsub("\\{\\s*(.*)}", "\\1", AllPMCIDS2)
AllPMCIDS3 = unique(AllPMCIDS3)
#print(AllPMCIDS3)
write(AllPMCIDS3, file = "PMCIDS.txt", ncolumns = 1, append = TRUE)

##PMCID
Temp5 = unlist(str_extract_all(rf1, "\\{PMCID:.*?\\}"))
AllPMCIDS = gsub("\\{PMCID:\\s*(.*)}", "\\1",Temp5)
AllPMCIDS = unique(AllPMCIDS)
#print(AllPMCIDS)
write(AllPMCIDS, file = "PMCIDS.txt", ncolumns = 1, append = TRUE)

 
###Note###
###Before proceeding, check these text files to make sure the proper IDS are in the proper places. Fix files as necessary.


print("DOI")
###Citations For DOIS
fdoi = "DOIS.txt"
rfdoi = readLines(fdoi)
#cit1 = GetBibEntryWithDOI(rfdoi)
file = "citations1.bib"
empty_lines = grepl('^\\s*$', rfdoi)
rfdoi = rfdoi[! empty_lines]

for(id1 in rfdoi){
  #id1 = rfdoi[i]
  print(id1)
  cit1 = GetBibEntryWithDOI(id1, temp.file = "citations1.bib", delete.file = FALSE)
  print(cit1)
  a = WriteBib(cit1, file = "citations1.bib")
  #print(a)
  try(system("bib2xml citations1.bib | xml2end > citations1.end", ignore.stderr=T))
  try(system("cat citations1.end >> export1.end",ignore.stderr=T))
}



###Citations For PMIDS
print("PMID")
fpmid = "PMIDS.txt"
rfpmid = readLines(fpmid)
empty_lines = grepl('^\\s*$', rfpmid)
rfpmid = rfpmid[! empty_lines]
#print(rfpmid)
file.create = "citations2.bib"
for(i in rfpmid){
  print(i)
  cit2 = tryCatch( {GetPubMedByID(i, db = "pubmed")}, error = function(err) {
  # error handler picks up where error was generated
  print(paste("MY_ERROR:  ",err))
})
  print(cit2)
  b = tryCatch({WriteBib(cit2, file = "citations2.bib")}, error = function(err) {
  # error handler picks up where error was generated
  print(paste("MY_ERROR:  ",err))
})
  #print(b)
  try(system("bib2xml citations2.bib | xml2end > citations2.end", ignore.stderr=T))
  #file2 = "citations2.end"
  #write(x = paste("%M", i, sep = " "), file = "export2.end", append = TRUE)
  try(system("cat citations2.end >> export2.end", ignore.stderr=T))
}

system("cat export1.end export2.end > finalexport.end")

###Citations For PMCIDS
print("PMCID")
fpmcid = "PMCIDS.txt"
rfpmcid = readLines(fpmcid)
empty_lines3 = grepl('^\\s*$', rfpmcid)
rfpmcid = rfpmcid[! empty_lines3]
if(length(rfpmcid) > 0)
{
 for(i in rfpmcid){
  fileURL = gsub(x = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=PMC14901", pattern = "PMC14901", replacement = i)
  #fileURL = paste0("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=", i)     ###########################################
  print(url)
  xData <- getURL(fileURL)
  xmldoc <- xmlParse(xData)
  xml_data = xmlToList(xmldoc)
  pmidinxml = as.list(xml_data[["record"]])
  finalpmid = pmidinxml$.attrs["pmid"]
  print(finalpmid)
  cit3 = GetPubMedByID(finalpmid, db = "pubmed")
  print(cit3)
  c = WriteBib(cit3, file = "citations3.bib")
  print(c)
  system("bib2xml citations3.bib | xml2end > citations3.end")
  file3 = "citations3.end"
  write(x = paste("%M", finalpmid, sep = " "), file = "export3.end", append = TRUE)
  write(x = paste("%2", i, sep = " "), file = "export3.end", append = TRUE)
  
  system("cat citations3.end >> export3.end")
  
  write(x = finalpmid, file = "PMIDS2.txt", append = TRUE)
 }
 system("cat export3.end >> finalexport.end")
}
#system("cat export1.end export2.end export3.end > finalexport.end")

system("sed -i 's/ %0/%0/g' finalexport.end")

