args = commandArgs(trailingOnly=TRUE)

#usage:
# cd /nfs3/PHARM/Morgun_Lab/richrr/db/amiran_mothur_format
#Rscript blast_hit_for_qiime_otus.R /nfs3/PHARM/Morgun_Lab/richrr/db/amiran_mothur_format/ncbi20.fasta.ncbi.format.fa /nfs3/PHARM/Morgun_Lab/richrr/db/amiran_mothur_format/ncbi20.tax test.qiime_13.8_97_otuid_repseq.txt qiime97_otus.txt

library(Biostrings)
library(seqRFLP)
#library(devtools) # sudo apt-get install libcurl4-openssl-dev libssl-dev
# install_github("mhahsler/rBLAST", force = TRUE)
library(rBLAST)
library(seqinr)

blastdbfile = args[1]
ncbitaxafile = args[2]
qiimequeryotus = args[3]
outstr = args[4]

## check if makeblastdb is correctly installed
Sys.which("makeblastdb")

## see possible arguments
blast_help("makeblastdb")

# the following created three files ending in extensions: nhr, nin, nsq
makeblastdb(blastdbfile, dbtype = "nucl")

# Load blast database
blast_db <- blast(db = blastdbfile)
print(blast_db, info=TRUE)


##############
if(FALSE){
processFile = function(filepath) {
  con = file(filepath, "r")
  id = ''
  seq = ''
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(grepl('>', line)){
        id = line
    } else {
        seq = line
        rep = paste(id, seq, sep='\t')
        print(rep)
        id = ''
        seq = ''
    }
    #print(line)
  }

  close(con)
}

processFile("test.qiime_13.8_97_otuid_repseq.txt")
}
##############

# id and taxonomy file
ncbi20.tax = read.delim(ncbitaxafile, header=F, sep='\t')
rownames(ncbi20.tax) = as.character(ncbi20.tax[,1])
colnames(ncbi20.tax) = c("id", "taxa")
#head(ncbi20.tax)

infile = readDNAStringSet(qiimequeryotus)
#print(infile)

# Function for selected otu and blast database
blast_otu <- function(otu_select, blast_db){
  compare <- predict(blast_db, otu_select)
  res10 = compare[1:10,]
  names = ncbi20.tax[as.character(res10$SubjectID), ]
  res10 = cbind(res10, names)
  return(res10)
}

big_res =  data.frame(matrix(NA, ncol=14, nrow=1))
colnames(big_res) = c("QueryID" , "SubjectID" , "Perc.Ident" , "Alignment.Length" , "Mismatches",  "Gap.Openings" , "Q.start"          , "Q.end" , "S.start" ,  "S.end" ,  "E"   , "Bits" , "id" , "taxa")

top_res = data.frame(matrix(NA, ncol=14, nrow=1))
colnames(top_res) = colnames(big_res)

# get the top 10 and top blast hit for each query
for(idx in 1:length(infile)){
    query = infile[idx]
    res = blast_otu(query,  blast_db)
    big_res = rbind(big_res, res)
    top_res = rbind(top_res, res[1,])
}



write.table(big_res[-1,] , paste0("top10_blast_hits_for_",outstr) , row.names=F, quote=F, sep='\t')
write.table(top_res[-1,] , paste0("top_blast_hits_for_", outstr), row.names=F, quote=F, sep='\t')





