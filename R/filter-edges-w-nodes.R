library(stringr)

args = commandArgs(trailingOnly=TRUE)

# remove pairs if they contain certain nodes. The output of this code can be used in --consistent to create networks. This is different than filter-edges.R which keeps edges that have the nodes



read_file_as_vector = function(infile){
	indata = read.csv(infile, header=FALSE)
	indata = as.vector(indata[!duplicated(indata[,1]),1] )	#delete duplicated array probes
	indata = indata[order(indata)] # ascending sort
	return(indata)
}



remove = read_file_as_vector(args[1])
#head(remove)


consist_elems = read_file_as_vector(args[2])
consist_elems = grep("<==>", consist_elems, value=TRUE) # keep pairs
   
pair = str_split( consist_elems ,"<==>") 
pairs = t(as.data.frame(pair))
colnames(pairs) = c("p1", "p2")
#head(pairs)
print(nrow(pairs))


# keep pairs where both partners ABSENT in vector
if(length(remove) > 0){
    pairs = pairs[ !pairs[,"p1"] %in% remove & !pairs[,"p2"] %in% remove, ]
    print(nrow(pairs))
}

good_pairs = paste(pairs[,"p1"], pairs[,"p2"], sep='<==>')
print(length(good_pairs))
#head(good_pairs)
#write.csv(good_pairs, outputFile, quote=F, row.names = F)
write(good_pairs, paste0(args[2], ".filtered-edges.txt"), sep="\n")
