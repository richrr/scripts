# csv files only
# first column has to be ids


args = commandArgs(trailingOnly=TRUE)
#print(length(args))


relativize = function(infile){
	df = read.csv(infile, header=T, check.names=F, row.names=1)
	reldf = df/sum(df)
	ID = rownames(df)
	reldf = cbind(ID, reldf)
	write.csv(reldf, paste(infile, ".rel.csv", sep=''), quote=F, row.names=F)
}


relativize(args[1])






