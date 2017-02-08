# csv files only
# first column has to be ids


# http://stackoverflow.com/questions/9447801/dividing-columns-by-colsums-in-r

args = commandArgs(trailingOnly=TRUE)
#print(length(args))



pow <- function(x=10, y=6) {
   # function to print x raised to the power y
   result <- x^y
   return(result)
}


relativize = function(infile){
	df = read.csv(infile, header=T, check.names=F, row.names=1)
	#reldf = df/sum(df) # this is incorrect since it divides by the sum of the matrix
	reldf = prop.table(as.table(as.matrix(df)), 2)
	ID = rownames(df)
	reldf = cbind(ID, reldf)
	write.csv(reldf, paste(infile, ".rel.csv", sep=''), quote=F, row.names=F)
}


relativize_multply_million = function(infile){
	df = read.csv(infile, header=T, check.names=F, row.names=1)
	#reldf = df/sum(df) # this is incorrect since it divides by the sum of the matrix
	reldf = prop.table(as.table(as.matrix(df)), 2)
	mill_reldf = reldf * pow(10, 6)
	ID = rownames(df)
	mill_reldf = cbind(ID, mill_reldf)
	write.csv(mill_reldf, paste(infile, ".million.rel.csv", sep=''), quote=F, row.names=F)
}


relativize(args[1])

relativize_multply_million(args[1])





