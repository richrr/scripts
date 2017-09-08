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
	# this and most of the code in TransNet keeps quote=F, so if the id names have
	# comma it could be problematic.
	# this part removes ',' from the id names, so that the csv files can be used
	# later without running into issues.
	ID = gsub(",", ".", rownames(df), fixed=TRUE)
	reldf = cbind(ID, reldf)
	write.csv(reldf, paste(infile, ".rel.csv", sep=''), quote=F, row.names=F)
}


relativize_multply_million = function(infile){
	df = read.csv(infile, header=T, check.names=F, row.names=1)
	#reldf = df/sum(df) # this is incorrect since it divides by the sum of the matrix
	reldf = prop.table(as.table(as.matrix(df)), 2)
	mill_reldf = reldf * pow(10, 6)
	# this and most of the code in TransNet keeps quote=F, so if the id names have
	# comma it could be problematic.
	# this part removes ',' from the id names, so that the csv files can be used
	# later without running into issues.
	ID = gsub(",", ".", rownames(df), fixed=TRUE)
	mill_reldf = cbind(ID, mill_reldf)
	write.csv(mill_reldf, paste(infile, ".rel.million.csv", sep=''), quote=F, row.names=F)
}


relativize(args[1])

relativize_multply_million(args[1])





