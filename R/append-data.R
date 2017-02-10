# needs two csv files

args = commandArgs(trailingOnly=TRUE)
#print(length(args))

infile1 = args[1]

df1 = read.csv(infile1, header=T, check.names=F, row.names=1)
rownames(df1) <- paste("OTU", rownames(df1), sep = "_")

infile2 = args[2]

df2 = read.csv(infile2, header=T, check.names=F, row.names=1)

outdf = rbind(df1, df2)


ID = rownames(outdf)
outdf = cbind(ID, outdf)
write.csv(outdf, paste(infile1, ".pheno.csv", sep=''), quote=F, row.names=F)

