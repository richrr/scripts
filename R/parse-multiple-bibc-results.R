args = commandArgs(trailingOnly=TRUE)


# always run it in the directory where you have the bibc results.

patern = args[1] # e.g. "-results.txt"
outfile = args[2] # e.g. "all-bibc-results.csv"

files = list.files(path = ".", pattern = patern)

df = read.csv(files[1], header=F, skip=7)
colnames(df) = c("ID", files[1])
#head(df)

for(f in files[-1]){
    #print(f)
    tmp_df = read.csv(f, header=F, skip=7)
    colnames(tmp_df) = c("ID", f)
    #head(tmp_df)
    df = merge(df, tmp_df, by="ID")
}

write.csv(df, outfile, row.names=F)


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
        return(reldf)
}


reldf = relativize(outfile)

# higher values will give lower ranks with the method (alternatively directly get rank(-x))
# the ties method of min matches excel
reldf[,-1] <- apply(reldf[,-1], 2, function (x) {rank(1/rank(as.numeric(x) , ties.method = "min") , ties.method = "min")})
write.csv(reldf, paste0(outfile,".rel.ranked.csv"), quote=F, row.names=F)


# convert the ranks to perc by dividing rank by max rank
percrankdf = reldf
percrankdf[,-1] = apply(reldf[,-1], 2, function (x) {(as.numeric(x)*100)/max(as.numeric(x))})

MinPercRank = apply(percrankdf[,-1], 1, function (x) min(as.numeric(x)))

FreqRankTop25Perc = apply(percrankdf[,-1], 1, function (x) sum(as.numeric(x) < 25))

percrankdf = cbind(percrankdf, MinPercRank, FreqRankTop25Perc)

head(percrankdf)

write.csv(percrankdf, paste0(outfile,".perc.ranked.csv"), quote=F, row.names=F)



