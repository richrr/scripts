

library(stringr)

args = commandArgs(trailingOnly=TRUE)

indir = args[1]

files = list.files(indir, pattern = ".txt",  full.names = TRUE)
print(paste(c("Processing ", length(files), " files"), collapse=''))

print(files[1])
f1 = read.delim(files[1], header=F, sep="\t")

# get sample name
group=str_split(files[1],"/")[[1]]
sample_name = tail(group, n=1)

colnames(f1) = c("ID", sample_name)

out=f1

# file 1 is already read above
for(f in files[2:length(files)]){
    print(f)
    df = read.delim(f, header=F, sep="\t")

    # get sample name 
    group=str_split(f,"/")[[1]]  
    sample_name = tail(group, n=1) 

    colnames(df) = c("ID", sample_name) 

    out=merge(out, df, by="ID", all=T)
}

write.csv(out, "suumarized-htseq-results.csv", row.names=FALSE, quote=F)

