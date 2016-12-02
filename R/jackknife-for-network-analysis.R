



library(argparser)
library(hash)
library(psych)
library(gtools)
library(corpcor)
library(reshape)




#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('perform jackknife analyses.R for calculates comparison and correlation analysis using random resampling without replacing in each iteration')
p <- add_argument(p, "expressionDataFile", help="file containing values (measurements) to be used for correlation") # required

# truly optional
p <- add_argument(p, "--output", help="output file", default="./result/network_")
p <- add_argument(p, "--symbolColumnName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID") # 'UniqueID'
p <- add_argument(p, "--warnings", help="print warnings", flag=TRUE)
p <- add_argument(p, "--noCorrelationsRequested", help="This is a preliminary analysis and no correlations are requested. Use this if you have a huge (e.g. >10K input elements as rows)", flag=TRUE)

# # required; but written as optional format so I explicitly mention the following in cmd line
p <- add_argument(p, "--lists", help="lists (of measurement labels e.g. taxa names, phenotypics tests) to generate pairs for calculate correlation", nargs=2) 
p <- add_argument(p, "--AnalysToDoList", help="AnalysToDoList", default="AnalysToDoList.txt") # tab-delimited; perform correlations and/or comparisons as per input
p <- add_argument(p, "--mapFile", help="tab delimited mapFile", default="mapFile.txt") # use the mapping file to identify the samples group affiliation to perform analysis
p <- add_argument(p, "--mapColumns", nargs=2, help="Columns in map file in which to search the sample group associations", default=c("#SampleID","Experiment")) # use the columns in the mapping file to identify the samples group affiliation to perform analysis
p <- add_argument(p, "--pairedInfoColumn", help="Column in map file in which to search the mouse ids. This info will help to identify samples from same mice", default="Mouse") # Column in map file in which to search the mouse ids. This info will help to identify samples from same mice, useful when doing paired comparisons

## optional methods for t test and correlations
p <- add_argument(p, "--comparMethod", help="method to compa A vs B", default="tt") # allowed: "tt" (Welch t test), "mw" (Mann-Whitney U Test)
p <- add_argument(p, "--correlMethod", help="method to corre A vs B", default="pearson") # allowed: pearson, "kendall", or "spearman", can be abbreviated.

p <- add_argument(p, "--sampleSize", help="number of samples per group", default=25) 
p <- add_argument(p, "--numbjackknife", help="number of times to jackknife", default=1000) 

p <- add_argument(p, "--excludeGroup", help="the group to be excluded from jackknife") 

argv <- parse_args(p)
#print (p)


expressionDataFile = argv$expressionDataFile

if(length(argv$lists)!= 2)
{
  print("Lists are required. Give --lists list1 list2 in cmd")
  quit()
}
list1File = argv$lists[1]
list2File = argv$lists[2]
outputFile = argv$output
symbleColumnName =  argv$symbolColumnName
AnalysFile = argv$AnalysToDoList

# the samples in samplIdCol belong to the groups in exptCol
samplIdCol = argv$mapColumns[1]
exptCol = argv$mapColumns[2]
sampsize = argv$sampleSize
numbbootstr = argv$numbjackknife

#print(sampsize)
#print(numbbootstr)
#print(argv$excludeGroup)

#==================================================================================================================
# create a group_sampleid dictionary
#==================================================================================================================
mapFile=read.delim(argv$mapFile,header = TRUE,sep="\t",check.names=FALSE)
dict = hash()
# for unique items in the experiment column
for(k in levels(mapFile[, exptCol])){
    # extract the rows whose col value in exptCol matches k
    # assign the sample names (i.e. values in samplIdCol) to key k
    dict[k] = as.vector(mapFile[which(mapFile[,exptCol]==k), samplIdCol]) 
}

# create if outdir doesn't exist:
dir.create('./rnd-maps/', recursive = TRUE)


# create a string for the random mapping files
#rnd_map_file = paste(argv$mapFile, "rand", sep='.')
rnd_map_file = paste('./rnd-maps/', "rand", sep='')
rnd_map_file = gsub('../', '', rnd_map_file, fixed = T)

#print(dict)
#print(length(dict))



cmds = ''
for (numb in 1:numbbootstr) {

    # create random sample to group maps
    out_df = data.frame()
    for (name in names(dict)){
    
    	    if(!is.na(argv$excludeGroup)) { if(name == argv$excludeGroup) next }
    	    #print(name)
            values = dict[[name]]
            #print(values)

            rnd_samp = ''
	    if(length(values) <= sampsize){
	        rnd_samp = sample(values) # a random permutation
	    } else {
	        rnd_samp = sample(values, sampsize) # sample without replacement
	    }
	    out_df = rbind(out_df, data.frame(rnd_samp,name))
    }
    colnames(out_df) = c(samplIdCol, exptCol)
    #print(out_df)
    rnd_map_file_numb = paste(c(rnd_map_file, numb, "txt"), collapse='.')
    outputFile_rnd_numb = paste(c("sampSize-", sampsize, "-numbjackknife-", numbbootstr, "/", outputFile, "/rnd-", numb, "/results_"), collapse='')
    
    write.table(out_df, rnd_map_file_numb, sep="\t", row.names=FALSE, quote=FALSE)
    
    
    # no need to keep only the required samples in the file. the mapping file decides which samples will be used for analysis
    
    if(argv$noCorrelationsRequested){
    	str = paste(c('Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/perform-analyses.R', argv$expressionDataFile, '--lists', list1File, list2File, '--mapFile', rnd_map_file_numb, '--mapColumns', samplIdCol, exptCol, '--AnalysToDoList', AnalysFile, '--comparMethod', argv$comparMethod, '--correlMethod', argv$correlMethod, '--symbolColumnName', symbleColumnName, '--pairedInfoColumn', argv$pairedInfoColumn, '--output', outputFile_rnd_numb, '--noCorrelationsRequested'), collapse=' ')
    } else {
    	str = paste(c('Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/perform-analyses.R', argv$expressionDataFile, '--lists', list1File, list2File, '--mapFile', rnd_map_file_numb, '--mapColumns', samplIdCol, exptCol, '--AnalysToDoList', AnalysFile, '--comparMethod', argv$comparMethod, '--correlMethod', argv$correlMethod, '--symbolColumnName', symbleColumnName, '--pairedInfoColumn', argv$pairedInfoColumn, '--output', outputFile_rnd_numb), collapse=' ')    
    }
    
    cmds = paste(cmds, str, sep='\n')
    #print(str)

}

write(cmds, "cmds-to-execute.txt")

print(getwd())
out_str = 'Run the following command in the above directory on Needleman: SGE_Array -c cmds-to-execute.txt -q transkingdom -m 25G -f 750G -b 25'
print(out_str)




