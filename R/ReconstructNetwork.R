# function: calculate correlation between two lists
# input files: 
	#	1.  two lists to generate pair for calculate correlation
	#	2.	a "gene expression" file
# output file:
	#	correlation coefficient,pvalue, fdr ;   for each pair of "genes"
#install.packages("stringr")
# note:
#   if there are lines with the same "ID"(gene symble), the lines other than the first are deleted

# [rodrrich@transkingdom forMichael]$ Rscript ReconstructNetwork.R ./AllControl4Strains.csv --lists ./partner1.txt ./partner2.txt --foldchange ./FoldChange.csv --mapFile mapFile.txt --mapColumns #SampleID Experiment --AnalysToDoList AnalysToDoList.txt

# Rscript ReconstructNetwork.R ./AllControl4Strains.csv --lists ./partner1.txt ./partner2.txt --foldchange ./FoldChange.csv --mapFile mapFile.txt --mapColumns #SampleID Experiment --AnalysToDoList AnalysToDoList.txt --comparMethod mw --correlMethod sp --symbolColumnName #SampleID --genNetwork

library(argparser)
library(hash)
library(psych)
library(gtools)
library(corpcor)
library(reshape)

# load the functions from a different file to keep code neater
source("/nfs1/Morgun_Lab/richrr/scripts/R/comp-correl-delta-paired-util-functs.R")



#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Calculate comparison and correlation analysis between lists')
p <- add_argument(p, "expressionDataFile", help="file containing values (measurements) to be used for correlation") # required
p <- add_argument(p, "--lists", help="lists (of measurement labels e.g. taxa names, phenotypics tests) to generate pairs for calculate correlation", nargs=2) # required; but written as optional format so I explicitly mention the "--lists"
p <- add_argument(p, "--output", help="output file", default="./result/network_")

#==================================================================================================================
# network generation parameters
#==================================================================================================================
p <- add_argument(p, "--foldchange", help="file containing foldchange to be used for correlation")
p <- add_argument(p, "--indivPvalCutoff", help="individualPvalueCutoff", default=0.3, type="numeric")
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=1, type="numeric")
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.3, type="numeric")
p <- add_argument(p, "--symbolColumnName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID") # 'UniqueID'

# # required; but written as optional format so I explicitly mention the following in cmd line
p <- add_argument(p, "--AnalysToDoList", help="AnalysToDoList", default="AnalysToDoList.txt") # tab-delimited; perform correlations or comparisons as per inmput
p <- add_argument(p, "--mapFile", help="tab delimited mapFile", default="mapFile.txt") # use the mapping file to identify the samples group affiliation to perform correlations
p <- add_argument(p, "--mapColumns", nargs=2, help="Columns in map file in which to search the sample group associations", default=c("#SampleID","Experiment")) # use the columns in the mapping file to identify the samples group affiliation to perform correlations
p <- add_argument(p, "--pairedInfoColumn", help="Column in map file in which to search the mouse ids. This info will help to identify samples from same mice", default="Mouse") # Column in map file in which to search the mouse ids. This info will help to identify samples from same mice



## optional methods for t test and correlations
p <- add_argument(p, "--comparMethod", help="method to compa A vs B", default="tt") # allowed: "tt" (Welch t test), "mw" (Mann-Whitney U Test)
p <- add_argument(p, "--correlMethod", help="method to corre A vs B", default="pearson") # allowed: pearson, "kendall", or "spearman", can be abbreviated.


p <- add_argument(p, "--genNetwork", help="generate Network", flag=TRUE)
p <- add_argument(p, "--dontgenPUC", help="don't generate PUC", flag=TRUE)



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

# create if outdir doesn't exist:
odir =  strsplit(outputFile, "/")[[1]]
odir = paste(odir[1:(length(odir) -1)], collapse='/')
print(odir)
dir.create(odir, recursive = TRUE)


split= argv$split
NumberOfSplit= argv$numberofsplit


individualPvalueCutoff = argv$indivPvalCutoff
combinedPvalueCutoff = argv$combPvalCutoff
combinedFDRCutoff = argv$combFDRCutoff


##### to do: 
## partial correlation


symbleColumnName =  argv$symbolColumnName

AnalysToDoList = read.delim(argv$AnalysToDoList,header = FALSE)

# the samples in samplIdCol belong to the groups in exptCol
samplIdCol = argv$mapColumns[1]
exptCol = argv$mapColumns[2]

# create a group_sampleid dictionary
mapFile=read.delim(argv$mapFile,header = TRUE,sep="\t",check.names=FALSE)
dict = hash()
# for unique items in the experiment column
for(k in levels(mapFile[, exptCol])){
    # extract the rows whose col value in exptCol matches k
    # assign the sample names (i.e. values in samplIdCol) to key k
    dict[k] = as.vector(mapFile[which(mapFile[,exptCol]==k), samplIdCol]) 
}



# check the number and order of samples assigned to groups in dict. MAY be as per the order in which the samples are added (order in the mapping file)
for(k in keys(dict)){
    str_r_l = paste(k, length(dict[[k]]), "values", paste(as.vector(dict[[k]]), collapse=", "), sep="---") 
    print(str_r_l)
}


# check if the order of samples in the two categories is correct; affects any time there is delta, paired, timelag
checkPairing = function(catg1, catg2, dict, mapFile, pairedInfoColumn, samplIdCol){
	s1 = as.vector(dict[[catg1]])
        s2 = as.vector(dict[[catg2]])

        tmp_mapFile = as.matrix(mapFile)
        row.names(tmp_mapFile) = mapFile[,samplIdCol]
        #print(tmp_mapFile)
        m1 = as.vector(tmp_mapFile[s1, pairedInfoColumn]) 
        print(m1)
        m2 = as.vector(tmp_mapFile[s2, pairedInfoColumn]) 
        print(m2)

        bool = all.equal(m1,m2)

        if(bool == "TRUE"){
            print("Passed pairing check")
        }else{
            print("XXXXXXXXXXXXXX----Failed pairing check----XXXXXXXXXXXXXX")
            print(bool)
            quit()
        }
        return(1)
}


for(indx in 1:nrow(AnalysToDoList)){
    c = as.vector(AnalysToDoList[indx,])
    
    #print(c)
    if(c[3] == "paired" || c[3] == "delta" || c[3] == "timelag"){
    print(c)

      if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == ""){
        print("2")
        categA= toString(c[1,"V1"])
        categB= toString(c[1,"V2"])

        cA = strsplit(categA, ",")[[1]]
        categ1 = cA[1]
        categ2 = cA[2] 
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

        cB = strsplit(categB, ",")[[1]]
        categ3 = cB[1]
        categ4 = cB[2]
        checkPairing(categ3, categ4, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

      } # end delta if


      if(c[1] != "" && c[2] != "" && c[3] == "paired" && c[4] == ""){
        print("3")
        categ1= toString(c[1,"V1"])
        categ2= toString(c[1,"V2"])
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

      } # end paired if


      if(c[1] != "" && c[2] == "" && c[3] == "delta" && c[4] == "correlation"){
        print("5")
        categA= toString(c[1,"V1"])
        cA = strsplit(categA, ",")[[1]]
        categ1 = cA[1]
        categ2 = cA[2]
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

      } # end delta correl if


      if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == "correlation"){
        print("7")
        categA= toString(c[1,"V1"])
        categB= toString(c[1,"V2"])

        cA = strsplit(categA, ",")[[1]]
        categ1 = cA[1]
        categ2 = cA[2]
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)
 
        cB = strsplit(categB, ",")[[1]]
        categ3 = cB[1]
        categ4 = cB[2]
        checkPairing(categ3, categ4, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

      } # end delta correl comparison if


      if(c[1] != "" && c[2] != "" && c[3] == "timelag" && c[4] == "correlation"){
        print("9")
        # std. correlations: g1T1-g2T1 or g1T2-g2T2
        # calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc. 
        # for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene
        categ1= toString(c[1,"V1"])
        categ2= toString(c[1,"V2"])
        checkPairing(categ1, categ2, dict, mapFile, argv$pairedInfoColumn, samplIdCol)

      } # end timelag if

      
    } # end outer if

}# end for


#quit()


###########################################################################################################
#	output filename
networkFile = paste(outputFile,"individualPvalue:",individualPvalueCutoff,"_combinedPvalue:",combinedPvalueCutoff,"_combinedFDR:",combinedFDRCutoff,".csv",sep="")
pairFile = "./testpair.txt"


#==================================================================================================================
# generate or load existing pairlist to avoid re-generating pairs(needs time) (for correlations and comparison of correlations)
#==================================================================================================================
list1 = read.csv(list1File,header = FALSE,sep="\t")
list2 = read.csv(list2File,header = FALSE,sep="\t")

genes1 = as.vector(list1[!duplicated(list1[,1]),1] )	#delete duplicated array probes
genes1 = genes1[order(genes1)] # ascending sort
genes2 = as.vector(list2[!duplicated(list2[,1]),1] )
genes2 = genes2[order(genes2)]

if (!file.exists(pairFile)){
	if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then generate unique pairs
		#pairs = t(combn(genes1,2))[,2:1]
                # I am using this https://gist.github.com/randy3k/10015496
                pairs = combinations(length(unique(genes1)), 2, genes1, repeats.allowed=TRUE)
	}else{
		pairs = expand.grid(genes1,genes2)  # this has all the above & their opposite (e.g. gene1 gene2 and gene2 gene1) & same gene pairs (e.g. gene1 gene1)
	}	
	write.csv(pairs,pairFile,row.names=FALSE)
}else{
	pairs = read.csv(pairFile)
}


#==================================================================================================================
# load expression data
#==================================================================================================================
expressionData = read.csv(expressionDataFile,header = TRUE,check.names=FALSE,na.strings=c("","na","NA"))
expressionData = expressionData[!duplicated(expressionData[,symbleColumnName]),]    #delete duplicated measurements (genes, taxa labels, phenotypic labels, etc.)
rownames(expressionData) = expressionData[,symbleColumnName]


calcFC = function(v1 , v2){
    v1 = as.numeric(as.vector(v1))
    v2 = as.numeric(as.vector(v2))

    FC = c() # append did not work, don't know why
    for(i in 1:length(v1)){
        if(is.na(v1[i]) || is.na(v2[i])){
            FC = c(FC, NA)
        }else{
            fc = as.numeric(v1[i])/as.numeric(v2[i])
            FC = c(FC, fc)
        }
    }
    #print(FC)
            
    return(FC)
}


n=0 # total number of analysis

corrout = c()
compout = c()
folchout = c()

# conditions: 
# add 4th column with "correlation" to indicate whatever operation is to be done, needs to be done on correlations
# 1. A vs B: 2 non-empty columns, last two columns empty [A	B empty empty]
# 2. delta A vs delta B: 3 non-empty columns, 3rd column delta, 4th column empty [A	B	delta	empty]  # needs 4 groups: two comma separated from A and two comma separated from B
# 3. A vs B paired:  3 non-empty columns, 3rd column paired, 4th column empty [A	B	paired	empty]
# 4. Correl A: 2 non-empty columns, 4th column correlation [A	empty	empty	correlation]
# 5. Correl delta A: 3 non-empty columns, 2nd column empty, 3rd column delta, 4th column correlation [A	empty	delta	correlation] # needs 2 groups: two comma separated from A 
# 6. Correl A vs. Correl B: 3 non-empty columns, 4th column correlation [A	B	empty	correlation]
# 7. Correl delta A vs. Correl delta B: 4 non-empty columns, 3rd column delta, 4th column correlation [A	B	delta	correlation] # needs 4 groups: two comma separated from A and two comma separated from B
# 8. Correl A vs. Correl B paired: 4 non-empty columns, 3rd column paired, 4th column correlation [A	B	paired	correlation]
## no need to do paired in 8
# 9. Timelag Correl A and B (where A and B are different time points): 4 non-empty columns, 3rd column timelag, 4th column correlation [A	B	timelag	correlation]
## std. correlations: g1T1-g2T1 or g1T2-g2T2
## calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc. 
## for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene


for(indx in 1:nrow(AnalysToDoList)){
    c = as.vector(AnalysToDoList[indx,])
    
    #print(c)
    if(c[1] != "" && c[2] != "" && c[3] == "" && c[4] == ""){
        print("1")
        categ1= toString(c[1,"V1"])
        categ2= toString(c[1,"V2"])

        genes = ''
        if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
	    genes = genes1
	}else{                              # else take union
	    genes = union(genes1,genes2)  
	} # end else

        # calculate comparison
        outEachExperiment = calculateComparison(genes, expressionData, categ1, categ2, dict, argv$comparMethod, c[3], indx)

        if (length(compout)==0){
	    compout <<- outEachExperiment
        }else{
	    compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
        } # end else

        # calculate fold change
        FolCh =  calcFC(outEachExperiment[,2],outEachExperiment[,3])
        FolCh = t(rbind(outEachExperiment[,1], t(FolCh)))
        colnames(FolCh) = c(colnames(outEachExperiment)[1],paste(categ1,"vs",categ2,"FoldChange",sep=" "))

        if (length(folchout)==0){
	    folchout <<- FolCh
        }else{
	    folchout <<- merge(folchout,FolCh,sort=FALSE,by=1)
        } # end else

    }else if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == ""){
        print("2")
        categA= toString(c[1,"V1"])
        categB= toString(c[1,"V2"])

        cA = strsplit(categA, ",")[[1]]
        categ1 = cA[1]
        categ2 = cA[2]
 
        cB = strsplit(categB, ",")[[1]]
        categ3 = cB[1]
        categ4 = cB[2]
 
        if(categ1 == "" || categ2 == "" || categ3 == "" || categ4 == ""){
            print("Error in number of categories. Give categ1,categ2\tcateg3,categ4")
            next
        }

        genes = ''
        if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
	    genes = genes1
	}else{                              # else take union
	    genes = union(genes1,genes2)  
	} # end else

        # calculate comparison
        outEachExperiment = calculateComparisonDelta(genes, expressionData, categ1, categ2, categ3, categ4, dict, argv$comparMethod, indx)
        if (length(compout)==0){
	    compout <<- outEachExperiment
        }else{
	    compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
        } # end else

        # calculate fold change
        FolCh =  calcFC(outEachExperiment[,2],outEachExperiment[,3])
        FolCh = t(rbind(outEachExperiment[,1], t(FolCh)))
        colnames(FolCh) = c(colnames(outEachExperiment)[1],paste(categA,"vs",categB,"FoldChange",sep=" "))

        if (length(folchout)==0){
	    folchout <<- FolCh
        }else{
	    folchout <<- merge(folchout,FolCh,sort=FALSE,by=1)
        } # end else

    }else if(c[1] != "" && c[2] != "" && c[3] == "paired" && c[4] == ""){
        print("3")
        categ1= toString(c[1,"V1"])
        categ2= toString(c[1,"V2"])

        genes = ''
        if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
	    genes = genes1
	}else{                              # else take union
	    genes = union(genes1,genes2)  
	} # end else

        # calculate comparison
        outEachExperiment = calculateComparison(genes, expressionData, categ1, categ2, dict, argv$comparMethod, c[3], indx)
        if (length(compout)==0){
	    compout <<- outEachExperiment
        }else{
	    compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
        } # end else

        # calculate fold change
        FolCh =  calcFC(outEachExperiment[,2],outEachExperiment[,3])
        FolCh = t(rbind(outEachExperiment[,1], t(FolCh)))
        colnames(FolCh) = c(colnames(outEachExperiment)[1],paste(categ1,"vs",categ2,"FoldChange",sep=" "))

        if (length(folchout)==0){
	    folchout <<- FolCh
        }else{
	    folchout <<- merge(folchout,FolCh,sort=FALSE,by=1)
        } # end else

    }else if(c[1] != "" && c[2] == "" && c[3] == "" && c[4] == "correlation"){
        print("4")
        categ= toString(c[1,"V1"])
        # calculate correlations
        outEachExperiment = calculateCorrelation(pairs, expressionData, categ, dict, argv$correlMethod, indx)
        if (length(corrout)==0){
	    corrout <<- outEachExperiment
        }else{
	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
        } # end else
        
        genes = ''
        if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
	    genes = genes1
	}else{                              # else take union
	    genes = union(genes1,genes2)  
	} # end else

        # calculate partial correlations only for the given list of genes (all pairs of the list)
        # instead of calc PC with the 35K genes which will be lot of (unnecessary) computation
        if(FALSE){
        outEachExperiment = calculatePartialCorrelation(pairs, genes, expressionData, categ, dict, argv$correlMethod, indx)
        if (length(corrout)==0){
    	    corrout <<- outEachExperiment
        }else{
            corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
        } # end else
        # the above does not have gene1-gene2 and gene2-gene1
        }

    }else if(c[1] != "" && c[2] == "" && c[3] == "delta" && c[4] == "correlation"){
        print("5")
        categA= toString(c[1,"V1"])
        cA = strsplit(categA, ",")[[1]]
        categ1 = cA[1]
        categ2 = cA[2]
  
        if(categ1 == "" || categ2 == ""){
            print("Error in number of categories. Give categ1,categ2")
            next
        }

        # calculate correlations
        outEachExperiment = calculateDeltaCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, indx)
        if (length(corrout)==0){
	    corrout <<- outEachExperiment
        }else{
	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
        } # end else


    }else if(c[1] != "" && c[2] != "" && c[3] == "" && c[4] == "correlation"){
        print("6")
        # compare correlations of gene pairs
        categ1= toString(c[1,"V1"])
        categ2= toString(c[1,"V2"])

        # calculate correlations
        outEachExperiment = calculateComparisCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, c[3], indx)
        if (length(corrout)==0){
	    corrout <<- outEachExperiment
        }else{
	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
        } # end else



    }else if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == "correlation"){
        print("7")
        categA= toString(c[1,"V1"])
        categB= toString(c[1,"V2"])

        cA = strsplit(categA, ",")[[1]]
        categ1 = cA[1]
        categ2 = cA[2]
 
        cB = strsplit(categB, ",")[[1]]
        categ3 = cB[1]
        categ4 = cB[2]
 
        if(categ1 == "" || categ2 == "" || categ3 == "" || categ4 == ""){
            print("Error in number of categories. Give categ1,categ2\tcateg3,categ4")
            next
        }

        # calculate correlations
        outEachExperiment = calculateComparisDeltaCorrelation(pairs, expressionData, categ1, categ2, categ3, categ4, dict, argv$correlMethod, indx)
        if (length(corrout)==0){
	    corrout <<- outEachExperiment
        }else{
	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
        } # end else


    }else if(c[1] != "" && c[2] != "" && c[3] == "paired" && c[4] == "correlation"){
        print("8")

    }else if(c[1] != "" && c[2] != "" && c[3] == "timelag" && c[4] == "correlation"){
        print("9")
        # std. correlations: g1T1-g2T1 or g1T2-g2T2
        # calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc. 
        # for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene
        categ1= toString(c[1,"V1"])
        categ2= toString(c[1,"V2"])
        # calculate correlations
        outEachExperiment = calculateTwoCategCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, indx)
        if (length(corrout)==0){
	    corrout <<- outEachExperiment
        }else{
	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
        } # end else


    }else{
        print("Unknown")
    }


}
print(warnings())

FoldChangeFile = argv$foldchange
if(is.na(FoldChangeFile)){
    FoldChangeFile =  paste(outputFile,"output.csv",sep='folch-')
}
write.csv(corrout,paste(outputFile,"output.csv",sep='corr-'),row.names=FALSE)
write.csv(compout,paste(outputFile,"output.csv",sep='comp-'),row.names=FALSE)
write.csv(folchout,FoldChangeFile,row.names=FALSE)


# if do NOT calc PUC is requested, then quit, default it will calc PUC
if(argv$dontgenPUC){
    quit()
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# calculate PUC 
	# 1. calculate PUC from correlation coefficient and foldchange information
# input File:
	# 1. output from 20130417CalculationCorrelation: correlation coefficient, pvalue, fdr for each strain(datasets)
# output:
	#  combined pvalue, combined fdr,PUC added to (correlation coefficient, pvalue, fdr for each dataset,)

# a fold change file:
# which has genes (phenotypes) as rows
# columns are the foldchange for each comparison in the analysis file

# the Data file which has the correlation coeff, pval, etc. for multiple correlations
# you need to grep the required columns


Data = corrout
dataName = ''
# calculate PUC
#==================================================================================================================
forPUC = function(FoldChangeFile, categ1, categ2){
	FoldChangeCol = paste(categ1,"vs",categ2,"FoldChange",sep=" ")
        geneIDCol = "geneName"

	#change out format for PUC
	library(stringr)
	pair = str_split(Data[,"pairName"],"<==>")
	pairs = t(as.data.frame(pair))
	#pairs = str_replace_all(pairs,perl("\\.\\d"),"")
	colnames(pairs) = c("partner1","partner2")
	#outForPUC = apply(outForPUC, 2, function(x) as.numeric(as.character(x))) # convert the chars to numeric
        Data = apply(Data, 2, function(x) as.numeric(as.character(x))) # convert the chars to numeric
        outForPUC = cbind(pairs,Data)
	
        g_grep_cols = c("partner1","partner2")
        for(categ in c(categ1, categ2)){
            grep_cols_c = grep(paste(categ," ",sep="") , colnames(outForPUC), value=TRUE, fixed=TRUE)
            for(suff_str in c("Coefficient" , "pvalue")){
                #regex_patt = paste(categ, suff_str , sep = ' .*')
                grep_cols = grep(suff_str , grep_cols_c, value=TRUE)
                # remove the 'pvalue' from comparisons:
                grep_cols = grep(" vs " , grep_cols, value=TRUE, invert=TRUE)
                # remove the 'timelag' correlations:
                grep_cols = grep(" timelag " , grep_cols, value=TRUE, invert=TRUE)
                # remove the 'delta' correlations if not present in input:
                if( !(grepl(" minus ", categ)) ){
                    grep_cols = grep(" minus " , grep_cols, value=TRUE, invert=TRUE)
                }
                #print(grep_cols)
                
                l_identical_inp = outForPUC[,grep_cols[1]]
                if(length(grep_cols)>1){
                    for(z in 2:length(grep_cols)){
                        l_identical_inp = cbind(l_identical_inp, outForPUC[,grep_cols[z]])
                    }

                    
                    #if(all(sapply(list(l_identical_inp), FUN = identical, outForPUC[,grep_cols[1]]))){
                    if(all(apply(l_identical_inp, 2, identical, outForPUC[,grep_cols[1]]))){
                        #print("All are equal")
                    } else {
                        print(grep_cols)
                        print(head(l_identical_inp))
                        print("All are NOT equal")
                        quit()
                    }
                }
                # double check that the vectors are the same, if true, add the first one
                g_grep_cols = c(g_grep_cols, grep_cols[1]) # just add the first one, since the others are just duplicates
            }
        }   

        #print(g_grep_cols)

        outForPUC = outForPUC[,g_grep_cols]
        #print(head(outForPUC))
        #quit()
#==================================================================================================================
#    						generate foldChange

        if( grepl(" minus ", FoldChangeCol) ){
                    FoldChangeCol = gsub(" minus " , ",", FoldChangeCol)
        }

	FoldChangeMetabolic =  read.csv(FoldChangeFile,header = TRUE,check.names=FALSE)
	#print(head(FoldChangeMetabolic))
 	FoldChangeMetabolic = FoldChangeMetabolic[!duplicated(FoldChangeMetabolic[,geneIDCol]),]
	rownames(FoldChangeMetabolic) = FoldChangeMetabolic[,"geneName"]
	
	FoldMetab1_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner1"]),c(geneIDCol,FoldChangeCol)]
	colnames(FoldMetab1_InPair) = c("partner1InFold","partner1_FoldChange")
	#print(head(FoldMetab1_InPair))
	
	FoldMetab2_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner2"]),c(geneIDCol,FoldChangeCol)]
	colnames(FoldMetab2_InPair) = c("partner2InFold","partner2_FoldChange")
	#print(head(FoldMetab2_InPair))

	
	outForPUC = cbind(outForPUC,FoldMetab1_InPair,FoldMetab2_InPair)
#==================================================================================================================
#    					calculate PUC
	coefficientColnames = colnames(outForPUC)[grep("Coefficient",colnames(outForPUC))]                                  ############ to be fixed
	interestedCoefficientColnames = coefficientColnames[sapply(dataName,grep,coefficientColnames)]
	#print(interestedCoefficientColnames)
	#print (head(outForPUC))
        #quit()
	#error
	interestedCorrelationData = outForPUC[,interestedCoefficientColnames]
	# calculate correlation Direction For correlation coefficient of interest
	interestedCorrelationData = as.matrix(interestedCorrelationData)
	interestedCorrelationData = apply(interestedCorrelationData,2,function(x){as.numeric(as.vector(x))})
	#print(dim(interestedCorrelationData))
	#error
	signOfInterestedCorrelationData = interestedCorrelationData/abs(interestedCorrelationData)					# forOutput
	rownames(signOfInterestedCorrelationData) = c()
	colnames(signOfInterestedCorrelationData) = paste(colnames(interestedCorrelationData),"correlationDirection",sep=".")

        #print(signOfInterestedCorrelationData)
	# see if direction match in interest
	ifMatch = apply(signOfInterestedCorrelationData,1
					,function(tt){
								if (length(tt)==1){
									return (1)
								}
								if (length(tt[!is.na(tt)])<2){
									return (NA)
								}
								uniqueValues = unique(tt[!is.na(tt)])
								NumberOfUniqueValues = length(uniqueValues)
								if (NumberOfUniqueValues>1){
									return (0)
								}else if(NumberOfUniqueValues==0){
									return (NA)
								}else{
									return (1*uniqueValues[1])
								}
							}
					)	# if each row of signOfInterestedCorrelationData is like -1,1,-1; if match, "ifMatch" will be 3 or -3

	names(ifMatch) = c()
        #print(ifMatch)  # forOutput
	
	# calculate matched correlation direction
	matchedExpressionDirection = ifMatch*signOfInterestedCorrelationData[,1]
	names(matchedExpressionDirection) = c()
	matchedExpressionDirection	= ifMatch	# forOutput	
        #print(colnames(outForPUC))
	# calculate fold change direction	
	FoldChangeColnames = colnames(outForPUC)[grep("FoldChange",colnames(outForPUC))]
	FoldChangeData = outForPUC[,FoldChangeColnames]
	FoldChangeDirection = (FoldChangeData-1)/abs(FoldChangeData-1)
	names(FoldChangeDirection) = c()
	colnames(FoldChangeDirection) = paste(colnames(FoldChangeData),"FoldChangeDirection",sep=".")
	FoldChangeDirection				## forOutput	
	# calculate if fold change direction are the same for the two partners
	IfFoldChangeDirectionMatch = apply(FoldChangeDirection,1,prod)
	names(IfFoldChangeDirectionMatch) = c()
	IfFoldChangeDirectionMatch			## forOutput	
	# use "matched correlation direction" and "if fold change direction are the same" to calculate PUC
	PUC = IfFoldChangeDirectionMatch * matchedExpressionDirection    				## forOutput	
	outForPUC = cbind(outForPUC,signOfInterestedCorrelationData,ifMatch,matchedExpressionDirection,FoldChangeDirection,IfFoldChangeDirectionMatch,PUC)

	if(FALSE){
	# calculate combined Pvalue for interest
	PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
	interestedPvalueColnames = PvalueColnames[sapply(dataName,grep,PvalueColnames)]
	interestedPvalueData = outForPUC[,interestedPvalueColnames]
	interestedPvalueData = as.matrix(interestedPvalueData)
	interestedPvalueData = apply(interestedPvalueData,2,function(x){as.numeric(as.vector(x))})
	combinedPvalue = apply(interestedPvalueData,1
							,function(pvalues){
										pvalues = pvalues[!is.na(pvalues)]
										statistics = -2*log(prod(pvalues))
										degreeOfFreedom = 2*length(pvalues)
										combined = 1-pchisq(statistics,degreeOfFreedom)
									}
							)
	outForPUC = cbind(outForPUC,combinedPvalue)
        }
                
	return(outForPUC)
		
} 

#=============================================================================================
zzzz = 0
PUCoutfile = paste(outputFile,"output.csv",sep='PUC-')

if(file.exists(PUCoutfile))
{
    print("Deleted old PUC output")
    file.remove(PUCoutfile)
}

for(indx in 1:nrow(AnalysToDoList)){
    c = as.vector(AnalysToDoList[indx,])
    zzzz = zzzz+1
    if(c[1] != "" && c[2] != "" && c[4] != "correlation"){
        if(c[1] == "DUMMY" || c[2] == "DUMMY" || c[3] == "DUMMY" || c[4] == "DUMMY"){
            next
        }
        print(c)
    
        categ1= toString(c[1,"V1"])
        categ2= toString(c[1,"V2"])

        # if the analysis is delta:
        if(grepl("," , categ1)){
           categ1 = gsub(",", " minus ", categ1)
        }
        if(grepl("," , categ2)){
           categ2 = gsub(",", " minus ", categ2)
        }

        result = forPUC(FoldChangeFile, categ1, categ2)

        # calculate percentage of unexpected correlations
        print("Total items")
        DEN = length(result[,"PUC"])
        print(DEN)

        print("Items not 1")
        NUM = DEN - length(result[which(result[,"PUC"]==1), "PUC"])       
        print(NUM)
        
        PUC_Prop = as.numeric(NUM*100/DEN)
        print(PUC_Prop)
        PUC_Prop = paste(c("Percentage of Unexpected Correlation=" , PUC_Prop, "%") , collapse='')
       


        #print(result[,c(3,5)])
        #write.csv(result,paste(c(outputFile,zzzz,"output.csv"),collapse='PUC-'),row.names=FALSE,append=TRUE)
        if(file.exists(PUCoutfile)){
          out = file(PUCoutfile ,'a')
          write.csv(result, file=out, row.names=FALSE)
          close(out)
        } else{
          write.csv(result, PUCoutfile ,row.names=FALSE)
        }

        if(file.exists(PUCoutfile)){
          out = file(PUCoutfile ,'a')
          write.csv(PUC_Prop, file=out, row.names=FALSE)
          close(out)
        }

    }

}
file.remove("./testpair.txt")

print("Done with PUC")


# if generate network is not requested, then quit
if(!argv$genNetwork){
    quit()
}




#calculate FDR for combined pvalue
combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
result = cbind(result,combinedFDR)


###########################################################################################################
#generate network
###########################################################################################################
data = result
#print(head(data))
#error
generateNetwork = function(){
	#generate network After Filter

	# find PUC expected and significant in combinedPvalue and significant in combinedFDR
	out = data[data[,"PUC"]==1 ,]
	out = out[(out[,"combinedFDR"]<combinedFDRCutoff)==1,]
	# find significant in individual pvalue: all of the pvalues for all of the datasets for each pair must be smaller than threshold
	# find pvalue data
	pvalueData = out[,grep("pvalue",colnames(out))]
	pvalueData = as.matrix(pvalueData)
	# calculate the largest pvalue among all datasets for each gene, this smallest pvalue must be smaller than threshold
	passIndevidualPvalue = apply(pvalueData,1,max)<individualPvalueCutoff
	outNetwork = out[passIndevidualPvalue,]


	write.csv (outNetwork,networkFile,row.names=FALSE)
}
###########################################################################################################


generateNetwork()




























