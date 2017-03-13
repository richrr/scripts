##### to do: 



# Usage: Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/calc-frequency.R ./AllControl4Strains.csv --lists ./partner1.txt ./partner2.txt --mapFile mapFile.txt --mapColumns #SampleID Experiment --AnalysToDoList AnalysToDoList.txt --comparMethod mw --correlMethod sp --symbolColumnName #SampleID 



# Allowed analyses (conditions): 
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



library(argparser)
library(hash)
library(psych)
library(gtools)
library(corpcor)
library(reshape)

# load the functions from a different file to keep code neater
source("/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/comp-correl-delta-paired-util-freq-functs.R")
#source("/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/get_base_stats_low_express_probes.R")


#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('perform-analyses.R Calculates comparison and correlation analysis between lists')
p <- add_argument(p, "expressionDataFile", help="file containing values (measurements) to be used for correlation") # required

# truly optional
p <- add_argument(p, "--output", help="output file", default="./result/network_")
p <- add_argument(p, "--symbolColumnName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID") # 'UniqueID'
p <- add_argument(p, "--warnings", help="print warnings", flag=TRUE)
p <- add_argument(p, "--noCorrelationsRequested", help="This is a preliminary analysis and no correlations are requested. Use this if you have a huge (e.g. >10K input elements as rows)", flag=TRUE)

p <- add_argument(p, "--skipDiffInputsCheck", help="does not quit if diff input check fails.", flag=TRUE) # avoid using this

# # required; but written as optional format so I explicitly mention the following in cmd line
p <- add_argument(p, "--lists", help="lists (of measurement labels e.g. taxa names, phenotypics tests) to generate pairs for calculate correlation", nargs=2) 
p <- add_argument(p, "--AnalysToDoList", help="AnalysToDoList", default="AnalysToDoList.txt") # tab-delimited; perform correlations and/or comparisons as per input
p <- add_argument(p, "--mapFile", help="tab delimited mapFile", default="mapFile.txt") # use the mapping file to identify the samples group affiliation to perform analysis
p <- add_argument(p, "--mapColumns", nargs=2, help="Columns in map file in which to search the sample group associations", default=c("#SampleID","Experiment")) # use the columns in the mapping file to identify the samples group affiliation to perform analysis
p <- add_argument(p, "--pairedInfoColumn", help="Column in map file in which to search the mouse ids. This info will help to identify samples from same mice", default="Mouse") # Column in map file in which to search the mouse ids. This info will help to identify samples from same mice, useful when doing paired comparisons

## optional methods for t test and correlations
p <- add_argument(p, "--comparMethod", help="method to compa A vs B", default="limma") # allowed: "tt" (Welch t test), "mw" (Mann-Whitney U Test), limma
p <- add_argument(p, "--correlMethod", help="method to corre A vs B", default="pearson") # allowed: pearson, "kendall", or "spearman", can be abbreviated.
p <- add_argument(p, "--logbase", help="calc log using the base", default=2) # allowed: 0 (no log), 1 (e), 2, 10
p <- add_argument(p, "--background", help="the background level usually obtained from the raw files. this is used to find low expression probes", default=0)


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

logbase = argv$logbase
backgrond = argv$background

outputFile = argv$output

# create if outdir doesn't exist:
odir =  strsplit(outputFile, "/")[[1]]
odir = paste(odir[1:(length(odir) -1)], collapse='/')
print(odir)
dir.create(odir, recursive = TRUE)


symbleColumnName =  argv$symbolColumnName

AnalysToDoList = read.delim(argv$AnalysToDoList,header = FALSE)

# the samples in samplIdCol belong to the groups in exptCol
samplIdCol = argv$mapColumns[1]
exptCol = argv$mapColumns[2]

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
#print(dict)

#==================================================================================================================
# check the number and order of samples assigned to groups in dict. MAY be as per the order in which the samples are added (order in the mapping file)
#==================================================================================================================
for(k in keys(dict)){
    str_r_l = paste(k, length(dict[[k]]), "values", paste(as.vector(dict[[k]]), collapse=", "), sep="---") 
    print(str_r_l)
}

#==================================================================================================================
# check if the order of samples in the two categories is correct; affects any time there is delta, paired, timelag
#==================================================================================================================
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



#	output filename
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

#print(genes1)
#print(genes2)

pairs = ''
if(!argv$noCorrelationsRequested){
  if (!file.exists(pairFile)){
	#if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then generate unique pairs
	if(identical(genes1, genes2)){
		pairs = t(combn(genes1,2))[,2:1]
		# use below for testing purposes and/or when calculating partial correlation # has same gene pairs (e.g. gene1 gene1)
                #pairs = combinations(length(unique(genes1)), 2, genes1, repeats.allowed=TRUE) # from https://gist.github.com/randy3k/10015496 for testing purposes 
	}else{
		pairs = expand.grid(genes1,genes2)  # this has all the above & their opposite (e.g. gene1 gene2 and gene2 gene1) & same gene pairs (e.g. gene1 gene1)
	}	
	#write.csv(pairs,pairFile,row.names=FALSE)
	#file.remove("./testpair.txt")
  }else{
	pairs = read.csv(pairFile)
  }
}


#==================================================================================================================
# load expression data
#==================================================================================================================
expressionData = read.csv(expressionDataFile,header = TRUE,check.names=FALSE,na.strings=c("","na","NA", "Na", "NaN"))
expressionData = expressionData[!duplicated(expressionData[,symbleColumnName]),]    #delete duplicated measurements (genes, taxa labels, phenotypic labels, etc.)
rownames(expressionData) = expressionData[,symbleColumnName]

## keep only OTUs
 oturows = grep("^OTU", rownames(expressionData), value = T)
 expressionData =  expressionData[oturows ,]

#if(logbase != 0){
#    if(logbase == 1) {
#          expressionData[, !names(expressionData) %in% c(symbleColumnName)] = log(expressionData[, !names(expressionData) %in% c(symbleColumnName)] + 1)  # using the default base e i.e. exp(1)
#    } else {
#          expressionData[, !names(expressionData) %in% c(symbleColumnName)] = log(expressionData[, !names(expressionData) %in% c(symbleColumnName)] + 1, logbase)
#    }
#}

#write.csv(expressionData,paste(c(outputFile,"log", logbase, "transform-output.csv"),collapse='-'),row.names=FALSE)


#==================================================================================================================
# relativize the entire expression data
#==================================================================================================================
#get_base_stats_exp_data(expressionDataFile, backgrond, outputFile, sample_perc_threshold=10)

pow <- function(x=10, y=6) {
  # function to print x raised to the power y
  result <- x^y
  return(result)
}


relativize = function(df){
  reldf = prop.table(as.table(as.matrix(df)), 2)
  ID = rownames(df)
  reldf = cbind(ID, reldf)
  write.csv(reldf, paste(outputFile, ".rel.csv", sep=''), quote=F, row.names=F)
  symbleColumnName <<- "ID"
  return(reldf)
}


relativize_multply_million = function(df){
  reldf = prop.table(as.table(as.matrix(df)), 2)
  mill_reldf = reldf * pow(10, 6)
  ID = rownames(df)
  mill_reldf = cbind(ID, mill_reldf)
  write.csv(mill_reldf, paste(outputFile, ".rel.million.csv", sep=''), quote=F, row.names=F)
  symbleColumnName <<- "ID"
  return(mill_reldf)
}

expressionData = relativize(expressionData[,-1])
#expressionData = relativize_multply_million(expressionData[,1])




# which elements from list1 are not present in list2
checkDiffInputs = function(list1, list2, search, inn){
    diff = setdiff(list1, list2)
    diff = diff[!diff %in% c('DUMMY')]  # remove known tmp values
    if(length(diff) > 0){
           strr = paste("\n\nThe following ids in the", search ,"file are not present in the", inn, "file\n", sep=' ')
           cat(strr)
           print(diff)
           if(!argv$skipDiffInputsCheck){
               quit()
           }
    }
    #return(diff)
}


#==================================================================================================================
# check if all ids from the different files are consistent. This is especially to check typos or different cases ("ABX" vs "Abx")
#==================================================================================================================
# partners file and expt file to check gene name
checkDiffInputs(genes1, as.vector(expressionData[,symbleColumnName]), "partners", "expression")
checkDiffInputs(genes2, as.vector(expressionData[,symbleColumnName]), "partners", "expression")

# mapping file and expt file to check sample name
checkDiffInputs(as.vector(mapFile[, samplIdCol]), as.vector(colnames(expressionData)), "mapping", "expression")

# analysis file and mapping file to check group names
tmp_a1 = sort(as.vector(AnalysToDoList[,1]))
tmp_a2 = sort(as.vector(AnalysToDoList[,2]))
tmp_a12 = unique(sort(append(tmp_a1, tmp_a2)))

# elements contain ','
csv_tmp_a12 = grep("," , tmp_a12, value = TRUE)
# split based on ','
no_csv_tmp_a12 = unlist(strsplit(csv_tmp_a12, ','))
# remove elements with ',' and replace with the vector obtained by splitting on ','
uniq_analys_elems = append(setdiff( tmp_a12, csv_tmp_a12 ), no_csv_tmp_a12) 

# elements contains ';' in case of anova like analysis
semcol_tmp_a12 = grep(";" , uniq_analys_elems, value = TRUE)
# split based on ';'
no_semcol_tmp_a12 = unlist(strsplit(semcol_tmp_a12, ';'))

# the above may contain elements with "_vs_"
vs_no_semcol_tmp_a12 = grep("_vs_" , no_semcol_tmp_a12, value = TRUE)
# split based on '_vs_'
no_vs_no_semcol_tmp_a12 = unlist(strsplit(vs_no_semcol_tmp_a12, '_vs_'))
# remove elements with '_vs_' and replace with the vector obtained by splitting on '_vs_'
uniq_analys_elems_tmp = append(setdiff( no_semcol_tmp_a12, vs_no_semcol_tmp_a12 ), no_vs_no_semcol_tmp_a12)


# remove elements with ';' and replace with the vector obtained by splitting on ';'
uniq_analys_elems = append(setdiff( uniq_analys_elems, semcol_tmp_a12 ), uniq_analys_elems_tmp) 
#print(uniq_analys_elems)

checkDiffInputs(unique(uniq_analys_elems[uniq_analys_elems != ""]), as.vector(mapFile[, exptCol]), "analysis", "mapping")




n=0 # total number of analysis

corrout = c()
compout = c()
folchout = c()

#==================================================================================================================
# perform analysis
#==================================================================================================================
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


#     }else if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == ""){
#         print("2")
#         categA= toString(c[1,"V1"])
#         categB= toString(c[1,"V2"])
# 
#         cA = strsplit(categA, ",")[[1]]
#         categ1 = cA[1]
#         categ2 = cA[2]
#  
#         cB = strsplit(categB, ",")[[1]]
#         categ3 = cB[1]
#         categ4 = cB[2]
#  
#         if(categ1 == "" || categ2 == "" || categ3 == "" || categ4 == ""){
#             print("Error in number of categories. Give categ1,categ2\tcateg3,categ4")
#             next
#         }
# 
#     genes = ''
# 	#if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then generate unique pairs
# 	if(identical(genes1, genes2)){
# 	    genes = genes1
# 	}else{                              # else take union
# 	    #print("Reached here")
# 	    genes = union(genes1,genes2)  
# 	} # end else
# 
#         # calculate comparison
#         outEachExperiment = calculateComparisonDelta(genes, expressionData, categ1, categ2, categ3, categ4, dict, argv$comparMethod, indx)
#         if (length(compout)==0){
# 	    compout <<- outEachExperiment
#         }else{
# 	    compout <<- merge(compout,outEachExperiment,sort=FALSE,by=1)
#         } # end else
# 
# 
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


#     }else if(c[1] != "" && c[2] == "" && c[3] == "" && c[4] == "correlation"){
#         print("4")
#         categ= toString(c[1,"V1"])
#         # calculate correlations
#         outEachExperiment = calculateCorrelation(pairs, expressionData, categ, dict, argv$correlMethod, indx)
#         if (length(corrout)==0){
# 	    corrout <<- outEachExperiment
#         }else{
# 	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
#         } # end else
#         
#         genes = ''
#         if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then use genes1
# 	    genes = genes1
# 	}else{                              # else take union
# 	    genes = union(genes1,genes2)  
# 	} # end else
# 
#         # calculate partial correlations only for the given list of genes (all pairs of the list)
#         # instead of calc PC with the 35K genes which will be lot of (unnecessary) computation
#         if(FALSE){
#         outEachExperiment = calculatePartialCorrelation(pairs, genes, expressionData, categ, dict, argv$correlMethod, indx)
#         if (length(corrout)==0){
#     	    corrout <<- outEachExperiment
#         }else{
#             corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
#         } # end else
#         # the above does not have gene1-gene2 and gene2-gene1
#         }
# 
#     }else if(c[1] != "" && c[2] == "" && c[3] == "delta" && c[4] == "correlation"){
#         print("5")
#         categA= toString(c[1,"V1"])
#         cA = strsplit(categA, ",")[[1]]
#         categ1 = cA[1]
#         categ2 = cA[2]
#   
#         if(categ1 == "" || categ2 == ""){
#             print("Error in number of categories. Give categ1,categ2")
#             next
#         }
# 
#         # calculate correlations
#         outEachExperiment = calculateDeltaCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, indx)
#         if (length(corrout)==0){
# 	    corrout <<- outEachExperiment
#         }else{
# 	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
#         } # end else
# 
# 
#     }else if(c[1] != "" && c[2] != "" && c[3] == "" && c[4] == "correlation"){
#         print("6")
#         # compare correlations of gene pairs
#         categ1= toString(c[1,"V1"])
#         categ2= toString(c[1,"V2"])
# 
#         # calculate correlations
#         outEachExperiment = calculateComparisCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, c[3], indx)
#         if (length(corrout)==0){
# 	    corrout <<- outEachExperiment
#         }else{
# 	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
#         } # end else
# 
# 
# 
#     }else if(c[1] != "" && c[2] != "" && c[3] == "delta" && c[4] == "correlation"){
#         print("7")
#         categA= toString(c[1,"V1"])
#         categB= toString(c[1,"V2"])
# 
#         cA = strsplit(categA, ",")[[1]]
#         categ1 = cA[1]
#         categ2 = cA[2]
#  
#         cB = strsplit(categB, ",")[[1]]
#         categ3 = cB[1]
#         categ4 = cB[2]
#  
#         if(categ1 == "" || categ2 == "" || categ3 == "" || categ4 == ""){
#             print("Error in number of categories. Give categ1,categ2\tcateg3,categ4")
#             next
#         }
# 
#         # calculate correlations
#         outEachExperiment = calculateComparisDeltaCorrelation(pairs, expressionData, categ1, categ2, categ3, categ4, dict, argv$correlMethod, indx)
#         if (length(corrout)==0){
# 	    corrout <<- outEachExperiment
#         }else{
# 	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
#         } # end else
# 
# 
#     }else if(c[1] != "" && c[2] != "" && c[3] == "paired" && c[4] == "correlation"){
#         print("8")
# 
#     }else if(c[1] != "" && c[2] != "" && c[3] == "timelag" && c[4] == "correlation"){
#         print("9")
#         # std. correlations: g1T1-g2T1 or g1T2-g2T2
#         # calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc. 
#         # for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene
#         categ1= toString(c[1,"V1"])
#         categ2= toString(c[1,"V2"])
#         # calculate correlations
#         outEachExperiment = calculateTwoCategCorrelation(pairs, expressionData, categ1, categ2, dict, argv$correlMethod, indx)
#         if (length(corrout)==0){
# 	    corrout <<- outEachExperiment
#         }else{
# 	    corrout <<- merge(corrout,outEachExperiment,sort=FALSE,by=1)
#         } # end else
# 

    }else{
        print("Unknown")
    }


}

if(argv$warnings){print(warnings())}


#write.csv(corrout,paste(outputFile,"output.csv",sep='corr-'),row.names=FALSE)
write.csv(compout,paste(outputFile,"output.csv",sep='comp-'),row.names=FALSE)
#write.csv(folchout,paste(outputFile,"output.csv",sep='folch-'),row.names=FALSE)

print("Finished performing the requested analyses.")

print("You may want to run merge-comp-correl-results-files.R next.")


























