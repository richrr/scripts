library(argparser)
#library(hash)
#library(psych)
#library(gtools)
#library(corpcor)
#library(reshape)
library(stringr)

# Rscript /nfs1/Morgun_Lab/richrr/scripts/R/generate-network.R --file testPUC.csv --combine 1-2 --group KO --indivPvalCutoff 0.05 --combPvalCutoff 0.1 --combFDRCutoff 0.3
# Rscript /nfs1/Morgun_Lab/richrr/scripts/R/generate-network.R --file testPUC.csv --combine 1 2 3 4 --group KO

#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Extract pairs which pass the FDR and PUC criteria across datasets for a single correlation (group)')
p <- add_argument(p, "--file", help="file which has PUC information", nargs=1) # required; but written as optional format so I explicitly mention the "--file"
p <- add_argument(p, "--consistent", help="file which has list of consistent elements (mostly pairs)", nargs=1) # required; but written as optional format so I explicitly mention the "--consistent"

p <- add_argument(p, "--group", help="use the correlation for 'group'. Generate the network for this group") # required; but written as optional format so I explicitly mention the arg
# OR
p <- add_argument(p, "--groups", help="use the correlations for 'groups'. Calc. PUC using both these groups (states). Default only uses one state as mentioned in --group", nargs=2) # required; but written as optional format so I explicitly mention the arg

p <- add_argument(p, "--foldchange", help="file containing foldchange to be used for correlation")

p <- add_argument(p, "--output", help="output file", default="./build_netw_")
p <- add_argument(p, "--combine", help="combine different expts. before generating networks", nargs=Inf)
  # if you do not provide the combine option, it defaults to only using the first analysis in PUC
  # first analysis (contains correl in A, correl in B for each gene pair; fold change in A/B for each gene; PUC result) 
  # the above ends in 
  # "x
  # Percentage of Unexpected Correlation.....%"
  # the next analysis has the same setup for something else
  
  # the arguments can be analysis numbers (delimited by the above breakpoint string), NOT the "Analys " in the column header
  # e.g. 1-5 or 1 2 3 4 5 or 1:5
  

p <- add_argument(p, "--indivPvalCutoff", help="individualPvalueCutoff", default=0.005, type="numeric") # 0.05
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.1, type="numeric") # 0.05
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.3, type="numeric") # 


argv <- parse_args(p)
#print (p)


if(length(argv$file) < 1)
{
  print("At least 2 files are required. Give --file file ... in cmd")
  quit()
}

outputFile = argv$output

individualPvalueCutoff = argv$indivPvalCutoff
combinedPvalueCutoff = argv$combPvalCutoff
combinedFDRCutoff = argv$combFDRCutoff
FoldChangeFile = argv$foldchange
search_group = argv$group

networkFile = paste(outputFile, search_group, "indiv-pval", individualPvalueCutoff ,"comb-pval", combinedPvalueCutoff, "comb-fdr", combinedFDRCutoff, ".csv", collapse='_')




data = read.csv( argv$file, header=TRUE, check.names=FALSE, row.names=1)

   consist_elems = read.csv( argv$consistent, header=FALSE)

   consist_elems = as.vector(consist_elems[!duplicated(consist_elems[,1]),1] )	#delete duplicated array probes
   consist_elems = consist_elems[order(consist_elems)] # ascending sort
   consist_elems = grep("<==>", consist_elems, value=TRUE)
   #head(consist_elems)
   
   subdata = data[consist_elems,-c(1)] 
   
   #head(subdata)
   
   #print(head(subdata)[1:5])
   #print(dim(subdata))

   outForPUC = subdata
   
   # calculate combined Pvalue for interest
   PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
   PvalueColnames = PvalueColnames[grep(search_group,PvalueColnames)]
   interestedPvalueData = outForPUC[,PvalueColnames]
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

   # calculate median coefficient for interest
   CoefficientColnames = colnames(outForPUC)[grep("Coefficient",colnames(outForPUC))]
   CoefficientColnames = CoefficientColnames[grep(search_group,CoefficientColnames)]
   interestedCoefficientData = outForPUC[,CoefficientColnames]
   interestedCoefficientData = as.matrix(interestedCoefficientData)
   interestedCoefficientData = apply(interestedCoefficientData,2,function(x){as.numeric(as.vector(x))})
   combinedCoefficient = apply(interestedCoefficientData,1, function(x){median(x, na.rm = TRUE)})
   outForPUC = cbind(outForPUC,combinedCoefficient)



   result = outForPUC 


#calculate FDR for combined pvalue
combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
result = cbind(result,combinedFDR)

   #head(result)
   
   write.csv (result,"testout2.txt")
   

  #sorted_result = result[with(result, order(combinedFDR)), ] 
 



############## to do:  calculate PUC #######################

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

####################### XD code below ##############################
### is this your way of checking consistency between three datasets?
if(FALSE){
	# see if direction match in interest
	ifMatch = apply(signOfInterestedCorrelationData,1
					,function(tt){
								uniqueValues = unique(tt[!is.na(tt)])
								NumberOfUniqueValues = length(uniqueValues)
								if (length(tt[!is.na(tt)])<2){
									return (uniqueValues[1])
								}
								if (NumberOfUniqueValues>1){
									return (0)
								}else if(NumberOfUniqueValues==0){
									return (NA)
								}else{
									return (1*uniqueValues[1])
								}
							}
					)	# if each row of signOfInterestedCorrelationData is like -1,1,-1; if match, "ifMatch" will be 3 or -3

}
###############################################################




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


### at this point you have multiple datasets, which dataset to use:
### for (i) FC and (ii) correlation, we can use any dataset since we (i) calc. correlations only on consistent genes and (ii) calc PUC only on consistent pairs 
### think about using median values

### how to get the FC values based on the analysis file?
  ### use the merged file from the check consistent genes

### use single state to calc network. Allow option to use both states if requested to identify regulatory genes?



#result = forPUC(FoldChangeFile, categ1, categ2)



#for(indx in 1:nrow(AnalysToDoList)){
if(FALSE){
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



###########################################################################################################
#generate network
###########################################################################################################
data = result
#print(head(data))


calc_stats = function(inNet, correlThreshold=0){

        # if it is a single dataset it becomes a vector so adding the drop=FALSE
        # get the correlation coeff for the group of interest
        coefficientData = inNet[,grep("Coefficient",colnames(inNet)),drop=FALSE]
        coefficientData = coefficientData[,grep("correlationDirection",colnames(coefficientData), invert=T),drop=FALSE]  # remove this column
	coefficientData = coefficientData[,grep(search_group,colnames(coefficientData)),drop=FALSE] # keep group of interest
	        
        #print(head(coefficientData))
        res_pos = apply(coefficientData, 2, function(x) sum(x > correlThreshold))
        res_neg = apply(coefficientData, 2, function(x) sum(x < correlThreshold))
        res = cbind(res_pos, res_neg)
        print(res)
        
        #print(head(inNet))

        print(paste("Total number of edges: ", nrow(inNet), collapse="")) 

        setpartner1 = unique(as.vector(inNet[,"1"]))
        #print(length(setpartner1))
        
        setpartner2 = unique(as.vector(inNet[,"2"]))
        #print(length(setpartner2))

        nodes = union(setpartner1, setpartner2)
        print(paste("Number of unique nodes: ", length(nodes), collapse=""))


        edgesDistinctNodes = inNet[which(inNet$"1"!=inNet$"2"), ] 
        #print(head(edgesDistinctNodes))
        print(paste("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes), collapse=""))

}

generateNetwork = function(){
	#generate network After Filter

	# find PUC expected and significant in combinedPvalue and significant in combinedFDR
	out = data[data[,"PUC"]==1 ,]
	out = out[!is.na(out$"PUC"),] # remove the rows with 'NA' in PUC columns
	#print(head(out))
	#print(dim(out))
	out = out[(out[,"combinedFDR"]<combinedFDRCutoff)==1,]
	# find significant in individual pvalue: all of the pvalues for all of the datasets for each pair must be smaller than threshold
	# find pvalue data

        # if it is a single dataset it becomes a vector so adding the drop=FALSE
	pvalueData = out[,grep("pvalue",colnames(out)), drop=FALSE]
	pvalueData = pvalueData[,grep(search_group,colnames(pvalueData)), drop=FALSE]
	
	pvalueData = as.matrix(pvalueData)
	# calculate the largest pvalue among all datasets for each gene, this smallest pvalue must be smaller than threshold
	passIndevidualPvalue = apply(pvalueData,1,max)<individualPvalueCutoff
	outNetwork = out[passIndevidualPvalue,]

        outNetwork$names<-rownames(outNetwork)
        splt = str_split_fixed(outNetwork$names, "<-->", 2)
        outNetwork = cbind(outNetwork, splt)
        
        # select the group of interest
        coefficientData = outNetwork[,grep("Coefficient",colnames(outNetwork)), drop=FALSE]
        coefficientData = coefficientData[,grep("correlationDirection",colnames(coefficientData), invert=T), drop=FALSE]  # remove this column
	coefficientData = coefficientData[,grep(search_group,colnames(coefficientData)), drop=FALSE] # keep group of interest
	
	# check consistency only if more than one dataset
	print(ncol(coefficientData))
	if(ncol(coefficientData) > 1){
            res = checkConsistency_across_expts(coefficientData, "Coefficient", total_numb_input_files)
            #print(res)

            resvec = as.vector(unlist(res))
            #print(resvec)
            outNetwork = outNetwork[resvec,]
        }
	write.csv (outNetwork,networkFile)

        calc_stats(outNetwork)

}
###########################################################################################################


generateNetwork()

