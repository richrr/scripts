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
   print(dim(subdata))
   write.csv (subdata,"testout2.txt")
   

   outForPUC = subdata
   
   # calculate combined Pvalue for interest
   PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
   #print(PvalueColnames)
   PvalueColnames = PvalueColnames[grep(search_group,PvalueColnames)]
   #print(PvalueColnames)
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
   result = outForPUC 



#calculate FDR for combined pvalue
combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
result = cbind(result,combinedFDR)

   head(result)
   quit()


############## to do:  calculate PUC #######################



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

