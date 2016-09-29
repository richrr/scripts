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


expand_combine_args = function(pattern, str)
{
    range = strsplit(argv$combine, pattern)
    #print(lengths(range))
    if(lengths(range) > 2) { print("Error in --combine argument. More than one delimiter detected") ; quit()}
    analys_vector = seq(as.numeric(range[[1]][1]), as.numeric(range[[1]][2]), by=1)
    return(analys_vector)   
}


# find where the "PUC results" end for each analysis 
indx = c()
fileName <- argv$file
conn <- file(fileName,open="r")
linn <-readLines(conn)
for (i in 1:length(linn)){
   #print(linn[i])
    #if((linn[i] == "\"x\"" && grepl("Percentage of Unexpected Correlation", linn[i+1])) || (linn[i] == "x" && grepl("Percentage of Unexpected Correlation", linn[i+1])))
    if(grepl("Percentage of Unexpected Correlation", linn[i]))
      {
            indx = append(indx, i)
      }
}
close(conn)
linn = ''
print("Index where PUC results end per analysis:")
print(indx)


# re-read as data frame
#data = read.csv( argv$file, header=TRUE, check.names=FALSE)
data = read.csv( argv$file, header=FALSE, check.names=FALSE)

#print(data)
result = ''
total_numb_input_files = 1

if(is.na(argv$combine) && length(argv$combine))
{
    # remove the last two lines from the first analysis and use the remaining
    #breakpt = indx[1]-3  # accounting for one line of header in the above calculation & read with header
    breakpt = indx[1]-2  # accounting for one line of header in the above calculation & read without header
    subdata = data[1:breakpt,]  
    #print(subdata)
    
    names <- paste(subdata[,1], subdata[,2], sep="<-->")
    rownames(subdata) = names
    subdata = subdata[,-c(1,2)] 

    print(dim(subdata))
    colnames(subdata) = as.character(unlist(subdata[1, ])) # the first row will be the header
    subdata = subdata[-1, ] 
    write.csv (subdata,"testout2.txt")
    #quit()
    outForPUC = subdata
   
    # calculate combined Pvalue for interest
    PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
    #print(PvalueColnames)
    PvalueColnames = PvalueColnames[grep(search_group,PvalueColnames)]
    #print(PvalueColnames)
    interestedPvalueData = outForPUC[,PvalueColnames]
    interestedPvalueData = as.matrix(interestedPvalueData)
    interestedPvalueData = apply(interestedPvalueData,2,function(x){as.numeric(as.vector(x))})
    combinedPvalue = interestedPvalueData
    outForPUC = cbind(outForPUC,combinedPvalue)
    result = outForPUC 

        
} else {

   # use some system to combine the results
   analys_vector = argv$combine
   # if range is given 1:5 or 1-5
   if(length(argv$combine) == 1)
   {
       if(grepl('-', argv$combine))
        {
            analys_vector = expand_combine_args('-', argv$combine)
        } else if(grepl(':', argv$combine))
        {
            analys_vector = expand_combine_args(':', argv$combine)
        } else {print("Error in --combine argument. No known delimiter detected. Use '-' or ':'") ; quit()}

   } else { analys_vector = as.numeric(argv$combine) }  # if separated by space (1 2 3 4 5)
   
   print(analys_vector)
   total_numb_input_files = length(analys_vector)
   
   append_numb = 1
   indx_counter = 1
   breakpt = indx[1]-2  # accounting for one line of header in the above calculation & read without header
   subdata = data[1:breakpt,]  
   
   names <- paste(subdata[,1], subdata[,2], sep="<-->")
   rownames(subdata) = names
   subdata = subdata[,-c(1,2)] 
   
   
   while(append_numb < max(analys_vector))
   {
       append_numb = append_numb + 1
       if( append_numb %in% analys_vector){
           # bind
           newstartpt = indx[indx_counter] + 1 # header of next analysis
           indx_counter = indx_counter + 1
           breakpt = indx[indx_counter]-2  # accounting for one line of header in the above calculation & read without header
           print(newstartpt)
           print(breakpt)
           newsubdata = data[newstartpt:breakpt,] 
           names <- paste(newsubdata[,1], newsubdata[,2], sep="<-->")
           rownames(newsubdata) = names 
           newsubdata = newsubdata[,-c(1,2)] 
           #print(dim(newsubdata))
           subdata = cbind(subdata, newsubdata)
           #print(head(subdata)[1:5])
           #print(head(newsubdata)[1:5])
       }
   
   }
   #print(head(subdata)[1:5])
   print(dim(subdata))
   colnames(subdata) = as.character(unlist(subdata[1, ])) # the first row will be the header
   subdata = subdata[-1, ] 
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
}


##### check which measurement ("gene") pairs have the same trend across experiments ############
checkConsistency_across_expts = function(s_df, condition, total_numb_input_files, correlThreshold=0 , pvalThreshold=0.1, foldchThreshold=1){

    
    list_rows_passing_consistency = list()
    
    if(condition == "Coefficient"){
        res_pos = apply(s_df, 1, function(x) sum(x > correlThreshold))
        res_neg = apply(s_df, 1, function(x) sum(x < correlThreshold))
        rows_passing_consistency = c()
        #print("Positive correlation")
        #print(head(res_pos))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_pos[[idx]])) && res_pos[[idx]] >= (total_numb_input_files-1)){               #---------------- (1)
            # temporary hack (2) to remove cases where they show one dir in expt1 & other in expt2
            # once I have three expts, this hack is no longer needed, use (i.e. uncomment) --------(1) and comment (2)
            #if((!is.na(res_pos[[idx]])) && res_pos[[idx]] > (total_numb_input_files-1)){               #---------------- (2)
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[1]] = rows_passing_consistency
        list_rows_passing_consistency$PosCorrel = rows_passing_consistency

        rows_passing_consistency = c()
        #print("Negative correlation")
        #print(head(res_neg))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_neg[[idx]])) && res_neg[[idx]] >= (total_numb_input_files-1)){               #---------------- (1)
            # temporary hack to remove cases where they show one dir in expt1 & other in expt2
            # once I have three expts, this hack is no longer needed, use (i.e. uncomment) --------(1) and comment (2)
            #if((!is.na(res_neg[[idx]])) && res_neg[[idx]] > (total_numb_input_files-1)){               #---------------- (2)
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[2]] = rows_passing_consistency
        list_rows_passing_consistency$NegCorrel = rows_passing_consistency

    } else if(condition == "pvalue"){
        res_pos = apply(s_df, 1, function(x) sum(x < pvalThreshold))
        rows_passing_consistency = c()
        #print("Pass pvalue cutoff")
        #print(head(res_pos))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_pos[[idx]])) && res_pos[[idx]] >= (total_numb_input_files-1)){               #---------------- (1)
            # temporary hack to remove cases where they show one dir in expt1 & other in expt2
            # once I have three expts, this hack is no longer needed, use (i.e. uncomment) --------(1) and comment (2)
            #if((!is.na(res_pos[[idx]])) && res_pos[[idx]] > (total_numb_input_files-1)){                #---------------- (2)
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[1]] = rows_passing_consistency
        list_rows_passing_consistency$PvalThresh = rows_passing_consistency

    } else if(condition == "FoldChange"){
        res_pos = apply(s_df, 1, function(x) sum(x > foldchThreshold))
        res_neg = apply(s_df, 1, function(x) sum(x < foldchThreshold))
        rows_passing_consistency = c()
        #print("Up regulation")
        #print(head(res_pos))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_pos[[idx]])) && res_pos[[idx]] >= (total_numb_input_files-1)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[1]] = rows_passing_consistency
        list_rows_passing_consistency$Upreg = rows_passing_consistency

        rows_passing_consistency = c()
        #print("Down regulation")
        #print(head(res_neg))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_neg[[idx]])) && res_neg[[idx]] >= (total_numb_input_files-1)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[2]] = rows_passing_consistency
        list_rows_passing_consistency$Dwnreg = rows_passing_consistency

    } 
    #print(list_rows_passing_consistency)
    return(list_rows_passing_consistency)
}










#calculate FDR for combined pvalue
combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
result = cbind(result,combinedFDR)

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

