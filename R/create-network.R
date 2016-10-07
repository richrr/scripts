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



 # check whether the first "number of files" columns are the same and remove them after setting row names
remove_redundant_columns = function(this_df, total_numb_input_files){
         l_identical_inp = as.character(this_df[, 1])
        
         
         # if there are more than 1 datasets
         if(total_numb_input_files>1){
             # create a data frame of all datasets
             for(z in 2:total_numb_input_files){
                    l_identical_inp = cbind(l_identical_inp, as.character(this_df[,z]))
             }

             if(all(apply(l_identical_inp, 2, identical, l_identical_inp[, 1]))){
                    print("All are equal")
                    rownames(this_df) <- this_df[, 1]                           ## set column 1 as rownames
                    this_df <- this_df[, -c(1:total_numb_input_files)]          ## remove the first "total numb of input files" columns      
             } else {
                    print(l_identical_inp)
                    print("All are NOT equal")
                    quit()
             }
          }

      return(this_df)
}


calc_PUC_at_thresholds = function(df){
    
    combined_fdr = c(0)
    puc_percent = c(0)
    df = as.matrix(df)
    df = apply(df[,c("combinedFDR", "PUC")],2,function(x){as.numeric(as.vector(x))})
    #fdr_range = range(as.vector(df[,"combinedFDR"]))
    fdr_min = min(as.vector(df[,"combinedFDR"]))
    fdr_max = max(as.vector(df[,"combinedFDR"]))
    print(fdr_min)
    print(fdr_max)
    m_df = as.matrix(df[,c("combinedFDR", "PUC")])
    m_df = apply(m_df,2,function(x){as.numeric(as.vector(x))})
    
    for(r in sort(seq(as.numeric(fdr_min), as.numeric(fdr_max), length.out=100))){  ##length.out: desired length of the sequence
        combined_fdr = c(combined_fdr, r)
        
        #nrow(m_df[which(m_df[,"combinedFDR"] <= r & m_df[,"PUC"] == 1), , drop=FALSE])
        df_under_fdr = m_df[which(m_df[,"combinedFDR"] <= r), , drop=FALSE]
        DEN = nrow(df_under_fdr)
        
        goog_puc_in_df_under_fdr =  df_under_fdr[which(df_under_fdr[,"PUC"] == 1), , drop=FALSE]
        NUM = DEN - nrow(goog_puc_in_df_under_fdr)
        
        puc_percent = c(puc_percent, as.numeric(NUM*100/DEN))
    }
    #print(head(combined_fdr))
    #print(head(puc_percent))
    
    plt = cbind(combined_fdr, puc_percent)
    #print(head(plt))
    
    pdf(paste(search_group, 'FDRvsPUC.pdf', sep='-'))
    plot(combined_fdr, puc_percent, type="o", ylab="PUC(%)", xlab="combined FDR", pch=10, cex=.2 ) 
    #plot(plt, type="o", ylab="PUC(%)", xlab="combined FDR", pch=10, cex=.2 ) 
    
    dev.off()

    return(plt)

}


   # consistent pairs
   consist_elems = read.csv( argv$consistent, header=FALSE)
   consist_elems = as.vector(consist_elems[!duplicated(consist_elems[,1]),1] )	#delete duplicated array probes
   consist_elems = consist_elems[order(consist_elems)] # ascending sort
   consist_elems = grep("<==>", consist_elems, value=TRUE) # keep pairs
   
   # keep consistent pairs only
   #data = read.csv( argv$file, header=TRUE, check.names=FALSE, row.names=1)
   #subdata = data[consist_elems,-c(1)] 
   data = read.csv( argv$file, header=TRUE, check.names=FALSE)
   subdata = data[data[,1] %in% consist_elems,] 
   total_numb_input_files = length(grep("pairName", colnames(subdata)))
   print(total_numb_input_files)
   print(dim(subdata))

   subdata = remove_redundant_columns(subdata, total_numb_input_files)
   print(dim(subdata))
   
   #head(subdata)
   
   #print(head(subdata)[1:5])
   #print(dim(subdata))

   outForPUC = subdata
   
   # calculate combined Pvalue for interest
   PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
   PvalueColnames = PvalueColnames[grep(search_group,PvalueColnames)]
   total_numb_input_files = length(PvalueColnames)
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
   #colnames(combinedPvalue) = paste("combinedPvalue", search_group, sep='_')
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
 

#==================================================================================================================
#    						generate foldChange


	FoldChangeMetabolic =  read.csv(FoldChangeFile,header = TRUE,check.names=FALSE)
	#print(head(FoldChangeMetabolic))

        FoldChangeMetabolic = remove_redundant_columns(FoldChangeMetabolic, total_numb_input_files)
        print(dim(FoldChangeMetabolic))
	
        #head(FoldChangeMetabolic)

   # calculate median FoldChange
   FoldChangeColnames = colnames(FoldChangeMetabolic)[grep("FoldChange",colnames(FoldChangeMetabolic))]
   #FoldChangeColnames = FoldChangeColnames[grep(search_group,FoldChangeColnames)] # you may not need this
   interestedFoldChangeData = FoldChangeMetabolic[,FoldChangeColnames]
   interestedFoldChangeData = as.matrix(interestedFoldChangeData)
   interestedFoldChangeData = apply(interestedFoldChangeData,2,function(x){as.numeric(as.vector(x))})
   combinedFoldChange = apply(interestedFoldChangeData,1, function(x){median(x, na.rm = TRUE)})
   FoldChangeMetabolic$geneName = rownames(FoldChangeMetabolic)
   FoldChangeMetabolic = cbind(FoldChangeMetabolic,combinedFoldChange)

   #head(FoldChangeMetabolic)





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


#Data = read.csv("testout2.txt",header = TRUE,check.names=FALSE,stringsAsFactors = FALSE)
Data = result


#head(Data)

# calculate PUC
#==================================================================================================================
forPUC = function(FoldChangeMetabolic){
	#FoldChangeCol = paste(categ1,"vs",categ2,"FoldChange",sep=" ")
        #geneIDCol = "geneName"

	#change out format for PUC
	row_names_Data = rownames(Data)
	pair = str_split( row_names_Data ,"<==>") 
	pairs = t(as.data.frame(pair))
	
	colnames(pairs) = c("partner1","partner2")
	Data = apply(Data, 2, function(x) as.numeric(as.character(x))) # convert the chars to numeric
	rownames(pairs) = row_names_Data # remove this if you do not want row names
        outForPUC = cbind(pairs,Data)
	#print(head(outForPUC))
	
        grep_cols_c = grep("pvalue", colnames(outForPUC), ignore.case = TRUE , value=TRUE, fixed=TRUE) 
        grep_cols_c = append(grep_cols_c, grep("combined" , colnames(outForPUC), value=TRUE, fixed=TRUE) )
 
        g_grep_cols = c("partner1","partner2", grep_cols_c) 
           

        #print(g_grep_cols)

        outForPUC = outForPUC[,g_grep_cols]
        
#==================================================================================================================
#    						generate foldChange

        FoldChangeCol = grep("combined" , colnames(FoldChangeMetabolic), value=TRUE, fixed=TRUE)
	FoldMetab1_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner1"]), c("geneName", FoldChangeCol)]
	colnames(FoldMetab1_InPair) = c("partner1InFold","partner1_FoldChange")
	
	FoldMetab2_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner2"]),c("geneName", FoldChangeCol)]
	colnames(FoldMetab2_InPair) = c("partner2InFold","partner2_FoldChange")

	outForPUC = cbind(outForPUC,FoldMetab1_InPair,FoldMetab2_InPair)
	
	
#==================================================================================================================
#    					calculate PUC

        
        # calculate correlation Direction For correlation coefficient of interest
	interestedCoefficientColnames = grep("Coefficient",colnames(outForPUC), value=TRUE, fixed=TRUE)     
	print(interestedCoefficientColnames)
	interestedCorrelationData = outForPUC[,interestedCoefficientColnames, drop=FALSE]
	interestedCorrelationData = apply(interestedCorrelationData,2,function(x){as.numeric(as.vector(x))})
	signOfInterestedCorrelationData = interestedCorrelationData/abs(interestedCorrelationData)					# forOutput
	rownames(signOfInterestedCorrelationData) = c()
	colnames(signOfInterestedCorrelationData) = paste(colnames(interestedCorrelationData),"correlationDirection",sep=".")
	matchedExpressionDirection = signOfInterestedCorrelationData
	
	# calculate fold change direction	
	FoldChangeColnames = colnames(outForPUC)[grep("FoldChange",colnames(outForPUC))]
	FoldChangeData = outForPUC[,FoldChangeColnames]
	FoldChangeDirection = (FoldChangeData-1)/abs(FoldChangeData-1)
	names(FoldChangeDirection) = c()
	
	colnames(FoldChangeDirection) = paste(colnames(FoldChangeData),"FoldChangeDirection",sep=".")
		
	# calculate if fold change direction are the same for the two partners
	IfFoldChangeDirectionMatch = apply(FoldChangeDirection,1,prod)
	names(IfFoldChangeDirectionMatch) = c()
	colnames(matchedExpressionDirection) = c()
	
	# use "signOfInterestedCorrelationData" and "if fold change direction are the same" to calculate PUC
	PUC = IfFoldChangeDirectionMatch * matchedExpressionDirection    				## forOutput	
	outForPUC = cbind(outForPUC,signOfInterestedCorrelationData,FoldChangeDirection,IfFoldChangeDirectionMatch,PUC)
	#print(head(outForPUC))

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
### using median values

### how to get the FC values based on the analysis file?
  ### use the merged file from the check consistent genes

### use single state to calc network. Allow option to use both states if requested to identify regulatory genes?

### add proper labels to the combined p value, coeff, to keep track of which is the group of interest

result = forPUC(FoldChangeMetabolic)

#print(head(result))


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

print("Done with PUC")

sorted_result = result[order(result$"combinedFDR"), ]
print(head(sorted_result))

  plt = calc_PUC_at_thresholds(sorted_result) # walk along fdr and calc puc at diff fdr





###########################################################################################################
#generate network
###########################################################################################################
data = sorted_result
#print(head(data))




calc_stats = function(inNet, correlThreshold=0){

        
        # if it is a single dataset it becomes a vector so adding the drop=FALSE
        # get the correlation coeff for the group of interest
        coefficientData = inNet[,grep("correlationDirection",colnames(inNet)),drop=FALSE]
	
	coefficientData = apply(coefficientData,2,function(x){as.numeric(as.vector(x))})
	
        res_pos = apply(coefficientData, 2, function(x) sum(x == 1))
        res_neg = apply(coefficientData, 2, function(x) sum(x == -1))
        res = cbind(res_pos, res_neg)
        print(res)
        
        print(paste("Total number of edges: ", nrow(inNet), collapse="")) 

        setpartner1 = unique(as.vector(inNet[,"partner1"]))
        #print(length(setpartner1))
        
        setpartner2 = unique(as.vector(inNet[,"partner2"]))
        #print(length(setpartner2))

        nodes = union(setpartner1, setpartner2)
        print(paste("Number of unique nodes: ", length(nodes), collapse=""))


        edgesDistinctNodes = inNet[which(inNet[,"partner1"]!=inNet[,"partner2"]), ] 
        #print(head(edgesDistinctNodes))
        print(paste("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes), collapse=""))

}


### for a given fdr threshold, genr network
generateNetwork = function(){
	#generate network After Filter

	# find PUC expected and significant in combinedPvalue and significant in combinedFDR
	out = data[data[,"PUC"]==1 ,]
	out = out[!is.na(out$"PUC"),] # remove the rows with 'NA' in PUC columns
	
	out = as.matrix(out)
        #df = apply(df[,c("combinedFDR", "PUC")],2,function(x){as.numeric(as.vector(x))})

	
	out = out[(as.numeric(out[,"combinedFDR"])<combinedFDRCutoff)==1,]
	# find significant in individual pvalue: all of the pvalues for all of the datasets for each pair must be smaller than threshold
	# find pvalue data
        # if it is a single dataset it becomes a vector so adding the drop=FALSE
	pvalueData = out[,grep("pvalue",colnames(out)), drop=FALSE]
	pvalueData = pvalueData[,grep(search_group,colnames(pvalueData)), drop=FALSE]
	
	pvalueData = as.matrix(pvalueData)
	pvalueData = apply(pvalueData,2,function(x){as.numeric(as.vector(x))})
	
	# calculate the largest pvalue among all datasets for each gene, this smallest pvalue must be smaller than threshold
	passIndevidualPvalue = apply(pvalueData,1,max)<individualPvalueCutoff
	outNetwork = out[passIndevidualPvalue,]
        
        #outNetwork$names<-rownames(outNetwork)
        #splt = str_split_fixed(outNetwork$names, "<==>", 2)
        #outNetwork = cbind(outNetwork, splt)
        
	write.csv (outNetwork,networkFile)

        calc_stats(outNetwork)

}
###########################################################################################################


generateNetwork()

