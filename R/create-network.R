library(argparser)
library(stringr)


### for (i) FC and (ii) correlation, we can use any dataset (of the multiple datasets) since we
### (i) calc. correlations only on consistent genes and (ii) calc PUC only on consistent pairs.
### using median values "USES" all datasets.


#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Extract pairs which pass the FDR and PUC criteria across datasets for a single correlation (group)')
p <- add_argument(p, "--file", help="file which has merged correlations from multiple expts", nargs=1) # required; but written as optional format so I explicitly mention the "--file"
p <- add_argument(p, "--consistent", help="file which has list of consistent elements (mostly pairs)", nargs=1) # required; but written as optional format so I explicitly mention the "--consistent"
p <- add_argument(p, "--group", help="use the correlation for 'group'. Generate the network for this group") # required; but written as optional format so I explicitly mention the arg
#p <- add_argument(p, "--groups", help="use the correlations for 'groups'. Calc. PUC using both these groups (states). Default only uses one state as mentioned in --group", nargs=2) # not implemented
p <- add_argument(p, "--foldchange", help="file containing foldchange to be used in analysis") ### use the merged file with consistent genes across datasets
p <- add_argument(p, "--output", help="output file", default="./build_netw_")
p <- add_argument(p, "--indivPvalCutoff", help="individualPvalueCutoff", default=0.005, type="numeric") # 0.05
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.1, type="numeric") # 0.05
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.3, type="numeric") # 

p <- add_argument(p, "--foldchMean", help="use fold change from mean", flag=TRUE)  # default fold change is median
#p <- add_argument(p, "--transposeOutput", help="tranpose the results so the rows are analysis and columns are gene or pairs", flag=TRUE)  # easier for downstream grep of required analysis, do not use when you expect many genes or pairs since it might truncate when you open in excel or librecalc due ot limited columns

p <- add_argument(p, "--noPUC", help="do not calc. PUC, so you do not need to use fold change", flag=TRUE) 

p <- add_argument(p, "--analysisfc", help="Partial header (Analysis number) to be used for median fold change calc", default="ALL") #e.g. "Analys 1 " ; avoids cutting the required columns as input for the script  # make sure there is space after the number to avoid selecting analysis 11, etc.

p <- add_argument(p, "--analysiscorr", help="Partial header (Analysis number) to be used for selecting the columns for combined pval and median correlations", default="ALL") #e.g. "Analys 1 " ; avoids cutting the required columns as input for the script   # make sure there is space after the number to avoid selecting analysis 11, etc.


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


foldchVar = 'FolChMedian'
if(argv$foldchMean){foldchVar = "FoldChange"}
outputFile = paste(outputFile , foldchVar, '_', sep='')

noPUC = argv$noPUC
analysisfc = argv$analysisfc
analysiscorr = argv$analysiscorr

if(analysisfc != "ALL" || analysiscorr != "ALL"){
  outputFile = paste(outputFile , "FC", analysisfc, "Corr", analysiscorr, '_', sep='')
  outputFile = gsub(' ', '_', outputFile)
}

networkFile = paste(c(outputFile, search_group, "indiv-pval", individualPvalueCutoff ,"comb-pval", combinedPvalueCutoff, "comb-fdr", combinedFDRCutoff, ".csv"), collapse='_')


#------------------------------------------------------------------------------------------------
# check whether the first "number of files" columns are the same and remove them after setting row names
#------------------------------------------------------------------------------------------------
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


#------------------------------------------------------------------------------------------------
# calc PUC (%) at each FDR threshold and plot
#------------------------------------------------------------------------------------------------
calc_PUC_at_thresholds = function(df, str='all'){
    
    write.csv(df, paste("tmp", str, "csv", sep='.'))
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
    
    plt = cbind(combined_fdr, puc_percent)
    
    pdf(paste(search_group, str, 'FDRvsPUC.pdf', sep='-'))
    plot(combined_fdr, puc_percent, type="o", ylab="PUC(%)", xlab="combined FDR", pch=10, cex=.2, ylim=c(0, 100) ) 
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
   data = read.csv( argv$file, header=TRUE, check.names=FALSE)
   subdata = data[data[,1] %in% consist_elems,] 
   total_numb_input_files = length(grep("pairName", colnames(subdata)))
   subdata = remove_redundant_columns(subdata, total_numb_input_files)
   

   outForPUC = subdata
   
   CorrAnalysColnames = colnames(outForPUC)
   # only keep the requested analysis number
   if(analysiscorr != "ALL"){
       CorrAnalysColnames = CorrAnalysColnames[grep(analysiscorr,CorrAnalysColnames)]
   }
   print(CorrAnalysColnames)
   outForPUC = outForPUC[,CorrAnalysColnames]
   
   # calculate combined Pvalue for interest group
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
   outForPUC = cbind(outForPUC,combinedPvalue)

   # calculate median coefficient for interest
   CoefficientColnames = colnames(outForPUC)[grep("Coefficient",colnames(outForPUC))]
   CoefficientColnames = CoefficientColnames[grep(search_group,CoefficientColnames)]
   interestedCoefficientData = outForPUC[,CoefficientColnames]
   interestedCoefficientData = as.matrix(interestedCoefficientData)
   interestedCoefficientData = apply(interestedCoefficientData,2,function(x){as.numeric(as.vector(x))})
   combinedCoefficient = apply(interestedCoefficientData,1, function(x){round(median(x, na.rm = TRUE), 3)})
   outForPUC = cbind(outForPUC,combinedCoefficient)
   
   result = outForPUC 
   #calculate FDR for combined pvalue
   combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
   result = cbind(result,combinedFDR)
   
   write.csv (result,paste( outputFile, search_group, "tmp-out.csv", sep='-'))
   
 
   # calculate median FoldChange
   FoldChangeMetabolic =  ''
   if(noPUC){
           # no need to loop up the fold change
        } else {
   FoldChangeMetabolic =  read.csv(FoldChangeFile,header = TRUE,check.names=FALSE)
   FoldChangeMetabolic = remove_redundant_columns(FoldChangeMetabolic, total_numb_input_files)
   FoldChangeColnames = colnames(FoldChangeMetabolic)[grep(foldchVar,colnames(FoldChangeMetabolic))]
   # if an Analysis number is specified then select only that analysis number
   if(analysisfc != "ALL"){
       FoldChangeColnames = FoldChangeColnames[grep(analysisfc,FoldChangeColnames)]
   }
   #print(FoldChangeColnames)

   interestedFoldChangeData = FoldChangeMetabolic[,FoldChangeColnames]
   interestedFoldChangeData = as.matrix(interestedFoldChangeData)
   interestedFoldChangeData = apply(interestedFoldChangeData,2,function(x){as.numeric(as.vector(x))})
   combinedFoldChange = apply(interestedFoldChangeData,1, function(x){round(median(x, na.rm = TRUE), 3)})
   FoldChangeMetabolic$geneName = rownames(FoldChangeMetabolic)
   FoldChangeMetabolic = cbind(FoldChangeMetabolic,combinedFoldChange)
   }

###########################################################################################################
#calc PUC
###########################################################################################################

   Data = result
#--------------------------------------------------------------------------------
# calculate PUC
#--------------------------------------------------------------------------------
forPUC = function(FoldChangeMetabolic,noPUC){

	#change out format to partner1 partner2 for PUC
	row_names_Data = rownames(Data)
	pair = str_split( row_names_Data ,"<==>") 
	pairs = t(as.data.frame(pair))
	
	colnames(pairs) = c("partner1","partner2")
	Data = apply(Data, 2, function(x) as.numeric(as.character(x))) # convert the chars to numeric
	rownames(pairs) = row_names_Data # remove this if you do not want row names
        outForPUC = cbind(pairs,Data)
	
        grep_cols_c = grep("pvalue", colnames(outForPUC), ignore.case = TRUE , value=TRUE, fixed=TRUE) 
        grep_cols_c = append(grep_cols_c, grep("combined" , colnames(outForPUC), value=TRUE, fixed=TRUE) )
 
        g_grep_cols = c("partner1","partner2", grep_cols_c) 
           
        outForPUC = outForPUC[,g_grep_cols]
        
        if(noPUC){
           # no need to loop up the fold change
        } else {
        # attach the foldChange information for each partner
        FoldChangeCol = grep("combined" , colnames(FoldChangeMetabolic), value=TRUE, fixed=TRUE)
	FoldMetab1_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner1"]), c("geneName", FoldChangeCol)]
	colnames(FoldMetab1_InPair) = c("partner1InFold","partner1_FoldChange")
	FoldMetab2_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner2"]),c("geneName", FoldChangeCol)]
	colnames(FoldMetab2_InPair) = c("partner2InFold","partner2_FoldChange")

	outForPUC = cbind(outForPUC,FoldMetab1_InPair,FoldMetab2_InPair)
	}
	
        # calculate correlation Direction For combined correlation coefficient of interest 
        # at this point we only have the consistent pairs left, so the value of combined corr coeff is ok to use
        interestedCoefficientColnames = grep("Coefficient",colnames(outForPUC), value=TRUE, fixed=TRUE)     
	print(interestedCoefficientColnames)
	interestedCorrelationData = outForPUC[,interestedCoefficientColnames, drop=FALSE]
	interestedCorrelationData = apply(interestedCorrelationData,2,function(x){as.numeric(as.vector(x))})
	signOfInterestedCorrelationData = interestedCorrelationData/abs(interestedCorrelationData)
	rownames(signOfInterestedCorrelationData) = c()
	colnames(signOfInterestedCorrelationData) = paste(colnames(interestedCorrelationData),"correlationDirection",sep=".")
	matchedExpressionDirection = signOfInterestedCorrelationData
	
	
        if(noPUC){
           # no need to loop up the fold change direction
           outForPUC = cbind(outForPUC,signOfInterestedCorrelationData)
        } else {
	# calculate fold change direction for each partner
	FoldChangeColnames = colnames(outForPUC)[grep("FoldChange",colnames(outForPUC))] # since this is using the combined fold change calculated above, you do not need the foldchVar variable
	FoldChangeData = outForPUC[,FoldChangeColnames]
	FoldChangeDirection = (FoldChangeData-1)/abs(FoldChangeData-1)
	names(FoldChangeDirection) = c()
	colnames(FoldChangeDirection) = paste(colnames(FoldChangeData),"FoldChangeDirection",sep=".")
	
	# calculate if fold change direction are the same for the two partners
	IfFoldChangeDirectionMatch = apply(FoldChangeDirection,1,prod)
	names(IfFoldChangeDirectionMatch) = c()
	colnames(matchedExpressionDirection) = c()
	
	# use "matchedExpressionDirection" and "IfFoldChangeDirectionMatch" to calc PUC, i.e. if these two are the same PUC=1 (good)
	PUC = IfFoldChangeDirectionMatch * matchedExpressionDirection 
	outForPUC = cbind(outForPUC,signOfInterestedCorrelationData,FoldChangeDirection,IfFoldChangeDirectionMatch,PUC)
        }
        
	return(outForPUC)
		
} 

    PUCoutfile = ''
    if(noPUC){
           PUCoutfile = paste(outputFile, search_group, "noPUC-output.csv",sep='-')
    } else {
         PUCoutfile = paste(outputFile, search_group, "PUC-output.csv",sep='-')
    }

    if(file.exists(PUCoutfile))
    {
        print("Deleted old PUC output")
        file.remove(PUCoutfile)
    }

    result = forPUC(FoldChangeMetabolic,noPUC)
    data = result
    
    if(noPUC){
           # no need to calc stats
    } else {

    # calculate percentage of unexpected correlations in the entire file
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

    # sort as per combined fdr
    sorted_result = result[order(result$"combinedFDR"), ]
    plt = calc_PUC_at_thresholds(sorted_result) # walk along fdr and calc puc at diff fdr
    write.csv(plt, paste(PUCoutfile, "all-edges.csv", sep='.'), row.names=FALSE)

    pos_plt = calc_PUC_at_thresholds(sorted_result[which(sorted_result$combinedCoefficient.correlationDirection == 1),], "pos") # walk along fdr and calc puc at diff fdr
    write.csv(pos_plt, paste(PUCoutfile, "pos-edges.csv", sep='.'), row.names=FALSE)

    neg_plt = calc_PUC_at_thresholds(sorted_result[which(sorted_result$combinedCoefficient.correlationDirection == -1),], "neg") # walk along fdr and calc puc at diff fdr
    write.csv(neg_plt, paste(PUCoutfile, "neg-edges.csv", sep='.'), row.names=FALSE)

    data = sorted_result
    }



###########################################################################################################
#generate network
###########################################################################################################


#----------------------------------------------------------------------
# calc stats of the generated network
#----------------------------------------------------------------------
calc_stats = function(inNet, correlThreshold=0){

        out = file(paste(networkFile,'-stats-log.txt',sep='') ,'a')
        
        combcoeffcol = as.vector(inNet[,"combinedCoefficient"])
        # not sure why it behaves this way!
        lowest_neg_corrl = min(combcoeffcol[combcoeffcol < 0])
        highest_neg_corrl = max(combcoeffcol[combcoeffcol < 0])

        lowest_pos_corrl = min(combcoeffcol[combcoeffcol > 0])
        highest_pos_corrl = max(combcoeffcol[combcoeffcol > 0])

        
        coefficientData = inNet[,grep("correlationDirection",colnames(inNet)),drop=FALSE]
	coefficientData = apply(coefficientData,2,function(x){as.numeric(as.vector(x))})
        res_pos = apply(coefficientData, 2, function(x) sum(x == 1))
        res_neg = apply(coefficientData, 2, function(x) sum(x == -1))
        res = cbind(res_pos, res_neg)
        #print(res)
        
        #print(paste(c("Total number of edges: ", nrow(inNet)), collapse="")) 

        setpartner1 = unique(as.vector(inNet[,"partner1"]))
        setpartner2 = unique(as.vector(inNet[,"partner2"]))

        nodes = union(setpartner1, setpartner2)
        #print(paste(c("Number of unique nodes: ", length(nodes)), collapse=""))

        edgesDistinctNodes = inNet[which(inNet[,"partner1"]!=inNet[,"partner2"]), ] 
        #print(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse=""))

        write.csv(paste(c("Lowest pos corr:" , lowest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Highest pos corr:" , highest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Lowest neg corr:" , lowest_neg_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Highest neg corr:" , highest_neg_corrl), collapse='') , file=out, row.names=FALSE)

        write.csv(paste(c("Pos corr edges:" , toString(res_pos)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Neg corr edges:" , toString(res_neg)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Ratio of Pos to Neg corr edges:" , toString(res_pos/res_neg)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Total number of edges: ", nrow(inNet)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of unique nodes: ", length(nodes)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Ratio of edges to nodes: ", nrow(inNet)/length(nodes)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse='') , file=out, row.names=FALSE)
        close(out)

}


#----------------------------------------------------------------------
# for a given combinedPvalue and combinedfdr threshold, generate network
#----------------------------------------------------------------------
generateNetwork = function(){
       #generate network After Filter

       out = ''
       if(noPUC){
	    # at this point you only have pairs with consistent correlations 
	    out = data
	    #out = out[!is.na(out$"PUC"),] # remove the rows with 'NA' in PUC columns  # decide whether to use this combinedCoefficient.correlationDirection
        } else {
	    # find PUC expected
	    out = data[data[,"PUC"]==1 ,]
	    out = out[!is.na(out$"PUC"),] # remove the rows with 'NA' in PUC columns
        }

	out = as.matrix(out)
        #df = apply(df[,c("combinedFDR", "PUC")],2,function(x){as.numeric(as.vector(x))})

	
	out = out[(as.numeric(out[,"combinedPvalue"])<combinedPvalueCutoff)==1,]
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
        
	write.csv (outNetwork,networkFile, quote=FALSE)

        calc_stats(outNetwork)

        print("Done!")
}


    generateNetwork()

