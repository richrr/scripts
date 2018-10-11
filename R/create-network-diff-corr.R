library(argparser)
library(stringr)


### for (i) FC and (ii) correlation, we can use any dataset (of the multiple datasets) since we
### (i) calc. correlations only on consistent genes and (ii) calc PUC only on consistent pairs.
### using median values "USES" all datasets.


#usage:
# cd /nfs3/PHARM/Morgun_Lab/richrr/test_diff_corr/css.quantif_quantile/prelim_result/merged/corr/mw_sp_t1hfhs-more-samples/p1/per_analysis/stool4wk_dietHFHS/ip_0.3_cp_0.1_cfdr_0.15
# Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/create-network-diff-corr.R --file ../../../merged_mic_mw_sp_t1hfhs-more-samples_corr_p1_FolChMedian_merged-parallel-output.csv --group stool4wk\ dietHFHS stool4wk\ dietNCD --consistent ../../analysis-3-consistent-corr-dir.txt --foldchange /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/microbe-pheno-analyses/analysis/microbes/css.quantif_quantile/prelim_result/merged/comp/mw_sp_t1hfhs-more-samples/p0.2/merged_mic_mw_sp_t1hfhs-more-samples_comp_p0.2_FolChMedian_merged-parallel-output.csv --indivPvalCutoff 0.3 --combPvalCutoff 0.1 --analysisfc Analys\ 1\  --analysiscorr Analys\ 4\  --combFDRCutoff 0.15


#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Create network of specified pairs and keep track of PUC for each group. Use this for differential correlation network')
p <- add_argument(p, "--file", help="file which has merged correlations from multiple expts", nargs=1) # required; but written as optional format so I explicitly mention the "--file"
p <- add_argument(p, "--consistent", help="file which has list of consistent elements (mostly pairs)", nargs=1) # required; but written as optional format so I explicitly mention the "--consistent"
p <- add_argument(p, "--group", help="use the correlation for 'group'. Generate the network for this group", nargs=2) # required; but written as optional format so I explicitly mention the args
p <- add_argument(p, "--foldchange", help="file containing foldchange to be used in analysis") ### use the merged file with consistent genes across datasets
p <- add_argument(p, "--output", help="output file", default="./build_netw_")
p <- add_argument(p, "--indivPvalCutoff", help="individualPvalueCutoff", default=0.005, type="numeric") # 0.05
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.1, type="numeric") # 0.05
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.3, type="numeric") # 
p <- add_argument(p, "--foldchthresh", help="fold change threshold", default=0, type="numeric") # default is logged data so check whether greater than 0. use 1 for unlog data
p <- add_argument(p, "--logbase", help="calc log using the base", default=2) # allowed: 0 (no log), 1 (e), 2, 10

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

if(length(argv$group) != 2)
{
  print("Exactly 2 groups are required. Give --group group1 group2 ... in cmd")
  quit()
}

correlThreshold = 0

outputFile = argv$output

individualPvalueCutoff = argv$indivPvalCutoff
combinedPvalueCutoff = argv$combPvalCutoff
combinedFDRCutoff = argv$combFDRCutoff
FoldChangeFile = argv$foldchange

search_group = argv$group # this contains two values.
search_group_santized = gsub('+', '_', search_group, fixed=TRUE)
print(search_group_santized)
search_group_santized = gsub("\\", '', search_group_santized, fixed=TRUE)
print(search_group_santized)


foldchthresh = argv$foldchthresh
logbase = argv$logbase

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

networkFilewcuts = paste(c(outputFile, search_group_santized, "indiv-pval", individualPvalueCutoff ,"comb-pval", combinedPvalueCutoff, "comb-fdr", combinedFDRCutoff, ".csv"), collapse='_')
networkFile = paste(c(outputFile, search_group_santized, ".csv"), collapse='_')



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
# calc combined coeff, pval and fdr per group.
#------------------------------------------------------------------------------------------------
calc_combined_coeff_pval_fdr = function(ldf){
	
   #print("here")
   for(search_gr in search_group){
   

   # calculate combined Pvalue for individual group
   PvalueColnames = colnames(ldf)[grep("CategCorPVal",colnames(ldf))]
   PvalueColnames = PvalueColnames[grep(paste(search_gr,collapse=" vs "),PvalueColnames)]
   #print(PvalueColnames)
   
   total_numb_input_files = length(PvalueColnames)
   interestedPvalueData = ldf[,PvalueColnames]
   interestedPvalueData = as.matrix(interestedPvalueData)
   interestedPvalueData = apply(interestedPvalueData,2,function(x){as.numeric(as.vector(x))})
   combinedCategCorPVal = apply(interestedPvalueData,1
							,function(pvalues){
										pvalues = pvalues[!is.na(pvalues)]
										statistics = -2*log(prod(pvalues))
										degreeOfFreedom = 2*length(pvalues)
										combined = 1-pchisq(statistics,degreeOfFreedom)
									}
							)
   ldf = cbind(ldf,combinedCategCorPVal)
   # add the group name to the header
   colnames(ldf)[colnames(ldf)=="combinedCategCorPVal"] <- paste("combCategCorPVal", search_gr, sep='_')
   #print(head(ldf))
   

   # calculate median coefficient
   CoefficientColnames = colnames(ldf)[grep("Coefficient",colnames(ldf))]
   CoefficientColnames = CoefficientColnames[grep(search_gr,CoefficientColnames)]
   #print(CoefficientColnames)
   
   interestedCoefficientData = ldf[,CoefficientColnames]
   interestedCoefficientData = as.matrix(interestedCoefficientData)
   interestedCoefficientData = apply(interestedCoefficientData,2,function(x){as.numeric(as.vector(x))})
   combinedCoefficient = apply(interestedCoefficientData,1, function(x){round(median(x, na.rm = TRUE), 3)})
   
   ##########
   # find the pairs which are consistent across cohorts per group
   #   the pairs can be inconsistent in cohorts even though they have the same direction of differential correlation between cohorts:
		#		cohort		CoeffgroupA		CoeffgroupB		CoeffgroupA-CoeffgroupB
		#		c1				0.4				0.1				0.3
		#		c2				-0.1			-0.8			0.7
		#
		#		c1				-0.4			-0.1			-0.3
		#		c2				0.1				0.8				-0.7
	#########

   numb_datast = ncol(interestedCoefficientData)
   res_pos = apply(interestedCoefficientData, 1, function(x) sum(x > correlThreshold))
   res_neg = apply(interestedCoefficientData, 1, function(x) sum(x < correlThreshold))
   res_corr = cbind(res_pos, res_neg)
   
   pass_consis = apply(res_corr, 1, function(x) {ifelse(any(x == numb_datast), 1, 0)} )   #  numb_datast %in% x
   #print(head(res_corr))
   #print(head(pass_consis))
   
   
   ldf = cbind(ldf,combinedCoefficient, pass_consis)
   # add the group name to the header
   colnames(ldf)[colnames(ldf)=="combinedCoefficient"] <- paste("combCategCoeff", search_gr, sep='_')
   colnames(ldf)[colnames(ldf)=="pass_consis"] <- paste("pass_consis", search_gr, sep='_')
   #print(head(ldf))
   

   #calculate FDR for combined pvalue
   combinedCategFDR = p.adjust(ldf[,paste("combCategCorPVal", search_gr, sep='_')],method="fdr")
   ldf = cbind(ldf,combinedCategFDR)
   # add the group name to the header
   colnames(ldf)[colnames(ldf)=="combinedCategFDR"] <- paste("combCategFDR", search_gr, sep='_')
   #print(head(ldf))
   
   }
   return(ldf)
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
    
    pdf(paste(search_group_santized, str, 'FDRvsPUC.pdf', sep='-'))
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
   #print(head(outForPUC))
   
   # calculate combined Pvalue for interest group
   PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
   PvalueColnames = PvalueColnames[grep(paste(search_group,collapse=" vs "),PvalueColnames)]
   #print(PvalueColnames)
   
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
   #print(head(outForPUC))
   
   
   # calculate median (difference of coefficient)
   CoefficientColnames = colnames(outForPUC)[grep("DiffCorrs",colnames(outForPUC))]
   CoefficientColnames = CoefficientColnames[grep(paste(search_group,collapse=" vs "),CoefficientColnames)]
   #print(CoefficientColnames)
   
   interestedCoefficientData = outForPUC[,CoefficientColnames]
   interestedCoefficientData = as.matrix(interestedCoefficientData)
   interestedCoefficientData = apply(interestedCoefficientData,2,function(x){as.numeric(as.vector(x))})
   combinedDiffCoefficient = apply(interestedCoefficientData,1, function(x){round(median(x, na.rm = TRUE), 3)})
   outForPUC = cbind(outForPUC,combinedDiffCoefficient)
   
   result = outForPUC 
   #calculate FDR for combined pvalue
   combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
   result = cbind(result,combinedFDR)
   #print(head(result))
   #write.csv (result,paste( outputFile, paste(search_group_santized,collapse='_'), "tmp-out.csv", sep='-'))
   
   
 
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


pow <- function(x=10, y=6) {
   # function to print x raised to the power y
   result <- x^y
   return(result)
}


unlog <- function(FoldChangeData, base){
	if(logbase != 0){
	    if(logbase == 1) {
	          FoldChangeData = apply(FoldChangeData, 2, function(x) {exp(x)})    # using the base e
	    } else {
	          FoldChangeData = apply(FoldChangeData, 2, function(x) {pow(base, x)})
	    }
	}
	return(FoldChangeData)
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
	
	#print(colnames(outForPUC))
	outForPUC = calc_combined_coeff_pval_fdr(outForPUC)
	#print(head(outForPUC))
	
	IfFoldChangeDirectionMatch = ''	
    if(noPUC){
           # no need to look up the fold change
    } else {
        # attach the foldChange information for each partner
        FoldChangeCol = grep("combined" , colnames(FoldChangeMetabolic), value=TRUE, fixed=TRUE)
		FoldMetab1_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner1"]), c("geneName", FoldChangeCol)]
		colnames(FoldMetab1_InPair) = c("partner1InFold","partner1_FoldChange")
		FoldMetab2_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner2"]),c("geneName", FoldChangeCol)]
		colnames(FoldMetab2_InPair) = c("partner2InFold","partner2_FoldChange")

        FoldMetab1_InPair = cbind(FoldMetab1_InPair[, 1, drop=F], unlog(FoldMetab1_InPair[, 2, drop=F], logbase))
		FoldMetab2_InPair = cbind(FoldMetab2_InPair[, 1, drop=F], unlog(FoldMetab2_InPair[, 2, drop=F], logbase))
		#print(head(FoldMetab1_InPair))
		#print(head(FoldMetab2_InPair))
		
		outForPUC = cbind(outForPUC,FoldMetab1_InPair,FoldMetab2_InPair)
		
		# calculate fold change direction for each partner
		FoldChangeColnames = colnames(outForPUC)[grep("FoldChange",colnames(outForPUC))] # since this is using the combined fold change calculated above, you do not need the foldchVar variable
		FoldChangeData = outForPUC[,FoldChangeColnames]
		

		#print(head(FoldChangeData))
		
		FoldChangeDirection = NA

		#if (foldchthresh==0){
		#    FoldChangeDirection <- as.matrix(FoldChangeData)
		#    FoldChangeDirection[FoldChangeDirection<0] <- -1
		#    FoldChangeDirection[FoldChangeDirection>=0] <- 1
		#    FoldChangeDirection = as.data.frame(FoldChangeDirection)
		#    #print(head(FoldChangeDirection))
		#} else { # for unlog data
			FoldChangeDirection = (FoldChangeData-1)/abs(FoldChangeData-1)
		#}
		names(FoldChangeDirection) = c()
		colnames(FoldChangeDirection) = paste(colnames(FoldChangeData),"FoldChangeDirection",sep=".")
		
		# calculate if fold change direction are the same for the two partners
		IfFoldChangeDirectionMatch = apply(FoldChangeDirection,1,prod)
		names(IfFoldChangeDirectionMatch) = c()
		
		outForPUC = cbind(outForPUC,FoldChangeDirection,IfFoldChangeDirectionMatch)
	}
	
	#print(head(outForPUC))
	

	#########
	# calc puc using the normal correlation, instead of diff correlation
	#
    # calculate correlation Direction For combined correlation coefficient for each group of interest (combCategCoeff) 
	###### IF inconsistent just treat the puc as 0
		## check if the pass consis column is 1 and continue with the rest else set the puc column to 0
	
	
	#print("there")
	print(colnames(outForPUC))
	
	#write.csv(outForPUC, "tmp-nopuc-file.csv", quote=F)
	
    for(search_gr in search_group){
	    #print(search_gr)
		pattrn = paste("combCategCoeff", search_gr, sep='_')
		
		# why not directly build the required string instead of using grep to get the same string. did this to get around the case where one gr is substring of the other, e.g. --group stool4wk AbxHFHS\+nor stool4wk AbxHFHS
		# while the \\b before and after, or ^ and $ before and after worked for normal strings, it seems it does not work for when we have special chars (\\, + , etc)
		#interestedCoefficientColnames = grep(pattrn,colnames(outForPUC), value=TRUE, fixed=TRUE)
		interestedCoefficientColnames = pattrn
		print(interestedCoefficientColnames)
		
		interestedCorrelationData = outForPUC[,interestedCoefficientColnames, drop=FALSE]
		interestedCorrelationData = apply(interestedCorrelationData,2,function(x){as.numeric(as.vector(x))})
		signOfInterestedCorrelationData = interestedCorrelationData/abs(interestedCorrelationData)
		rownames(signOfInterestedCorrelationData) = c()
		colnames(signOfInterestedCorrelationData) = paste(colnames(interestedCorrelationData),"correlationDirection",sep=".")
		matchedExpressionDirection = signOfInterestedCorrelationData
		colnames(matchedExpressionDirection) = c()

		
		outForPUC = cbind(outForPUC,signOfInterestedCorrelationData)
		
		if(!noPUC){
			
			local_consis_search = paste("pass_consis", search_gr, sep='_')
			#print(head(outForPUC[,local_consis_search,drop=F]))
			
			# use "matchedExpressionDirection" and "IfFoldChangeDirectionMatch" to calc PUC, i.e. if these two are the same PUC=1 (good)
			pdff = outForPUC[, local_consis_search,drop=F]
			
			prod_fc_corr_dir = IfFoldChangeDirectionMatch * matchedExpressionDirection
			
			pdff = cbind(pdff, IfFoldChangeDirectionMatch, matchedExpressionDirection, prod_fc_corr_dir)
			
			#print(head(pdff, 6))
			# if the corr is not consistent in cohorts, make puc 0
			PUC = ifelse(pdff[, local_consis_search] == 1, pdff$prod_fc_corr_dir,0 )
			
			outForPUC = cbind(outForPUC,PUC)
			colnames(outForPUC)[colnames(outForPUC)=="PUC"] <- paste("PUC", search_gr, sep='_')
			#print(head(outForPUC))
        }
		
    }   
	return(outForPUC)
		
} 

    PUCoutfile = ''
    if(noPUC){
           PUCoutfile = paste(outputFile, paste(search_group_santized,collapse=" vs "), "noPUC-output.csv",sep='-')
    } else {
         PUCoutfile = paste(outputFile, paste(search_group_santized,collapse=" vs "), "PUC-output.csv",sep='-')
    }  

    if(file.exists(PUCoutfile))
    {
        print("Deleted old PUC output")
        file.remove(PUCoutfile)
    }

	
	#print(head(result))
	
    result = forPUC(FoldChangeMetabolic,noPUC)
	write.csv (result,networkFile, quote=FALSE)

    #q()
	data = result
	
	if(FALSE){
	
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

		pos_plt = calc_PUC_at_thresholds(sorted_result[which(sorted_result$combinedDiffCoefficient.correlationDirection == 1),], "pos") # walk along fdr and calc puc at diff fdr
		write.csv(pos_plt, paste(PUCoutfile, "pos-edges.csv", sep='.'), row.names=FALSE)

		neg_plt = calc_PUC_at_thresholds(sorted_result[which(sorted_result$combinedDiffCoefficient.correlationDirection == -1),], "neg") # walk along fdr and calc puc at diff fdr
		write.csv(neg_plt, paste(PUCoutfile, "neg-edges.csv", sep='.'), row.names=FALSE)

		data = sorted_result
		}
	}



###########################################################################################################
#generate network
###########################################################################################################


#----------------------------------------------------------------------
# calc stats of the generated network
#----------------------------------------------------------------------
calc_stats = function(inNet, correlThreshold=0){

        out = file(paste(networkFilewcuts,'-stats-log.txt',sep='') ,'a')
        
		if(FALSE){
			combcoeffcol = as.vector(inNet[,"combinedDiffCoefficient"])
			# not sure why it behaves this way!
			lowest_neg_corrl = min(combcoeffcol[combcoeffcol < 0])
			highest_neg_corrl = max(combcoeffcol[combcoeffcol < 0])

			lowest_pos_corrl = min(combcoeffcol[combcoeffcol > 0])
			highest_pos_corrl = max(combcoeffcol[combcoeffcol > 0])

			
			coefficientData = inNet[,grep("correlationDirection",colnames(inNet)),drop=FALSE]
			
			## artifically add row (with 0) incase it is only 1 row
			if(nrow(coefficientData) == 1){
				tmp_v = c("0")
				tmp_artif = data.frame(tmp_v)
				colnames(tmp_artif) = c("combinedDiffCoefficient.correlationDirection")
				coefficientData = rbind(coefficientData, tmp_artif)
			}
			
			coefficientData = apply(coefficientData,2,function(x){as.numeric(as.vector(x))})
			
			res_pos = apply(coefficientData, 2, function(x) sum(x == 1))
			res_neg = apply(coefficientData, 2, function(x) sum(x == -1))
			res = cbind(res_pos, res_neg)
			#print(res)
        }
        #print(paste(c("Total number of edges: ", nrow(inNet)), collapse="")) 

        setpartner1 = unique(as.vector(inNet[,"partner1"]))
        setpartner2 = unique(as.vector(inNet[,"partner2"]))

        nodes = union(setpartner1, setpartner2)
        #print(paste(c("Number of unique nodes: ", length(nodes)), collapse=""))

        edgesDistinctNodes = inNet[which(inNet[,"partner1"]!=inNet[,"partner2"]), ] 
        #print(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse=""))

        #write.csv(paste(c("Lowest pos corr:" , lowest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        #write.csv(paste(c("Highest pos corr:" , highest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        #write.csv(paste(c("Lowest neg corr:" , lowest_neg_corrl), collapse='') , file=out, row.names=FALSE)
        #write.csv(paste(c("Highest neg corr:" , highest_neg_corrl), collapse='') , file=out, row.names=FALSE)

        #write.csv(paste(c("Pos corr edges:" , toString(res_pos)), collapse='') , file=out, row.names=FALSE)
        #write.csv(paste(c("Neg corr edges:" , toString(res_neg)), collapse='') , file=out, row.names=FALSE)
        #write.csv(paste(c("Ratio of Pos to Neg corr edges:" , toString(res_pos/res_neg)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Total number of edges: ", nrow(inNet)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of unique nodes: ", length(nodes)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Ratio of edges to nodes: ", nrow(inNet)/length(nodes)), collapse='') , file=out, row.names=FALSE)
        #write.csv(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse='') , file=out, row.names=FALSE)
        close(out)

}


#----------------------------------------------------------------------
# for a given combinedPvalue and combinedfdr threshold, generate network
#----------------------------------------------------------------------
generateNetwork = function(){
       #generate network After Filter

    out = data
	if(FALSE){
		if(noPUC){
	    # at this point you only have pairs with consistent correlations 
	    out = data
	    #out = out[!is.na(out$"PUC"),] # remove the rows with 'NA' in PUC columns  # decide whether to use this combinedDiffCoefficient.correlationDirection
        } else {
	    # find PUC expected
	    out = data[data[,"PUC"]==1 ,]
	    out = out[!is.na(out$"PUC"),] # remove the rows with 'NA' in PUC columns
        }
	}
	out = as.matrix(out)
    	
	out = out[(as.numeric(out[,"combinedPvalue"])<combinedPvalueCutoff)==1,]
	out = out[(as.numeric(out[,"combinedFDR"])<combinedFDRCutoff)==1,]
	# find significant in individual pvalue: all of the pvalues for all of the datasets for each pair must be smaller than threshold
	# find pvalue data
  # if it is a single dataset it becomes a vector so adding the drop=FALSE
	pvalueData = out[,grep("pvalue",colnames(out)), drop=FALSE]
	pvalueData = pvalueData[,grep(paste(search_group,collapse=" vs "),colnames(pvalueData)), drop=FALSE]
		
	pvalueData = as.matrix(pvalueData)
	pvalueData = apply(pvalueData,2,function(x){as.numeric(as.vector(x))})
	
	# added this on july 17 2018
	if(is.null(nrow(pvalueData))){
		print("Qutting since no edges passed the fisher pval or fdr cuts.")
		q()
	}
	# calculate the largest pvalue among all datasets for each gene, this smallest pvalue must be smaller than threshold
	passIndevidualPvalue = apply(pvalueData,1,max)<individualPvalueCutoff
	#print(passIndevidualPvalue)
	outNetwork = out[passIndevidualPvalue, , drop=FALSE]
	#print(outNetwork)
        
        #outNetwork$names<-rownames(outNetwork)
        #splt = str_split_fixed(outNetwork$names, "<==>", 2)
        #outNetwork = cbind(outNetwork, splt)
        
	write.csv (outNetwork,networkFilewcuts, quote=FALSE)

	# the code does not handle pos and neg corrrelated stats for diff corr result files correctly
    calc_stats(outNetwork)

    print("Done!")
}


   generateNetwork()

