
args = commandArgs(trailingOnly=TRUE)


# usage: # /nfs3/PHARM/Morgun_Lab/richrr/Lean_Abx_data/analysis/merged/corr/microb-pheno_corr_all_samps/p1/per_analysis/all_samples/ip_0.2_cp_0.05_otu-fdr_0.01_non-otu-fdr_0.1/non OTU edges network/selected_metab_tk_net
# Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/calc-metapval.R corr-w-gr-info.csv combinedPvalue _combinedCoefficient one

infile = args[1]
search_pval = args[2]
search_coeff = args[3]
sided=args[4]

outForPUC = read.csv(infile, header=T, check.names=F)
#head(outForPUC)
dim(outForPUC)

# keep only rows that have consistent correlations
  CoeffColnames = colnames(outForPUC)[grep(search_coeff,colnames(outForPUC))]
  print(CoeffColnames)

  coeffDf = outForPUC[,CoeffColnames]
  head(coeffDf)

  coeffDf = outForPUC[,CoeffColnames]
  head(coeffDf)

  out = apply(coeffDf, 1, function(coeffs){
										notNAcoeffs = coeffs[!is.na(coeffs)]
										all_pos = sum(coeffs > 0, na.rm=TRUE)
										all_neg = sum(coeffs < 0, na.rm=TRUE)
										
										res = 0
										if(all_pos==length(notNAcoeffs) || all_neg==length(notNAcoeffs))
										{
										    res=TRUE
										} else{
										    res=FALSE
										}
									})
  
  out = unlist(out)
  #out
  
  # these are consistent correlations
  outForPUC = outForPUC[out, ]
  

  
# calculate combined Pvalue for interest group
   PvalueColnames = colnames(outForPUC)[grep(search_pval,colnames(outForPUC))]
   print(PvalueColnames)
   
   
   interestedPvalueData = outForPUC[,PvalueColnames]
   #head(interestedPvalueData)
   # divide the pvalues by 2 to make it one-sided
   if(sided=="one"){
       interestedPvalueData = interestedPvalueData/2
   }
   head(interestedPvalueData)
   

   
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
   head(outForPUC)
   
write.csv(outForPUC, paste0(infile, "consisCoeff-metap.csv"), row.names=FALSE, quote=FALSE)
   