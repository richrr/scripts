library(stringr)


# usage:
#cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/microbe-pheno-analyses/analysis/m-p-pairs/rel.quantif_quantile/prelim_result/merged/corr/mw_sp
#Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/find-shortlisted-otu-pheno-pairs.R ./p1/per_analysis/stool-comp-out-no-cutoffs-freq-taxa-contrastfc.csv ./p1/per_analysis/stool-corr-out.csv stool-no-cutoff-sel-otu-pheno-pairs

args = commandArgs(trailingOnly=TRUE)
#print(length(args))


#------------------------------------------------------------------
# this calculates fisher pvalue based on different analysis
#------------------------------------------------------------------
calcFisherPval = function(outForPUC, searchterm){
   # calculate combined Pvalue for interest group
   PvalueColnames = colnames(outForPUC)[grep(searchterm,colnames(outForPUC))]
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
	return(combinedPvalue)
}


#------------------------------------------------------------------
# this calcs FDR
#------------------------------------------------------------------
calcFDR = function(result){
   #calculate FDR for combined pvalue
   combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
   return(combinedFDR)
}



#------------------------------------------------------------------
# this finds the analysis number along with the analys string
#------------------------------------------------------------------
find_analysis_number_w_analys_and_group = function(x){
  #res = as.numeric(str_match(x, "Analys ([0-9]+) .*")[,2])
  res = str_match(x, "(Analys [0-9]+-?[0-9]* .+) Spearman.*")[,2]
  res
}


#------------------------------------------------------------------
# this finds the analysis number along with the analys string
#------------------------------------------------------------------
find_analysis_number_w_analys = function(x){
  #res = as.numeric(str_match(x, "Analys ([0-9]+) .*")[,2])
  res = str_match(x, "(Analys [0-9]+-?[0-9]*) .*")[,2]
  res
}


#------------------------------------------------------------------
# this finds the analysis name along with the expt number
#------------------------------------------------------------------
find_analysis_name_w_expt_numb = function(x){
  res = str_match(x, "Analys [0-9]+-?[0-9]* .* Median (.*)")[,2]
  res
}




add_analysis_group_info = function(colheaders){
	#print(colheaders)
	previous_info = NA
	newcols = c()
	for(c in colheaders){
		newc = c
		#print(newc)
		
		info = find_analysis_number_w_analys_and_group(c)
		#print(info)
		
		# reset previous_info if info is not NA
		if(!is.na(info)){
			previous_info = info
			#print(paste("previous_info is now" , info, sep=" "))
		} 
		
		# if previous_info is not NA but info is then join the previous info to the colname
		if(!is.na(previous_info) & is.na(info)) {
			newc = paste(previous_info, c, sep=' ')
			#print(newc)
		}
		newcols = c(newcols, newc)
	}
	
	#print(newcols)
	return(newcols)
}


round_df <- function(df, digits=3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}

add_taxa_freq_cols = function(ids, cmpdf){
	rdf = NA
	flag = TRUE
	
	counter = 0
	for(r in ids){
		otu_ph = unlist(str_split(r , "<==>"))
		res = cmpdf[otu_ph[1], grep("FOLDCHANGE", colnames(cmpdf), value = T, invert=T)]
		#print(res)
		
		newcs = c()
		for(c in colnames(res)){
			
			if(c == 'ID' || grepl("df\\$taxa" , c) || grepl("contrastFC", c) || grepl("freqplayer", c)) {
				newc = c
			} else {
				numb = find_analysis_number_w_analys(c)
				name = find_analysis_name_w_expt_numb(c)
				newc = paste(numb, name, sep=" ")
			}
			newcs = c(newcs, newc)
		}
		
		colnames(res) = newcs
		
		if(flag){
			rdf = res
			flag = FALSE
		} else{
			rdf = rbind(rdf, res)
		}
		
		#if(counter == 2){
		#quit()
		#}
		#counter = counter + 1
	}
	
	rdf = round_df(rdf, 5) # since you can convert to percentage in excel
	return(rdf)
}

compfile = args[1]
corrfile = args[2]
outstr = args[3]

# merged consistent comp file which has the info about contrastfc and freq players
cmpdf = read.csv( compfile, header=TRUE,check.names=FALSE, stringsAsFactors=FALSE)
rownames(cmpdf) = cmpdf[,1]
#colnames(cmpdf)


# merged consistent corr file which has the info about the otu pheno pairs
corrdf = read.csv( corrfile, header=TRUE,check.names=FALSE, stringsAsFactors=FALSE)
rownames(corrdf) = corrdf[,1]


# only drop these since we still need to assign the analysis and group info
# remcols1= grep( "partner1[^_]" , colnames( corrdf ), fixed=TRUE, value=TRUE ) # this keeps "_"
# remcols2= grep( "partner2[^_]" , colnames( corrdf ), fixed=TRUE, value=TRUE ) # this keeps "_"
# remcols3= grep( "partner1InFold" , colnames( corrdf ), fixed=TRUE, value=TRUE )
# remcols4= grep( "partner2InFold" , colnames( corrdf ), fixed=TRUE, value=TRUE )
# remcols5= grep( "Direction" , colnames( corrdf ), fixed=TRUE, value=TRUE )
# remcols6= grep( "PUC" , colnames( corrdf ), fixed=TRUE, value=TRUE )
# remcols7= grep( "FDR" , colnames( corrdf ), fixed=TRUE, value=TRUE )
# remcols = c(remcols1, remcols2, remcols3, remcols4, remcols5, remcols6, remcols7)
# corrdf = corrdf[ , -which(colnames(corrdf) %in% remcols)]  


# only drop these since we still need to assign the analysis and group info
remcols1= grep( "partner[1,2][^_].*|partner[1,2]$" , colnames( corrdf ), value=TRUE ) # this keeps "_"
remcols2= grep( "Direction" , colnames( corrdf ), fixed=TRUE, value=TRUE )
remcols3= grep( "PUC" , colnames( corrdf ), fixed=TRUE, value=TRUE )
remcols4= grep( "FDR" , colnames( corrdf ), fixed=TRUE, value=TRUE )
remcols = c(remcols1, remcols2, remcols3, remcols4)
corrdf = corrdf[ , -which(colnames(corrdf) %in% remcols)]  


# assign the analysis and group info
colnames(corrdf) = add_analysis_group_info(colnames(corrdf))

# now drop the individual pval cols since we don't need them and also have the analysis and group info
remcols8= grep( "pvalue_" , colnames( corrdf ), fixed=TRUE, value=TRUE )
corrdf = corrdf[ , -which(colnames(corrdf) %in% remcols8)]  

contrstfcCol = "contrastFC"
if("contrastFC.x" %in% colnames(cmpdf)){
	contrstfcCol = "contrastFC.x"
}

selOTUs = ''
if(length(args) == 3){
	# select the otus that have 1 in contrastFC and freqplayer
	selOTUs = cmpdf[cmpdf[,contrstfcCol]==1 & cmpdf[,"freqplayer"]==1,"ID"]
} else if (length(args) > 3 && args[4] == "nocontrastfc"){
# select the otus that have 1 in freqplayer
	selOTUs = cmpdf[cmpdf[,"freqplayer"]==1,"ID"]
}
# add the "<==>" string to avoid selecting other otus that are partial match of the shortlisted otus
otu_pairs = paste(selOTUs, "<==>", sep='')


resdf = subset(corrdf , grepl(paste(otu_pairs, collapse= "|"), ID))


combinedPvalue = calcFisherPval(resdf, "combinedPvalue")
resdf_w_Fisherp = cbind(resdf, combinedPvalue)
combinedFDR = calcFDR(resdf_w_Fisherp)
resdf = cbind(resdf_w_Fisherp, combinedFDR)

resdf = round_df(resdf, 2)

cols_to_add = add_taxa_freq_cols(resdf$ID, cmpdf)

#### full table ####
resdf_full = cbind(resdf, cols_to_add)
write.csv(resdf_full, paste(outstr , ".csv", sep=""), row.names=FALSE, quote=F)


# note that for correl, we might want to use something different.
DirectionFunc = function(col){
	
    newcol = tryCatch(
			{
				as.numeric(col)/abs(as.numeric(col))
			},
			error=function(cond) {
            # Choose a return value in case of error
            return(col)
			},
			warning=function(cond) {
            # Choose a return value in case of warning
            return(col)
			}
			)
	newcol
}

# check if not na and assign 1 to numeric
KeepNumericAs1 = function(x){
	out = x
	if(!is.na(x) & suppressWarnings(all(!is.na(as.numeric(as.character(x))))))
	{
	   		out = 1
	}
	# otherwise do nothing as out stays the string or na
	out

}


summarize_table = function(resdf){
	# keep only the first column and combined coeff col
	keepCols = grep("combinedCoefficient" , colnames(resdf), fixed=T, value=T)
	keepCols = c("ID", keepCols)
	
	resdf = resdf[, keepCols]
	resdf = apply(resdf,c(1,2), KeepNumericAs1) 
	return(resdf)
}

### summary table ###
keepasis = resdf[, c("combinedPvalue", "combinedFDR")]
resdf = summarize_table(resdf)
resdf = cbind(resdf, keepasis)
resdf_summary = cbind(resdf, cols_to_add)
write.csv(resdf_summary, paste(outstr , "-summary.csv", sep=""), row.names=FALSE, quote=F)

#warnings()