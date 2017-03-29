#==================================================================================================================
#		calculate all gene or pairs, comparison or correlation,pvalue,fdr
#==================================================================================================================

# conditions: 
# add 4th column with "correlation" in analysis file to indicate whatever operation is to be done, needs to be done on correlations
# 1. A vs B: 2 non-empty columns, last two columns empty [A	B empty empty]
# 		# for an anova like analysis give the comparison as 
#		 	# ";" separated comparison as column1 e.g. "A_vs_B;C_vs_D" in decreasing order of priority
#				since the venn diagram only plots the top three comparisons
#  				## note that in each comparison "_vs_" separates the two categories being compared
# 				## code replaces "_vs_" by "-" and then use to make contrasts
#			# ";" separated categories as column2 e.g. "A;B;D;C"
# 				## these are used to create the levels and column names


# 2. delta A vs delta B: 3 non-empty columns, 3rd column delta, 4th column empty [A	B	delta	empty]  # needs 4 groups: two comma separated from A and two comma separated from B

# 3. A vs B paired:  3 non-empty columns, 3rd column paired, 4th column empty [A	B	paired	empty]

# 4. Correl A: 2 non-empty columns, 4th column correlation [A	empty	empty	correlation]

# 5. Correl delta A: 3 non-empty columns, 2nd column empty, 3rd column delta, 4th column correlation [A	empty	delta	correlation] # needs 2 groups: two comma separated from A 

# 6. Correl A vs. Correl B: 3 non-empty columns, 4th column correlation [A	B	empty	correlation]

# 7. Correl delta A vs. Correl delta B: 4 non-empty columns, 3rd column delta, 4th column correlation [A	B	delta	correlation] # needs 4 groups: two comma separated from A and two comma separated from B

# 8. Correl A vs. Correl B paired: 4 non-empty columns, 3rd column paired, 4th column correlation [A	B	paired	correlation]
## no need to do paired in 8, so skipping this, because without paired, this defaults to condition 6

# 9. Timelag Correl A and B (where A and B are different time points): 4 non-empty columns, 3rd column timelag, 4th column correlation [A	B	timelag	correlation]
## std. correlations: g1T1-g2T1 or g1T2-g2T2
## calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc. 
## for all gene pairs: give group 1 (T1) to first gene, group 2 (T2) to second gene


# the indexes allow creating a unique column which does not interfere with merge. merge likes to have unique column names.


#http://stackoverflow.com/questions/9455167/converting-data-frame-factors-with-data-matrix-weird-results
#else as.numeric gives weird numbers instead of actual values
numericizeVector = function(c){
    nc = sapply(c, function(x) as.numeric(as.character(x)))
    return(nc)
}


#==================================================================================================================
#		Condition 4
#==================================================================================================================
#     function: calculate correlation coefficient of a pair of genes
# input:
#         1.pair: a pair of "gene"s:
#                 value: a vector of length 2, containing names of the two "genes"
#         2. expression data
#	  3. the category on which to perform correlations between measurements
#	  4. the dictionary containing mapping of category and samples which belong to it
#         5. method for correlation
# output:
#	 vector of length 3: method, "correlation Coefficient","pvalue"

CalcCor = function(pair, edata, c, dict, correlMethod){
    idxs = as.vector(dict[[c]])

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste(rev(pair),collapse=" ___ ")
    cat (calculating)
    cat ("\n")

    c1 = numericizeVector(as.vector(edata[pair[1], idxs]))
    c2 = numericizeVector(as.vector(edata[pair[2], idxs]))

    outLine = ''
    # there are less than 3 samples with values (i.e. without NA)
    if (sum(!is.na(c1+c2))<3){
        outLine = as.matrix(c(NA,NA,NA))
    }else{
    p = cor.test(c1,c2,method=correlMethod)
    outLine = as.matrix(c(p$method,p$estimate,p$p.value))
    }
    # no point in adding row names since it is not retained when it returns from the apply method
    #rownames(outLine) = c("correlation.Coefficient", paste(p$method,"pvalue",sep=' '))
    outLine
}


#==================================================================================================================
#		Get number of samples, used for Condition 6
#==================================================================================================================
CalcNCor = function(pair, edata, c, dict, correlMethod){
    idxs = as.vector(dict[[c]])

    c1 = numericizeVector(as.vector(edata[pair[1], idxs]))
    c2 = numericizeVector(as.vector(edata[pair[2], idxs]))

    # number of non-na samples
    N = sum(!is.na(c1+c2))
    return(N)
}

#==================================================================================================================
#		Get number of samples, used for Condition 7
#==================================================================================================================
CalcNDeltaCor = function(pair, edata, c1, c2, dict, correlMethod){
    idxs1 = as.vector(dict[[c1]])
    idxs2 = as.vector(dict[[c2]])

    c1 = numericizeVector(as.vector(edata[pair[1], idxs1]))
    c2 = numericizeVector(as.vector(edata[pair[1], idxs2]))

    c3 = numericizeVector(as.vector(edata[pair[2], idxs1]))
    c4 = numericizeVector(as.vector(edata[pair[2], idxs2]))

    # calculate delta 
    #ca = c1-c2
    #cb = c3-c4

    # calculate delta normalized to starting value
    ca = (c1-c2)/c2
    cb = (c3-c4)/c4

    # number of non-na samples
    N = sum(!is.na(ca+cb))
    return(N)
}


#==================================================================================================================
#		Condition 4
#==================================================================================================================
# function : calculate correlation,pvalue,fdr for all gene pairs
# input:
#    	  1. pairs: for which to calculate the correlation
#         2. expression data
#	  3. the category on which to perform correlations between measurements
#	  4. the dictionary containing mapping of category and samples which belong to it
#         5. method for correlation
# output : 
#	a table containing pair information, correlation,pvalue,fdr

calculateCorrelation = function (pairs, expressionData, c, dict, correlMethod, indxg){

    # calculate for each pair
    out = apply(pairs, 1, CalcCor, expressionData, c, dict, correlMethod)

    # add the pair as the column name
    colnames(out) = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")

    # the pairs become row names; and the corr coeff and pvalue are the column names
    out = t(out)

    # method used; # just take the first non-NA element
    tmp_uniq = unique(out[,1]) 
    NonNAindex <- which(!is.na(tmp_uniq))[1]
    method_label = tmp_uniq[NonNAindex]

    # discards the first column of method used
    out = out[,-1]    
    # append method used to column labels
    colnames(out) = c(paste("Analys", indxg, c, method_label,"Coefficient",sep=" "), paste("Analys", indxg, c,method_label,"pvalue",sep=" "))

    # calculate FDR
    FDR = p.adjust(out[,colnames(out)[grep("pvalue",colnames(out))]],method="fdr")
    oldColnames = colnames(out)
    out = cbind(out, FDR)
    colnames(out) = c(oldColnames,paste("Analys", indxg, c,method_label,"FDR",sep=" "))
    pairName = rownames(out)
    out = cbind(pairName,out)
    out
}


#==================================================================================================================
#		Partial correlation without pvalue and FDR
#==================================================================================================================
# function : calculate Partial correlation for all gene pairs. this does NOT give you pvalue and FDR
# input:
#    	  1. genes: for whose pairs to calculate the partial correlation
#         2. expression data
#	  3. the category on which to perform correlations between measurements
#	  4. the dictionary containing mapping of category and samples which belong to it
#         5. method for correlation
# output : 
#	a table containing pair information, correlation

calculatePartialCorrelation = function (pairs, genes, edata, c, dict, correlMethod, indxg){
    # extract the required rows and columns using the gene list and the categ
    idxs = as.vector(dict[[c]])
    subset_edata = as.matrix(edata[genes, idxs])
    subset_edata[subset_edata == ""] <- NA
    subset_edata[subset_edata == "na"] <- NA
    class(subset_edata) <- "numeric"
    #print(subset_edata)

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste("Calculating Partial Correl", " all pairs of specified gene list", sep=" ___ ")
    cat (calculating)
    cat ("\n")
 
    Corr_ = cor(t(subset_edata), method=correlMethod, use="pairwise.complete.obs") # this calculates correlation between genes  "pairwise.complete.obs" "na.or.complete"
    # this will be used to set the row names
    tmp_MC_Corr_ = melt(Corr_)

    C = make.positive.definite(Corr_)
    C = cor2pcor(C)
    # The output from below might be 3 column (X1, X2, value), with rows labeled lost
    tmp_MC = melt(C)

    out = matrix(tmp_MC[,"value"])
    # the pairs become row names
    row.names(out) = paste(tmp_MC_Corr_[,"X1"],tmp_MC_Corr_[,"X2"],sep="<==>")

    # method used
    method_label = correlMethod

    # append method used to column labels
    colnames(out) = c(paste("Analys", indxg, c, method_label,"Partial Correl Coefficient",sep=" "))

    # select only the pairs that we need: i.e. if you have gene1-gene2, you don't need gene2-gene1
    select_pairs = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")
    tttt = t(subset( t(as.matrix(out)) , select = c(select_pairs)))
    pairName = rownames(tttt)
    out = cbind(pairName,tttt)
    out
}


#==================================================================================================================
#  DO NOT USE. DOESN'T WORK: Partial correlation with pvalue and FDR # http://www.yilab.gatech.edu/pcor.html
#==================================================================================================================
# function : calculate Partial correlation for all gene pairs. this gives you pvalue and FDR
# input:
#    	  1. genes: for whose pairs to calculate the partial correlation
#         2. expression data
#	  3. the category on which to perform correlations between measurements
#	  4. the dictionary containing mapping of category and samples which belong to it
#         5. method for correlation
# output : 
#	a table containing pair information, correlation, pvalue and FDR

calculatePartialCorrelationWFDR = function (genes, edata, c, dict, correlMethod, indxg){
    # load the R code "pcor.test"
    source("pcor.R")

    # extract the required rows and columns using the gene list and the categ
    idxs = as.vector(dict[[c]])
    subset_edata = edata[genes, idxs]

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste("Calculating Partial Correl", " all pairs of specified gene list", sep=" ___ ")
    cat (calculating)
    cat ("\n")

    # you might have to create three items: vector gene1, vector gene2 and matrix( other genes)
    c1 = numericizeVector(as.vector(subset_edata[pair[1], ]))
    c2 = numericizeVector(as.vector(subset_edata[pair[2], ]))
    
    subset_edata_mat = as.matrix(t(subset_edata))
    pair_subset_edata = subset( subset_edata_mat , select = -c(get(pair[1]),get(pair[2])))

    pcor.test("Analys", indxg, c1,c2,pair_subset_edata, method=correlMethod) # this calculates correlation between genes

    
    #Corr_ = cor(t(subset_edata), method=correlMethod) 
    # this will be used to set the row names
    #tmp_MC_Corr_ = melt(Corr_)


    #C = make.positive.definite(Corr_)
    #C = cor2pcor(C)
    # The output from below might be 3 column (X1, X2, value), with rows labeled lost
    #tmp_MC = melt(C)
    #out = matrix(tmp_MC[,"value"])
    # the pairs become row names
    #row.names(out) = paste(tmp_MC_Corr_[,"X1"],tmp_MC_Corr_[,"X2"],sep="<==>")

    # method used
    #method_label = correlMethod

    # append method used to column labels
    #colnames(out) = c(paste(c, method_label,"Partial Correl Coefficient",sep=" "))


    #####
    # if you need the pvalue to calc FDR for the significance of the partial correlation use pcor.test
    # http://www.yilab.gatech.edu/pcor.html
    
    ####
    # calculate FDR
    #FDR = p.adjust(out[,colnames(out)[grep("pvalue",colnames(out))]],method="fdr")
    #oldColnames = colnames(out)
    #out = cbind(out, FDR)
    #colnames(out) = c(oldColnames,paste(c,method_label,"FDR",sep=" "))
    #pairName = rownames(out)
    #out = cbind(pairName,out)
    #out
}


#==================================================================================================================
#		Condition 9
#==================================================================================================================
#     function: calculate (timelag) correlation coefficient of a pair of genes (across two categories)
#               calculate timelag correlations for g?T1-g?T2 : will take care of g1T1-g1T2, g1T1-g2T2, etc. 
# input:
#         1.pair: a pair of "gene"s:
#                 value: a vector of length 2, containing names of the two "genes"
#         2. expression data
#	  3 & 4. the categories across which to perform timelag correlations between measurements
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for correlation
# output:
#	 vector of length 3: method, "correlation Coefficient","pvalue"

CalcTwoCategCor = function(pair, edata, c1, c2, dict, correlMethod){
    idxs1 = as.vector(dict[[c1]])
    idxs2 = as.vector(dict[[c2]])

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste(rev(pair),collapse=" ___ ")
    calculating = paste(calculating, c1, "timelag" , c2, sep=" ___ ")
    cat (calculating)
    cat ("\n")

    c1 = numericizeVector(as.vector(edata[pair[1], idxs1]))
    c2 = numericizeVector(as.vector(edata[pair[2], idxs2]))

    outLine = ''
    # there are less than 3 samples with values (i.e. without NA)
    if (sum(!is.na(c1+c2))<3){
        outLine = as.matrix(c(NA,NA,NA))
    }else{
    p = cor.test(c1,c2,method=correlMethod)
    outLine = as.matrix(c(p$method,p$estimate,p$p.value))
    }
    # no point in adding row names since it is not retained when it returns from the apply method
    #rownames(outLine) = c("correlation.Coefficient", paste(p$method,"pvalue",sep=' '))
    outLine
}


#==================================================================================================================
#		Condition 9
#==================================================================================================================
# function : calculate (timelag) correlation,pvalue,fdr for all pairs
# input:
#    	  1. pairs: for which to calculate the correlation
#         2. expression data
#	  3 & 4. the categories across which to perform correlations between measurements
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for correlation
# output : 
#	a table containing pair information, correlation,pvalue,fdr

calculateTwoCategCorrelation = function (pairs, expressionData, c1, c2, dict, correlMethod, indxg){

    # calculate for each pair
    out = apply(pairs, 1, CalcTwoCategCor, expressionData, c1, c2, dict, correlMethod)

    # add the pair as the column name
    colnames(out) = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")

    # the pairs become row names; and the corr coeff and pvalue are the column names
    out = t(out)
    
    # method used; # just take the first non-NA element
    tmp_uniq = unique(out[,1]) 
    NonNAindex <- which(!is.na(tmp_uniq))[1]
    method_label = tmp_uniq[NonNAindex]


    # discards the first column of method used
    out = out[,-1]    
    # append method used to column labels
    colnames(out) = c(paste("Analys", indxg, c1, "timelag" , c2, method_label,"Coefficient",sep=" "), paste("Analys", indxg, c1, "timelag" , c2,method_label,"pvalue",sep=" "))

    # calculate FDR
    FDR = p.adjust(out[,colnames(out)[grep("pvalue",colnames(out))]],method="fdr")
    oldColnames = colnames(out)
    out = cbind(out, FDR)
    colnames(out) = c(oldColnames,paste("Analys", indxg, c1, "timelag" , c2,method_label,"FDR",sep=" "))
    pairName = rownames(out)
    out = cbind(pairName,out)
    out
}


#==================================================================================================================
#		Condition 5
#==================================================================================================================
#     function: calculate (delta) correlation coefficient of a pair of genes using difference between two categories
# input:
#         1.pair: a pair of "gene"s:
#                 value: a vector of length 2, containing names of the two "genes"
#         2. expression data
#	  3 & 4. the categories whose difference is to be used to perform correlations between measurements
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for correlation
# output:
#	 vector of length 3: method, "correlation Coefficient","pvalue"

CalcDeltaCor = function(pair, edata, c1, c2, dict, correlMethod){
    idxs1 = as.vector(dict[[c1]])
    idxs2 = as.vector(dict[[c2]])

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste(rev(pair),collapse=" ___ ")
    calculating = paste(calculating, c1, "minus" , c2, sep=" ___ ")
    cat (calculating)
    cat ("\n")

    c1 = numericizeVector(as.vector(edata[pair[1], idxs1]))
    c2 = numericizeVector(as.vector(edata[pair[1], idxs2]))

    c3 = numericizeVector(as.vector(edata[pair[2], idxs1]))
    c4 = numericizeVector(as.vector(edata[pair[2], idxs2]))

    # calculate delta 
    #ca = c1-c2
    #cb = c3-c4

    # calculate delta normalized to starting value
    ca = (c1-c2)/c2
    cb = (c3-c4)/c4

    outLine = ''
    # there are less than 3 samples with values (i.e. without NA)
    if (sum(!is.na(ca+cb))<3){
        outLine = as.matrix(c(NA,NA,NA))
    }else{
    p = cor.test(ca,cb,method=correlMethod)
    outLine = as.matrix(c(p$method,p$estimate,p$p.value))
    }
    # no point in adding row names since it is not retained when it returns from the apply method
    #rownames(outLine) = c("correlation.Coefficient", paste(p$method,"pvalue",sep=' '))
    outLine
}


#==================================================================================================================
#		Condition 5
#==================================================================================================================
# function : calculate (delta) correlation,pvalue,fdr for all pairs using difference between two categories
# input:
#    	  1. pairs: for which to calculate the correlation
#         2. expression data
#	  3 & 4. the categories whose difference is to be used to perform correlations between measurements
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for correlation
# output : 
#	a table containing pair information, correlation,pvalue,fdr

calculateDeltaCorrelation = function (pairs, expressionData, c1, c2, dict, correlMethod, indxg){

    # calculate for each pair
    out = apply(pairs, 1, CalcDeltaCor, expressionData, c1, c2, dict, correlMethod)

    # add the pair as the column name
    colnames(out) = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")

    # the pairs become row names; and the corr coeff and pvalue are the column names
    out = t(out)

    # method used; # just take the first non-NA element
    tmp_uniq = unique(out[,1]) 
    NonNAindex <- which(!is.na(tmp_uniq))[1]
    method_label = tmp_uniq[NonNAindex]


    # discards the first column of method used
    out = out[,-1]    
    # append method used to column labels
    colnames(out) = c(paste("Analys", indxg, c1, "minus" , c2, method_label,"Coefficient",sep=" "), paste("Analys", indxg, c1, "minus" , c2, method_label,"pvalue",sep=" "))

    # calculate FDR
    FDR = p.adjust(out[,colnames(out)[grep("pvalue",colnames(out))]],method="fdr")
    oldColnames = colnames(out)
    out = cbind(out, FDR)
    colnames(out) = c(oldColnames,paste("Analys", indxg, c1, "minus" , c2, method_label,"FDR",sep=" "))
    pairName = rownames(out)
    out = cbind(pairName,out)
    out
}


library(Biobase)
library(limma)

# http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  gm = exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  return(gm)
}



###################
# returns the mean, meadian, and fold change for the elements in the two categories
###################
GetSimpleStats = function(lgene, edata, c1, c2, dict){

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste(lgene, c1, "vs" , c2, sep=" ___ ")
    cat (calculating)
    cat ("\n")

    idxs1 = as.vector(dict[[c1]])
    idxs2 = as.vector(dict[[c2]])

    #print(idxs1)
    #print(idxs2)

    c1 = as.vector(edata[lgene, idxs1])
    c2 = as.vector(edata[lgene, idxs2])

    c1 =  numericizeVector(c1) 
    #print(c1)
    c2 = numericizeVector(c2)  
    #print(c2)

    
    outLine = ''
    p = ''

 
    mean_c1 = NA
    mean_c2 = NA
    
    median_c1 = NA
    median_c2 = NA
    
    gm_c1 = NA
    gm_c2 = NA

    # there are at least 3 samples with values (i.e. without NA)
    if(sum(!is.na(c1)) >= 3){
        mean_c1 = mean(c1, na.rm=TRUE)
        median_c1 = median(c1, na.rm=TRUE)
        gm_c1 = gm_mean(c1, na.rm=TRUE)
    } 

    if(sum(!is.na(c2)) >= 3){
        mean_c2 = mean(c2, na.rm=TRUE)
        median_c2 = median(c2, na.rm=TRUE)
        gm_c2 = gm_mean(c2, na.rm=TRUE)
    } 

    
    outLine = as.matrix(c(mean_c1, mean_c2, median_c1 , median_c2 , gm_c1 , gm_c2))
    
    outLine
}

#==================================================================================================================
#		Conditions 1 and 3
#==================================================================================================================
# function : calculate comparison for a gene using two categories
# input:
#    	  1. genes: for which to calculate the comparison
#         2. expression data
#	  3 & 4. the categories to compare
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for comparison (t test, Man Whitney, etc.)
#         7. paired or unpaired
# output : 
#	vector of length 4: method, mean in c1, mean in c2,pvalue

CalcCom = function(lgene, edata, c1, c2, dict, comparMethod, pairedd){

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste(lgene, c1, "vs" , c2, sep=" ___ ")
    cat (calculating)
    cat ("\n")

    idxs1 = as.vector(dict[[c1]])
    idxs2 = as.vector(dict[[c2]])

    #print(idxs1)
    #print(idxs2)

    c1 = as.vector(edata[lgene, idxs1])
    c2 = as.vector(edata[lgene, idxs2])
        
    c1 =  numericizeVector(c1) 
    print(c1)
    c2 = numericizeVector(c2)  
    print(c2)

    
    outLine = ''
    p = ''

 
    mean_c1 = NA
    mean_c2 = NA
    
    median_c1 = NA
    median_c2 = NA

    gm_c1 = NA
    gm_c2 = NA

    # there are at least 3 samples with values (i.e. without NA)
    if(sum(!is.na(c1)) >= 3){
        mean_c1 = mean(c1, na.rm=TRUE)
        median_c1 = median(c1, na.rm=TRUE)
        gm_c1 = gm_mean(c1, na.rm=TRUE)
    } 

    if(sum(!is.na(c2)) >= 3){
        mean_c2 = mean(c2, na.rm=TRUE)
        median_c2 = median(c2, na.rm=TRUE)
        gm_c2 = gm_mean(c2, na.rm=TRUE)
    } 

    
    outLine = as.matrix(c(mean_c1, mean_c2, median_c1 , median_c2 , gm_c1 , gm_c2))

    outLine

}



#==================================================================================================================
#		Conditions 1 and 3
#==================================================================================================================


# function : calculate comparison,pvalue,fdr for all genes using two categories
# input:
#    	  1. list of genes: for which to calculate the comparison
#         2. expression data
#	  3 & 4. the categories to compare
#		# note that if there are multiple categories being comapred like an ANOVA
#		# give the comparison as 
#		 	# ";" separated comparison as column1 e.g. "A_vs_B;C_vs_D"
#  				## note that in each comparison "_vs_" separates the two categories being compared
# 				## code replaces "_vs_" by "-" and then use to make contrasts
#			# ";" separated categories as column2 e.g. "A;B;D;C"
# 				## these are used to create the levels and column names
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for comparison (t test, Man Whitney, etc.)
#         7. paired or unpaired
# output : 
#	a table containing gene information, means in c1 and c2,pvalue,fdr


calculateComparison = function (lgenes, expressionData, c1, c2, dict, comparMethod, pairedd, indxg){

    out = ''

    if(comparMethod == 'limma'){
            
            all_categs = c("MU","WT")   # use vector of categories to create colnames(design)
            all_idxs = c()
            all_idxs_lengths = c()
            Group = ''
            #comparisonsargs = ''
            
            ### check if c1 AND c2 have ";"
            # for anova like behavior, multi-categories
            if( (grepl(";", c1))&(grepl(";", c2)) ) {
                
                    # c2
                    all_categs = unlist(strsplit(c2, ";")) # split the c2 and create a vector of categories
                    for(catg in all_categs){              # loop over the categories to create all_idxs and all_idxs_lengths
                        tmp_v = as.vector(dict[[catg]])
                        all_idxs = c(all_idxs, tmp_v)
                        all_idxs_lengths = c(all_idxs_lengths, length(tmp_v))
                    }
                    #print(all_categs)
                    #print(all_idxs_lengths)
                    
                    targets=rep(all_categs, all_idxs_lengths) # use all_categs and all_idxs_lengths to create targets
                    Group <- factor(targets, levels=all_categs)
		    #print(Group)
		    
                    
            } else {
	    
		    # create the design
		    idxs1 = as.vector(dict[[c1]])
		    idxs2 = as.vector(dict[[c2]])
		    all_idxs = c(idxs1, idxs2)
		    
		    # see approach 2 on page 41 of https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
		    numbSamplc1 = length(idxs1)
		    numbSamplc2 = length(idxs2)
		    targets=rep(c(c1,c2), c(numbSamplc1, numbSamplc2))
		    Group <- factor(targets, levels=c(c1,c2))
		    print(Group)
		    
	    }
	    
	    
	    design <- model.matrix(~0+Group)
	    colnames(design) <- all_categs
	    rownames(design) <- all_idxs
	    #print(design)
	    
	    ### what about give all the list of genes and select only required later ###
	    # extract the required expression data
	    #subsetexpressionData = expressionData[lgenes, all_idxs]
	    subsetexpressionData = expressionData[, all_idxs]
	    #print(head(subsetexpressionData))
	    print(ncol(subsetexpressionData))
	    
	    
	    eset = ExpressionSet(data.matrix(subsetexpressionData))
	    #print(head(eset))
	    
	    fit <- lmFit(eset, design)

        if( (grepl(";", c1))&(grepl(";", c2)) ) {
            
			# c1
			res = strsplit(c1, ";")          # split the c1
			unlistedres = unlist(res)
			comparisonsargs = gsub("_vs_", "-" , unlistedres)     # create limma compatiable comparisons
			#print(comparisonsargs)
            
			# build the makecontrasts
			# https://support.bioconductor.org/p/27900/
			astr=paste(comparisonsargs, collapse=",")
			#Then add correct strings before and after.
			prestr="makeContrasts("
			poststr=",levels=design)"
			commandstr=paste(prestr,astr,poststr,sep="")
			#Now evaluate the command string
			contrast.matrix <- eval(parse(text=commandstr))
			#print(contrast.matrix)
			
			fit2 <- contrasts.fit(fit, contrast.matrix)
			fit2 <- eBayes(fit2)
			#print(fit2)
			
			#https://www.bioconductor.org/help/course-materials/2009/BioC2009/labs/limma/limma.pdf
			# https://stat.ethz.ch/pipermail/bioconductor/2012-November/049385.html
			#What limma does is always equivalent to Anova.  Of course it is not 
			#classical Anova, but rather an empirical Bayes extension of Anova that is 
			#appropriate for microarray data.

			#Post hoc tests are done in limma using decideTests(), and many options are 
			#offered.  You won't find classical methods like TukeyHSD though, because 
			#limma isn't doing classical Anova and because methods like TukeyHSD don't 
			#generalize well to high-dimensional datasets like microarrays.
			
			#The statistic fit2$F and the corresponding fit2$F.p.value combine the three pair-wise comparisons into one F
			#-test.  This is equivalent to a one-way ANOVA for each gene except that the residual mean squares have been moderated between genes.
			#To find genes which vary between any comparisons, look for genes with small p-values.  To find the top 50 genes:
			#topdeg = topTableF(fit2, number=50)
			#print(head(topdeg))
			
			#Note that toptable technically does not filter, it only sorts as per p-value
			#The outcome of each hypothesis test can be assigned using
			resultsdt <- decideTests(fit2, adjust.method = "none" )
			#print(resultsdt)
			if(ncol(resultsdt) > 3){
			    resultsdt = resultsdt[, 1:3]
			}
			
			#A Venn diagram showing numbers of genes significant in each comparison can be obtained from
			pdf(paste("Analys", indxg, "no-fdr-limma-first-3-comp-venn.pdf" , sep='-'))
			vennDiagram(resultsdt)
			dev.off()
			

			# get the results for the required input gene list
			outf1 = topTable(fit2[lgenes,],sort="none",n=Inf, adjust="BH") # in order of lgenes
			#print(head(outf1))
			out = outf1[,c("F" , "P.Value")]
			colnames(out) = c(paste("Analys", indxg, "limma F" , sep=' ') , paste("Analys", indxg, "limma F.P.Value" , sep=' '))

			geneName = rownames(out)
			out = cbind(geneName, out)
            
            # add the base statistics and the t test results		
			for (contarst in 1:length(comparisonsargs)){
			    indxg_ = paste(c(indxg, contarst), collapse='-') # to keep track of the comparisons in each analysis
		    	
			    tmpout = topTable(fit2[lgenes,],coef=contarst,sort="none",n=Inf, adjust="BH") # in order of lgenes
			    tmpout = tmpout[, -2] # get rid of the Average express across ALL samples (c1 and c2 category)
			    colnames(tmpout) = c(paste("Analys", indxg_, comparMethod,comparisonsargs[contarst],"limmaslog2FC=meanA-meanB",sep=" "), paste("Analys", indxg_, comparMethod,comparisonsargs[contarst],"limma t",sep=" "), paste("Analys", indxg_, comparMethod,comparisonsargs[contarst],"pvalue",sep=" "), paste("Analys", indxg_, comparMethod,comparisonsargs[contarst],"FDR",sep=" "), paste("Analys", indxg_, comparMethod,comparisonsargs[contarst],"log_odds_deg",sep=" "))
			    #print(head(tmpout))
			    
			    c1_c2 = unlist(strsplit(unlistedres[contarst], "_vs_"))
			    c1_ = c1_c2[1]
			    c2_ = c1_c2[2]
			    out2 = sapply(lgenes, GetSimpleStats, subsetexpressionData, c1_, c2_, dict)
			    out2 = t(out2)
			    colnames(out2) = c(paste("Analys", indxg_, comparMethod,"Mean", c1_, sep=" "), paste("Analys", indxg_, comparMethod,"Mean", c2_, sep=" "), paste("Analys", indxg_, comparMethod,"Median", c1_, sep=" "), paste("Analys", indxg_, comparMethod,"Median", c2_, sep=" "), paste("Analys", indxg_, comparMethod,"GeoMean", c1_, sep=" "), paste("Analys", indxg_, comparMethod,"GeoMean", c2_, sep=" "))
			    #print(head(out2))
		    	
			    out = cbind(out, out2)
			}
			#print(head(out))

	    
	    } else { # for two categories
		    cont.matrix <- makeContrasts(MUvsWT=MU-WT, levels=design)
		    print(cont.matrix)
		    fit2 <- contrasts.fit(fit, cont.matrix)
		    fit2 <- eBayes(fit2)
		    ## keep only the genes of interest # https://support.bioconductor.org/p/23611/
		    ##out1 = topTable(fit2,sort="none",n=Inf, adjust="BH")
		    out1 = topTable(fit2[lgenes,],sort="none",n=Inf, adjust="BH") # in order of lgenes
		    
		    out1 = out1[, -2] # get rid of the Average express across ALL samples (c1 and c2 category)
		    # https://www.biostars.org/p/100460/
		    colnames(out1) = c(paste("Analys", indxg, comparMethod,c1,"vs",c2,"limmaslog2FC=meanA-meanB",sep=" "), paste("Analys", indxg, comparMethod,c1,"vs",c2,"limma t",sep=" "), paste("Analys", indxg, comparMethod,c1,"vs",c2,"pvalue",sep=" "), paste("Analys", indxg, comparMethod,c1,"vs",c2,"FDR",sep=" "), paste("Analys", indxg, comparMethod,c1,"vs",c2,"log_odds_deg",sep=" "))
		    #print(head(out1))
		    
		    out2 = sapply(lgenes, GetSimpleStats, subsetexpressionData, c1, c2, dict)
		    out2 = t(out2)
		    
		    colnames(out2) = c(paste("Analys", indxg, comparMethod,"Mean", c1, sep=" "), paste("Analys", indxg, comparMethod,"Mean", c2, sep=" "), paste("Analys", indxg, comparMethod,"Median", c1, sep=" "), paste("Analys", indxg, comparMethod,"Median", c2, sep=" "), paste("Analys", indxg, comparMethod,"GeoMean", c1, sep=" "), paste("Analys", indxg, comparMethod,"GeoMean", c2, sep=" "))
		    #print(head(out2))
		    geneName = rownames(out2)
		    out = cbind(geneName, out2)
		    
	    }
	    #print(head(out))
        
    } else {
	    # calculate for each gene
	    out = sapply(lgenes, CalcCom, expressionData, c1, c2, dict, comparMethod, pairedd)

	    # the genes become row names; and the method, meanx, meany, and pvalue are the column names
	    out = t(out)

	    colnames(out) = c(paste("Analys", indxg, "Mean", c1, sep=" "), paste("Analys", indxg, "Mean", c2, sep=" "), paste("Analys", indxg, "Median", c1, sep=" "), paste("Analys", indxg, "Median", c2, sep=" "), paste("Analys", indxg, "GeoMean", c1, sep=" "), paste("Analys", indxg, "GeoMean", c2, sep=" "))

	    geneName = rownames(out)
	    out = cbind(geneName,out)

    }
    
    out
}



#==================================================================================================================
#		Condition 6
#==================================================================================================================
# function : calculate comparison of correlations between categories for a gene pair
# input:
#    	  1. pair: for which to calculate the comparison of correlation
#         2. expression data
#	  3 & 4. the category on which to perform comparison of correlations between measurements
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for correlation
#         7. paired or unpaired. For now only unpaired is used (independent groups)
#
# output : 
#	vector of length 3: correlation in category 1, correlation in category 2, pvalue of comparing correlations in categories

CalcComparisCor = function(pair, edata, c1, c2, dict, correlMethod, pairedd){

    # calculate correlation for categ1
    ca = CalcCor(pair, edata, c1, dict, correlMethod)
    N_a = CalcNCor(pair, edata, c1, dict, correlMethod)
    # calculate correlation for categ2
    cb = CalcCor(pair, edata, c2, dict, correlMethod)
    N_b = CalcNCor(pair, edata, c2, dict, correlMethod)

    searchString_pattern = ''
    if(grepl('pe',as.vector(correlMethod), ignore.case = TRUE)){
        searchString_pattern = 'cor'
    }else if(grepl('sp',as.vector(correlMethod), ignore.case = TRUE)){
        searchString_pattern = 'rho'
    }

    #print(ca)
    #print(cb)
    r_12 = NA
    r_34 = NA

    
    if(grepl(searchString_pattern, ca)){
        r_12 = as.numeric(as.character(ca[searchString_pattern,1]))
    }

    if(grepl(searchString_pattern, cb)){
        r_34 = as.numeric(as.character(cb[searchString_pattern,1]))
    }

    #print(N_a)
    #print(N_b)
    #print("R12")
    #print(r_12)
    #print("R34")
    #print(r_34)
    
    # made the results in vector format
    fca = as.vector(t(ca)) 
    fcb = as.vector(t(cb))


    if(grepl(searchString_pattern, ca) && grepl(searchString_pattern, cb)){

      # compare correlations # sometimes this test gives NA when the two correlations are same and also obtained from the same sample size
      res = r.test(n = N_a, r12 = r_12, r34 = r_34, n2 = N_b)
      #print(res)
      #print(toString(res))
      #print("\n")
      #print(toString(res[4]))
    
      # same results as the methods from the DAPfinder and Reverse engg paper
      # Convert correlations to z-scores 
      #z1 = fisherz(r_12)
      #z2 = fisherz(r_34)
      # Calculate vector of t-tests to compare correlations between classes
      #fisher = (z1 - z2) / sqrt((1/(N_a - 3)) + (1/(N_b - 3))) 
      # Calculate raw p-values
      #pv.dif.cor = 2*pt(-abs(fisher),Inf)
      #print(pv.dif.cor)

      #ret_result = c(r_12, r_34, res$p)  # OR
      ret_result = c( as.vector(fca[-1]) , as.vector(fcb[-1]) , res$p)
      #print(ret_result)
    } else {
        # because you still expect the correlations to be the same
        #ret_result = c(r_12, r_34, NA)    # OR
        ret_result = c( as.vector(fca[-1]) , as.vector(fcb[-1]) , NA)
        #print("here")
    }
  
    ret_result

}


#==================================================================================================================
#		Condition 6
#==================================================================================================================
# function : calculate comparison of correlations between categories,pvalue,fdr for all pairs
# input:
#    	  1. pairs: for which to calculate the comparison of correlation
#         2. expression data
#	  3 & 4. the category on which to perform comparison of correlations between measurements
#	  5. the dictionary containing mapping of category and samples which belong to it
#         6. method for correlation
#         7. paired or unpaired. For now only unpaired is used (independent groups)
#
# output : 
#	a table containing pair information, correlation in category 1, correlation in category 2, pvalue of comparing correlations in categories,fdr

calculateComparisCorrelation = function (pairs, expressionData, c1, c2, dict, correlMethod, pairedd, indxg){

    # calculate for each pair
    out = apply(pairs, 1, CalcComparisCor, expressionData, c1, c2, dict, correlMethod, pairedd)

    # add the pair as the column name
    colnames(out) = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")

    # the pairs become row names; and the corr coeff and pvalue are the column names
    out = t(out)

    # method used
    method_label = correlMethod
    # append method used to column labels
    colnames(out) = c( paste("Analys", indxg, c1, method_label,"Coefficient",sep=" "), paste("Analys", indxg, c1, method_label, "CategCorPVal",sep=" "), 
    		       paste("Analys", indxg, c2, method_label,"Coefficient",sep=" "), paste("Analys", indxg, c2, method_label, "CategCorPVal",sep=" "), 
    		       paste("Analys", indxg, method_label,c1,"vs",c2,"pvalue",sep=" "))

    # calculate FDR
    FDR = p.adjust(out[,colnames(out)[grep("pvalue",colnames(out))]],method="fdr")
    oldColnames = colnames(out)
    out = cbind(out, FDR)
    colnames(out) = c(oldColnames,paste("Analys", indxg, method_label,c1,"vs",c2,"FDR",sep=" "))
    pairName = rownames(out)
    out = cbind(pairName,out)
    out
}


#==================================================================================================================
#		Condition 7
#==================================================================================================================
# function : calculate comparison of (delta) correlation (correl(c1-c2) vs correl(c3-c4)) for gene pair using 4 categories
# input:
#    	  1. pairs: for which to calculate the correlation
#         2. expression data
#	  3 - 6. the category on which to perform comparison of delta correlations (correl(c1-c2) vs correl(c3-c4)) between measurements
#	  7. the dictionary containing mapping of category and samples which belong to it
#         8. method for correlation
# output : 
#	vector of length 3: correlation in delta1 (i.e. category1-category2), correlation in delta2 (i.e. category3-category4), pvalue of comparing delta correlations

CalcDeltaComparisCor = function(pair, edata, c1, c2, c3, c4, dict, correlMethod){

    #calculate (delta) correlation coefficient of a pair of genes for c1 and c2
    ca = CalcDeltaCor(pair, edata, c1, c2, dict, correlMethod)
    N_a = CalcNDeltaCor(pair, edata, c1, c2, dict, correlMethod)
    #calculate (delta) correlation coefficient of a pair of genes for c3 and c4
    cb = CalcDeltaCor(pair, edata, c3, c4, dict, correlMethod)
    N_b = CalcNDeltaCor(pair, edata, c3, c4, dict, correlMethod)

    searchString_pattern = ''
    if(grepl('pe',as.vector(correlMethod), ignore.case = TRUE)){
        searchString_pattern = 'cor'
    }else if(grepl('sp',as.vector(correlMethod), ignore.case = TRUE)){
        searchString_pattern = 'rho'
    }

    r_12 = NA
    r_34 = NA

    if(grepl(searchString_pattern, ca)){
        r_12 = as.numeric(as.character(ca[searchString_pattern,1]))
    }

    if(grepl(searchString_pattern, cb)){
        r_34 = as.numeric(as.character(cb[searchString_pattern,1]))
    }

    #print(N_a)
    #print(N_b)
    #print("R12")
    #print(r_12)
    #print("R34")
    #print(r_34)

    # made the results in vector format
    fca = as.vector(t(ca)) 
    fcb = as.vector(t(cb))

    if(grepl(searchString_pattern, ca) && grepl(searchString_pattern, cb)){

      # compare correlations
      res = r.test(n = N_a, r12 = r_12, r34 = r_34, n2 = N_b)
      #print(toString(res))
      #print("\n")
      #print(toString(res[4]))
    
      # same results as the methods from the DAPfinder and Reverse engg paper
      # Convert correlations to z-scores 
      #z1 = fisherz(r_12)
      #z2 = fisherz(r_34)
      # Calculate vector of t-tests to compare correlations between classes
      #fisher = (z1 - z2) / sqrt((1/(N_a - 3)) + (1/(N_b - 3))) 
      # Calculate raw p-values
      #pv.dif.cor = 2*pt(-abs(fisher),Inf)
      #print(pv.dif.cor)

      #ret_result = c(r_12, r_34, res$p) # OR
      ret_result = c( as.vector(fca[-1]) , as.vector(fcb[-1]) , res$p)
    } else {
        # because you still expect the correlations to be the same
        #ret_result = c(r_12, r_34, NA) # OR
        ret_result = c( as.vector(fca[-1]) , as.vector(fcb[-1]) , res$p)
    }

    ret_result

}


#==================================================================================================================
#		Condition 7
#==================================================================================================================
# function : calculate comparison of (delta) correlation (correl(c1-c2) vs correl(c3-c4)),pvalue,fdr for all pairs using 4 categories
# input:
#    	  1. pairs: for which to calculate the correlation
#         2. expression data
#	  3 - 6. the category on which to perform comparison of delta correlations (correl(c1-c2) vs correl(c3-c4)) between measurements
#	  7. the dictionary containing mapping of category and samples which belong to it
#         8. method for correlation
# output : 
#	a table containing pair information, correlation(c1-c2), correlation(c3-c4), comparison of delta correlations, pvalue,fdr

calculateComparisDeltaCorrelation = function (pairs, expressionData, c1, c2, c3, c4, dict, correlMethod, indxg){

    # calculate for each pair
    out = apply(pairs, 1, CalcDeltaComparisCor, expressionData, c1, c2, c3, c4, dict, correlMethod)

    # add the pair as the column name
    colnames(out) = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")

    # the pairs become row names; and the corr coeff and pvalue are the column names
    out = t(out)


    # method used
    method_label = correlMethod

    # append method used to column labels
    colnames(out) = c(paste("Analys", indxg, c1, "minus" , c2, method_label,"Coefficient",sep=" "), paste("Analys", indxg, c1, "minus" , c2, method_label, "CategCorPVal",sep=" "), 
    		      paste("Analys", indxg, c3, "minus" , c4, method_label,"Coefficient",sep=" "), paste("Analys", indxg, c3, "minus" , c4, method_label, "CategCorPVal",sep=" "), 
    		      paste("Analys", indxg, method_label, c1, "minus" , c2, "vs", c3, "minus" , c4, "pvalue",sep=" "))

    # calculate FDR
    FDR = p.adjust(out[,colnames(out)[grep("pvalue",colnames(out))]],method="fdr")
    oldColnames = colnames(out)
    out = cbind(out, FDR)
    colnames(out) = c(oldColnames,paste("Analys", indxg, method_label,c1, "minus" , c2, "vs", c3, "minus" , c4,"FDR",sep=" "))
    pairName = rownames(out)
    out = cbind(pairName,out)
    out
}


#==================================================================================================================
#		Condition 2
#==================================================================================================================
# function : calculate comparison of deltas (i.e. c1-c2 vs c3-c4) for gene using 4 categories
# input:
#    	  1. gene: for which to calculate the comparison of deltas
#         2. expression data
#	  3 - 6. the categories to perform delta and then compare
#	  7. the dictionary containing mapping of category and samples which belong to it
#         8. method for comparison (t test, Man Whitney, etc.)
# output : 
#	vector of length 3: mean of c1-c2, mean of c3-c4, pvalue

CalcComDelta = function(lgene, edata, c1, c2, c3, c4, dict, comparMethod){
    idxs1 = as.vector(dict[[c1]])
    idxs2 = as.vector(dict[[c2]])

    idxs3 = as.vector(dict[[c3]])
    idxs4 = as.vector(dict[[c4]])

    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste(lgene, c1, "-", c2, "vs" , c3, "-", c4, sep=" ___ ")
    cat (calculating)
    cat ("\n")

    c1 = numericizeVector(as.vector(edata[lgene, idxs1]))
    c2 = numericizeVector(as.vector(edata[lgene, idxs2]))

    c3 = numericizeVector(as.vector(edata[lgene, idxs3]))
    c4 = numericizeVector(as.vector(edata[lgene, idxs4]))

    
    ca = ''
    cb = ''
    if((sum(!is.na(c1)) >= 3) & (sum(!is.na(c2)) >= 3)){
    	# calculate delta
    	#ca = c1-c2

    	# calculate delta normalized to starting value
    	ca = (c1-c2)/c2

    } else {
        ca[1:3] <- NA  # since we expect this to NA
    }

    if((sum(!is.na(c3)) >= 3) & (sum(!is.na(c4)) >= 3)){
    	# calculate delta
    	#cb = c3-c4

    	# calculate delta normalized to starting value
    	cb = (c3-c4)/c4

    } else {
        cb[1:3] <- NA  # since we expect this to NA
    }


    outLine = ''
    p = ''


    mean_ca = NA
    mean_cb = NA
    
    median_ca = NA
    median_cb = NA

    # there are at least 3 samples with values (i.e. without NA)
    if(sum(!is.na(ca)) >= 3){
        # on deltas
        mean_ca = mean(ca, na.rm=TRUE)
        median_ca = median(ca, na.rm=TRUE)
    } 

    if(sum(!is.na(cb)) >= 3){
        # on deltas
        mean_cb = mean(cb, na.rm=TRUE)
        median_cb = median(cb, na.rm=TRUE)
    } 
   
    outLine = as.matrix(c(mean_ca, mean_cb, median_ca, median_cb)) 
    outLine

}


#==================================================================================================================
#		Condition 2
#==================================================================================================================
# function : calculate comparison of deltas (i.e. c1-c2 vs c3-c4),pvalue,fdr for all genes using 4 categories
# input:
#    	  1. list of genes: for which to calculate the comparison of deltas
#         2. expression data
#	  3 - 6. the categories to perform delta and then compare
#	  7. the dictionary containing mapping of category and samples which belong to it
#         8. method for comparison (t test, Man Whitney, etc.)
# output : 
#	a table containing gene information, means & medians of normalized deltas for (c1-c2)/c2 and (c3-c4)/c4, fold change ((c1/c2):(c3/c4)), pvalue,fdr

calculateComparisonDelta = function (lgenes, expressionData, c1, c2, c3, c4, dict, comparMethod, indxg){
    # calculate for each gene
    out = sapply(lgenes, CalcComDelta, expressionData, c1, c2, c3, c4, dict, comparMethod)

    #print(as.data.frame(out))
       
    # the genes become row names; and the method, meanx, meany, and pvalue are the column names
    out = t(out)
    
    colnames(out) = c(paste("Analys", indxg, "Mean", c1, "minus", c2, sep=" "), paste("Analys", indxg, "Mean", c3, "minus" , c4, sep=" "), paste("Analys", indxg, "Median", c1, "minus", c2, sep=" "), paste("Analys", indxg, "Median", c3, "minus" , c4, sep=" "))

    geneName = rownames(out)
    out = cbind(geneName,out)
    out
}



