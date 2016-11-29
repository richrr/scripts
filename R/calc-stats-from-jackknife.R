# usage: Rscript ~/Morgun_Lab/richrr/scripts/R/calc-stats-from-jackknife.R ./ corr 4,5 parallel


args = commandArgs(trailingOnly=TRUE)
path = args[1]
patt = args[2]

analys_to_do = unlist(strsplit(args[3], ","))

files = list.files(path = path, pattern = patt, recursive = T)

# get the pairs from any file
fi = read.csv(files[1], header=T, check.names=F)
pairs = as.vector(fi$pairName)


median_indx <- function(x) {
  x <- sort(x, na.last = NA)
  n <- length(x)
  m <- (n+1)/2
  if (floor(m) != m) {  # when n is even
    l <- m-1/2; u <- m+1/2
  } else {  # when n is odd
    l <- m; u <- m
  }
  return(c(l,u))
}


# for each analysis, for each pair, loop over the files and get the correl and pvalue

if(args[4] == "serial"){
	final_out_df = data.frame(pairs)
	colnames(final_out_df) = c("pairName")

	for(a in analys_to_do){

		analys = paste("Analys ", a, " ", sep='')
		out_df = data.frame()
		for (p in pairs){

		    print(p)
		    df = data.frame()
		    for (f in files){
			lines = read.csv(f, header=T, row.names=1, check.names=F) # read file
			selectCols = grep( analys, colnames(lines), value=T) # select the required analysis
			df_ = lines[p, selectCols] # keep the required row and columns
			df = rbind(df, df_)
		
		    }
		    #print(head(df))
		    
		    pvalCol = grep("pvalue" , colnames(df), value=T)
		    median_pval = median(df[,pvalCol], na.rm = T) # get the median p-value
		    
		    # remove NA from ordering
		    sorted_df = df[order(df[,pvalCol], na.last = NA),]
		    #print(nrow(sorted_df))

		    # get the coefficient corresponding to the median p value
		    med_indx = median_indx(sorted_df[,pvalCol])
		    l = med_indx[1]
		    u = med_indx[2]
		    coeffCol = grep("Coefficient" , colnames(df), value=T)

		    median_coeff = NA
		    if( l != u){
			#print(l)
			#print(u)
			v1 = NA
			v2 = NA
			if(l > 0) {v1 = sorted_df[,coeffCol][l]} # this could happen when all the vectors in that list were NAs and so we got empty vector in sorted_df
			if(u > 0) {v2 = sorted_df[,coeffCol][u]}
			#print(v1)
			#print(v2)
			median_coeff = (v1+v2)/2.0
		    } else {
			median_coeff = sorted_df[,coeffCol][l]
		    }
		    
		    # why do we need confidence interval? #print(t.test(sorted_df[,coeffCol]))
		    
		    out_df = rbind(out_df, data.frame(p , median_pval, median_coeff))
		    
		    
		    
		}
		colnames(out_df) = c("pairName", paste(analys, "pvalue", sep='') , paste(analys, "Coefficient", sep=''))
		#print(out_df)
		final_out_df = merge(final_out_df , out_df)
	}
	write.csv(final_out_df, "median-across-jackknife.csv", row.names=FALSE, quote=FALSE)
} else {

	library(foreach)
	library(doParallel)

	cores=detectCores()
	cl <- makeCluster(50) # or cores[1]-1
	registerDoParallel(cl)

	final_out_df = data.frame(pairs)
	colnames(final_out_df) = c("pairName")

	for(a in analys_to_do){  # Although this can be parallelized, it doesn't save much time since it still might need the pairs to be analyzed sequentially

		analys = paste("Analys ", a, " ", sep='')

		# the BEST loop to parallelize
                out_df = foreach(p=pairs, .combine=rbind) %dopar% { 
		    print(p)
		    df = data.frame()
		    for (f in files){ # SHOULD NOT be parallelized since you need all the information from all the files to calculate the median
			lines = read.csv(f, header=T, row.names=1, check.names=F) # read file
			selectCols = grep( analys, colnames(lines), value=T) # select the required analysis
			df_ = lines[p, selectCols] # keep the required row and columns
			df = rbind(df, df_)		
		    }
		    #print(head(df))
		    
		    pvalCol = grep("pvalue" , colnames(df), value=T)
		    median_pval = median(df[,pvalCol], na.rm = T) # get the median p-value
		    
		    # remove NA from ordering
		    sorted_df = df[order(df[,pvalCol], na.last = NA),]
		    #print(nrow(sorted_df))

		    # get the coefficient corresponding to the median p value
		    med_indx = median_indx(sorted_df[,pvalCol])
		    l = med_indx[1]
		    u = med_indx[2]
		    coeffCol = grep("Coefficient" , colnames(df), value=T)

		    median_coeff = NA
		    if( l != u){
			#print(l)
			#print(u)
			v1 = NA
			v2 = NA
			if(l > 0) {v1 = sorted_df[,coeffCol][l]} # this could happen when all the vectors in that list were NAs and so we got empty vector in sorted_df
			if(u > 0) {v2 = sorted_df[,coeffCol][u]}
			#print(v1)
			#print(v2)
			median_coeff = (v1+v2)/2.0
		    } else {
			median_coeff = sorted_df[,coeffCol][l]
		    }
		    
		    # why do we need confidence interval? #print(t.test(sorted_df[,coeffCol]))
		    
		    res = data.frame(p , median_pval, median_coeff)
		    res
		    
		}

		colnames(out_df) = c("pairName", paste(analys, "pvalue", sep='') , paste(analys, "Coefficient", sep=''))
		#print(out_df)
		final_out_df = merge(final_out_df , out_df)
	}
	write.csv(final_out_df, "median-across-jackknife-parallel.csv", row.names=FALSE, quote=FALSE)

	#stop cluster
	stopCluster(cl)

}
