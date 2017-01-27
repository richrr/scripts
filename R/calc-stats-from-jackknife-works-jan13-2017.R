# usage: Rscript ~/Morgun_Lab/richrr/scripts/R/calc-stats-from-jackknife.R ./ corr 4,5 parallel
# /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/jackknife-25-samples-comparisons/cohort1/test
# Rscript ~/Morgun_Lab/richrr/scripts/R/calc-stats-from-jackknife.R ./ comp 1 parallel

library(data.table)

args = commandArgs(trailingOnly=TRUE)
path = args[1]
patt = args[2]

uniq_col_name = "pairName"
metric_col = "Coefficient"

if(patt == "comp") { 
    uniq_col_name = "geneName" 
    metric_col = "FolChMedian"
}

analys_to_do = unlist(strsplit(args[3], ","))

ofname = paste(c("median-across-jackknife", patt, analys_to_do), collapse='_')


files = list.files(path = path, pattern = patt, recursive = T)


# get the pairs from any file
fi = read.csv(files[1], header=T, check.names=F)
pairs = as.vector(fi[,uniq_col_name])


median_indx <- function(x) {
  # if you resort it here, and for some reason if the sorting does not match the input,
  # then the indexes returned do not match with the input order
  #x <- sort(x, na.last = NA)  
  n <- length(x)
  m <- (n+1)/2
  if (floor(m) != m) {  # when n is even
    l <- m-1/2; u <- m+1/2
  } else {  # when n is odd
    l <- m; u <- m
  }
  return(c(l,u))
}

# custom function that performs rbind, then deletes the pairs from the wholedf to reduce memory
cfun <- function(...){
    arguments <- list(...)
    res = rbindlist(arguments)
    print(res)

    # delete the finished pairs from the wholedf to clear memory    
    for(e in res$p){
        wholedf <<- wholedf[!wholedf[, uniq_col_name] == e, ]
    }

    res

}

# for each analysis, for each pair, loop over the files and get the correl and pvalue

if(args[4] == "serial"){
	final_out_df = data.frame(pairs)
	colnames(final_out_df) = c(uniq_col_name)

	for(a in analys_to_do){
                #print(a)
                #quit()
		analys = paste("Analys ", a, " ", sep='')

                # extract the correct column names from the first file
                tmp_read = read.csv(files[1], header=T, row.names=1, check.names=F) # read file
                tmp_selectanalyCols = grep( analys, colnames(tmp_read), value=T) # select the required analysis
		tmp_selectpvalCols = grep( "pvalue", tmp_selectanalyCols, value=T) # select the required pval columns
		tmp_selectmetricCols = grep(metric_col , tmp_selectanalyCols, value=T) # select the required metric columns


		out_df = data.frame()
		for (p in pairs){
		    print(p)
		    df = data.frame()
		    for (f in files){
			lines = read.csv(f, header=T, row.names=1, check.names=F) # read file
			#print(colnames(lines))
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

		    # get the coefficient or foldchange corresponding to the median p value
		    med_indx = median_indx(sorted_df[,pvalCol])
		    l = med_indx[1]
		    u = med_indx[2]
		    coeff_or_fc_col = grep(metric_col , colnames(df), value=T)

                    # in case there is single fc or correlation column
                    if(length(coeff_or_fc_col) == 1){  
	    		    median_coeff = NA
			    if( l != u){
				#print(l)
				#print(u)
				v1 = NA
				v2 = NA
				if(l > 0) {v1 = sorted_df[,coeff_or_fc_col][l]} # this could happen when all the vectors in that list were NAs and so we got empty vector in sorted_df
				if(u > 0) {v2 = sorted_df[,coeff_or_fc_col][u]}
				#print(v1)
				#print(v2)
				median_coeff = (v1+v2)/2.0
			    } else {
				median_coeff = sorted_df[,coeff_or_fc_col][l]
			    }
			    
			    # why do we need confidence interval? #print(t.test(sorted_df[,coeff_or_fc_col]))
			    out_df = rbind(out_df, data.frame(p , median_coeff, median_pval))
		    } else{ # in case there is more than one fc or correlation column. mainly written for compare correlations. not tested for other conditions
			    median_coeff = data.frame("NA")
			    for(group in coeff_or_fc_col)
			    {
			            median_coeff_tmp = NA
				    if( l != u){
					v1 = NA
					v2 = NA
					if(l > 0) {v1 = sorted_df[,group][l]} # this could happen when all the vectors in that list were NAs and so we got empty vector in sorted_df
					if(u > 0) {v2 = sorted_df[,group][u]}
					median_coeff_tmp = (v1+v2)/2.0					
				    } else {
					median_coeff_tmp = sorted_df[,group][l]
				    }
				    median_coeff = cbind(median_coeff, median_coeff_tmp)
			    }
			    median_coeff = median_coeff[, -1] # remove the temporarily added NA during intialization
			    
			    out_df = data.frame(p , median_coeff, median_pval)		    
		   
		    }
		    
		    
		}
		#colnames(out_df) = c(uniq_col_name , paste(analys, metric_col, sep=''), paste(analys, "pvalue", sep=''))
		colnames(out_df) = c(uniq_col_name , tmp_selectmetricCols, tmp_selectpvalCols)
		#print(out_df)
		final_out_df = merge(final_out_df , out_df, by=uniq_col_name)
	}
	write.csv(final_out_df, paste(ofname, ".csv", sep=''), row.names=FALSE, quote=FALSE)
} else {

	library(foreach)
	library(doParallel)
	library(doSNOW) # print output on screen
        
             
	#cores=detectCores()
	#cl <- makeCluster(50) # or cores[1]-1
	#registerDoParallel(cl)
	cl <- makeCluster(50, outfile="")  # print output on screen
	registerDoSNOW(cl) # print output on screen
	
        

	final_out_df = data.frame(pairs)
	colnames(final_out_df) = c(uniq_col_name)
	
	# grab all the needed files once:

	wholedf = data.frame()


	pattrn = paste( uniq_col_name , metric_col, "pvalue" , sep='|')

	#need all the information from all the files to calculate the median
	wholedf = foreach(f=files, .combine=rbind) %dopar% { 
			print(f)
			alllines = read.csv(f, header=T, check.names=F) # read file
			
			selCols = grep( pattrn, colnames(alllines), value=T) # select the required cols
			#print(selCols)
			lines = alllines[, selCols] 
			
			lines
	}
	#print(wholedf)
        
        # sort it before you print it
        wholedf = wholedf[with(wholedf, order(wholedf[,1])), ]
	write.csv(wholedf, paste(ofname, "-wholedf-parallel.csv", sep=''), row.names=FALSE, quote=FALSE)
        print("Starting median calculation") 
      ### read this file completely, analyze pair one by one and delete the rows every time you are done with it (longest time and reducing memory)
      ### read this file completely, analyze pair in parallel and delete the rows every time you are done with it (not sure how multiple it works with multiple processors deleting simultaneosuly)
          ### how about one processor deletes after the combine function
      ### read this file in chunks to reduce memory and analyze each pair one by one (longest time, medium memory)
      ### do not read the file, parallelize unix's grep to get the rows you want from a file (medium time, lowest memory)
      
	#if(FALSE){
	for(a in analys_to_do){  # Although this can be parallelized, it doesn't save much time since it still might need the pairs to be analyzed sequentially

		analys = paste("Analys ", a, " ", sep='')

		# extract the correct column names from the first file
		tmp_read = read.csv(files[1], header=T, row.names=1, check.names=F) # read file
		tmp_selectanalyCols = grep( analys, colnames(tmp_read), value=T) # select the required analysis
		tmp_selectpvalCols = grep( "pvalue", tmp_selectanalyCols, value=T) # select the required pval columns
		tmp_selectmetricCols = grep(metric_col , tmp_selectanalyCols, value=T) # select the required metric columns
                tmp_read = NULL

		# the BEST loop to parallelize
		#out_df = foreach(p=pairs, .combine=rbind) %dopar% { 
		out_df = foreach(p=pairs, .combine='cfun') %dopar% { 
		    print(p)
		    
		    selectCols = grep( analys, colnames(wholedf), value=T) # select the required analysis
		    
		    df = wholedf[which(wholedf[, uniq_col_name] == p), c( uniq_col_name, selectCols)]
		    #print(df)
		    
		    # if you edit the wholedf (delete rows) in parallel the row numbers might be off
		    # but it might be possible to delete the rows using the above condition instead
		    # of using specific row numbers
		    
		    # alternative is create the wholedf using only the cols you need
		    # this would save memory
		    
		    pvalCol = grep("pvalue" , colnames(df), value=T)
		    #print(df[,pvalCol])
		    median_pval = median(df[,pvalCol], na.rm = T) # get the median p-value
		    #print(median_pval)


		    # remove NA from ordering
		    sorted_df = df[order(df[,pvalCol], na.last = NA),]
		    #print(sorted_df)

		    # get the coefficient or foldchange corresponding to the median p value
		    med_indx = median_indx(sorted_df[,pvalCol])
		    l = med_indx[1]
		    u = med_indx[2]
		    coeff_or_fc_col = grep(metric_col , colnames(df), value=T)

		    res = ''
		    
            	    # in case there is single fc or correlation column
            	    if(length(coeff_or_fc_col) == 1){  
			    median_coeff = NA
			    if( l != u){
				#print(l)
				#print(u)
				v1 = NA
				v2 = NA
				if(l > 0) {v1 = sorted_df[,coeff_or_fc_col][l]} # this could happen when all the vectors in that list were NAs and so we got empty vector in sorted_df
				if(u > 0) {v2 = sorted_df[,coeff_or_fc_col][u]}
				#print(v1)
				#print(v2)
				median_coeff = (v1+v2)/2.0
			    } else {
				median_coeff = sorted_df[,coeff_or_fc_col][l]
			    }
			    
			    # why do we need confidence interval? #print(t.test(sorted_df[,coeff_or_fc_col]))
			    
			    res = data.frame(p , median_coeff, median_pval)
		    
		    } else{ # in case there is more than one fc or correlation column. mainly written for compare correlations. not tested for other conditions
			    median_coeff = data.frame("NA")
			    for(group in coeff_or_fc_col)
			    {
			            median_coeff_tmp = NA
				    if( l != u){
					v1 = NA
					v2 = NA
					if(l > 0) {v1 = sorted_df[,group][l]} # this could happen when all the vectors in that list were NAs and so we got empty vector in sorted_df
					if(u > 0) {v2 = sorted_df[,group][u]}
					median_coeff_tmp = (v1+v2)/2.0					
				    } else {
					median_coeff_tmp = sorted_df[,group][l]
				    }
				    median_coeff = cbind(median_coeff, median_coeff_tmp)
			    }
			    median_coeff = median_coeff[, -1] # remove the temporarily added NA during intialization
			    
			    res = data.frame(p , median_coeff, median_pval)		    
		   
		    }
		    res
		}
                
        #print(out_df)
		#colnames(out_df) = c(uniq_col_name , paste(analys, metric_col, sep=''), paste(analys, "pvalue", sep=''))
		colnames(out_df) = c(uniq_col_name , tmp_selectmetricCols, tmp_selectpvalCols)
		#print(out_df)
		#print(final_out_df)
		final_out_df = merge(final_out_df , out_df, by=uniq_col_name)
	}
	write.csv(final_out_df, paste(ofname, "-parallel.csv", sep=''), row.names=FALSE, quote=FALSE)
        #}
	#stop cluster
	stopCluster(cl)

}
