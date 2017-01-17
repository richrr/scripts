# usage: Rscript ~/Morgun_Lab/richrr/scripts/R/calc-stats-from-jackknife.R ./ corr 4,5 parallel
# /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/jackknife-25-samples-comparisons/cohort1/test
# Rscript ~/Morgun_Lab/richrr/scripts/R/calc-stats-from-jackknife.R ./ comp 1 parallel

library(data.table)

args = commandArgs(trailingOnly=TRUE)
path = args[1]
patt = args[2]

grepFileName = ''

scratchDir = './scratch/'

uniq_col_name = "pairName"
metric_col = "Coefficient"

if(patt == "comp") { 
    uniq_col_name = "geneName" 
    metric_col = "FolChMedian"
}

analys_to_do = unlist(strsplit(args[3], ","))

ofname = paste(c("median-across-jackknife", patt, analys_to_do), collapse='_')



	cmd=paste('rm' , '-rf', scratchDir, sep=' ')
	print(cmd)
    d1 <- system(cmd, intern = TRUE)


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
# not needed if using grep
# this also does not work for multiple analysis, since you no longer have the columns for the next analysis
cfun <- function(...){
    arguments <- list(...)
    res = rbindlist(arguments)
    #print(res)

    # this increases with every process
    #pattrn = paste(res$p, collapse='|')
    # use the last element
    pattrn = res$p[length(res$p)] 
    
    pattrn = paste("'", pattrn, "'", sep='')
    cmd=paste('LC_ALL=C', 'fgrep', "-v", pattrn , grepFileName, "> filename2 ; mv filename2", grepFileName, sep=' ')
	#print(cmd)
    d1 <- system(cmd, intern = TRUE)
    
    # delete the finished pairs from the big file to save read time    
    # http://stackoverflow.com/questions/5410757/delete-lines-in-a-text-file-that-containing-a-specific-string
    #for(e in res$p){
    #   cmd=paste('LC_ALL=C', 'grep', "-v", e , grepFileName, "> filename2 ; mv filename2", grepFileName, sep=' ')
	#   d1 <- system(cmd, intern = TRUE)
    # }

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
        library(reshape)
             
	#cores=detectCores()
	#cl <- makeCluster(50) # or cores[1]-1
	#registerDoParallel(cl)
	cl <- makeCluster(50, outfile="")  # print output on screen
	registerDoSNOW(cl) # print output on screen
	
        

	final_out_df = data.frame(pairs)
	colnames(final_out_df) = c(uniq_col_name)
	
	# grab all the needed files once:
	grepFileName <<- paste(ofname, "-wholedf-parallel.csv", sep='')
	
    if(FALSE){
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
	write.csv(wholedf, grepFileName, row.names=FALSE, quote=FALSE)
	rm(wholedf)  # empty the variable to reduce size
	}
	
	# keep the header from first file for later use
	CMD = paste('head', '-1', files[1], sep=' ')
    HEADER <- system(CMD, intern = TRUE)
    
    # using cat to save time   
	CMD = paste(c('cat', files), collapse=' ')
	CMD = paste(CMD, '|' , 'sort', '>' , grepFileName, sep=' ')
	CREATEFILE <- system(CMD, intern = TRUE)
    print("Starting median calculation") 
        
    # the number of columns in this file
    cmD = paste("awk", "-F','", "'{ print NF }'" , grepFileName, '|' , 'sort', '|' , 'uniq',  sep=' ')
    numbCols = as.integer(system(cmD, intern = TRUE))
    print(numbCols)
    
    
    # create a tmp folder to split the big file into small files
    cmd = paste('mkdir' , scratchDir , sep=' ')
    CREATEDIR <- system(cmd, intern = TRUE)
    cmd = paste('split -a 5 -d -l 100000', grepFileName, scratchDir, sep=' ')
    SPLITCMD <- system(cmd, intern = TRUE)
    
    splitfiles = list.files(path = scratchDir)
    #print(splitfiles)
    #print(typeof(splitfiles))
    #print(class(splitfiles))
    
    mapdf = data.frame(c("NA"), c("NA"))
    colnames(mapdf) = c("vals", "fname")
    for(f in splitfiles){
        #print(f)
        tmp = paste(scratchDir, f, sep='')
        cmd=paste('cut -d, -f1', tmp, '|' , 'sort', '|' , 'uniq' , sep=' ')
        vals = system(cmd, intern = TRUE)
        vals <- gsub('\"', '', vals)
        tmpdf = as.data.frame(vals)
        tmpdf$fname = f
        mapdf = rbind(mapdf, tmpdf)
    }
    #print(mapdf)
    

        
      ### read this file completely, analyze pair one by one and delete the rows every time you are done with it (longest time and reducing memory)
      ### read this file completely, analyze pair in parallel and delete the rows every time you are done with it (not sure how multiple it works with multiple processors deleting simultaneosuly)
          ### how about one processor deletes after the combine function
          ### using this in memory needs Rmpi and couldn't get it to work
      ### read this file in chunks to reduce memory and analyze each pair one by one (longest time, medium memory)
      ### do not read the file, parallelize unix's grep to get the rows you want from a file (medium time, lowest memory)
          ### this works
      
	  ### to do ###
	  ## make the code faster
	     ### split into files each with 100K lines, use numeric suffix, allows upto 100000 files, and keep track of the unique gene names in each file
	        #### split -a 5 -d -l 100000 filename prefix
	     ### only search the files which has the pairs
	  ## continue where left off
	  
	for(a in analys_to_do){  # Although this can be parallelized, it doesn't save much time since it still might need the pairs to be analyzed sequentially

		analys = paste("Analys ", a, " ", sep='')

		# extract the correct column names from the first file
		#print(files[1])
		tmp_read = read.csv(files[1], header=T, row.names=1, check.names=F) # read file
		
		tmp_selectanalyCols = grep( analys, colnames(tmp_read), value=T) # select the required analysis
		tmp_selectpvalCols = grep( "pvalue", tmp_selectanalyCols, value=T) # select the required pval columns
		tmp_selectmetricCols = grep(metric_col , tmp_selectanalyCols, value=T) # select the required metric columns
        rm(tmp_read)

		# the BEST loop to parallelize
		out_df = foreach(p=sort(pairs), .combine=rbind) %dopar% { 
		#out_df = foreach(p=sort(pairs), .combine='cfun') %dopar% { 
		    print(p)
		    
		    #### in case there is more than one analysis in the big file pick as per analys number
		    ### try fgrep to exit after 1000 matches; plus it may be faster
		    ### http://stackoverflow.com/questions/13913014/grepping-a-huge-file-80gb-any-way-to-speed-it-up https://groups.google.com/forum/#!topic/comp.lang.awk/gngk-epT9Ug
		    #cmd=paste('LC_ALL=C', 'fgrep', "-m 1000", p , grepFileName, sep=' ')

                    scractfilename = ''
                    if(length(as.vector(mapdf[which(mapdf[, "vals"] == p), "fname"])) == 1){
		        scractfilename = paste(c(scratchDir, as.vector(mapdf[which(mapdf[, "vals"] == p), "fname"])), collapse="")
		    } else {
		        print("Redo the calculation for this id since there seem to be one or more files for this id.")
		    }
		    
		    #print(scractfilename)
		    pattren = paste("'", p, "'", sep='')		    
		    cmd=paste('LC_ALL=C', 'fgrep', "-m 1000", pattren , scractfilename, sep=' ')
		    t1 <- system(cmd, intern = TRUE)
		    
		    # add header to the results vector
		    t1 <- c(HEADER, t1)
		    t1 <- gsub('\"', '', t1) 
		    
		    
		    # conver the char vector to data frame
		    df = reshape::colsplit(t1, split=",", c(1:numbCols))  # or # #df <- data.frame(do.call(rbind, strsplit(t1, ",", fixed=TRUE)))
		    colnames(df) = unlist(df[1, ]) # the first row will be the header
		    df = df[-1, ] # drop the row that you don't need
		    
		    # convert the non-name cols to numeric
		    dfCol1 = as.data.frame(df[,1])
		    colnames(dfCol1) = uniq_col_name
		    newdf = as.data.frame(lapply(df[ , -1], function(x) as.numeric(as.character(x))), check.names=F)
		    df = cbind(dfCol1, newdf)
		    #print(df)
		    #print(colnames(df))
		    
            # keep columns for the required analysis
		    selectCols = grep( analys, colnames(df), value=T) # select the required analysis
		    df = df[which(df[, uniq_col_name] == p), c( uniq_col_name, selectCols)]
		    
		    # OR if everything was loaded in wholedf
		    #selectCols = grep( analys, colnames(wholedf), value=T) # select the required analysis
		    #df = wholedf[which(wholedf[, uniq_col_name] == p), c( uniq_col_name, selectCols)]
	    
		        
		    pvalCol = grep("pvalue" , colnames(df), value=T)
		    median_pval = median(df[,pvalCol], na.rm = T) # get the median p-value
		    
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
		    print(res)
		    res
		    
		}
                
		colnames(out_df) = c(uniq_col_name , tmp_selectmetricCols, tmp_selectpvalCols)
		final_out_df = merge(final_out_df , out_df, by=uniq_col_name)
	}
	write.csv(final_out_df, paste(ofname, "-parallel.csv", sep=''), row.names=FALSE, quote=FALSE)
        
	#stop cluster
	stopCluster(cl)
	
	cmd=paste('rm' , '-f', grepFileName, sep=' ')
	#print(cmd)
    d1 <- system(cmd, intern = TRUE)

	cmd=paste('rm' , '-rf', scratchDir, sep=' ')
	print(cmd)
    d1 <- system(cmd, intern = TRUE)

}
