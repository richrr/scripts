# this files takes as input the "consistent" files from the per_analysis folder for deltas.
# the input file is assumed to be consistent for fc>0 and pval, so just 
# use this to process the file for delta consistency using the correct fc (1).


args = commandArgs(trailingOnly=TRUE)
#print(length(args))


check_delta_consistency = function(infile, total_numb_input_files, foldchVar = 'FolChMedian'){
	print(infile)
	df = read.csv( infile, header=TRUE,check.names=FALSE)
	rownames(df) = df[,1]
	outdf = df
	
	select_analys_cols = grep(foldchVar, colnames(df), value=T)
	df = df[,select_analys_cols]
	#print(head(df))
	cols_processed = 1
	rows_pass_consistency = c()
	while(cols_processed < ncol(df))  # loop over the merged dataset in steps of (number of input files)
	{
		subset_df = df[,cols_processed:(cols_processed+total_numb_input_files-1)]
		cols_processed = cols_processed + total_numb_input_files
		print(colnames(subset_df))
		rows_pass_consistency = c(rows_pass_consistency, checkConsistency_across_expts(subset_df, total_numb_input_files))
	}
	outdf = outdf[rows_pass_consistency,]
	# write the detailed output
	outputFile = paste(infile,"delta-fc-fixed.csv",sep='-')
    write.csv(outdf, outputFile, row.names=FALSE, quote=F)
	# write the names of the consistent genes
	outputFile_consis = paste(outputFile,"consis_genes.csv",sep='-')
    write.table(row.names(outdf), outputFile_consis, row.names=FALSE, quote=F, col.names=FALSE)
	#s_df = 
}


checkConsistency_across_expts = function(s_df, total_numb_input_files, foldchThreshold=1){

		res_pos = apply(s_df, 1, function(x) sum(x > foldchThreshold))
        res_neg = apply(s_df, 1, function(x) sum(x < foldchThreshold))
        rows_passing_consistency = c()
        #print("Up regulation")
        #print(head(res_pos))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_pos[[idx]])) && res_pos[[idx]] == (total_numb_input_files)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency$Upreg = rows_passing_consistency

        #rows_passing_consistency = c()
        #print("Down regulation")
        #print(head(res_neg))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_neg[[idx]])) && res_neg[[idx]] == (total_numb_input_files)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency$Dwnreg = rows_passing_consistency
		return(rows_passing_consistency)
}

if(length(args)<2){
	print("Script needs atleast 2 args: one file and one number indicating how many files(expt) were merged\n")
	print("Usage: Rscript ~/Morgun_Lab/richrr/scripts/R/process-delta-consistency.R Analys-pX-consis.csv-comb-pval-output.csv.combpval.Y.combfdr.Z.cutoff.csv ... 3\n")
	quit()
}

total_numb_input_files = as.integer(args[length(args)])
counter=1
while(counter <= length(args)-1){
	check_delta_consistency(args[counter],total_numb_input_files)
	counter = counter+1
}