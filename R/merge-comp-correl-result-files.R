## to do:
# implement code for majority of datasets as consistency

# Usage: Rscript /nfs1/Morgun_Lab/richrr/scripts/R/merge-comp-correl-result-files.R --files network_spearm-manwhit_expt1_comp-output.csv network_spearm-manwhit_expt2_comp-output.csv --parallel --output merged_two_expts_sp_mw_files_comp_all_



library(argparser)
#library(hash)
#library(psych)
#library(gtools)
#library(corpcor)
#library(reshape)
library(stringr)



#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Calculate which gene or pairs show consistency across experiments for same (comparison and/or correlation) analysis')
p <- add_argument(p, "--files", help="files where (corresponding (for --parallel)) columns are pasted one after other", nargs=Inf) # required; but written as optional format so I explicitly mention the "--files"
p <- add_argument(p, "--output", help="output file", default="./merged_files_")
p <- add_argument(p, "--parallel", help="merge file1-col1 file2-col1 ... filen-col1 file1-col2 file2-col2 ... filen-col2", flag=TRUE) # only option; # else prints all cols from file 1 followed by all from file 2                                                                                   
p <- add_argument(p, "--pvals", help="merge only columns with pvalues. Default ALL columns", flag=TRUE) # default all columns
p <- add_argument(p, "--nofileinfo", help="do not keep track of which column came from which file", flag=TRUE)  # do not use this flag. It can get confusing to know which sample came from which file if this flag is used.
p <- add_argument(p, "--majority", help="the element is consistent across MAJORITY of the datasets. Default is consistent across ALL datasets.", flag=TRUE)
     # Note that to use this option you need to implement the consistency check differently. Just switching the condition FROM (= numb. files) TO (>= numb. files/2) is not sufficient
     # Remember that a gene may be up (sig.), up (not-sig.), down (sig.) and still show up as consistent if you check for pvalue and FC/co-eff separately
     # You need to implement code so that pvalue AND FC/co-eff is simultanesouly checked if you plan to use the majority option 
     # this is not a problem if you use the default of consistent across ALL datasets
# thresholds to check consistency
p <- add_argument(p, "--correlThreshold", help="correlThreshold", default=0, type="numeric")
p <- add_argument(p, "--pvalThreshold", help="pvalThreshold", default=0.03, type="numeric") # individual Pvalue Cutoff
p <- add_argument(p, "--foldchThreshold", help="foldchThreshold", default=1, type="numeric")

# use these only in case of genes (not pairs)
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.05, type="numeric") 
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.1, type="numeric") 


p <- add_argument(p, "--warnings", help="print warnings", flag=TRUE)
p <- add_argument(p, "--foldchMean", help="use fold change from mean", flag=TRUE)  # default fold change is median
p <- add_argument(p, "--transposeOutput", help="tranpose the results so the rows are analysis and columns are gene or pairs", flag=TRUE)  # easier for downstream grep of required analysis, do not use when you expect many genes or pairs since it might truncate when you open in excel or librecalc due ot limited columns


 
argv <- parse_args(p)
#print (p)


if(length(argv$files) < 2)
{
  print("At least 2 files are required. Give --files file1 file2 ... in cmd")
  quit()
}

if(argv$majority){
  print("This option has not yet been implemented. Bye!")
  quit()
}

outputFile = argv$output
correlThreshold = argv$correlThreshold 
pvalThreshold = argv$pvalThreshold
foldchThreshold = argv$foldchThreshold

combinedPvalueCutoff = argv$combPvalCutoff
combinedFDRCutoff = argv$combFDRCutoff

# create if outdir doesn't exist:
res_directory = paste(c("./p", pvalThreshold, "/"), collapse='')
dir.create(paste(res_directory, "per_analysis", sep=''), recursive = TRUE)
outputFile = paste(res_directory, outputFile, sep='')


foldchVar = 'FolChMedian'
if(argv$foldchMean){foldchVar = "FoldChange"}

outputFile = paste(outputFile , foldchVar, '_', sep='')

# number of columns
max_numb_cols = 0

df = list()
#------------------------------------------------------------------
# read the files one after the other:
#------------------------------------------------------------------
for( numb in 1:length(argv$files)){
   df[numb] = list(read.csv( argv$files[numb], header=TRUE,check.names=FALSE))

   # set the max_numb_cols
   if(length(df[[numb]]) > 0 && max_numb_cols == 0){
       max_numb_cols = length(df[[numb]])
   } else if(max_numb_cols != length(df[[numb]])){ # this has been set before
       print("The number of columns in the input files are different")
       quit()
   } 
}


#------------------------------------
# keep only the common rows since only they can be consistent
#------------------------------------
all_present_rows = ''
for(file in 1:length(df)){
    rown = df[[file]][1]
    #print(typeof(rown))
    #print(class(rown))

    if(file == 1)
      {
         all_present_rows = rown
      }
    else
    {
        all_present_rows = merge(all_present_rows, rown)
    }
    
}
print(paste(c("Using", nrow(all_present_rows), "common", "rows"), collapse=' '))
#print(head(all_present_rows))

# keep only the common rows present in all files
intersect_rows_df = list()
for(file in 1:length(df)){
    mergd_df = merge(all_present_rows, df[[file]])
    intersect_rows_df[[file]] = mergd_df
}


df = intersect_rows_df

# add this arg in name of output file
if(argv$pvals){
outputFile = paste(outputFile, "only-pvals-", sep = '_')
}


#------------------------------------------------------------------
# merge the files:
#------------------------------------------------------------------
if(!argv$parallel){
    outputFile = paste(outputFile,"output.csv",sep='merged-serial-')
    write.csv(as.data.frame(df), outputFile, row.names=FALSE)
} else {

    para_df = c()

    for(col in 1:max_numb_cols){
        
        for(file in 1:length(df)){
            
            if(file == 1 && col == 1){
                para_df = df[[file]][col]
            } else {
                # if only pvalue columns requested
                if(argv$pvals){
                    if(grepl("pvalue" , colnames(df[[file]][col]))){
                        para_df = cbind(para_df, df[[file]][col])
                    }
                    else{
                        next # skip to next iteration, don't change the column label
                    }
                } else {
                # add all columns
                    para_df = cbind(para_df, df[[file]][col])
                }
            }
            if(!argv$nofileinfo){  # keep track of which file the column came from
                new_col_name = paste(c(toString(colnames(df[[file]][col])), "Expt" , file), collapse="_")
                colnames(para_df)[ncol(para_df)] = new_col_name            
            }
        }
    }
   outputFile = paste(outputFile,"output.csv",sep='merged-parallel-')
   write.csv(as.data.frame(para_df), outputFile, row.names=FALSE)

   if(argv$transposeOutput){
       #outputFile_tr = paste(outputFile,"transposed.csv",sep='-')
       #write.csv(t(as.data.frame(para_df)), outputFile_tr, col.names=FALSE)
   }


}



if(argv$pvals){
    print("Cannot check for change in direction without FoldChange or Coefficient columns. Do NOT use '--pvals'. I quit!")
    quit()
}



#------------------------------------------------------------------
# calculates the combined fisher p value for each gene for each comparison
# the last arg lets you control whether to apply the combined pvalue and fdr cutoff to the results
#------------------------------------------------------------------
calc_comb_pval_fdr = function(in_df, total_numb_input_files, l_outputFile, applyCombCutoff = FALSE){

  cols_processed = 1
  comb_in_df = data.frame(matrix(NA, nrow = nrow(in_df), ncol = 2)) # ceate dataframe of "NA" columns
  rownames(comb_in_df) = rownames(in_df)
  
  while(cols_processed < ncol(in_df))  # loop over the merged dataset in steps of (number of input files)
  {
      subset_df = in_df[,cols_processed:(cols_processed+total_numb_input_files-1)]
      cols_processed = cols_processed + total_numb_input_files
      comb_in_df = cbind(comb_in_df, subset_df)
      
      
      if(grepl("pvalue", colnames(subset_df)[1])){

          outForPUC = subset_df
          # calculate combined Pvalue for interest group
          PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
          total_numb_input_files = length(PvalueColnames)
          interestedPvalueData = outForPUC[,PvalueColnames]
          #print(interestedPvalueData)
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
	   #calculate FDR for combined pvalue
           combinedFDR = p.adjust(combinedPvalue,method="fdr")
  

           comb_in_df = cbind(comb_in_df,combinedPvalue,combinedFDR)

      } # end if for p value
   }
   comb_in_df = comb_in_df[,-c(1,2)] # remove the tmp "NA" column added during initialization
   outputFile1 = paste(l_outputFile,"output.csv",sep='-comb-pval-')
   write.csv(as.data.frame(comb_in_df), outputFile1)
   
   if(argv$transposeOutput){
       outputFile_tr = paste(outputFile1,"transposed.csv",sep='-')
       write.csv(t(as.data.frame(comb_in_df)), outputFile_tr) #, col.names=FALSE) # since the col.names is ignored anyways but unnecessarily throws a warning about it
   }

   # apply the cutoff
   if(applyCombCutoff){
       outputFile_tr = paste(c(outputFile1,"combpval",combinedPvalueCutoff, "combfdr",combinedFDRCutoff,"cutoff.csv"),collapse='.')
       #print(head(as.data.frame(comb_in_df)))
       #quit()
       new.data <- comb_in_df[comb_in_df$"combinedPvalue" < combinedPvalueCutoff  &  comb_in_df$"combinedFDR" < combinedFDRCutoff, ]
       #new.data <- subset(as.data.frame(comb_in_df), "combinedPvalue" < combinedPvalueCutoff  &  "combinedFDR" < combinedFDRCutoff)
       write.csv(new.data, outputFile_tr)
       
       # this prints only the row names incase you want to use this for calc correl of only these genes
       out_vec = rownames(new.data)
       write(out_vec, paste(outputFile_tr, "consis_genes.csv", sep='.'))
       
   }


}



#------------------------------------------------------------------
# re-read the merged file as a dataframe 
#------------------------------------------------------------------
total_numb_input_files = length(argv$files)
merged_df = read.csv( outputFile, header=TRUE,check.names=FALSE)
#print(head(merged_df)[1:5])
if(argv$pvals){
    rownames(merged_df) <- merged_df[, 1] ## set rownames
    merged_df <- merged_df[, -1]          ## remove the first column
} else {
  # check whether the first "number of files" columns are the same and remove them after setting row names
         l_identical_inp = as.character(merged_df[, 1])
                
         if(total_numb_input_files>1){
             for(z in 2:total_numb_input_files){
                    l_identical_inp = cbind(l_identical_inp, as.character(merged_df[,z]))
             }

             if(all(apply(l_identical_inp, 2, identical, l_identical_inp[, 1]))){
                    print("All are equal")
                    rownames(merged_df) <- merged_df[, 1]                           ## set rownames
                    merged_df <- merged_df[, -c(1:total_numb_input_files)]          ## remove the first "total numb of input files" columns      
             } else {
                    print(l_identical_inp)
                    print("All are NOT equal")
                    quit()
             }
          }
}


#------------------------------------------------------------------
# calculates the combined fisher p value for each gene for each comparison
# it does not care about the gene being consistent in either FC dir or pval significance. 
# calcs. FDR. Of course once you keep consistent genes, you may have to recalc fdr but NOT combined pvalue
#------------------------------------------------------------------
calc_comb_pval_fdr(merged_df, total_numb_input_files, outputFile)



#------------------------------------------------------------------
# check which measurement ("gene") or pairs have the same trend across experiments
# identify the rows where values are >/< than threshold
#------------------------------------------------------------------
checkConsistency_across_expts = function(s_df, condition, total_numb_input_files, correlThreshold , pvalThreshold, foldchThreshold){

    
    list_rows_passing_consistency = list()
    
    if(condition == "Coefficient"){
        res_pos = apply(s_df, 1, function(x) sum(x > correlThreshold))
        res_neg = apply(s_df, 1, function(x) sum(x < correlThreshold))
        rows_passing_consistency = c()
        #print("Positive correlation")
        #print(head(res_pos))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_pos[[idx]])) && res_pos[[idx]] == (total_numb_input_files)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[1]] = rows_passing_consistency
        list_rows_passing_consistency$PosCorrel = rows_passing_consistency

        rows_passing_consistency = c()
        #print("Negative correlation")
        #print(head(res_neg))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_neg[[idx]])) && res_neg[[idx]] == (total_numb_input_files)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[2]] = rows_passing_consistency
        list_rows_passing_consistency$NegCorrel = rows_passing_consistency

    } else if(condition == "pvalue"){
        res_pos = apply(s_df, 1, function(x) sum(x < pvalThreshold))
        rows_passing_consistency = c()
        #print("Pass pvalue cutoff")
        #print(head(res_pos))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_pos[[idx]])) && res_pos[[idx]] == (total_numb_input_files)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[1]] = rows_passing_consistency
        list_rows_passing_consistency$PvalThresh = rows_passing_consistency

    } else if(condition == foldchVar){
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
        #list_rows_passing_consistency[[1]] = rows_passing_consistency
        list_rows_passing_consistency$Upreg = rows_passing_consistency

        rows_passing_consistency = c()
        #print("Down regulation")
        #print(head(res_neg))
        for(idx in 1:nrow(s_df)){
            if((!is.na(res_neg[[idx]])) && res_neg[[idx]] == (total_numb_input_files)){
                rows_passing_consistency = append(rows_passing_consistency,rownames(s_df)[idx])
            }
        }
        #list_rows_passing_consistency[[2]] = rows_passing_consistency
        list_rows_passing_consistency$Dwnreg = rows_passing_consistency

    } 
    #print(list_rows_passing_consistency)
    return(list_rows_passing_consistency)
}


#------------------------------------------------------------------
# check Consistency across ALL expts using the thresholds 
#------------------------------------------------------------------
result_dumper = list()  # to dump the consistent results in output files
cols_processed = 1
while(cols_processed < ncol(merged_df))  # loop over the merged dataset in steps of (number of input files)
{
   subset_df = merged_df[,cols_processed:(cols_processed+total_numb_input_files-1)]
   cols_processed = cols_processed + total_numb_input_files

   
   if(grepl("Coefficient", colnames(subset_df)[1])){
       out_res = list()
       out_res$name = colnames(subset_df) #paste(colnames(subset_df), collapse='<-->')
       out_res = append(out_res, checkConsistency_across_expts(subset_df, "Coefficient", total_numb_input_files, correlThreshold , pvalThreshold, foldchThreshold))
       #result_dumper = append(result_dumper, in_cols)
       #result_dumper = append(result_dumper, out_res)
       result_dumper[[length(result_dumper)+1]] <- out_res
       
   } else if(grepl("pvalue", colnames(subset_df)[1])){
       out_res = list()
       out_res$name = colnames(subset_df) #paste(colnames(subset_df), collapse='<-->')
       out_res = append(out_res, checkConsistency_across_expts(subset_df, "pvalue", total_numb_input_files, correlThreshold , pvalThreshold, foldchThreshold))
       #result_dumper = append(result_dumper, in_cols)
       #result_dumper = append(result_dumper, out_res)
       result_dumper[[length(result_dumper)+1]] <- out_res
       
   } else if(grepl(foldchVar, colnames(subset_df)[1])){
       out_res = list()
       out_res$name = colnames(subset_df) #paste(colnames(subset_df), collapse='<-->')
       out_res = append(out_res, checkConsistency_across_expts(subset_df, foldchVar, total_numb_input_files, correlThreshold , pvalThreshold, foldchThreshold))
       #result_dumper = append(result_dumper, in_cols)
       #result_dumper = append(result_dumper, out_res)
       result_dumper[[length(result_dumper)+1]] <- out_res
       
   }
   
}

         
#------------------------------------------------------------------
# Print the consistent results for each aspect (e.g. pos correl, neg correl, pvalue, FC, etc.) to file 
#print(result_dumper)
#dput(result_dumper, file = paste(outputFile,"dump_consistent_output.csv",sep='-'))
#print(length(result_dumper))
#------------------------------------------------------------------
outputFile1 = paste(outputFile,"consistent_across_expts_output.csv",sep='-')
strr = '' 
for(i in 1:length(result_dumper)){
    p_elem = result_dumper[[i]]
    l_strr = ''
    for(x in 1:length(p_elem)){ # number of child elements (e.g. names, pval, pcorr, ncorr, etc.)
        length_c_elem = toString(lengths(p_elem)[x]) # gives the number of elements in the child elements 
        name_c_elem = names(lengths(p_elem)[x]) # name of the child element
        elem = p_elem[x]
        c_elem = elem[[1]] # this returns a vector of all the elements in the child element
        c_elemts_str = paste(c_elem, collapse=',') # convert the elements to string
        c_elemts_str = paste(name_c_elem, c_elemts_str , sep=':\t') # add the name of the child element
        l_strr = paste(l_strr, c_elemts_str , sep='\n')
    }
    strr = paste(strr, l_strr, sep='\n')   
}
write(strr, file=outputFile1)


#------------------------------------------------------------------
# this finds the analysis number 
#------------------------------------------------------------------
find_analysis_number = function(x){
   #res = as.numeric(str_match(x, "Analys ([0-9]+) .*")[,2])
   res = str_match(x, "Analys ([0-9]+-?[0-9]*) .*")[,2]
   res
}


#------------------------------------------------------------------
# finds which gene (or pairs) show consistency (across expts) in multiple aspects of the same analysis (e.g. direction of comparison (or correlation) and sign of pval)
#------------------------------------------------------------------
extract_common_vals = function(nam, searchIn, searchThese){
   res = intersect(searchIn[[nam]],searchThese)
   res
}


#------------------------------------------------------------------
# finds which gene pairs show consistency (across expts) in multiple aspects of the same analysis (e.g. direction of correlations and sign of pval) 
# (e.g. pos correl of first, neg correl of second and pval  i.e. pos correl of first but not second and pval
#       neg correl of first, pos correl of second and pval  i.e. neg correl of first but not second and pval)
#------------------------------------------------------------------
extract_common_vals_3_vectors = function(names, searchIn1, searchIn2, searchThese){
 res = list()
 for(nam in names){
   resint = intersect(searchIn1[[nam]],searchThese)
   res = append(res, setdiff(resint, searchIn1[[nam]])) #produces all elements of the first input vector without any matching elements from the second input vector
 }
 res
}


#------------------------------------------------------------------
# Calc comb p val and fdr for consistent genes/pairs (based on FC/Coeff and pvalue) per analysis 
#------------------------------------------------------------------
CalcCombPvalFdrPerAnalysis = function(l_strr, consis_elem, merged_df){

  analyses = paste(c("Analys ", l_strr, " "), collapse='')
  
  # paste(c(res_directory, "per_analysis/","Analys ", l_strr, " "), collapse='')
  #print(analyses)
  #print(length(consis_elem))
  
  #print(rownames(merged_df))
  
  select_analys_cols = grep(analyses, colnames(merged_df), value=T)
  
  tmp_out_df = merged_df[consis_elem, select_analys_cols]
  #print(nrow(tmp_out_df))
  
  o_f_name_cmob_pval_fdr = gsub(' ' , '', paste(res_directory, "per_analysis/",analyses, paste(c("-p", pvalThreshold), collapse=''), "-consis.csv", sep=''))
  
  #print(head(tmp_out_df))
  if(nrow(tmp_out_df) > 0){ # if no genes have reached uptil this point, no point to call the function
      calc_comb_pval_fdr(tmp_out_df, total_numb_input_files, o_f_name_cmob_pval_fdr, TRUE)
  }

}



#------------------------------------------------------------------
# Identify which "genes" or "pairs" are present in FC/Coeff and pvalue
# the last method arg runs the combined pval and fdr output per analysis code for only the genes and not for correlations
#------------------------------------------------------------------
SameDirectionAndPvalue = function(analy_name_c_elem1, analy_name_c_elem2, p_elem1, p_elem2, merged_df, isFCanalys){

    # this returns the analysis number    
    res1=lapply(analy_name_c_elem1, find_analysis_number)
    res2=lapply(analy_name_c_elem2, find_analysis_number)
    
    # if the columns are from different analysis
    if(length(union(res1, res2)) != 1){
        print(analys_processed)
        print(p_elem1)
        print(p_elem2)
        print(union(res1, res2))
        print("A: something ridiculous happened here, so I quit")
        quit() 
    } 
    
    
    names_elem1 = names(lengths(p_elem1)[-1])
    #print(p_elem1[[names_elem1[1]]])   
    names_elem2 = names(lengths(p_elem2)[-1])
    #print(names_elem2)
    
    res=lapply(names_elem1, extract_common_vals, p_elem1, p_elem2[[names_elem2[1]]])
    analysis_names = paste(paste(analy_name_c_elem1, collapse=',') , paste(analy_name_c_elem2, collapse=',') , sep='\nAND\n')
    
    l_strr = paste('\n', analysis_names, sep='')
    
    #print("THE LENGTH IS")
    #print(length(res))
    
    #print(res)
    
    # this happens when the fold change is not consistent but p value is
    if(length(res) == 0){ res = list(NULL) }

    all_consis_elements = c()    # this compiles the conistent upreg and dwnreg (or pos and neg pairs) before sending them to calc comb fdr
    for(x in 1:length(res)){ # number of child elements (e.g. names, pval, pcorr, ncorr, upreg, dwnreg, etc.)
        consis_elem = res[[x]]
        if(length(consis_elem) > 0){
            #print(length(consis_elem))
            all_consis_elements = c(all_consis_elements, consis_elem)
            #print(length(all_consis_elements))
        }
    

        if(length(consis_elem) > 0){
          l_strr = paste(l_strr, paste(consis_elem, collapse=',') , sep='\n')
        }
    }
    
    if(isFCanalys == 'FC'){
    # this part is useful to calc. comb FDR on only the consistent genes
        CalcCombPvalFdrPerAnalysis(union(res1, res2), all_consis_elements, merged_df)
    }

    return(l_strr)

}


#------------------------------------------------------------------

#------------------------------------------------------------------
CompareCorrelations = function(analy_name_c_elem1, analy_name_c_elem2, analy_name_c_elem3, p_elem1, p_elem2, p_elem3){

    # this returns the analysis number    
    res1=lapply(analy_name_c_elem1, find_analysis_number)
    res2=lapply(analy_name_c_elem2, find_analysis_number)
    res3=lapply(analy_name_c_elem3, find_analysis_number)

    # if the columns are from different analysis
    if(length(union(union(res1, res2), res3)) != 1){
        print(analys_processed)
        print(p_elem1)
        print(p_elem2)
        print(p_elem3)
        print(union(union(res1, res2), res3))
        print("B: something ridiculous happened here, so I quit")
        quit() 
    }
    
    # names of child element except the child element "name"
    names_elem1 = names(lengths(p_elem1)[-1])
    names_elem2 = names(lengths(p_elem2)[-1])
    names_elem3 = names(lengths(p_elem3)[-1])
    
    res=extract_common_vals_3_vectors( names_elem1, p_elem1, p_elem2, p_elem3[[names_elem3[1]]])
        
    analysis_names = paste(paste(analy_name_c_elem1, collapse=',') , paste(analy_name_c_elem2, collapse=',') , paste(analy_name_c_elem3, collapse=',') , sep='\nAND\n')
       
    l_strr = paste('\n', analysis_names, sep='')
    if(length(res) > 0) {
      for(x in 1:length(res)){ # number of child elements (e.g. names, pval, pcorr, ncorr, etc.)
        consis_elem = res[[x]]
        if(length(consis_elem) > 0){
          l_strr = paste(l_strr, paste(consis_elem, collapse=',') , sep='\n')
        }
      }
    }
    return(l_strr)

}




#------------------------------------------------------------------
# Check which measurement ("gene") (or pairs) show same trends in multiple aspects of analysis across expts: same direction of comparison (or correlation) and sig pval across expts
#------------------------------------------------------------------
outputFile2 = paste(outputFile,"consistent_across_expts_and_analysis_output.csv",sep='-')
analys_processed = 1
consis_strr = '' 

while(analys_processed < length(result_dumper)){
    p_elem1 = result_dumper[[analys_processed]]
    p_elem2 = result_dumper[[analys_processed+1]]
            
    elem1 = p_elem1[1] # child element "name"
    analy_name_c_elem1 = elem1[[1]] # this returns a vector of all the elements in the child element "name"
    elem2 = p_elem2[1] # child element "name"
    analy_name_c_elem2 = elem2[[1]] # this returns a vector of all the elements in the child element "name"

    # contains foldchange
    res_foldch = lapply(analy_name_c_elem1, function(x) grepl(foldchVar, x))
    # contains coefficient
    res_coeff = lapply(analy_name_c_elem1, function(x) grepl("Coefficient", x))
    # contains pvalue
    res_pvalu = lapply(analy_name_c_elem2, function(x) grepl("pvalue", x))

    check_pvalue_in_elem1 = lapply(analy_name_c_elem1, function(x) grepl("pvalue", x))
    
    #print(analys_processed)
    #print(unique(check_pvalue_in_elem1))
    if(unique(check_pvalue_in_elem1)[1] == "TRUE") { 
        # pvalue only; because fold change was not calculated for the delta(a-b) vs delta(c-d)
            analys_processed = analys_processed + 1 # do nothing and go to next
    } else if(identical(unique(res_coeff) , unique(res_pvalu))){
        # correlation coeff and pvalue
        l_strr = SameDirectionAndPvalue(analy_name_c_elem1, analy_name_c_elem2, p_elem1, p_elem2, merged_df, "")
        consis_strr = paste(consis_strr, l_strr, sep='\n')
        analys_processed = analys_processed + 2 # coefficient and pvalue   
    } else if(identical(unique(res_foldch) , unique(res_pvalu)) && unique(res_foldch)[1] == "TRUE" ){  
        # because in the compare correl, the first condition is satisfied (FALSE=FALSE) but not second 
        # foldch and pvalue
        l_strr = SameDirectionAndPvalue(analy_name_c_elem1, analy_name_c_elem2, p_elem1, p_elem2, merged_df, "FC")
        consis_strr = paste(consis_strr, l_strr, sep='\n')
        analys_processed = analys_processed + 2 # foldch and pvalue   
    } else { 
        # coefficient, coefficient and pvalue
        p_elem3 = result_dumper[[analys_processed+2]]
        elem3 = p_elem3[1] # child element "name"
        analy_name_c_elem3 = elem3[[1]] # this returns a vector of all the elements in the child element "name"
        res_pvalu = lapply(analy_name_c_elem3, function(x) grepl("pvalue", x))
        
        if(  identical(unique(res_coeff) , unique(res_pvalu))  ){
            l_strr = CompareCorrelations(analy_name_c_elem1, analy_name_c_elem2, analy_name_c_elem3, p_elem1, p_elem2, p_elem3)
            consis_strr = paste(consis_strr, l_strr, sep='\n')
            analys_processed = analys_processed + 3 # coefficient, coefficient and pvalue
        } else { 
            print("C: something ridiculous happened here, so I quit"); quit() 
        }  
    } 
}


print(consis_strr)
write(consis_strr, file=outputFile2)

finaloutputFile = paste(c(outputFile,  "pval", pvalThreshold ,"corr", correlThreshold, "foldch", foldchThreshold, "consistent_results.csv"), collapse='_')
write(gsub(",", "\n", consis_strr), file=finaloutputFile)


if(argv$warnings){print(warnings())}

print("Finished performing the requested analyses.")

print("You may want to run create-network.R next.")













