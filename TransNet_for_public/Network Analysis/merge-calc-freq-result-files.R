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
p <- add_argument(p, "--output", help="output file", default="./merged_freq_files_")
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
p <- add_argument(p, "--foldchThreshold", help="foldchThreshold", default=0, type="numeric") # if the data is log transformed, we calc fold change by log(A/B) which is log(A) - log (B), so the difference is compared against 0. since log transform is default the arg defaults to 0. use 1 if using non log transformed data
p <- add_argument(p, "--foldchMean", help="use fold change from mean", flag=TRUE)  # default fold change is median

# use these only in case of genes (not pairs)
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.05, type="numeric") 
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.1, type="numeric") 


p <- add_argument(p, "--warnings", help="print warnings", flag=TRUE)
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
res_directory = paste(c("./p", pvalThreshold, "/freq/"), collapse='')
dir.create(res_directory, recursive = TRUE)
outputFile = paste(res_directory, outputFile, sep='')


foldchVar = 'Median'
if(argv$foldchMean){foldchVar = " Mean"} # to separate from GeoMean
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
# this finds the analysis number 
#------------------------------------------------------------------
find_analysis_number = function(x){
   #res = as.numeric(str_match(x, "Analys ([0-9]+) .*")[,2])
   res = str_match(x, "Analys ([0-9]+-?[0-9]*) .*")[,2]
   res
}

#------------------------------------------------------------------
# this finds the analysis name
#------------------------------------------------------------------
find_analysis_name = function(x){
   res = str_match(x, "Analys [0-9]+-?[0-9]* (.*)_Expt_.*")[,2]
   res
}


#print(head(merged_df))
#------------------------------------------------------------------
# calc median freq across expts.
# Note: at this point we do not care about consistency in fold change and/or p-value
#------------------------------------------------------------------
ID=rownames(merged_df)
result_dumper = merged_df[,1,drop=FALSE]  
result_dumper = cbind(result_dumper,ID)
result_dumper = result_dumper[,-1,drop=FALSE]
#print(head(result_dumper))
cols_processed = 1
while(cols_processed < ncol(merged_df))  # loop over the merged dataset in steps of (number of input files)
{
   subset_df = merged_df[,cols_processed:(cols_processed+total_numb_input_files-1)]
   cols_processed = cols_processed + total_numb_input_files
   result_dumper = cbind(result_dumper, subset_df)
   
   	if(grepl(foldchVar, colnames(subset_df)[1])){

	SelecColnames = colnames(subset_df)[grep(foldchVar,colnames(subset_df))]
	#print(SelecColnames)

	# this returns the analysis number and name
    res1=lapply(SelecColnames, find_analysis_number)
    analys_numb = unique(unlist(res1))[1] 
	#print(analys_numb)
    res2=lapply(SelecColnames, find_analysis_name)
    analys_name = unique(unlist(res2))[1] 
	#print(analys_name)

	interestedSelecData = subset_df[,SelecColnames]
	interestedSelecData = as.matrix(interestedSelecData)
	interestedSelecData = apply(interestedSelecData,2,function(x){as.numeric(as.vector(x))})
	combinedSelec = apply(interestedSelecData,1, function(x){round(median(x, na.rm = TRUE), 3)})
	result_dumper = cbind(result_dumper,combinedSelec)
	colnames(result_dumper)[length(result_dumper)] = paste(analys_numb, analys_name, sep='_')
	
   }
   
}

outputFile2 = paste(outputFile,"med.csv",sep='-')

write.csv(result_dumper, outputFile2, row.names=FALSE)

if(argv$warnings){print(warnings())}

print("Finished performing the requested analyses.")

print("You may want to run create-network.R next.")


