library(argparser)


# usage: Rscript ~/Morgun_Lab/richrr/scripts/R/convert-groups-to-datasets.R --files /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-liver-ileum/rel_quantile/merged/corr/sp/p1_cp1_cfdr1/per_analysis/Analys1-consis.csv-comb-pval-output.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-liver-ileum/rel_quantile/merged/corr/sp/p1_cp1_cfdr1/per_analysis/Analys2-consis.csv-comb-pval-output.csv --output test-Analysis-1-2-consis-corr-dir --outstr Analys\ 12\ li_il_T2_HFHSNCD


#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('convert the merged files per group to merged file of datasets. merged file of 2 hfhs, 2 ncd becomes merged file of 4 dataset')
p <- add_argument(p, "--files", help="files where (corresponding (for --parallel)) columns are pasted one after other", nargs=Inf) # required; but written as optional format so I explicitly mention the "--files"
p <- add_argument(p, "--output", help="output file", default="merged_files_")
p <- add_argument(p, "--pergrdataset", help="datasets (expts) per group", default=2, type="numeric")
p <- add_argument(p, "--corrmethod", help="correlation method used", default="Spearman")
p <- add_argument(p, "--outstr", help="analysis number and name for new group", default="Analys 12 newgroup")



if(FALSE){
p <- add_argument(p, "--correlThreshold", help="correlThreshold", default=0, type="numeric")
p <- add_argument(p, "--pvalThreshold", help="pvalThreshold", default=0.03, type="numeric") # individual Pvalue Cutoff
p <- add_argument(p, "--foldchThreshold", help="foldchThreshold", default=0, type="numeric") # if the data is log transformed, we calc fold change by log(A/B) which is log(A) - log (B), so the difference is compared against 0. since log transform is default the arg defaults to 0. use 1 if using non log transformed data

# use these only in case of genes (not pairs)
p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.05, type="numeric") 
p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.1, type="numeric") 


p <- add_argument(p, "--warnings", help="print warnings", flag=TRUE)
p <- add_argument(p, "--foldchMean", help="use fold change from mean", flag=TRUE)  # default fold change is median
}
 
argv <- parse_args(p)
#print (p)



files = argv$files
pergrdataset = argv$pergrdataset
outputFile = argv$output
columns_to_keep = 1 + (2*pergrdataset) # id + (coeff, pvalue)*number of datasets per group
number_of_new_datasets = length(argv$files) * pergrdataset

corrmethod = argv$corrmethod
outstr = argv$outstr


#correlThreshold = argv$correlThreshold 
#pvalThreshold = argv$pvalThreshold
#foldchThreshold = argv$foldchThreshold

#combinedPvalueCutoff = argv$combPvalCutoff
#combinedFDRCutoff = argv$combFDRCutoff

# create if outdir doesn't exist:
#res_directory = paste(c("./p", pvalThreshold, "_cp", combinedPvalueCutoff, "_cfdr", combinedFDRCutoff, "/"), collapse='')
#dir.create(paste(res_directory, "per_analysis", sep=''), recursive = TRUE)
#outputFile = paste(res_directory, outputFile, sep='')


#foldchVar = 'FolChMedian'
#if(argv$foldchMean){foldchVar = "FoldChange"}

#outputFile = paste(outputFile , foldchVar, '_', sep='')






# first file taken as is for making a data frame
df = read.csv(files[1], header=T)
head(df)

# keep columns depending on number of expts per group
df = df[,c(1:columns_to_keep)]
head(df)

for (f in files[2:length(files)]){
    df_tmp = read.csv(f, header=T)
	df_tmp = df_tmp[,c(1:columns_to_keep)]
	df = merge(df, df_tmp, by=1)
}
# at this point df has common pairs but may not be consistent between hfhs and ncd.
# since these are already consistent in two expts, we can use only the first expt to
# check consistency between hfhs and ncd.

rownames(df) = df[,1]
head(df)



selcols = grep("Coefficient_Expt_1", colnames(df), value=T)
df_selcols = df[, selcols]
head(df_selcols)



# either 1, -1, 0 or a combination thereof in unique. length of 1 indicates consistency in groups
consistency_check = function(x){
    out = length(unique(sign(x)))
    out
}

len_uniq_sign = apply(df_selcols, 1, consistency_check)
df = cbind(df, len_uniq_sign)

head(df)

df_consis = df[df$len_uniq_sign == 1, ]
dim(df)
dim(df_consis)


write.table(df_consis, paste0(outputFile, "-info.txt"), quote=F, row.names=F, sep='\t')
write.csv(rownames(df_consis), paste0(outputFile, ".txt"), quote=F, row.names=F)


# drop the first and last column, replicate row names for first 4 columns
df_consis = df_consis[, -c(1, ncol(df_consis))]
df_new = cbind(replicate(number_of_new_datasets, rownames(df_consis)), df_consis)
head(df_new)


# edit the group and analysis
# names such that they seem as 4 expts of the same group. add the pairName cols to indicate 4 expts.

colnames(df_new) = gsub(paste0("Analys.(.+).", corrmethod), paste0(outstr, " ", corrmethod), colnames(df_new))

colnames(df_new) = gsub("Expt_\\d+", "Expt_", colnames(df_new))
df_new = df_new[ , order(names(df_new))]
colnames(df_new)[1:number_of_new_datasets] = paste0("pairName_Expt_" , colnames(df_new)[1:number_of_new_datasets])

write.csv(df_new, paste0("../merged-", outputFile, number_of_new_datasets, "-datasets.csv"), quote=F, row.names=F)

