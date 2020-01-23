library(argparser)

# uage: cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/summarize_per_sample/relativised/brb-16s-improv-expts/9-11-2018/16s-improv-expts-rnaseq-network
# Rscript ~/Morgun_Lab/richrr/scripts/R/split-analys-file-as-separate-datasets.R --file RQ/poolexpt/ile/sp_corr-output.csv --outstr poolexpt_meta_microbe_ile
#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('convert the analysis file (obtained from perform_analysis) containing several analysis to several datasets of one analysis. analysis file of 1 hfhs, 1 ncd becomes 2 separate datasets')
p <- add_argument(p, "--file", help="files to be split per analysis into datasets", nargs=1) # required; but written as optional format so I explicitly mention the "--files"
p <- add_argument(p, "--output", help="output file", default="_dataset_")
#p <- add_argument(p, "--pergrdataset", help="datasets (expts) per group", default=2, type="numeric")
#p <- add_argument(p, "--corrmethod", help="correlation method used", default="Spearman")
p <- add_argument(p, "--outstr", help="analysis number and name for new group", default="newgroup")


argv <- parse_args(p)
files = argv$file
#print(files)

#pergrdataset = argv$pergrdataset
outputFile = argv$output
# = 1 + (2*pergrdataset) # id + (coeff, pvalue)*number of datasets per group
#number_of_new_datasets = length(argv$files) * pergrdataset

#corrmethod = argv$corrmethod
outstr = argv$outstr


df = read.csv(files[1], header=T, check.names=F, row.names=1)
#head(df)
#print(ncol(df))

number_of_new_datasets = length(grep('Coefficient', colnames(df), value=T))  # corr file
columns_to_keep = 3
Id = 'pairName'

if(number_of_new_datasets == 0){ # comp file
	number_of_new_datasets = length(grep('FolChMedian', colnames(df), value=T))
	columns_to_keep = 11
	Id = 'geneName'
}

# keep columns depending on number of expts per group
endcol = 0
for(col in c(1:number_of_new_datasets)){
	print(col)
	df_tmp = df[,c((endcol+1):(endcol + columns_to_keep))]
	endcol = endcol + columns_to_keep
	#print(endcol)
	colnames(df_tmp) = gsub("(Analys .+ )", paste0("Analys 1 ", outstr, " "), colnames(df_tmp))
	IdName = rownames(df_tmp)
	
	df_tmp = cbind(IdName, df_tmp)
	colnames(df_tmp)[1] = Id
	
	print(head(df_tmp))
	write.csv(df_tmp, paste0(c(files,outputFile,col,".csv"), collapse=''),quote=F, row.names=F)
	
}


q()