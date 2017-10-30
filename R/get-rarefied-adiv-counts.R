# based on the documentation of http://qiime.org/scripts/compare_alpha_diversity.html , qiime averages the values of multiple iterations at the maximum depth of the collated file and uses that as a single number per sample. to give the reviewer rarefied counts, I will do the same and create a table similar to the counts obtained from unrarefied data.


library(stringr)

# cd /nfs3/PHARM/Morgun_Lab/richrr/Lean_Abx_data/analysis/microbes/post-biom-analysis/otu_tables/99perc_cum_abund/quantile/core-div-analysis

# Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/get-rarefied-adiv-counts.R /nfs3/PHARM/Morgun_Lab/richrr/Lean_Abx_data/analysis/microbes/post-biom-analysis/otu_tables/99perc_cum_abund/quantile/core-div-analysis/corediv_c1_e200000/adiv_mc200000/arare_max200000/alpha_div_collated c1_rarefied_200K_avg.csv

# Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/get-rarefied-adiv-counts.R /nfs3/PHARM/Morgun_Lab/richrr/Lean_Abx_data/analysis/microbes/post-biom-analysis/otu_tables/99perc_cum_abund/quantile/core-div-analysis/corediv_c2_e200000/adiv_mc200000/arare_max200000/alpha_div_collated c2_rarefied_200K_avg.csv

args = commandArgs(trailingOnly=TRUE)
base = args[1]
outfile = args[2]
files = list.files(base, full.names=TRUE)

big_res = ''

for(f in files){
	
	print(f)
	df = read.csv(f ,header = TRUE,check.names=FALSE,na.strings=c("","na","NA", "Na", "NaN"), row.names=1, sep="\t")
	#print(head(df))
	df_max = tail(df, 10) # get the iterations at the maximum depth. I know I had 10 iterations.
	#print(df_max)
	
	fname = gsub(base, "", f)
	fname = gsub("/", "", fname)
	fname = gsub(".txt", "", fname)
	
	out = t(colMeans(df_max))
	rownames(out) = c(fname)
	print(out)
	
	big_res = rbind(big_res, out)
	
}

# delete the first dummy empty row
big_res = big_res[-1, ]
head(big_res)

write.csv(t(big_res), outfile, row.names=TRUE, quote=F)
