library(ggplot2)


get_base_stats_exp_data = function (infile, backgrond, outputFile, sample_perc_threshold=10){
 
	#infile = "~/Morgun_Lab/richrr/Cervical_Cancer/Data/gene-expression/HumanWG6-IlluminaGeneExpr_QuantileNormalized_Cohort1.csv"
	#backgrond = 65 #100 #40

	#infile = "~/Morgun_Lab/richrr/Cervical_Cancer/Data/gene-expression/HumanHT12v4_IlluminaGeneExpr_QuantileNorm_Cohort2.csv"
	#backgrond = 110 #150 #70

	df = read.csv(infile,row.names = 1, header=T, check.names=F,na.strings=c("","na","NA", "Na", "NaN"))

	#head(df)

	# overall stats
	mindf = min(df) 
	maxdf = max(df)

	ldf = sort(unique(round(as.vector(t(df)))))
	print(ldf[1:10])

	#print(mindf)
	#print(maxdf)
	tot_num_samples = ncol(df)

	# per row stats
	meandf = apply(df, 1, mean)
	meandf = data.frame(meandf)
	#head(meandf)

	mediandf = apply(df, 1, median)
	mediandf = data.frame(mediandf)
	#head(mediandf)

	sddf = apply(df, 1, sd)
	sddf = data.frame(sddf)
	#head(sddf)

	covdf = sddf/meandf
	colnames(covdf) = c("covdf")
	#head(covdf)

	# number of samples with value greater than background
	tot_samps_more_than_bkg = apply(df, 1, function(x) sum(abs(x) > backgrond))
	tot_samps_more_than_bkg = data.frame(tot_samps_more_than_bkg)
	#head(tot_samps_more_than_bkg)

	# percent of samples with value greater than background
	perc_samps_more_than_bkg = (100*tot_samps_more_than_bkg) / tot_num_samples
	colnames(perc_samps_more_than_bkg) = c("%_samps_more_than_bkg")
	#head(perc_samps_more_than_bkg)
	probes_to_remove = sum(perc_samps_more_than_bkg < sample_perc_threshold)  # in librecalc =COUNTIF(F2:F48804, "<=10")
	print(paste(c("Numb. of probes with <",sample_perc_threshold,"% samples with value greater than raw file's background", " (", backgrond, "): ", probes_to_remove), collapse=""))


	resultdf =  cbind(meandf, mediandf, covdf, tot_samps_more_than_bkg, perc_samps_more_than_bkg)
	#head(resultdf)
	write.csv(resultdf, paste(outputFile, "base_stats_output.csv" ,sep = "_"))
	
	probes_to_remove_names = resultdf[resultdf[,"%_samps_more_than_bkg"] < sample_perc_threshold, ]
	write.csv(probes_to_remove_names, paste(outputFile, "remove_probes.csv" ,sep = "_"))

}



