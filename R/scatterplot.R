library(stringr)


# cd ~/Morgun_Lab/richrr/Type2_Diabetes/microbe-pheno-analyses/analysis/m-p-pairs/rel.quantif_quantile/prelim_result/merged/comp/mw_sp/p1/
# Rscript scatterplot.R 3 merged_mp_mw_sp_comp_p1_FolChMedian_merged-parallel-output.csv
args = commandArgs(trailingOnly=TRUE)

# read csv file
#------------------------------------------------------------------
# read the merged file as a dataframe 
#------------------------------------------------------------------
total_numb_input_files = as.numeric(args[1])
outputFile = args[2]
foldchVar = ' Median '

ODIR = "scatterplots/"
dir.create(ODIR)

merged_df = read.csv( outputFile, header=TRUE,check.names=FALSE)
#print(head(merged_df)[1:5])
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
   res = str_match(x, "Analys [0-9]+-?[0-9]* .* Median (.*)Expt.*")[,2]
   res
}

#------------------------------------------------------------------
# fancy scatter plot matrix with R and p value
# https://www.r-bloggers.com/scatter-plot-matrices-in-r/
#------------------------------------------------------------------
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)

  # r squared #https://help.gooddata.com/display/doc/Covariance+and+Correlation+and+R-Squared
  r2 = r^2
  txt <- format(c(r2, 0.123456789), digits = digits)[1]
  txt <- paste("r-squared= ", txt, sep = "")
  text(0.5, 0.4, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  if(is.na(p)) p=1
  #print(p)
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.2, txt2)
}


# remove the pheno rows
SelecRownames = rownames(merged_df)[grep("^OTU_",rownames(merged_df))]
merged_df = merged_df[ SelecRownames ,]
cols_processed = 1
while(cols_processed < ncol(merged_df))  # loop over the merged dataset in steps of (number of input files)
{
   subset_df = merged_df[,cols_processed:(cols_processed+total_numb_input_files-1)]
   cols_processed = cols_processed + total_numb_input_files
   
   
   # keep the cols you need
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
	colnames(interestedSelecData) = c(1:total_numb_input_files)
	#print(head(interestedSelecData))
	#print("Next")
	ofile = paste(analys_numb, analys_name, sep='_')
	pdf(paste(ODIR, ofile, ".pdf", sep=''))
	#pairs(interestedSelecData)
	pairs(interestedSelecData, upper.panel = panel.cor)
	dev.off()
   }
   
}






#with(DF, plot(VAR1, VAR2))
#abline(fit <- lm(VAR2 ~ VAR1, data=DF), col='red')
#legend("topright", bty="n", legend=paste("R2 is",  format(summary(fit)$adj.r.squared, digits=4)))