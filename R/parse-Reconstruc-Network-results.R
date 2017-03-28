##### to do: 
## how to merge multiple miseq expts?

library(argparser)
library(stringr)
library(hash)


# /nfs1/Morgun_Lab/richrr/Abx/MiSeq1/analysis_read1/OTU_ClosedRef/Abx/network_analysis/result
 
# Rscript /nfs1/Morgun_Lab/richrr/scripts/R/parse-Reconstruc-Network-results.R --comparison network_comp-output.csv --correlation network_corr-output.csv --AnalysToDoList ../AnalysToDoList.txt

# Rscript /nfs1/Morgun_Lab/richrr/scripts/R/parse-Reconstruc-Network-results.R ../miseq10-pheno-measurm.csv --correlation network_corr-output.csv --AnalysToDoList ../AnalysToDoList.txt --mapFile ../miseq10-mapping-file.txt --mapColumns #SampleID Sample_typeTime_pointTreatment_per_Abx --symbolColumnName #SampleID --pval 0.005

# Rscript /nfs1/Morgun_Lab/richrr/scripts/R/parse-Reconstruc-Network-results.R ../miseq10-pheno-measurm.csv --correlation network_corr-output.csv --AnalysToDoList ../AnalysToDoList.txt --mapFile ../miseq10-mapping-file.txt --mapColumns #SampleID Sample_typeTime_pointTreatment_per_Abx --symbolColumnName #SampleID --pval 0.05 --excludesamepartnerpairs

# /nfs1/Morgun_Lab/richrr/Abx/MiSeq2/analysis_read1/OTU_ClosedRef/Abx/network_analysis/result
# Rscript /nfs1/Morgun_Lab/richrr/scripts/R/parse-Reconstruc-Network-results.R ../miseq11-pheno-measurm.csv --correlation network_corr-output.csv --AnalysToDoList ../AnalysToDoList.txt --mapFile ../miseq11-mapping-file.txt --mapColumns #SampleID Sample_typeTime_pointTreatment_per_Abx --symbolColumnName #SampleID --pval 0.05 --excludesamepartnerpairs


#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Parses individual result files and plots correlations of pairs that pass the pvalue cutoff. Mainly useful for correlation results.')
p <- add_argument(p, "expressionDataFile", help="file containing values (measurements) to be used for correlation") # required

## comparison or correlations file
p <- add_argument(p, "--comparison", help="results are from comparison") 
p <- add_argument(p, "--correlation", help="results are from correlation") 
p <- add_argument(p, "--pval", help="pvalue threshold", default=0.005, type="numeric") 
p <- add_argument(p, "--AnalysToDoList", help="AnalysToDoList", default="AnalysToDoList.txt") # tab-delimited; perform correlations or comparisons as per inmput
p <- add_argument(p, "--symbolColumnName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID") 
p <- add_argument(p, "--mapFile", help="tab delimited mapFile", default="mapFile.txt") # use the mapping file to identify the samples group affiliation to perform correlations
p <- add_argument(p, "--mapColumns", nargs=2, help="Columns in map file in which to search the sample group associations", default=c("#SampleID","Experiment")) # use the columns in the mapping file to identify the samples group affiliation to perform correlations
p <- add_argument(p, "--excludesamepartnerpairs", help="exclude same partner pairs", flag=TRUE)


argv <- parse_args(p)
#print (p)


sig_thresh = argv$pval

if(is.na(argv$comparison) && is.na(argv$correlation))
{
  print("Either comparison or correlation result file is required. Parses individual result files and plots correlations of pairs that pass the pvalue cutoff. Mainly useful for correlation results.")
  quit()
}


symbleColumnName =  argv$symbolColumnName
AnalysToDoList = read.delim(argv$AnalysToDoList,header = FALSE)

#==================================================================================================================
# load expression data
#==================================================================================================================
expressionDataFile = argv$expressionDataFile
expressionData = read.csv(expressionDataFile,header = TRUE,check.names=FALSE)
expressionData = expressionData[!duplicated(expressionData[,symbleColumnName]),]    #delete duplicated measurements (genes, taxa labels, phenotypic labels, etc.)
rownames(expressionData) = expressionData[,symbleColumnName]



# the samples in samplIdCol belong to the groups in exptCol
samplIdCol = argv$mapColumns[1]
exptCol = argv$mapColumns[2]

# create a group_sampleid dictionary
mapFile=read.delim(argv$mapFile,header = TRUE,sep="\t",check.names=FALSE)
dict = hash()
# for unique items in the experiment column
for(k in levels(mapFile[, exptCol])){
    # extract the rows whose col value in exptCol matches k
    # assign the sample names (i.e. values in samplIdCol) to key k
    dict[k] = as.vector(mapFile[which(mapFile[,exptCol]==k), samplIdCol]) 
}



# check the number and order of samples assigned to groups in dict. MAY be as per the order in which the samples are added (order in the mapping file)
for(k in keys(dict)){
    str_r_l = paste(k, length(dict[[k]]), "values", paste(as.vector(dict[[k]]), collapse=", "), sep="---") 
    print(str_r_l)
}



if(FALSE){
for(indx in 1:nrow(AnalysToDoList)){
    c = as.vector(AnalysToDoList[indx,])
    #print(c)
    if(c[1] != "" && c[2] == "" && c[4] == "correlation"){
        # conditions 6,7,9
        if(c[3] == "delta"){
            categA= toString(c[1,"V1"])

            cA = strsplit(categA, ",")[[1]]
            categ1 = cA[1]
            categ2 = cA[2]
           
            interest_column = paste(categ1, "minus", categ2, sep=" ")
            print(interest_column)
        } else{
            categ = toString(c[1,"V1"])
            interest_column = categ
            print(interest_column)

        }
    }
}
}# end false


if(!is.na(argv$comparison)){
    inFileComp = argv$comparison

    inFileCompData = read.csv(inFileComp,header = TRUE,check.names=FALSE, row.names=1)
    #print(head(inFileCompData))

    # selecting only the pvalue columns
    subset_data = inFileCompData[, grep("pvalue", colnames(inFileCompData))]
    #print(head(subset_data))

    for(column  in 1:ncol(subset_data)){
       print("---------Analysis---------")
       print(colnames(subset_data)[column])
       sign_df  = subset_data[which(subset_data[,column] < sig_thresh), column, drop=FALSE]
       #print(sign_df)
       print(rownames(sign_df))

    }

}


#http://stackoverflow.com/questions/9455167/converting-data-frame-factors-with-data-matrix-weird-results
#else as.numeric gives weird numbers instead of actual values
numericizeVector = function(c){
    nc = sapply(c, function(x) as.numeric(as.character(x)))
    return(nc)
}


extract_raw_vals = function(pair, edata, c, dict, pVal, cOrr){

    idxs = as.vector(dict[[c]])

    #n<<- n+1
    #cat (n)
    #cat (":\t")
    #calculating = paste(rev(pair),collapse=" ___ ")
    #cat (calculating)
    #cat ("\n")

    c1 = numericizeVector(as.vector(edata[pair[1], idxs]))
    c2 = numericizeVector(as.vector(edata[pair[2], idxs]))

    #print(c1)
    #print(c2)

    cOrr = signif(cOrr,3)
    pVal = signif(pVal,3)

    plot(c1, c2, xlab=pair[1], ylab=pair[2])
    label = paste(c("Correl=", cOrr, "p-value=", pVal), collapse=" ")
    mtext(label, side = 3, adj = 0.05, line = -1.3)

}


extract_raw_vals_delta = function(pair, edata, c, dict, pVal, cOrr){

    cA = strsplit(c, ",")[[1]]
    c_1 = cA[1]
    c_2 = cA[2]


    #n<<- n+1
    #cat (n)
    #cat (":\t")
    #calculating = paste(rev(pair),collapse=" ___ ")
    #cat (calculating)
    #cat ("\n")

    idxs1 = as.vector(dict[[c_1]])
    idxs2 = as.vector(dict[[c_2]])

    c1 = numericizeVector(as.vector(edata[pair[1], idxs1]))
    c2 = numericizeVector(as.vector(edata[pair[1], idxs2]))

    c3 = numericizeVector(as.vector(edata[pair[2], idxs1]))
    c4 = numericizeVector(as.vector(edata[pair[2], idxs2]))

    ca = c1-c2
    cb = c3-c4

    #print(ca)
    #print(cb)

    cOrr = signif(cOrr,3)
    pVal = signif(pVal,3)

    plot(ca, cb, xlab=pair[1], ylab=pair[2])
    label = paste(c("Correl=", cOrr, "p-value=", pVal), collapse=" ")
    mtext(label, side = 3, adj = 0.05, line = -1.3)

}


if(!is.na(argv$correlation)){
    inFileCorr = argv$correlation

    inFileCorrData = read.csv(inFileCorr,header = TRUE,check.names=FALSE, row.names=1)
    #print(head(inFileCorrData))

    # selecting only the pvalue columns
    subset_data = inFileCorrData[, grep("pvalue", colnames(inFileCorrData)), drop=FALSE]
    # unselect comparisons of delta or normal correlations
    subset_data = subset_data[, grep(" vs ", colnames(subset_data), invert=TRUE), drop=FALSE]
    # unselect comparisons of timelag correlations
    #print(head(subset_data))
    subset_data = subset_data[, grep(" timelag ", colnames(subset_data), invert=TRUE), drop=FALSE]
    #print(head(subset_data))

    for(column  in 1:ncol(subset_data)){
       # column name
       print("---------Analysis---------")
       analy_colname = toString(colnames(subset_data)[column])
       print(analy_colname)

       # the analysis for which the above column was generated
       res=as.numeric(str_match(analy_colname, "Analys ([0-9]+) .*")[,2])
       print(AnalysToDoList[res,])
       categ_name = toString(AnalysToDoList[res,"V1"])
       print(categ_name)

       # identify which "gene"-"gene" pairs are significant for the given analysis
       #sign_df  = subset_data[which(subset_data[,column] < sig_thresh), column, drop=FALSE]
       sign_df  = subset_data[which(subset_data[,analy_colname] < sig_thresh), analy_colname, drop=FALSE]
       
       print(rownames(sign_df))
       #print(sign_df)


       total_numb_sig_pairs = length(rownames(sign_df))
       #print(total_numb_sig_pairs)

       # discard cases where partner1 = partner2 (same gene pair)
       if(argv$excludesamepartnerpairs){
           diff_partners = c()
           for(r in 1: total_numb_sig_pairs){
               print("Row number")
               print(r)
               rown=rownames(sign_df)[r]
               partner1 = str_split(rown,"<==>")[[1]][1]
               partner2 = str_split(rown,"<==>")[[1]][2]
               print(partner1)
               print(partner2)
               if(partner1 != partner2)
               {
                   diff_partners = c(diff_partners, rown)
               }
           }
           print(diff_partners)
           sign_df  = sign_df[diff_partners, , drop=FALSE]
           total_numb_sig_pairs = length(rownames(sign_df))
       }

       print("Total number of subfigures to be generated for the category")
       print(total_numb_sig_pairs)
       #print(sign_df)
       # no need to generate the figure
       if(total_numb_sig_pairs==0 || total_numb_sig_pairs>10){
           next
       }

       #jpeg('rplot.jpg') 
       #png(paste(categ_name , '.png' , sep=''), width=1000, height=1200) 
       pdf(paste(categ_name , '.pdf' , sep=''), width=10, height=20) 

		
	
       if(total_numb_sig_pairs>1){

         # arrange plots in 2 column plots
         cols_in_fig = 2
         if(total_numb_sig_pairs>10) # arrange plots in 4 column plots
         {
             cols_in_fig = 4
         }
         # arrange plots in 2 column plots with rows_in_fig rows
         rows_in_fig = ceiling(total_numb_sig_pairs/cols_in_fig)
         par(mfrow=c(rows_in_fig,cols_in_fig))
       } else {
         par(mfrow=c(1,1))
       }

       # sets the margin sizes in the following order: bottom, left, top, and right. default c(5.1, 4.1, 4.1, 2.1)
       #par(mar=c(1.1, 1.1, 1.1, 1.1)) # reduces margins
       
       for(r in 1: total_numb_sig_pairs){
           #print("Row number")
           #print(r)
           rown=rownames(sign_df)[r]
           # pvalue
           pVal = sign_df[rown,analy_colname]
           # correlation
           analy_colname_correl_coeff = sub("pvalue", "Coefficient", analy_colname)
           cOrr = inFileCorrData[rown, analy_colname_correl_coeff]
           

           partner1 = str_split(rown,"<==>")[[1]][1]
           partner2 = str_split(rown,"<==>")[[1]][2]
           #print(partner1)
           #print(partner2)
           if(grepl(',', categ_name)){
               extract_raw_vals_delta(c(partner1, partner2), expressionData, categ_name, dict, pVal, cOrr)
           } else{
               extract_raw_vals(c(partner1, partner2), expressionData, categ_name, dict, pVal, cOrr)
           }
       }
       dev.off()
       
    }

}



