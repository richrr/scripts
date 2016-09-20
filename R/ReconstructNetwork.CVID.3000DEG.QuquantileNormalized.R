# function: calculate correlation between two lists
# input files: 
	#	1.  two lists to generate pair for calculate correlation
	#	2.	a "gene expression" file
# output file:
	#	correlation coefficient,pvalue, fdr ;   for each pair of "genes"

# note:
#   if there are lines with the same "ID"(gene symble), the lines other than the first are deleted
# careful, BRB will change the name of the ids when doing class comparison, such as " <alpha subdivision>;uncultured alpha proteobacterium HF0070_14E07;"  "<alpha subdivision>" will be deleted. this will lead to the inconsistancy of the name in partner and expression file.
#==================================================================================================================
													#parameters
#==================================================================================================================
split= as.numeric(commandArgs(T)[1])		# which subset to run(input from perl)
NumberOfSplit= as.numeric(commandArgs(T)[2])	# split the whole pair combinations to smaller splits to run simutaniously when there are too many "genes", use with a perl script

#==================================================================================================================
#  network generation parameters
individualPvalueCutoff = 0.2
combinedPvalueCutoff = 2
combinedFDRCutoff = 2
#==================================================================================================================
# data import parameters
symbleColumnName = "UniqueID"     # the name of the column containing the "gene symble"
strainToDoList = c("CVID.GI","CVID.noGI","CVID.HealthyVolunteer")	# which column to analyze
#ExperimentToDoList =  c("control")
ExperimentToDoList =  c("CVID")
#==================================================================================================================
													#input Files 
#==================================================================================================================
expressionDataFile = "/capecchi/pharmacy/morgunlab/Dong/projects/CVID/network/hostGeneNetwork/input/expressionTable/GEO.RNASEQ.quantileNormalized.expression.txt"
list1File = "/capecchi/pharmacy/morgunlab/Dong/projects/CVID/network/hostGeneNetwork/input/geneList/geneList.3000DEGs.RNASEQ.quantileNormalized.txt"
list2File = "/capecchi/pharmacy/morgunlab/Dong/projects/CVID/network/hostGeneNetwork/input/geneList/geneList.3000DEGs.RNASEQ.quantileNormalized.txt"
FoldChangeFile = "/capecchi/pharmacy/morgunlab/Dong/projects/CVID/network/hostGeneNetwork/input/FoldChange/foldchange.3000DEGs.RNASEQ.quantileNormalized.txt"
annotationFile = "/capecchi/pharmacy/morgunlab/Dong/projects/CVID/network/hostGeneNetwork/input/geneAnnotations/geneDescription3000.unique.geneSymbol.txt"
#==================================================================================================================
outputFile = "/capecchi/pharmacy/morgunlab/Dong/projects/CVID/network/hostGeneNetwork/result/3000DEG.quantileNormalized."
#==================================================================================================================
#  R -q --slave <ReconstructNetwork.2014.12.06.Akk-pos-vs-Akk-neg-blocked.fromAndrey12-5-14.R --args 1 50000
#  SGE_Batch -c 'R -q --slave </capecchi/pharmacy/morgunlab/Dong/projects/diabetes/microbe-microbe-network/scripts/ReconstructNetwork.2014.12.19.Akk-pos-vs-Akk-neg-blocked.fromAndrey12-5-14.unionExpression.R --args $SGE_TASK_ID 50' -t 1-50 -f 10G -F 10G -r rnaSeqAkk

###########################################################################################################
#	output filename
networkFile = paste(outputFile,"individualPvalue:",individualPvalueCutoff,"_combinedPvalue:",combinedPvalueCutoff,"_combinedFDR:",combinedFDRCutoff,"-",split,"-in-",NumberOfSplit,".csv",sep="")
pairFile = "./pairs.3000DEG.quantileNormalized"


# generate pairlist to avoid re-generating pairs(needs time)
# load pair list data
list1 = read.csv(list1File,header = FALSE,sep="\t")
list2 = read.csv(list2File,header = FALSE,sep="\t")

genes1 = as.vector(list1[!duplicated(list1[,1]),1] )	#delete duplicated array probes
genes1 = as.character(genes1[order(genes1)])
genes2 = as.vector(list2[!duplicated(list2[,1]),1] )
genes2 = as.character(genes2[order(genes2)])

if (!file.exists(pairFile)){
	if ( sum(genes1 != genes2)==0){     # if gene list1 is the same as gene list2 then generate unique pairs
		pairs = t(combn(genes1,2))[,2:1]
	}else{
		pairs = expand.grid(genes1,genes2)  
	}	
	write.csv(pairs,pairFile,row.names=FALSE)
}else{
	pairs = read.csv(pairFile)
}

#==================================================================================================================
# generate subset start and end position
numberOfPairs = c()
splitSize = 1
begin = 1
end = 1
selectSubset = function(){
	numberOfPairs <<- nrow(pairs)
	splitSize <<- round(numberOfPairs/NumberOfSplit,0)+1
	begin <<- (split-1)*splitSize+1
	end <<- split*splitSize
}
#==================================================================================================================
													# done   parameters
#==================================================================================================================

#==================================================================================================================
# load expression data
#==================================================================================================================
expressionData = read.table(expressionDataFile,header = TRUE,stringsAsFactors=FALSE,sep="\t",quote='"',comment.char = "",check.names=F)
expressionData = expressionData[!duplicated(expressionData[,symbleColumnName])&!is.na(expressionData[,symbleColumnName]),]    # delete duplicated "genes"
rownames(expressionData) = as.character(expressionData[,symbleColumnName])
head(expressionData)[1:5]

#==================================================================================================================
#								calculate all pair correlation,pvalue,fdr for GF,Control in balbc/b10a
#==================================================================================================================

# function : calculate correlation,pvalue,fdr for all pairs between list 1 and list 2,  for strain and experiment
# input:
#     1.strain : a string
#     2. experiment : a string
#     3. list 1 of "genes" : a list of "gene symble" character
#     4. list 2 of "genes": a list of "gene symble" character
#       5. data : a dataframe containing the data    
# output : a table containing pair information, correlation,pvalue,fdr
calculateCorrelation = function (strain,experiment,pairs,data){

    # calculate for each pair
#	cat ("pair is \t",pairs[1,],"\n\n")   
	print (paste(">>>>>>",pairs[1,1],sep=""))
    out = apply(pairs,1,ForCor,strain,experiment,data)  
    #cat("here")
    #print ("\n")
    colnames(out) = paste(as.vector(pairs[,1]),as.vector(pairs[,2]),sep="<==>")
    out = t(out)
    colnames(out) = c(paste(strain,experiment,"correlation.Coefficient",sep="."),paste(strain,experiment,"pvalue",sep="."))
    
	#calculate FDR
	FDR = p.adjust(out[,colnames(out)[grep("pvalue",colnames(out))]],method="fdr")
	oldColnames = colnames(out)
	out = cbind(out, FDR)
	colnames(out) = c(oldColnames,paste(strain,experiment,"FDR",sep="."))
	pairName = rownames(out)
	out = cbind(pairName,out)
    out
}
combineEachExperiment = function(strainExperimentPair,pairs){
	#outEachExperiment = calculateCorrelation(strain=strainExperimentPair[1],experiment=strainExperimentPair[2],list1=head(list1),list2=head(list2),data=expressionData)
	outEachExperiment = calculateCorrelation(strain=strainExperimentPair[1],experiment=strainExperimentPair[2],pairs=pairs,data=expressionData)
	if (length(out)==0){
		out <<- outEachExperiment
	}else{
		out <<- merge(out,outEachExperiment)
	}
}

main = function(strainToDoList,ExperimentToDoList,pairs){
			experimentCombination = expand.grid(strainToDoList,ExperimentToDoList)
			t=apply(experimentCombination,1,combineEachExperiment,pairs)
			#out	
}

#     function: calculate correlation coefficient of a pair of genes
# input:
#         1.pair: a pair of "gene"s:
#                 value: a vector of length 2, containing names of the two "genes"
#         2. strain,experiment,data
# output correlation coefficient: vector of length 2: "correlation Coefficient","pvalue"
ForCor = function(pair,strain,experiment,data){    
	#print (paste(">>",pair[1],sep=""))
	
    pair = as.vector(pair)
    pair = gsub(" ","",pair)  # when using "apply" function, R will insert space in front of the shorter names, introducing error, so make this correction
	#print (paste(">>",pair[1],sep=""))
    n<<- n+1
    cat (n)
    cat (":\t")
    calculating = paste(strain,experiment,paste(rev(pair),collapse=" ___ "),sep=" --> ")
	cat (calculating)
    cat ("\n")
   #########################################
   # find which column to calculate
   ##########################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    columnsForStrain = colnames(data)[grep(strain,colnames(data))]
    collumnsForStrainAndExperiment = columnsForStrain[grep(experiment,columnsForStrain)]
    #print(strain)
   # print(colnames(data))
   # print(columnsForStrain)
   # print(collumnsForStrainAndExperiment)
    c1 = as.numeric(
                    as.vector(data[as.character(pair[1]),
                                   collumnsForStrainAndExperiment
                                   ]
                    )
                )
    c2 = as.numeric(
                    as.vector(data[as.character(pair[2]),
                                   collumnsForStrainAndExperiment
                                   ]
                    )
                )
  #  print(row.names(head(data)))
   # print(colnames(head(data)))
  #  print(collumnsForStrainAndExperiment)
    
  #  print(as.character(pair[1]))
   
    if (sum(!is.na(c1+c2))<3){
    	cat ("not enough observation\n\n\n>>>>>>>>>>>>>>>>>>>\n",pair[1],":\n",c1,"\n",pair[2],":\n",c2,"\n>>>>>>>>>>>>>>>>\n")
        return (c(NA,NA))
    }
    #error()
    p = cor.test(c1,c2)
    outLine = as.matrix(c(p$estimate,p$p.value))
    rownames(outLine) = c("correlation.Coefficient","pvalue")
    #colnames(outLine) = paste(x,collapse=",")
    outLine
}


#####################################################################################################################################################
#test a small dataset
#==================================================================================================================
test = function (){	
	split <<- 1
	NumberOfSplit <<- 100000
}
#test()
selectSubset()
# pick subset to calculate	
pairsTodo = pairs[begin:min(end,nrow(pairs)),]	
#print ("row names!!!!!!!!!!!!!!!")
#==================================================================================================================
n=0
out = c()
main(strainToDoList,ExperimentToDoList,pairsTodo)

#write.csv(out, paste(networkFile,"correlation",sep="-"),row.names=FALSE)
#write.csv(out,paste(outfile,"_",begin,"_",end,".csv",sep=""),row.names=FALSE)
#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# calculate PUC 
	# 1. calculate PUC from correlation coefficient and foldchange information
	# 2. calculate combined pvalue from different datasets
	# 3. generate Combined FDR from combined pvalue
# input File:
	# 1. output from 20130417CalculationCorrelation: correlation coefficient, pvalue, fdr for each strain(datasets)
# output:
	#  combined pvalue, combined fdr,PUC added to (correlation coefficient, pvalue, fdr for each dataset,)
#strainToDoList = c("SW","B6","B10A")
#dataName = paste(strainToDoList,ExperimentToDoList,sep=".")


#ResultFile = "/projects1/antibiotic_mice/dong/scripts/network/result/forMatt/test_1_92-PUC.csv"
Data = out
# calculate PUC
#==================================================================================================================
forPUC = function(FoldChangeFile){
	
	#change out format for PUC
	library(stringr)
	pair = str_split(Data[,"pairName"],"<==>")
	pairs = t(as.data.frame(pair))
	#pairs = str_replace_all(pairs,perl("\\.\\d"),"")
	colnames(pairs) = c("partner1","partner2")
	outForPUC = cbind(pairs,Data)
	head(outForPUC)
#==================================================================================================================
#    						generate metabolite foldChange
	if (FoldChangeFile == "NotUsingPUC"){
		
	}else {
		FoldChangeMetabolic =  read.table(FoldChangeFile,header = TRUE,stringsAsFactors=FALSE,sep="\t",quote="",comment.char = "")
		head(FoldChangeMetabolic)
	 	FoldChangeMetabolic = FoldChangeMetabolic[!duplicated(FoldChangeMetabolic[,"ID"]),]
		rownames(FoldChangeMetabolic) = FoldChangeMetabolic[,"ID"]
	
		FoldMetab1_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner1"]),c("ID","FoldChange")]
		colnames(FoldMetab1_InPair) = c("partner1InFold","partner1_FoldChange")
		head(FoldMetab1_InPair)
	
		FoldMetab2_InPair = FoldChangeMetabolic[as.vector(outForPUC[,"partner2"]),c("ID","FoldChange")]
		colnames(FoldMetab2_InPair) = c("partner2InFold","partner2_FoldChange")
		head(FoldMetab2_InPair)
	
		outForPUC = cbind(outForPUC,FoldMetab1_InPair,FoldMetab2_InPair)
	}
#==================================================================================================================
#    					calculate PUC
	coefficientColnames = colnames(outForPUC)[grep("Coefficient",colnames(outForPUC))]
#	interestedCoefficientColnames = coefficientColnames[sapply(dataName,grep,coefficientColnames)]
	interestedCoefficientColnames = coefficientColnames
	interestedCorrelationData = outForPUC[,interestedCoefficientColnames]
	interestedCorrelationData = as.matrix(interestedCorrelationData)
	# calculate correlation Direction For correlation coefficient of interest
	
	interestedCorrelationData = apply(interestedCorrelationData,2,function(x){as.numeric(as.vector(x))})
	signOfInterestedCorrelationData = interestedCorrelationData/abs(interestedCorrelationData)					# forOutput
	rownames(signOfInterestedCorrelationData) = c()
	colnames(signOfInterestedCorrelationData) = paste(colnames(interestedCorrelationData),"correlationDirection",sep=".")
	# see if direction match in interest
	ifMatch = apply(signOfInterestedCorrelationData,1
					,function(tt){
								uniqueValues = unique(tt[!is.na(tt)])
								NumberOfUniqueValues = length(uniqueValues)
								if (length(tt[!is.na(tt)])<2){
									return (uniqueValues[1])
								}
								if (NumberOfUniqueValues>1){
									return (0)
								}else if(NumberOfUniqueValues==0){
									return (NA)
								}else{
									return (1*uniqueValues[1])
								}
							}
					)	# if each row of signOfInterestedCorrelationData is like -1,1,-1; if match, "ifMatch" will be 3 or -3
	names(ifMatch) = c()
	head(ifMatch)																									# forOutput	
	# calculate matched correlation direction
	matchedExpressionDirection = ifMatch*signOfInterestedCorrelationData[,1]
	names(matchedExpressionDirection) = c()
	matchedExpressionDirection	= ifMatch																			# forOutput	
	# calculate fold change direction	
	if (FoldChangeFile == "NotUsingPUC"){
		outForPUC = cbind(outForPUC,signOfInterestedCorrelationData,ifMatch,matchedExpressionDirection)
	}else {
		FoldChangeColnames = colnames(outForPUC)[grep("FoldChange",colnames(outForPUC))]
		FoldChangeData = outForPUC[,FoldChangeColnames]
		FoldChangeDirection = (FoldChangeData-1)/abs(FoldChangeData-1)
		names(FoldChangeDirection) = c()
		colnames(FoldChangeDirection) = paste(colnames(FoldChangeData),"FoldChangeDirection",sep=".")
		FoldChangeDirection																							## forOutput	
		# calculate if fold change direction are the same for the two partners
		IfFoldChangeDirectionMatch = apply(FoldChangeDirection,1,prod)
		names(IfFoldChangeDirectionMatch) = c()
		IfFoldChangeDirectionMatch																					## forOutput	
		# use "matched correlation direction" and "if fold change direction are the same" to calculate PUC
		PUC = IfFoldChangeDirectionMatch * matchedExpressionDirection												## forOutput	
		outForPUC = cbind(outForPUC,signOfInterestedCorrelationData,ifMatch,matchedExpressionDirection,FoldChangeDirection,IfFoldChangeDirectionMatch,PUC)
	}
	head(outForPUC)

	
	# calculate combined Pvalue for interest
	PvalueColnames = colnames(outForPUC)[grep("pvalue",colnames(outForPUC))]
#	interestedPvalueColnames = PvalueColnames[sapply(dataName,grep,PvalueColnames)]
	interestedPvalueColnames = PvalueColnames
	interestedPvalueData = outForPUC[,interestedPvalueColnames]
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
	outForPUC = cbind(outForPUC,combinedPvalue)
	outForPUC
		
} 

#=============================================================================================

result = forPUC(FoldChangeFile)

#calculate FDR for combined pvalue
combinedFDR = p.adjust(result[,"combinedPvalue"],method="fdr")
result = cbind(result,combinedFDR)
annotation = read.table(annotationFile,header=TRUE,stringsAsFactors=FALSE,sep="\t",quote="",comment.char = "")
result = merge(result,annotation,by.x = "partner1",by.y = "UniqueID",suffixes = c("",".partner1"))
result = merge(result,annotation,by.x = "partner2",by.y = "UniqueID",suffixes = c("",".partner2"))

write.csv(result, paste(networkFile,"PUC",sep="-"),row.names=FALSE)

###########################################################################################################
#generate network
###########################################################################################################
data = result
generateNetwork = function(){
	#generate network After Filter

	# find PUC expected and significant in combinedPvalue and significant in combinedFDR
	if (FoldChangeFile == "NotUsingPUC"){
		out = data
	}else {
		out = data
#		out = data[data[,"PUC"]==1 ,]
	}
	out = out[abs(out[,"ifMatch"]) ==1 & !is.na(out[,"ifMatch"]),]
	out = out[(out[,"combinedFDR"]<combinedFDRCutoff)==1,]
	# find significant in individual pvalue: all of the pvalues for all of the datasets for each pair must be smaller than threshold
	# find pvalue data
	pvalueData = out[,grep("pvalue",colnames(out))]
	pvalueData = as.matrix(pvalueData)
	# calculate the largest pvalue among all datasets for each gene, this smallest pvalue must be smaller than threshold
	passIndevidualPvalue = apply(pvalueData,1,max,na.rm=TRUE)<individualPvalueCutoff
	outNetwork = out[passIndevidualPvalue,]


	write.csv (outNetwork,networkFile,row.names=FALSE)
}
###########################################################################################################


generateNetwork()


#file.remove(pairFile)

addAnnotation = function(){
	if(!exists("result")){
		#result = read.csv(paste(networkFile,"PUC",sep="-"),header=TRUE,strinsAsFactors=FALSE)
		result = read.csv("/capecchi/pharmacy/morgunlab/Dong/projects/diabetes/microbe-microbe-network/result/OverlapindividualPvalue:2_combinedPvalue:2_combinedFDR:2-1-in-1.csv-PUC",header=TRUE,stringsAsFactors=FALSE)
	}
	head(result)
	annotation = read.table("/capecchi/pharmacy/morgunlab/Dong/projects/diabetes/microbe-microbe-network/input/geneAnnotations/annotations-Uid_symbl.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t",quote="",comment.char = "")
	head(annotation)
	result.annotation = merge(result,annotation,by.x = "partner1",by.y = "UniqueID",suffixes = c("",".partner1"))
	head(result.annotation)
	result.annotation = merge(result.annotation,annotation,by.x = "partner2",by.y = "UniqueID",suffixes = c("",".partner2"))
	head(result.annotation)
	write.csv(result.annotation,"/capecchi/pharmacy/morgunlab/Dong/projects/diabetes/microbe-microbe-network/result/OverlapindividualPvalue:2_combinedPvalue:2_combinedFDR:2-1-in-1.csv-PUC-annotation",row.names=FALSE)

	#write.csv(result,"/capecchi/pharmacy/morgunlab/Dong/projects/diabetes/microbe-microbe-network/result/OverlapindividualPvalue:2_combinedPvalue:2_combinedFDR:2-1-in-1.csv-PUC",row.names=FALSE)
}




#  SGE_Batch -c "R -q --slave < ReconstructNetwork.2014.09.20.illumina7new.33new.222old.R --args 1 100000" -f 10G -F 10G -r OMEcolor



















