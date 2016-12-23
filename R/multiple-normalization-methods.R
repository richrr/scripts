######Supplemental Material to 
######State-of-the-art normalizations improve NMR-based metabolomics analysis
######Stefanie M. Kohl, Matthias S. Klein, Peter J. Oefner, Rainer Spang, Wolfram Gronwald 
###### 11306_2011_350_MOESM2_ESM

######In the Following, the R code for the normalizations used in the paper is given.


######The non-normalized data matrix of feature intensities is stored in a matrix called "data" 
######Each row represents a feature, each column a sample


library(argparser)
library(affy)

#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Perform all or the requested normalization.')
p <- add_argument(p, "--file", help="file which has non-normalized data of feature intensities. Each row represents a feature, each column a sample", nargs=1) # required; but written as optional format so I explicitly mention the "--file"
p <- add_argument(p, "--output", help="output file", default="./normalizations/res_")
p <- add_argument(p, "--symbolColumnName", help="the name of the column containing the measurement labels in expression data file", default="#SampleID") # 'UniqueID'

# explicitly mention atleast one of the following
p <- add_argument(p, "--all", help="Run all Normalizations", flag=TRUE)
# OR
p <- add_argument(p, "--pq", help="Probabilistic Quotient Normalization", flag=TRUE)
p <- add_argument(p, "--cl", help="Cyclic Loess Normalization", flag=TRUE)
p <- add_argument(p, "--c", help="Contrast Normalization", flag=TRUE)
p <- add_argument(p, "--q", help="Quantile Normalization", flag=TRUE)
p <- add_argument(p, "--lb", help="Linear Baseline Normalization", flag=TRUE)
p <- add_argument(p, "--lw", help="Li-Wong Normalization", flag=TRUE)
p <- add_argument(p, "--cs", help="Cubic Spline Normalization", flag=TRUE)
p <- add_argument(p, "--a", help="Auto Scaling", flag=TRUE)
p <- add_argument(p, "--p", help="Pareto Scaling", flag=TRUE)
p <- add_argument(p, "--v", help="Variance Stabilization Normalization", flag=TRUE)




argv <- parse_args(p)
#print (p)

outputFile = argv$output
symbleColumnName =  argv$symbolColumnName

# set all normalization flags to TRUE
if(argv$all){
    print("Requested ALL normalizations!")
    argv$pq = T
    argv$cl = T
    argv$c = T
    argv$q = T
    argv$lb = T
    argv$lw = T
    argv$cs = T
    argv$a = T
    argv$p = T
    argv$v = T
}


# create if outdir doesn't exist:
odir =  strsplit(outputFile, "/")[[1]]
odir = paste(odir[1:(length(odir) -1)], collapse='/')
print(odir)
dir.create(odir, recursive = TRUE)


# read the file and create a matrix
expressionData = read.csv( argv$file, header=TRUE, check.names=FALSE, na.strings=c("","na","NA", "Na", "NaN"), row.names= symbleColumnName)
data = data.matrix(expressionData)


# output to file
send_to_write = function(matrix, ofile, delim=','){
    write.table(matrix, ofile, sep=delim)
}


#####---Probabilistic Quotient Normalization-----##########################
if(argv$pq){
        print("Running PQN")
	reference <- apply(data,1,median)
	quotient <- data/reference
	quotient.median <- apply(quotient,2,median)
	pqn.data <- t(t(data)/quotient.median)
	send_to_write(pqn.data, paste(outputFile, ".txt", sep='pqn'))
}



#####---Cyclic Loess Normalization-----####################################
if(argv$cl){
        print("Running Loess")
	loess.data <- normalize.loess(data, 
						subset=1:nrow(data), 
						epsilon=10^-2, 
						maxit=2, 
						log.it=FALSE, 
						verbose=TRUE, 
						span=0.75, 
						family.loess="gaussian")
	send_to_write(loess.data, paste(outputFile, ".txt", sep='loess'))
}



#####---Contrast Normalization-----########################################
if(argv$c){
        print("Running Contrast")
	#---First adaption: Make the data matrix non-negative
	smallvalue <- function(x, threshold=1e-11){				#threshold was chosen such that it is sufficiently smaller than the data
		for(i in 1:length(x)){
			if(!x[i]>0)	x[i] <- threshold
		}
	}

	nonnegative.data=smallvalue(data)

	#---Apply normalization
	maffy.data <- maffy.normalize(nonnegative.data,
						subset=1:nrow(nonnegative.data),
						span=0.75,
						verbose=TRUE,
						family="gaussian",
						log.it=FALSE)

	#---Second adaption: Subtract 10% Quantile from each sample
	subtract <- function(x){
		t(t(x)-apply(x,2,quantile,0.1))
	}

	contrast.data <- subtract(maffy.data)
	send_to_write(contrast.data, paste(outputFile, ".txt", sep='contrast'))
}



#####---Quantile Normalization-----########################################
if(argv$q){
        print("Running Quantile")
	normalize.quantile <- get("normalize.quantiles", en=asNamespace("affy"))
	quantile.data <- normalize.quantile(data)
	send_to_write(quantile.data, paste(outputFile, ".txt", sep='quantile'))
}



#####---Linear Baseline Normalization-----#################################
if(argv$lb){
        print("Running Linear Baseline")
	linear.baseline <- apply(data,1,median)
	baseline.mean <- mean(linear.baseline)
	sample.means <- apply(data,2,mean)
	linear.scaling <- baseline.mean/sample.means
	linear.baseline.data <- t(t(data)*linear.scaling)
	send_to_write(linear.baseline.data, paste(outputFile, ".txt", sep='linear.baseline'))
}



#####---Li-Wong Normalization-----#########################################
if(argv$lw){
        print("Running Li-Wong")
	#---First step: Find baseline sample
	average.intensity <- apply(data,2,mean)
	median.number <- round(ncol(data)/2 + 0.1)				#R has an add way of rounding. 
												#the additional 0.1 ensures that it rounds properly
	ordering <- order(average.intensity)
	median.sample.number <- ordering[median.number]
	median.sample <- data[,median.sample.number]

	#---Apply normalization
	liwong.data=vector()
	for(i in 1:ncol(data)){
		liwong.model <- normalize.invariantset(data=data[,i],
						ref=median.sample,
						prd.td=c(0.003,0.007))		#the threshold of the rank-invariant set might need to be adjusted from case to case
		liwong.sample <- predict(liwong.model$n.curve$fit,		#chosen such that the rank-invariant set it sufficiently large
						data[,i])
		liwong.data <- cbind(liwong.data,liwong.sample$y)
	}
	send_to_write(liwong.data, paste(outputFile, ".txt", sep='liwong'))
}



#####---Cubic Spline Normalization-----####################################
if(argv$cs){
        print("Running Spline")
	spline.data <- normalize.qspline(data,
					samples=0.02,
					target=apply(data,1,mean))
	send_to_write(spline.data, paste(outputFile, ".txt", sep='spline'))
}



#####---Auto Scaling-----##################################################
if(argv$a){
        print("Running Auto")
	centered.data <- data - apply(data,1,mean)
	scaling.auto <- apply(data,1,sd)
	auto.data <- centered.data/scaling.auto
	send_to_write(auto.data, paste(outputFile, ".txt", sep='auto'))
}



#####---Pareto Scaling-----################################################
if(argv$p){
        print("Running Pareto")
	centered.data <- data - apply(data,1,mean)
	scaling.pareto <- sqrt(apply(data,1,sd))
	pareto.data <- centered.data/scaling.pareto
	send_to_write(pareto.data, paste(outputFile, ".txt", sep='pareto'))
}



#####---Variance Stabilization Normalization (VSN)-----####################
if(argv$v){
        print("Running VSN")
	library(vsn)									#load package unless it is already loaded
	vsn.model <- vsn2(data)
	vsn.data <- predict(vsn.model,data)
	send_to_write(vsn.data, paste(outputFile, ".txt", sep='vsn'))
}


print("DONE")
