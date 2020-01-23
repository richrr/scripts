library(argparser)

p <- arg_parser('Calc basic stats about the network')
p <- add_argument(p, "--file", help="file which has merged correlations from multiple expts", nargs=1) # required; but written as optional format so I explicitly mention the "--file"
p <- add_argument(p, "--output", help="output file", default='-stats-log.txt')
p <- add_argument(p, "--noPUC", help="did not calc. PUC", flag=TRUE) 

argv <- parse_args(p)
#print (p)


if(length(argv$file) < 1)
{
  print("At least 1 file is required. Give --file file ... in cmd")
  quit()
}

outputFile = argv$output

noPUC = argv$noPUC

data = read.csv( argv$file, header=TRUE, check.names=FALSE, row.names=1)
head(data)



#----------------------------------------------------------------------
# calc stats of the generated network
#----------------------------------------------------------------------
calc_stats = function(inNet, noPUC, correlThreshold=0){

        
        combcoeffcol = as.numeric(as.vector(inNet[,"combinedCoefficient"]))
        # not sure why it behaves this way!
        lowest_neg_corrl = min(combcoeffcol[combcoeffcol < 0])
        highest_neg_corrl = max(combcoeffcol[combcoeffcol < 0])

        lowest_pos_corrl = min(combcoeffcol[combcoeffcol > 0])
        highest_pos_corrl = max(combcoeffcol[combcoeffcol > 0])
        
        print(lowest_neg_corrl)
        print(highest_neg_corrl)
        print(lowest_pos_corrl)
        print(highest_pos_corrl)


        
        coefficientData = inNet[,grep("correlationDirection",colnames(inNet)),drop=FALSE]
		
		## artifically add row (with 0) incase it is only 1 row
		if(nrow(coefficientData) == 1){
			tmp_v = c("0")
			tmp_artif = data.frame(tmp_v)
			colnames(tmp_artif) = c("combinedCoefficient.correlationDirection")
			coefficientData = rbind(coefficientData, tmp_artif)
		}
		
		coefficientData = apply(coefficientData,2,function(x){as.numeric(as.vector(x))})
		
        res_pos = apply(coefficientData, 2, function(x) sum(x == 1))
        res_neg = apply(coefficientData, 2, function(x) sum(x == -1))
        res = cbind(res_pos, res_neg)
        #print(res)
        
        
        print(paste(c("Total number of edges: ", nrow(inNet)), collapse="")) 

        setpartner1 = unique(as.vector(inNet[,"partner1"]))
        setpartner2 = unique(as.vector(inNet[,"partner2"]))

        nodes = union(setpartner1, setpartner2)
        print(paste(c("Number of unique nodes: ", length(nodes)), collapse=""))
        
        dump = cbind(highest_neg_corrl, lowest_neg_corrl, lowest_pos_corrl, highest_pos_corrl, res_pos, res_neg, res_pos/res_neg, nrow(inNet), length(nodes), nrow(inNet)/length(nodes)) 
        colnames(dump) = c("Lowest neg corr", "Highest neg corr", "Lowest pos corr", "Highest pos corr", "Pos corr edges", "Neg corr edges", "Ratio of Pos to Neg corr edges" , "Total number of edges", "Number of unique nodes", "Ratio of edges to nodes")
        print(dump)
        
        
        
        if(!noPUC){
			df_a = inNet[,c("partner1InFold","partner1_FoldChange")]
			colnames(df_a) = c("partnerInFold","partner_FoldChange")
			rownames(df_a) <- NULL
		
			df_b = inNet[,c("partner2InFold","partner2_FoldChange")]
			colnames(df_b) = c("partnerInFold","partner_FoldChange")
			rownames(df_b) <- NULL
		
			df_ab = rbind(df_a, df_b)
			uniq_df_ab = unique(df_ab)
			print("Unique nodes (partner, fc) df:")
			print(dim(uniq_df_ab))
		
        
        	upreg = length(uniq_df_ab[as.numeric(uniq_df_ab[,"partner_FoldChange"]) > 1, "partner_FoldChange"])
        	dnreg = length(uniq_df_ab[as.numeric(uniq_df_ab[,"partner_FoldChange"]) < 1, "partner_FoldChange"])
        
        	potential_pos_corr_edges = choose(upreg, 2) + choose(dnreg, 2)
        	potential_neg_corr_edges = upreg * dnreg
		}
		
        #edgesDistinctNodes = inNet[which(inNet[,"partner1"]!=inNet[,"partner2"]), ] 
        #print(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse=""))

        
        write.csv(dump, paste0(argv$file, outputFile), row.names=F)
} 


calc_stats(data, noPUC, correlThreshold=0)