args = commandArgs(trailingOnly=TRUE)

poolnetfile = args[1]
metagrfile = args[2]
consis_metagrfile = args[3]
exptmwfile = args[4]
outputFile = args[5]  #"test"
exptanalys = args[6]
metagranalys = args[7]

# cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/pool/brb_deg_lipid_pheno_hfhs_sp/net2group-spfgf/

# poolnetfile:
#poolnetfile = "~/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/pool/brb_deg_lipid_pheno_hfhs_sp/net2group-spfgf/build_netw_FolChMedian_FCAnalys_1_CorrAnalys_1___ileum4wk AbxHFHS_indiv-pval_0.01_comb-pval_0.01_comb-fdr_0.05_.csv"


# edges to be kept from the metagr 2 hfhs file (consistent and dals)
#consis_metagrfile = "/nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/merged/corr/brb_deg_lipid_pheno_sp/p1_cp1_cfdr1/per_analysis/diet-spf-gf-degs-Analys1-consis.csv-comb-pval-output.csv.consis_genes.csv"

#
#metagrfile="/nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/merged/corr/brb_deg_lipid_pheno_sp/p1_cp1_cfdr1/merged_rq_sp_corr_FolChMedian_merged-parallel-output.csv"
 

# expt file:
#exptmwfile = "/nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/pool/expts_mw/mw_comp-output.csv" 

# output file
#outputFile = "test"

#exptanalys = "Analys 1 "
#metagranalys = "Analys 1 "

# usage:
#cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/pool/brb_deg_lipid_pheno_hfhs_sp/net2group-spfgf

#Rscript ~/Morgun_Lab/richrr/scripts/R/parse_pool_metagr_nets_onegr_2datasets.R ~/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/pool/brb_deg_lipid_pheno_hfhs_sp/net2group-spfgf/build_netw_FolChMedian_FCAnalys_1_CorrAnalys_1___ileum4wk\ AbxHFHS_indiv-pval_0.01_comb-pval_0.01_comb-fdr_0.05_.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/merged/corr/brb_deg_lipid_pheno_sp/p1_cp1_cfdr1/merged_rq_sp_corr_FolChMedian_merged-parallel-output.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/merged/corr/brb_deg_lipid_pheno_sp/p1_cp1_cfdr1/per_analysis/diet-spf-gf-degs-Analys1-consis.csv-comb-pval-output.csv.consis_genes.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/pool/expts_mw/mw_comp-output.csv lipid_pheno_net-ip0.01_cp0.01_fdr0.05 Analys\ 1\  Analys\ 1\ 


if(FALSE){



# need to compare edges from pool net
# with the 4 group networks 
# use the expt wise analysis and extract the fc and pval info for each partner of a pair

# to check whether the improved fdr and puc is due to any bias of pooling different expts, we will check whether
#the two genes with correlations are:
#- both increased in ncd of two expts; both increased in hfhs of two expts; positive correlated in pooled analysis
#- both decreased in ncd of two expts; both decreased in hfhs of two expts; positive correlated in pooled analysis
#- 1 incre vs other decr in ncd of two expts; similarly 1 incre vs other decr in hfhs of two expts; negative correlated in pooled analysis


# remove the FP and pairs that are not same dir as of dataset correlations
}

print(args)

poolnet = read.csv(poolnetfile, header=T, check.names=F, sep=',')
#head(poolnet)


# to check whether the improved fdr and puc is due to any bias of pooling different expts, we will check whether
#the two genes with correlations are:
#- both increased in ncd of two expts; both increased in hfhs of two expts; positive correlated in pooled analysis
#- both decreased in ncd of two expts; both decreased in hfhs of two expts; positive correlated in pooled analysis
#- 1 incre vs other decr in ncd of two expts; similarly 1 incre vs other decr in hfhs of two expts; negative correlated in pooled analysis


expt = read.csv(exptmwfile, header=T, check.names=F, sep=',')
rownames(expt) = expt[,"geneName"]
#head(expt)

# keep only the required cols from the particular analysis (group)
ExptCols = grep(exptanalys, colnames(expt), value=T)
exptfc = grep("FolChMedian", ExptCols, value=T)
exptpv = grep("pvalue", ExptCols, value=T)


selexpt = expt[, c("geneName", exptfc, exptpv)] 
#head(selexpt)


p1 = selexpt[ as.vector(poolnet$partner1) , ]
colnames(p1) = paste0("p1_" , colnames(p1))
#head(p1)

p2 = selexpt[ as.vector(poolnet$partner2) , ]
colnames(p2) = paste0("p2_" , colnames(p2))
#head(p2)


# https://stackoverflow.com/questions/46173023/combine-two-data-frames-of-the-same-size-one-column-after-each-other?rq=1
p = cbind(p1, p2)[order(c(seq_along(p1), seq_along(p2)))]
rownames(p) = NULL
#head(p)

write.table(p, paste0(outputFile, "-exptfcpv-tmp-out.txt"), quote=F, row.names=F, sep='\t')



#https://stackoverflow.com/questions/39582031/r-multiple-conditions-in-if-statement
#if (any(i <- (tester$V3> 200 & tester$V4>250))) {tester$V5[i] <- "one"} else {tester$V5[i] <-NA}

# fc dir for genes betwn expts for hfhs
p$fc_dir_of_genes_betn_expts <- 
          ifelse (p[,3] > 0 & p[,4] > 0 ,"up",
			ifelse (p[,3] < 0 & p[,4] < 0 , "down",
			  ifelse (p[,3] > 0 & p[,4] < 0 , "opposite",
			    ifelse (p[,3] < 0 & p[,4] > 0 , "opposite",
			       	"inconsistent"		
			))))


# significant pvals for genes betwn expts for particular group
p$sig_pval_of_genes_betn_expts <- 
            ifelse (p[,5] < 0.05 & p[,6] < 0.05 , 1, 0)


poolnetExpt = cbind(poolnet, p)

poolnetExpt$EdgeType =  
  ifelse( poolnetExpt$fc_dir_of_genes_betn_expts == "up" & poolnetExpt$sig_pval_of_genes_betn_expts == 1 & poolnetExpt$combinedCoefficient > 0, "FP",
ifelse( poolnetExpt$fc_dir_of_genes_betn_expts == "down" & poolnetExpt$sig_pval_of_genes_betn_expts == 1 & poolnetExpt$combinedCoefficient > 0, "FP",
ifelse( poolnetExpt$fc_dir_of_genes_betn_expts == "opposite" & poolnetExpt$sig_pval_of_genes_betn_expts == 1 & poolnetExpt$combinedCoefficient < 0, "FP",
"ok"
)))


head(poolnetExpt)
write.table(poolnetExpt, paste0(outputFile, "-edgetype-out.txt"), quote=F, row.names=F, sep='\t')




# has only the correlation with consistent dir in different groups.(2 hfhs)
metagr = read.csv(metagrfile, header=T, check.names=F, sep=',')
#head(metagr)
#dim(metagr)


# keep only the consis edges
consis_metagr = as.vector(unlist(read.csv(consis_metagrfile, header=F)))
#head(consis_metagr)
#length(consis_metagr)

metagr = metagr[metagr$pairName_Expt_1 %in% consis_metagr, ]
#head(metagr)
#dim(metagr)


# keep only the required cols from the particular analysis (group)
# since at this point the rows are already consistent between different expts
metagrCols = grep(metagranalys, colnames(metagr), value=T)
selcorrcol = grep("Coefficient_Expt_1", metagrCols, value=T)
selmetagr = metagr[, c("pairName_Expt_1", selcorrcol), drop=F]
#head(selmetagr)


# keep edges present in metagr net:
df = merge(poolnetExpt, selmetagr, by=1)
head(df)

dim(poolnetExpt)
dim(selmetagr)
dim(df)

if(nrow(df) == 0){
 print("Quitting since there are no edges common in the pooled and meta 4 gr nets")
 q()
}

# at this point df has common pairs between pool and metagr net 
# these may not be consistent between pool and metagr net.
# since these are already consistent in two expts, we can use only the first expt to
# check consistency between hfhs and ncd.

rownames(df) = df[,1]
#head(df)

selcols = grep("Coefficient", colnames(df), value=T)
selcols = grep("Direction", selcols, value=T, invert = T) # exclude the direction column from pooled net
df_selcols = df[, selcols]
#head(df_selcols)


# either 1, -1, 0 or a combination thereof in unique. length of 1 indicates consistency in groups
consistency_check = function(x){
    out = length(unique(sign(x)))
    out
}


pool_metagr_consis_corr_dir = apply(df_selcols, 1, consistency_check)
df = cbind(df, pool_metagr_consis_corr_dir)

#print("for Consistency check")
#head(df)

# keep pool edges that are consistent dir with metagr and not false positive
df_consis = df[df$pool_metagr_consis_corr_dir == 1 & df$EdgeType == "ok", ]
dim(df)
dim(df_consis)

head(df_consis)


write.table(df_consis, paste0(outputFile, "-metagrconsis_nonfp-out.txt"), quote=F, row.names=F, sep='\t')
#write.csv(rownames(df_consis), paste0(outputFile, ".txt"), quote=F, row.names=F)


# keep only the rows from the df_consis net:
rownames(poolnet) = poolnet[,1]
inNet = poolnet[rownames(df_consis), ]
head(inNet)
dim(df_consis)
dim(inNet)


#----------------------------------------------------------------------
# calc stats of the generated network
#----------------------------------------------------------------------
calc_stats = function(inNet, correlThreshold=0){

        out = file(paste(outputFile,'-stats-log.txt',sep='') ,'a')
        
        combcoeffcol = as.vector(inNet[,"combinedCoefficient"])
        # not sure why it behaves this way!
        lowest_neg_corrl = min(combcoeffcol[combcoeffcol < 0])
        highest_neg_corrl = max(combcoeffcol[combcoeffcol < 0])

        lowest_pos_corrl = min(combcoeffcol[combcoeffcol > 0])
        highest_pos_corrl = max(combcoeffcol[combcoeffcol > 0])

        
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
        
        #print(paste(c("Total number of edges: ", nrow(inNet)), collapse="")) 

        setpartner1 = unique(as.vector(inNet[,"partner1"]))
        setpartner2 = unique(as.vector(inNet[,"partner2"]))

        nodes = union(setpartner1, setpartner2)
        #print(paste(c("Number of unique nodes: ", length(nodes)), collapse=""))
        
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

        #edgesDistinctNodes = inNet[which(inNet[,"partner1"]!=inNet[,"partner2"]), ] 
        #print(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse=""))

        write.csv(paste(c("Lowest pos corr:" , lowest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Highest pos corr:" , highest_pos_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Lowest neg corr:" , lowest_neg_corrl), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Highest neg corr:" , highest_neg_corrl), collapse='') , file=out, row.names=FALSE)

        write.csv(paste(c("Pos corr edges:" , toString(res_pos)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Neg corr edges:" , toString(res_neg)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Ratio of Pos to Neg corr edges:" , toString(res_pos/res_neg)), collapse='') , file=out, row.names=FALSE)
        
        write.csv(paste(c("Potential Pos corr edges:" , toString(potential_pos_corr_edges)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Potential Neg corr edges:" , toString(potential_neg_corr_edges)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Potential Ratio of Pos to Neg corr edges:" , toString(potential_pos_corr_edges/potential_neg_corr_edges)), collapse='') , file=out, row.names=FALSE)
        
        
        write.csv(paste(c("Total number of edges: ", nrow(inNet)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of unique nodes: ", length(nodes)), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of upreg nodes: ", upreg), collapse='') , file=out, row.names=FALSE)
        write.csv(paste(c("Number of dwnreg nodes: ", dnreg), collapse='') , file=out, row.names=FALSE)
        
        write.csv(paste(c("Ratio of edges to nodes: ", nrow(inNet)/length(nodes)), collapse='') , file=out, row.names=FALSE)
        #write.csv(paste(c("Number of unique edges (without self loops): ", nrow(edgesDistinctNodes)), collapse='') , file=out, row.names=FALSE)
        
        close(out)

}

calc_stats(inNet)

print("Done")
