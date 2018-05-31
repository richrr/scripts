args = commandArgs(trailingOnly=TRUE)

poolnetfile = args[1]
metagrfile = args[2]
exptmwfile = args[3]

outputFile = args[4]  #"test"

# usage:
#cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-limf/rel_quantile/merged/pool/brb_deg_f_f_hfhsncd/p1_cp1_cfdr1/per_analysis/net2group
#Rscript ~/Morgun_Lab/richrr/scripts/R/parse_pool_metagr_nets.R /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-limf/rel_quantile/merged/pool/brb_deg_f_f_hfhsncd/p1_cp1_cfdr1/per_analysis/net2group/build_netw_FolChMedian_FCAnalys_12_CorrAnalys_1___ileum4wk\ Abx_indiv-pval_0.3_comb-pval_0.05_comb-fdr_0.15_.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-limf/rel_quantile/merged/corr/brb_deg_f_f_sp/p1_cp1_cfdr1/merged-Analysis-1-2-consis-corr-dir4-datasets.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-limf/rel_quantile/pool/expts_mw/mw_comp-output.csv f_f_net-ip0.3_cp0.05_fdr0.15

### code tested using ###
#sftp://files.cgrb.oregonstate.edu:732/raid1/home/pharmacy/rodrrich/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-liver-ileum/rel_quantile/merged/pool/brb_deg_il_li_hfhsncd/p1_cp1_cfdr1/per_analysis/net2group/pool-li-il-net.xlsx

#Rscript ~/Morgun_Lab/richrr/scripts/R/parse_pool_metagr_nets.R /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-liver-ileum/rel_quantile/merged/pool/brb_deg_il_li_hfhsncd/p1_cp1_cfdr1/per_analysis/net2group/build_netw_FolChMedian_FCAnalys_12_CorrAnalys_1___ileum4wk\ Abx_indiv-pval_0.3_comb-pval_0.05_comb-fdr_0.15_.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-liver-ileum/rel_quantile/merged/corr/brb_deg_sp/p1_cp1_cfdr1/merged-Analysis-1-2-consis-corr-dir4-datasets.csv /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-liver-ileum/rel_quantile/pool/expts_mw_brb_degs_ph/mw_comp-output.csv l_i_net-ip0.3_cp0.05_fdr0.15
###


if(FALSE){



# need to compare edges from pool net
#/nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-limf/rel_quantile/merged/pool/brb_deg_f_f_hfhsncd/p1_cp1_cfdr1/per_analysis/net2group/build_netw_FolChMedian_FCAnalys_12_CorrAnalys_1___ileum4wk\ Abx_indiv-pval_0.3_comb-pval_0.05_comb-fdr_0.15_.csv
# with the 4 group networks 
#/nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-limf/rel_quantile/merged/corr/brb_deg_f_f_sp/p1_cp1_cfdr1/merged-Analysis-1-2-consis-corr-dir4-datasets.csv


# use the analysis from below and extract the fc and pval info for each partner of a pair
#/nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-limf/rel_quantile/pool/expts_mw/mw_comp-output.csv



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


#Mat[which(Mat[,'A'] == 10), ]


exptfc = grep("FolChMedian", colnames(expt), value=T)
exptpv = grep("pvalue", colnames(expt), value=T)

selexpt = expt[, c("geneName", exptfc, exptpv)] 
#head(selexpt)
#selexpt["ENSMUSG00000003949_L",]



p1 = selexpt[ as.vector(poolnet$partner1) , ]
colnames(p1) = paste0("p1_" , colnames(p1))
#head(p1)

p2 = selexpt[ as.vector(poolnet$partner2) , ]
colnames(p2) = paste0("p2_" , colnames(p2))
#head(p2)

#p = cbind(poolnet, p1, p2)
# https://stackoverflow.com/questions/46173023/combine-two-data-frames-of-the-same-size-one-column-after-each-other?rq=1
p = cbind(p1, p2)[order(c(seq_along(p1), seq_along(p2)))]
rownames(p) = NULL
#head(p)

write.table(p, paste0(outputFile, "-exptfcpv-tmp-out.txt"), quote=F, row.names=F, sep='\t')



#https://stackoverflow.com/questions/39582031/r-multiple-conditions-in-if-statement
#if (any(i <- (tester$V3> 200 & tester$V4>250))) {tester$V5[i] <- "one"} else {tester$V5[i] <-NA}

# fc dir for genes betwn expts for hfhs and ncd
p$fc_dir_of_genes_betn_expts <- 
          ifelse (p[,3] > 0 & p[,4] > 0 & p[,5] > 0 & p[,6] > 0,"up",
			ifelse (p[,3] < 0 & p[,4] < 0 & p[,5] < 0 & p[,6] < 0, "down",
			  ifelse (p[,3] > 0 & p[,4] < 0 & p[,5] > 0 & p[,6] < 0 , "opposite",
			    ifelse (p[,3] < 0 & p[,4] > 0 & p[,5] < 0 & p[,6] > 0, "opposite",
			       	"inconsistent"		
			))))


# significant pvals for genes betwn expts for hfhs and ncd
p$sig_pval_of_genes_betn_expts <- 
            ifelse (p[,7] < 0.05 & p[,8] < 0.05 & p[,9] < 0.05 & p[,10] < 0.05, 1, 0)


poolnetExpt = cbind(poolnet, p)

poolnetExpt$EdgeType =  
  ifelse( poolnetExpt$fc_dir_of_genes_betn_expts == "up" & poolnetExpt$sig_pval_of_genes_betn_expts == 1 & poolnetExpt$combinedCoefficient > 0, "FP",
ifelse( poolnetExpt$fc_dir_of_genes_betn_expts == "down" & poolnetExpt$sig_pval_of_genes_betn_expts == 1 & poolnetExpt$combinedCoefficient > 0, "FP",
ifelse( poolnetExpt$fc_dir_of_genes_betn_expts == "opposite" & poolnetExpt$sig_pval_of_genes_betn_expts == 1 & poolnetExpt$combinedCoefficient < 0, "FP",
"ok"
)))


head(poolnetExpt)

write.table(poolnetExpt, paste0(outputFile, "-edgetype-out.txt"), quote=F, row.names=F, sep='\t')




# has only the correlation with consistent dir in different groups.(2 hfhs and 2 ncd)
metagr = read.csv(metagrfile, header=T, check.names=F, sep=',')
#head(metagr)

selcorrcol = grep(".Coefficient_Expt_.1", colnames(metagr), value=T)
selmetagr = metagr[, c("pairName_Expt_1", selcorrcol), drop=F]
#head(selmetagr)

# keep edges present in metagr net:

df = merge(poolnetExpt, selmetagr, by=1)
#head(df)

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

