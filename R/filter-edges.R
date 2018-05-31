library(argparser)
library(stringr)

# usage: 
#cd ~/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-liver-ileum/rel_quantile/merged/pool/il_li_ph_hfhsncd/p1_cp1_cfdr1/per_analysis
# Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/filter-edges.R --file Analys1-consis.csv-comb-pval-output.csv.consis_genes.csv --keeppartners nodes_to_keep.txt --splitfile _I _L

#cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/merged/pool/brb_deg_lipid_pheno_hfhsncd/p1_cp1_cfdr1/per_analysis/
#Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/filter-edges.R --file Analys1-consis.csv-comb-pval-output.csv.consis_genes.csv --keeppartner1 ~/Morgun_Lab/richrr/Type2_Diabetes/phenotypes/liverlipids/brb-2expt-diet-lipid-fdr15.txt --output filtered_edges.csv --nosplit

#cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-lipids/rel_quantile/merged/pool/brb_deg_lipid_lipid_hfhsncd/p1_cp1_cfdr1/per_analysis/
#Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/filter-edges.R --file Analys1-consis.csv-comb-pval-output.csv.consis_genes.csv --keeppartners ~/Morgun_Lab/richrr/Type2_Diabetes/phenotypes/liverlipids/brb-2expt-diet-lipid-fdr15-union-1gf-diet-lipid-fdr15.txt --output diet-spf-gf-degs-Analys1-consis.csv-comb-pval-output.csv.consis_genes.csv --nosplit



#==================================================================================================================
# parameters
#==================================================================================================================
p <- arg_parser('Keep pairs if they contain certain nodes. The output of this code can be used in --consistent to create networks.')
p <- add_argument(p, "--file", help="consistent edges", nargs=1) # required; but written as optional format so I explicitly mention the arg
# the following keeps edges if both partner1 and partner2 have allowed nodes. should work for most of the cases.
p <- add_argument(p, "--keeppartners", help="list of nodes allowed in partners", nargs=1)
# use the following args if you specifically want to keep edges based on nodes from either partner1 or partner2
p <- add_argument(p, "--keeppartner1", help="list of nodes allowed in partner1", nargs=1)
p <- add_argument(p, "--keeppartner2", help="list of nodes allowed in partner2", nargs=1)
# be careful when using --keeppartners with --keeppartner1 and/or --keeppartner2. 

p <- add_argument(p, "--splitfile", help="string in partner1 on which the allowed pairs (say from --file) will be split into multiple files", nargs=Inf) # this specific case can be used if gene_I-pheno and gene_L-pheno edges were present in the same file
# and now you want to create separate networks of only gene_I-pheno and gene_L-pheno edges.
# so give --splitfile _I _L
p <- add_argument(p, "--nosplit", help="no need to split file.", flag=TRUE)  # default allows split

p <- add_argument(p, "--output", help="output file", default="./filt_netw.csv")

#p <- add_argument(p, "--indivPvalCutoff", help="individualPvalueCutoff", default=0.005, type="numeric") # 0.05
#p <- add_argument(p, "--combPvalCutoff", help="combinedPvalueCutoff", default=0.1, type="numeric") # 0.05
#p <- add_argument(p, "--combFDRCutoff", help="combinedFDRCutoff", default=0.3, type="numeric") # 
#p <- add_argument(p, "--foldchthresh", help="fold change threshold", default=0, type="numeric") # default is logged data so check whether greater than 0. use 1 for unlog data

#p <- add_argument(p, "--foldchMean", help="use fold change from mean", flag=TRUE)  # default fold change is median


argv <- parse_args(p)
#print (p)


if(is.na(argv$file) || (is.na(argv$keeppartners) && is.na(argv$keeppartner1) && is.na(argv$keeppartner2)) )
{
  print (p)
  print("**At least 2 files are required. Give --file file --keeppartner(s/1/2) keeppartner... in cmd**")
  quit()
}


outputFile = argv$output


read_file_as_vector = function(infile){
	indata = read.csv(infile, header=FALSE)
	indata = as.vector(indata[!duplicated(indata[,1]),1] )	#delete duplicated array probes
	indata = indata[order(indata)] # ascending sort
	return(indata)
}


good_nodes = c()
if(!is.na(argv$keeppartners)){
    good_nodes = read_file_as_vector(argv$keeppartners)
    print(length(good_nodes))
}

good_nodes1 = c()
if(!is.na(argv$keeppartner1)){
    good_nodes1 = read_file_as_vector(argv$keeppartner1)
    print(length(good_nodes1))
}

good_nodes2 = c()
if(!is.na(argv$keeppartner2)){
    good_nodes2 = read_file_as_vector(argv$keeppartner2)
    print(length(good_nodes2))
}


consist_elems = read_file_as_vector(argv$file)
consist_elems = grep("<==>", consist_elems, value=TRUE) # keep pairs
   
pair = str_split( consist_elems ,"<==>") 
pairs = t(as.data.frame(pair))
colnames(pairs) = c("p1", "p2")
#head(pairs)
print(nrow(pairs))

# keep pairs where both partners present in vector
if(length(good_nodes) > 0){
    pairs = pairs[ pairs[,"p1"] %in% good_nodes & pairs[,"p2"] %in% good_nodes, ]
    print(nrow(pairs))
}


# keep pairs where partner1 present in vector
if(length(good_nodes1) > 0){
    pairs = pairs[ pairs[,"p1"] %in% good_nodes1 , ]
    print(nrow(pairs))
}


# keep pairs where partner2 present in vector
if(length(good_nodes2) > 0){
    pairs = pairs[ pairs[,"p2"] %in% good_nodes2 , ]
    print(nrow(pairs))
}

good_pairs = paste(pairs[,"p1"], pairs[,"p2"], sep='<==>')
print(length(good_pairs))
#head(good_pairs)
#write.csv(good_pairs, outputFile, quote=F, row.names = F)
write(good_pairs, outputFile, sep="\n")


if(argv$nosplit){
    print("Exiting since split not requested")
    q()
}

if(length(argv$splitfile) > 0){
    
    for(patt in argv$splitfile){
        print(patt)
        tmp = pairs[grep( patt , pairs[,"p1"] ), ] 
        tmp_pairs = paste(tmp[,"p1"], tmp[,"p2"], sep='<==>')
        print(length(tmp_pairs))
        write.csv(tmp_pairs, paste0(outputFile,patt,".csv"), quote=F, row.names = F)
    }

}






