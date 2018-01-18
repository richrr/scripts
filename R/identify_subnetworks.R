library("ProNet")

# /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/DNAShotgunMetagenome/analysis/t2-corr/merged/seed_seed_sp_otu_seed_corr_0.6/p1_cp0.05_cfdr0.1/per_analysis
# [Linux@waterman per_analysis]$ Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/identify_subnetworks.R Analys1-consis.csv-comb-pval-output-combcoeff.csv 2


#' Identify top subnetwork
#'
#' This function returns the induced subgraph of the nodes in the top cluster using mcode
#'
#'
#' @param network in which to identify subnetworks
#' @return A data frame (network) that is Puc compatiable
#' @export


args = commandArgs(trailingOnly=TRUE)

outNetwork = read.csv(args[1], header=T, check.names=F, row.names=1)
#head(outNetwork)

numb = as.numeric(args[2])

  # identify subnetworks
  dfNetwork = outNetwork[,c("partner1", "partner2", "combinedCoefficient")]
  #print(head(dfNetwork))

  #?graph_from_data_frame
  # http://kateto.net/networks-r-igraph
  g = graph_from_data_frame(dfNetwork,directed = F, vertices = NULL)
  #print(g, e=TRUE, v=TRUE)
  # see the edges and vertices
  #E(g) ; V(g) ; edge_attr(g) ; vertex_attr(g)

  # plots the networks
  #cluster(g,method="MCODE",layout="fruchterman.reingold")
  #cluster(g,method="FN",layout="fruchterman.reingold")

  #https://cran.r-project.org/web/packages/ProNet/vignettes/Tutorial.pdf
  # identify clusters
  result <- mcode(g,vwp=0.5,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE)

  
print(summary(result$COMPLEX))
clustrs = nrow(summary(result$COMPLEX))

if(numb > clustrs){
	print(paste("There are only", clustrs, "clusters", collapse=''))
	numb = clustrs
}


counter = 1
while(counter <= numb){
print(counter)
# plot the cluster
cluster1<-induced.subgraph(g,result$COMPLEX[[counter]])
outfile = paste0("mcode_cluster", counter, "_edges.txt")
write_graph(cluster1, outfile, "ncol")
pdf(paste0(outfile,".pdf"))
visualization(cluster1,node.size=4,node.label=V(cluster1)$name,node.label.color="blue")
dev.off()
counter = counter + 1
}


