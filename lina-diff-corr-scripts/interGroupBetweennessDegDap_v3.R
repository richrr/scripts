
## SGE_Batch -c 'R < ~/scripts/interGroupBetweennessDegDap_v3.R --no-save "/capecchi/pharmacy/morgunlab/lina/imsys/Network/Correlation/Merged_correlation_network_state0.txt" "/capecchi/pharmacy/morgunlab/lina/imsys/Network/Correlation/Supertable_first1000daps.txt" "/capecchi/pharmacy/morgunlab/lina/imsys/Network/Correlation/Betweeness.overall.state0.txt" /capecchi/pharmacy/morgunlab/lina/imsys/Network/Correlation/Norm_betweeness_state0.txt' -r cor_difcor_mine -f 128G -m 128G -q samwise

#############################

args = commandArgs(trailingOnly = FALSE)

network_table = args[3]
nodeIdMappingFile = args[4]
out_table = args[5]
out_table2 = args[6]

############################
#functions

source("scripts/betweeness_centrality_function.R")


pairs.from.vectors = function(group1, group2)
{
  
  allPairs = as.matrix(expand.grid(group1, group2))
  allPairs = as.data.frame(allPairs, stringsAsFactors=F)
  if(length(which(allPairs$Var1 == allPairs$Var2))>0)
    allPairs = allPairs[-which(allPairs$Var1 == allPairs$Var2),]
  
  Pairs = paste(as.vector(t(allPairs[1])),as.vector(t(allPairs[2])),sep="_")  
  
  return(Pairs)
}



###########
#main

idNetworkData = read.table( network_table, stringsAsFactors=F, header=T, sep = "\t")[,c("Gene.symbol.x", "Gene.symbol.y")]
idNetworkData=idNetworkData[!duplicated(paste(idNetworkData$Gene.symbol.x, idNetworkData$Gene.symbol.y)),]

nodeInNetwork = c(idNetworkData$Gene.symbol.x, idNetworkData$Gene.symbol.y)
nodeInNetwork = nodeInNetwork[!duplicated(nodeInNetwork)]

nodeIdMappingData = read.table(nodeIdMappingFile,header=TRUE,stringsAsFactors=FALSE,sep="\t")[,1:4]

nodeIdMappingData = nodeIdMappingData[ nodeIdMappingData$Gene.symbol %in% nodeInNetwork,]

nodeIdMappingData$node.type.deg = ifelse(nodeIdMappingData$node.type.deg == "nDEG",0,1)
nodeIdMappingData$node.type.dap = ifelse(nodeIdMappingData$node.type.dap == "not in DAPs",0,1)
nodeIdMappingData$node.type = ifelse(apply(nodeIdMappingData[3:4],1,sum) == 2, "both", ifelse(nodeIdMappingData$node.type.deg == 1, "DEG", "inDAP"))

genesymboltable = as.data.frame(table(nodeIdMappingData$Gene.symbol))
genesymboldupl = genesymboltable[ genesymboltable[,"Freq"] >1, ]
genesymboldupl = merge(genesymboldupl, nodeIdMappingData[,c("Gene.symbol","node.type")], by.x = "Var1", by.y = "Gene.symbol")
genesymboldupl = genesymboldupl[!duplicated(paste( genesymboldupl$Var1, genesymboldupl$node.type)),]
genesymboldupl = genesymboldupl[duplicated(genesymboldupl$Var1),]

nodeIdMappingData$node.type = ifelse( nodeIdMappingData$Gene.symbol %in% genesymboldupl$Var1, "both", nodeIdMappingData$node.type ) 
nodeIdMappingData = nodeIdMappingData[!duplicated(nodeIdMappingData$Gene.symbol),]
nodeIdMappingData$id.for.igraph = c(1:nrow(nodeIdMappingData))

# View(nodeIdMappingData[, c("Gene.symbol", "id.for.igraph")])
# setwd("/media/Storage/imsys/Project_03mp/Correlation/Betweeness/")
# write.table(nodeIdMappingData[, c("Gene.symbol", "id.for.igraph")], "id.for.igraph.state1.txt", row.names = F, sep = "\t")

idNetworkData = merge(idNetworkData,nodeIdMappingData[,c("Gene.symbol","id.for.igraph")], by.x = "Gene.symbol.x", by.y = "Gene.symbol")
idNetworkData = merge(idNetworkData,nodeIdMappingData[,c("Gene.symbol","id.for.igraph")], by.x = "Gene.symbol.y", by.y = "Gene.symbol")

library(igraph)
myNetwork = graph.data.frame(idNetworkData[,3:4],directed = FALSE)

# DegNodes = nodeIdMappingData[nodeIdMappingData[,"node.type"] %in% c("DEG","both"),"id.for.igraph"]
# DapNodes = nodeIdMappingData[nodeIdMappingData[,"node.type"] %in% c("inDAP","both"),"id.for.igraph"]
allNodes = nodeIdMappingData$id.for.igraph

# OnlyDapNodes = nodeIdMappingData[nodeIdMappingData[,"node.type"] %in% c("inDAP"), "id.for.igraph"]
# OnlyDegNodes = nodeIdMappingData[nodeIdMappingData[,"node.type"] %in% c("DEG"), "id.for.igraph"]

# betweeness.list = betweeness.centrality(myNetwork, DegNodes, DapNodes, selection = nodeIdMappingData$id.for.igraph)
betweeness.list = betweeness.centrality(myNetwork, allNodes, allNodes, selection = allNodes)

sumFractions = colSums(betweeness.list[[1]])

norm.betweeness = sumFractions/nrow(betweeness.list[[1]])

write.table(norm.betweeness, out_table2 , sep = "\t")


write.table(betweeness.list[[1]], out_table , sep = "\t")



# Degs_subtable = pairs.by.nodes(DegNodes, OnlyDapNodes, DegNodes)
# write.table(Degs_subtable, "all_Degs_by_only_Daps_table_for_betweeness_state0.txt")
# 
# Daps_subtable = pairs.by.nodes(OnlyDegNodes, DapNodes, DapNodes)
# write.table(Daps_subtable, "all_Daps_by_only_Degs_table_for_betweeness_state0.txt")
# 
# Degs_subtable = pairs.by.nodes(OnlyDegNodes, DapNodes, DegNodes)
# write.table(Degs_subtable, "only_Degs_by_all_Daps_table_for_betweeness_state0.txt")
# 
# Daps_subtable = pairs.by.nodes(DegNodes, OnlyDapNodes, DapNodes)
# write.table(Daps_subtable, "only_Daps_by_all_Degs_table_for_betweeness_state0.txt")
