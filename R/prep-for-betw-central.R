# make graphml file from dataframe with vertex attributes
library(igraph)
args = commandArgs(trailingOnly=TRUE)

edge_file = args[1]
attrib_file = args[2] # two column file containing "name" and "Class" (otu or gene)

outstr = args[3]

network <- read.csv(edge_file, header = TRUE, stringsAsFactors = FALSE)
networkgraph <- graph.data.frame(network, directed = FALSE)

attributesgraph <- read.csv(attrib_file, header = TRUE, stringsAsFactors = FALSE)

genenode <- attributesgraph$name[attributesgraph$Class == "gene"]
otunode <- attributesgraph$name[attributesgraph$Class == "otu"]


networkgraph <- set.vertex.attribute(networkgraph, "node_type", 
                     index = V(networkgraph)[V(networkgraph)$name %in% genenode],
                     value = c("degs"))

networkgraph <- set.vertex.attribute(networkgraph, "node_type", 
                                     index = V(networkgraph)[V(networkgraph)$name %in% otunode],
                                     value = c("microbe"))

out_file = paste(outstr, "graphml", sep = ".")
write.graph(networkgraph, out_file,format = "graphml")
