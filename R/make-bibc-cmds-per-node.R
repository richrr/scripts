library(igraph)
args = commandArgs(trailingOnly=TRUE)

# this file makes a different node file for every node for which you want to run bibc
# this is useful when you want to edit the node file of all phenos to run bibc for each pheno

### see if/how this code is different than the make_per_pheno_cmds.py ###

node_file = args[1] # this file contains the node for which you want to make new node files
attrib_file = args[2] # two column file containing "name" and "Class" (otu or gene)

outstr = args[3]
edge_file = args[4] # full path

nodes <- as.vector(unlist(read.csv(node_file, header = FALSE, stringsAsFactors = FALSE)))
print(nodes)

attributesgraph <- read.csv(attrib_file, header = TRUE, stringsAsFactors = FALSE)
head(attributesgraph)
dim(attributesgraph)

#genenode <- attributesgraph$name[attributesgraph$Class == "gene"]
#head(genenode)
#otunode <- attributesgraph$name[attributesgraph$Class == "otu"]
#head(otunode)

cmds = c()

for( node in nodes){
    print(node)
    tmp = attributesgraph
    sanitize_node = gsub("[[:space:]]",'-',node)
    sanitize_node = gsub("[[:punct:]]",'-',sanitize_node)
    
    dir.create(sanitize_node)
    
    #outfile = paste0(sanitize_node, "-", outstr, "-node.csv")
    outfile = paste0(sanitize_node, "/", outstr, "-node.csv")
    
    for (r in 1:nrow(tmp)){
        # change all items from nodes file except the node in iteration to "ignore"
        if(tmp[r, "name"] %in% nodes){
            if(tmp[r, "name"] != node){
            	tmp[r, "Class"] = "ignore"
        	}
        }   
    }
    
    write.csv(tmp, outfile, row.names=F, quote=F)
    print(paste0("printed to ", outfile))


    cmds = c(cmds, paste0("cd ", sanitize_node))  
    
    cmds = c(cmds, paste("Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/prep-for-betw-central.R", args[4], paste0(outstr, "-node.csv"), paste0(outstr,"-",sanitize_node), sep=' '))

    if(args[5] == 'sge'){
		cmds = c(cmds, paste0('SGE_Batch -c "/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/betweennessadapter node_type microbe degs ', paste0(outstr,"-",sanitize_node), '.graphml > bibc-', paste0(outstr,"-",sanitize_node), '-results.txt" -m 50G -F 50G -q transkingdom -M rodrrich@oregonstate.edu -r ', paste0("log-", outstr,"-",sanitize_node)))}
	else{
		cmds = c(cmds, paste0('/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/betweennessadapter node_type microbe degs ', paste0(outstr,"-",sanitize_node), '.graphml > bibc-', paste0(outstr,"-",sanitize_node), '-results.txt'))
	}
	
	cmds = c(cmds, "cd ..")

}

write(cmds, paste0(outstr,"-","bibc-cmds.txt"), sep = "\n")

