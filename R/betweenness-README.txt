betweennessadapter [attribute] [group1] [group2] [input_filepath]

the input file is in the form .graphml that I created using igraph in R, you need a node attribute such as "node_type" to differentiate between different groups

example:
./betweennessadapter node_type degs microbe ./transking_uncoll_gene.graphml

it outputs onto the command line though, perhaps you know a way to place the output into a file using bash


