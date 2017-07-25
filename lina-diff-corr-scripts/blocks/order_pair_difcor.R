# SGE_Batch -c 'R < scripts/blocks/order_pair_difcor.R --no-save "/capecchi/pharmacy/morgunlab/lina/CervicalCancer/Project_03mp_alias/Merged_difcor_table_Biewenga_Igen_smp.txt" "/capecchi/pharmacy/morgunlab/lina/CervicalCancer/Project_03mp_alias/Merged_difcor_table_Biewenga_Igen_smp2.txt"' -r order_pair_difcor -f 1G -m 1G -q samwise

# arguments from terminal

args = commandArgs(trailingOnly = FALSE)

args

data.file = args[3]
out.file = args[4]

data = read.table(data.file, header=T, sep="", stringsAsFactors=F)

separated_pair = strsplit(data$pair, split = " ")

separated_pair = matrix(unlist(separated_pair), ncol = 2,byrow = T)
colnames(separated_pair) = c("Gene.symbol.x","Gene.symbol.y")
data[, c("Gene.symbol.x", "Gene.symbol.y")] = separated_pair[, c("Gene.symbol.x", "Gene.symbol.y")]


data = data[ order( data$Gene.symbol.x, data$Gene.symbol.y, abs(data$difcor_Biewenga),abs(data$difcor_Igen), decreasing = T),]

data = data[!duplicated(data$pair),]

write.table(data[,-grep("Gene.symbol", colnames(data))], out.file, row.names=F, sep="\t")
