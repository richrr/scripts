# arguments from terminal
args = commandArgs(trailingOnly = FALSE)

table.file = args[3]
annotation.file = args[4]

table = read.table(table.file, header = T, stringsAsFactor = F, sep = "\t")
annotation = read.table(annotation.file, header = T, stringsAsFactor = F, sep = "\t")
annotation$probe.id = paste("X", annotation$probe.id, sep = "")

table = merge(table, annotation, by.x = "gene.x", by.y = "probe.id", all.x=T)
table = merge(table, annotation, by.x = "gene.y", by.y = "probe.id", all.x=T)

table$pair = ifelse(table$Gene.symbol.x < table$Gene.symbol.y, paste(table$Gene.symbol.x , table$Gene.symbol.y), paste(table$Gene.symbol.y , table$Gene.symbol.x))

write.table(table, paste(table.file, "_withpair.txt", sep = ""), sep = "\t", row.names = F)
