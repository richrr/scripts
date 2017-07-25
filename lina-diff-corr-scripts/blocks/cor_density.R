
# table1 = read.table( "/media/Storage/imsys/cor_dircor_table_byblocks_dados1.txt", header = T)
table1 = read.table( "/media/Storage/imsys/cor0_dados1.txt", header = T)
cor0.group1 = table1[,1]
plot(density(cor0.group1))

table1 = read.table( "/media/Storage/imsys/cor1_dados1.txt", header = T)
cor1.group1 = table1[,1]
plot(density(cor1.group1))

rm(table1)

# table2 = read.table( "/media/Storage/imsys/cor_dircor_table_byblocks_dados2.txt", header = T)
table2 = read.table( "/media/Storage/imsys/cor0_dados2.txt", header = T)

cor0.group2 = table2[,1]
plot(density(cor0.group2))

cor1.group2 = table2[,4]
plot(density(cor0.group2))

rm(table2)

# table3 = read.table( "/media/Storage/imsys/cor_dircor_table_byblocks_dados3.txt", header = T)
table3 = read.table( "/media/Storage/imsys/cor0_dados3.txt", header = T)
cor0.group3 = table3[,1]
plot(density(cor0.group3))


table3 = read.table( "/media/Storage/imsys/cor1_dados3.txt", header = T)
cor1.group3 = table3[,1]
plot(density(cor1.group3))

