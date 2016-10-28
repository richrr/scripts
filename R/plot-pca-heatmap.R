#install.packages('ggfortify', repos='http://cran.rstudio.com/')
library(ggfortify)
library(gplots)

plot_pca = function(df, idata, outfile, color){


pdf(paste("pca-", outfile, ".pdf", sep=''))
print(autoplot(prcomp(na.omit(df)), data = idata, colour = color)) # if you use this outside a loop or function, you do not need to say print()
dev.off()

pdf(paste("pca", outfile, "-w-loading.pdf", sep='-'))
print(autoplot(prcomp(na.omit(df)), data = idata, colour = color,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3))
dev.off()

}


#https://learningomics.wordpress.com/tag/heatmap-2/
color.map <- function(group) { 

if (group=="amp") "red" 
else if(group=="cef") "cyan"
 else if(group=="nor") "blue"
 else if(group=="neo") "green"
 else if(group=="met") "yellow"
 else if(group=="vanc") "purple"

}


plot_heatmap = function(m_t1, scaleV, data, outfile, colorbar){

abxcolors <- unlist(lapply(as.vector(data[,colorbar]), color.map))
pdf(paste("htmap-", outfile, ".pdf", sep=''))
print(heatmap.2(m_t1,trace="none",col=greenred(30),scale=scaleV,margins=c(12,8),srtCol=45, RowSideColors = abxcolors))
dev.off()

}


args = commandArgs(trailingOnly=TRUE)

infile = args[1]
outfile = args[2]

data = read.delim(infile, header=T, row.names=1, check.names=F)
dim(data)
#sub_data = data[,-c(18,19)]
sub_data = subset(data, select = -c(Abx,geneName_Expt)) #variable names would NOT be specified in quotes when using subset() function

# na.omit removes rows with NA
# transpose to it removes columns with NA and then transpose it back to original form
sub_data = t(na.omit(t(sub_data)))
#print(dim(sub_data))

plot_pca(sub_data, data, outfile, 'Abx')

m_t1 = as.matrix(sub_data)
plot_heatmap(m_t1, "row", data, outfile, 'Abx')



