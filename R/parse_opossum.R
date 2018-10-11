# cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/t2-rna-seq-cyto/pooled-net/brb-deg/enrichment/TF/Opossum


args = commandArgs(trailingOnly=TRUE)


infile = args[1]
clean = args[2]

input = read.table(infile, header=T, check.names=T, row.names=1, sep='\t')
#head(input)
dim(input)

# filtered input data as per suggestions in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3429929/
if(clean == "yes"){
   input = input[ which(input$IC >= 9 & input$IC <= 19 & input$GC.Content >= 0.33 & input$GC.Content <= 0.66), ]
   print(dim(input))
   infile = paste0(infile,".filtered")
}

fish.cutoff = median(input$Fisher.score) + (sd(input$Fisher.score))
z.cutoff = median(input$Z.score) + (2*sd(input$Z.score))
input$fpass_cuts = factor( ifelse((input$Fisher.score > fish.cutoff & input$Z.score > z.cutoff), TRUE, FALSE) )
ids = rownames(input)

write.table(input, paste0(infile, ".opossum.parsed.out.txt"), sep='\t',row.names=T, quote=F)


pdf(paste0(infile,".plot.pdf"))
#pairs(~GC.Content+Z.score+Fisher.score,data = input,  main = "Scatterplot Matrix")


library(ggplot2)


source("http://peterhaschke.com/Code/multiplot.R")  # https://stackoverflow.com/questions/18024810/why-am-i-getting-error-could-not-find-function-multiplot-in-rshiny-app?rq=1

p1 = ggplot(input, aes(x=GC.Content, y=Z.score)) + geom_point()

p2 = ggplot(input, aes(x=GC.Content, y=Fisher.score)) + geom_point()

#p3 = ggplot(input, aes(x=Z.score, y=Fisher.score, color = fpass_cuts)) + geom_point()  + geom_hline(yintercept=fish.cutoff, linetype="dashed", color = "blue")  + geom_vline(xintercept=z.cutoff, linetype="dashed", color = "blue")

# https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
# https://stackoverflow.com/questions/8750871/ggplot2-reverse-order-of-scale-brewer/30224727
#boolColors <- as.character(c("TRUE"="#FF0000", "FALSE"="#000000"))
#p3 = ggplot(input, aes(x=Z.score, y=Fisher.score)) + geom_point(aes(colour = fpass_cuts)) + scale_colour_manual(name="passed", values=rev(boolColors))  + geom_hline(yintercept=fish.cutoff, linetype="dashed", color = "blue")  + geom_vline(xintercept=z.cutoff, linetype="dashed", color = "blue")

#http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
# https://stackoverflow.com/questions/15624656/label-points-in-geom-point
#p3 = ggplot(input, aes(x=Z.score, y=Fisher.score, colour = fpass_cuts, label= ids)) + geom_point() + geom_text(aes(label=ifelse((Fisher.score > fish.cutoff & Z.score > z.cutoff), as.character(ids),'')),hjust=-0.1,vjust=1,  size=2) + scale_colour_manual(values=c("#000000", "#FF0000")) + geom_hline(yintercept=fish.cutoff, linetype="dashed", color = "blue")  + geom_vline(xintercept=z.cutoff, linetype="dashed", color = "blue") + theme(legend.position = "none")
#https://stackoverflow.com/questions/17241182/how-to-make-geom-text-plot-within-the-canvass-bounds?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
p3 = ggplot(input, aes(x=Z.score, y=Fisher.score, colour = fpass_cuts, label= ids)) + geom_point() + geom_text(aes(label=ifelse((Fisher.score > fish.cutoff & Z.score > z.cutoff), as.character(ids),'')),hjust="inward",vjust="inward",  size=2) + scale_colour_manual(values=c("#000000", "#FF0000")) + geom_hline(yintercept=fish.cutoff, linetype="dashed", color = "blue")  + geom_vline(xintercept=z.cutoff, linetype="dashed", color = "blue") + theme(legend.position = "none")


multiplot(p1, p2, p3, cols=1)

dev.off()
