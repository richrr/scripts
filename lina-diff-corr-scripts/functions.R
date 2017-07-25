########################################################################################################################################
#function that calculates de coefficient of variation

#cv = function(x)
# { 
#  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
#}

#cv2 = function(x)
#{ 
# sd(x, na.rm=TRUE)/median(x, na.rm=TRUE)
#}

quartis=function(x)
{
  quantile(x,probs=c(0.25,0.5,0.75),na.rm=TRUE)
}

cv3 = function(x)
{ 
  #quart=quantile(x,probs=c(0.25,0.5,0.75),na.rm=TRUE)
  quart=apply(x,2,quartis)
  abs(quart[1,]-quart[3,])/quart[2,]
}

########################################################################################################################################

# Define a function to calculate sample sizes for gene-gene correlations -- courtesy of Vivek Gopalan

# x           is the sequence from 1 to G (i.e. 1:nrow(aa)), where G = number of genes
# aa          is a G x S matrix of gene expression values, where G = number of genes and S = number of samples

# WARNING!!!  The input aa must be a matrix -- 'aa' cannot be a data frame.

min_size_cal         <- function(x,aa){rowSums(!is.na(t(t(aa)*aa[x,])))}

#################################################################################
#########################################################################################################################

##function which orders the pair of genes. The matrix A must have two columns named "Gene.symbol.x" and "Gene.symbol.y"

ordpair.x=function(A){
  A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
  genecol1=unique(c(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]))
  genecol1=sort(genecol1)
  
  for(k in 1:length(genecol1))
  {
    
    if(k == 1)  where = which(A[,"Gene.symbol.y"] == genecol1[k]) else
      where = which(A[,"Gene.symbol.y"] == genecol1[k] & !(A[,"Gene.symbol.x"] %in% genecol1[1:k]))
    
    who = A[where,"Gene.symbol.x"]
    A[where, "Gene.symbol.x"]= c(rep(genecol1[k], length(who)))
    A[where,"Gene.symbol.y"] = who
    
  } 
  A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
  
  return (A)
}
#############################################################################################################################

##function which orders the pair of genes. The matrix A must have two columns named "Gene.symbol.x" and "Gene.symbol.y"
## improvement, if probe.id is true, also switch the probe.ids

ordpair.x_v2=function(A, probe.id){
  A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
  genecol1 = unique(c(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]))
  genecol1 = sort(genecol1)
  
  for(k in 1:length(genecol1))
  {
    
    if(k == 1)  where = which(A[,"Gene.symbol.y"] == genecol1[k]) else
      where = which(A[,"Gene.symbol.y"] == genecol1[k] & !(A[,"Gene.symbol.x"] %in% genecol1[1:k]))
    
    who = A[where,"Gene.symbol.x"]
    A[where, "Gene.symbol.x"]= c(rep(genecol1[k], length(who)))
    A[where,"Gene.symbol.y"] = who
    
    if(probe.id == T)
    {
      who = A[where,"probe.id.x"]
      A[where, "probe.id.x"] = A[where, "probe.id.y"]
      A[where,"probe.id.y"] = who
    }
    
  } 
  A=A[order(A[,"Gene.symbol.x"],A[,"Gene.symbol.y"]),]
  
  return (A)
}

#######################3
median.normalization = function(data)
{
  median = apply(data, 2, median)
  
  for(i in 1:ncol(data))
    data[,i] = data[,i]/median[i]
  
  return(data)
  
}
###########################
# This function replaces all alias names by its official name 
# input data: data frame with "Gene.symbol" column
#       alias.data: data frame with columns "Gene.symbol.official" and "alias"
#       column.name: name or number of the column with gene symbol
replace.alias.by.official.gene.symbol = function(data,alias.data, column.name)
{
  data[,"alias"] = ifelse(data[,column.name] %in% alias.data[,3],1,0)
  data[,"official"] = ifelse(data[,column.name] %in% alias.data[,2],1,0)
  data[,"total"] = data[,"alias"]+data[,"official"]
  # write.table(data[data$total == 2,],"/media/Storage/CervicalCancer/cervicalCancerdataOfficialAndAlias.txt", row.names=F, sep="\t")
  
  #only considered Gene.symbols present in taxID file that are either official or alias
  data = data[data[,"total"]!=0,]
  
  data = merge(data, alias.data[,2:3], by.x=column.name, by.y="alias", all.x=T)
  # removing data that are alias with duplicated official Gene.symbol
  data = data[! (data[,column.name] %in% data[data$alias==1,][which(duplicated(data[data$alias==1,column.name])),column.name]),]
  
  data$Gene.symbol.official = ifelse(data$official==1, data[,column.name], data[,paste(column.name,"official",sep=".")])
  
  return(data)
}

#################################################################################
# this fuction separates de column pair and removes duplicates according to highest difference of correlation
# Requirement
remove.duplicates.pair.by.difcor = function(data, difcor.column)
{
  separated_pair = strsplit(data$pair, split = " ")
  
  separated_pair = matrix(unlist(separated_pair), ncol = 2,byrow = T)
  colnames(separated_pair) = c("Gene.symbol.x","Gene.symbol.y")
  data[, c("Gene.symbol.x", "Gene.symbol.y")] = separated_pair[, c("Gene.symbol.x", "Gene.symbol.y")]
  
  
  data = data[ order( data$Gene.symbol.x, data$Gene.symbol.y, abs(data[, difcor.column]), decreasing = T),]
  
  data = data[!duplicated(data$pair),]
  return(data)
  
}


#################################################################################

# this function returns all genes that are partners of a specific pairs table
# pairs table must have columns Gene.symbol.x and Gene.symbol.y

partners = function(gene, pairs_table)
{
  pairs = pairs_table[(pairs_table$Gene.symbol.x == gene) | (pairs_table$Gene.symbol.y == gene),]
  partners = ifelse(pairs$Gene.symbol.y == gene, pairs$Gene.symbol.x, pairs$Gene.symbol.y)
  partners = paste(partners, collapse = " ")
  return(partners)
}

#################################################################################

# this function returns all genes that are partners of a specific pairs table
# pairs table must have columns Gene.symbol.x and Gene.symbol.y
# input: if selection == TRUE it picks only partners that are in selection
#        if selection == FALSE, it picks all partners                        

partners.selection = function(gene, pairs_table, selection = FALSE, selection.vector = NULL)
{
  pairs = pairs_table[(pairs_table$Gene.symbol.x == gene) | (pairs_table$Gene.symbol.y == gene),]
  partners = ifelse(pairs$Gene.symbol.y == gene, pairs$Gene.symbol.x, pairs$Gene.symbol.y)
  
  if(selection == TRUE)
    partners = partners[which(partners %in% selection.vector)]
  
  degree = length(partners)
  if(degree == 0) partners = NA else
    partners = paste(partners, collapse = " ")
  
  return(list(Freq = degree, Partners = partners))
}


###############################################################################33
#
# This function filters the expression data according to 
# filtro_BRB=TRUE: filtering done in BRB or
# 
# filtro_BRB=FALSE: mp which is the max missing percentage allowed 

filtering_missing=function(dados,filtro_BRB=FALSE,mp=0.3){
  if(filtro_BRB=="yes") filter=read.table("FILTER.txt") else
    filter=as.data.frame(ifelse(colSums(is.na(dados))<nrow(dados)*mp,1,0)) #FALSE=1 e TRUE = 0
  
  #filtrando os dados 
  #filtering data 
  dados=dados[,which(filter==1)]
  
  return(dados)
  
}

###################################################################################

separate.classes=function(data,file.name,class_type=c(0,1)){
  
  if(!is.list(data)) data=as.list(data)
  
  class1= lapply(data, function(x) which(x[,"class"] %in% class_type[1]))
  class2= lapply(data, function(x) which(x[,"class"] %in% class_type[2]))
  
  dados_list=list()
  for(j in 1:length(data))
    dados_list[[substr(file.name[j],1,nchar(study.file[j])-4)]]=list(class1 = as.matrix(data[[j]][class1[[j]],-ncol(data[[j]])]), 
                                                                     class2 = as.matrix(data[[j]][class2[[j]],-ncol(data[[j]])]))
  return(dados_list)
}
#############################################################################################
separate.classes2=function(data,study.names,class_type=c(0,1)){
  
  if(!is.list(data)) data=as.list(data)
  
  class1= lapply(data, function(x) which(x[,"class"] %in% class_type[1]))
  class2= lapply(data, function(x) which(x[,"class"] %in% class_type[2]))
  
  dados_list=list()
  for(j in 1:length(data))
    dados_list[[study.names[j]]]=list(class1 = as.matrix(data[[j]][class1[[j]],-ncol(data[[j]])]), 
                                      class2 = as.matrix(data[[j]][class2[[j]],-ncol(data[[j]])]))
  return(dados_list)
}

##########################################################################################
# this funtion returns the different of median : class1 - class2
# data: a list with case and control separated. it can be a list of lists
#

dif.median.classes = function(data){
  
  if(!is.list(data)) print("data must be a list with class and control data separate") else
   if (length(data[[1]]) == 1) data=list(data)
  
  dif_med=list()
  
  med = lapply(data, function(x) lapply(x, function(x) apply(x, 2, median, na.rm = T)))
  dif_degs=lapply(med, function(x) t(x[[1]]-x[[2]]))
  
  return(dif_degs)
}

##########################################################################################
# this funtion returns the different of median : class2 - class1
# input data: a list with case and control separated. it can be a list of lists
#

dif.median.classes.v2 = function(data){
  
  if(!is.list(data)) print("data must be a list with class and control data separate") else
    if (length(data[[1]]) == 1) data=list(data)
  
  dif_med=list()
  
  med = lapply(data, function(x) lapply(x, function(x) apply(x, 2, median, na.rm = T)))
  dif_degs=lapply(med, function(x) x[[2]]-x[[1]])
  
  dif_degs = lapply(dif_degs, as.matrix)
  
  dif_degs = as.data.frame(dif_degs)
  
  colnames(dif_degs) = paste("dif_median",names(data), sep = ".")
  
  dif_degs[, "probe.id"] = substr(row.names(dif_degs), 2 , nchar(row.names(dif_degs)))
 
  row.names(dif_degs) = NULL
  
  return(dif_degs[,c(length(data)+1,1:length(data))])
}

#################################################################################
# This function returns the genes with the n_f highest coefficient of variance
# in each state

# n_f: number of genes with highest cv in each state

cvfilter=function(dados_list,n_f)
{
  cvlist=lapply(dados_list, cv3)
  
  # sorting both lists separately
  cvlist= lapply(cvlist, sort,decreasing=TRUE)
  
  #pegando os n genes com os n maiores coeficientes de variacao no caso e controle
  #getting n genes with greater variance coefficients in both conditions
  geneidfilter=c(names(cvlist[[1]][1:n_f]),names(cvlist[[2]][1:n_f]))
  
  #removing the duplicates.
  geneidfilter=geneidfilter[!duplicated(geneidfilter)]
  
  return(geneidfilter)
  
}
############################################################################################################################

subsetting=function(studies,p,pi){
  
  if(!is.vector(studies)) return("object is not a vector!!")
  
  if(p>length(studies)) return("p is greater than length of studies!!")
  
  if(pi>p) return("impossible to get a pith part greater than the total number of parts!!")
  
  #divisao inteira
  n=length(studies)%/%p
  # #resto da divisao
  f=(length(studies)%%p)/p
  
  if (pi <= f*p) studies=studies[((pi-1)*n+pi):(pi*(n+1))] else
    studies=studies[((pi-1)*n+1+f*p):(pi*n+f*p)]
  
  return(studies)
  
}

####################################################################################################
# this function selects the degs expression and returns a list with expression of case and control 
# separately
# 
# input: data = file with expression data with one column named "class" containing the class information
# degs = vector cantaining the degs
# case = character of number specifying samples from case studies (default = 0)
# control = character of number specifying samples from control studies (default = 1)
# total_studies = total number of studies (default = 1)

degs.expression =function(data,degs,case=0,control=1,total_studies=1){
  
  data_case=list()
  data_control=list()
  
  
  for(i in 1:total_studies)
  {
    #SUPONDO QUE OS GENES ESTAO NA MESMA ORDEM EM CADA TABELA e j? com os geneids!!!!!!!!!!!!!!!!
    
    #transposing in order to filter
    
    
    data_case[[i]]=data[[i]][data[[i]][,"class"]==control,]
    data_case[[i]]=data_case[[i]][,-ncol(data_case[[i]])]
    
    data_case[[i]]=data_case[[i]][,colnames(data_case[[i]]) %in% degs]
    
    #################################################################################################################
    data_control[[i]]=data[[i]][data[[i]][,"class"]==case,]
    data_control[[i]]=data_control[[i]][,-ncol(data_control[[i]])]
    
    data_control[[i]]=data_control[[i]][,colnames(data_control[[i]]) %in% degs]
    
  }
  
  dl=list(case=data_case,control=data_control)
  
  return(dl)
}
###########################################################################################


###########################################################################################
pairwise_num_obs=function(dados_degs,total_studies=1){
  
  min.mat=list()
  
  for (j in 1:total_studies)
    min.mat[[j]]   <- sapply(1:nrow(t(dados_degs[[j]])),FUN=min_size_cal,as.matrix(t(dados_degs[[j]])))
  
  return(min.mat)
}

#################################################################################################################
#     function: calculate correlation coefficient of a pair of genes
# input:
#         1.pair: a pair of "gene"s:
#                 value: a vector of length 2, containing names of the two "genes"
#         2. strain,method (default: "pearson") and use (default: )
# output correlation coefficient: vector of length 2: "correlation Coefficient","pvalue"

cor.perpair=function(data, pair, method = "pearson", use = "pairwise.complete.obs"){
  
  pair = as.vector(pair)
  
#   print(pair)
  
  c1 = as.numeric( as.vector(data[, grep(pair[1], colnames(data))]))
  c2 = as.numeric( as.vector(data[, grep(pair[2], colnames(data))]))
  
  if (sum(! is.na(c1  +c2)) < 3){
    return (matrix(c(NA, NA, NA)))
  }
  p = cor.test(c1, c2, method = method, use = use)
  outLine = data.frame(c(p$estimate,p$p.value,as.integer(p$parameter)))
  
  if(method == "pearson") rownames(outLine) = c(paste(method,"correlation"),"pvalue","df") else
    rownames(outLine) = c(paste(method,"correlation"),"pvalue") 
  
  #colnames(outLine) = paste(x,collapse=",")
  return(outLine)
  
}


# cor.perpair(pair,data,method="spearman")

#     function: calculate correlation coefficient of a pair of genes
# input:
#         1.pair: a pair of "gene"s:
#                 value: a vector of length 2, containing names of the two "genes"
#         2. strain,method (default: "pearson") and use (default: )
# output correlation coefficient: vector of length 2: "correlation Coefficient","pvalue"

t.test.pergene=function(data){
  
  counter = 1
  outLine = data.frame(mean.class1=NA, mean.class2=NA, dif.mean.pvalue = NA, dif.mean = NA)
  
  if (sum(colnames(data[[1]]) == colnames(data[[2]])) != ncol(data[[1]]))
    stop("colnames must match in both classes")
  
  for(i in 1:ncol(data[[1]]))
  {
#     print(i)
    c1 = data$class1[,i]
    c2 = data$class2[,i]
    if(sum(!is.na(c1)) < 2 | sum(!is.na(c2)) < 2 )
      outLine[counter, 1:4] = NA else
      {
        t = t.test(c1, c2)
        x = data.frame(c(t$estimate, t$p.value))
        outLine[counter, 1:3] = x[, 1]
        outLine[counter, "dif.mean"] = t$estimate[1] - t$estimate[2]
      }
    counter = counter + 1
  }
  
  row.names(outLine) = colnames(data[[1]])
  nrow(outLine)
  return(outLine[,c(1,2,4,3)])
  
}



#################################################################################################################
# this function returns a vector with pvalue for correlation
# input: a vector with correlation and a vector with number of samples

cor.pv = function(cor,num_smp, method = "pearson")
{
  if(method == "pearson")
  {
    # t statistics
    t <- (sqrt(num_smp-2)*cor)/sqrt(1-cor^2)
    #  pvalue
    pvalue <- 2 * pt( - abs(t), num_smp - 2)
  }
  
  return(pvalue)
}

  
#############################################################################################  
# argumentos dessa funcao sao listas - calcular para varias bases de uma vez (as bases estao na listas)
cor.pv.pairwise=function(cor,min.mat,total_studies=1){
  
  pvalue=list()
  t=list()
  
  for (j in 1:total_studies)
  {
    t[[j]] <- (sqrt(min.mat[[j]]-2)*cor[[j]])/sqrt(1-cor[[j]]^2)
    #  pvalue
    pvalue[[j]] <- 2*pt(-abs(t[[j]]), min.mat[[j]]-2)
  }
  
  return(pvalue)
}

# calcula para matrizes e retorna uma matriz simetrica
cor.pv.pairwise1=function(cor,min.mat){
  t <- (sqrt(min.mat-2)*cor)/sqrt(1-cor^2)
  #  pvalue
  pvalue <- 2*pt(-abs(t), min.mat-2)
  
  return(pvalue)
}

cor.pv.pairwise2=function(cor,min.mat){
  
  if(!is.list(cor)) cor = list(cor)
  if(!is.list(min.mat)) min.mat = list(min.mat)
  
  if( length(cor) != length(min.mat) )
  {
    return(print("ERROR: length of both lists differ"))
    stop()
  }
  
  pvalue=list()
    
  for (j in 1:length(cor))
    pvalue[[j]] = cor.pv.pairwise1(cor[[j]],min.mat[[j]])
  
  if(length(names(cor)) != 0) names(pvalue) = names(cor)
   
  return(pvalue)
}

cor.pv.pairwise3=function(cor,min.mat){
  if(!is.list(cor)) cor=list(cor)
  if(!is.list(min.mat)) cor=list(min.mat)
  
  if( length(cor) != length(min.mat) )
  {
    return(print("ERROR: length of both lists differ"))
    stop()
  }
  
  pvalue=list()
    
  for (j in 1:length(cor))
      pvalue[[j]]=cor.pv.pairwise2(cor[[j]],min.mat[[j]])
    
  if( length(names(cor)) != 0 ) names(pvalue) = names(cor)
        
  return(pvalue)
  
}
#############################################################################################
# This function gives a table with values of correlation, pvalue and number or 
# samples used to calculate the pairwise correlation (just for Pearson) and direction of correlation
# 1 : positive correlation
# -1: negative correlation
# 
# Inputs
# It uses function cor.perpair written above which uses cor.test from base package
# data: expression values
# pairs: the pairs of genes which you want to calculate it. Must be organized in a matrix with two
#         columns. the rows represent the pairs and the columns the component of the pairs
# method: method to calculate the correlation coefficient (default = pearson)
# use: how to deal with missing values (default: pairwise.complete.obs) for more info check the cor function from base

table.perpair.cor.pv = function(data, pairs, method = "pearson", use = "pairwise.complete.obs"){
  
  final_table = data.frame(pairs) 
  colnames(final_table) = c("ID.x","ID.y")
  
  for(pair_index in 1:nrow(pairs))
  {
    
    temp = cor.perpair(data, pair = pairs[pair_index ,], method = method, use = use)
    
    for(j in 1: total_studies)
      final_table[pair_index, "cor"] = temp[1, 1]
    for(j in 1:total_studies)
      final_table[pair_index, "pv"] = temp[2, 1]
    
    if(method == "pearson")
      for(j in 1: total_studies)
        final_table[pair_index, "pairwise_num_smp"] = temp[3,1] + 2
  }
  
  for(j in 1: total_studies)
    final_table[, "dircor"] = ifelse( final_table[, "cor"] > 0, 1, -1)   
  
  return(final_table)
  
}
###########################################################################################



join.pv=function(dados_finais,total_studies=1)
{
  pvalue=matrix(ncol=total_studies,nrow=nrow(dados_finais))
  for(j in 1:total_studies)
    pvalue[,j]=dados_finais[,paste("pv_gr",j,sep="")]
  return(pvalue)
} 

#####################33

Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p),na.rm = T)
  p.val <- 1-pchisq(Xsq, df = 2*length(!is.na(p)))
  return(c(Xsq = Xsq, p.value = p.val))
}

#################################################################################################
# this function returns fisher pvalue for a meta analysis
# ONLY SIZE 3 
# input: pvalue: a matrix where each column represents the pvalue of a study. 
# The number of columns must be that same as the number of studies.

pv.meta.analysis=function(pvalue){
  
#   print(class(pvalue))
  
  if(!is.list(pvalue))
  {
    temp=pvalue
    pvalue=list()
    
    print(class(temp))
    if(is.matrix(temp))
      for(i in 1:ncol(temp)){
        pvalue[[i]]=temp[,i]
        print(class(pvalue[[i]]))
  }}
  
  
  t=-2*(log(pvalue[[1]])+log(pvalue[[2]])+log(pvalue[[3]]))
  
  pv_fish=1-pchisq(t,6)
    
  return(pv_fish)
}

################################################################################################
# this function returns fisher pvalue for a meta analysis
# ANY SIZE 
# input: pvalue: a matrix where each column represents the pvalue of a study. 
# The number of columns must be that same as the number of studies.


pv.meta.analysis.v2 = function(pvalue){
  
  if(!is.list(pvalue))
  {
    temp=pvalue
    pvalue=list()
    
    print(class(temp))
    if(is.matrix(temp))
      for(i in 1:ncol(temp)){
        pvalue[[i]]=temp[,i]
        print(class(pvalue[[i]]))
      }}
  
  sum = log(pvalue[[1]])
  for(j in 2 : length(pvalue))
    sum = sum + log(pvalue[[j]])
  
  t=-2*(sum)
  
  pv_fish=1-pchisq(t, 2*length(pvalue))
  
  print(paste("total studies ",length(pvalue)))
  
  return(pv_fish)
}
##################################################################################

pv.meta.analysis.v3=function(pvalue){
  
  #   print(class(pvalue))
  
  if(!is.list(pvalue))
  {
    temp=pvalue
    pvalue=list()
    
    print(class(temp))
    if(is.matrix(temp))
      for(i in 1:ncol(temp)){
        pvalue[[i]]=temp[,i]
        print(class(pvalue[[i]]))
      }
  } 
        
 
  logpvalue = log(pvalue)
	
  sum = apply(logpvalue, 1,sum, na.rm = T)

  sum.not.na = function(vector) {
    sum.not.na = sum(!is.na(vector))
    return(sum.not.na)
  }
  
  num.datasets = apply(pvalue, 1, sum.not.na)
 
  t=-2*(sum)

  pv_fish=1-pchisq(t, 2*num.datasets)
  
  print(paste("total studies ",length(pvalue)))
  
  return(pv_fish)
}

########################################
pv.meta.analysis.v4=function(pvalue){
  
  
  logpvalue = log(pvalue)
  
  sum = apply(logpvalue, 1,sum, na.rm = T)
  
  sum.not.na = function(vector) {
    sum.not.na = sum(!is.na(vector))
    return(sum.not.na)
  }
  
  num.datasets = apply(pvalue, 1, sum.not.na)
  
  t=-2*(sum)
  
  pv_fish=1-pchisq(t, 2*num.datasets)
  
  print(paste("total studies ",length(pvalue)))
  
  return(pv_fish)
}


##############################################################################################
# date: 12/13/13
# this function filters a table with only rows that have pvalues higher than a certain threshold
# inputs
# data = table with your data
# total_studies = total number of studies (default = 1)
# pv_thrs = threshold (default = 0.2)
# pv_gr = name patern in your data table fot the pvalue you want to filter (default = "pv_gr")
# ID.x = name patern of first gene id (default = "ID.x")
# ID.y = name patern of second gene id (default = "ID.y")
# meta = boolean variable indicating if you want to run pvalue metanalysis on the filtered pvalues (default = FALSE)
# fdr = boolean variable indicating if you want to run fdr on the filtered pvalues (default = FALSE)

pv.filter = function(data, total_studies = 1, pv_thrs = 0.2, pv_gr = "pv_gr", ID.x = "ID.x", ID.y = "ID.y", meta = FALSE, fdr = FALSE){
  
  pvalue=list()
  
  for(j in 1:total_studies)
  {      
    if(j == 1) 
    {
      pv_filter = data.frame(ifelse(data[, paste(pv_gr, j, sep = "")] < pv_thrs,1,0), row.names = paste(data[, ID.x], data[, ID.y]), stringsAsFactors = F) 
      colnames(pv_filter)[j] = paste(pv_gr, j, sep= "")
    } else
      pv_filter[, paste(pv_gr, j, sep = "")] = ifelse(data[, paste(pv_gr, j, sep = "")] < pv_thrs, 1, 0)
  }
  
  data = data[apply(pv_filter, 1, sum, na.rm = T) == total_studies,]
  
  for(j in 1:total_studies)
    pvalue[[j]] = data[, paste("pv_gr", j, sep = "")]
  
  if(meta == TRUE)
    data[,"pv_fisher"] = pv.meta.analysis(pvalue)
  
  if(fdr == TRUE)
    data[,"fdr_fisher"] = p.adjust(data[, "pv_fisher"], method = "fdr")
  
  return(data)
}






###########################################################################################
# This function is to be used when we have more than one probe.ID per Gene.symbol and we    # 
# would like to compare by Gene.symbol.                                                     #
#                                                                                           #       
# It selects the Gene.symbol with higher comparison value: correlation difference in DAPs   #
# for example                                                                               #
# INPUT: Table with at least two columns (Gene.symbol and difference column)                #
# OUTPUT: Table with all columns (unique Gene.symbol and higher difference) with row.names  #
###########################################################################################

select.gene.symbol=function(dif_degs, gene.symbol.col=1, comp.value.col=2,row.names=TRUE,rm.na=FALSE){
  
  if(rm.na) dif_degs = dif_degs[!is.na(dif_degs[,comp.value.col]),]
  
  dif_degs = dif_degs[order(dif_degs[,gene.symbol.col],- abs(dif_degs[,comp.value.col]) ), ] #sort by id and abs(value)
  dif_degs = dif_degs[ !duplicated(dif_degs[,gene.symbol.col]), ]   
  
  if(!row.names) row.names(dif_degs)=NULL
  
  return(dif_degs)
}


#####################################################################################################

# function that calculates the pvalue for difference of correlation test from a list

pv.dif.cor = function(data1){
      
    if(!"psych" %in% installed.packages()) install.packages("psych")
    library(psych)
  
    if(is.list(data1))
    {
    # convert correlations to z-scores 
    z.class1            = fisherz(data1$class1[grep("correlation",row.names(data1$class1)),1])
    z.class2            = fisherz(data1$class2[grep("correlation",row.names(data1$class2)),1])
    
    df.class1=data1$class1["df",1]
    df.class2=data1$class2["df",1]
    }
    
    # Calculate vector of t-tests to compare correlations between classes
    fisher          = (z.class1 - z.class2) / sqrt((1/((df.class1+2) - 3)) + (1/((df.class2+2) - 3)))
    
    # Calculate raw p-values
    pv    = 2*pt(-abs(fisher),Inf)
  
  
  return(pv)
}


#####################################################################################################

# function that calculates the pvalue for difference of correlation test from vectors

pv.dif.cor2 = function(cor1,cor2,num_smp1,num_smp2){
  
  if(!"psych" %in% installed.packages()) install.packages("psych")
  library(psych)
  
  # convert correlations to z-scores 
  z.class1            = fisherz(cor1)
  z.class2            = fisherz(cor2)
  
  # Calculate vector of t-tests to compare correlations between classes
  fisher          = (z.class1 - z.class2) / sqrt((1/(num_smp1 - 3)) + (1/(num_smp2 - 3)))
  
  # Calculate raw p-values
  pv.dif.cor    = 2*pt(-abs(fisher),Inf)
  
  pv.dif.cor = ifelse(num_smp1 < 3 | num_smp2 < 3, NA ,pv.dif.cor)
    
  return(pv.dif.cor)
}
#########################################################################################
# this function returns the pvalue of difference of correlation
pv.dif.cor.table = function(final_table, study.names, class_type)
{
  
  for(j in 1:length(study.names))
  {
    final_table[,paste("pv_difcor_",study.names[j],sep="")] = 
      pv.dif.cor2(cor1 = final_table[,paste("cor",class_type[1],"_",study.names[j],sep = "")],
                  cor2 = final_table[,paste("cor",class_type[2],"_",study.names[j],sep = "")],
                  num_smp1= final_table[,paste("num_smp_cor",class_type[1],"_",study.names[j],sep = "")],
                  num_smp2= final_table[,paste("num_smp_cor",class_type[2],"_",study.names[j],sep = "")])
    
  }
    
  return(final_table)
}

##############################################################################

fdr.pv.dif.cor.table = function(final_table, study.names)
{
  
  for(j in 1:length(study.names))
  {
    final_table[,paste("fdr_pv_difcor_",study.names[j],sep="")] = 
      p.adjust(final_table[,paste("pv_difcor_",study.names[j],sep="")],method="fdr")  
  }
  
  return(final_table)
}


##############################################################################
#returns a matrix with pairwise number of samber

num.pairwise.sample = function(data) {
  
  data=sapply(1:nrow(t(data)),FUN=min_size_cal,as.matrix(t(data)))
  colnames(data)=row.names(data)

  return(data)
}
#############################################################################################
# this function filters a table by pairs of genes who have the same direction in difference of correlation for all studies
# inputs: 
#       final_table.list: either one table containing all studies columns named difcor_blablabla or a list of table containing difcor separately
#       study_names: a vector containing the names that difference difcor amongst each other - for example: for "difcor_study1" and "difcor_study2", study.names would be c("study1","study2")
# output: one filtered single table containing all info about all studies in study.names 

filter.dircor = function(final_table, study.names)
{
  total_studies = length(study.names)
  
  for(j in 2:total_studies)
  {
    cond1 = (final_table[,paste("difcor_",study.names[1],sep="")] > 0) & (final_table[,paste("difcor_",study.names[j],sep="")] > 0)
    cond2 = (final_table[,paste("difcor_",study.names[1],sep="")] < 0) & (final_table[,paste("difcor_",study.names[j],sep="")] < 0)
    final_table = final_table[cond1 | cond2,]
  }
  return(final_table)
}
######################################################################################################################
# this function is an extension of the previous one but also considering cases when we have missing is some tables
#       final_table.list: either one table containing all studies columns named difcor_blablabla or a list of table containing difcor separately
#       study_names: a vector containing the names that difference difcor amongst each other - for example: for "difcor_study1" and "difcor_study2", study.names would be c("study1","study2")
# output: one filtered single table containing all info about all studies in study.names 

filter.dircor.v2 = function(final_table, study.names )
{

  total_studies = length(study.names)
  
  difcor = final_table[, c(paste("difcor", study.names, sep = "_"))]
  
  for(j in 1:total_studies)
  {
    if(j ==1)
      dirdifcor = data.frame(ifelse(difcor[, j] > 0, 1, -1)) else
        dirdifcor[, study.names[j]] = ifelse(difcor[, j] > 0, 1, -1)
    
    colnames(dirdifcor)[1] = study.names[1]
  }
  
  final_table[, "number.no.na"] = apply(!is.na(dirdifcor), 1, sum , na.rm = T)
  
  final_table = final_table[abs(apply(dirdifcor, 1, sum , na.rm = T)) ==  final_table[, "number.no.na"], ]
  
  
  return(final_table)
  
}

##########################################################################################
# this function returns a table with columns indicating the direction of correlation 
# it also calculates the correlation pvalue in case there is no columns beginning with "pv_cor"
#
# directions under pv_thrsd
#    1: signifcant positive correlation 
#    0: insignifcant correlation  
#    -1: signifcant negative correlation
#
# WARNINGS: if there is columns with correlation pvalue, name it pv_cor_<study.name>
#           don't start any other column with "pv_cor"\
#
# inputs: 
#

cor.directions = function(final_table, pv_thrsd, total_studies=1, study.names = "gr1",class_type=c(0,1), method = "pearson")
{
  if(is.list(class_type)) class_type.list = class_type
  
  #calculating pvalues in case they haven't been already calculated
  pv_calculated = grep("pv_cor",colnames(final_table))
  if(length(pv_calculated) == 0)
  {
    for( j in 1:length(study.names))
    {
      if(exists("class_type.list")) class_type = class_type.list[[j]]
      
      for(i in 1:length(class_type))
        final_table[,paste("pv_cor",class_type[i],"_",study.names[j],sep = "")] = cor.pv(final_table[,paste("cor",class_type[i],"_",study.names[j],sep = "")], final_table[,paste("num_smp_cor",class_type[i],"_",study.names[j],sep = "")], method = method)
    
    }
  }
    
  for( j in 1:length(study.names))
  {
    if(exists("class_type.list")) class_type = class_type.list[[j]]    
    for(i in 1:length(class_type))
    {
      cond1 = final_table[,paste("pv_cor",class_type[i],"_",study.names[j],sep = "")] > pv_thrsd
      cond2 = final_table[,paste("cor", class_type[i],"_",study.names[j],sep = "")] > 0
      final_table[,paste("dircor",class_type[i],"_",study.names[j],sep = "")] = ifelse(cond1,0,ifelse(cond2,1,-1))
    }  
  }
  return(final_table)
}

#############################################################################################

filter.insign.cor.both.states = function(final_table, study.names = "gr1", class_type = c(0,1))
{
  sum = rep(0, nrow(final_table))
  
  for(j in 1:length(study.names))
  {
    sum1 = sum + final_table[,paste("dircor",class_type[1],"_",study.names[j],sep = "")]
    sum2 = sum + final_table[,paste("dircor",class_type[2],"_",study.names[j],sep = "")]
  }
  
  final_table = final_table[sum1 != 0 | sum2 != 0,]
  
  return(final_table)

}


#############################################################################################
# this function returns a table with only the pairs with coherent correlation directions.
# which means which have the same direction in all studies for all states
# 
# inputs:
# 
# WARNING: class_type names must be the same to all studies
# 

coherent.cor.directions=function(final_table,pv_thrsd = 0.01, total_studies=1, study.names = "gr1" ,class_type=c(0,1), insign.cor = TRUE, method = "pearson"){
   
  #checking if final_table already has cor directions assigned
  dir_assigned = grep("dircor",colnames(final_table))
  
  if(length(dir_assigned) == 0)
    final_table = cor.directions(final_table, pv_thrsd, total_studies, study.names , class_type, method = method)
       
  # only keeping the pairs with the same correlation directions in all stydies for both classes
  
  if(insign.cor == FALSE)
    for(i in 1:length(class_type))
      for(j in 2:total_studies)
      {
        cond1 = final_table[, paste("dircor",class_type[i],"_",study.names[j],sep = "")] == 0
        final_table = final_table[cond1,]
      }
    
  for(i in 1:length(class_type))
    for(j in 2:total_studies)
    {
      cond1 = final_table[, paste("dircor",class_type[i],"_",study.names[1],sep = "")] == final_table[,paste("dircor",class_type[i],"_",study.names[j],sep = "")]
      final_table = final_table[cond1,]
    }
  
  
  return(final_table)
}

################################################################################################

coherent.cor.directions1=function(final_table,pv_thrsd,total_studies){
  
  cor.class1=data.frame(cor1=ifelse(final_table[,paste("pv_cor",class_type[1],"_gr1",sep="")]>pv_thrsd,0,final_table[,paste("cor",class_type[1],"_gr1",sep="")]))
  cor.class2=data.frame(cor1=ifelse(final_table[,paste("pv_cor",class_type[2],"_gr1",sep="")]>pv_thrsd,0,final_table[,paste("cor",class_type[2],"_gr1",sep="")]))
  
  dir.class1=data.frame(dir1=ifelse(cor.class1[,"cor1"]==0,0,final_table[,paste("dir_cor",class_type[1],"_gr1",sep="")]))
  dir.class2=data.frame(dir1=ifelse(cor.class2[,"cor1"]==0,0,final_table[,paste("dir_cor",class_type[2],"_gr1",sep="")]))
  
  if(total_studies>1)
    for(j in 2:total_studies)
    {
      # assigning 0 to the statistically insignificant correlations
      cor.class1[,paste("cor",j,sep="")]=
        ifelse(final_table[,paste("pv_cor",class_type[1],"_gr",j,sep="")]>pv_thrsd,0,final_table[,paste("cor",class_type[1],"_gr",j,sep="")])
      cor.class2[,paste("cor",j,sep="")]=
        ifelse(final_table[,paste("pv_cor",class_type[2],"_gr",j,sep="")]>pv_thrsd,0,final_table[,paste("cor",class_type[2],"_gr",j,sep="")])
      
      # assigning 0 to the direction of statistically insignificant correlations
      dir.class1[,paste("cor",j,sep="")]= ifelse(cor.class1[,paste("cor",j,sep="")]==0,0,final_table[,paste("dir_cor",class_type[1],"_gr",j,sep="")])
      dir.class2[,paste("cor",j,sep="")]= ifelse(cor.class2[,paste("cor",j,sep="")]==0,0,final_table[,paste("dir_cor",class_type[2],"_gr",j,sep="")])
    }
  
  for(j in 1:total_studies)
  {
    final_table[,paste("dir_cor",class_type[1],"_gr",j,sep="")]=dir.class1[,j]
    final_table[,paste("dir_cor",class_type[2],"_gr",j,sep="")]=dir.class2[,j]
  }
  
  # only keeping the pairs with the same correlation directions in both classes
  
  # 1: signifcant positive correlation in all studies -----------------------
  # 0: insignifcant correlation in all studies ------------------------------
  # -1: signifcant negative correlation in all studies ----------------------
  
  dir.cor.filter.class1=apply(apply(dir.class1,1,duplicated),2,sum)+1==total_studies
  dir.cor.filter.class2=apply(apply(dir.class2,1,duplicated),2,sum)+1==total_studies
  
  final_table=final_table[dir.cor.filter.class1 & dir.cor.filter.class2,]
  
  return(final_table)
}
##############################################################################################

# This function returns a vector with the direction coherent in a minimum amount of studies    
# 
# input: cor_table: a matrix with correlation. Each row is a pair and each column is a study.
#       min_studies: the minimum amount of studies required 

coherent.cor.dir.vector = function(cor_table, min_studies = 2)
{
  if(min_studies < ncol(cor_table)/2)
    print("minimum studies required less than halph the number of studies. Does not make sense!")
  
  countPos = rep(0,nrow(cor_table))
  countNeg = rep(0,nrow(cor_table))
  
  for(j in 1 : ncol(cor_table))
  {
    countPos = ifelse(cor_table[, j] > 0, countPos + 1, countPos)
    countNeg = ifelse(cor_table[, j] < 0, countNeg + 1, countNeg)
  }
  
  coherent.cor.dir = ifelse(countPos >= min_studies , 1,
                            ifelse(countNeg >= min_studies, -1, NA))
  return(coherent.cor.dir)
}

#############################################################################################
# this function returns a table with only the pairs with coherent correlation directions in
# at least one class.
# which means which have the same direction in all studies for at least one state
# 
# inputs:
# 
# WARNING: class_type names must be the same to all studies
# 

coherent.cor.directions.one.state=function(final_table,pv_thrsd_cor = 0.2, total_studies=1, study.names = "gr1" ,class_type=c(0,1), method = "pearson"){
  
  #checking if final_table already has cor directions assigned
  dir_assigned = grep("dircor",colnames(final_table))
  
  if(length(dir_assigned) == 0)
    final_table = cor.directions(final_table, pv_thrsd_cor, total_studies, study.names , class_type, method = method)
  
  # only keeping the pairs with the same correlation directions in all stydies for both classes
  
  dircor0_index = grep("dircor0", colnames(final_table))
  sum_dircor0 = apply(final_table[,dircor0_index],1,sum)
  
  dircor1_index = grep("dircor1", colnames(final_table))
  sum_dircor1 = apply(final_table[,dircor1_index],1,sum)
  
  cond = abs(sum_dircor0) == total_studies | abs(sum_dircor1) == total_studies 
  final_table = final_table[cond,]
  
  return(final_table)
}

# ###################################################################3
# this function returns a list, where each component is related to one class and
# has only the genes with coherent directions in all studies for that particular class

coherent.cor.directions.list=function(final_table,pv_thrsd_cor = 0.2, total_studies=1, study.names = "gr1" ,class_type=c(0,1), method = "pearson"){
  #checking if final_table already has cor directions assigned
  dir_assigned = grep("dircor",colnames(final_table))
  
  if(length(dir_assigned) == 0)
    final_table = cor.directions(final_table, pv_thrsd_cor, total_studies, study.names , class_type, method = method)

  final_table.list = list()
  
  for(i in 1:length(class_type))
  {
    dircor_index = grep(paste("dircor", class_type[i], sep = ""), colnames(final_table))
    sum_dircor = apply(final_table[,dircor_index],1,sum)
    cond = abs(sum_dircor) == total_studies 
    
    final_table.list[[i]] = final_table[cond, c(grep("\\.x", colnames(final_table)), grep("\\.y", colnames(final_table)), grep(paste("cor", class_type[i], sep = ""), colnames(final_table)))]
  }
  
  names(final_table.list) = c("class1", "class2")
  return(final_table.list)
  
}
################################################################################################

coherent.dircor.pvalue = function(final_table, pv_thrsd_difcor = 0.2)
{
  pv.dircor.index = grep("pv_difcor", colnames(final_table))
  
  for( i in 1: length(pv.dircor.index))
  {
    cond = final_table[, pv.dircor.index[i]] < pv_thrsd_difcor | is.na(final_table[, pv.dircor.index[i]])
    final_table = final_table[cond,]
  }
  return(final_table)
}

##############################################################################################
coherent.directions=function(final_table,pv_thrsd,total_studies){
  
  final_table = coherent.cor.directions1(final_table,pv_thrsd,total_studies)
    
  if(nrow(final_table)==0) return(final_table) else
  {
    dir_difcor=data.frame(dir_difcor_gr1=final_table[,"dir_difcor_gr1"])
    if(total_studies > 1)
      for(j in 2: total_studies)
        dir_difcor[,paste("dir_difcor_gr",j,sep="")]=final_table[,paste("dir_difcor_gr",j,sep="")]
    
    dir.difcor.filter=apply(apply(dir_difcor,1,duplicated),2,sum)+1==total_studies
    
    final_table = final_table[dir.difcor.filter,]
    
    return(final_table)
    
  }
  
}
################################################################################
# columns must be a vector with the number of columns corresponding to the different studies

coherent.directions.general = function(final_table, columns){
  
  total_studies = length(columns)
  
  dir = matrix(NA, nrow=nrow(final_table), ncol=total_studies)
  for(j in 1:total_studies)
  {
    dir[,j] = ifelse(final_table[,columns[j]] > 0, 1, -1)
  }
  
  sum_dir = apply(dir,1,sum)
  
  final_table = final_table[abs(sum_dir) == total_studies,]
    
    
  return(final_table)
}
##############################################################################################

# this function receives a list of expression of different states 
# it outputs a list of expression of the genes that have the highest coefficient of variation 
# in at least one of the states based on the median

filter.by.cv = function(data_list, max_genes)
{
  if(ncol(data_list[[1]]) != ncol(data_list[[2]]))
    return(print("ERROR: You must have the same number of genes in both states/conditions"))
  
  if(max_genes >= ncol(data_list[[1]]))
    print(paste("WARNING: less genes than specified - considering total number of genes:",ncol(data_list[[1]])))

  if(max_genes < ncol(data_list[[1]]))
  {
    geneidfilter=cvfilter(data_list, max_genes)
    
    #pegando os dados dos genes filtrados para caso na posicao 1 e controle na 2
    #getting the data for the genes in last step for each condition
  
    data_list=lapply(data_list, function(x) x[,geneidfilter])
    
  }
      
  return(data_list)
}


#####################################################################################################
# function that find DAPs
# r a list with two numeric matrices or data frames. Each one related to one state

dapfinder=function(dados_list,n_f)
{
  if(n_f > ncol(dados_list[[1]]) | n_f > ncol(dados_list[[2]]))
    return(print("ERROR: less rows than specified")) else
      
      if(n_f < ncol(dados_list[[1]]) | n_f < ncol(dados_list[[2]]))
      {
        geneidfilter=cvfilter(dados_list,n_f)
        n=length(geneidfilter)
        
        #pegando os dados dos genes filtrados para caso na posicao 1 e controle na 2
		# Taking the data of the filtered genes for case in position 1 and control in the 2
        #getting the data for the genes in last step for each condition
        dadosf=lapply(dados_list, function(x) x[,geneidfilter])
        
        rm(geneidfilter)
        gc()
        
      }else if(n_f == ncol(dados_list[[1]]) & n_f == ncol(dados_list[[2]]))
      {
        dadosf=dados_list
        n=ncol(dados_list[[1]])
        
      }else
        return(print("ERROR: Check the conditions... there might be some mistake..."))
  
  
  remove(dados_list)
  gc()
  
  #separating into 2 matrices
  dados_caso=dadosf[[1]]
  dados_controle=dadosf[[2]]
  
  #cleaning up memory
  remove(dadosf)
  gc()
  
  # tabela de dados com as seguintes colunas: todos os possiveis pares de gens, correlacao no caso,
  # correlacao no controle, pvalor p caso, pvalor p controle, z, F, pvalorF, FDR
  
  #table with the following columns: every possibel pair of genes, correlation in condition 1 (caso), 
  #correlation in condition 2 (controle), pvalue for caso, pvalue for controle, z, F, Fpvalue, FDR
  #
  
  
  dados_finais=as.data.frame(matrix(data=NA, (n*(n-1)/2),1))
  colnames(dados_finais)=c("gene1")
  
  # colnames(dados_finais)=c("gene1", "gene2", "r_caso", "r_controle","r_caso-r_controle", "pv_dif","direction","filtro")
  
  #Matriz de correlacao de pearson de caso e controle separadamente
  #Pearson correlation matrix of caso and controle separately
  cor_caso=cor(dados_caso,use="pairwise.complete.obs")
  cor_controle=cor(dados_controle,use="pairwise.complete.obs")
  
  #organizing the names of the genes. Is there a better way to do this?
  
  genes=colnames(cor_caso)
  
  comb=combn(rev(genes),2)
  
  dados_finais[,"gene1"]=rev(comb[2,])
  dados_finais[,"gene2"]=rev(comb[1,])
  
  # dados_finais[,1]=rev(comb[2,])
  # dados_finais[,2]=rev(comb[1,])
  
  rm(comb,genes)
  gc()
  
  #arranging correlation in table
  dados_finais[,"r_caso"]=cor_caso[upper.tri(cor_caso)]
  dados_finais[,"r_controle"]=cor_controle[upper.tri(cor_caso)]
  
  rm(cor_caso,cor_controle)
  gc()
  
  # test H0: rho = 0 for all pairwise correlations
  
  min.mat_caso          <- sapply(1:nrow(t(dados_caso)),FUN=min_size_cal,as.matrix(t(dados_caso)))
  min.mat_controle      <- sapply(1:nrow(t(dados_controle)),FUN=min_size_cal,as.matrix(t(dados_controle)))
  
  min_caso              <- min.mat_caso[upper.tri(min.mat_caso)]
  rm(min.mat_caso)
  gc()
  
  min_controle          <- min.mat_controle[upper.tri(min.mat_controle)]
  rm(min.mat_controle)
  gc()
  
  #removendo os pares tem menos obs que 50%)
  #removing pairs with more than 50% missing)
  
  missing=c(which(min_caso<nrow(dados_caso)/2 | min_caso<=3),which(min_controle<nrow(dados_controle)/2 | min_controle<=3))
  missing=sort(missing[!duplicated(missing)])
  
  rm(dados_caso,dados_controle)
  gc()
  
  if(length(missing)!=0)
  {
    dados_finais=dados_finais[-missing,]
    min_caso=min_caso[-missing]
    min_controle=min_controle[-missing]
  }
  
  rm(missing)
  gc()
  
  #absolute value of the difference of correlations in both conditions
  dados_finais[,"r_caso-r_controle"]=dados_finais[,"r_caso"]-dados_finais[,"r_controle"]
  
  # convert correlations to z-scores (requires psych package)
  z_caso                 <- fisherz(dados_finais[,"r_caso"])
  z_controle             <- fisherz(dados_finais[,"r_controle"])
  
  # Calculate vector of t-tests to compare correlations between classes
  
  fisher          <- (z_caso - z_controle) / sqrt((1/(min_caso - 3)) + (1/(min_controle - 3)))
  
  rm(min_caso,z_caso)
  rm(min_controle,z_controle)
  gc()
  
  
  # Calculate raw p-values
  
  dados_finais[,"pv_dif"]    <- 2*pt(-abs(fisher),Inf)
  rm(fisher)
  gc()
  
  # FDR for t-tests (requires nFDR package)
  
  dados_finais[,"fdr"]=p.adjust(dados_finais[,"pv_dif"],method="fdr")
  
  #verificar pq fdr esta dando prox de 1 e depois mudar p fdr!!!
  #dados_finais=subset(dados_finais,(pv_dif<0.2)|(fdr<0.2))
  gc()
  
  
  #direcao
  #1:dif cor > 0
  #-1:dif cor < 0
  
  dados_finais[,"direction"]=ifelse(dados_finais[,"r_caso-r_controle"]>0 ,1,-1)
  
  dados_finais[,"filtro"]=paste(dados_finais[,"gene1"],dados_finais[,"gene2"],dados_finais[,"direction"])
  
  return(dados_finais)
}


