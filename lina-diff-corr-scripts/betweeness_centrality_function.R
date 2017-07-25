betweeness.centrality = function(myNetwork, group1, group2, selection)
{
  
  allPairs = as.matrix(expand.grid(group1, group2))
  allPairs = as.data.frame(allPairs, stringsAsFactors=F)
  if(length(which(allPairs$Var1 == allPairs$Var2))>0)
    allPairs = allPairs[-which(allPairs$Var1 == allPairs$Var2),]
    
  sumAllFractionsForAllNodes = matrix(0,ncol=length(selection)
                                      ,nrow=nrow(allPairs)
                                      ,dimnames= list(
                                        paste(as.vector(t(allPairs[1])),as.vector(t(allPairs[2])),sep="_")
                                        , as.character(selection)
                                      )
  )
  
  
  getShortestPaths = function(pair){
    
    allShortestPaths = get.all.shortest.paths(myNetwork, from= pair[1], to= pair[2], mode = "all", weights=NULL)
    
    getFractionForAllNodes(allShortestPaths,pair)
    
  }
  findExistanceInShortestPath = function (nodeID,allShortestPaths){
    nodesInPath = as.vector(allShortestPaths[,c(-1,-ncol(allShortestPaths))])
    count = table(as.factor(nodesInPath))[as.character(nodeID)]
    if (is.na(count)){
      count=0
    }
    fractionForThisNodeInThisPair =count/nrow(allShortestPaths)
    sumAllFractionsForAllNodes[as.character(nodeID)] <<- sumAllFractionsForAllNodes[as.character(nodeID)] + fractionForThisNodeInThisPair
    
  }
  getFractionForAllNodes = function(allShortestPaths,pair){
    if (length(allShortestPaths$res)==0){
      fractions = cbind(selection,rep(0,times=length(selection)))
      
    }else{
      
      allShortestPaths = do.call(rbind,allShortestPaths$res)
      
      nodesInPath = as.vector(allShortestPaths[,c(-1,-ncol(allShortestPaths))])
      count = table(as.factor(nodesInPath))/nrow(allShortestPaths)
      
      sumAllFractionsForAllNodes[paste(pair,collapse="_"),as.character(intersect(names(count),selection))] <<- count[as.character(intersect(names(count),selection))]
      
    }
    
  }
  
  t = apply(allPairs,1,getShortestPaths)
  sumFractions = colSums(sumAllFractionsForAllNodes)
  
  return(list(fraction.by.pair = sumAllFractionsForAllNodes, sum.fraction.by.node = sumFractions))
}  