library(Rcpp)
source('TSPfunctions.R')
sourceCpp('Updates.cpp')

initialRing <- function(n, x_scale, y_scale, x_bias, y_bias){
  matrix(c(x_scale * cos(2 * (1:n) * pi / n) + x_bias,
           y_scale * sin(2 * (1:n) * pi / n) + y_bias),
           n, 2)
}

plotRing <- function(nodes, ring){
  ring <- rbind(ring, ring[1,])
  plot(nodes[,1], nodes[,2], cex = 0.5)
  lines(ring[,1], ring[,2])
}

getK <- function(i, amp = 2, scale = 10){
  amp * (1 / (1 + exp(i / scale)))
}

runElasticNet <- function(nodes, distMat, n_iterations = 50, alpha = 1, beta = 2){
  
  # save mean and sd of columns
  colmeans <- apply(nodes, 2, mean)
  colsds <- apply(nodes, 2, sd)
  
  # scale node coordinates to prevent numerical overflow
  nodes <- scale(nodes)
  
  # initialise ring
  ring <- initialRing(nrow(nodes) * 2, 
                      sd(nodes[,1]) / 2, 
                      sd(nodes[,2]) / 2, 
                      mean(nodes[,1]), 
                      mean(nodes[,2]))
  
  for(i in 1:n_iterations){
    plotRing(nodes, ring)
    weights <- weightUpdate(nodes, ring, getK(i))
    delta_y <- yUpdate(nodes, ring, weights, getK(i), alpha, beta)
    ring <- ring + delta_y
    Sys.sleep(2)
  }
  
  ring <- apply(ring, 1, function(row) row + colmeans)
  ring <- apply(ring, 1, function(row) row + colsds)
    
  return(ring)
  
}

nodes <- readNodeData('data/berlin52.tsp')

distanceMatrix <- as.matrix(dist(nodes))

best_tour <- runElasticNet(nodes, distanceMatrix)
