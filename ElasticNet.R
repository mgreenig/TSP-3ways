library(Rcpp)
source('TSPfunctions.R')
sourceCpp('Updates.cpp')

# function to get the initial path 
initialPath <- function(n, x_scale, y_scale, x_bias, y_bias){
  matrix(c(x_scale * cos(2 * (1:n) * pi / n) + x_bias,
           y_scale * sin(2 * (1:n) * pi / n) + y_bias),
           n, 2)
}

# plot a path overlayed over the nodes
plotPath <- function(nodes, path){
  path <- rbind(path, path[1,])
  plot(nodes[,1], nodes[,2], cex = 0.5)
  lines(path[,1], path[,2], col = 'red', cex = 0.5)
}

# get a value of K for the given iteration
getK <- function(i, scale = 10){
  1 / exp(i / scale)
}

# run the elastic net for a set of nodes
runElasticNet <- function(nodes, Mfactor = 1.5, n_iterations = 50, alpha = 1, beta = 1, plot = T){
  
  # save mean and sd of columns
  colmeans <- apply(nodes, 2, mean)
  colsds <- apply(nodes, 2, sd)
  
  # scale node coordinates to prevent numerical errors
  nodes <- scale(nodes)
  
  # initialise path as a ellipse
  ring <- initialPath(ceiling(nrow(nodes) * Mfactor), 
                      sd(nodes[,1]) / 2, 
                      sd(nodes[,2]) / 2, 
                      mean(nodes[,1]), 
                      mean(nodes[,2]))
  
  for(i in 1:n_iterations){
    if(plot){
      plotPath(nodes, ring)
      Sys.sleep(1)
    }
    weights <- weightUpdate(nodes, ring, getK(i))
    delta_y <- yUpdate(nodes, ring, weights, getK(i), alpha, beta)
    ring <- ring + delta_y
  }
  
  ring <- t(apply(ring, 1, function(row) row * colsds))
  ring <- t(apply(ring, 1, function(row) row + colmeans))
  
  return(ring)
  
}

# data on 52 locations in berlin
b52 <- readNodeData('data/berlin52.tsp', 'data/berlin52.opt.tour')
b52_distanceMatrix <- as.matrix(dist(b52$nodes))

# run elastic net
b52_results <- runElasticNet(b52$nodes)

print(paste('Final tour distance:', 
            round(getTourDistance(as.matrix(dist(b52_results)), 1:nrow(b52_results)))))
print(paste('Optimal tour distance:', 
            round(getTourDistance(b52_distanceMatrix, b52$opt_tour))))

# data on 76 cities in Germany
pr76 <- readNodeData('data/pr76.tsp', 'data/pr76.opt.tour')
pr76_distanceMatrix <- as.matrix(dist(pr76$nodes))

# run elastic net
pr76_results <- runElasticNet(pr76$nodes)

print(paste('Final tour distance:', 
            round(getTourDistance(as.matrix(dist(pr76_results)), 1:nrow(pr76_results)))))
print(paste('Optimal tour distance:', 
            round(getTourDistance(pr76_distanceMatrix, pr76$opt_tour))))

# data on 48 state capitals in US
att48 <- readNodeData('data/att48.tsp', 'data/att48.opt.tour')
att48_distanceMatrix <- as.matrix(dist(att48$nodes))

# run elastic net
att48_results <- runElasticNet(att48$nodes)

print(paste('Final tour distance:', 
            round(getTourDistance(as.matrix(dist(att48_results)), 1:nrow(att48_results)))))
print(paste('Optimal tour distance:', 
            round(getTourDistance(att48_distanceMatrix, att48$opt_tour))))

