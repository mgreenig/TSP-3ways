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
runElasticNet <- function(nodes, n_iterations = 50, alpha = 1, beta = 2, plot = T){
  
  # save mean and sd of columns
  colmeans <- apply(nodes, 2, mean)
  colsds <- apply(nodes, 2, sd)
  
  # scale node coordinates to prevent numerical errors
  nodes <- scale(nodes)
  
  # initialise path as a ellipse
  ring <- initialPath(nrow(nodes) * 2, 
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

nodes <- readNodeData('data/berlin52.tsp')

best_tour <- runElasticNet(nodes, plot = F)
node_distances <- as.matrix(dist(best_tour))
best_tour_distance <- round(getTourDistance(node_distances, 1:nrow(best_tour)), 2)

print(paste('Final tour distance:', best_tour_distance))

