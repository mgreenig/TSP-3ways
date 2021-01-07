library(Rcpp)

# function for selecting candidate nodes to swap
selectCandidates <- function(n_nodes, s1, s2){
  # continuity correction for beta distribution, to sample a distance between points
  distance <- ceiling(n_nodes * rbeta(1, s1, s2))
  # number of possible pairs of nodes with that distance
  possible_pairs <- n_nodes - distance
  # sample from the possible nodes
  candidate1 <- sample(possible_pairs, 1)
  candidate2 <- candidate1 + distance
  return(c(candidate1, candidate2))
}

# function for reading in a TSP data set
readNodeData <- function(filepath){
  
  lines <- readLines(filepath)
  table_lines <- lines[7:(length(lines)-2)]
  entries <- strsplit(table_lines, ' ') 
  entries <- lapply(entries, function(entry) entry[2:3])
  table <- do.call(rbind, entries)
  colnames(table) <- c('c1', 'c2')
  table <- as.matrix(apply(table, 2, as.numeric))
  
  return(table)
}


cppFunction('double getTourDistance(NumericMatrix d, NumericVector tour){
  tour.push_back(tour[0]);
  int n = tour.size() - 1;
  double tourDistance = 0;
  int currentNode;
  int nextNode;
  for(int i = 0; i < n; i++){
    currentNode = tour[i] - 1;
    nextNode = tour[i+1] - 1;;
    tourDistance += d(currentNode, nextNode);
  }
  return(tourDistance);
}')

plotTour <- function(nodes, tour){
  
  plot(nodes[,1], nodes[,2], cex = 0.5)
  tour <- c(tour, tour[1])
  
  for(i in 1:(length(tour)-1)){
    lines(c(nodes[tour[i],1], nodes[tour[i+1],1]),
          c(nodes[tour[i],2], nodes[tour[i+1],2]))
  }
}