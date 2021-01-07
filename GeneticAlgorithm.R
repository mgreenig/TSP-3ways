library(Rcpp)

source('TSPfunctions.R')
sourceCpp('Crossover.cpp')

# selection process for tours
selectTours <- function(tours, distanceMatrix, proportion = 0.5, state = 123){
  
  set.seed(state)
  # get distance of each tour, scale for numerical efficiency
  distances <- apply(tours, 1, getTourDistance, d = distanceMatrix)
  distances <- scale(distances)
  
  # probability of a tour surviving is inversely proportional to its distance
  probs <- exp(-distances)/sum(exp(-distances))
  keep <- sample(nrow(tours), floor(nrow(tours) * proportion), prob = probs)
  new_tours <- tours[keep,]
  return(new_tours)
}

# function for mutating individual tours
insertionMutation <- function(tour){
  i <- sample(length(tour), 2)
  node_to_move <- tour[i[1]]
  mutated <- tour[-i[1]]
  mutated <- append(mutated, node_to_move, after = i[2])
  return(mutated)
}

runEvolution <- function(nodes, distMat, population = 1000, 
                         mutation_rate = 0.2, n_iterations = 5000, 
                         state = 123){
  
  set.seed(state)
  n_nodes <- nrow(nodes)
  # generate random initial population
  tours <- t(replicate(population, sample(n_nodes)))
  # keep sampling until no duplicates are present
  while(any(duplicated(tours))){
    tours[duplicated(tours),] <- replicate(sum(duplicated(tours)), sample(n_nodes))
  }
  
  for(i in 1:n_iterations){
    tours <- selectTours(tours, distMat)
    crossover_order <- sample(nrow(tours))
    # sections to swap in each pair of tours
    swap_sections <- t(replicate(floor(nrow(tours) / 2), sample(ncol(tours), 2)))
    offspring <- runCrossover(tours, swap_sections, crossover_order)
    tours <- rbind(tours, offspring)
    # add mutations with specified frequency
    mutate <- runif(nrow(tours)) < mutation_rate
    tours[mutate,] <- t(apply(tours[mutate,], 1, insertionMutation))
  }
  
  distances <- apply(tours, 1, getTourDistance, d = distMat)
  best_tour <- tours[which.min(distances),]
  results <- list('best_tour' = best_tour, 
                  'distance' = distances[which.min(distances)])
  
  return(results)
  
}

nodes <- readNodeData('data/berlin52.tsp')

distanceMatrix <- as.matrix(dist(nodes))

best_tour <- runEvolution(nodes, distanceMatrix, population = 2000, mutation_rate = 0.2)
