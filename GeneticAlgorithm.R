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

# function for mutating individual tours via insertion
insertionMutation <- function(tour){
  i <- sample(length(tour), 2)
  node_to_move <- tour[i[1]]
  mutated <- tour[-i[1]]
  mutated <- append(mutated, node_to_move, after = i[2])
  return(mutated)
}

# function for mutating individual tours via inversion
simpleInversionMutation <- function(tour){
  i <- sample(length(tour), 2)
  mutated <- tour
  mutated[i[1]:i[2]] <- rev(tour[i[1]:i[2]])
  return(mutated)
}

runEvolution <- function(nodes, distMat, mutate_func, population = 2000, 
                         mutation_rate = 0.2, n_iterations = 2000, 
                         state = 123, plot = T, n_checkpoints = 100){
  
  if(plot){
    checkpoints <- unique(round(seq(1, population, length.out = n_checkpoints)))
  }
  
  set.seed(state)
  n_nodes <- nrow(nodes)
  # generate random initial population
  tours <- t(replicate(population, sample(n_nodes)))
  # keep sampling until no duplicates are present
  while(any(duplicated(tours))){
    tours[duplicated(tours),] <- replicate(sum(duplicated(tours)), sample(n_nodes))
  }
  
  for(i in 1:n_iterations){
    if(plot & i %in% checkpoints){
      distances <- apply(tours, 1, getTourDistance, d = distMat)
      best_tour <- tours[which.min(distances),]
      plotTour(nodes, best_tour)
    }
    tours <- selectTours(tours, distMat)
    crossover_order <- sample(nrow(tours))
    # sections to swap in each pair of tours
    swap_sections <- t(replicate(floor(nrow(tours) / 2), sample(ncol(tours), 2)))
    offspring <- runCrossover(tours, swap_sections, crossover_order)
    tours <- rbind(tours, offspring)
    # add insertion and inversion mutations with specified frequency
    mutate <- runif(nrow(tours)) < mutation_rate
    tours[mutate,] <- t(apply(tours[mutate,], 1, mutate_func))
  }
  
  distances <- apply(tours, 1, getTourDistance, d = distMat)
  best_tour <- tours[which.min(distances),]
  results <- list('best_tour' = best_tour, 
                  'distance' = distances[which.min(distances)])
  
  return(results)
  
}

# data on 52 locations in berlin
b52 <- readNodeData('data/berlin52.tsp', 'data/berlin52.opt.tour')
b52_distanceMatrix <- as.matrix(dist(b52$nodes))

# compare insertion mutation with simple inversion mutation
b52_results_ins <- runEvolution(b52$nodes, b52_distanceMatrix, 
                                mutate_func = insertionMutation)
b52_results_inv <- runEvolution(b52$nodes, b52_distanceMatrix, 
                                mutate_func = simpleInversionMutation)
b52_best_result <- which.min(c(b52_results_ins$best_distance, 
                               b52_results_inv$best_distance))
b52_results <- c(b52_results_ins, b52_results_inv)[b52_best_result]

# compare to optimal tour
b52_comparison <- compareOptimal(berlin52$nodes, b52_distanceMatrix, 
                                 b52_results$best_tour, berlin52$opt_tour)

# data on 76 cities in Germany
pr76 <- readNodeData('data/pr76.tsp', 'data/pr76.opt.tour')
pr76_distanceMatrix <- as.matrix(dist(pr76$nodes))

# compare insertion mutation with simple inversion mutation
pr76_results_ins <- runEvolution(pr76$nodes, pr76_distanceMatrix, 
                                 mutate_func = insertionMutation)
pr76_results_inv <- runEvolution(pr76$nodes, pr76_distanceMatrix, 
                                 mutate_func = simpleInversionMutation)
pr76_best_result <- which.min(c(pr76_results_ins$best_distance, 
                                pr76_results_inv$best_distance))
pr76_results <- c(pr76_results_ins, pr76_results_inv)[pr76_best_result]

# compare to optimal tour
pr76_comparison <- compareOptimal(pr76$nodes, pr76_distanceMatrix, 
                                  pr76_results$best_tour, pr76$opt_tour)

# data on 48 state capitals in US
att48 <- readNodeData('data/att48.tsp', 'data/att48.opt.tour')
att48_distanceMatrix <- as.matrix(dist(att48$nodes))

# compare insertion mutation with simple inversion mutation
att48_results_ins <- runEvolution(att48$nodes, att48_distanceMatrix, 
                                  mutate_func = insertionMutation)
att48_results_inv <- runEvolution(att48$nodes, att48_distanceMatrix, 
                                  mutate_func = simpleInversionMutation)
att48_best_result <- which.min(c(att48_results_ins$best_distance, 
                                 att48_results_inv$best_distance))
att48_results <- c(att48_results_ins, att48_results_inv)[att48_best_result]

# compare to optimal tour
att48_comparison <- compareOptimal(att48$nodes, att48_distanceMatrix, 
                                   att48_results$best_tour, att48$opt_tour)

