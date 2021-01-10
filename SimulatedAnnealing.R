source('TSPfunctions.R')

# get the temperature for iteration x
getTemp <- function(x, amp, loc, scale){
  amp * (1 / (1 + exp((x - loc) / scale)))
}

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

# get probability of accepting a new tour, given two distances and temperature
getProb <- function(old_dist, new_dist, temp){
  min(1, exp((old_dist - new_dist) / temp))
}


# function for running simulated annealing
runAnnealing <- function(nodes, distMat, n_iterations = 50000, 
                         state = 123, tempAmplitude = 4000, tempScale = 3000,
                         preferNeighbors = T, s1 = 2, s2 = 12,
                         plot = T, n_checkpoints = 100){
  
  if(plot){
    checkpoints <- unique(round(seq(1, n_iterations, length.out = n_checkpoints)))
  }
  
  set.seed(state)
  n_nodes <- nrow(nodes)
  # start with a random tour
  tour <- sample(n_nodes)
  tour_dist <- getTourDistance(distanceMatrix, tour)
  best_tour <- tour
  best_distance <- Inf
  # initialise vectors for tracking probabilities and distances
  probs <- numeric(length = n_iterations)
  distances <- numeric(length = n_iterations)
  # getCandidates function uses either the custom function or random sampling
  if(preferNeighbors){
    getCandidates <- function(n) selectCandidates(n, s1 = s1, s2 = s2)
  } else {
    getCandidates <- function(n) sample(n, size = 2, replace = F)
  }
  # iterate the annealing process
  for(i in 1:n_iterations){
    if(plot & i %in% checkpoints){
      plotTour(nodes, tour)
      Sys.sleep(1)
    }
    temp <- getTemp(i, amp = tempAmplitude, loc = 0, scale = tempScale)
    # get candidates, drawn randomly depending on the current iteration
    swapped <- getCandidates(n_nodes)
    candidate_tour <- tour
    # swap the selected points in the candidate tour
    candidate_tour[swapped[1]:swapped[2]]<- rev(tour[swapped[1]:swapped[2]])
    candidate_dist <- getTourDistance(distanceMatrix, candidate_tour)
    prob <- getProb(tour_dist, candidate_dist, temp)
    probs[i] <- prob
    if(runif(1) < prob){
      tour <- candidate_tour
      tour_dist <- candidate_dist
      if(tour_dist < best_distance){
        best_tour <- tour
        best_distance <- tour_dist
      }
    }
    distances[i] <- tour_dist
  }
  results <- list('distances' = distances, 'probs' = probs, 
                  'best_tour' = best_tour, 'best_distance' = best_distance)
  
  return(results)
}

# function for plotting temperature change over iterations
plotTemp <- function(iterations, tempFunc, amp, loc, scale){
  newTempFunc <- function(x) tempFunc(x, amp, loc, scale)
  curve(newTempFunc, 1, iterations)
}

# data on 52 locations in berlin
b52 <- readNodeData('data/berlin52.tsp', 'data/berlin52.opt.tour')
b52_distanceMatrix <- as.matrix(dist(b52$nodes))

# compare models with random swap and neighbor preference swap
b52_results_random <- runAnnealing(b52$nodes, b52_distanceMatrix)
b52_results_prefN <- runAnnealing(b52$nodes, b52_distanceMatrix, preferNeighbors = T)
b52_best_result <- which.min(c(b52_results_random$best_distance,
                               b52_results_prefN$best_distance))
b52_results <- c(b52_results_random, b52_results_prefN)[b52_best_result]

# compare to optimal tour
b52_comparison <- compareOptimal(b52$nodes, b52_distanceMatrix, 
                                 b52_results$best_tour, b52$opt_tour)

# data on 76 cities in Germany
pr76 <- readNodeData('data/pr76.tsp', 'data/pr76.opt.tour')
pr76_distanceMatrix <- as.matrix(dist(pr76$nodes))

# compare models with and without biased swap 
pr76_results_random <- runAnnealing(pr76$nodes, pr76_distanceMatrix, 
                                    tempScale = 5000, tempAmplitude = 6000)
pr76_results_prefN <- runAnnealing(pr76$nodes, pr76_distanceMatrix, 
                                   tempScale = 5000, preferNeighbors = T)
pr76_best_result <- which.min(c(pr76_results_random$best_distance,
                                pr76_results_prefN$best_distance))
pr76_results <- c(pr76_results_random, pr76_results_prefN)[pr76_best_result]

# compare to optimal tour
pr76_comparison <- compareOptimal(pr76$nodes, pr76_distanceMatrix, 
                                  pr76_results$best_tour, pr76$opt_tour)

# data on 48 state capitals in US
att48 <- readNodeData('data/att48.tsp', 'data/att48.opt.tour')
att48_distanceMatrix <- as.matrix(dist(att48$nodes))

# compare models with and without biased swap 
att48_results_random <- runAnnealing(att48$nodes, att48_distanceMatrix)
att48_results_prefN <- runAnnealing(att48$nodes, att48_distanceMatrix)
att48_best_result <- which.min(c(att48_results_random$best_distance,
                                 att48_results_prefN$best_distance))
att48_results <- c(att48_results_random, att48_results_prefN)[att48_best_result]

# compare to optimal tour
att48_comparison <- compareOptimal(att48$nodes, att48_distanceMatrix, 
                                   att48_results$best_tour, att48$opt_tour)

