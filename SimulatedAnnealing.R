source('TSPfunctions.R')

# get the temperature for iteration x
getTemp <- function(x, amp = 4000, loc = 0, scale = 3000){
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
                         maxBetaShape2 = 20, preferNeighbors = T, s1 = 2, s2 = 12,
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
    temp <- getTemp(i, amp = tempAmplitude, scale = tempScale)
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

plotTemp <- function(iterations, tempFunc){
  curve(tempFunc, 1, iterations)
}

# berlin52 data 
berlin52 <- readNodeData('data/berlin52.tsp')
b52_opt_tour <- readLines('data/berlin52.opt.tour')
b52_opt_tour <- as.integer(b52_opt_tour)
b52_opt_tour <- b52_opt_tour[!is.na(b52_opt_tour) & b52_opt_tour != -1]

b52_distanceMatrix <- as.matrix(dist(berlin52))

b52_results <- runAnnealing(berlin52, b52_distanceMatrix)

plotTour(berlin52, b52_results$best_tour)

b52_comparison <- compareOptimal(berlin52, b52_distanceMatrix, 
                                 b52_results$best_tour, b52_opt_tour)

# pr76 data
pr76 <- readNodeData('data/pr76.tsp')
pr76_opt_tour <- readLines('data/pr76.opt.tour')
pr76_opt_tour <- as.integer(pr76_opt_tour)
pr76_opt_tour <- pr76_opt_tour[!is.na(pr76_opt_tour) & pr76_opt_tour != -1]

pr76_distanceMatrix <- as.matrix(dist(pr76))

pr76_results <- runAnnealing(pr76, pr76_distanceMatrix)

plotTour(pr76, pr76_results$best_tour)

pr76_comparison <- compareOptimal(pr76, pr76_distanceMatrix, 
                                  pr76_results$best_tour, pr76_opt_tour)

