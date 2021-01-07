source('TSPfunctions.R')

# get the temperature for iteration x
getTemp <- function(x, amp = 4000, loc = 0, scale = 3000){
  amp * (1 / (1 + exp((x - loc) / scale)))
}

# get probability of accepting a new tour, given two distances and temperature
getProb <- function(old_dist, new_dist, temp){
  min(1, exp((old_dist - new_dist) / temp))
}


# function for running simulated annealing
runAnnealing <- function(nodes, distMat, n_iterations = 50000, 
                         state = 123, tempAmplitude = 4000, tempScale = 3000,
                         maxBetaShape2 = 20, preferNeighbors = T, s1 = 2, s2 = 12){
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
  # getCandidates function uses either adaptive selection or random sampling
  if(preferNeighbors){
    getCandidates <- function(n) selectCandidates(n, s1 = s1, s2 = s2)
  } else {
    getCandidates <- function(n) sample(n, size = 2, replace = F)
  }
  # iterate the annealing process
  for(i in 1:n_iterations){
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

plotNeighborDistances <- function(nodes, exact = T){
  if(exact){
    n <- nrow(nodes)
    xs <- 1:(n-1)
    probs <- 2*(n - xs) / (n * (n - 1))
    plot(xs, probs)
    lines(xs, probs)
  }
}

nodes <- readNodeData('data/berlin52.tsp')

distanceMatrix <- as.matrix(dist(nodes))

results <- runAnnealing(nodes, distanceMatrix, preferNeighbors = F)

plotTour(nodes, results$best_tour)
