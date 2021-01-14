library(Rcpp)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dichromat)

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
                         mutation_rate = 0.2, n_iterations = 1000, 
                         state = 123, plot = T, n_checkpoints = 100,
                         track_distance = F){
  
  if(plot){
    checkpoints <- unique(round(seq(1, n_iterations, length.out = n_checkpoints)))
  }
  
  set.seed(state)
  n_nodes <- nrow(nodes)
  # generate random initial population
  tours <- t(replicate(population, sample(n_nodes)))
  # keep sampling until no duplicates are present
  while(any(duplicated(tours))){
    tours[duplicated(tours),] <- replicate(sum(duplicated(tours)), sample(n_nodes))
  }
  d <- numeric(length = n_iterations)
  for(i in 1:n_iterations){
    if(plot){
      if(i %in% checkpoints){
        distances <- apply(tours, 1, getTourDistance, d = distMat)
        best_tour <- tours[which.min(distances),]
        plotTour(nodes, best_tour)
      }
    }
    if(track_distance){
      distances <- apply(tours, 1, getTourDistance, d = distMat)
      d[i] <- min(distances)
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
                  'best_distance' = distances[which.min(distances)])
  if(track_distance){
    results[['distances']] <- d
  }
  
  return(results)
  
}

# function for searching for optimal values for amplitude and scale
gridSearch <- function(nodes, distMat, populations, mutation_rates, n_iterations = 1000){
  
  best_dist <- Inf
  heatmap <- matrix(0, nrow = length(populations), ncol = length(mutation_rates),
                    dimnames = list(populations, mutation_rates))
  for(pop in populations){
    for(m in mutation_rates){
      results_ins <- runEvolution(nodes, distMat, population = pop, mutation_rate = m,
                                  mutate_func = insertionMutation, plot = F)
      results_inv <- runEvolution(nodes, distMat, population = pop, mutation_rate = m,
                                  mutate_func = simpleInversionMutation, plot = F)
      best_result <- which.min(c(results_ins$best_distance,
                                 results_inv$best_distance))
      if(best_result == 1){ best_results <- results_ins } else { best_results <- results_inv }
      heatmap[rownames(heatmap) == pop, colnames(heatmap) == m] <- best_results$best_distance
      if(best_results$best_distance < best_dist){
        optimal_params <- list('mutation_rate' = m, 'population' = pop, 
                               'mutation' = if(best_result == 1) 'insertion' else 'inversion')
        optimal <- best_results
        best_dist <- best_results$best_distance
      } 
    }
  }
  
  return(list('results' = optimal, 'params' = optimal_params, 'heatmap' = heatmap))
}

# function for running grid search and making figures
runPipeline <- function(dataset_name, node_filepath, opt_tour_filepath, 
                        populations, mutation_rates, n_iterations){
  
  data <- readNodeData(node_filepath, opt_tour_filepath)
  distanceMatrix <- as.matrix(dist(data$nodes))
  
  # run grid search
  results <- gridSearch(data$nodes, distanceMatrix, populations, mutation_rates, n_iterations)
  
  # time the run
  pt <- proc.time()
  runEvolution(data$nodes, distanceMatrix, 
               if(results$params$mutation == 'insertion') insertionMutation 
               else simpleInversionMutation,
               mutation_rate = results$params$mutation_rate, 
               population = results$params$population)
  time <- proc.time() - pt
  writeLines(paste('Time elapsed finding best path:', round(time['elapsed'], 3)))
  
  # plot best tour
  tourplot <- plotTour(data$nodes, results$results$best_tour)
  ggsave(paste('figures/GeneticAlgorithm/', dataset_name, '_tour.png', sep = ''), tourplot)
  
  # plot distances
  best_run <- runEvolution(data$nodes, distanceMatrix, 
                           if(results$params$mutation == 'insertion') insertionMutation 
                           else simpleInversionMutation,
                           mutation_rate = results$params$mutation_rate, 
                           population = results$params$population, 
                           track_distance = T)
  distplot <- plotDistances(best_run$distances)
  ggsave(paste('figures/GeneticAlgorithm/', dataset_name, '_distances.png', sep = ''), distplot)
  
  # normalise heatmap
  results$heatmap <- ((results$heatmap / getTourDistance(distanceMatrix, data$opt_tour)) - 1) * 100
  # plot heatmap from grid search
  heatmap <- pheatmap::pheatmap(results$heatmap, cluster_rows = F, 
                                cluster_cols = F, angle_col = 0,
                                color = colorRampPalette(rev(colorschemes$BluetoGreen.14))(100),
                                filename = paste('figures/GeneticAlgorithm/', dataset_name, 
                                                 '_gridsearch.png', sep = ''), 
                                cellwidth = 30, cellheight = 30)
  
  writeLines(paste('Tours for the', dataset_name, 
                   'data set completed, optimal model found using the following parameters:',                             '\nPopulation size:', results$params$population,
                   '\nMutation rate:', results$params$mutation_rate,
                   '\nMutation operator:', results$params$mutation))
  
  # compare to optimal tour
  compareOptimal(data$nodes, distanceMatrix, 
                 results$results$best_tour, data$opt_tour)
  
}

# population values and mutation rates to test for each data set
populations <- seq(1000, 3000, 500)
mutation_rates <- seq(0.1, 0.4, 0.05)
n_iterations <- 1000

# data on 52 locations in berlin
writeLines(paste('Loading b52 data set, running genetic algorithm for', 
                 n_iterations, 'iterations...'))
runPipeline('b52', 'data/berlin52.tsp', 'data/berlin52.opt.tour',
            populations, mutation_rates, n_iterations)

# data on 76 cities in Germany
writeLines(paste('Loading pr76 data set, running genetic algorithm for', 
                 n_iterations, 'iterations...'))
runPipeline('pr76', 'data/pr76.tsp', 'data/pr76.opt.tour',
            populations, mutation_rates, n_iterations)

# data on 48 state capitals in US
writeLines(paste('Loading att48 data set, running genetic algorithm for', 
                 n_iterations, 'iterations...'))
runPipeline('att48', 'data/att48.tsp', 'data/att48.opt.tour',
            populations, mutation_rates, n_iterations)