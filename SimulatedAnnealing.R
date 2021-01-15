library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dichromat)

source('TSPfunctions.R')

# get the temperature for iteration x
getTemp <- function(x, amp, scale){
  amp * (1 / (1 + exp(x / scale)))
}

# get probability of accepting a new tour, given two distances and temperature
getProb <- function(old_dist, new_dist, temp){
  min(1, exp((old_dist - new_dist) / temp))
}


# function for running simulated annealing
runAnnealing <- function(nodes, distMat, tempAmplitude, tempScale, 
                         n_iterations = 50000, state = 123,
                         plot = T, n_checkpoints = 100){
  
  if(plot){
    checkpoints <- unique(round(seq(1, n_iterations, length.out = n_checkpoints)))
  }
  
  set.seed(state)
  n_nodes <- nrow(nodes)
  # start with a random tour
  tour <- sample(n_nodes)
  tour_dist <- getTourDistance(distMat, tour)
  best_tour <- tour
  best_distance <- Inf
  # initialise vectors for tracking probabilities and distances
  probs <- numeric(length = n_iterations)
  distances <- numeric(length = n_iterations)
  # iterate the annealing process
  for(i in 1:n_iterations){
    if(plot & interactive()){
      if(i %in% checkpoints){
        plotTour(nodes, tour)
      }
    }
    temp <- getTemp(i, amp = tempAmplitude, scale = tempScale)
    # get candidates, drawn randomly depending on the current iteration
    swapped <- sample(n_nodes, size = 2, replace = F)
    candidate_tour <- tour
    # swap the selected points in the candidate tour
    candidate_tour[swapped[1]:swapped[2]]<- rev(tour[swapped[1]:swapped[2]])
    candidate_dist <- getTourDistance(distMat, candidate_tour)
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
plotTemp <- function(iterations, tempFunc, amp, scale){
  newTempFunc <- function(x) tempFunc(x, amp, scale)
  temp_plot <- ggplot(data.frame('Iteration' = 1:iterations), aes(x = Iteration)) +
    labs(x = '\nIteration', y = 'Temperature\n') +
    stat_function(fun = newTempFunc) + 
    theme(text = element_text(size = 26),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(color = 'black'))
  return(temp_plot)
}

# function for searching for optimal values for amplitude and scale
gridSearch <- function(nodes, distMat, amplitudes, scales, n_iterations = 50000){
  
  best_dist <- Inf
  heatmap <- matrix(Inf, nrow = length(amplitudes), ncol = length(scales),
                    dimnames = list(amplitudes, scales))
  for(amp in amplitudes){
    for(scale in scales){
      results <- runAnnealing(nodes, distMat, tempAmplitude = amp,
                                     n_iterations = n_iterations,
                                     tempScale = scale, plot = F)
      heatmap[rownames(heatmap) == amp, colnames(heatmap) == scale] <- results$best_distance
      if(results$best_distance < best_dist){
        optimal_params <- list('amp' = amp, 'scale' = scale)
        optimal <- results
        best_dist <- results$best_distance
      } 
    }
  }
  
  return(list('results' = optimal, 'params' = optimal_params, 'heatmap' = heatmap))
}

# function for running grid search and making figures
runPipeline <- function(dataset_name, node_filepath, opt_tour_filepath, 
                        amplitudes, scales, n_iterations){
  
  data <- readNodeData(node_filepath, opt_tour_filepath)
  distanceMatrix <- as.matrix(dist(data$nodes))
  
  # run elastic net
  results <- gridSearch(data$nodes, distanceMatrix, 
                        amplitudes, scales, n_iterations = n_iterations)
  
  # time the run
  pt <- proc.time()
  runAnnealing(data$nodes, distanceMatrix, results$params$amp, results$params$scale)
  time <- proc.time() - pt
  writeLines(paste('Time elapsed finding best path:', round(time['elapsed'], 3)))
    
  # plot tour
  tourplot <- plotTour(data$nodes, results$results$best_tour)
  ggsave(paste('figures/SimulatedAnnealing/', dataset_name, '_tour.png', sep = ''), tourplot)
  
  # plot optimal tour
  tourplot <- plotTour(data$nodes, data$opt_tour)
  ggsave(paste('figures/', dataset_name, '_opt_tour.png', sep = ''), tourplot)
  
  # plot temperature
  tempplot <- plotTemp(50000, getTemp, results$params$amp, results$params$scale)
  tempplot <- tempplot + scale_x_continuous(labels = c(0, paste(1:(n_iterations / 10000), '0k', sep = '')))
  ggsave(paste('figures/SimulatedAnnealing/', dataset_name, '_temperature.png', sep = ''))
  
  # plot distances 
  distplot <- plotDistances(results$results$distances)
  distplot <- distplot + scale_x_continuous(labels = c(0, paste(1:(n_iterations / 10000), '0k', sep = '')))
  ggsave(paste('figures/SimulatedAnnealing/', dataset_name, '_distance.png', sep = ''), distplot)
  
  # normalise heatmap
  results$heatmap <- ((results$heatmap / getTourDistance(distanceMatrix, data$opt_tour)) - 1) * 100
  # plot heatmap from grid search
  heatmap <- pheatmap::pheatmap(results$heatmap, cluster_rows = F, 
                                cluster_cols = F, angle_col = 0,
                                color = colorRampPalette(rev(colorschemes$BluetoGreen.14))(100),
                                filename = paste('figures/SimulatedAnnealing/', dataset_name, 
                                                 '_gridsearch.png', sep = ''), 
                                cellwidth = 30, cellheight = 30)
  
  writeLines(paste('Tours for the', dataset_name, 
                   'data set completed, optimal model found using the following parameters:',
                   '\nTemperature function amplitude:', results$params$amp,
                   '\nTemperature function scale parameter:', results$params$scale,
                   '\nRandom selection?', results$params$prefN))
  
  compareOptimal(data$nodes, distanceMatrix, 
                 results$results$best_tour, data$opt_tour)
  
}

amplitudes <- seq(1000, 10000, by = 1000)
scales <- seq(1000, 10000, by = 1000)
n_iterations <- 50000

# data on 52 locations in berlin
writeLines(paste('Loading b52 data set, running simulated annealing for', 
                 n_iterations, 'iterations...'))
runPipeline('b52', 'data/berlin52.tsp', 'data/berlin52.opt.tour',
            amplitudes, scales, n_iterations)

# data on 76 cities in Germany
writeLines(paste('Loading pr76 data set, running simulated annealing for', 
                 n_iterations, 'iterations...'))
runPipeline('pr76', 'data/pr76.tsp', 'data/pr76.opt.tour',
            amplitudes, scales, n_iterations)

# data on 48 state capitals in US
writeLines(paste('Loading att48 data set, running simulated annealing for', 
                 n_iterations, 'iterations...'))
runPipeline('att48', 'data/att48.tsp', 'data/att48.opt.tour',
            amplitudes, scales, n_iterations)


