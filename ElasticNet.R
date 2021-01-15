library(Rcpp)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dichromat)

source('TSPfunctions.R')
sourceCpp('Updates.cpp')

# function to get the initial path 
initialPath <- function(n, x_scale, y_scale, x_bias, y_bias){
  matrix(c(x_scale * cos(2 * (1:n) * pi / n) + x_bias,
           y_scale * sin(2 * (1:n) * pi / n) + y_bias),
           n, 2)
}

# get a value of K for the given iteration
getK <- function(i, scale = 10){
  1 / exp(i / scale)
}

# function for plotting temperature change over iterations
plotK <- function(iterations, KFunc, scale){
  newKFunc <- function(x) KFunc(x, scale)
  K_plot <- ggplot(data.frame('Iteration' = 1:iterations), aes(x = Iteration)) +
    labs(x = '\nIteration', y = 'K\n') +
    stat_function(fun = newKFunc) + 
    theme(text = element_text(size = 26),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.ticks = element_blank())
  return(K_plot)
}

# run the elastic net for a set of nodes
runElasticNet <- function(nodes, Mfactor = 1.5, n_iterations = 100, 
                          alpha = 1, beta = 1, Kscale = 10, plot = T,
                          n_checkpoints = 10){
  
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
  
  if(plot){
    checkpoints <- unique(round(seq(1, n_iterations, length.out = n_checkpoints)))
  }
  
  # list for saving intermediate paths
  paths <- list()
  for(i in 1:n_iterations){
    if(plot){
      if(i %in% checkpoints){
        idx <- which(checkpoints == i)
        paths[[idx]] <- plotElasticNet(nodes, ring)
        if(interactive()) print(paths[[idx]])
      }
    }
    weights <- weightUpdate(nodes, ring, getK(i, scale = Kscale))
    delta_y <- yUpdate(nodes, ring, weights, getK(i, scale = Kscale), alpha, beta)
    ring <- ring + delta_y
  }
  
  ring <- t(apply(ring, 1, function(row) row * colsds))
  ring <- t(apply(ring, 1, function(row) row + colmeans))
  
  return(list('final_path' = ring, 'paths' = paths))
  
}

plotElasticNet <- function(nodes, path){
  
  path <- rbind(path, path[1,])
  
  ggplot(data.frame('c1' = nodes[,1], 'c2' = nodes[,2]), aes(x = c1, y = c2)) + 
    geom_point(size = 2) + 
    geom_path(data = data.frame('c1' = path[,1],
                                'c2' = path[,2]),
              aes(x = c1, c2), col = 'red') + 
    theme(text = element_text(size = 26),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
}

# function for searching for optimal values for amplitude and scale
gridSearch <- function(nodes, Mfactors, betas, Kscales, n_iterations = 100){
  
  best_dist <- Inf
  heatmap <- matrix(Inf, nrow = length(Kscales), ncol = length(betas),
                    dimnames = list(Kscales, betas))
  for(b in betas){
    for(scale in Kscales){
      for(Mf in Mfactors){
        results <- runElasticNet(nodes, Mfactor = Mf, beta = b, Kscale = scale, 
                                 n_iterations = n_iterations)
        distance <- getTourDistance(as.matrix(dist(results$final_path)), 1:nrow(results$final_path))
        current_distance <- heatmap[rownames(heatmap) == scale, colnames(heatmap) == b]
        if(distance < current_distance){
          heatmap[rownames(heatmap) == scale, colnames(heatmap) == b] <- distance
        }
        if(distance < best_dist){
          optimal_params <- list('beta' = b, 'Kscale' = scale, 'Mfactor' = Mf)
          best_dist <- distance
          optimal <- list('best_tour' = results$final_path, 
                          'best_distance' = round(distance, 2),
                          'paths' = results$paths)
        } 
      }
    }
  }
  
  return(list('results' = optimal, 'params' = optimal_params, 'heatmap' = heatmap))
}

# function for running grid search and making figures
runPipeline <- function(dataset_name, node_filepath, opt_tour_filepath, 
                        Mfactors, betas, Kscales, n_iterations){
  
  data <- readNodeData(node_filepath, opt_tour_filepath)
  distanceMatrix <- as.matrix(dist(data$nodes))
  
  # run elastic net
  results <- gridSearch(data$nodes, Mfactors, betas, Kscales, n_iterations = n_iterations)
  
  # time the run
  pt <- proc.time()
  runElasticNet(data$nodes, Mfactor = results$params$Mfactor,
                beta = results$params$beta, 
                Kscale = results$params$Kscale)
  time <- proc.time() - pt
  writeLines(paste('Time elapsed finding best path:', round(time['elapsed'], 3)))
  
  # plot tour
  tourplot <- plotElasticNet(data$nodes, results$results$best_tour)
  ggsave(paste('figures/ElasticNet/', dataset_name, '_tour.png', sep = ''), tourplot)
  
  # plot intermediate paths
  for(i in 1:length(results$results$paths)){
    ggsave(paste('figures/ElasticNet/', dataset_name, '_', i, '.png', sep = ''), results$results$paths[[i]])
  }
  
  # plot K
  Kplot <- plotK(n_iterations, getK, results$params$Kscale)
  ggsave(paste('figures/ElasticNet/', dataset_name, '_K.png', sep = ''))
  
  # normalise heatmap
  results$heatmap <- ((results$heatmap / getTourDistance(distanceMatrix, data$opt_tour)) - 1) * 100
  # plot heatmap from grid search
  heatmap <- pheatmap::pheatmap(results$heatmap, cluster_rows = F, 
                                cluster_cols = F, angle_col = 0,
                                color = colorRampPalette(rev(colorschemes$BluetoGreen.14))(100),
                                filename = paste('figures/ElasticNet/', dataset_name, 
                                                 '_gridsearch.png', sep = ''), 
                                cellwidth = 30, cellheight = 30)
  
  writeLines(paste('Tours for the', dataset_name, 
                   'data set completed, optimal model found using the following parameters:',
                   '\nBeta:', results$params$beta,
                   '\nK scale:', results$params$Kscale,
                   '\nM factor', results$params$Mf))
  
  writeLines(paste('Final tour distance, ', dataset_name, ': ', 
                   results$results$best_distance, sep = ''))
  writeLines(paste('Optimal tour distance, ', dataset_name, ': ', 
                   round(getTourDistance(distanceMatrix, data$opt_tour)), sep = ''))
}

Mfactors <- c(1.5, 1.75, 2)
betas <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
Kscales <- seq(10, 20, 2)
n_iterations <- 100

# data on 52 locations in berlin
writeLines(paste('Loading b52 data set, running elastic net for', n_iterations, 'iterations...'))
runPipeline('b52', 'data/berlin52.tsp', 'data/berlin52.opt.tour', 
            Mfactors, betas, Kscales, n_iterations)

# data on 76 cities in Germany
writeLines(paste('Loading pr76 data set, running elastic net for', n_iterations, 'iterations...'))
runPipeline('pr76', 'data/pr76.tsp', 'data/pr76.opt.tour', 
            Mfactors, betas, Kscales, n_iterations)

# data on 48 state capitals in US
writeLines(paste('Loading att48 data set, running elastic net for', n_iterations, 'iterations...'))
runPipeline('att48', 'data/att48.tsp', 'data/att48.opt.tour',
            Mfactors, betas, Kscales, n_iterations)
