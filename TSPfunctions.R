library(Rcpp)
library(ggplot2)

# function for reading in a TSP data set
readNodeData <- function(node_filepath, opt_tour_filepath){
  
  options(warn = -1)
  
  lines <- readLines(node_filepath)
  entries <- strsplit(lines, '\\s+') 
  start <- which(sapply(entries, function(e) e[e != ''][1] == '1'))
  end <- which(sapply(entries, function(e) 'EOF' %in% e)) - 1
  entries <- entries[start:end]
  entries <- lapply(entries, function(entry) entry[2:3])
  table <- do.call(rbind, entries)
  colnames(table) <- c('c1', 'c2')
  nodes <- as.matrix(apply(table, 2, as.numeric))

  opt_tour <- readLines(opt_tour_filepath)
  opt_tour <- as.integer(opt_tour)
  opt_tour <- opt_tour[!is.na(opt_tour) & opt_tour != -1]
  
  data <- list('nodes' = nodes, 'opt_tour' = opt_tour)
  
  return(data)
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

# function for plotting a tour
plotTour <- function(nodes, tour){
  
  tour_nodes <- nodes[tour,]
  tour_nodes <- rbind(tour_nodes, tour_nodes[1,])
  
  ggplot(data.frame('c1' = nodes[,1], 'c2' = nodes[,2]), aes(x = c1, y = c2)) + 
    geom_point(size = 2) + 
    geom_path(data = data.frame('c1' = tour_nodes[,1],
                                'c2' = tour_nodes[,2]),
              aes(x = c1, y = c2), col = 'red') + 
    theme(text = element_text(size = 26),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
}

# function for comparing a tour to the optimal tour
compareOptimal <- function(nodes, distMat, tour, opt_tour){
  
  best_tour_distance <- round(getTourDistance(distMat, tour), 2)
  opt_tour_distance <- round(getTourDistance(distMat, opt_tour), 2)
  
  writeLines(paste('Final tour distance:', best_tour_distance))
  writeLines(paste('Optimal tour distance:', opt_tour_distance))
  
  return(list('best' = best_tour_distance, 'opt' = opt_tour_distance))
}

# function for plotting distances from a TSP optimization
plotDistances <- function(distances){
  distance_df <- data.frame('Iteration' = 1:length(distances),
                            'Distance' = distances)
  dist_plot <- ggplot(distance_df, aes(x = Iteration, y = Distance)) +
    labs(x = '\nIteration', y = 'Distance\n') +
    geom_line() + theme(text = element_text(size = 26),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line = element_line(color = 'black'))
  
  return(dist_plot)
}

