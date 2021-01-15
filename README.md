# Travelling salesman problem 

This repository contains code for three algorithms implemented to solve the TSP:
- [Simulated annealing](https://toddwschneider.com/posts/traveling-salesman-with-simulated-annealing-r-and-shiny/)
- [Genetic algorithms](https://link.springer.com/article/10.1023/A:1006529012972)
- [Elastic net](https://pubmed.ncbi.nlm.nih.gov/3561510/)

Data can be found [here](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/). The performance of the three algorithms was tested on three data sets, included in the `data/` directory: att48 (48 state capitals in the US), berlin52 (52 locations in Berlin), and pr76 (76 cities around the world).

The reproduce the analysis, run the scripts from the command line: 

```
Rscript SimulatedAnnealing.R
```
```
Rscript GeneticAlgorithm.R
```
```
Rscript ElasticNet.R
```

The scripts automatically save all relevant figures and a .txt file containing a summary of the results for each data set. 

If in RStudio, the progression of each algorithm can be visualised by running each function (`runAnnealing()`, `runEvolution()`, and `runElasticNet()`) with `plot = T`, e.g.

``` {R}
source('ElasticNet.R')
berlin52 <- readNodeData('data/berlin52.tsp', 'data/berlin52.opt.tour')
results <- runElasticNet(berlin52$nodes, plot = T)
```