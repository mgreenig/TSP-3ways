#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix partiallyMappedCrossover(IntegerVector t1, IntegerVector t2, int swap1, int swap2){
  
  // convert swap indices to 0 index
  int s1 = swap1 - 1;
  int s2 = swap2 - 1;
  
  IntegerVector offspring1 (t1.size());
  IntegerVector offspring2 (t2.size());
  IntegerMatrix offspring (2, t1.size());  
  
  std::unordered_map<int, int> offspring1_mapping;
  std::unordered_map<int, int> offspring2_mapping;
  
  // insert elements from the swap into the other offspring
  for(int i = s1; i <= s2; i++){
    offspring1[i] = t2[i];
    offspring1_mapping.insert(std::make_pair(t2[i], t1[i]));
    offspring2[i] = t1[i];
    offspring2_mapping.insert(std::make_pair(t1[i], t2[i]));
  }
  
  int repl1;
  int repl2; 
  
  // for other elements, add using the mapping
  for(int i = 0; i < offspring1.size(); i++){
    // if in the swapped region, don't make any changes
    if(i >= s1 && i <= s2){
      continue;
    }
    // start with the node from tour 1 as the default
    repl1 = t1[i];
    // if the node is already in the offspring tour, use the mapping to find a node that is not
    while(std::find(offspring1.begin(), offspring1.end(), repl1) != offspring1.end()){
      repl1 = offspring1_mapping[repl1];
    }
    offspring1[i] = repl1;
    // same for tour 2
    repl2 = t2[i];
    while(std::find(offspring2.begin(), offspring2.end(), repl2) != offspring2.end()){
      repl2 = offspring2_mapping[repl2];
    }
    offspring2[i] = repl2;
  }
  
  offspring.row(0) = offspring1;
  offspring.row(1) = offspring2;
  
  return(offspring);
  
}

// [[Rcpp::export]]
IntegerMatrix runCrossover(IntegerMatrix tours, IntegerMatrix swaps, IntegerVector ordering){
  
  IntegerVector t1;
  IntegerVector t2;
  IntegerMatrix offspring (tours.nrow(), tours.ncol());
  IntegerMatrix new_offspring;
  int swap1;
  int swap2;
  
  // loop through pairs of tours and run crossover
  for(int i = 0; i < tours.nrow(); i = i + 2){
    t1 = tours(ordering[i] - 1, _);
    t2 = tours(ordering[i+1] - 1, _);
    swap1 = swaps((i / 2), 0) - 1;
    swap2 = swaps((i / 2), 1) - 1;
    new_offspring = partiallyMappedCrossover(t1, t2, 
                                             std::min(swap1, swap2),
                                             std::max(swap1, swap2));
    offspring(i, _) = new_offspring(0, _);
    offspring(i+1, _) = new_offspring(1, _);
  }
  
  return(offspring);
  
}

