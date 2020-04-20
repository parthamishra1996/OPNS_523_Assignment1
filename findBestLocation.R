find_best_location <- function(num.sites, assemblers, competetor, union.rate.fn, beta, num.tries){
  # Function to find best location for foreign firm, given the mentioned parameters. All the location co-ordinates must be 2D. 
  euclid_dist <- function(d1,d2){
    #'''Function to compute Euclidean distance betweeen two 2D vectors d1 and d2'''
    return(sqrt((d1[1] - d2[1])^2 + (d1[2] - d2[2])^2));
  }
  
  union_cost <- function(plant){
    #'''Function to compute the union component of cost'''
    return(union.rate.fn(plant[1], plant[2]));
  }
  
  Cbar_ri <- function(plant, assembler){
    #'''Function estimates the total cost of a plant to supply a part to an assembler  '''
    return(beta["distance"]*euclid_dist(plant, assembler) + beta["union"]*union_cost(plant));
  }
  
  
  #init = runif(2*num.sites)                    # Initialize the 6 locations of plants uniformly in [0,1]^2 for the foreign firm because the assemblers seems to be within this compact space
  solution = readRDS(paste0(var_save, 'solution.rds'))
  init = split(solution, seq(nrow(solution))) # [ An alternate initialization given]
  sol1 = unlist(init, use.names = TRUE)        # Converts list to a vector of rows(co-ordinates)
  
  #Based on the data, it is assumed that there is a single competing domestic firm with multiple plants at diff. locations
  objective <- function(sol1){
    sol2 = matrix(sol1, nrow= nrow(solution), ncol = ncol(solution))     # Create a matrix of size of input  
    
    tot = 0.0                                                            #Final obj value is accumulated here
    for (y in 1:nrow(assemblers)){
      var1 = 0.0                                                         #Cost component of the foreign firm
      var2 = 0.0                                                         #Cost component of the domestic competing firm
      
      for (z in 1:nrow(solution)){
        var1 = var1 + exp(-Cbar_ri(sol2[z,],assemblers[y,]))
      }
      
      for (w in 1:nrow(competetor)){
        var2 = var2 + exp(-Cbar_ri(competetor[w,],assemblers[y,]))
      }
      tot = tot - log(var1 +var2)
    }
    
    return(tot)  
  }
  
  res = optim(sol1, objective, control = list(maxit=10000))
  
  print("Optimal location: ")
  m = matrix(res$par, nrow=6, ncol=2)
  print(as.data.frame(m,row.names = NULL))
  print("Minimum cost: ")
  print(res$value)
  print("Convergence achieved?(Yes(1)/No(0): ")
  print(res$convergence)
}