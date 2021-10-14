#######################################################
## Table 3 (Small Graphs, Nom vs Ord, Results) ####
#######################################################

if(F){ # for debugging, delete when submitting
  
  rm(list=ls())
  source("./codes/hhmc_00_load_ftns.R")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  
  queue_arg = read.table("./codes/queue_list_mmhc", sep =",", strip.white = T)
  queue_args = queue_arg[1,]
  graph_type = as.character(queue_args[1])
  method = as.character(queue_args[2])
  iter = as.integer(queue_args[3])

  save_file_name = paste0("./results/simu2_small/simu2_small_",
                          graph_type,"_",method,"_",iter,"_mmhc.rds")
  
}

if(T){ # for server running
  
  # Clear the memory
  rm(list=ls())
  
  # Load required packages and functions
  source("./hhmc_00_load_ftns.R")
  
  # Set seed generating system
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  
  #!/usr/bin/env Rscript
  queue_args = commandArgs(trailingOnly=TRUE)
  graph_type = as.character(queue_args[1])
  method = as.character(queue_args[2])
  iter = as.integer(queue_args[3])
  save_file_name = paste0("simu2_small_", graph_type, "_", method, "_", 
                          iter, "_mmhc.rds")

}


## Set up parameters
n_obs = 50 # the number of observations
n_conti = 6 # the number of continuous nodes
n_multi = 4 # the number of multinomial nodes
n_nodes = n_conti + n_multi # the number of total nodes


## Generate graphs, parameters and datapoints 

# Generate the true graph and obtain corresponding adjacency matrix
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = 1)
graph_true = graph_set$graph_true
A_true = graph_set$A_true

# Randomly assign continuous or multinomial nodes
types_by_node = gen_node_types(n_conti, n_multi, 0, seed = 1)

# Generate the number of levels for each node randomly 
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = 1)
#.. it would be 1 if the node is continuous  

# Generate the grand parameter matrix W 
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = 1)

# Generate datapoints
data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
#.. Note that the seed set with the 'iter' number


## Discretizing continuous nodes either 2 or 4 levels
conti_nodes = which(types_by_node == "c")
two_level_nodes = conti_nodes[1:(n_conti/2)]
four_level_nodes = conti_nodes[(n_conti/2+1):n_conti]

# Discretizing into 2 levels
temp = data_input[two_level_nodes]
temp1 = lapply(temp, function(x) 
  cut(x, c(min(x)-1, median(x), max(x)), labels = 1:2))

# Discretizing into 4 levels
temp = data_input[four_level_nodes]
temp2 = lapply(temp, function(x) 
  cut(x, c(min(x)-1, quantile(as.matrix(x), probs = c(0.25, 0.5, 0.75, 1))), 
      labels = 1:4))

discretized_data = data.frame(temp1, temp2)


## Generate a different type of dataset according to the method
data_input = 
  switch(method,
         # Oringinal data
         "mc" = { data_input },
         # Consider discretized nodes as 'nominal-level' nodes
         "mm" = { data_input[conti_nodes] <- discretized_data; data_input},
         # Consider discretized nodes as 'ordinal-level' nodes
         "mo" = {
           ordered = lapply(discretized_data, function(x) factor(x, ordered = T))
           data_input[conti_nodes] <- ordered; data_input
         } )


## Run main functions
A_est_hc = hc(data_input)
A_est_mmhc = mmhc(data_input)
result = list(hc = A_est_hc, mmhc = A_est_mmhc)

## Save the result
saveRDS(result, file = save_file_name)
