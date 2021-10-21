# Debugging Block
if (F){
  if (Sys.info()['sysname'] == "Darwin") source("./master_script.R")
  seed_para = 1
  n_obs = 50 # the number of observations
  n_nodes = 30 # the number of total (multinomial) node
  n_lams = 30
  
}

# Randomly assign continuous or multinomial nodes
types_by_node = gen_node_types(n_conti + n_ordin, n_multi, 0, seed = seed_para)

# Generate the number of levels for each node randomly 
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
#.. it would be 1 if the node is continuous  

# Generate the grand parameter matrix W 
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)

# Generate datapoints
data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
#.. Note that the seed set with the 'iter' number

if (size == "small"){   # .. All the continuous nodes are discretized

  ## Discretizing continuous nodes either 3 or 5 levels
  # disc_nodes = which(types_by_node == "c")
  # n_disc = length(disc_nodes)
  # three_level_nodes = disc_nodes[1:(n_disc/2)]
  # five_level_nodes = disc_nodes[(n_disc/2+1):n_disc]
  
  ## Discretizing continuous nodes 5 levels
  disc_nodes = which(types_by_node == "c")
  n_disc = length(disc_nodes)
  five_level_nodes = disc_nodes[1:n_disc]
  
} else if (size == "large") { # .. Some of them are discretized
  
  ## Discretizing continuous nodes either 3 or 5 levels
  disc_nodes = sample(which(types_by_node == "c"), n_ordin, replace = F)
  n_disc = length(disc_nodes)
  three_level_nodes = disc_nodes[1:(n_disc/2)]
  five_level_nodes = disc_nodes[(n_disc/2+1):n_disc]
  
}

## Discretizing into 3 levels
# temp = data_input[three_level_nodes]
# temp1 = lapply(temp, function(x) 
#   cut(x, c(min(x)-1, quantile(as.matrix(x), probs = c(0.34, 0.67, 1))), 
#       labels = 1:3))

## Discretizing into 5 levels
temp = data_input[five_level_nodes]
temp2 = lapply(temp, function(x) 
  cut(x, c(min(x)-1, quantile(as.matrix(x), probs = c(0.2, 0.4, 0.6, 0.8, 1))), 
      labels = 1:5))

discretized_data = data.frame(temp2)

## Generate a different type of dataset according to the method
data_input = 
  switch(method,
         # Oringinal data
         "mc" = { data_input },
         # Consider discretized nodes as 'nominal-level' nodes
         "mm" = { data_input[disc_nodes] <- discretized_data; data_input},
         # Consider discretized nodes as 'ordinal-level' nodes
         "mo" = {
           ordered = lapply(discretized_data, function(x) factor(x, ordered = T))
           data_input[disc_nodes] <- ordered; data_input
         } )

## Run main functions

result = glmdag(data_input, n_lams = n_lams, eps_lam = eps_lam,  
                path_par = path_par, path_par_num = ell, verbose = T)
