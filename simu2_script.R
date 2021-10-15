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

if (size == "small"){
  
  ## Discretizing continuous nodes either 2 or 4 levels
  conti_nodes = which(types_by_node == "c")
  two_level_nodes = conti_nodes[1:(n_conti/2)]
  four_level_nodes = conti_nodes[(n_conti/2+1):n_conti]
  
} else if (size == "large") {
  
  ## Discretizing continuous nodes either 2 or 4 levels
  
  # Select will-be discretizing nodes (20 nodes)
  disc_nodes = sample(which(types_by_node == "c"), n_ordin, replace = F)
  n_disc = length(disc_nodes)
  conti_raw_nodes = which(types_by_node == "c")
  conti_nodes = setdiff(conti_raw_nodes, disc_nodes) 
  two_level_nodes = disc_nodes[1:(n_disc/2)]
  four_level_nodes = disc_nodes[(n_disc/2+1):n_disc]
  
}

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
result = glmdag(data_input, n_lams = n_lams, eps_lam = eps_lam,  
                path_par = path_par, path_par_num = ell, verbose = T)
