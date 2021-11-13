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
  
  ## Discretizing continuous nodes 5 levels
  disc_nodes = which(types_by_node == "c")
  n_disc = length(disc_nodes)
  five_level_nodes = disc_nodes[1:n_disc] # not necessary for five-level only
  
} else if (size == "large") { # .. Some of them are discretized
  
  ## Discretizing continuous nodes either 5 levels
  disc_nodes = sample(which(types_by_node == "c"), n_ordin, replace = F)
  n_disc = length(disc_nodes)
  five_level_nodes = disc_nodes[1:n_disc] # not necessary for five-level only
  
}

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

# glmdag
if (simu_case == "simu2") {
  result = glmdag(data_input, n_lams = n_lams, eps_lam = eps_lam,  
                  path_par = path_par, path_par_num = ell, verbose = T)
  saveRDS(result, file = paste0(name_temp,".rds"))
}

# mmhc for "mm" type data
if (simu_case %in% c("simu2", "simu2_mmhc_only")) 
{
  if (method == "mm") {
    if (path_par == F | (path_par == T & ell == 1)) {
      A_est_hc = hc(data_input)
      A_est_mmhc = mmhc(data_input)
      result.mmhc = list(hc = A_est_hc, mmhc = A_est_mmhc)
      saveRDS(result.mmhc, file = paste0(name_temp, "_mmhc.rds") )
    }
  }
}