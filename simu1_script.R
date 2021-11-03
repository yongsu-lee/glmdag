# Debugging Block
if (F){
  if (Sys.info()['sysname'] == "Darwin") source("./master_script.R")
  seed_para = 1
  n_obs = 50 # the number of observations
  n_nodes = 30 # the number of total (multinomial) node
  n_lams = 30
  
}

# Set type of each nodes
types_by_node = rep("m", n_nodes)

# Generate the number of levels for each node randomly 
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
#.. it would be 1 if the node is continuous  

# Generate the grand parameter matrix W 
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)

# Generate datapoints
data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
#.. Note that the seed set with the 'iter' number

# Convert dataset for cd algorithm
data_input_cd = conv_to_cd_data(data_input) 

## Generate result containers

## Run main functions - GLMDAG
result_glmdag = glmdag(data_input, n_lams = n_lams, eps_lam = eps_lam, 
                       verbose = T)

## Run main functions - CD
result_cd = cd.run(indata = data_input_cd)

# Combine two results
result <- list(glmdag = result_glmdag, cd = result_cd)

## Save the result
saveRDS(result, file = paste0(name_temp,".rds"))
