## Set up ####

rm(list=ls())


# Set some global parameters ++++
seed_para = 1
n_obs = 20 # the number of observations
n_nodes = 10 # the number of total (multinomial) nodes
types_by_node = rep("m", n_nodes)
graph_type = "bi"
seed_para = 1
seed_iter = 1

# Load functions 
source("./codes/00_load_ftns.R")

# Set some parameters 
types_by_node = rep("m", n_nodes)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)

# Generate the true graph and obtain corresponding adjacency matrix
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
(A_true = graph_set$A_true)
n_edge = sum(A_true) 


(W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para))
data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = seed_iter)
data_input
