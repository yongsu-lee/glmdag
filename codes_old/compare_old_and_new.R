rm(list=ls())

source("load_r_pkgs.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") }
BiocManager::install()
BiocManager::install(c("graph", "RBGL"))

if(!require(pcalg)){
  install.packages("pcalg", quiet = T)
  library("pcalg", quietly = T)
}


# Set some global parameters ++++
seed_para = 1
n_obs = 20 # the number of observations
n_nodes = 10 # the number of total (multinomial) nodes
types_by_node = rep("m", n_nodes)
graph_type = "bi"
seed_para = 1
seed_iter = 1

## Old functions ####

# Load functions 
source("gen_adj_mat.R")
source("shared_fun.R")
source("gen_para.R")
source("gen_data.R")

source("multi.gen_data.R")

(A_true_old = gen_adj_mat(n_nodes, graph_type, seed = seed_para))
(n_edges_old = sum(A_true_old))

n_resps_by_node = gen_node_levels(types_by_node, max_n_levels = 4, 
                                  seed = seed_para)

(W_true_old = gen_para_old(A_true_old, types_by_node,
                          n_resps_by_node, seed = seed_para))

data_input_old = gen_data(n_obs, A_true, W_true, seed = seed_data)

