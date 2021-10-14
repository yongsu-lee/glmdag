setwd("~/iCloud/glmdag")

rm(list=ls())
source("./codes/00_load_ftns.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")


g = "rand"

iter = 1
seed_para = 1
n_obs = 50 # the number of observations
n_conti = 3 # the number of continuous nodes
n_ordin = 3
n_multi = 4 # the number of multinomial nodes
n_nodes = n_conti + n_multi + n_ordin # the number of total nodes
types_by_node = gen_node_types(n_conti, n_multi, n_ordin, seed = seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
n_levels_by_node[1] <- 4

graph_set = gen_graph_adj_mat(n_nodes, g, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)

# Generate an original dataset
data_input_origin = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
data_input_nom = data_input_ord <- data_input_origin
