if(!require(bnlearn)){
  install.packages("bnlearn", quiet = T)
  library("bnlearn", quietly = T)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz")

require(Rgraphviz)

seed_para = 1
seed_iter = 1
graph_type = "rand"
n_conti = 3
n_multi = 3
n_ordin = 4
n_obs = 50
n_nodes = n_conti + n_multi + n_ordin
types_by_node = gen_node_types(n_conti, n_multi, n_ordin, seed = seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)

graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)
n_edges = sum(A_true) 

data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = seed_iter)


result.hc = hc(data_input)
graphviz.plot(result.hc)
A_est_hc = amat(result.hc)

result.mmhc = mmhc(data_input)
graphviz.plot(result.mmhc)
A_est_mmhc = amat(result.mmhc)

str(data_input)
data_input_convert = data_input
ord_nodes = sapply(data_input, function(x) any( class(x) == "ordered"))
data_input_convert[, ord_nodes] <- 
  lapply(data_input[, ord_nodes], function(x) factor(x, ordered = F)) 

str(data_input_convert)

result.hc_convert = hc(data_input_convert)
graphviz.plot(result.hc_convert)
A_est_hc_convert = amat(result.hc_convert)

result.mmhc_convert = mmhc(data_input_convert)
graphviz.plot(result.mmhc_convert)
A_est_mmhc_convert = amat(result.mmhc_convert)
