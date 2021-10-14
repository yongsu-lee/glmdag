setwd("~/iCloud/paper_a/paper_a_submit/biom_j/biom_j_code_submit/")

rm(list=ls())
source("./codes/00_load_ftns.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

eps_lambda = 0.4
g = "bi"
m = "mc"
iter = 67

n_lams = 30

seed_para = 1
n_obs = 50 # the number of observations
n_conti = 6 # the number of continuous nodes
n_multi = 4 # the number of multinomial nodes
n_nodes = n_conti + n_multi # the number of total nodes
types_by_node = gen_node_types(n_conti, n_multi, 0, seed = seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)


graph_set = gen_graph_adj_mat(n_nodes, g, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true
W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)

# Generate an original dataset
data_input_origin = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
data_input_nom = data_input_ord <- data_input_origin

# Ready for discretizing
conti_nodes = which(types_by_node == "c")
two_level_nodes = conti_nodes[1:(n_conti/2)]
four_level_nodes = conti_nodes[(n_conti/2+1):n_conti]

# Discretizing into 2 levels
temp = data_input_origin[two_level_nodes]
temp1 = lapply(temp, function(x) 
  cut(x, c(min(x)-1, median(x), max(x)), labels = 1:2))

# Discretizing into 4 levels
temp = data_input_origin[four_level_nodes]
temp2 = lapply(temp, function(x) 
  cut(x, c(min(x)-1, quantile(as.matrix(x), probs=c(0.25,0.5,0.75,1))), 
      labels = 1:4))

discretized_data = data.frame(temp1, temp2)

# Replace continuous nodes with discretized nominal level
data_input_nom[conti_nodes] <- discretized_data

# Replace continuous nodes with discretized ordinal level
ordered_data = lapply(discretized_data, function(x) factor(x,ordered = T))
data_input_ord[conti_nodes] <- ordered_data

data_input = switch(m, "mc" = data_input_origin, "mo" = data_input_ord,
                    "mm" = data_input_nom)

Sys.sleep(1)
result = glmdag(data_input, n_lams = n_lams, verbose = T, 
                eps_lambda = eps_lambda, fit_hist = T)


file_name = paste0("simu2_small_",g,"_",m,"_",iter,".rds")
saveRDS(result, file_name)


