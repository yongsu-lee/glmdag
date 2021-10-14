## Load necessary packages

if(!require(igraph)){
  install.packages("igraph", quiet = T)
  library("igraph", quietly = T)
}

if(!require(bnlearn)){
  install.packages("bnlearn", quiet = T)
  library("bnlearn", quietly = T)
}


## Load commonly-used functions (short-length functions)
source("./codes/shared_fun.R")

## Cpp-enhanced functions to deal with bottlenecks
# sourceCpp("./codes/arma_trexpm.cpp")
# sourceCpp("./codes/obj0_mix.cpp")
# sourceCpp("./codes/grad0_mix.cpp")
# sourceCpp("./codes/lkhd_mix.cpp")
# sourceCpp("./codes/lkhd_multi.cpp")
# sourceCpp("./codes/lkhd_ordin.cpp")
# sourceCpp("./codes/lkhd_conti.cpp")

## Generating data source
source("./codes/gen_graph_adj_mat.R")
source("./codes/gen_para.R")
source("./codes/gen_data.R")

## Main algorithm
# source("./codes/glmdag.R")
# source("./codes/initialize_W.R")
# source("./codes/gen_lambdas.R")
# source("./codes/admm_loop.R")
# source("./codes/mix.W_update.R")
# source("./codes/Beta_update.R")

## Evaluation source
# source("./codes/dag2cpdag.R")
# source("./codes/push_dag.R")
# source("./codes/tun_sel_lkhd.R")
# source("./codes/eval_metrics.R")


## Load old function for comparison








