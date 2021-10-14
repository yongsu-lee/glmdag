## Set subdirectory name according to the system
subdir = ifelse(Sys.info()["sysname"] == "Linux", "./", "./codes/")
names(subdir) <- NULL

## Load necessary packages
source(paste0(subdir, "load_r_pkgs.R"))

## Load commonly-used functions (short-length functions)
source(paste0(subdir, "shared_fun.R"))

## Cpp-enhanced functions to deal with bottlenecks
sourceCpp(paste0(subdir, "arma_trexpm.cpp"))
sourceCpp(paste0(subdir, "obj0_mix.cpp"))
sourceCpp(paste0(subdir, "grad0_mix.cpp"))
sourceCpp(paste0(subdir, "lkhd_mix.cpp"))
sourceCpp(paste0(subdir, "lkhd_multi.cpp"))
sourceCpp(paste0(subdir, "lkhd_ordin.cpp"))
sourceCpp(paste0(subdir, "lkhd_conti.cpp"))

## Generating data source
source(paste0(subdir, "gen_graph_adj_mat.R"))
source(paste0(subdir, "gen_para.R"))
source(paste0(subdir, "gen_data.R"))

## Main algorithm
source(paste0(subdir, "glmdag.R"))
source(paste0(subdir, "initialize_W.R"))
source(paste0(subdir, "gen_lambdas.R"))
source(paste0(subdir, "admm_loop.R"))
source(paste0(subdir, "mix.W_update.R"))
source(paste0(subdir, "Beta_update.R"))

## Evaluation source
source(paste0(subdir, "dag2cpdag.R"))
source(paste0(subdir, "push_dag.R"))
source(paste0(subdir, "tun_sel_lkhd.R"))
# source(paste0(subdir, "eval_metrics.R"))


## Load old function for comparison








