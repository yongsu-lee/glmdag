## Load queue arguments ####
rm(list=ls())
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

## Set queue number accoring to system ####
sysname = Sys.info()['sysname']
if (sysname == "Linux"){ # for CHTC server
  
  simu_info = read.table("./simu_info.txt", header = T)
  list2env(as.list(simu_info), environment())
  
  simu_info_common = read.table("./simu_info_common.txt", header = T)
  list2env(as.list(simu_info_common), environment())
  
  eps_lam = as.double(eps_lam)
  simu_script = ifelse(simu_case == "simu2_mmhc_only", "simu2", simu_case)
  
  if (path_par == F) ell = NULL
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Check this block carefully! ++++++++++++++++++++++++++++++++++++++++++++++++
  
  #!/usr/bin/env Rscript
  queue_args = commandArgs(trailingOnly=TRUE)
  graph_type = as.character(queue_args[1])
  method = as.character(queue_args[2])
  iter = as.integer(queue_args[3])
  # ell = as.integer(queue_args[2])
  
  name_temp = paste(simu_case, size,  
                    graph_type, iter, method, sep = "_")
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Pull glmdag codes from github
  url = "https://github.com/yongsu-lee/glmdag/archive/refs/heads/master.zip"
  download.file(url = url, destfile = "glmdag.zip")
  Sys.sleep(10)
  unzip("glmdag.zip")
  Sys.sleep(10)
  init_dir = "./glmdag-master/"
  
} else { # for local macOS ####
  
  simu_info = read.table("./simu_info.txt", header = T)
  list2env(as.list(simu_info), environment())
  
  simu_info_common = read.table("./simu_info_common.txt", header = T)
  list2env(as.list(simu_info_common), environment())
  
  eps_lam = as.double(eps_lam)
  simu_script = ifelse(simu_case == "simu2_mmhc_only", "simu2", simu_case)
  
  setwd("~/Dropbox/glmdag/")
  init_dir = "./"
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Check this block carefully! ++++++++++++++++++++++++++++++++++++++++++++++++
  # queue = read.table("queue_list", sep =",", strip.white = T)
  # queue_args = queue[1, ]
  graph_type = "rand"
  method ="mm"
  iter = 1
  # ell = 15
  if (path_par == F) ell = NULL
  name_temp = paste(simu_case, size,  
                    graph_type, iter, method, sep = "_")
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
}

## Sourcing files ####
Sys.sleep(10)
source(paste0(init_dir,"codes/00_load_ftns.R"))

## Generate the true graph and obtain corresponding adjacency matrix
n_nodes = n_conti + n_multi + n_ordin
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true
n_edge = sum(A_true) 

## Run simulation
source(paste0(init_dir, simu_script, "_script.R"))

file.remove("glmdag.zip")