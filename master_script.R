## Load queue arguments ####
rm(list=ls())
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

## Set queue number accoring to system ####
sysname = Sys.info()['sysname']
if (sysname == "Linux"){ # for CHTC server
  
  simu_info = readRDS("./simu_info.rds")
  list2env(simu_info, globalenv())
  if (path_par == F) ell = NULL
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Check this block carefully! ++++++++++++++++++++++++++++++++++++++++++++++++
  
  #!/usr/bin/env Rscript
  queue_args = commandArgs(trailingOnly=TRUE)
  graph_type = as.character(queue_args[1])
  method = as.character(queue_args[2])
  iter = as.integer(queue_args[3])
  ell = as.integer(queue_args[4])
  
  name_temp = paste(simu_case, size,"n", n_obs, 
                    graph_type, method, iter, ell, sep = "_")
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Pull glmdag codes from github
  url = "https://github.com/yongsu-lee/glmdag/archive/refs/heads/master.zip"
  download.file(url = url, destfile = "glmdag.zip")
  Sys.sleep(10)
  unzip("glmdag.zip")
  Sys.sleep(10)
  init_dir = "./glmdag-master/"
  
} else { # for local macOS
  
  setwd("~/iCloud/glmdag/")
  init_dir = "./"

  file.edit("~/iCloud/glmdag/gen_input.R")
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Check this block carefully! ++++++++++++++++++++++++++++++++++++++++++++++++
  # queue = read.table("queue_list", sep =",", strip.white = T)
  # queue_args = queue[1, ]
  graph_type = "rand"
  method = "mo"
  iter = 1
  if (path_par == F) ell = NULL
  name_temp = paste(simu_case, size,"n", n_obs, 
                    graph_type, method, iter, ell, sep = "_")
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
}

## Sourcing files
Sys.sleep(10)
source(paste0(init_dir,"codes/00_load_ftns.R"))

## Generate the true graph and obtain corresponding adjacency matrix
n_nodes = n_conti + n_multi + n_ordin
graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
graph_true = graph_set$graph_true
A_true = graph_set$A_true
n_edge = sum(A_true) 

## Run simulation
source(paste0(init_dir,simu_case,"_script.R"))

## Save the result
saveRDS(result, file = paste0(name_temp,".rds"))

