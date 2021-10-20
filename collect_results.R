## Setup parameters and containers ####
rm(list=ls())
require(tidyverse)

source("~/iCloud/glmdag/codes/00_load_ftns.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# graph_types_long = c("Bipartite","Random DAG","Scale-Free","Small-World","Tree")
# methods_long = c("glmdag-orginal", "glmdag-ordianl", "glmdag-nominal", 
#                  "hc", "mmhc")
graph_types = c("bi","rand","sf","sw","tree")
# graph_types = "rand"
methods =  c("mc","mo", "mm")
thre_list = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
eps_list = c("005", "01", "02", "03", "04", "05")
# eps_list = c("005", "01", "02", "03", "04")

# Setup global table
global_table <- data.frame(
  P = double(), E = double(), TPR = double(), JI = double(), 
  R = double(), M = double(), FP = double(), FDP = double(), SHD = double(),
  three = double(),  lam_eps = character(), 
  method = character(), graph = character())

# Set error flag 
err_table = data.frame(file_name = character(), eps = character())


## GLMDAG results ####
for (eps in eps_list){ # eps = "005"
  
  init_dir = "~/iCloud/glmdag_result/small_eps_"
  setwd(paste0(init_dir, eps))
  
  for (graph_type in graph_types){ # graph_type = "rand"
    
    simu_info = readRDS("simu_info.rds")
    list2env(simu_info, globalenv())
    
    # Setup evaluation table 
    empty_eval_table <- array(NA, dim = c(7, 9, n_iter))
    
    n_nodes = n_conti + n_multi + n_ordin
    types_by_node = gen_node_types(n_conti + n_ordin, n_multi, 0, seed = seed_para)
    n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
    
    graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
    graph_true = graph_set$graph_true
    A_true = graph_set$A_true
    n_edge = sum(A_true) 
    W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)
    
    for (m in methods){ # m = "hc"
      
      m_table <- empty_eval_table
      
      for (iter in 1:n_iter){ # iter = 1
        
        cat(paste0("******* iter = ", iter, "********\n"))
        init_file_name = paste0("simu2_small_n_50_", graph_type, "_", 
                                m, "_", iter, ".rds")
        
        result = readRDS(init_file_name)
        A_est = result$pushed_A_est_by_lam 
        sel_lam = result$sel_tun
        
        if (is.na(A_est[1]) == T) { 
          err_table_current = data.frame(file_name = init_file_name, eps = eps)
          err_table = rbind(err_table_current, err_table)
          next
        }
        
        temp = eval_table(A_est, A_true, "obs", sel_lam)[,-c(10,11)]
        row.names(temp) <- NULL
        
        m_table[, , iter] <- as.matrix(temp)
        
      } # end of iter
      
      mean_table = as.data.frame(apply(m_table, c(1,2), 
                                       function(x) mean(x, na.rm =T)))
      colnames(mean_table) <- c("P", "E", "TPR", "JI", "R", "M", "FP", "FDP", "SHD")
      mean_table2 = data.frame(mean_table, thre = thre_list, lam_eps = eps, 
                               method = m, graph = graph_type)
      
      global_table = rbind(global_table, mean_table2)
      
    }
    
  }
  
}

saveRDS(global_table, "~/iCloud/glmdag_result/global_table.rds")

err_table


## HC case ####
for (eps in eps_list){ # eps = "005"
  
  init_dir = "~/iCloud/glmdag_result/small_eps_"
  setwd(paste0(init_dir, eps))
  
  for (graph_type in graph_types){ # graph_type = "rand"
    
    simu_info = readRDS("simu_info.rds")
    list2env(simu_info, globalenv())
    
    # Setup evaluation table 
    empty_eval_table <- array(NA, dim = c(7, 9, n_iter))
    
    n_nodes = n_conti + n_multi + n_ordin
    types_by_node = gen_node_types(n_conti + n_ordin, n_multi, 0, seed = seed_para)
    n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
    
    graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
    graph_true = graph_set$graph_true
    A_true = graph_set$A_true
    n_edge = sum(A_true) 
    W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)
    
    for (m in methods){ # m = "hc"
      
      m_table <- empty_eval_table
      
      for (iter in 1:n_iter){ # iter = 1
        
        cat(paste0("******* iter = ", iter, "********\n"))
        init_file_name = paste0("simu2_small_n_50_", graph_type, "_", 
                                m, "_", iter, ".rds")
        
        result = readRDS(init_file_name)
        A_est = result$pushed_A_est_by_lam 
        sel_lam = result$sel_tun
        
        if (is.na(A_est[1]) == T) { 
          err_table = rbind(init_file_name, err_table)
          next
        }
        
        temp = eval_table(A_est, A_true, "obs", sel_lam)[,-c(10,11)]
        row.names(temp) <- NULL
        
        m_table[, , iter] <- as.matrix(temp)
        
      } # end of iter
      
      mean_table = as.data.frame(apply(m_table, c(1,2), 
                                       function(x) mean(x, na.rm =T)))
      colnames(mean_table) <- c("P", "E", "TPR", "JI", "R", "M", "FP", "FDP", "SHD")
      mean_table2 = data.frame(mean_table, thre = thre_list, lam_eps = eps, 
                               method = m, graph = graph_type)
      
      global_table = rbind(global_table, mean_table2)
      
    }
    
  }
  
}



## MMHC case ####






## Find the error table ####
err_table = data.frame(file_name = character(), eps = character())

for (eps in eps_list){ # eps = "005"
  
  init_dir = "~/iCloud/glmdag_result/small_eps_"
  setwd(paste0(init_dir, eps))
  
  for (graph_type in graph_types){ # graph_type = "rand"
    
    simu_info = readRDS("simu_info.rds")
    list2env(simu_info, globalenv())

    for (m in methods){ # m = "hc"
      
      # m_table <- empty_eval_table
      
      for (iter in 1:n_iter){ # iter = 1
        
        cat(paste0("******* iter = ", iter, "********\n"))
        init_file_name = paste0("simu2_small_n_50_", graph_type, "_", 
                                m, "_", iter, ".rds")
        
        result = readRDS(init_file_name)
        A_est = result$pushed_A_est_by_lam 
        
        if (is.na(A_est[1]) == T) { 
          err_table_current = data.frame(file_name = init_file_name, eps = eps)
          err_table = rbind(err_table_current, err_table)
          next
        }
      } 
    }
  }
}

saveRDS(err_table,"~/iCloud/glmdag_result/err_report.rds")


## Debugging Block ####

# ## Simulation Setup +++++++++++++++++++++++++++
# seed_para = 1
# n_iter = 100
# n_conti = 3
# n_ordin = 3
# n_multi = 4
# n_ordin = 0
# n_nodes = n_conti + n_multi + n_ordin
# types_by_node = gen_node_types(n_conti + n_ordin, n_multi, 0, seed = seed_para)
# n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
# # +++++++++++++++++++++++++++++++++++++++++++++


## Simulation Setup +++++++++++++++++++++++++++
# seed_para = 1
# n_iter = 100
# n_conti = 6
# n_multi = 4
# n_ordin = 0
# n_nodes = n_conti + n_multi + n_ordin
# types_by_node = gen_node_types(n_conti + n_ordin, n_multi, 0, seed = seed_para)
# n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
# +++++++++++++++++++++++++++++++++++++++++++++
