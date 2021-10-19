rm(list=ls())
require(tidyverse)

source("~/iCloud/glmdag/codes/00_load_ftns.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# graph_types_long = c("Bipartite","Random DAG","Scale-Free","Small-World","Tree")
# graph_types = c("bi","rand","sf","sw","tree")
graph_types = "rand"
methods =  c("mc","mo", "mm")
eps_list = c("005", "01", "02", "03", "04", "05")
# thre_list =  c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
eps_list = c("005", "01", "02", "03", "04")


## Simulation Setup +++++++++++++++++++++++++++
seed_para = 1
n_iter = 100
n_conti = 6
n_multi = 4
n_ordin = 0
n_nodes = n_conti + n_multi + n_ordin
types_by_node = gen_node_types(n_conti + n_ordin, n_multi, 0, seed = seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)
# +++++++++++++++++++++++++++++++++++++++++++++


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

# Setup global table
global_table <- data.frame(
  P = double(), E = double(), TPR = double(), JI = double(), 
  R = double(), M = double(), FP = double(), FDP = double(), SHD = double(),
  three = double(),  lam_eps = character(), 
  method = character(), graph = character())


# Setup evaluation table 
empty_eval_table <- array(NA, dim = c(7, 9, n_iter))

# empty_eval_table <- data.frame(
#   P = double(), E = double(), TPR = double(), JI = double(), 
#   R = double(), M = double(), FP = double(), FDP = double(), SHD = double())

for (eps in eps_list){ # eps = "005"
  
  init_dir = "~/iCloud/glmdag_result/small_eps_"
  setwd(paste0(init_dir, eps))
  
  for (graph_type in graph_types){ # graph_type = "rand"
    
    graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
    graph_true = graph_set$graph_true
    A_true = graph_set$A_true
    n_edge = sum(A_true) 
    W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)
    
    for (m in methods){ # m = "mm"
      
      m_table = hc_table = mmhc_table <- empty_eval_table
      
      for (iter in 1:n_iter){ # iter = 1
        
        cat(paste0("******* iter = ", iter, "\n"))
        init_file_name = paste0("simu2_small_n_50_", graph_type, "_", 
                                m, "_", iter, ".rds")
        
        result = readRDS(init_file_name)
        A_est = result$pushed_A_est_by_lam 
        sel_lam = result$sel_tun
        
        if (is.na(A_est[1]) == T) next
        
        temp = eval_table(A_est, A_true, "obs", sel_lam)[,-c(10,11)]
        row.names(temp) <- NULL
        
        m_table[, , iter] <- as.matrix(temp)
        
        # if (m %in% c("hc", "mmhc")){ # m = "hc"
        #   
        #   mm = "mmhc"
        #   mmhc_init_file_name = "~/iCloud/glmdag_result/small_mmhc/mmhc_small_n_50_"
        # 
        #   mmhc_file_name = paste0(mmhc_init_file_name, graph_type, "_",
        #                           m, "_", iter, ".rds")
        #   
        #   result = readRDS(mmhc_file_name)
        #   
        #   result.hc <- result$hc
        #   result.mmhc <-  result$mmhc
        # 
        #   temp.hc = eval_table(amat(result.hc), A_true, "obs", single = T)
        #   temp.mmhc = eval_table(amat(result.mmhc), A_true,"obs", single=T)
        # 
        #   temp2.hc = data.frame(temp.hc, thre = thre_list, eps = eps,
        #                         method = m, graph_type = graph_type)
        #   temp2.mmhc = data.frame(temp.mmhc, thre = thre_list, eps = eps,
        #                           method = m, graph_type = graph_type)
        # 
        #   hc_table = rbind(hc_table, temp2.hc)
        #   mmhc_table = rbind(mmhc_table, temp2.mmhc)
        # 
        # }
        
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

