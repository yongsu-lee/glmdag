## Delete this region when you submit +++++++++++++++++
setwd("~/iCloud/paper_a/paper_a_submit/biom_j/biom_j_code_submit/")
use_saved_data = 1
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Clear the memory
rm(list=ls())

# Load required packages and functions
source("./codes/00_load_ftns.R")

# Set seed generating system
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Download saved result to regenerate the results
use_saved_data = 1
if (use_saved_data == 1){
  #.. Download dataset from github and saved under "./results/"
  #.. There will be five folders
  #.. simu1_n_50, simu1_n_200, simu2_small, simu2_large, hcc
}

graph_types_long = c("Bipartite", "Random DAG", "Scale-Free", 
                     "Small-World", "Tree")

graph_types = c("bi","rand","sf","sw","tree")


## Large Case +++++++++++++++++


# Set some parameters
methods =  c("mc","mo", "mm")
n_iter = 100
seed_para = 1
n_obs = 200 # the number of observations
n_conti = 35 # the number of continuous nodes
n_multi = 15 # the number of multinomial nodes
n_nodes = n_conti + n_multi # the number of total nodes

types_by_node = gen_node_types(n_conti, n_multi, 0, seed = seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)

# Setup table 

write(x=NULL, file ="./results/Table_4.txt", append = F)
empty_eval_table <- data.frame(P = double(), E = double(), TPR = double(),
                               JI = double(), R = double(), M = double(), 
                               FP = double(), FDP = double(), SHD = double())

# Simulations starts and generates Tabl3 3
#.. We parallelizing following nested loop in the original simulation.
#.. If you run this loop using a single core, it will take a long time.
#.. As a default, we are going to use saved results. 
#.. If you want to run simulations directly, set 'use_saved_data = 0'
for (g in graph_types){ # g = "bi"
  
  graph_set = gen_graph_adj_mat(n_nodes, g, seed = seed_para)
  graph_true = graph_set$graph_true
  A_true = graph_set$A_true
  W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)
  
  mc_table = mo_table = mm_table = hc_table = mmhc_table <- empty_eval_table
  
  for (iter in 1:n_iter){ # iter = 1
    
    # Generate an original dataset
    data_input_origin = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
    data_input_nom = data_input_ord <- data_input_origin
    
    # Select will-be discretizing nodes (20 nodes)
    disc_nodes = sample(which(types_by_node == "c"), 20, replace = F)
    n_disc = length(disc_nodes)
    conti_raw_nodes = which(types_by_node == "c")
    conti_nodes = setdiff(conti_raw_nodes, disc_nodes) 
    two_level_nodes = disc_nodes[1:(n_disc/2)]
    four_level_nodes = disc_nodes[(n_disc/2+1):n_disc]
    
    # Discretizing into 2 levels
    temp = data_input_origin[two_level_nodes]
    temp1 = lapply(temp, function(x) 
      cut(x, c(min(x)-1, median(x), max(x)), labels = 1:2))
    
    # Discretizing into 4 levels
    temp = data_input_origin[four_level_nodes]
    temp2 = lapply(temp, function(x) 
      cut(x, c(min(x)-1, quantile(as.matrix(x), probs = c(0.25, 0.5, 0.75, 1))), 
          labels = 1:4))
    
    discretized_data = data.frame(temp1, temp2)
    
    # Replace continuous nodes with discretized nominal level
    data_input_nom[disc_nodes] <- discretized_data
    
    # Replace continuous nodes with discretized ordinal level
    ordered_data = lapply(discretized_data, function(x) factor(x,ordered = T))
    data_input_ord[disc_nodes] <- ordered_data
    
    for (m in methods){ # m = "mo"
      
      data_input = switch(m, "mc" = data_input_origin, "mo" = data_input_ord,
                          "mm" = data_input_nom)
      
      if (use_saved_data != 1){
        
        result = glmdag(data_input, n_lams = n_lams, verbose = T,
                        admm_args = admm_arg_ctrl(warm_start = F))
        
        if (m == "mm") { 
          result.hc <- hc(data_input)
          result.mmhc <- mmhc(data_input)
        }
        
      } else { # when using stored results
        
        temp = paste("./results/simu2_large/simu2_large",g, m, iter, sep="_")
        result = readRDS(file = paste0(temp, ".rds"))
        
        if (m == "mm") {
          result.temp = readRDS(file = paste0(temp,"_mmhc.rds"))
          result.hc <- result.temp$hc
          result.mmhc <-  result.temp$mmhc
        }
        
      } # end of use_saved_data
      # end of methods
      
      A_est = result$pushed_A_est_by_lam
      sel_lam = result$sel_tun
      
    
      Beta_est = result$Beta_est_by_lam
      
      sel_tun = tun_sel_lkhd(A_est, Beta_est, result, thre_list = thre_list)
      
      
      sel_lam = result$sel_tun
      
      if (m == "mc") {
        mc_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[5,-c(10,11)]
      } else if (m == "mo"){
        mo_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[5,-c(10,11)]
      } else {
        mm_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[5,-c(10,11)]
        hc_table[iter,] <- eval_table(amat(result.hc), A_true, "obs", single = T)
        mmhc_table[iter,] <- eval_table(amat(result.mmhc),A_true,"obs",single=T)
      }
    }
  } # end of iter
  
  graph_name = graph_types_long[which(graph_types == g)]
  
  write(paste(graph_name, ": n_edges = ", sum(A_true), "\n"),
        file ="./results/Table_4.txt", append = T)
  
  write.table(rbind(c("","","","P","E","TPR","JI","R*","M","FP","FDP","SHD")),
              file ="./results/Table_4.txt", append = T, sep ="\t",
              col.names = F, row.names = F, quote = F)
  
  write(paste0("--------------------------------------------------------------",
               "-----------------------------"),  
        file ="./results/Table_4.txt", append = T)
  
  ave_mc_table = as.data.frame(rbind(colMeans(mc_table) %>% round(digits = 2)))
  print_mc_table = data.frame(m = "GLMDAG: Original ", ave_mc_table )
  write.table(print_mc_table, file ="./results/Table_4.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mo_table = as.data.frame(rbind(colMeans(mo_table) %>% round(digits = 2)))
  print_mo_table = data.frame(m = "GLMDAG: Ordinal  ", ave_mo_table )
  write.table(print_mo_table, file ="./results/Table_4.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mm_table = as.data.frame(rbind(colMeans(mm_table) %>% round(digits = 2)))
  print_mm_table = data.frame(m = "GLMDAG: Nominal  ", ave_mm_table )
  write.table(print_mm_table, file ="./results/Table_4.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  write(paste0("--------------------------------------------------------------",
               "-----------------------------"),  
        file ="./results/Table_4.txt", append = T)
  
  ave_hc_table = as.data.frame(rbind(colMeans(hc_table) %>% round(digits = 2)))
  print_hc_table = data.frame(m = "HC               ", ave_hc_table )
  write.table(print_hc_table, file ="./results/Table_4.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mmhc_table = as.data.frame(rbind(colMeans(mmhc_table) %>% round(digits = 2)))
  print_mmhc_table = data.frame(m = "MMHC             ", ave_mmhc_table )
  write.table(print_mmhc_table, file ="./results/Table_4.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  write("\n\n",  file ="./results/Table_4.txt", append = T)
  
  print(paste0("Graph: ", graph_name, " has been generated! ..."))
  
  
} # end of g in graph_types

# .. Check "./results/Table_4.txt". 





## Small Case +++++++++++++++++


# Set some parameters
methods =  c("mc","mo", "mm")
n_iter = 100
seed_para = 1
n_obs = 50 # the number of observations
n_conti = 6 # the number of continuous nodes
n_multi = 4 # the number of multinomial nodes
n_nodes = n_conti + n_multi # the number of total nodes
types_by_node = gen_node_types(n_conti, n_multi, 0, seed = seed_para)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)

# Setup table 

write(x=NULL, file ="./results/Table_3_thre.txt", append = F)
empty_eval_table <- data.frame(P = double(), E = double(), TPR = double(),
                               JI = double(), R = double(), M = double(), 
                               FP = double(), FDP = double(), SHD = double())

# Simulations starts and generates Tabl3 3
#.. We parallelizing following nested loop in the original simulation.
#.. If you run this loop using a single core, it will take a long time.
#.. As a default, we are going to use saved results. 
#.. If you want to run simulations directly, set 'use_saved_data = 0'
for (g in graph_types){ # g = "bi"
  
  graph_set = gen_graph_adj_mat(n_nodes, g, seed = seed_para)
  graph_true = graph_set$graph_true
  A_true = graph_set$A_true
  W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)
  
  mc_table = mo_table = mm_table = hc_table = mmhc_table <- empty_eval_table
  
  for (iter in 1:n_iter){ # iter = 1
    
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
    
    for (m in methods){ # m = "mc"
      
      data_input = switch(m, "mc" = data_input_origin, "mo" = data_input_ord,
                          "mm" = data_input_nom)
      
      if (use_saved_data != 1){
        
        result = glmdag(data_input, n_lams = n_lams, verbose = T,
                        admm_args = admm_arg_ctrl(warm_start = F))
        
        if (m == "mm") { 
          result.hc <- hc(data_input)
          result.mmhc <- mmhc(data_input)
        }
        
      } else { # when using stored results
        
        temp = paste("./results/simu2_small/simu2_small",g, m, iter, sep="_")
        result = readRDS(file = paste0(temp, ".rds"))
        
        if (m == "mm") {
          result.temp = readRDS(file = paste0(temp,"_mmhc.rds"))
          result.hc <- result.temp$hc
          result.mmhc <-  result.temp$mmhc
        }
        
      } # end of use_saved_data
      # end of methods
      
      A_est = result$pushed_A_est_by_lam
      Beta_est = result$Beta_est_by_lam

      sel_tun = tun_sel_lkhd(A_est, Beta_est, result, thre_list = thre_list)

      sel_lam = result$sel_tun
      
      if (m == "mc") {
        mc_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[7,-c(10,11)]
      } else if (m == "mo"){
        mo_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[7,-c(10,11)]
      } else {
        mm_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[7,-c(10,11)]
        hc_table[iter,] <- eval_table(amat(result.hc), A_true, "obs", single = T)
        mmhc_table[iter,] <- eval_table(amat(result.mmhc),A_true,"obs",single=T)
      }
    }
  } # end of iter
  
  graph_name = graph_types_long[which(graph_types == g)]
  
  write(paste(graph_name, ": n_edges = ", sum(A_true), "\n"),
        file ="./results/Table_3_thre.txt", append = T)
  
  write.table(rbind(c("","","","P","E","TPR","JI","R*","M","FP","FDP","SHD")),
              file ="./results/Table_3_thre.txt", append = T, sep ="\t",
              col.names = F, row.names = F, quote = F)
  
  write(paste0("-----------------------------------------------------------",
               "----------------------------------------"), 
        file ="./results/Table_3_thre.txt", append = T)
  
  ave_mc_table = as.data.frame(rbind(colMeans(mc_table) %>% round(digits = 2)))
  print_mc_table = data.frame(m = "GLMDAG: Original ", ave_mc_table )
  write.table(print_mc_table, file ="./results/Table_3_thre.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mo_table = as.data.frame(rbind(colMeans(mo_table) %>% round(digits = 2)))
  print_mo_table = data.frame(m = "GLMDAG: Ordinal  ", ave_mo_table )
  write.table(print_mo_table, file ="./results/Table_3_thre.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mm_table = as.data.frame(rbind(colMeans(mm_table) %>% round(digits = 2)))
  print_mm_table = data.frame(m = "GLMDAG: Nominal  ", ave_mm_table )
  write.table(print_mm_table, file ="./results/Table_3_thre.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  write(paste0("-----------------------------------------------------------",
               "----------------------------------------"), 
        file ="./results/Table_3_thre.txt", append = T)
  
  ave_hc_table = as.data.frame(rbind(colMeans(hc_table) %>% round(digits = 2)))
  print_hc_table = data.frame(m = "HC               ", ave_hc_table )
  write.table(print_hc_table, file ="./results/Table_3_thre.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mmhc_table = as.data.frame(rbind(colMeans(mmhc_table) %>% round(digits = 2)))
  print_mmhc_table = data.frame(m = "MMHC             ", ave_mmhc_table )
  write.table(print_mmhc_table, file ="./results/Table_3_thre.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  write("\n\n",  file ="./results/Table_3_thre.txt", append = T)
  
  print(paste0("Graph: ", graph_name, " has been generated! ..."))
  
  
} # end of g in graph_types

# .. Check "./results/Table_3_thre.txt". 


