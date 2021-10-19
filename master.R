# Last Updated: Sep 23, 2021
# Set working direction to "codes_biomj_yongsu" folder
# Outputs will be generated in the "results" folder

## Delete this region when you submit +++++++++++++++++
setwd("~/iCloud/glmdag")
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Figure 2 (Graphs, CD vs GLMDAG) ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Set some parameters
seed_para = 10
n_nodes = 30
# graph_types = c("bi", "rand", "sf", "sw","tree")
# graph_names = c("Bipartite", "Random DAG", "Scale-Free", "Small-World", "Tree")
types_by_node = gen_node_types(0, n_multi = n_nodes, 0, seed = seed_para)

# Set plot layout 
m = c(1,1,2,2,3,3,1,1,2,2,3,3,0,4,4,5,5,0,0,4,4,5,5,0)
M = matrix(m, nrow =4, ncol = 6, byrow = T); layout(mat = M)

# Draw graphs
pdf("./results/Figure_2.pdf", width = 10)
layout(mat = M); par(mar=c(0,0.5,3,0.5)); par(oma=c(0,0,0,0))
for (j in 1:5){
  
  graph_type = graph_types[j]
  graph_name = graph_types_long[j]
  
  graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
  graph_true = graph_set$graph_true
  A_true = graph_set$A_true
  
  fancy_dag_plot(A_true, types_by_node, title = graph_name, legend = F, 
                 edge.arrow.size = 0.3,
                 vertex.size = 10,
                 layout = layout_set(graph_true, graph_type))
}
dev.off()
# .. Check "./results/Figure_2.pdf".


#++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Table 2 (CD vs GLMDAG results) ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

# graph_types = c("bi","rand")
# n_iter = 5

# Set some parameters
n_iter = 100
seed_para = 10
n_obs_list = c(50, 200)
n_nodes = 30 
types_by_node = rep("m", n_nodes)
n_levels_by_node = gen_node_levels(types_by_node, 4, seed = seed_para)

# Setup table 
gg = 0 
write(x = NULL, file ="./results/Table_2.txt", append = F)
empty_eval_table <- data.frame(P = double(), E = double(), TPR = double(),
                               JI = double(), R = double(), M = double(), 
                               FP = double(), FDP = double(), SHD = double())

# Simulations starts and generates Tabel 2
#.. We parallelizing following nested loop in the original simulation.
#.. If you run this loop using a single core, it will take a long time.
#.. As a default, we are going to use saved results. 
#.. If you want to run simulations directly, set 'use_saved_data = 0'
for (g in graph_types){ # g = "bi"
  
  gg = gg + 1
  
  graph_set = gen_graph_adj_mat(n_nodes, g, seed = seed_para)
  graph_true = graph_set$graph_true
  A_true = graph_set$A_true
  W_true = gen_para(A_true, types_by_node, n_levels_by_node, seed = seed_para)
  n_edges = sum(A_true) 
  
  eval_table_glmdag <- empty_eval_table
  eval_table_cd <- empty_eval_table
  
  for (n_obs in n_obs_list){ # n_obs = 50
    
    iter_table_glmdag <- empty_eval_table
    iter_table_cd <- empty_eval_table
    
    ## !! Check old results (should be deleted) !!!
    # table.num = ifelse(n_obs == 50, "table3", "table4")
    
    for (iter in 1:n_iter){ # iter = 1
      
      # Generate a dataset
      data_input = gen_data(n_obs, A_true, graph_true, W_true, seed = iter)
      
      # Convert dataset for cd algorithm
      
      suppressMessages(data_input_cd <- conv_to_cd_data(data_input)) 
      
      if (use_saved_data != 1){
        
        result_glmdag = glmdag(data_input, n_lams = n_lams, verbose = T)
        result_cd = cd.run(indata = data_input_cd)
        
      } else {
        
        ## Original
        subdir = paste0("./results/simu1_n_",n_obs,"/simu1_n_",n_obs)
        temp = paste(subdir,g,iter, sep="_")
        result = readRDS(file = paste0(temp, ".rds"))
        result_glmdag = result$glmdag
        result_cd = result$cd
        
        ## !!! Test previous results !!!
        # subdir = paste0("./results/simu1_n_",n_obs,"_old/")
        # glmdag.name = paste0(subdir, "noteargis_", g, "_data_", iter, "_", 
        #                      table.num, ".rds")
        # cd.name = paste0(subdir, "cd_", g, "_data_", iter, "_", table.num, ".rds")
        # result_glmdag = readRDS(glmdag.name)
        # result_cd = readRDS(cd.name)
        
      }
      
      ## Origianl 
      temp_eval_glmdag = eval_by_lam(result_glmdag$pushed_A_est_by_lam,
                                     A_true, n_edges, type = "obs")
      
      ## !!! Test previous resutls !!!
      # temp_eval_glmdag = eval_by_lam(result_glmdag$A_est_by_lam,
      #                                A_true, n_edges, type = "obs")
      
      temp_sel_glmdag = mod_sel(temp_eval_glmdag, n_edges, sel_crit = "JI")
      iter_table_glmdag = rbind(iter_table_glmdag, temp_sel_glmdag)
      
      A_est_by_lam_cd = calc_A_est_by_lam_cd(result_cd)
      temp_eval_cd = eval_by_lam(A_est_by_lam_cd, A_true, n_edges, type = "obs")
      temp_sel_cd = mod_sel(temp_eval_cd, n_edges, sel_crit = "JI")
      iter_table_cd = rbind(iter_table_cd, temp_sel_cd)
      
    } # end of iter in 1:n_iter
    
    eval_table_glmdag = rbind(eval_table_glmdag, 
                              colMeans(iter_table_glmdag) 
                              %>% round(digits = 2))
    eval_table_cd = rbind(eval_table_cd, 
                          colMeans(iter_table_cd) %>% round(digits = 2))
    
  } # end of n_obs in n_obs_list
  
  # c("","P","E","TPR","JI","R*","M","FP","FDP","SHD")
  
  graph_name = graph_types_long[which(graph_types == g)]
  
  write(paste(graph_name, ": n_edges = ", sum(A_true), "\n"),
        file ="./results/Table_2.txt", append = T)
  
  write.table(rbind(c("","","","P","E","TPR","JI","R*","M","FP","FDP","SHD")),
              file ="./results/Table_2.txt", append = T, sep ="\t",
              col.names = F, row.names = F, quote = F)
  
  write(paste0("-----------------------------------------------------------",
               "------------------------------"), 
        file ="./results/Table_2.txt", append = T)
  
  eval_table_glmdag = data.frame(m = c("GLMDAG: n =  50 ", "        n = 200 "),
                                 eval_table_glmdag)
  write.table(eval_table_glmdag, file ="./results/Table_2.txt", sep ="\t",
              quote = F, append = T, row.names = F, col.names =F)
  
  write(paste0("-----------------------------------------------------------",
               "------------------------------"), 
        file ="./results/Table_2.txt", append = T)
  
  eval_table_cd = data.frame(m = c("    CD: n =  50 ", "        n = 200 "),
                             eval_table_cd)
  write.table(eval_table_cd, file ="./results/Table_2.txt", sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  write("\n\n",  file ="./results/Table_2.txt", append = T)
  
  print(paste0("Graph: ", graph_name, " has been generated! ..."))
  
  
} # end of g in graph_types

# .. Check "./results/Table_2.txt". 



#++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Figure 3 (Small Graphs, Nom vs Ord) ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Set some parameters
seed_para = 1
n_conti = 6
n_multi = 4
n_nodes = n_conti + n_multi 
graph_types = c("bi", "rand", "sf", "sw", "tree")
graph_names = c("Bipartite", "Random DAG", "Scale-Free", "Small-World", "Tree")
types_by_node = gen_node_types(n_conti, n_multi, 0, seed = seed_para)

# Set plot layout 
m = c(1,1,2,2,3,3,1,1,2,2,3,3,0,4,4,5,5,0,0,4,4,5,5,0)
M = matrix(m, nrow =4, ncol = 6, byrow = T); layout(mat = M)

pdf("./results/Figure_3.pdf", width = 10)
layout(mat = M); par(mar=c(0,0.5,3,0.5)); par(oma=c(0,0,0,0))

# Draw graphs
for (j in 1:5){
  
  graph_type = graph_types[j]
  graph_name = graph_names[j]
  
  graph_set = gen_graph_adj_mat(n_nodes, graph_type, seed = seed_para)
  graph_true = graph_set$graph_true
  A_true = graph_set$A_true
  
  fancy_dag_plot(A_true, types_by_node, title = graph_name, legend = F, 
                 edge.arrow.size = 0.5, 
                 layout = layout_set(graph_true, graph_type))
}
dev.off()
# .. Check "./results/Figure_3.pdf".


#++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Table 3 (Small Graphs, Nom vs Ord, Results) ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

## !!! Debugging Block !!! ++++++++++++++++++++
graph_types = c("bi", "rand","sf","sw","tree")
# +++++++++++++++++++++++++++++++++++++++++++++

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
thre_idx = 7
table_file_name = paste0("./results/Table_3_eps_03_using_",thre_idx,".txt")
write(x=NULL, file = table_file_name, append = F)
empty_eval_table <- data.frame(P = double(), E = double(), TPR = double(),
                               JI = double(), R = double(), M = double(), 
                               FP = double(), FDP = double(), SHD = double())

# Simulations starts and generates Tabl3 3
#.. We parallelizing following nested loop in the original simulation.
#.. If you run this loop using a single core, it will take a long time.
#.. As a default, we are going to use saved results. 
#.. If you want to run simulations directly, set 'use_saved_data = 0'
for (g in graph_types){ # g = "rand"
  
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
    
    for (m in methods){ # m = "mm"
      
      data_input = switch(m, "mc" = data_input_origin, "mo" = data_input_ord,
                          "mm" = data_input_nom)
      
      if (use_saved_data != 1){
        
        result = glmdag(data_input, n_lams = n_lams, verbose = T)
        
        if (m == "mm") { 
          result.hc <- bnlearn::hc(data_input)
          result.mmhc <- bnlearn::mmhc(data_input)
        }
        
      } else { # when using stored results
        
        temp = paste("./results/simu2_small/simu2_small",g, m, iter, sep="_")
        result = readRDS(file = paste0(temp, ".rds"))
        
        if (m == "mm") {
          temp2 = paste("./results/simu2_small_mmhc/simu2_small",
                        g, m, iter, sep="_")
          result.temp = readRDS(file = paste0(temp2,"_mmhc.rds"))
          result.hc <- result.temp$hc
          result.mmhc <-  result.temp$mmhc
        }
        
      } # end of use_saved_data
      # end of methods
      
      A_est = result$pushed_A_est_by_lam
      sel_lam = result$sel_tun
      
      if (m == "mc") {
        mc_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[thre_idx,-c(10,11)]
      } else if (m == "mo"){
        mo_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[thre_idx,-c(10,11)]
      } else {
        mm_table[iter,] <- eval_table(A_est, A_true,"obs",sel_lam)[thre_idx,-c(10,11)]
        hc_table[iter,] <- eval_table(amat(result.hc), A_true, "obs", single = T)
        mmhc_table[iter,] <- eval_table(amat(result.mmhc),A_true,"obs",single=T)
      }
    }
  } # end of iter
  
  graph_name = graph_types_long[which(graph_types == g)]
  
  write(paste(graph_name, ": n_edges = ", sum(A_true), "\n"),
        file = table_file_name, append = T)
  
  write.table(rbind(c("","","","P","E","TPR","JI","R*","M","FP","FDP","SHD")),
              file = table_file_name, append = T, sep ="\t",
              col.names = F, row.names = F, quote = F)
  
  write(paste0("-----------------------------------------------------------",
               "----------------------------------------"), 
        file = table_file_name, append = T)
  
  ave_mc_table = as.data.frame(rbind(colMeans(mc_table) %>% round(digits = 2)))
  print_mc_table = data.frame(m = "GLMDAG: Original ", ave_mc_table )
  write.table(print_mc_table, file = table_file_name, sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mo_table = as.data.frame(rbind(colMeans(mo_table) %>% round(digits = 2)))
  print_mo_table = data.frame(m = "GLMDAG: Ordinal  ", ave_mo_table )
  write.table(print_mo_table, file = table_file_name, sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mm_table = as.data.frame(rbind(colMeans(mm_table) %>% round(digits = 2)))
  print_mm_table = data.frame(m = "GLMDAG: Nominal  ", ave_mm_table )
  write.table(print_mm_table, file = table_file_name, sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  write(paste0("-----------------------------------------------------------",
               "----------------------------------------"), 
        file = table_file_name, append = T)
  
  ave_hc_table = as.data.frame(rbind(colMeans(hc_table) %>% round(digits = 2)))
  print_hc_table = data.frame(m = "HC               ", ave_hc_table )
  write.table(print_hc_table, file = table_file_name, sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  ave_mmhc_table = as.data.frame(rbind(colMeans(mmhc_table) %>% round(digits = 2)))
  print_mmhc_table = data.frame(m = "MMHC             ", ave_mmhc_table )
  write.table(print_mmhc_table, file = table_file_name, sep ="\t",
              quote = F, append = T, row.names = F,
              col.names =F)
  
  write("\n\n",  file = table_file_name, append = T)
  
  print(paste0("Graph: ", graph_name, " has been generated! ..."))
  
  
} # end of g in graph_types

# .. Check "./results/Table_3.txt". 



#++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Table 4 (Large Graphs, Nom vs Ord) ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

## !!! Debugging Block !!! ++++++++++++++++++++
graph_types = c("bi", "rand","sf","sw","tree")
# +++++++++++++++++++++++++++++++++++++++++++++

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
  
  write(paste0("-----------------------------------------------------------",
               "----------------------------------------"), 
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
  
  write(paste0("-----------------------------------------------------------",
               "----------------------------------------"), 
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Figure 6 (Original Est Graph with HCC) ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Import dataset

# Get hcc dataset from UCI ML Repository
url = paste0("https://archive.ics.uci.edu/ml/machine-learning-databases/00423/",
             "hcc-survival.zip")
download.file(url, destfile = file.path("./data/.", basename(url)))
unzip(zipfile = "./data/hcc-survival.zip", exdir = "./data/.")
file.copy("./data/hcc-survival/hcc-data.txt", "./data/.")

# Delete 
unlink("./data/hcc-survival.zip")
unlink("./data/__MACOSX/", recursive = T) # if your system is macOS
unlink("./data/hcc-survival/", recursive = T)

# read original data
hcc_temp <- read.table("./data/hcc-data.txt", sep=",", na.strings = "?")
n = nrow(hcc_temp)
hcc_temp_b4_impute = as.matrix(hcc_temp)

# Export the original data with xlsx format
hcc_info <- read.csv("./data/hcc_info.csv") 
var_names = hcc_info$var_names
names(hcc_temp) <- var_names


## Data cleaning

source("./codes/hcc_cleaning.R") # See details in "./codes/hcc_cleaning.R"
data_input = hcc_imputed
var_names_modified = names(data_input)
write.csv(hcc_imputed, "./data/hcc_imputed.csv") # save the final data


## Run GLMDAG

# Set some parameters and containers
terminal_nodes = c("class")
root_nodes = c("gend", "ende", "age")
n_lams = 30 
graph_est_by_lam = as.list(rep(NA,n_lams))
pushed_A_est_by_lam = as.list(rep(NA,n_lams))
attr(pushed_A_est_by_lam, "n_edges") <- c()
pushed_Beta_est_by_lam = pushed_A_est_by_lam

# Learning GLMDAG using hcc data
# .. Note that we parallelzing the algorithm along with paths.
for (ell in 1:n_lams){ # ell = 1
  
  if(use_saved_data != 1){
    
    result = glmdag(data_input, n_lams = n_lams, 
                    root_nodes = root_nodes, terminal_nodes = terminal_nodes,
                    path_par = T, path_par_num = ell, verbose = T)
    
  } else {
    
    result_file_name = paste0("./results/hcc/hcc_path_",ell,".rds")
    result = readRDS(file = result_file_name)
    
  }
  
  graph_est_by_lam[[ell]] <- result$graph_est_by_lam[[1]]
  pushed_A_est_by_lam[[ell]] <- result$pushed_A_est_by_lam[[1]]
  attr(pushed_A_est_by_lam, "n_edges")[ell] <-  
    attr(result$pushed_A_est_by_lam, "n_edges")
  pushed_Beta_est_by_lam[[ell]] <- result$pushed_Beta_est_by_lam[[1]]
  
}

tun_sel_lkhd(pushed_A_est_by_lam, pushed_Beta_est_by_lam, result)
# .. We choose 17th path according to the threshold (alpha) = 0.3


## Drawing Graphs

# Excluding isolated nodes (no child and no parent)
A_est_reduced = pushed_A_est_by_lam[[17]]
combined = cbind(colSums(A_est_reduced), rowSums(A_est_reduced))
sig_nodes = rowSums(combined) !=0

# Save the modified graph
A_est_reduced = A_est_reduced[sig_nodes, sig_nodes]
graph_est_reduced = graph_from_adjacency_matrix(A_est_reduced)
var_names_reduced = var_names_modified[sig_nodes]

# Set colors of nodes according to data types
node_types = sapply(hcc_imputed, class)
node_types_reduced = node_types[var_names_reduced]
nodes_colors_reduced = sapply(node_types_reduced, function(x) {
  if ( any(x %in% "numeric") ) {
    "chartreuse4" # green
  } else if (any(x %in% "ordered")) {
    "dodgerblue3" # blue
  } else{
    "firebrick1" # red
  } })

# Set root nodes as "square" nodes (otherwise "circle")
nodes_shapes = ifelse(colSums(A_est_reduced) == 0, "square", "circle")

# Set direct factors' edges as "red" (otherwise "black")
edges_info = as_edgelist(graph_est_reduced)
target.num = which(var_names_reduced == "class")
direct.factors <- which(edges_info[,2] == target.num)
edges_colors = rep("black", nrow(edges_info))
edges_colors[direct.factors] <- "red"

# Generate graph 
pdf("./results/Figure_6.pdf")
par(mar=c(0,0,0,0))
plot(graph_est_reduced, 
     edge.arrow.size=0.4,
     edge.color = edges_colors,
     edge.width = 1, 
     layout = layout_on_grid, 
     vertex.label = var_names_reduced, 
     vertex.label.cex = rep(0.8,50),
     vertex.label.family = "Helvetica",
     vertex.color = nodes_colors_reduced,
     vertex.shape = nodes_shapes,
     vertex.size = 13,
     vertex.label.color = "white")
dev.off()
# .. Note that it is the 'Figure_6' because it will be shown in Appendix.
# .. Check "./results/Figure_6.pdf".



#++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Figure 5 (Trimmed Est Graph with HCC) ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Modified graph for the trimmed version

# Find all the non-descedents of the target (class)
A.df = as.data.frame(A_est_reduced)
non.desc = 30 # class node
temp = 30
while (length(temp) > 0){
  temp1 = lapply(A.df[,temp, drop = F], function(x) which(x != 0) )
  temp = unique(unlist(temp1))
  non.desc <- unique(c(temp, non.desc))
}
var_names_trimmed <- var_names_reduced[non.desc]

# Set the position of 16 nodes
## > var_names_trimmed
## [1] "nodules" "afp"     "ps"      "ggt"     "hemo"    "major"   "albu"   
## [8] "hbcab"   "hcvab"   "age"     "alch"    "asc"     "alt"     "ast"    
## [15] "ferr"    "class"  

# Set nodes' coordinates for better illustration
layout_info <- c(
  -5.2, 3  , # nodules
  -2.5, 1.2, # afp
  5  , 0.7, # ps
  0.7, 2.2, # ggt
  5.5, 1.5, # hemo
  6  , 2  , # major
  0.5, 1.7, # albu
  4.5, 3  , # hbcab
  6.5, 3  , # hcvab
  2.5, 3  , # age
  -2  , 3  , # alch
  -4.8, 1.9, # asc
  4.5, 2.5, # alt
  -1  , 2.5, # ast
  -6  , 1  , # ferr
  0  , 0.5  # class
)
layout_info <- matrix(layout_info, ncol = 2, byrow = T)
## > data.frame (var = var_names_trimmed, layout_info) # check coordinates
A_est_trimmed <- A_est_reduced[non.desc, non.desc]
graph_est_trimmed <-  graph_from_adjacency_matrix(A_est_trimmed )

# Set colors of nodes according to data types
node_types = sapply(hcc_imputed, class)
node_types_trimmed = node_types[var_names_trimmed]
nodes_colors = sapply(node_types_trimmed, function(x) {
  if ( any(x %in% "numeric") ) {
    "chartreuse4" # green
  } else if (any(x %in% "ordered")) {
    "dodgerblue3" # blue
  } else{
    "firebrick1" # red
  } })

# Set root nodes as "square" nodes (otherwise "circle")
nodes_shapes = ifelse(colSums(A_est_trimmed) == 0, "square", "circle")

# Set direct factors' edges as "red" (otherwise "black")
edges_info = as_edgelist(graph_est_trimmed)
target.num = which(var_names_trimmed == "class")
direct.factors <- which(edges_info[,2] == target.num)
edges_colors = rep("black", nrow(edges_info))
edges_colors[direct.factors] <- "red"

# Generate graph
pdf("./results/Figure_5.pdf")
par(mar=c(0,0,0,0))
plot(graph_est_trimmed,
     edge.arrow.size = 0.4, 
     edge.color = edges_colors,
     edge.width = 1.5, 
     layout = layout_info,
     vertex.color = nodes_colors,
     vertex.label.family = "Helvetica",
     vertex.shape = nodes_shapes,
     vertex.label = var_names_trimmed, 
     vertex.label.cex= rep(1.2,length(var_names_trimmed)),
     vertex.label.color = "white",
     vertex.size = 23 )
dev.off()
# .. Check "./results/Figure_5.pdf".





