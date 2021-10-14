


for (g in graph_types){ # g = "rand"
  
  graph_set = gen_graph_adj_mat(n_nodes, g, seed = seed_para)
  A_true = graph_set$A_true
  output_name = paste0(simu_title, "_", g, ".txt")
  write(x=NULL, file= output_name, append = F)
  
  for (m in methods){ # m = "mc"
    
    result_collect = array(dim = c(length(thre_list), length(eval_crit), n_iter))
    time_table = c()
    cat(paste0("Method = ", m), sep = "\n", append = T, file = output_name)
    
    for (i in 1:n_iter){ # i = 1
      
      file_name = paste0(subdir, simu_title, "_", g, "_", m, "_", i, ".rds")
      result = readRDS(file_name)
      
      A_est = result$pushed_A_est_by_lam
      sel_lam = result$sel_tun
      
      result_collect[,,i] <- as.matrix(
        eval_table(A_est, A_true, "obs", sel_lam, thre_list)[,-c(10,11)] )
      
    }
    
    ave_result = apply(result_collect, c(1,2), function(x) mean(x))
    colnames(ave_result) <- eval_crit
    ave_result_with_thre = cbind(thre = thre_list, ave_result)
    capture.output( print(ave_result_with_thre, print.gap = 3), 
                    append = T, file = output_name )
    # cat(paste0("Time: ", round(mean(time_table),2), "sec."), sep = "\n\n", 
    #     append = T,  file = output_name)
  }
  
}
