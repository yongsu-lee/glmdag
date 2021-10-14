setwd("~/iCloud/paper_a/paper_a_submit/biom_j/biom_j_code_submit/results/simu2_large_makeup/")

# subdir = "simu2_large_sw_"
subdir = "simu2_large_bi_"
n_lams = 30 
# name.list = c("mm_79", "mm_91", "mo_9", "mo_31")
name.list = c("mc_100")

for (k in seq_along(name.list)){ # k = 1
  
  graph_est_by_lam = as.list(rep(NA,n_lams))
  pushed_A_est_by_lam = as.list(rep(NA,n_lams))
  attr(pushed_A_est_by_lam, "n_edges") <- c()
  pushed_Beta_est_by_lam = pushed_A_est_by_lam
  
  for (ell in 1:n_lams){ # ell = 1
    
    file.name.temp = paste0(subdir, name.list[k])  
    file.name.path = paste0(file.name.temp, "_", ell,".rds")
    result = readRDS(file = file.name.path)
    
    graph_est_by_lam[[ell]] <- result$graph_est_by_lam[[1]]
    pushed_A_est_by_lam[[ell]] <- result$pushed_A_est_by_lam[[1]]
    attr(pushed_A_est_by_lam, "n_edges")[ell] <-  
      attr(result$pushed_A_est_by_lam, "n_edges")
    pushed_Beta_est_by_lam[[ell]] <- result$pushed_Beta_est_by_lam[[1]]
    
  }
  
  sel_tun = tun_sel_lkhd(pushed_A_est_by_lam, pushed_Beta_est_by_lam, result)
  
  fit = list(graph_est_by_lam = graph_est_by_lam,
             sel_tun = sel_tun,
             pushed_A_est_by_lam = pushed_A_est_by_lam,
             pushed_Beta_est_by_lam = pushed_Beta_est_by_lam)
  
  file.name.combined = paste0(file.name.temp, ".rds")
  saveRDS(fit, file.name.combined )
  
  
}

