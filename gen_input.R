## Generate Simulation Information ####
simu_info = list(
  simu_case = "simu2",
  size = "large",
  seed_para = 1,
  n_obs = 50,
  n_multi = 4,
  n_conti = 4,
  n_ordin = 4,
  n_lams = 3,
  eps_lam = 0.4,
  path_par = TRUE
)

saveRDS(simu_info, "simu_info.rds")

## Generate queue list ####
queue_list = function(args_list, file_name){
  
  n_args = length(args_list)
  combn = expand.grid(args_list)
  write.table(combn, file = file_name, sep=", ", 
              quote = F, row.names = F, col.names = F)
}

args_list = list(c("tree"), c("mc","mo","mm"), 1:50, 1:30)
queue_list(args_list, "queue_list")