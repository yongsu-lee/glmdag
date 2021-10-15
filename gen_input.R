## Generate Simulation Information ####
simu_info = list(
  simu_case = "mmhc",
  size = "small",
  seed_para = 1,
  n_obs = 50,
  n_multi = 4,
  n_conti = 3,
  n_ordin = 3,
  n_lams = 30,
  eps_lam = 0.05,
  path_par = FALSE
)
saveRDS(simu_info, "simu_info.rds")

## Generate queue list ####
queue_list = function(args_list, file_name){
  
  n_args = length(args_list)
  combn = expand.grid(args_list)
  write.table(combn, file = file_name, sep=", ", 
              quote = F, row.names = F, col.names = F)
}

args_list = list(c("bi", "rand", "sf", "sw", "tree" ), c("mc", "mo", "mm"), 1:100)
queue_list(args_list, "queue_list")