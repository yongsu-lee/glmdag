## Generate Simulation Information ####
simu_info = list(
  simu_case = "simu2",
  size = "small",
  method = "mc",
  seed_para = 1,
  n_obs = 50,
  n_multi <- 4,
  n_conti <- 3,
  n_ordin <- 3,
  n_nodes = n_multi + n_conti + n_ordin,
  n_lams = 3,
  eps_lam = 0.3,
  path_par = FALSE
)

saveRDS(simu_info, "simu_info_test.rds")

## Generate queue list ####
queue_list = function(args_list){
  
  n_args = length(args_list)
  combn = expand.grid(args_list)
  write.table(combn, file = "queue_list", sep=", ", 
              quote = F, row.names = F, col.names = F)
}

args_list = list(c("tree"), c("mc","mo","mm"), 1:50, 1:30)

# args_list = list(1:20, c("alarm"), c("FALSE"))

queue_list(args_list)