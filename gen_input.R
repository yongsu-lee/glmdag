rm(list=ls())

# ## Generate Simulation Information (Large) #####

# simu_info = list(
#   simu_case = "simu2",
#   size = "large",
#   seed_para = 1,
#   n_obs = 200,
#   n_multi = 15,
#   n_conti = 15,
#   n_ordin = 20,
#   n_lams = 30,
#   eps_lam = 0.2,
#   path_par = TRUE
# )


## Generate Simulation Information (Small) #####
simu_info = list(
  simu_case = "mmhc",
  size = "small",
  n_iter = 100,
  seed_para = 100,
  n_obs = 50,
  n_multi = 4,
  n_conti = 3,
  n_ordin = 3,
  n_lams = 30,
  eps_lam = 0.5,
  path_par = FALSE
)
list2env(simu_info, globalenv())
saveRDS(simu_info, "~/Dropbox/glmdag/simu_info_queue/simu_info.rds")

## Generate queue list ####
queue_list = function(args_list, file_name){
  
  n_args = length(args_list)
  combn = expand.grid(args_list)
  write.table(combn, file = file_name, sep=", ", 
              quote = F, row.names = F, col.names = F)
}

args_list = list(c("bi", "rand", "sf", "sw", "tree"), "mm", 1:n_iter)

queue_list(args_list, "~/iCloud/glmdag/simu_info_queue/queue_list")
