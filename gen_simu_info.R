simu_info = list(
  simu_case = "simu2",
  size = "small",
  method = "mc",
  seed_para = 1,
  n_obs = 50,
  n_multi = 4,
  n_conti = 3,
  n_ordin = 3,
  n_nodes = n_multi + n_conti + n_ordin,
  n_lams = 3,
  eps_lam = 0.3
)

saveRDS(simu_info, "simu_info_test.rds")
