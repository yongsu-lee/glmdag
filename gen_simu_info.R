simu_info = list(
  simu_case = "simu1",
  size = "small",
  method = NULL,
  seed_para = 1,
  n_obs = 50,
  n_nodes = 10,
  n_lams = 3,
  eps_lam = 0.3
)

saveRDS(simu_info, "simu_info_test.rds")
