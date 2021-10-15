simu_info = list(
  simu_case = "simu2",
  size = "small",
  method = NULL,
  seed_para = 1,
  n_obs = 50,
  n_nodes = 30,
  n_lams = 30,
  eps_lam = 0.3
)

saveRDS(simu_info, "simu_info_test.rds")
