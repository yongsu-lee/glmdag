## Import dataset

# Get hcc dataset from UCI ML Repository
url = paste0("https://archive.ics.uci.edu/ml/machine-learning-databases/00423/",
             "hcc-survival.zip")
download.file(url, destfile = file.path("./data/.", basename(url)))
unzip(zipfile = "./data/hcc-survival.zip", exdir = "./data/.")
file.copy("./data/hcc-survival/hcc-data.txt", "./data/.")

# Delete 
unlink("./data/hcc-survival.zip")
if (sysname == "Darwin") { # if your system is macOS
  unlink("./data/__MACOSX/", recursive = T) 
}
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
graph_est_by_lam = as.list(rep(NA,n_lams))
pushed_A_est_by_lam = as.list(rep(NA,n_lams))
attr(pushed_A_est_by_lam, "n_edges") <- c()
pushed_Beta_est_by_lam = pushed_A_est_by_lam

# Learning GLMDAG using hcc data

result = glmdag(data_input, n_lams = n_lams, eps_lam = eps_lam, 
                root_nodes = root_nodes, terminal_nodes = terminal_nodes,
                path_par = T, path_par_num = ell, verbose = T)

saveRDS(result, file = paste0(name_temp,".rds"))