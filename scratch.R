rm(list=ls())
test <- read.table("test.txt", header = T)
list2env(as.list(test), environment())
