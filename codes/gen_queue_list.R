queue_list = function(args_list){
  
  n_args = length(args_list)
  combn = expand.grid(args_list)
  write.table(combn, file = "queue_list", sep=", ", 
              quote = F, row.names = F, col.names = F)
}

args_list = list(c("tree"), c("mc","mo","mm"), 1:50, 1:30)


# args_list = list(1:20, c("alarm"), c("FALSE"))
queue_list(args_list)