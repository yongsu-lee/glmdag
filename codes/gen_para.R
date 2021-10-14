#' Generate true parameters using an adjacency matrix
#'
#' @param A_true binary matrix, a true adjacecny matrix from a desired graph.
#' @param types_by_node 
#' @param n_levels_by_node 
#' @param ordin_pred_as blah
#' @param intercept
#' @param range a possible range of (absolute value of) parameters.
#' @param prob_sparse double scalar, define a sparseness of parameters within each parameter sub-matrix with \eqn{0\le prob_sparse \le 1} - default is 0. If \eqn{prob_sparse = 0}, all the parameter for each submatrix will be non-zeros. Not recommnedable, but, if \eqn{prob_sparse = 1}, all the parameter for each submatrix will be zeros. See details below for the parameter sub-matrix.
#' @param seed 
#'
#' @return double matrix, true parameter matrix, \code{W_true} is returend with \code{n_levels} or \code{n_group_elmts} as an attribute.
#' @export
#' @importFrom stats runif
#' @importFrom igraph graph_from_adjacency_matrix is_dag
#' @details The parameter matrix can be divided into submatrices. For example, a submatrix (i,j) consists of parameters related to (node i) -> (node j). The \code{prob_sparse} controls the sparseness of parameters within each submatrix.
#'
#' @examples
#' n_nodes = 10
#' graph_set = gen_graph_adj_mat(n_nodes = n_nodes, graph_type = "bi")
#' A_true = graph_set$A_true
#'
#' ## Discrete Case ####
#' gen_para(A_true = A_true, n_levels = rep(2, n_nodes), intercept = TRUE)
#'
#' ## Continuous Case ####
#' gen_para(A_true = A_true, n_grp_elmts = rep(2, n_nodes), intercept = FALSE)
gen_para=function(A_true, types_by_node=NULL, n_levels_by_node=NULL,
                  ordin_pred_as = "num",
                  intercept=c("always","for_child_only","none"),
                  range=c(0.5,1), prob_sparse=0, seed=1){
  
  if (FALSE){
    ordin_pred_as = "num"
    intercept = "always"
    range = c(0.5,1)
    prob_sparse = 0
    seed = 1
  }
  
  if (missing(intercept)) {
    intercept = "always"
  } else {
    intercept = match.arg(intercept)
  }
  
  graph_temp = igraph::graph_from_adjacency_matrix( A_true )
  
  if (!igraph::is_dag(graph_temp)) {
    stop("Provided adjacency matrix does not represent a DAG. Check the A_true.")
  }
  
  if (is.null(n_levels_by_node)){
    stop("n_levels_by_node is required!")
  }
  
  if (any(n_levels_by_node[(types_by_node == "m")] == 1)){
    stop("The number of levels for categorical nodes should be more than 1")
  }
  
  if (any(n_levels_by_node[(types_by_node == "o")] == 1)){
    stop("The number of levels for categorical nodes should be more than 1")
  }
  
  n_nodes = length(types_by_node)
  n_resps_by_node = ifelse(types_by_node == "c", 
                           n_levels_by_node, n_levels_by_node-1)
  
  idx_resps_by_node = gen_indices_by_node(n_resps_by_node)
  
  n_expls_by_node = n_resps_by_node
  if (ordin_pred_as == "num"){
    n_expls_by_node[types_by_node == "o"] = 1
  }
  
  idx_expls_by_node = gen_indices_by_node(n_expls_by_node)
  
  n_row = sum(n_expls_by_node)
  n_col = sum(n_resps_by_node)
  
  W_temp = matrix(0, n_row, n_col)
  
  set.seed(seed)
  for (j in 1:n_nodes){ #j = 1
    
    idx_resps_j = idx_resps_by_node[[j]]
    
    for (i in 1:n_nodes){ # i = 2
      
      idx_expls_i = idx_expls_by_node[[i]]
      
      if (A_true[i,j] == 1){
        
        n_cell = length(idx_expls_i) * length(idx_resps_j)
        W_temp[idx_expls_i, idx_resps_j] = sample(c(0, 1), n_cell, replace = T,
                                                  prob = c(prob_sparse, 1-prob_sparse))
      }
    }
  }
  
  
  # Adding intercept term ----
  if (intercept == "for_child_only") {
    
    idx_intercepts = which(apply(W_temp, 2, sum) != 0)
    W_true = rbind(W_temp, rep(0,n_col))
    
    for (j in idx_intercepts) W_true[n_row + 1, j] = 1
    
  } else if (intercept == "always") {
    
    W_true = rbind(W_temp, rep(1,n_col))
    
  } else {
    
    W_true = W_temp
  }
  
  # Assign generated coefficients  ----
  min_para = min(range); max_para = max(range)
  idx_nnz = which(W_true  == 1)
  n_nnz = sum(W_true)
  coefs = round(runif(n_nnz, min_para, max_para), 4)
  
  # Flipping the signs ----
  signs = sample(c(-1,1), n_nnz, replace = T)
  W_true[idx_nnz] = coefs * signs
  
  # Attaching types of nodes to W_true ----
  attr(W_true, "types_by_node") = types_by_node
  attr(W_true, "n_levels_by_node") = n_levels_by_node
  if (any(types_by_node == "o")){
    attr(W_true, "ordin_pred_as") = ordin_pred_as 
  }
  attr(W_true, "n_resps_by_node") = n_resps_by_node
  
  colnames(W_true) <- resp_names(n_nodes, types_by_node, n_resps_by_node)
  
  if (intercept == "none") {
    rownames(W_true) <- c(expl_names(n_nodes, types_by_node, n_expls_by_node,
                                     ordin_pred_as = ordin_pred_as))
    
  } else {
    rownames(W_true) <- c(expl_names(n_nodes, types_by_node, n_expls_by_node,
                                     ordin_pred_as = ordin_pred_as), "(Intcpt)")}
  
  
  return(W_true)
  
}
