#' Generate Variable and Edge Names for Graph Representation
#'
#' Creates node and edge tags for the adjacency matrix and identifies category structure for mixed graphical models.
#'
#' @param p Integer. Number of original predictors.
#' @param q Integer. Number of expanded predictors after encoding categorical variables.
#' @param categories Integer vector. Number of categories per original predictor.
#' @param var_names Integer vector. Maps expanded variables to original predictors.
#'
#' @return A list with:
#'   \item{tags}{Character vector of node names.}
#'   \item{tag_Beta}{Matrix of edge names.}
#'   \item{indicators_noedge}{Logical matrix; TRUE if edge is allowed, FALSE if not.}
#' @export
#'
get_names_graph <- function(p, q, categories, var_names){
  tag <- NULL
  for(i in 1:p){
    if(categories[i] == 1){
      tag <- c(tag, paste0(i))
    } else {
      length <- categories[i]
      for(j in 1:length){
        tag <- c(tag, paste0(i, "_", j))
      }
    }
  }
  tag_Beta <- matrix(nrow = q, ncol = q)
  ind_noedge <- matrix(TRUE, nrow = q, ncol = q)
  for(i in 1:q){
    for(j in 1:q){
      tag_Beta[i,j] <- paste0(tag[i],"-",tag[j])
      if(var_names[i] == var_names[j]){
        ind_noedge[i,j] <- FALSE
      }
    }
  }
  diag(ind_noedge) <- TRUE
  return(list("tags" = tag, "tag_Beta" = tag_Beta, "indicators_noedge" = ind_noedge))
}

#' Context-specific Graph Output
#'
#' Processes posterior Beta and inclusion samples into estimated context-specific adjacency matrices.
#'
#' @param q Number of nodes.
#' @param post_Bta Posterior Beta matrix.
#' @param post_Z Posterior inclusion matrix.
#' @param tags Node tags.
#' @param bfdr Desired Bayesian FDR.
#'
#' @return List with estimated adjacency matrices.
#' @export
context_spec_graph <- function(q, post_Bta, post_Z, tags, bfdr){
  esti_Beta <- matrix(rep(0, q*q), nrow = q, ncol = q)
  esti_Beta[upper.tri(esti_Beta)] <- colMeans(post_Bta)
  esti_Beta <- esti_Beta + t(esti_Beta)
  colnames(esti_Beta) <- tags
  rownames(esti_Beta) <- tags

  esti_Z <- matrix(rep(0, q*q), nrow = q, ncol = q)
  post_inclusion <- colMeans(post_Z)
  fdr_c <- function(c){
    E_fdr <- sum((1 - post_inclusion)*(post_inclusion > c))/(sum(post_inclusion > c) +
                                                               rnorm(1, mean = 0, sd = 0.001))
    return(E_fdr)
  }
  pos_c <- seq(0, 1, by = 0.01)
  expected_fdr <- sapply(pos_c, fdr_c)
  pos <- pos_c[expected_fdr < bfdr]
  cutoff <- min(pos)

  esti_Z[upper.tri(esti_Z)] <- colMeans(post_Z)
  esti_Z <- esti_Z + t(esti_Z)
  esti_Z <- esti_Z*(esti_Z > cutoff)
  colnames(esti_Z) <- tags
  rownames(esti_Z) <- tags

  esti_Beta <- esti_Beta*esti_Z
  return(list(ce_esti_Beta = esti_Beta, ce_esti_Z = esti_Z))
}

#' Category-based Graph Output
#'
#' Aggregates context-specific adjacency matrices into variable-level adjacency matrices.
#'
#' @param q Number of nodes after expansion.
#' @param p Number of original variables.
#' @param var_names Mapping from nodes to original variables.
#' @param esti_Z Estimated context-specific inclusion matrix.
#' @param esti_Beta Estimated context-specific Beta matrix.
#' @param categories Vector of categories per variable.
#'
#' @return List with adjacency matrices at variable level.
#' @export
categories_graph <- function(q, p, var_names, esti_Z, esti_Beta, categories){
  esti_Z_gen <- matrix(nrow = q, ncol = p)
  diag(esti_Z_gen) <- 0

  for(s in 1:p){
    if(categories[s] == 1){
      esti_Z_gen[,s] <- esti_Z[,which(var_names == s)]
    } else {
      esti_Z_gen[,s] <- pmin(1,rowSums(esti_Z[,which(var_names == s)]))
    }
  }

  esti_Z_gen_mod <- matrix(nrow = p, ncol = p)
  for(l in 1:p){
    if(categories[l] == 1){
      esti_Z_gen_mod[l,] <- esti_Z_gen[which(var_names == l),]
    } else {
      esti_Z_gen_mod[l,] <- pmin(1, colSums(esti_Z_gen[which(var_names == l),]))
    }
  }

  esti_Beta_gen <- matrix(nrow = q, ncol = p)
  diag(esti_Beta_gen) <- 0

  for(s in 1:p){
    if(categories[s] == 1){
      esti_Beta_gen[,s] <- esti_Beta[,which(var_names == s)]
    } else {
      esti_Beta_gen[,s] <- rowSums(esti_Beta[,which(var_names == s)])
    }
  }

  esti_Beta_gen_mod <- matrix(nrow = p, ncol = p)
  for(l in 1:p){
    if(categories[l ] == 1){
      esti_Beta_gen_mod[l,] <- esti_Beta_gen[which(var_names == l),]
    } else {
      esti_Beta_gen_mod[l,] <- colSums(esti_Beta_gen[which(var_names == l),])
    }
  }

  esti_Beta_gen_mod <- esti_Beta_gen_mod*lower.tri(esti_Beta_gen_mod)
  esti_Beta_gen_mod <- esti_Beta_gen_mod + t(esti_Beta_gen_mod)

  return(list("Adj_Beta" = esti_Beta_gen_mod, "Adj_Z" = esti_Z_gen_mod))
}
