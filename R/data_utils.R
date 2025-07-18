#' Split Categorical Variables for Design Matrix
#'
#' Expands categorical variables in X and returns the new design matrix and categories per variable.
#'
#' @param X Matrix of predictors.
#' @param type Vector of variable types (e.g. "c", "d", "z", "m").
#'
#' @return List with expanded matrix and number of categories per variable.
#' @export
split_X_cat <- function(X, type){
  X_design <- c()
  categories <- c()
  p <- ncol(X)
  for(s in 1:p){
    if(type[s] == "m"){
      x <- factor(X[,s])
      X_design <- cbind(X_design, model.matrix(~x)[,-1])
      categories[s] <- length(levels(x)) - 1
    } else {
      X_design <- cbind(X_design, X[,s])
      categories[s] <- 1
    }
  }
  return(list("matrix" = X_design, "categories" = categories))
}

#' Find Lambda for F Transformation
#'
#' Computes the lambda parameter for the F-transformation based on Arkaprava & Dunson (2020).
#'
#' @param X Matrix of predictors.
#' @param type Vector of variable types.
#'
#' @return Estimated lambda.
#' @export
find_lambda <- function(X, type){
  d <- which(type != "m")
  X_d <- X[, d]
  p0 <- max(X_d)
  p_v <- seq(0, 2*p0, by = 0.1)
  tes <- c()
  for(i in p_v){
    F_X <- F_transformation(X = X_d, type = type[d], parameter = i, cont = FALSE)
    sq_error <- sqrt(sum((cov(X_d) - cov(F_X))^2))
    tes <- c(tes, sq_error)
  }
  lambda_est <- p_v[which.min(tes)]
  return(lambda_est)
}

#' F Transformation for Predictors
#'
#' Applies transformation F to predictors depending on type and lambda.
#'
#' @param X Predictor matrix or vector.
#' @param type Vector of variable types.
#' @param parameter Lambda parameter.
#' @param cont Logical; if TRUE, transforms continuous variables as well.
#'
#' @return Transformed predictors.
#' @export
F_transformation <- function(X, type, parameter, cont = FALSE){
  X_t <- X
  transform <- function(x, typ, param){
    switch(typ,
           m = x,
           d = atan(x)^param,
           z = atan(x)^param,
           c = if(cont) sign(x)*atan(abs(x))^param else x)
  }
  if(!(is.matrix(X) | is.data.frame(X))) {
    X_t <- transform(X, type, parameter)
  } else {
    for(c in 1:ncol(X)) X_t[,c] <- transform(X[,c], type[c], parameter)
  }
  return(X_t)
}
