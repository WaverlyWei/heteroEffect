#' Fit propensity score and outcome model on all observations
#'
#' @param X Covariates
#' @param D Treatment
#' @param Y Outcome
#' @param weight case weight
#'
#' @return Fitted model objects
#'
ModelFit <- function(X, D, Y, weight=NULL){

  if(!is.null(weight)){
    # propensity score estimates
    psEst <- ps(X = X, D = D, weight = weight)
    # outcome estimates
    outEst <- outcomeModel(D = D, X = X, Y = Y, weight = weight)
  }

  else{
    # propensity score estimates
    psEst <- ps(X = X, D = D)
    # outcome estimates
    outEst <- outcomeModel(D = D, X = X, Y = Y)
  }

  return(list(psEst = psEst, outEst = outEst,
              D = D, Y=Y, weight = weight))

}
