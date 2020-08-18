
#'
#' This function computes the propensity score
#'
#' @param X Covariates
#' @param D Treatment
#' @param weight Observation weights
#' @return e1: propensity score of the treated; e0: propensity score of the control
#' @import ranger
#' @export
#'
ps <- function(X, D, weight = NULL) {
    if(!is.null(weight)){

      psModel <- ranger(x = X,
                        y = D,
                        case.weights = weight)

    }
  else{
    psModel <- ranger(x = X,
                      y = D)

  }

  e1 <- predict(psModel,
                data = X,
                type = "response")

  # propensity score the treated
  e1 <- e1$predictions

  # propensity score the control
  e0 <- 1 - e1

  return(list(e1 = e1, e0 = e0))
}
