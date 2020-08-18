#'
#' This function computes the propensity score
#'
#' @param X Covariates
#' @param D Treatment
#' @param Y observed outcome
#' @return initiMu: initial estimation of the expected outcome Y
#' @import ranger
#' @export
#'

outcomeModel <- function(D, X, Y, weight=NULL) {
  if(!is.null(weight)){
    muModel <- ranger(x = cbind(D, X),
                      y = Y,
                      case.weights = weight)
  }
  else{
    muModel <- ranger(x = cbind(D, X),
                      y = Y)
  }

  mu <- predict(muModel,
                data = cbind(D, X),
                type = "response")

  mat1 <- mat0 <- cbind(D, X)

  mat1[, 1] <- rep(1, nrow(mat1))

  mat0[, 1] <- rep(0, nrow(mat0))

  mu1 <- predict(muModel,
                 data = mat1,
                 type = "response")

  mu0 <- predict(muModel,
                 data = mat0,
                 type = "response")

  initMu <- cbind(mu$predictions,
                  mu0$predictions,
                  mu1$predictions)

  colnames(initMu) <- c("mu", "mu0", "mu1")

  return(initMu)
}
