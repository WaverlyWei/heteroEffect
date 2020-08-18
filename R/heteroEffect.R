#'
#' This function implements the heteroEffect estimator
#'
#' @param dat dataset to analyze
#' @param X Covariates
#' @param D Treatment
#' @param Y observed outcome
#' @param sublabel Label of the subgroup
#' @return psi: heteroEffect estimate of the target parameter
#' @import ranger
#' @export
#'
heteroEffect <- function(dat,X, D, Y, sublabel) {

  P <- groupWeights(sublabel = sublabel)

  # propensity score estimates
  psEst <- ps(X = X,
              D = D)

  e1 <- psEst$e1

  e0 <- psEst$e0

  # outcome estimates
  outEst <- outcomeModel(D = D,
                         X = X,
                         Y = Y)

  m1 <- mean(outEst[, "mu1"])

  m0 <- mean(outEst[, "mu0"])

  # Compute Perturbing Directions for the treated and control
  S1 <- (dat$IssuingBankName == sublabel)/P * D/e1

  S0 <- (dat$IssuingBankName == sublabel)/P * (1 - D)/e0

  # compute step size epsilon

  # Note: this step can be coded iterative optimization
  # use glm to compute optimal epsilon is a more efficient
  # one step method

  epsilon <- coef(glm(Y ~ -1 + offset(outEst[, "mu"]) + S0 + S1,
                      family = "binomial"))

  # Update the initial potential outcome model
  mustar <- outEst + c((epsilon[1] * S0 + epsilon[2] * S1), epsilon[1]/e0,
                       epsilon[2]/e1)

  colnames(mustar) <- c("mu", "mu0", "mu1")

  mustar <- plogis(mustar)

  # calculate parameters
  # u1: empirical average of E[Y|D=1,X]
  u1 <- mean(mustar[, "mu1"])

  # u0: empirical average of E[Y|D=0,X]
  u0 <- mean(mustar[, "mu0"])

  # ATE
  psi <- u1 - u0


  return(psi)
}
