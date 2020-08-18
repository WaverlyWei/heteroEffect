#'
#' This function implements the doubly robust estimator
#'
#' @param X Covariates
#' @param D Treatment
#' @param Y observed outcome
#' @param sublabel Label of the subgroup
#' @return drPsi: doubly robust estimate of the target parameter
#' @import ranger
#' @export
#'
#' Doubly Robust Estimator Implementation
#' Answer: This function receives the fitted objects returned from ModelFit
#' ModelFit takes care of weighting. No need to weight in this function
#'
#' @param  sublabel Boolean vector of the sublable
#' @param modFit Fiited model object from ModelFit()
#'
#' @retun doubly robust estimate
CausalEffect <- function(type = c("general", "finite"),
                         sublabel = NULL,
                         weighted = FALSE,
                         modFit) {

# Initialization
  weight <- drPsi <- NA


# extract model objects
  psEst <- modFit$psEst
  outEst <- modFit$outEst
  D <- modFit$D
  Y <- modFit$Y

  # propensity score
  e1 <- psEst$e1
  e0 <- psEst$e0

  # extract weight
  if(weighted){

    weight <- modFit$weight
    # outcome estimates
    u1 <- weighted.mean(sublabel * outEst[, "mu1"], na.rm = TRUE, w = weight)
    u0 <- weighted.mean(sublabel * outEst[, "mu0"], na.rm = TRUE, w = weight)
  }

  else{
    u1 <- mean(sublabel * outEst[, "mu1"], na.rm = TRUE)
    u0 <- mean(sublabel * outEst[, "mu0"], na.rm = TRUE)
  }


  init <- u1 - u0

  # Doubly Robust Estimator
  if(type == "general"){
    if(!is.null(sublabel)){

      if(weighted){
        drPsi <- init + weighted.mean(sublabel * (D/e1 - (1 - D)/e0) *
                                        (Y - outEst[, "mu"]) + outEst[, "mu1"] -
                                        outEst[, "mu0"] - (u1 - u0), w = weight)
      }
      else{
        drPsi <- init + mean(sublabel * (D/e1 - (1 - D)/e0) *
                                        (Y - outEst[, "mu"]) + outEst[, "mu1"] -
                                        outEst[, "mu0"] - (u1 - u0))
      }
    }
    else{
      drPsi <- init + mean((D/e1 - (1 - D)/e0) *
                             (Y - outEst[, "mu"]) + outEst[, "mu1"] -
                             outEst[, "mu0"] - (u1 - u0))
    }
  }

  # doubly robust estimator
  drPsi <- init + weighted.mean(sublabel * (D/e1 - (1 - D)/e0) *
                                  (Y - outEst[, "mu"]) + outEst[, "mu1"] -
                                  outEst[, "mu0"] - (u1 - u0), w = weight)
  return(drPsi)
}

