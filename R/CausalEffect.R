#'
#' This function implements the doubly robust estimator
#'
#' @param type Use "general" when the sample size is large, use "finite" when the sample size is limited
#' @param sublabel boolean vector, default value is NULL
#' @return drPsi: doubly robust estimate of the target parameter
#' @import ranger
#' @export
#'
CausalEffect <- function(type = c("general", "finite"),
                         sublabel = NULL,
                         weighted = FALSE,
                         modFit) {

# Initialization
  weight <- drPsi <- subPsi <- psi <-  NA


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
  # Only estimate general case NOT subgroup case
  if(type == "general"){
      if(weighted){
        drPsi <- init + weighted.mean((D/e1 - (1 - D)/e0) *
                                        (Y - outEst[, "mu"]) + outEst[, "mu1"] -
                                        outEst[, "mu0"] - (u1 - u0), w = weight)

      }else{
        drPsi <- init + mean((D/e1 - (1 - D)/e0) *
                               (Y - outEst[, "mu"]) + outEst[, "mu1"] -
                               outEst[, "mu0"] - (u1 - u0))
      }
    psi <- drPsi

  }else if(type == "finite"){
    # subgroup case
    if(!is.null(sublabel)){
      # Compute Perturbing Directions for the treated and control
      S1 <- sublabel/P * D/e1
      S0 <- sublabel/P * (1 - D)/e0

      if(weighted){

        weight <- modFit$weight
        # estimate of the group weight
        P <- groupWeights(sublabel, weight)

        # Compute step size epsilon
        epsilon <- coef(glm(Y ~ -1 + offset(outEst[, "mu"]) + S0 + S1,
                            family = "binomial",
                            weights = weight))

        # Update the initial potential outcome model
        mustar <- outEst + c((epsilon[1] * S0 + epsilon[2] * S1),
                             epsilon[1]/e0,
                             epsilon[2]/e1)

        colnames(mustar) <- c("mu", "mu0", "mu1")

        mustar <- plogis(mustar)

        # calc parameters u1: empirical average of E[Y|D=1,X]
        u1 <- weighted.mean(mustar[, "mu1"],w = weight)
        u0 <- weighted.mean(mustar[, "mu0"],w = weight)
      }else{

        P <- groupWeights(sublabel)

        # Compute step size epsilon
        epsilon <- coef(glm(Y ~ -1 + offset(outEst[, "mu"]) + S0 + S1,
                            family = "binomial"))

        # Update the initial potential outcome model
        mustar <- outEst + c((epsilon[1] * S0 + epsilon[2] * S1),
                             epsilon[1]/e0,
                             epsilon[2]/e1)

        colnames(mustar) <- c("mu", "mu0", "mu1")

        mustar <- plogis(mustar)


        u1 <- mean(mustar[, "mu1"],na.rm = TRUE)
        u0 <- mean(mustar[, "mu0"],na.rm = TRUE)
      }



    }else{ #Finite sample, no subgroup

      S1 <-  D/e1
      S0 <- (1 - D)/e0

      if(weighted){

        weight <- modFit$weight

        # Compute step size epsilon
        epsilon <- coef(glm(Y ~ -1 + offset(outEst[, "mu"]) + S0 + S1,
                            family = "binomial",
                            weights = weight))

        # Update the initial potential outcome model
        mustar <- outEst + c((epsilon[1] * S0 + epsilon[2] * S1),
                             epsilon[1]/e0,
                             epsilon[2]/e1)

        colnames(mustar) <- c("mu", "mu0", "mu1")

        mustar <- plogis(mustar)

        # calc parameters u1: empirical average of E[Y|D=1,X]
        u1 <- weighted.mean(mustar[, "mu1"],w = weight)
        u0 <- weighted.mean(mustar[, "mu0"],w = weight)

      }else{  # finite, no subgroup, no weight

        epsilon <- coef(glm(Y ~ -1 + offset(outEst[, "mu"]) + S0 + S1,
                            family = "binomial"))

        # Update the initial potential outcome model
        mustar <- outEst + c((epsilon[1] * S0 + epsilon[2] * S1),
                             epsilon[1]/e0,
                             epsilon[2]/e1)

        colnames(mustar) <- c("mu", "mu0", "mu1")

        mustar <- plogis(mustar)

        u1 <- mean(mustar[, "mu1"],na.rm = TRUE)
        u0 <- mean(mustar[, "mu0"],na.rm = TRUE)
    }
    }

    subPsi <- u1 - u0
    psi <- subPsi

  }else{
    stop("Tyep is specified incorrectly")
  }

  return(psi)
}

