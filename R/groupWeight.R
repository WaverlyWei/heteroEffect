#' Get Group Weight
#'
#' This function computes the subgroup weight.
#'
#' @param sublabel label of the subgroup
#' @param weight observation weights
#' @return subgroup weights, a numeric value
#' @export


groupWeight <- function(sublabel,weight = NULL) {

  if(!is.null(weight)){
    return(sum(weight[sublabel]) / sum(weight))
  }

  else{
    return(sum(sublabel)/length(sublabel))
  }

}
