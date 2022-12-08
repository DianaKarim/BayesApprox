#' Ridge prior model using a fixed shrinkage parameter
#'
#' This function builds a Bayesian ridge model with fixed lambda parameter.
#' The posterior of the coefficients in this case is normal.
#'
#' @param beta_MLE Initial MLE estimates.
#' @param lambda2 Value of the shrinkage parameter, with default of 1.
#' @param Nsample Number of iterations of the Gibbs Sample with the default of 500.
#' @return A matrix with Nsample rows and number of columns equal to the number of coefficients.
#' @export
ridge.fixed <- function(beta_MLE, Sigma_MLE,
                        lambda2 = 1,
                        Nsample = 500){

  #print("Define pre-iteration variables ...")
  ## Define dimensions
  P <- length(beta_MLE)
  ## Define the vectors to store the results

  if(!is.matrix(Sigma_MLE)){
    Sigma_MLE <- as.matrix(Sigma_MLE)
  }

  lambda2Input <- lambda2 # value of the shrinkage parameter


  D <- diag(rep(lambda2Input, P), nrow = P)
  D_inv <- diag(rep(lambda2Input, P)^(-1), nrow = P)
  beta <- rep(0,P)


  Sigma_MLE_inv <- solve(Sigma_MLE)
  var_beta <- solve(Sigma_MLE_inv + D_inv)
  mu_beta  <- var_beta%*%Sigma_MLE_inv%*%beta_MLE

  beta <- rmvnorm(Nsample, mean = mu_beta, sigma = var_beta)

  if (P==1) {
    beta <- c(beta)
  }


  return(list(beta = beta))
}






