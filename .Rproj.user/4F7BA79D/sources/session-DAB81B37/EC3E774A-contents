#' Lasso prior model.
#'
#' This function builds a Bayesian lasso model.
#'
#' @param beta_MLE Initial MLE estimates.
#' @param Sigma_MLE Initial variance-covariance matrix.
#' @param Nsample Number of iterations of the Gibbs Sample with the default of 500.
#' @param store Fraction of the iteration to save from the sampler to improve quality of mixing, with default of 10.
#' @param burnin Number of first iterations to leave out, with default of 100.
#' @param a1 Hyper-parameter of mixing distribution with the default value of 0.5.
#' @param a2 Hyper-parameter of mixing distribution with the default value of 0.5.
#' @param b1 Hyper-parameter of mixing distribution with the default value of 1.
#' @param variableNames Character vector of the effects' names.
#' @return A list with a matrix of sampled coefficients, and a matrix of sampled hyper-parameter.
#' @export
lasso <- function(beta_MLE, Sigma_MLE,
                  Nsample = 500, store = 10, burnin = 100,
                  a1 = 0.5, a2 = 0.5, b1 = 1,
                  variableNames){
  # Define dimensions
  P <- length(beta_MLE)
  # Check the covariance matrix:
  if(!is.matrix(Sigma_MLE)){
    Sigma_MLE <- as.matrix(Sigma_MLE)
  }

  Sigma_MLE_inv <- solve(Sigma_MLE)

  #Initial values of parameters:
  lambda2 <- delta <- 1
  tau2 <- rep(1, P)
  D_inv <- diag(1/(tau2*lambda2),nrow=P)


  # Define the vectors to store the results
  beta_STORE    <- matrix(0,nrow = Nsample/store, ncol = P)
  lambda2_STORE <- delta_STORE <- matrix(0,nrow = Nsample/store, ncol = 1)
  tau2_STORE    <- matrix(0,nrow = Nsample/store, ncol = P)


  for (t in 1:burnin){
    # 1. sample beta
    # beta updates through tau
    var_beta <- solve(Sigma_MLE_inv + D_inv)
    mu_beta  <- var_beta%*%Sigma_MLE_inv%*%beta_MLE

    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 2. Sample 1/tau^2
    mu_tau <- sqrt(lambda2/beta^2)
    tau2 <- 1/rinvgauss(P, mean = mu_tau, shape = 1)
    D_inv <- diag(1/(tau2 * lambda2),nrow=length(tau2))

    # 3. Sample lambda2 & delta
    lambda2 <- rinvgamma(1, a1 + P/2, delta + sum(beta^2/(2*tau2)) )
    delta <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2)

  }

  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = Nsample, clear = F, width = 80)
  storecount <- 0
  for (t in 1:Nsample){

    # 1. sample beta
    # beta updates through tau
    var_beta <- solve(Sigma_MLE_inv + D_inv)
    mu_beta  <- var_beta%*%Sigma_MLE_inv%*%beta_MLE

    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 2. Sample 1/tau^2
    mu_tau <- sqrt(lambda2/beta^2)
    tau2 <- 1/rinvgauss(P, mean = mu_tau, shape = 1)
    D_inv <- diag(1/(tau2 * lambda2),nrow=length(tau2))

    # 3. Sample lambda2 & delta
    lambda2 <- rinvgamma(1, a1 + P/2, delta + sum(beta^2/(2*tau2)) )
    delta <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2)

    # save each 'store' iterations
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      tau2_STORE[storecount,] <- tau2
      lambda2_STORE[storecount, ] <- lambda2
      delta_STORE[storecount, ] <- delta

    }

    ###################

    pb$tick()
   Sys.sleep(1/Nsample)
  }

  colnames(beta_STORE) <- variableNames
  return(list(beta = beta_STORE,
              tau2 = tau2_STORE,
              lambda2 = lambda2_STORE,
              delta = delta_STORE))
}


