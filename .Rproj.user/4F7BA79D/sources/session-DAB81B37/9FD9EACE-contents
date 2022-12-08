#' Horseshoe prior model using a fixed shrinkage parameter
#'
#' This function builds a Bayesian horseshoe model with fixed lambda parameter.
#'
#' @param beta_MLE Initial MLE estimates.
#' @param Sigma_MLE Initial variance-covariance matrix.
#' @param Nsample Number of iterations of the Gibbs Sample with the default of 500.
#' @param store Fraction of the iteration to save from the sampler to improve quality of mixing, with default of 10.
#' @param burnin Number of first iterations to leave out, with default of 100.
#' @param a1 Hyper-parameter of mixing distribution with the default value of 0.5.
#' @param a2 Hyper-parameter of mixing distribution with the default value of 0.5.
#' @param a3 Hyper-parameter of mixing distribution with the default value of 0.5.
#' @param a4 Hyper-parameter of mixing distribution with the default value of 0.5.
#' @param b1 Hyper-parameter of mixing distribution with the default value of 1.
#' @param b2 Hyper-parameter of mixing distribution with the default value of 1.
#' @param variableNames Character vector of the effects' names.
#' @return A list with a matrix of sampled coefficients, and a matrix of sampled hyper-parameter.
#' @export
horseshoe <- function(beta_MLE, Sigma_MLE,
                            a1 = .5, a2 = 0.5, a3 = .5, a4 = 0.5,
                            b1 = 1, b2 = 1,
                            Nsample = 500, burnin = 1000,
                            store = 10, variableNames){
  # Define dimensions
  P <- length(beta_MLE)
  # Check the covariance matrix:
  if(!is.matrix(Sigma_MLE)){
    Sigma_MLE <- as.matrix(Sigma_MLE)
  }

  Sigma_MLE_inv <- solve(Sigma_MLE)

  #Initial values of parameters:
  tau2 <- psi2 <- rep(1, P)
  lambda2 <- gamma2 <- 1
  D <- diag(tau2*lambda2)
  D_inv <- diag((tau2*lambda2)^(-1))
  beta <- rep(0,P)


  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  lambda2_STORE <- gamma2_STORE <- matrix(0,nrow = Nsample/store, ncol = 1)
  tau2_STORE <- psi2_STORE <- matrix(0,nrow = Nsample/store, ncol = P)


  # Burn-in

  for (t in 1:burnin){
    # 1. sample beta
    var_beta <- solve(Sigma_MLE_inv + D_inv)
    mu_beta  <- var_beta%*%Sigma_MLE_inv%*%beta_MLE

    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 2. Sample tau2 (vector)
    tau2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*lambda2))

    # 3. Sample lambda2 (scalar)
    lambda2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*tau2)))
    D_inv <- diag(1/(tau2*lambda2))

    # 4.  Sample gamma2
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2)

    # 5. Sample psi2
    psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 +1/tau2)

  }

  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                       total = Nsample, clear = F, width = 80)
  storecount <- 0
  for (t in 1:Nsample){
    # 1. sample beta
    var_beta <- solve(Sigma_MLE_inv + D_inv)
    mu_beta  <- var_beta%*%Sigma_MLE_inv%*%beta_MLE

    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))

    # 2. Sample tau2 (vector)
    tau2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*lambda2))

    # 3. Sample lambda2 (scalar)
    lambda2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*tau2)))
    D_inv <- diag(1/(tau2*lambda2))

    # 4.  Sample gamma2
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2)

    # 5. Sample psi2
    psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 +1/tau2)

    # save each 'store' iterations
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      lambda2_STORE[storecount,] <- lambda2
      tau2_STORE[storecount,] <- tau2
      gamma2_STORE[storecount,] <- gamma2
      psi2_STORE[storecount,] <- psi2

    }
    ###################

    pb$tick()
    Sys.sleep(1/Nsample)
  }

  colnames(beta_STORE) <- variableNames


  return(list(beta = beta_STORE,
              lambda2 = lambda2_STORE, gamma2 = gamma2_STORE, psi2 = psi2_STORE,
              tau2 = tau2_STORE))
}











