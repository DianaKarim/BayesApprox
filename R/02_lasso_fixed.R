#' Lasso prior model using a fixed shrinkage parameter
#'
#' This function builds a Bayesian lasso model with fixed lambda parameter.
#'
#' @param beta_MLE Initial MLE estimates.
#' @param Sigma_MLE Initial variance-covariance matrix.
#' @param lambda2 Value of the shrinkage parameter, with default of 1.
#' @param Nsample Number of iterations of the Gibbs Sample with the default of 500.
#' @param store Fraction of the iteration to save from the sampler to improve quality of mixing, with default of 10.
#' @param burnin Number of first iterations to leave out, with default of 100.
#' @return A list with a matrix of sampled coefficients, and a matrix of sampled hyper-parameter.
#' @export
lasso.fixed <- function(beta_MLE, Sigma_MLE,
                        Nsample = 500,
                        store = 10, burnin = 100,
                        lambda2 = 1){
  # Define dimensions
  P <- length(beta_MLE)

  if(!is.matrix(Sigma_MLE)){
    Sigma_MLE <- as.matrix(Sigma_MLE)
  }

  # Define the vectors to store the results
  beta_STORE    <- matrix(0,nrow = Nsample/store, ncol = P)
 # lambda2_STORE <- gamma2_STORE <- matrix(0,nrow = Nsample/store, ncol = 1)
  tau2_STORE    <- matrix(0,nrow = Nsample/store, ncol = P)

  # initial values
  tau2 <- rep(1, P)
  D <- diag(tau2*lambda2,nrow=P)
  D_inv <- diag(1/(tau2*lambda2),nrow=P)
  beta <- rep(0,P)

  # Posterior mean and variance:
  Sigma_MLE_inv <- solve(Sigma_MLE)

  # Burn-in   ----
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
  }

  # Sampler iterations ----
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

    # save each 'store' iterations
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      tau2_STORE[storecount,] <- tau2
    }

    ###################

    pb$tick()
    Sys.sleep(1/Nsample)
  }


  return(list(beta = beta_STORE,
              tau2 = tau2_STORE))
}

