#' This is a Demo Example of the AVEmGLMM package
#'
#' It sets up an estimation problem, and performs the estimation with methods AVE and WAVE
#'@return Returns the demo example output. It runs as a script with no input.
#'@export
demoExample <- function() {

  library(MASS)
    
  n <- 100  # sample size
  Q <- 4 # number of longitudinal binary outcomes
  P <- choose(Q,2) # number of pairs of items
  times <- seq(0, 5, length.out = 11)
  times <- times - (0+times[length(times)])/2    # center around zero
  id <- rep(seq_len(n), each = length(times))
  betas <- c(-1, 0.5)
  D <- bdiag(rep(list(cbind(c(1, 0.2), c(0.2, 0.25))), Q))
  D[D == 0] <- c(0.1)
  D[3,1]=0.15
  D[5,1]=0.15
  D[7,1]=0.15
  D[5,3]=0.15
  D[7,3]=0.15
  D[7,5]=0.15
  
  D[1,3]=0.15
  D[1,5]=0.15
  D[1,7]=0.15
  D[3,5]=0.15
  D[3,7]=0.15
  D[5,7]=0.15
  
  
  X <- cbind(1, rep(times, n))
  Z <- cbind(1, rep(times, n))
  ncz <- ncol(Z)
  
  m <- 1
  set.seed(1000 + m)
  
  b <- mvrnorm(n, rep(0, Q*2), D)
  
  Data <- generateData(id,times,n,X,Z,betas,b,Q)
  
  modelFit <- estimateModelFit(Data,Q,n)
  
  betasN=c(-1,0.5,1,0.25,0.2,  -1,0.5,1,0.25,0.2, -1,0.5,1,0.25,0.2, -1,0.5,1,0.25,0.2)
  
  modelFitAndTarget <- as.data.frame(cbind("True" = betasN, modelFit))
  
  return(modelFitAndTarget)
  
}