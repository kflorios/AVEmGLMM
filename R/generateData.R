#' This generates the Data
#'
#' It acts on the primitives of the estimation problem and returns the Data
#'@param id
#'@param times
#'@param n
#'@param X
#'@param Z
#'@param betas
#'@param b
#'@param Q
#'@return Returns the generated Data as data.frame
#'@export
generateData <- function(id,times,n,X,Z,betas,b,Q) {
  
  ncz <- ncol(Z)
  Data <- data.frame(id = id, time = rep(times, n))
  for (q in seq_len(Q)) {
    eta <- c(X %*% betas) + rowSums(Z * b[id, c(q*ncz - 1, q*ncz)])
    ###eta <- c(X %*% betas) + rowSums(Z * b[id, c(q*ncz)])
    Data[[paste("y", q, sep = "")]] <- rbinom(length(eta), 1, plogis(eta))
  }
  
  
  
  return(Data)
}