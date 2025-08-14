## Set of R functions for the optimal block length computation for the
## Sharpe ratio

# #' @name .sharpeBlockSize
# #' @title See sharpeBlockSize
# #' @importFrom stats rgeom
# #' @import compiler
.sharpeBlockSize <- function(x, y, control = list(), b.vec = c(1, 3, 6, 10), 
                             alpha = 0.05, M = 199, K = 500, b.av = 5, 
                             T.start = 50) {
  
  sb.sequence <- function(T, b.av, length = T) {
    
    index.sequence <- c(1:T, 1:T)
    sequence <- rep.int(0, length + T)
    current <- 0
    while (current < length) {
      start <- sample(1:T, 1)
      b <- stats::rgeom(1, 1/b.av) + 1
      sequence[(current + 1):(current + b)] <- index.sequence[start:(start + 
                                                                       b - 1)]
      current <- current + b
    }
    out <- sequence[1:length]
    return(out)
    
  }
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  rets <- cbind(x, y)
  
  # process control parameters
  ctr <- processControl(control)
  
  b.len <- length(b.vec)
  emp.reject.probs <- rep.int(0, b.len)
  d <- sharpe.ratio.diff(x, y, ctr$ttype)
  T <- length(x)
  Var.data <- matrix(data = 0, nrow = T.start + T, ncol = 2)
  Var.data[1, ] <- rets[1, ]
  fit1 <- stats::lm(x[2:T] ~ x[1:(T - 1)] + y[1:(T - 1)])
  fit2 <- stats::lm(y[2:T] ~ x[1:(T - 1)] + y[1:(T - 1)])
  coef1 <- as.numeric(fit1$coef)
  coef2 <- as.numeric(fit2$coef)
  resid.mat <- cbind(as.numeric(fit1$resid), as.numeric(fit2$resid))
  for (k in 1:K) {
    ids <- sb.sequence(T - 1, b.av, T.start + T - 1)
    resid.mat.star <- rbind(c(0, 0), resid.mat[ids, ])
    for (t in 2:(T.start + T)) {
      Var.data[t, 1] <- coef1[1] + coef1[2] * Var.data[t - 1, 1] + 
        coef1[3] * Var.data[t - 1, 2] + resid.mat.star[t, 1]
      Var.data[t, 2] <- coef2[1] + coef2[2] * Var.data[t - 1, 2] + 
        coef2[3] * Var.data[t - 1, 2] + resid.mat.star[t, 2]
    }
    Var.data.trunc <- Var.data[(T.start + 1):(T.start + T), ]
    for (j in 1:b.len) {
      bsids <- bootIndices(T, M, b.vec[j])
      tmp <- sharpeTestBootstrap(Var.data.trunc, bsids, b.vec[j], 
                                 ctr$ttype, pBoot = 1, d)
      if (tmp$pval <= alpha) {
        emp.reject.probs[j] <- emp.reject.probs[j] + 1
      }
    }
  }
  emp.reject.probs <- emp.reject.probs/K
  b.order <- order(abs(emp.reject.probs - alpha))
  b.opt <- b.vec[b.order[1]]
  return(b.opt)
}
sharpeBlockSize <- compiler::cmpfun(.sharpeBlockSize)
