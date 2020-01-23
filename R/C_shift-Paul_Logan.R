orig_renorm <- function(Cov_obs, g, n, eps = 0.001) {
  cntr <- 0
  r <- rep(0, 2)
  for (i in 1:(g - 1)) {
    for (j in (i + 1):g) {
      cntr = cntr + 1
      temp <- (Cov_obs[i, i] * Cov_obs[j, j] - Cov_obs[j, i] ^ 2) / 
        (Cov_obs[j, j] + Cov_obs[i, i] - 2 * Cov_obs[j, i])
      r[cntr] <- ifelse(abs(Cov_obs[i, j] / sqrt(Cov_obs[i, i] * Cov_obs[j, j])) < 1 - eps, temp, NA)
    }
  }
  r0 <- 1 / (t(rep(1, g)) %*% solve(Cov_obs) %*% rep(1, g))
  r <- c(r, diag(Cov_obs), r0)
  w <- min(r)
  Cov_new <- Cov_obs - w
  return(Cov_new)
}

proj_renorm <- function(Cov_obs, g, n, eps = 0.001) {
  temp <- eigen(Cov_obs)
  r <- rankMatrix(Cov_obs)
  Dhat <- temp$values[1:r]
  Vhat <- temp$vectors[, 1:r]
  Den <- sum((colSums(Vhat)) ^ 2 / Dhat)
  w <- 1 / Den
  Cov_new <- Cov_obs - w * (Vhat %*% t(Vhat)) %*% matrix(1, g, g) %*% (Vhat %*% t(Vhat))
  return(Cov_new)
}

new_norm <- function(Cov_obs, g, n, eps = 0.001) {
  if(g <= n) {
    Cov_new <- orig_renorm(Cov_obs, g, n, eps)
  }
  if(g > n) {
    Cov_new <- proj_renorm(Cov_obs, g, n, eps)
  }
  return(Cov_new)
}