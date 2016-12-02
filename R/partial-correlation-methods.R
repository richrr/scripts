

# this implements different methods for partial correlations


## ggmridge

library(GGMridge)
p <- 20 # number of variables
n <- 10 # sample size
###############################
# Simulate data
###############################
simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
data <- simulation$data[[1]]
stddata <- scale(x = data, center = TRUE, scale = TRUE)
###############################
# estimate ridge parameter
###############################
lambda.array <- seq(from = 0.1, to = 20, by = 0.1) * (n - 1.0)
fit <- lambda.cv(x = stddata, lambda = lambda.array, fold = 10L)
lambda <- fit$lambda[which.min(fit$spe)] / (n - 1.0)
###############################
# calculate partial correlation
# using ridge inverse
###############################
w.upper <- which(upper.tri(diag(p)))
partial <- solve(lambda * diag(p) + cor(data))
partial <- (-scaledMat(x = partial))[w.upper]
###############################
# get p-values from empirical
# null distribution
###############################
efron.fit <- getEfronp(z = transFisher(x = partial))


###############################
# Simulate data
###############################
simulation <- simulateData(G = p, etaA = 0.02, n = n, r = 1)
dat <- simulation$data[[1L]]
stddat <- scale(x = dat, center = TRUE, scale = TRUE)
shrinkage.lambda <- lambda.TargetD(x = stddat)
###############################
# the ridge parameter
###############################
ridge.lambda <- shrinkage.lambda / (1.0 - shrinkage.lambda)
###############################
# partial correlation matrix
###############################
partial <- solve(cor(dat) + ridge.lambda * diag(ncol(dat)))
partial <- -scaledMat(x = partial)
###############################
# Fisher s Z transformation of
# upper diagonal of the partial
# correlation matrix
###############################
w.upper <- which(upper.tri(diag(nrow(dat))))
psi <- transFisher(x = partial[w.upper])
