#==================================================
# Import packages
#==================================================
library(nimble)



#==================================================
# Model 1 (Normal)
#==================================================
Code <- nimbleCode({
  for(i in 1:n){
    y[i] ~ dbeta(phi[i]*mu[i], phi[i]*(1 - mu[i]))
    logit(mu[i]) <- beta0 + beta1*x1[i] + beta2*x2[i] + psi[i]
    log(phi[i]) <- gamma0 + gamma1*x1[i]
    psi[i] ~ dnorm(0, var = sigma2)
  }
  
  beta0 ~ dnorm(0, var = 100)
  beta1 ~ dnorm(0, var = 100)
  beta2 ~ dnorm(0, var = 100)
  gamma0 ~ dnorm(0, var = 100)
  gamma1 ~ dnorm(0, var = 100)
  sigma2 ~ dinvgamma(0.1, 0.1)
})


Data <- list(y = y, x1 = x1, x2 = x2)
Consts <- list(n = length(y))
Inits <- list(list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10),
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1)),
              list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10),
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0.1, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1)))
mcmc.out1 <- nimbleMCMC(code = Code, data = Data, constants = Consts, inits = Inits,
                        niter = 300000, nburnin = 100000, thin = 20, nchains = 2,
                        monitors = c("beta0", "beta1", "beta2", "gamma0", "gamma1",
                                     "sigma2", "phi", "mu"),
                        summary = TRUE, WAIC = TRUE)

bayes.est1 = mcmc.out1$summary$all.chains[c(1:5, 2*length(y) + 6), ]
bayes.est1





#==================================================
# Model 2 (SN)
#==================================================
dSN <- nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma2 = double(0), 
                 alpha = double(0), log = integer(0, default = rnorm(1, 0, 10))){
    returnType(double(0))
    logProb <- log((2/sqrt(sigma2))*dnorm((x - mu)/sqrt(sigma2))*pnorm(alpha*((x - mu)/sqrt(sigma2))))
    if(log) return(logProb)
    else return(exp(logProb))
  })


rSn = function(n, mu, sigma2, alpha){
  X = rnorm(n)
  Y = rnorm(n)
  Z = rep(NA, n)
  for(i in 1:n){
    if(Y[i] <= alpha*X[i]){
      Z[i] = mu + sqrt(sigma2)*X[i]
    }else{
      Z[i] = mu - sqrt(sigma2)*X[i]
    }
  }
  return(Z)
}

rSN <- nimbleRcall(function(n = integer(0), mu = double(0), sigma2 = double(0),
                            alpha = double(0)){}, 
                   Rfun = 'rSn', returnType = double(0))


registerDistributions(list(
  dSN = list(
    BUGSdist = "dSN(mu, sigma2, alpha)",
    types = c('mu = double(0)', 'sigma2 = double(0)', 'alpha = double(0)')
  )
))


Code <- nimbleCode({
  for(i in 1:n){
    y[i] ~ dbeta(phi[i]*mu[i], phi[i]*(1 - mu[i]))
    logit(mu[i]) <- beta0 + beta1*x1[i] + beta2*x2[i] + psi[i]
    log(phi[i]) <- gamma0 + gamma1*x1[i]
    psi[i] ~ dSN(0, sigma2, alpha)
  }
  
  beta0 ~ dnorm(0, var = 100)
  beta1 ~ dnorm(0, var = 100)
  beta2 ~ dnorm(0, var = 100)
  gamma0 ~ dnorm(0, var = 100)
  gamma1 ~ dnorm(0, var = 100)
  sigma2 ~ dinvgamma(0.1, 0.1)
  alpha ~ dnorm(0, var = 100)
})


Data <- list(y = y, x1 = x1, x2 = x2)
Consts <- list(n = length(y))
Inits <- list(list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10),
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1), alpha = rnorm(1, 0, 10)),
              list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10),
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0.1, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1), alpha = rnorm(1, 0, 10)))
mcmc.out2 <- nimbleMCMC(code = Code, data = Data, constants = Consts, inits = Inits, 
                        niter = 300000, nburnin = 100000, thin = 20, nchains = 2,
                        monitors = c("beta0", "beta1", "beta2", "gamma0", "gamma1", 
                                     "sigma2", "alpha", "phi", "mu"),
                        summary = TRUE, WAIC = TRUE)

bayes.est2 = mcmc.out2$summary$all.chains[c(1:6, 2*length(y) + 7), ]
bayes.est2





#==================================================
# Model 3 (SFN)
#==================================================
dSFN <- nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma2 = double(0), alpha = double(0),
                 delta = double(0), log = integer(0, default = rnorm(1, 0, 10))){
    returnType(double(0))
    logProb <- log((1/(sqrt(sigma2) - sqrt(sigma2)*pnorm(delta)))*dnorm((abs(x - mu)/sqrt(sigma2)) + delta)*
                     pnorm(alpha*((x - mu)/sqrt(sigma2))))
    if(log) return(logProb)
    else return(exp(logProb))
  })


library(TruncatedNormal)
rSfn = function(n, mu, sigma2, alpha, delta){
  Y = S = Z = rep(NA, n)
  for(i in 1:n){
    Y[i] = rtnorm(1, -delta, sd = 1, lb = rnorm(1, 0, 10), ub = Inf)
    S[i] = sample(c(1, -1), size = 1, prob = c(pnorm(alpha*Y[i]), 1-pnorm(alpha*Y[i])))
    Z[i] = mu + sqrt(sigma2)*(S[i]*Y[i])
  }
  return(Z)
}

rSFN <- nimbleRcall(function(n = integer(0), mu = double(0), sigma2 = double(0),
                             alpha = double(0), delta = double(0)){}, 
                    Rfun = 'rSfn', returnType = double(0))


registerDistributions(list(
  dSFN = list(
    BUGSdist = "dSFN(mu, sigma2, alpha, delta)",
    types = c('mu = double(0)', 'sigma2 = double(0)', 'alpha = double(0)', 'delta = double(0)')
  )
))


Code <- nimbleCode({
  for(i in 1:n){
    y[i] ~ dbeta(phi[i]*mu[i], phi[i]*(1 - mu[i]))
    logit(mu[i]) <- beta0 + beta1*x1[i] + beta2*x2[i] + psi[i]
    log(phi[i]) <- gamma0 + gamma1*x1[i]
    psi[i] ~ dSFN(0, sigma2, alpha, delta)
  }
  
  beta0 ~ dnorm(0, var = 100)
  beta1 ~ dnorm(0, var = 100)
  beta2 ~ dnorm(0, var = 100)
  gamma0 ~ dnorm(0, var = 100)
  gamma1 ~ dnorm(0, var = 100)
  sigma2 ~ dinvgamma(0.1, 0.1)
  alpha ~ dnorm(0, var = 100)
  delta ~ dnorm(0, var = 100)
})


Data <- list(y = y, x1 = x1, x2 = x2)
Consts <- list(n = length(y))
Inits <- list(list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10),
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1), 
                   alpha = rnorm(1, 0, 10), delta = rnorm(1, 0, 10)),
              list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10),
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0.1, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1),
                   alpha = rnorm(1, 0, 10), delta = rnorm(1, 0, 10)))
mcmc.out3 <- nimbleMCMC(code = Code, data = Data, constants = Consts, inits = Inits, 
                        niter = 300000, nburnin = 100000, thin = 20, nchains = 2,
                        monitors = c("beta0", "beta1", "beta2", "gamma0", "gamma1",
                                     "sigma2", "alpha", "delta", "phi", "mu"),
                        summary = TRUE, WAIC = TRUE)

bayes.est3 = mcmc.out3$summary$all.chains[c(1:7, 2*length(y) + 8),]
bayes.est3





#==================================================
# Model 4 (BSN)
#==================================================

dBSN <- nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma2 = double(0), alpha = double(0),
                 delta = double(0), log = integer(0, default = rnorm(1, 0, 10))) {
    returnType(double(0))
    logProb <- log(2*((sigma2 + delta*((x - mu)^2))/((sqrt(sigma2)^3)*(1 + delta)))
                   *dnorm((x- mu)/sqrt(sigma2))*pnorm(alpha*((x- mu)/sqrt(sigma2))))
    if(log) return(logProb)
    else return(exp(logProb))
  })


rBsn = function(n, mu, sigma2, alpha, delta){
  X = rchisq(n, 3)
  U = sample(c(-1, 1), size = n, prob = c(0.5, 0.5), replace = TRUE)
  W = sqrt(X)*U
  Z1 = rnorm(n)
  G = sqrt(delta/(1 + delta))*W + sqrt(1/(1 + delta))*Z1
  Z2 = rnorm(n)
  T1 = rep(NA, n)
  for(i in 1:n){
    if (Z2[i] < alpha*G[i]) {
      T1[i] = mu + sqrt(sigma2)*G[i]
    } else {
      T1[i] = mu - sqrt(sigma2)*G[i]
    }
  }
  return(T1)
}

rBSN <- nimbleRcall(function(n = integer(0), mu = double(0), sigma2 = double(0),
                             alpha = double(0), delta = double(0)){}, 
                    Rfun = 'rBsn', returnType = double(0))


registerDistributions(list(
  dBSN = list(
    BUGSdist = "dBSN(mu, sigma2, alpha, delta)",
    types = c('mu = double(0)', 'sigma2 = double(0)', 'alpha = double(0)', 'delta = double(0)')
  )
))


Code <- nimbleCode({
  for(i in 1:n){
    y[i] ~ dbeta(phi[i]*mu[i], phi[i]*(1 - mu[i]))
    logit(mu[i]) <- beta0 + beta1*x1[i] + beta2*x2[i] + psi[i]
    log(phi[i]) <- gamma0 + gamma1*x1[i]
    psi[i] ~ dBSN(0, sigma2, alpha, delta)
  }
  
  beta0 ~ dnorm(0, var = 100)
  beta1 ~ dnorm(0, var = 100)
  beta2 ~ dnorm(0, var = 100)
  gamma0 ~ dnorm(0, var = 100)
  gamma1 ~ dnorm(0, var = 100)
  sigma2 ~ dinvgamma(0.1, 0.1)
  alpha ~ dnorm(0, var = 100)
  delta ~ dinvgamma(0.1, 0.1)
})


Data <- list(y = y, x1 = x1, x2 = x2)
Consts <- list(n = length(y))
Inits <- list(list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10),
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1), 
                   alpha = rnorm(1, 0, 10), delta = rinvgamma(1, 0.1, 0.1)),
              list(beta0 = rnorm(1, 0, 10), beta1 = rnorm(1, 0, 10), beta2 = rnorm(1, 0, 10), 
                   gamma0 = rnorm(1, 0, 10), gamma1 = rnorm(1, 0, 10),
                   psi = rep(0.1, Consts$n), sigma2 = rinvgamma(1, 0.1, 0.1), 
                   alpha = rnorm(1, 0, 10), delta = rinvgamma(1, 0.1, 0.1)))
mcmc.out4 <- nimbleMCMC(code = Code, data = Data, constants = Consts, inits = Inits, 
                        niter = 300000, nburnin = 100000, thin = 20, nchains = 2,
                        monitors = c("beta0", "beta1", "beta2", "gamma0", "gamma1", 
                                     "sigma2", "alpha", "delta", "phi", "mu"),
                        summary = TRUE, WAIC = TRUE)

bayes.est4 = mcmc.out4$summary$all.chains[c(1:7, 2*length(y) + 8),]
bayes.est4