#==================================================
# Import packages
#==================================================
library(openxlsx)



#==================================================
# Generating covariate
#==================================================
n = 100

x1 = rnorm(n)
hist(x1, nclass = 50, freq = FALSE)
lines(density(x1), col = "red", lwd = 2)

x2 = rnorm(n)
hist(x2, nclass = 50, freq = FALSE)
lines(density(x2), col = "red", lwd = 2)



#==================================================
# Standard Normal (Scenario 1)
#==================================================
psi = rnorm(n)
hist(psi, nclass = 50, freq = FALSE)
lines(density(psi), col = "red", lwd = 2)

beta0 = 0.1
beta1 = 0.3
beta2 = 0.2
logit_mu = beta0 + beta1*x1 + beta2*x2 + psi
mu = exp(logit_mu) / (1 + exp(logit_mu))

gamma0 = 5
gamma1 = 2
ln.phi = gamma0 + gamma1*x1
phi = exp(ln.phi)

y1 = c()
for (i in 1:n) {
  y1[i] = rbeta(1, phi[i]*mu[i], phi[i]*(1 - mu[i]))
}
hist(y1, nclass = 50, labels = T, freq = F)
lines(density(y1), col = "red", lwd = 2)


# Save data
data1 <- data.frame(y = y1, x1 = x1, x2 = x2)
write.xlsx(data1, file = "path", rowNames = FALSE)



#==================================================
# Mixture Normal (Scenario 2)
#==================================================
u = runif(n)
psi = c()
for(i in 1:n){
  if(u[i] < 0.3){
    psi[i] = rnorm(1, -2, 1)
  }else{
    psi[i] = rnorm(1, 2, 1)
  }
}
hist(psi, nclass = 50, freq = FALSE)
lines(density(psi), col = "red", lwd = 2)

beta0 = 0.1
beta1 = 0.3
beta2 = 0.2
logit_mu = beta0 + beta1*x1 + beta2*x2 + psi
mu = exp(logit_mu) / (1 + exp(logit_mu))

gamma0 = 5
gamma1 = 2
ln.phi = gamma0 + gamma1*x1
phi = exp(ln.phi)

y2 = c()
for (i in 1:n) {
  y1[i] = rbeta(1, phi[i]*mu[i], phi[i]*(1 - mu[i]))
}
hist(y2, nclass = 50, labels = T, freq = F)
lines(density(y2), col = "red", lwd = 2)


# Save data
data1 <- data.frame(y = y2, x1 = x1, x2 = x2)
write.xlsx(data1, file = "Path", rowNames = FALSE)