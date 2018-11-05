# What relative importance do the individual motivational domains have on PO?
# Model Individuals nested in Classes , with "alternative centering"
# Bayesian P-Value = 0.53
# DIC info: (pD = var(deviance)/2) 
# pD = 21.9 and DIC = 1187.722 

library(R2jags)
library(jagsUI)
library(rjags)

library(readr)

library(ggplot2)
library(gridExtra)
library(gridBase)
library(MCMCvis)
library(coda)

JAGS_Model_PO_class = '
model
{
  # Likelihood
  for (i in 1:N) {
  y[i] ~ dnorm(alphaClass[class[i]] + 
  betaMD[Cond[i],1]*x[i,1] + 
  betaMD[Cond[i],2]*x[i,2] + 
  betaMD[Cond[i],3]*x[i,3] + 
  betaMD[Cond[i],4]*x[i,4] + 
  betaGend*x4[i],
  sigma^-2)
  y_sim[i] ~ dnorm(alphaClass[class[i]] + 
  betaMD[Cond[i],1]*x[i,1] + 
  betaMD[Cond[i],2]*x[i,2] + 
  betaMD[Cond[i],3]*x[i,3] + 
  betaMD[Cond[i],4]*x[i,4] + 
  betaGend*x4[i],
  sigma^-2)
  
  res[i] <- y[i] - y_sim[i]   
  emp.new[i] ~ dnorm(y_sim[i], sigma^-2)
  res.new[i] <- emp.new[i] - y_sim[i]
  }
  # Priors
  for (j in 1:N_Class) {
  alphaClass[j] ~ dnorm(0, sigma_alpha^-2)
  
  }
  
  for (k in 1:4) {
  for (l in 1:N_Cond) {
  betaMD[l,k] ~ dnorm(mu_alpha, 1)
  
  }
  }
  
  mu_alpha ~ dnorm(0, 2^-2)
  sigma_alpha ~ dt(0, 5, 1)T(0, )
  sigma ~ dnorm(0, 1)T(0, )
  betaGend ~ dnorm(0, 1)
  
  #Derived parameters
  fit <- sum(res[])
  fit.new <- sum(res.new[])
}
'



JAGS_Run_PO_C_class = jags(data = list(N = nrow(Relevant_Long),
                                       N_Cond = length(unique(Relevant_Long$Cond)),
                                       N_Class = length(unique(Relevant_Long$Class)),
                                       y = Relevant_Long$PO,
                                       x = as.matrix(Relevant_Long[,c(12,13,14,15)]),
                                       Cond = Relevant_Long$Cond,
                                       class = Relevant_Long$Class,
                                       x4 = Relevant_Long$Gender),
                           parameters.to.save = c('betaMD',
                                                  'mu_alpha',
                                                  'sigma_alpha',
                                                  'betaGend',
                                                  'sigma',
                                                  'fit',
                                                  'fit.new'),
                           n.iter = 50000,
                           n.thin = 10,
                           n.chains = 3,
                           n.burnin = 25000,
                           model.file = textConnection(JAGS_Model_PO_class))

plot(JAGS_Run_PO_C_class)
print(JAGS_Run_PO_C_class)


library(jagsUI)
pp.check(JAGS_Run_PO_C_class, actual = 'fit', new = 'fit.new')

MCMCplot(JAGS_Run_PO_C_class, 
         params = c('betaMD', 'betaGend'),
         xlim = c(-1, 1),
         horiz = TRUE,
         rank = FALSE,
         ref_ovl = TRUE,
         xlab = 'Parameter Estimate', 
         main = 'Perceived Relevance of Contents', 
         labels = c('InRc Indoor', 'InRc Outdoor', 'IdRc Indoor', 
                    'IdRc Outdoor', 'IjRc Indoor', 'IjRc Outdoor', 'ExRc Indoor',
                    'ExRc Outdoor', 'Gender'),
         labels_sz = 1, med_sz = 2, thick_sz = 4, thin_sz = 3, 
         ax_sz = 3, main_text_sz = 1.5)
