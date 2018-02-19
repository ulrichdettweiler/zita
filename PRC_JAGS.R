library(R2jags)

# Data
HD_GCM_data <- read_excel("~/Desktop/Datasets/HD_GCM_data.xls")
HD_GCM_data <- as.data.frame((HD_GCM_data))
HD_GCM_dataNA = na.omit(HD_GCM_data)

# Centralize BN Variables
HD_GCM_dataNA$A_c = HD_GCM_dataNA$autonomy - mean(HD_GCM_dataNA$autonomy)
HD_GCM_dataNA$C_c = HD_GCM_dataNA$competence - mean(HD_GCM_dataNA$competence)
HD_GCM_dataNA$RT_c = HD_GCM_dataNA$RT - mean(HD_GCM_dataNA$RT)
HD_GCM_dataNA$RS_c = HD_GCM_dataNA$RS - mean(HD_GCM_dataNA$RS)

# Group, timepoint as numeric string
library(plyr)

HD_GCM_dataNA$grp <- revalue(HD_GCM_dataNA$Group,
                             c("int"="2", "ctr"="1"))

HD_GCM_dataNA$tpt <- revalue(HD_GCM_dataNA$TimePoint,
                             c("morning"="0", "midmorning"="1", "noon" = "2"))


HD_GCM_dataNA$ID <- revalue(as.factor(HD_GCM_dataNA$VP),
                             c("1"="1", "2"="2", "3" = "3", "4" = "4",
                               "5" = "5", "6" = "6", "7" = "7", "8" = "8",
                               "9" = "9", "10" = "10", "11" = "11", "12" = "12",
                               "13" = "13", "15" = "14", "16" = "15", "17" = "16",
                               "18" = "17", "19" = "18", "20" = "19", "21" = "20",
                               "22" = "21", "23" = "22", "24" = "23", "25" = "24",
                               "26" = "25", "27" = "26", "28" = "27", "29" = "28",
                               "30" = "29", "31" = "30", "32" = "31", "33" = "32",
                               "34" = "33", "35" = "34", "36" = "35", "37" = "36",
                               "38" = "37", "39" = "38", "40" = "39", "41" = "40",
                               "42" = "41", "43" = "42", "44" = "43", "45" = "44",
                               "46" = "45", "47" = "46", "48" = "47"))




# JAGS-Model with interaction grp x (BN + eqr)

model_PRC_HP  = '
model
{
  # Likelihood
  for (i in 1:N) {
  y[i] ~ dnorm(alpha_rGrp[grp[i]] + 
  beta_BN[grp[i],1]*x[i,1] + 
  beta_BN[grp[i],2]*x[i,2] + 
  beta_BN[grp[i],3]*x[i,3] + 
  beta_BN[grp[i],4]*x[i,4] + 
  beta_eqr[grp[i],1]*x[i,1] + 
  beta_eqr[grp[i],2]*x[i,2] + 
  beta_eqr[grp[i],3]*x[i,3] + 
  beta_gdr*x6[i],  
  sigma^-2)
  y_sim[i] ~ dnorm(alpha_rGrp[grp[i]] + 
    beta_BN[grp[i],1]*x[i,1] + 
    beta_BN[grp[i],2]*x[i,2] + 
    beta_BN[grp[i],3]*x[i,3] + 
    beta_BN[grp[i],4]*x[i,4] + 
    beta_eqr[grp[i],1]*x[i,1] + 
    beta_eqr[grp[i],2]*x[i,2] + 
    beta_eqr[grp[i],3]*x[i,3] + 
  beta_gdr*x6[i],  
  sigma^-2)
  }
  # Priors
  for (j in 1:N_grp) {
  alpha_rGrp[j] ~ dnorm(mu_alpha, sigma_alpha)
  }
  
  for (k in 1:4) {
  for (l in 1:N_grp) {
  beta_BN[l,k] ~ dnorm(0, 1)
  }
    }
  for (m in 1:3) {
  for (n in 1:N_eqr) {
  beta_eqr[m,n] ~ dnorm(0, 1)
  }
    }
  
  mu_alpha ~ dt(0, 25, 1)T(0, )
  #mu_alpha ~ dgamma(2,3)
  sigma_alpha ~ dt(0, 25, 1)T(0, )
  #sigma_alpha ~ dgamma(2,3)
  sigma ~ dt(0,25,1)T(0, )
  beta_Ssn ~ dt(0, 25, 1)T(0, )
  beta_gdr ~ dt(0, 25, 1)T(0, )
}
'

# Run it in R

library(R2jags)
JAGSrun_PRC_HP = jags(data = list(N = nrow(HD_GCM_dataNA),
                                      N_grp = length(unique(HD_GCM_dataNA$grp)),
                                      N_eqr = length(unique(HD_GCM_dataNA$enquiry)),
                                      y = log(HD_GCM_dataNA$PRC),
                                      x = as.matrix(HD_GCM_dataNA[,48:51]),
                                      grp = HD_GCM_dataNA$grp,
                                      x6 = HD_GCM_dataNA$gender),
                          parameters.to.save = c('alpha_rGrp',
                                                 'beta_BN',
                                                 'beta_eqr',
                                                 'mu_alpha',
                                                 'sigma_alpha',
                                                 'beta_gdr',
                                                 'sigma'),
                          n.iter = 10000,
                          n.thin = 5,
                          model.file = textConnection(model_PRC_HP))

# output
plot(JAGSrun_PRC_HP)
print(JAGSrun_PRC_HP)


# Model allowing for individual intercepts rather than group

model_PRC_HP_id  = '
model
{
  # Likelihood
  for (i in 1:N) {
  y[i] ~ dnorm(alpha_rID[id[i]] + 
  beta_BN[grp[i],1]*x[i,1] + 
  beta_BN[grp[i],2]*x[i,2] + 
  beta_BN[grp[i],3]*x[i,3] + 
  beta_BN[grp[i],4]*x[i,4] + 
  beta_eqr*x5[i] + 
  beta_gdr*x6[i],  
  sigma^-2)
  y_sim[i] ~ dnorm(alpha_rID[id[i]] + 
  beta_BN[grp[i],1]*x[i,1] + 
  beta_BN[grp[i],2]*x[i,2] + 
  beta_BN[grp[i],3]*x[i,3] + 
  beta_BN[grp[i],4]*x[i,4] + 
  beta_eqr*x5[i] + 
  beta_gdr*x6[i],  
  sigma^-2)
  }
  # Priors
  for (j in 1:N_id) {
  alpha_rID[j] ~ dnorm(mu_alpha, sigma_alpha)
  }
  
  for (k in 1:4) {
  for (l in 1:N_grp) {
  beta_BN[l,k] ~ dnorm(0, 1)
  }
  }
  
  beta_eqr ~ dnorm(0, 10^-2)
  mu_alpha ~ dnorm(0, 10^-2)
  sigma_alpha ~ dt(0,25,1)T(0, )
  sigma ~ dt(0,25,1)T(0, )
  beta_Ssn ~ dnorm(0, 10^-2)
  beta_gdr ~ dnorm(0, 10^-2)
}
'

# Run the model in R

JAGSrun_PRC_HP_id = jags(data = list(N = nrow(HD_GCM_dataNA),
                                      N_id = length(unique(HD_GCM_dataNA$ID)),
                                     N_grp = length(unique(HD_GCM_dataNA$grp)),
                                      y = log(HD_GCM_dataNA$PRC),
                                      #y = HD_GCM_dataNA$PRC,
                                      x = as.matrix(HD_GCM_dataNA[,48:51]),
                                      grp = HD_GCM_dataNA$grp,
                                      id = HD_GCM_dataNA$ID,
                                      x5 = HD_GCM_dataNA$enquiry,
                                      x6 = HD_GCM_dataNA$gender),
                          parameters.to.save = c('alpha_rID',
                                                 'beta_BN',
                                                 'beta_eqr',
                                                 'mu_alpha',
                                                 'sigma_alpha',
                                                 'beta_gdr',
                                                 'sigma',
                                                 'y_sim'),
                          n.iter = 50000,
                          n.thin = 10,
                          model.file = textConnection(model_PRC_HP_id))

# output, pD = 71.0 and DIC = -477.0
plot(JAGSrun_PRC_HP_id)
print(JAGSrun_PRC_HP_id)

# Posterior predictive check
pars_simID = JAGSrun_PRC_HP_id$BUGSoutput$sims.list$y_sim
plot(log(HD_GCM_dataNA$PRC), apply(pars_sim,2,'mean'))
abline(a=0, b = 1, col = 'red')

# CI's
postID = JAGSrun_PRC_HP_id$BUGSoutput$sims.matrix
apply(post,2, quantile, probs = c(0.025, 0.5, 0.975))

