# Code to reproduce simulation & analysis results in Performance Evaluation (Section 5)
# of 'Bayesian sample size determination using robust commensurate priors with
# interpretable discrepancy weights'

# Lou E. Whitehead 28th May 2024, l.whitehead2@newcastle.ac.uk
 
################################################################################
# Config A: settings (calculated in samplesizes_and_priors_simstudy.R) 
# NB - replace: n , CPmu_wdash, CPvar_wdash accordingly to replicate results  for other configs
################################################################################
n = 204 # sample size for Config A

# parameters for the collective prior for Config 1:
# mu_delta ~ N(CPmu_wdash*, CPvar_wdash*)
CPmu_wdash=0.1306122
CPvar_wdash=0.4047898

################################################################################
# Generate fake data in the case of efficacy, theta = 1 ** change to 0 to reproduce results for theta = 0 
################################################################################

true.trt = c(rep(1, n/2), rep(0, n/2)) # treatment assignment
Nobs = n # number of observations
sd = 3.69 # sd at the observation level, assume known
true.mu = 1 # true treatment effect ** change to 0 to reproduce results for mu = 0 
true.g0 = 0

f=function(true.g0,true.mu,true.trt, Nobs){true.g0 + c(true.mu*true.trt) + rnorm(Nobs, mean = 0, sd = sd)}

set.seed(452)
nsims=10000
y=matrix(nrow = n, ncol = nsims)
for (i in 1:nsims){
  y[,i] <- f(true.g0,true.mu,true.trt, Nobs) # outcome data
}

################################################################################
## Analyse fake data according to Eqns (8) and (9) in manuscript (using CP(wdash) for config A)
################################################################################

R=0.5 # proportion of participants randomized to treatment
var_outcome=sd^2 # variance of outcome (assume known)
MCID = 1 # min clinically important difference

mean_tx=vector()
mean_ctl=vector()
tx_eff=vector()
post_mean=vector()
post_var=vector()
prob_mu_geq_0=vector()
prob_mu_leq_d=vector()

for (i in 1:nsims){
  
  mean_tx[i]=mean(y[,i][which(true.trt==1)]) # mean of observed data treatment group
  mean_ctl[i]=mean(y[,i][which(true.trt==0)]) # mean of observed data control group
  tx_eff[i]=mean_tx[i]-mean_ctl[i] # observed difference in means
  
  post_mean[i]=(CPvar_wdash^-1*CPmu_wdash+tx_eff[i]*n*R*(1-R)/var_outcome)/(CPvar_wdash^-1+n*R*(1-R)/var_outcome) # posterior mean
  post_var[i]=(CPvar_wdash^-1+(n*R*(1-R)/var_outcome))^-1 # posterior variance
  
  prob_mu_geq_0[i]=pnorm(post_mean[i]/sqrt(post_var[i]))# p(mu>0)
  prob_mu_leq_d[i]=pnorm((MCID-post_mean[i])/sqrt(post_var[i])) #p(mu<1)
  
}

################################################################################
# Verify pre-specified statistical properties are upheld 
# i.e., decision making is possible in 100% of trials using following criteria:
# p(mu>0) >= 0.95 -> Effective, p(mu<=delta) >= 0.80 -> Futile
################################################################################

(a=sum(prob_mu_geq_0>=0.95)/nsims) # p(Effective|mu)
(b=sum(prob_mu_leq_d[which(prob_mu_geq_0<0.95)]>=0.80)/nsims) # p(Futile|mu)
a+b # = 100% as required
