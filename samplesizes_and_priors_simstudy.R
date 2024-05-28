# Code to reproduce sample size calculations used in Performance Evaluation (Section 5)
# in 'Bayesian sample size determination using robust commensurate priors with
# interpretable discrepancy weights'
#
# Lou E. Whitehead 28th May 2024, l.whitehead2@newcastle.ac.uk
# 

################################################################################
# Proposed sample size function borrowing from Q sources
################################################################################

ss.fun <- function(w,tau2){
  
  R = 0.5; # proportion of trial participants allocated to the new treatment
  dw =c(1.01, 1.01); br =c(1e6, 1); # downweighting / borrowing parameters for the Gamma mixture prior
  targEff = 1; # target MCID
  eta = 0.95;zeta = 0.80; # posterior decision thresholds
  newtrial_sig02=3.69^2; # assumed variance in outcomes in the new trial
  
  a=(qnorm(eta) + qnorm(zeta))^2/targEff^2
  b=sum(1/(tau2+w*dw[2]/(dw[1]-1) + (1-w)*br[2]/(br[1]-1)))
  c=newtrial_sig02/(R*(1-R))
  
  n=(a-b)*c # minimum sample size required in the new trial
}

################################################################################
# Config A sample size (replace theta_q , tau2_q, w_q accordingly for other configs)
################################################################################
# Historical data, summarized as lambda_q ~ N(theta_q, tau2_q)
# Description: (relatively) small/neutral historic treatment effects with larger variances
theta=c(0.10, 0.24, 0.37, 0, -0.05) # Historical theta_q
tau2=c(1.25,0.73,0.92,1.29,0.66) # Historical tau^2_q

w=c(0.2,0.4,0.8,0.6,0.7) # w_q values for discounting

# ################################################################################
# # Config B sample size (replace theta_q , tau2_q, w_q accordingly for other configs)
# ################################################################################
# # Historical data, summarized as lambda_q ~ N(theta_q, tau2_q)
# # Description: mixed historic treatment effects , weightings favour less positive trials
# theta=c(0, -0.05, 2.14, 0.37, 1.10) # Historical means, theta_q
# tau2=c(1.29,0.66,0.50,0.92, 0.75) # Historical variances , tau^2_q
# 
# w=c(0.2,0.4,0.8,0.6,0.7) # w_q values for discounting
# 
# ################################################################################
# # Config C sample size (replace theta_q , tau2_q, w_q accordingly for other configs)
# ################################################################################
# # Historical data, summarized as lambda_q ~ N(theta_q, tau2_q)
# # Description: mixed historic treatment effects , weightings favour more positive trials
# theta=c(1.10, 0.37, -0.05,2.14,0) # Historical means, theta_q
# tau2=c(0.75,0.92,0.66, 0.50,1.29) # Historical variances , tau^2_q
# 
# w=c(0.2,0.4,0.8,0.6,0.7) # w_q values for discounting
# 
# ################################################################################
# # Config D sample size (replace theta_q , tau2_q, w_q accordingly for other configs)
# ################################################################################
# # Historical data, summarized as lambda_q ~ N(theta_q, tau2_q)
# # Description: more positive historic treatment effects, lower variances
# theta=c(1.10, 2.14, 1.07, 0.6, 0.85) # Historical means, theta_q
# tau2=c(0.75,0.50,0.82,0.89,0.26) # Historical variances , tau^2_q
# 
# w=c(0.2,0.4,0.8,0.6,0.7) # w_q values for discounting

################################################################################
# Process to transform w -> wdash
################################################################################

m1=0.01/1.01;m2=(1e6-1)/1 # m1 = (dw[1]-1)/dw[2], m2 = (br[1]-1)/br[2]

# g = predictive precision, g=xi_q^-2 [Equation (14)]
g=function(w,tau2,m1,m2){(tau2+w/m1+(1-w)/m2)^-1} 

# h = linear interpolation on g=xi_q^-2 [Equation (15)]
h=function(w, tau2, m1, m2) {g(0,tau2,m1,m2)+w*(g(1,tau2,m1,m2)-g(0,tau2,m1,m2))}

# g_inv = inverse of g=xi_q^2 [Equation (16)]
g_inv=function(p, tau2, m1, m2){((1/p - tau2)*m1*m2 -m1)/(m2-m1)}

# Implementation (Transform w -> wdash)
p=h(w, tau2, m1, m2)
wdash = g_inv(p, tau2, m1, m2)

# Sample size for configuration A
(n=ss.fun(wdash, tau2))

################################################################################
# Calculation of CP*, mu_delta ~ N(mu_CP*, var_CP*)
# (used for analysis, see analysis_sim_study.R)
################################################################################
# synthesis weights
pstar=g(wdash,tau2,m1,m2)/sum(g(wdash,tau2,m1,m2)) 

# CP* used for analysis
(mu_CPwdash=sum(pstar*theta))
(var_CPwdash=1/sum((tau2+wdash/m1+(1-wdash)/m2)^-1))

