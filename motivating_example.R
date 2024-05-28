################################################################################
## Code to reproduce Motivating Example (Section 4) in 
# 'Bayesian sample size determination using robust commensurate priors with
# interpretable discrepancy weights'
#
# Lou E. Whitehead 28th May 2024, l.whitehead2@newcastle.ac.uk 
#
################################################################################

################################################################################
# First, calculate theta_q, tau_q^2 from the meta-analysis data provided in Du et al (2018):
################################################################################

######################################
# Source 1, Vreugdenhil et al., 2012 #
######################################
mean1=23.9;mean2=19 # mean1 = mean experimental group, mean2 = mean control group
(theta1 = mean1-mean2) # difference in means
n1=20;n2=20 # n1 = sample size experimental group, n1 = sample size control group
sd1=5;sd2=7.7 # sd1 = sd experimental group, sd2 = sd control group
sd_pooled1 = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2)) # pooled SD
(SEpooled1 = sd_pooled1*sqrt(1/n1+1/n2)) # pooled SE
tau1_2 = SEpooled1^2 # tau_q^2

###################################
# Source 2, Hoffmann et al., 2016 #
###################################
mean1=23.9;mean2=23.9
(theta2 = mean1-mean2)
n1=107;n2=93
sd1=3.4;sd2=3.9

sd_pooled2 = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))

SEpooled2 = sd_pooled2*sqrt(1/n1+1/n2)
tau2_2 = SEpooled2^2

#####################################
# Source 3, Ventuerlli et al., 2011 #
#####################################
mean1=12;mean2=6
(theta3 = mean1-mean2)
n1=11;n2=10
sd1=2;sd2=2

sd_pooled3 = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))

SEpooled3 = sd_pooled3*sqrt(1/n1+1/n2)
tau3_2 = SEpooled3^2 

##############################
# Source 4, Dky et al., 2008 #
##############################
m1=17.4;m2=19.2
(theta4 = mean1-mean2)
n1=24;n2=28
sd1=5.7;sd2=4.2

sd_pooled4 = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))

SEpooled4 = sd_pooled4*sqrt(1/n1+1/n2)
tau4_2 = SEpooled4^2 

#######################
# Source 5, Yang 2015 #
#######################
mean1=22.83;mean2=19.54
theta5 = mean1-mean2
n1=25;n2=25
sd1=2.75;sd2=3.43

sd_pooled5 = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))

SEpooled5 = sd_pooled5*sqrt(1/n1+1/n2)
tau5_2 = SEpooled5^2

###################################
# Source 6, Holthoff et al., 2015 #
###################################
mean1=22.11;mean2=20.72
theta6 = mean1-mean2
n1=15;n2=15
sd1=0.57;sd2=0.55

sd_pooled6 = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))

SEpooled6 = sd_pooled6*sqrt(1/n1+1/n2)
tau6_2 = SEpooled6^2

###############################
# Source 7, Kwak et al., 2008 #
###############################
mean1=19.1;mean2=12.3
theta7 = mean1-mean2
n1=15;n2=15
sd1=6.5;sd2=6.7

sd_pooled7 = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))

SEpooled7 = sd_pooled7*sqrt(1/n1+1/n2)
tau7_2 = SEpooled7^2

################################################################################
# Frequentist sample size calculation
################################################################################
(sample_size_frequentist = 3.69^2/(0.5*(1-0.5))*(qnorm(0.95)+qnorm(0.80))^2/1)

################################################################################
# Bayesian sample size function
################################################################################
ss.fun <- function(w,tau2){
  
  R = 0.5; # proportion of trial participants allocated to the new treatment
  dw =c(1.001, 1.01); br =c(1e6, 1); # downweighting / borrowing parameters for the Gamma mixture prior
  targEff = 1; # target MCID
  eta = 0.95;zeta = 0.80; # posterior decision thresholds
  newtrial_sig02=3.69^2; # assumed variance in outcomes in the new trial
  
  a=(qnorm(eta) + qnorm(zeta))^2/targEff^2
  b=sum(1/(tau2+w*dw[2]/(dw[1]-1) + (1-w)*br[2]/(br[1]-1)))
  c=newtrial_sig02/(R*(1-R))
  
  n=(a-b)*c # minimum sample size required in the new trial
}

################################################################################
# Summary prior historical information from all 7 studies
################################################################################
w = c(0.65, 0.90, 0.75, 0.75, 0.40, 0.95, 0.50)  # expert elicited w_q values for discounting
tau2 = c(tau1_2, tau2_2, tau3_2, tau4_2, tau5_2, tau6_2, tau7_2)
theta = c(theta1, theta2, theta3, theta4, theta5, theta6, theta7)

################################################################################
# Process to transform w -> wdash
################################################################################

m1=0.001/1.01;m2=(1e6-1)/1 # m1 = (dw[1]-1)/dw[2], m2 = (br[1]-1)/br[2]

# g = predictive precision, g=xi_q^-2 [Equation (14)]
g=function(w,tau2,m1,m2){(tau2+w/m1+(1-w)/m2)^-1} 

# h = linear interpolation on g=xi_q^-2 [Equation (15)]
h=function(w, tau2, m1, m2) {g(0,tau2,m1,m2)+w*(g(1,tau2,m1,m2)-g(0,tau2,m1,m2))}

# g_inv = inverse of g=xi_q^2 [Equation (16)]
g_inv=function(p, tau2, m1, m2){((1/p - tau2)*m1*m2 -m1)/(m2-m1)}

# Implementation (Transform w -> wdash)
p=h(w, tau2, m1, m2)
(wdash = g_inv(p, tau2, m1, m2))

################################################################################
# Bayesian sample size
################################################################################
(n=ss.fun(wdash, tau2))

################################################################################
# Formulation of CP*, mu_delta ~ N(mu_CP*, var_PP*)
################################################################################
# synthesis weights
pstar=g(wdash,tau2,m1,m2)/sum(g(wdash,tau2,m1,m2)) 

# CP*
(mu_CPwdash=sum(pstar*theta))
(var_CPwdash=1/sum((tau2+wdash/m1+(1-wdash)/m2)^-1))
