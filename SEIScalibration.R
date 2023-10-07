
#DOXY-PEP model calibration
# E Reichert, 02.2023

#Model Calibration
#see what parameters lead to under status quo conditions

# Change to your file path for folder where all doxy-PEP code is stored
filepath <- "~/Documents/2021 Grad lab research/DOXYPEP/"

source(paste0(filepath, "DoxyPEPcode/DOXYfunctions.R"))

#load necessary packages
library(deSolve)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

######## Calibrate a ceftriaxone-only model to 3.0% GC prevalence #######

#SET PARAM VALUES
#parameters that will stay constant throughout model
pop = 10^6                                           #pop size
pop.p = c(0.3, 0.6, 0.1)                             #relative size of each risk group; low, M, high

rho = 1/(20*365)     #model entry/exit rate
omega_a = 10^-8      #Pr of emergence of resistance on treatment with A (ceftriaxone)
prA = 1              #Pr of treatment with A
fA = 0.98            #relative fitness, resistant to A
pi_s = 0.9           #Pr of retreatment if initial treatment failure, symptomatic
resA = 0.0001        #initial prev of resistance to A
gamma = 1/1          #rate of removal from exposure compartment (1/day)

#Set model duration + initial conditions
years = 2
tstep = 1 #in days

#set N for each sexual risk group (1 = low risk, 2 = med risk, 3 = hi risk)
N <- pop
N1 <- N*pop.p[1]
N2 <- N*pop.p[2]
N3 <- N*pop.p[3]

#distribute GC cases to have overall 3% prevalence
x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)

gc_lo <- round(N1*x*0.029,0)
gc_md <- round(N2*x*0.154,0)
gc_hi <- round(N3*x*0.817,0)

#estimated proportion of symp. infections at cross-sectional point in time
prev_symp <- 0.089
Y_gc_lo <- gc_lo*prev_symp
Z_gc_lo <- gc_lo*(1-prev_symp)
Y_gc_md <- gc_md*prev_symp
Z_gc_md <- gc_md*(1-prev_symp)
Y_gc_hi <- gc_hi*prev_symp
Z_gc_hi <- gc_hi*(1-prev_symp)

Y0 <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * (1-resA)
Z0 <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * (1-resA)

YA <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resA
ZA <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resA

inits <- c(S1 = N1-gc_lo, S2 = N2-gc_md, S3 = N3-gc_hi, 
           E01 = 0, E02 = 0, E03 = 0,
           Y01 = Y0[1], Y02 = Y0[2], Y03 = Y0[3],
           Z01 = Z0[1], Z02 = Z0[2], Z03 = Z0[3], 
           Ea1 = 0, Ea2 = 0, Ea3 = 0,
           Ya1 = YA[1], Ya2 = YA[2], Ya3 = YA[3], 
           Za1 = ZA[1], Za2 = ZA[2], Za3 = ZA[3]) 

#create vector w/ all timepoints
dt <- seq(0, 365*years, tstep) 

# Calibration Model (Drug A Only)
calibration.SI <- function(t, x, parms){
  with(as.list(c(t, x, parms)),{
    N = c(N1, N2, N3) #total population
    S = c(S1, S2, S3) #not infected
    E0 = c(E01, E02, E03) #exposed, no resistance
    Ea = c(Ea1, Ea2, Ea3) #exposed, resistant to A
    Y0 = c(Y01, Y02, Y03) #symptomatic infected, no resistance
    Ya = c(Ya1, Ya2, Ya3) #symptomatic infected, resistant to A
    Z0 = c(Z01, Z02, Z03) #asymptomatic infected, no resistance
    Za = c(Za1, Za2, Za3) #asymptomatic infected, resistant to A
    #makes 3x3 contact matrix
    activities <- c(1*c_min/365, 5*c_min/365, 20*c_min/365)
    beta <- (1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + 
      epsilon*activities/(pop*pop.p)*diag(3)
    beta <- beta * b #contacts * transmission pr per partnership
    #susceptibles
    dS <- -beta %*% ((Y0 + Z0) + fA*(Ya+Za)) *S + 
      (1-omega_a*prA)*(Ts*Y0 + Tm*Z0) +
      pi_s*Tsr*Ya +
      g*(Y0 + Ya + Z0 + Za) + 
      rho*(S + Y0 + Ya + Z0 + Za) - rho*S
    #exposed + infected w/ no resistance
    dE0 <- (beta %*% (Y0+Z0)*S) - gamma*E0
    dY0 <- sigma*gamma*E0 - Ts*Y0 - g*Y0 - rho*Y0
    dZ0 <- (1-sigma)*gamma*E0- Tm*Z0 - g*Z0 - rho*Z0
    #exposed + infected w/ resistance to A
    dEa <- fA*(beta %*% (Ya+Za)*S) - gamma*Ea
    dYa <- sigma*gamma*Ea +
      omega_a*prA*Ts*Y0 - 
      pi_s*Tsr*Ya -
      g*Ya - rho*Ya
    dZa <- (1-sigma)*gamma*Ea +
      omega_a*prA*Tm*Z0 -
      g*Za - rho*Za
    der <- c(dS, dE0, dY0, dZ0, dEa, dYa, dZa)
    list(der)
  })
}


##########################################
### model parameters from MLE fitting  ###
##########################################

#start with parameter value estimates from the literature
# (Reichert et al. 2023 and other sources listed in Supp Table 1)
theta <- c(logit.b=logit(0.5),        #transmission probability
           log.c_min=log(1.22),       #min rate of partner change
           logit.epsilon=logit(0.24), #epsilon
           logit.sigma = logit(0.5),  #pr of incident symptomatic infection
           log.Ts=log(0.04*365),      #duration of infectiousness if symptomatic (days) = average time to treatment
           log.g=log(.288*365),       #duration of infectiousness if asymptomatic and untreated (days)
           logit.Tm =logit(0.4))      #screening rate per year


### set parameters ###
b <- ilogit(theta["logit.b"])
c_min <- exp(theta["log.c_min"]) 
epsilon <- ilogit(theta["logit.epsilon"])
sigma <- ilogit(theta["logit.sigma"])
Ts <- exp(theta["log.Ts"])/365
g <- exp(theta["log.g"])/365
Tm <- ilogit(theta["logit.Tm"])


#######################################
#### Maximum likelihood estimation ####
#######################################
#used for estimating parameters

model.epi.loglik <- function(theta) {
  b <- ilogit(theta["logit.b"])
  c_min <- exp(theta["log.c_min"]) 
  epsilon <- ilogit(theta["logit.epsilon"])
  sigma <- ilogit(theta["logit.sigma"])
  Ts <- 1/(exp(theta["log.Ts"]))
  g <- 1/(exp(theta["log.g"]))
  Tm <- ilogit(theta["logit.Tm"])/365
  params <-list(c_min = c_min, epsilon = epsilon, sigma = sigma, b=b,Ts = Ts,Tm = Tm, g = g, Tsr = Ts/3)
  pred <- pred_fun_er(params)
  beta.params.prev <- estBetaParams(mu = pred, var = 1.47e-5)
  ll <- sum(dbeta(x= 0.03, beta.params.prev$alpha, beta.params.prev$beta, log=TRUE)) #calculate likelihood
  ll[is.na(ll)]<-(-1e20)
  print(ll)
  c(prev=pred,ll=ll)
}


f.optim <-function(theta) {
  Res <-  model.epi.loglik(theta)
  LogLL <- Res["ll"] #Model is returning the Log likelihood
  return(-LogLL)
}

#### run MLE optimization ###
library(bbmle)
values.start=theta
parnames(f.optim)<-names(values.start)
fit0 <- bbmle::mle2(f.optim, start=values.start,  vecpar=TRUE,  optimizer="optim"); fit0
fit <- mle2(f.optim, start=coef(fit0),  vecpar=TRUE, optimizer="optim"); fit
theta.fit<-coef(fit)

exp(theta.fit)
ilogit(theta.fit)
################################################

#Check calibration model fit with these params

#SET PARAM VALUES
#parameters that will stay constant throughout model
pop = 10^6                      #pop size
pop.p = c(0.3, 0.6, 0.1)        #relative size of each risk group; low, M, high
c_min = exp(theta.fit[2])
activities = c(1*c_min/365, 5*c_min/365, 20*c_min/365)  #sexual contacts per day
epsilon =  ilogit(theta.fit[3])                                    #mixing parameter

b = ilogit(theta.fit[1])        #transmission Pr per partnership
sigma = ilogit(theta.fit[4])    #Pr of symptomatic infection
g = 1/(exp(theta.fit[6]))  #natural recovery rate from infection
Ts =1/(exp(theta.fit[5]))   #time to treatment for symptomatic infection
Tm = ilogit(theta.fit[7])/365   #screening rate (time to treatment for asymptomatic infection)
rho = 1/(20*365)     #model entry/exit rate
omega_a = 10^-8      #Pr of emergence of resistance on treatment with A (ceftriaxone)
prA = 1              #Pr of treatment with A
fA = 0.98            #relative fitness, resistant to A
pi_s = 0.90          #Pr of retreatment if initial treatment failure, symptomatic
Tsr = Ts/3           #time to retreatment for symptomatic infection, if failure
resA = 0.0001        #initial prev of resistance to A
gamma = 1            #rate of removal from exposure compartment (1/day)

#Set model duration + initial conditions
years = 2
tstep = 1 #in days

#set N for each sexual risk group (1 = low risk, 2 = med risk, 3 = hi risk)
N <- pop
N1 <- N*pop.p[1]
N2 <- N*pop.p[2]
N3 <- N*pop.p[3]

#distribute GC cases to have overall 3% prevalence
x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)

gc_lo <- round(N1*x*0.029,0)
gc_md <- round(N2*x*0.154,0)
gc_hi <- round(N3*x*0.817,0)

#estimated proportion of symp. infections at cross-sectional point in time
prev_symp <- 0.089
Y_gc_lo <- gc_lo*prev_symp
Z_gc_lo <- gc_lo*(1-prev_symp)
Y_gc_md <- gc_md*prev_symp
Z_gc_md <- gc_md*(1-prev_symp)
Y_gc_hi <- gc_hi*prev_symp
Z_gc_hi <- gc_hi*(1-prev_symp)

Y0 <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * (1-resA)
Z0 <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * (1-resA)

YA <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resA
ZA <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resA

#start with overall prevalence of 3% gonorrhea
inits <- c(S1 = N1-gc_lo, S2 = N2-gc_md, S3 = N3-gc_hi, 
           E01 = 0, E02 = 0, E03 = 0,
           Y01 = Y0[1], Y02 = Y0[2], Y03 = Y0[3],
           Z01 = Z0[1], Z02 = Z0[2], Z03 = Z0[3], 
           Ea1 = 0, Ea2 = 0, Ea3 = 0,
           Ya1 = YA[1], Ya2 = YA[2], Ya3 = YA[3], 
           Za1 = ZA[1], Za2 = ZA[2], Za3 = ZA[3]) 

#create vector w/ all timepoints
dt <- seq(0, 365*years, tstep) 

#list out parameters
parms <- list(g=g, pop=pop, pop.p=pop.p, epsilon=epsilon, b = b, 
              sigma = sigma, Ts = Ts, Tm = Tm, rho = rho, prA = prA, 
              omega_a = omega_a, fA = fA, pi_s = pi_s, Tsr = Tsr, c_min = c_min, gamma = gamma)

#Run the model
calibration_sim <- as.data.frame(ode(inits, dt, calibration.SI, parms=parms))   

#ensure that total pop size does not change
total_check <- calibration_sim %>% mutate(RowSum=rowSums(.[setdiff(names(.),"time")])) %>% select(RowSum)

#Data tidying
calibration_sim_long <- SI.clean(calibration_sim)

# Visualize some outputs
SI.visualize(data = calibration_sim_long)

end <- calibration_sim[nrow(calibration_sim),]

#Overall prevalence of gonorrhea at t=end years
prev_GC_calibration <- sum(end[,8:22])/sum(end[,2:22])
prev_GC_calibration

#overall prevalence of resistance among cases at t=end years
prev_resistance_calibration <- sum(end[,17:22])/(sum(end[,8:22]))
prev_resistance_calibration

#cross-section prev. of exposed at t=end years (for initiating transmission model)
(end$E01 + end$Ea1)/(end$E01 + end$Ea1 + end$Y01 + end$Z01+ end$Ya1 + end$Za1)
(end$E02 + end$Ea2)/(end$E02 + end$Ea2 + end$Y02 + end$Z02+ end$Ya2 + end$Za2)
(end$E03 + end$Ea3)/(end$E03 + end$Ea3 + end$Y03 + end$Z03+ end$Ya3 + end$Za3)

#cross-section prev. of symptomatic at t=end years (for initiating transmission model)
(end$Y01+end$Y02+end$Y03+end$Ya1+end$Ya2+end$Ya3)/
  (end$Y01+end$Y02+end$Y03+end$Ya1+end$Ya2+end$Ya3 + end$Z01+end$Z02+end$Z03+end$Za1+end$Za2+end$Za3)

#view prevalence distribution over 2-yr calibration period
calibration_sim <- calibration_sim %>%
  mutate(prev_GC = (Y01 + Y02 + Y03 + Z01 + Z02 + Z03 + Ya1 + Ya2 + Ya3 + Za1 + Za2 + Za3)/10^6)

ggplot(data = calibration_sim, aes(x = time, y = prev_GC)) + geom_point()
hist(calibration_sim$prev_GC)

######## APPENDIX -- explore params for prevalence distribution #######
#define range
p = seq(0, 0.10, length=100)
#create plot of Beta distribution with shape parameters 2 and 10
vars <- estBetaParams(mu = 0.03, var = 1.47e-5)
plot(p, dbeta(p, vars$alpha, vars$beta), type='l')
plot(p, pbeta(p, vars$alpha, vars$beta), type='l')

# Sample size
n = 1000
# Parameters of the beta distribution
alpha = vars$alpha
beta = vars$beta
# Simulate some data
set.seed(1)
x = rbeta(n, alpha, beta)

# Note that the distribution is not symmetrical
curve(dbeta(x,alpha,beta))


vars$alpha/(vars$alpha + vars$beta)

# Negative log likelihood for the beta distribution
nloglikbeta = function(mu, sig) {
  alpha = mu^2*(1-mu)/sig^2-mu
  beta = alpha*(1/mu-1)
  -sum(dbeta(x, alpha, beta, log=TRUE))
}

library(stats4)
est = mle(nloglikbeta, start=list(mu=mean(x), sig=sd(x)))
confint(est)
t.test(x)$conf.int

