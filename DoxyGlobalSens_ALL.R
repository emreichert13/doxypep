
# E Reichert 02.2024
# 
# Doxy-PEP Model Impact
# with doxy-PEP targeted to entire population (regardless of risk or sexual activity group)

# GLOBAL Sensitivity Analysis - Sample select model parameters from distributions

#load necessary packages
library(tidyverse)
library(deSolve)
library(RColorBrewer)
library(showtext)
library(reshape2)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Enter file path for folder where doxy-PEP code is stored
filepath <- "~/Documents/2021 Grad lab research/DOXYPEP/"

source(paste0(filepath, "DoxyPEPcode/Final_02.24/DOXYfunctions.R"))

# Set parameters that are fixed
pop = 10^6                    #pop size
pop.p = c(0.3, 0.6, 0.1)      #relative size of each risk group; low, M, high
rho = 1/(20*365)     #model entry/exit rate
prA = 1              #Pr of treatment with A
prB = 0              #Pr of DOXYPEP (vary this coverage %)
pi_s = 0.90          #Pr of retreatment if initial treatment failure, symptomatic
resA = 0.0001        #initial prev of resistance to A
resB = 0.109         #initial prev of resistance to doxy
gamma = 1/(1)        #rate of recovery from exposure
c_min = 1.4076374
activities = c(1*c_min/365, 
               5*c_min/365, 
               20*c_min/365)  #sexual contacts per day

# Define sampling distribution for params not fixed
vars_epsilon <- estBetaParams(mu = 0.25, var = 0.01)

# b, transmission Pr per partnership
vars_b <- estBetaParams(mu = 0.55, var = 0.01)

# sigma, Pr of symptomatic infection
vars_sigma <- estBetaParams(mu = 0.45, var = 0.01)

#g, natural recovery rate
vars_g <- estGammaParams(mu = 0.209*365.25, var = 15^2)

# Ts, treatment rate for symptomatic infection
vars_Ts <- estGammaParams(mu = 15.4, var = 3.5^2)

# Tm, screening rate for asymptomatic infection
vars_Tm <- estGammaParams(mu = 2.80*365, var = 15^2)

vars_f <- estBetaParams(mu = 0.98, var = 0.01)

vars_kappa <- estBetaParams(mu = 0.38, var = 0.01)

years = 20
tstep = 1 #in days

#set N for each sexual risk group (1 = lo risk, 2 = med risk, 3 = hi risk)
N <- pop
N1 <- N*pop.p[1]
N2 <- N*pop.p[2]
N3 <- N*pop.p[3]


#estimated prevalence of GC in each sexual activity group, @ equilibrium
# (mean from calibration model run; gives overall weighted GC prev of 3.0%)
gc_lo <- round(N1*0.00480,0)
gc_md <- round(N2*0.02576,0)
gc_hi <- round(N3*0.13719,0)

#estimated proportion of symp. infections at cross-sectional point in time, @ equilibrium
# (mean from calibration model run, 13%)
prev_symp <- 0.13
Y_gc_lo <- gc_lo*prev_symp
Z_gc_lo <- gc_lo*(1-prev_symp)
Y_gc_md <- gc_md*prev_symp
Z_gc_md <- gc_md*(1-prev_symp)
Y_gc_hi <- gc_hi*prev_symp
Z_gc_hi <- gc_hi*(1-prev_symp)

Y0 <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * (1-(resA + resB - resA*resB))
Z0 <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * (1-(resA + resB - resA*resB))

YA <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resA * (1-resB)
ZA <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resA * (1-resB)

YB <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resB * (1-resA)
ZB <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resB * (1-resA)

YAB <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * (resA*resB)
ZAB <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * (resA*resB)

#estimated proportion of individuals exposed at cross-sectional point in time, @ equilibrium
# (mean from calibration model run, 0.07%)
prev_exp <- 0.0007
exp_lo <- round(N1*prev_exp,0)
exp_md <- round(N2*prev_exp,0)
exp_hi <- round(N3*prev_exp,0)

E0 <- c(gc_lo, gc_md, gc_hi) * 0.016 * (1-(resA + resB - resA*resB))
Ea <- c(gc_lo, gc_md, gc_hi) * 0.016 * resA * (1-resB)
Eb <- c(gc_lo, gc_md, gc_hi) * 0.016 * resB * (1-resA)
Eab <- c(gc_lo, gc_md, gc_hi) * 0.016 * resA*resB

E0 <- c(exp_lo, exp_md, exp_hi) * (1-(resA + resB - resA*resB))
Ea <- c(exp_lo, exp_md, exp_hi) * resA * (1-resB)
Eb <- c(exp_lo, exp_md, exp_hi) * resB * (1-resA)
Eab<- c(exp_lo, exp_md, exp_hi) * resA * resB

#start with overall prevalence of 3% of gonorrhea
inits <- c(S1 = N1-gc_lo-exp_lo, S2 = N2-gc_md-exp_md, S3 = N3-gc_hi-exp_hi,
           E01 = E0[1], E02 = E0[2], E03 = E0[3], 
           Y01 = Y0[1], Y02 = Y0[2], Y03 = Y0[3],
           Z01 = Z0[1], Z02 = Z0[2], Z03 = Z0[3],
           Ea1 = Ea[1], Ea2 = Ea[2], Ea3 = Ea[3],
           Ya1 = YA[1], Ya2 = YA[2], Ya3 = YA[3],
           Za1 = ZA[1], Za2 = ZA[2], Za3 = ZA[3],
           Eb1 = Eb[1], Eb2 = Eb[2], Eb3 = Eb[3],
           Yb1 = YB[1], Yb2 = YB[2], Yb3 = YB[3], 
           Zb1 = ZB[1], Zb2 = ZB[2], Zb3 = ZB[3],
           Eab1 = Eab[1], Eab2 = Eab[2], Eab3 = Eab[3],
           Yab1 = YAB[1], Yab2 = YAB[2], Yab3 = YAB[3],
           Zab1 = ZAB[1], Zab2 = ZAB[2], Zab3 = ZAB[3])

dt <- seq(0, 365*years, tstep)

#Set up the ODE system of equations
doxyPEP.SI <- function(t, x, parms){
  with(as.list(c(t, x, parms)),{
    N = c(N1, N2, N3) #total population
    S = c(S1, S2, S3) #not infected
    E0 = c(E01, E02, E03)
    Ea = c(Ea1, Ea2, Ea3)
    Eb = c(Eb1, Eb2, Eb3)
    Eab = c(Eab1, Eab2, Eab3)
    Y0 = c(Y01, Y02, Y03) #symptomatic infected, no resistance
    Ya = c(Ya1, Ya2, Ya3) #symptomatic infected, resistant to A
    Yb = c(Yb1, Yb2, Yb3)
    Yab = c(Yab1, Yab2, Yab3)
    Z0 = c(Z01, Z02, Z03) #asymptomatic infected, no resistance
    Za = c(Za1, Za2, Za3) #asymptomatic infected, resistant to A
    Zb = c(Zb1, Zb2, Zb3)
    Zab = c(Zab1, Zab2, Zab3)
    
    #makes 3x3 contact matrix
    activities <- c(1*c_min/365, 5*c_min/365, 20*c_min/365)
    beta <- (1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + 
      epsilon*activities/(pop*pop.p)*diag(3)
    beta <- beta * b #contacts * transmission pr per partnership
    
    #keep track of prev resistance to B
    prevA = (sum(Ya) + sum(Za) + sum(Yab) + sum(Zab))/(sum(Y0) + sum(Z0) + sum(Ya) + sum(Za) + sum(Yb)+sum(Zb)+sum(Yab)+sum(Zab))
    prevB = (sum(Yb) + sum(Zb) + sum(Yab) + sum(Zab))/(sum(Y0) + sum(Z0) + sum(Ya) + sum(Za) + sum(Yb)+sum(Zb)+sum(Yab)+sum(Zab))
    prevAB = (sum(Yab) + sum(Zab))/(sum(Y0) + sum(Z0) + sum(Ya) + sum(Za) + sum(Yb)+sum(Zb)+sum(Yab)+sum(Zab))
    prev0 = (sum(Y0) + sum(Z0))/(sum(Y0) + sum(Z0) + sum(Ya) + sum(Za) + sum(Yb)+sum(Zb)+sum(Yab)+sum(Zab))
    #keep track of inc infections each day
    Inc = sum(beta %*% ((Y0 + Z0) + fA*(Ya+Za) + fB*(Yb+Zb) + fAB*(Yab+Zab))*S)
    prevGC = (sum(Y0) + sum(Z0) + sum(Ya) + sum(Za) + sum(Yb)+sum(Zb)+sum(Yab)+sum(Zab))/10^6
    ARx = sum(Ts*Y0 + Tm*Z0 + pi_s*Tsr*Ya + Ts*Yb + Tm*Zb + pi_s*Tsr*Yab)
    BRx = sum(prB*gamma*(E0 + Ea + Eb + Eab))
    
    #susceptibles
    dS <- -beta %*% ((Y0 + Z0) + fA*(Ya+Za) + fB*(Yb+Zb) + fAB*(Yab+Zab))*S +
      (1-omega_a)*(Ts*Y0 + Tm*Z0) +
      pi_s*(Tsr*Ya + Tsr*Yab) +
      (1-omega_a)*(Ts*Yb + Tm*Zb) +
      prB*(1-kappa)*(1-omega_b)*(gamma*E0 + gamma*Ea) +
      g*(Y0 + Ya + Z0 + Za + Yb + Zb + Yab + Zab) + 
      rho*(S + Y0 + Ya + Z0 + Za + Yb + Zb + Yab + Zab) - rho*S
    #infections w/ no resistance
    dE0 <- (beta %*% (Y0+Z0)*S) - gamma*E0
    dY0 <- sigma*(1-prB)*(gamma*E0) + sigma*prB*kappa*(1-omega_b)*(gamma*E0) - 
      g*Y0 - rho*Y0 - Ts*Y0
    dZ0 <- (1-sigma)*(1-prB)*(gamma*E0) + (1-sigma)*prB*kappa*(1-omega_b)*(gamma*E0) - 
      g*Z0 - rho*Z0 - Tm*Z0
    #infections w/ resistance to A
    dEa <- fA*(beta %*% (Ya+Za)*S) - gamma*Ea
    dYa <- sigma*(1-prB)*(gamma*Ea) + 
      sigma*prB*kappa*(1-omega_b)*(gamma*Ea) +
      omega_a*Ts*Y0 -
      g*Ya - rho*Ya - 
      pi_s*Tsr*Ya
    dZa <- (1-sigma)*(1-prB)*(gamma*Ea) + 
      (1-sigma)*prB*kappa*(1-omega_b)*(gamma*Ea) + 
      omega_a*Tm*Z0 -
      g*Za - rho*Za
    #infections w/ resistance to B
    dEb <- fB*(beta %*% (Yb+Zb)*S) - gamma*Eb
    dYb <- sigma*gamma*Eb + sigma*prB*omega_b*(gamma*E0) -
      g*Yb - rho*Yb - Ts*Yb
    dZb <- (1-sigma)*gamma*Eb + (1-sigma)*prB*omega_b*(gamma*E0)-
      g*Zb - rho*Zb - Tm*Zb
    #infections w/ dual resistance
    dEab <- fAB*(beta %*% (Yab+Zab)*S) - gamma*Eab
    dYab <- sigma*gamma*Eab + sigma*prB*omega_b*(gamma*Ea) +
      omega_a*Ts*Yb - 
      pi_s*Tsr*Yab -
      g*Yab - rho*Yab
    dZab <- (1-sigma)*gamma*Eab + (1-sigma)*prB*omega_b*(gamma*Ea) +
      omega_a*Tm*Zb -
      g*Zab - rho*Zab
    der <- c(dS, dE0, dY0, dZ0, dEa, dYa, dZa, dEb, dYb, dZb, dEab, dYab, dZab)
    list(der, Inc = Inc, prevA = prevA, prevB = prevB, prevAB = prevAB, prev0 = prev0, 
         prevGC = prevGC, ARx = ARx, BRx = BRx)
  })
}

#run ODE models over all combinations of these params
RunODE <- function (x) {
  #provide column names
  parms$prB <- x
  doxy_sim <- as.data.frame(ode(inits, dt, doxyPEP.SI, parms=parms))   
  return(doxy_sim)
}

res_global <- data.frame()

for (i in 1:1000) {
  
  set.seed(i)
  # Draw parameter values
  ## DRAW PARAMS FROM DISTRIBUTIONS ##
  epsilon <- rbeta(1, vars_epsilon$alpha, vars_epsilon$beta)
  b       <- rbeta(1, vars_b$alpha, vars_b$beta)
  sigma   <- rbeta(1, vars_sigma$alpha, vars_sigma$beta)
  x       <- rgamma(1, shape = vars_g$shape, scale = vars_g$scale)
  g       <- 1/x
  x       <- rgamma(1, shape = vars_Ts$shape, scale = vars_Ts$scale)
  Ts      <- 1/x
  x       <- rgamma(1, shape = vars_Tm$shape, scale = vars_Tm$scale)
  Tm      <- 1/x

  omega_a <- runif(1, min = 0, max = 2*10^-8)
  omega_b <- runif(1, min = 0, max = 2*10^-8)
  fA      <- rbeta(1, vars_f$alpha, vars_f$beta)
  fB      <- rbeta(1, vars_f$alpha, vars_f$beta)
  fAB     <- fA*fB
  Tsr     <- Ts/3 
  kappa   <- rbeta(1, vars_kappa$alpha, vars_kappa$beta)


#list out parameters defined above
  parms <- list(g=g, 
              pop=pop, 
              pop.p=pop.p, 
              epsilon=epsilon, 
              c_min = c_min,
              b = b, 
              sigma = sigma, 
              Ts = Ts, 
              Tm = Tm, 
              rho = rho, 
              prA = prA, 
              prB = prB, 
              omega_a = omega_a, 
              omega_b = omega_b, 
              fA = fA, fB = fB, 
              fAB = fAB, 
              pi_s = pi_s, 
              Tsr = Tsr,
              gamma = gamma,
              kappa = kappa)

  # RUN ODE function for various DoxyPEP uptake levels (0 to 90%)
  temp <- lapply(X = c(0,0.1,0.25,0.5,0.75,0.90), RunODE)
  doxy_sim0 <- temp[[1]]
  doxy_sim10 <- temp[[2]]
  doxy_sim25 <- temp[[3]]
  doxy_sim50 <- temp[[4]]
  doxy_sim75 <- temp[[5]]
  doxy_sim90 <- temp[[6]]

  # drop unnecessary cols
  # doxy_sim0 <- doxy_sim0 %>% select(-Doxy1, -Doxy2, -Doxy3)
  # doxy_sim10 <- doxy_sim10 %>% select(-Doxy1, -Doxy2, -Doxy3)
  # doxy_sim25 <- doxy_sim25 %>% select(-Doxy1, -Doxy2, -Doxy3)
  # doxy_sim50 <- doxy_sim50 %>% select(-Doxy1, -Doxy2, -Doxy3)
  # doxy_sim75 <- doxy_sim75 %>% select(-Doxy1, -Doxy2, -Doxy3)
  # doxy_sim90 <- doxy_sim90 %>% select(-Doxy1, -Doxy2, -Doxy3)
  
  # add indicator col for DoxyPEP uptake level
  doxy_sim0 <- doxy_sim0 %>% mutate(DoxyPEP = "0%")
  doxy_sim10 <- doxy_sim10 %>% mutate(DoxyPEP = "10%")
  doxy_sim25 <- doxy_sim25 %>% mutate(DoxyPEP = "25%")
  doxy_sim50 <- doxy_sim50 %>% mutate(DoxyPEP = "50%")
  doxy_sim75 <- doxy_sim75 %>% mutate(DoxyPEP = "75%")
  doxy_sim90 <- doxy_sim90 %>% mutate(DoxyPEP = "90%")
  
  doxy_sim_all <- rbind(doxy_sim0, doxy_sim10, doxy_sim25, doxy_sim50, doxy_sim75, doxy_sim90)
  
  #View model outcomes by DoxyPEP utilization
  doxy_sim_all <- doxy_sim_all %>%
    group_by(DoxyPEP) %>%
    mutate(CumInc = cumsum(Inc),
           CumARx = cumsum(ARx),
           CumBRx = cumsum(BRx)) %>%
    ungroup()
  
  baseline <- doxy_sim_all %>%
    filter(DoxyPEP == "0%") %>%
    select(time, 
           prev_0 = prevGC)
  
  #output key outcomes by DoxyPEP uptake level
  run_summary <- left_join(doxy_sim_all, baseline, by = "time") %>%
    group_by(DoxyPEP) %>%
    summarise(CI5 = CumInc[time == 1825],
              CI20 = CumInc[time == 7300],
              Prev5 = prevGC[time == 1825],
              Prev20 = prevGC[time == 7300],
              MinPR = min(prevGC/prev_0),
              TimeA = min(time[prevA >= 0.05])/365,
              TimeAB = min(time[prevAB >= 0.05])/365,
              TimeB = min(time[prevB >= 0.84])/365,
              ARx_5 = CumARx[time == 1825],
              ARx_20 = CumARx[time == 7300],
              BRx_5 = CumBRx[time == 1825],
              BRx_20 = CumBRx[time == 7300]) %>%
    mutate(Rel_CI5 = CI5/CI5[DoxyPEP == "0%"],
           Rel_CI20 = CI20/CI20[DoxyPEP == "0%"],
           Rel_PR5 = Prev5/Prev5[DoxyPEP == "0%"],
           Rel_PR20 = Prev20/Prev20[DoxyPEP == "0%"],
           Rel_A5 = ARx_5/ARx_5[DoxyPEP == "0%"],
           Rel_A20 = ARx_20/ARx_20[DoxyPEP == "0%"],
           iter = i)
  
  res_global <- rbind(res_global, run_summary)
  print(i)
}

write.csv(res_global, paste0(filepath, "GlobalSensResults_ALL.csv"), row.names=FALSE)
#res_global <- read.csv(paste0(filepath, "GlobalSensResults_ALL.csv"))

res_global %>%
  group_by(DoxyPEP) %>%
  summarise(TimeA_inf = sum(TimeA == Inf)/1000,
            TimeAB_inf = sum(TimeAB == Inf)/1000,
            TimeB_inf = sum(TimeB == Inf)/1000)

# replace Inf values with value > 20 yrs - will leave a note in table
res_global <- res_global %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 7301/365, x))

res_global %>%
  group_by(DoxyPEP) %>%
  summarise_each(funs(median), 
                 Prev5, Prev20, Rel_PR5, Rel_PR20, MinPR,
                 TimeA, TimeAB, TimeB, 
                 Rel_CI5, Rel_CI20, 
                 Rel_A5, Rel_A20,
                 BRx_5, BRx_20) %>%
  ungroup

res_global %>%
  group_by(DoxyPEP) %>%
  summarise_each(funs(quantile(., 0.25)), TimeA, TimeAB, TimeB, Rel_CI5, Rel_CI20, Rel_A5, Rel_A20) %>%
  ungroup

res_global %>%
  group_by(DoxyPEP) %>%
  summarise_each(funs(quantile(., 0.75)), TimeA, TimeAB, TimeB, Rel_CI5, Rel_CI20, Rel_A5, Rel_A20) %>%
  ungroup


