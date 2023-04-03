## List of Functions for DoxyPEP project ##
# E Reichert, 03.2023

# clean output of ODE solver into long format
SI.clean <- function(data){
  
  data_long <- data %>%
    gather(key = type, value = individuals, -time)
  
  #Classify infections by sexual risk group
  data_long <- data_long %>%
    mutate(RiskGroup = ifelse(type == "S1" | type == "Y01" | type == "Ya1" | type == "Yb1" | 
                                type == "Yab1" | type == "Z01" | type == "Za1" | type == "Zb1" |
                                type == "Zab1" | type == "E01" | type == "Ea1" | type == "Eb1" | type == "Eab1",
                              "Low",
                              ifelse(type == "S2" | type == "Y02" | type == "Ya2" | type == "Yb2" | 
                                       type == "Yab2" | type == "Z02" | type == "Za2" | type == "Zb2" | 
                                       type == "Zab2" | type == "E02" | type == "Ea2" | type == "Eb2" | type == "Eab2", 
                                     "Intermediate", 
                                     ifelse(type == "S3" | type == "Y03" | type == "Ya3" | type == "Yb3" | 
                                              type == "Yab3" | type == "Z03" | type == "Za3" | type == "Zb3" | 
                                              type == "Zab3" | type == "E03" | type == "Ea3" | type == "Eb3" | type == "Eab3",
                                            "High", NA))))
  #classify infections by resistance profile & symp/asymp status
  data_long <- data_long %>%
    mutate(Profile = str_sub(type, end=-2),
           InfectState = ifelse(grepl("S", type, fixed = TRUE), "Susceptible",
                                ifelse(grepl("Y", type, fixed = TRUE), "Sympomatic Infected",
                                       ifelse(grepl("Z", type, fixed = TRUE), "Asymptomatic Infected",
                                              ifelse(grepl("E", type, fixed = TRUE), "Exposed", NA)))),
           ResistState = ifelse(Profile == "Y0" | Profile == "Z0", 'Neither',
                                ifelse(Profile == "Ya" | Profile == "Za", "A only",
                                       ifelse(Profile == "Yb" | Profile == "Zb", "B only",
                                              ifelse(Profile == "Yab" | Profile == "Zab", "A and B", NA)))))
  return(data_long)
}

# quick visualization of results of ODE solver by Resistance Profile, sexual activity group, Infection status
SI.visualize <- function(data) { 
  
  p1 <- data %>%
    group_by(time, ResistState) %>%
    summarise(N = sum(individuals)) %>%
    ggplot() + 
    geom_line(aes(x = time, y = log10(N), col = ResistState), size = 1.5) + 
    theme_light() + 
    ggtitle("GC Infections by Resistance Profile over Time") + ylim(c(0,6))
  
  p2 <- data %>%
    group_by(time, ResistState, RiskGroup) %>%
    summarise(N = sum(individuals)) %>%
    ggplot() + 
    geom_line(aes(x = time, y = log10(N), col = ResistState), size = 1.5) + 
    theme_light() + 
    ggtitle("GC Infections by Resistance Profile over Time, by Risk") + ylim(c(0,6)) +
    facet_wrap(~RiskGroup, nrow = 3)
  
  
  p3 <- data %>% 
    group_by(time, InfectState) %>% 
    summarise(N = sum(individuals)) %>%
    ggplot() + 
    geom_line(aes(x = time, y = log10(N), col = InfectState), size = 1.5) + 
    theme_light() + 
    scale_color_manual(values = c("deeppink3","midnightblue", "orangered", "gold")) +
    ggtitle("GC Infections by Clinical Profile over Time")
  
  p4 <- data %>%
    group_by(time, InfectState, RiskGroup) %>%
    summarise(N = sum(individuals)) %>%
    ggplot() +
    geom_line(aes(x = time, y = log10(N), col = InfectState), size = 1.5) + 
    theme_light() + 
    facet_wrap(~RiskGroup) + 
    scale_color_manual(values = c("deeppink3","midnightblue", "orangered", "gold")) +
    ggtitle("GC Infections by Clinical Profile over Time, by Risk")
  list(p1, p2, p3, p4)
}

# calculate summary stats from ODE output
# including time to loss (5% resistance) of Drug A, drug B, dual resistance
# as well as prevalence of gonorrhea at equilibrium, and @ time of loss
SI.outputs <- function(data, data_long) {
  end <- data[nrow(data),]
  #overall prevalence of gonorrhea at t=end
  prev_GC <- rowSums(end[setdiff(names(end),c("time", "S1", "S2", "S3", "E01", "E02", "E03", "Ea1", "Ea2", "Ea3",
                                                  "Eb1", "Eb2", "Eb3", "Eab1", "Eab2", "Eab3", "DoxyPEP", 
                                              "Inc", "prevB", "prevA", "prevAB", "prevGC", "ARx"))])/10^6
  
  #proportion of cases w/ some resistance at t=end
  prev_resistance <- rowSums(end[setdiff(names(end),c("time", "S1", "S2", "S3", "E01", "E02", "E03", "Ea1", "Ea2", "Ea3",
                                                      "Eb1", "Eb2", "Eb3", "Eab1", "Eab2", "Eab3", "Y01", "Y02", "Y03", 
                                                      "Z01", "Z02", "Z03", "DoxyPEP", "Inc", "prevB", "prevA", "prevAB", 
                                                      "prevGC", "ARx"))])/
    rowSums(end[setdiff(names(end),c("time", "S1", "S2", "S3", "E01", "E02", "E03", "Ea1", "Ea2", "Ea3",
                                                                       "Eb1", "Eb2", "Eb3", "Eab1", "Eab2", "Eab3", 
                                     "DoxyPEP", "Inc", "prevB", "prevA", "prevAB", "prevGC", "ARx"))])
  
  data_cond <- data_long %>%
    group_by(time, Profile) %>%
    summarise(N = sum(individuals)) %>%
    spread(., key = Profile, value = N) %>%
    mutate(prevA = round((Ya + Za + Yab + Zab)/(Y0 + Z0 + Ya + Za + Yb + Zb + Yab + Zab),3),
           prevB = round((Yb + Zb + Yab + Zab)/(Y0 + Z0 + Ya + Za + Yb + Zb + Yab + Zab),3),
           prevAB = round((Yab + Zab)/(Y0 + Z0 + Ya + Za + Yb + Zb + Yab + Zab),3))
  
  #Time to loss of one or both drugs (>5% prevalence of resistance)
  LossA_5 = min(data_cond$time[data_cond$prevA >= 0.05])/365
  LossB_5 = min(data_cond$time[data_cond$prevB >= 0.05])/365
  LossAB_5 = min(data_cond$time[data_cond$prevAB >= 0.05])/365
  
  #Time to loss of one or both drugs (>1% prevalence of resistance)
  LossA_1 = min(data_cond$time[data_cond$prevA >= 0.01])/365
  LossB_1 = min(data_cond$time[data_cond$prevB >= 0.01])/365
  LossAB_1 = min(data_cond$time[data_cond$prevAB >= 0.01])/365
  
  Loss_Each5 = max(LossA_5[is.finite(LossA_5)], LossB_5[is.finite(LossB_5)])*365
  time_loss = data[data$time == Loss_Each5,]
  prev_GC_loss <- rowSums(time_loss[setdiff(names(time_loss),c("time", "S1", "S2", "S3", "E01", "E02", "E03", "Ea1", "Ea2", "Ea3",
                                                   "Eb1", "Eb2", "Eb3", "Eab1", "Eab2", "Eab3", "DoxyPEP", 
                                                   "Inc", "prevB", "prevA", "prevAB", "prevGC", "ARx"))])/10^6
  
  list(PrevGC = prev_GC, PrevRes = prev_resistance, PrevGCLoss = prev_GC_loss,
       LossA_5 = LossA_5, LossB_5 = LossB_5, LossAB_5 = LossAB_5,
       LossA_1 = LossA_1, LossB_1 = LossB_1, LossAB_1 = LossAB_1)
}

# calculate GC prevalence based on current parameters, for calibration
pred_fun_er <- function(params){
  params2 <- list(pop=10^6, pop.p=c(0.3, 0.6, 0.1), rho = 1/(20*365), 
                  prA = 1, omega_a = 10^-8, fA = 0.98, 
                  pi_s = 0.90)
  parms <- c(params, params2)
  years = 2
  dt <- seq(0, 365*years, 1) 
  calibration_sim <- as.data.frame(ode(inits, dt, calibration.SI, parms = parms))   
  calibration_sim <- round(calibration_sim, 0)
  end <- calibration_sim[nrow(calibration_sim),]
  #Overall prevalence of gonorrhea at t=end years
  prev_GC <- sum(end[,8:22])/sum(end[,2:22])
  prev_GC
  return(prev_GC)
}

# estimate alpha and beta for our beta distribution of prevalence of GC, for calibration
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

logit<-function(x) {log(x/(1-x))}

ilogit <-function(x) {1/(1+exp(-x))}
