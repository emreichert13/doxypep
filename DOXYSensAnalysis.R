
##Run Sens Analysis##
# E Reichert, 03.2023

library(knitr)
library(tidyr)
library(ggplot2)
library(readr)
library(dplyr)


##I. PROPERTIES OF DRUG B
#Set parameter values for fB and omegaB that we will explore
#Must first comment out values of fB and omega_b in "DOXYtransmission.Rmd" script
fB_range <- c(seq(0.8,1,0.05)) #ranges from 0-20% fitness cost
omega_b_range <- c(0, 10^-8, 10^-4) #ranges for Pr of resistance emergence upon treatment

#run ODE models over all combinations of these params
SensAnalyze <- function (fB_range, omega_b_range) {
  df <- data.frame(matrix(ncol = 11, nrow = 0))
  #provide column names
  colnames(df) <- c("time", "DoxyPEP", "Inc", "CumInc", "IR", "prevA", "prevB", "prevAB", "prevGC", "fB", "omega_b")
  for(i in omega_b_range)
  {
    for(j in fB_range)
    {
      fB <- j
      omega_b <- i
      knit("DOXYtransmission.Rmd")
      
      output <- doxy_sim_all %>% 
        select(time, DoxyPEP, Inc, CumInc, IR, prevA, prevB, prevAB, prevGC) %>%
        mutate(fB = j, omega_b = i)
      df <- rbind(df, output)
    }
  }
  return(df)
}

SensResults <- SensAnalyze(fB_range, omega_b_range)
#write.csv(SensResults, "~/Documents/2021 Grad lab research/DOXYPEP/SensResults.csv")
#SensResults <- read.csv("~/Documents/2021 Grad lab research/DOXYPEP/SensResults.csv")

#calculate model outcomes
table_res <- SensResults %>% group_by(DoxyPEP, fB, omega_b) %>%
  summarise(MinA = round(min(time[prevA >= 0.05])/365,3),
            MinAB = round(min(time[prevAB >= 0.05])/365,3),
            MinB = round(min(time[prevB >= 0.875])/365,3),
            Prev20 = round(prevGC[time == 7300],3),
            A_20 = round(prevA[time == 7300],3),
            B_20 = round(prevB[time == 7300],3),
            AB_20 = round(prevAB[time == 7300],3),
            Inc20 = CumInc[time==7300],
            Inc500 = round(min(time[CumInc >= 500000])/365,3))

resB <- SensResults %>%
  group_by(DoxyPEP, fB, omega_b) %>%
  summarise(MinT = min(time[prevB >= 0.87], na.rm = T),
            CumInc = CumInc[time == MinT],
            Inc = Inc[time == MinT],
            prevGC = prevGC[time == MinT])


SensResults$omega_b <- as.character(SensResults$omega_b)
SensResults$omega_b[SensResults$omega_b == "1e-04"] <- "1e-4"

resB$omega_b <- as.character(resB$omega_b)
resB$omega_b[resB$omega_b == "1e-04"] <- "1e-4"

#Visualize prevalence of GC over time, by properties of Drug B (doxycycline)
#pdf("~/Documents/2021 Grad lab research/DOXYPEP/DOXYPEP_SensAnalysisPrev.pdf", height = 6, width = 8, encoding = "Greek.enc")
ggplot() +
  geom_line(data = SensResults, aes(x = time/365, y = prevGC*100, col = factor(DoxyPEP)), size = 1) +
  geom_point(data = resB, aes(x = MinT/365, y = prevGC*100, col = factor(DoxyPEP)), size = 3, shape = 19) +
  theme_light() + xlab("Years") + ylab("Gonococcal Infection Prevalence (%)") + labs(col = "DoxyPEP\nUptake") + 
  theme(text = element_text(size=14)) +
  scale_color_manual(values = c("#172869", "#1BB6AF", "#A6E000", "#FC6882", "#C70E7B")) +
  #scale_y_continuous(breaks = c(0, 1000000, 2000000, 3000000), labels = c("0", "1 M", "2 M", "3 M")) +
  facet_grid(paste0("\u03C9   = ", omega_b) ~ paste0("fB = " , fB))
#scale_color_manual(values = c("#172869", "#1BB6AF", "#A6E000", "#FF5300", "#FC6882"))
#dev.off()

baseline <- SensResults %>%
  filter(DoxyPEP == "0%") %>%
  select(time, fB, omega_b, 
         IR_compare = IR,
         prev_compare = prevGC,
         CumInc_compare = CumInc)

cumcases <- left_join(SensResults, baseline) %>%
  mutate(PR = prevGC/prev_compare,
         IRR = IR/IR_compare,
         CasesAverted = CumInc_compare - CumInc)

resB <- cumcases %>%
  group_by(DoxyPEP, fB, omega_b) %>%
  summarise(MinT = min(time[prevB >= 0.87]),
            PR = PR[time == MinT],
            IRR = IRR[time == MinT],
            CasesAverted = CasesAverted[time == MinT]) 

resB_20 <- cumcases %>%
  group_by(DoxyPEP, fB, omega_b) %>%
  summarise(CI20 = CumInc[time == 7300])

#Visualize PR over time relative to baseline (0% DoxyPEP uptake), by properties of Drug B (Doxycycline)
#pdf("~/Documents/2021 Grad lab research/DOXYPEP/DOXYPEP_SensAnalysisPR.pdf", height = 6, width = 8, encoding = "Greek.enc")
ggplot() +
  geom_line(data = cumcases, aes(x = time/365, y = PR, col = factor(DoxyPEP)), size = 1) +
  geom_point(data = resB, aes(x = MinT/365, y = PR, col = factor(DoxyPEP)), size = 3, shape = 19) +
  theme_light() + xlab("Years") + ylab("Prevalence Ratio") + labs(col = "DoxyPEP\nUptake") + 
  theme(text = element_text(size=14)) +
  scale_color_manual(values = c("#172869", "#1BB6AF", "#A6E000", "#FC6882", "#C70E7B")) +
  #scale_y_continuous(breaks = c(0, 1000000, 2000000, 3000000), labels = c("0", "1 M", "2 M", "3 M")) +
  facet_grid(paste0("\u03C9   =" , omega_b) ~ paste0("fB = " , fB))
#dev.off()

LossA <- SensResults %>% 
  group_by(DoxyPEP, fB, omega_b) %>%
  summarise(LossA = min(time[prevA >= 0.05])/365)

### II. EFFECTIVENESS OF DOXYPEP
# Must first comment out baseline value for kappa in "DOXYtransmission.Rmd" script
kappa_range <- c(seq(0,0.8,0.2))

#run ODE models over all combinations of these params
SensAnalyze <- function (kappa) {
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  #provide column names
  colnames(df) <- c("time", "DoxyPEP", "Inc", "CumInc", "IR", "prevA", "prevB", "prevAB", "prevGC", "kappa")
  for(k in kappa)
  {
      kappa <- k
      knit("DOXYtransmission.Rmd")
      
      output <- doxy_sim_all %>% 
        select(time, DoxyPEP, Inc, CumInc,IR, prevA, prevB, prevAB, prevGC) %>%
        mutate(kappa = k)
      df <- rbind(df, output)
  }
  return(df)
}

SensResults2 <- SensAnalyze(kappa_range)
#write.csv(SensResults2, "~/Documents/2021 Grad lab research/DOXYPEP/SensResults2.csv")
#SensResults2 <- read.csv("~/Documents/2021 Grad lab research/DOXYPEP/SensResults2.csv")

#calculate model outcomes
table_res2 <- SensResults2 %>% group_by(DoxyPEP, kappa) %>%
  summarise(MinA = round(min(time[prevA >= 0.05])/365,3),
            MinAB = round(min(time[prevAB >= 0.05])/365,3),
            MinB = round(min(time[prevB >= 0.875])/365,3),
            Prev20 = round(prevGC[time == 7300],3),
            A_20 = round(prevA[time == 7300],3),
            B_20 = round(prevB[time == 7300],3),
            AB_20 = round(prevAB[time == 7300],3),
            Inc20 = CumInc[time==7300],
            Inc500 = round(min(time[CumInc >= 500000])/365,3))

resB2 <- SensResults2 %>%
  group_by(DoxyPEP, kappa) %>%
  summarise(MinT = min(time[prevB >= 0.87]),
            CumInc = CumInc[time == MinT],
            Inc = Inc[time == MinT],
            prevGC = prevGC[time == MinT])

#Visualize prevalence over time, by DoxyPEP effectiveness per exposure
g1 <- ggplot() +
  geom_line(data = SensResults2, aes(x = time/365, y = prevGC*100, col = factor(DoxyPEP)), size = 1) +
  geom_point(data = resB2, aes(x = MinT/365, y = prevGC*100, col = factor(DoxyPEP)), size = 3, shape = 19) +
  theme_classic() + xlab("Years") + ylab("Gonococcal Infection\nPrevalence (%)") + 
  labs(col = "DoxyPEP Uptake") + theme(text = element_text(size=12), legend.text = element_text(size=12), axis.text = element_text(size=12)) +
  scale_color_manual(values = c("#172869", "#1BB6AF", "#A6E000", "#FC6882", "#C70E7B")) +
  #scale_y_continuous(breaks = c(0, 1000000, 2000000, 3000000), labels = c("0", "1 M", "2 M", "3 M")) +
  facet_grid(~paste0("\u03BA  = " , kappa)) + ggtitle("A.")

baseline2 <- SensResults2 %>%
  filter(DoxyPEP == "0%") %>%
  select(time, kappa, 
         IR_compare = IR,
         prev_compare = prevGC,
         CumInc_compare = CumInc)

cumcases2 <- left_join(SensResults2, baseline2) %>%
  mutate(PR = prevGC/prev_compare,
         IRR = IR/IR_compare,
         CasesAverted = CumInc_compare - CumInc)

resB2 <- cumcases2 %>%
  group_by(DoxyPEP, kappa) %>%
  summarise(MinT = min(time[prevB >= 0.87]),
            PR = PR[time == MinT],
            IRR = IRR[time == MinT],
            CasesAverted = CasesAverted[time == MinT])

#Visualize PR over time relative to baseline (0% DoxyPEP uptake), by DoxyPEP effectiveness per exposure
g2 <- ggplot() +
  geom_line(data = cumcases2, aes(x = time/365, y = PR, col = factor(DoxyPEP)), size = 1) +
  geom_point(data = resB2, aes(x = MinT/365, y = PR, col = factor(DoxyPEP)), size = 3, shape = 19) +
  theme_classic() + xlab("Years") + ylab("Prevalence Ratio\n   ") + labs(col = "DoxyPEP Uptake") + 
  theme(text = element_text(size=12), legend.text = element_text(size=12), axis.text = element_text(size=12)) +
  scale_color_manual(values = c("#172869", "#1BB6AF", "#A6E000", "#FC6882", "#C70E7B")) +
  #scale_y_continuous(breaks = c(0, 250000, 500000, 750000), labels = c("0", "0.25 M", "0.5 M", "0.75 M")) +
  facet_grid(~paste0("\u03BA  = " , kappa)) + ggtitle("B.")

library(gridExtra)
#pdf("~/Documents/2021 Grad lab research/DOXYPEP/DOXYPEP_SensAnalysisKAPPA.pdf", height = 6, width = 8, encoding = "Greek.enc")
grid.arrange(g1, g2, nrow = 2)
#dev.off()

LossA <- SensResults2 %>% 
  group_by(DoxyPEP, kappa) %>%
  summarise(LossA = min(time[prevA >= 0.05])/365)

#### III. DOXY RESISTANCE AT MODEL START
# Must first comment out baseline value for resB in "DOXYtransmission.Rmd" script
resB_range <- c(seq(0.2,0.8,0.2))

#run ODE models over all combinations of these params
SensAnalyze <- function (resB) {
  df <- data.frame(matrix(ncol = 11, nrow = 0))
  #provide column names
  colnames(df) <- c("time", "DoxyPEP", "Inc", "CumInc", "IR", 
                    "prev0", "prevA", "prevB", "prevAB", "prevGC", "resB")
  for(i in resB)
  {
    resB <- i
    knit("DOXYtransmission.Rmd")
    
    output <- doxy_sim_all %>% 
      select(time, DoxyPEP, Inc, CumInc,IR, prev0, prevA, prevB, prevAB, prevGC) %>%
      mutate(resB = i)
    df <- rbind(df, output)
  }
  return(df)
}

SensResults3 <- SensAnalyze(resB_range)
#write.csv(SensResults3, "~/Documents/2021 Grad lab research/DOXYPEP/SensResults3.csv")
#SensResults3 <- read.csv("~/Documents/2021 Grad lab research/DOXYPEP/SensResults3.csv")

#calculate model outcomes
table_res3 <- SensResults3 %>% group_by(DoxyPEP, resB) %>%
  summarise(MinA = round(min(time[prevA >= 0.05])/365,3),
            MinAB = round(min(time[prevAB >= 0.05])/365,3),
            MinB = round(min(time[prevB >= 0.875])/365,3),
            Prev20 = round(prevGC[time == 7300],3),
            A_20 = round(prevA[time == 7300],3),
            B_20 = round(prevB[time == 7300],3),
            AB_20 = round(prevAB[time == 7300],3),
            Inc20 = CumInc[time==7300],
            Inc500 = round(min(time[CumInc >= 500000])/365,3))

SensResults3 <- SensResults3 %>%
  mutate(resB_cat = ifelse(resB == 0.20, "20%",
                           ifelse(resB == 0.40, "40%", 
                                  ifelse(resB == 0.60, "60%", "80%"))))

resB3 <- SensResults3 %>%
  group_by(DoxyPEP, resB, resB_cat) %>%
  summarise(MinT = min(time[prevB >= 0.87]),
            CumInc = CumInc[time == MinT],
            Inc = Inc[time == MinT],
            prevGC = prevGC[time == MinT])

#Visualize prevalence over time, by doxy resistance at model start
g3 <- ggplot() +
  geom_line(data = SensResults3, aes(x = time/365, y = prevGC*100, col = factor(DoxyPEP)), size = 1) +
  geom_point(data = resB3, aes(x = MinT/365, y = prevGC*100, col = factor(DoxyPEP)), size = 3, shape = 19) +
  theme_classic() + xlab("Years") + ylab("Gonococcal Infection\nPrevalence (%)") + 
  labs(col = "DoxyPEP Uptake", title = "A.", subtitle = "Initial Prevalence of Doxycycline Resistance") + 
  theme(text = element_text(size=12), legend.text = element_text(size=12), axis.text = element_text(size=12), 
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#172869", "#1BB6AF", "#A6E000", "#FC6882", "#C70E7B")) +
  #scale_y_continuous(breaks = c(0, 1000000, 2000000, 3000000), labels = c("0", "1 M", "2 M", "3 M")) +
  facet_grid(~paste0(resB_cat))

baseline3 <- SensResults3 %>%
  filter(DoxyPEP == "0%") %>%
  select(time, resB, resB_cat,
         IR_compare = IR,
         prev_compare = prevGC,
         CumInc_compare = CumInc)

cumcases3 <- left_join(SensResults3, baseline3) %>%
  mutate(PR = prevGC/prev_compare,
         IRR = IR/IR_compare,
         CasesAverted = CumInc_compare - CumInc)

resB3 <- cumcases3 %>%
  group_by(DoxyPEP, resB, resB_cat) %>%
  summarise(MinT = min(time[prevB >= 0.87]),
            PR = PR[time == MinT],
            IRR = IRR[time == MinT],
            CasesAverted = CasesAverted[time == MinT])

#Visualize PR over time relative to baseline (0% DoxyPEP uptake), by doxy resistance at model start
g4 <- ggplot() +
  geom_line(data = cumcases3, aes(x = time/365, y = PR, col = factor(DoxyPEP)), size = 1) +
  geom_point(data = resB3, aes(x = MinT/365, y = PR, col = factor(DoxyPEP)), size = 3, shape = 19) +
  theme_classic() + xlab("Years") + ylab("Prevalence Ratio\n   ") + 
  labs(col = "DoxyPEP Uptake", title = "B.", subtitle = "Initial Prevalence of Doxycycline Resistance") + 
  theme(text = element_text(size=12), legend.text = element_text(size=12), axis.text = element_text(size=12),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#172869", "#1BB6AF", "#A6E000", "#FC6882", "#C70E7B")) +
  #scale_y_continuous(breaks = c(0, 250000, 500000, 750000), labels = c("0", "0.25 M", "0.5 M", "0.75 M")) +
  facet_grid(~paste0(resB_cat)) 

library(gridExtra)
pdf("~/Documents/2021 Grad lab research/DOXYPEP/DOXYPEP_SensAnalysisRESB.pdf", height = 7, width = 9)
grid.arrange(g3, g4, nrow = 2)
dev.off()

# ##add resistance profiles over time plot
# infectiondat <- SensResults3 %>%
#   select(resB, DoxyPEP, time, Neither = prev0, prevA, prevB, Dual = prevAB) %>%
#   mutate(Ceftriaxone = prevA-Dual,
#          Doxy = prevB-Dual,
#          totalprev = Neither + Ceftriaxone + Doxy + Dual) %>%
#   select(-prevA, -prevB, -totalprev) %>%
#   gather(., ResistState, percent, 4:7)
# 
# ggplot(infectiondat,aes(x=time/365, y=percent*100, fill = factor(ResistState, levels = c("Neither", "Doxy", "Ceftriaxone", "Dual")))) +
#   geom_area() + scale_fill_manual(values = c("turquoise3", "#E9A17C", "mediumpurple", "deeppink2")) +
#   theme_classic() + facet_grid(resB~DoxyPEP) +
#   xlab("Years") + ylab("% of Gonococcal Infections") + labs(fill = "Resistance Profile") +
#   theme(legend.position = c(0.87,0.2)) +
#   geom_hline(yintercept = 5, col = "white", linetype = "dashed") + geom_hline(yintercept = 87, col = "white", linetype = "dotted") +
#   geom_vline(xintercept = 10, col = 'white')

