# Post processing of the results of the simulation of the model.
# Author YS
# Date: 1/21/2025
library(tidyverse)
library(emmeans)
library(broom)
library(splines)
library(mmrm)
library(truncnorm)
library(tmvtnorm)
library(lme4)
load(".")

names(sim_FTD_output)

# power evaluation
sim_FTD_output %>%
  group_by(dgm, model, misrate) %>%
  summarise(
    power = mean(pvalue < 0.05)
  )

# Getting the threshold for p-value to have < 5% under Null
load("./Outputs/tak594_FTD_SIM_02122025.RData")
sim_FTD_output %>% filter(misrate ==0, (dgm == "NULL" | dgm == "NULL_N"), model == "ncs-ranslp") %>%
  summarize(q = quantile(pvalue, 0.05))
# P-value cutoff should be 0.00446 for natural cubic spline


sim_FTD_output %>%
  group_by(b0, b1, b2, delta2, dgm, misrate, model) %>%
  summarise(
    power = mean(pvalue < 0.05),
                 estimate = mean(estimate)) %>%
  filter(dgm == "SP")


sim_FTD_output %>%
  group_by(b0, b1, b2, delta1 , delta2, delta3, sd1, sd2, sd3, sd4, sd5, cor, jitter_sd, missingPercentage,
           mtP1, mtP2, mtP3, mtP4, mtP5, dgm, model) %>%
  summarise(
    power = mean(pvalue < 0.05),
    estimate = mean(estimate)) %>% view()
  filter(dgm != "NULL_N", misrate ==.1)%>%view()
