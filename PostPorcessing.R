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
load("./Outputs/tak594_FTD_SIM_01312025.RData")

names(sim_FTD_output)

# power evaluation

sim_FTD_output %>%
  group_by(dgm, model, misrate) %>%
  summarise(
    power = mean(pvalue < 0.05)
  )

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
