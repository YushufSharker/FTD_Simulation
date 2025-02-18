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


# Getting the threshold for p-value to have < 5% under Null
load("./Outputs/tak594_FTD_SIM_02122025.RData")
sim_FTD_output %>% filter(misrate ==0, (dgm == "NULL" | dgm == "NULL_N"), model == "ncs-ranslp") %>%
  summarize(q = quantile(pvalue, 0.05))
# P-value cutoff should be 0.00446 for natural cubic spline
# datafiles
# load("./Outputs/tak594_FTD_SIM_02122025.RData")
#load("./Outputs/tak594_FTD_SIM_02132025.RData")
# write.csv(
# sim_FTD_output %>%
#   group_by(dgm, model, misrate) %>%
#   summarise(
#     power = mean(pvalue < 0.05)
#   ),
# "general12pvalu153.csv"
# )

# power evaluation
load("./Outputs/tak594_FTD_SIM_02132025.RData")
# for MMRM
write.csv(
sim_FTD_output %>%
  group_by(dgm, model, misrate) %>% filter(model == "MMRM")%>%
  summarise(
    power = mean(pvalue < 0.05)
  )
,"mmrm153.csv"
)

# for cubic spline
write.csv(
  sim_FTD_output %>%
  group_by(dgm, model, misrate) %>% filter(model == "ncs-ranslp")%>%
  summarise(
    power = mean(pvalue < 0.004)
  )
,"cs153.csv"
)


# getting the estimates or reduction
write.csv(
  sim_FTD_output %>%
    group_by(dgm, model, misrate) %>%
    summarise(
      reduction = mean(emmean), absreduction = mean(estimate)
    )
  ,"estimates153.csv"
)

# Comment: make a figure, histogram for the estimates
#load("./Outputs/tak594_FTD_SIM_02122025.RData")

# For presentation
load("./Outputs/tak594_FTD_SIM_02142025.RData") # random intercept + slope for each time point
sim_FTD_output %>%
  group_by(dgm, model, misrate) %>%
  summarise(
    power = mean(pvalue < 0.05)
  ) %>% pivot_wider(names_from = dgm, values_from = power) %>% group_by(misrate) %>%
  gt()

load("./Outputs/tak594_FTD_SIM_02132025.RData") # random intercept only
sim_FTD_output %>%
  group_by(dgm, model, misrate) %>%
  summarise(
    power = mean(pvalue < 0.05)
  ) %>% pivot_wider(names_from = dgm, values_from = power) %>% group_by(misrate) %>%
  gt()

sim_FTD_output %>%
  group_by(dgm, model, misrate) %>% filter(model == "ncs-ranslp")%>%
  summarise(
    power = mean(pvalue < 0.004)
  )%>% pivot_wider(names_from = dgm, values_from = power)%>%
  group_by(misrate) %>%
  gt()
