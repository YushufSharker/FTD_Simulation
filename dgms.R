# Author: Yushuf Sharker
# Date: 2024-12-16
# These are the functions to supplement or take the lead in the simulation for FTD
# to evaluate the power and type-1 error of CS2 compared to MMRM
# https://studio-insight.rda.onetakeda.com/
# Functions for generating placebo progression

# Packages
library(tidyverse)
library(emmeans)
library(broom)
library(splines)
library(mmrm)
library(truncnorm)
library(plyr)
library(tmvtnorm)


# library(mice) # Multiple imputation
# library(miceadds)
# library(pan)
# library(longpower) # For power/sample size calculation for longitudinal data

# Function to generate placebo progression
#inputs
source("placebo_model.R")
b0 = 3.4
b1=11.3
b2 = 3.7
# ncs time basis coefficients
#sds fro 5 time points
sd1 = 2
sd2 = 2
sd3 = 2
sd4 = 3
sd5 = 3

delta1 = 0.3 # 30% reduction for times >0
delta2 = 0.3 # 30% slower progression
delta3 = 4 # x unit gradual reduction

n_pbo = 40
n_act = 80

dgm <- function(delta1 = 0.3, delta2 = 0.3, delta3 = 4, b0 = 3.4, b1=11.3, b2 = 3.7,
                 n_pbo = 40, nact=80,
                sd1 = 2, sd2 = 2, sd3 = 2, sd4 = 2, sd5 = 2  ){

M = c(0,6,12,18,24)
jitter_sd = 0.8 ### patients visit windows
ncs_df = 2
beta = c (b0, b1, b2)
sd = c(sd1, sd2, sd3, sd4, sd5)

corr = matrix(c(1,    0.65, 0.40, 0.25, 0.15,
                0.65, 1,    0.65, 0.40, 0.25,
                0.40, 0.65, 1,    0.65, 0.40,
                0.25, 0.40, 0.65, 1,    0.65,
                0.15, 0.25, 0.40, 0.65, 1   ),  nrow = 5, byrow = TRUE)

cov = diag(sd) %*% corr %*% diag(sd)
m = length(M)
n = n_pbo + n_act

error <- as.vector(matrix(t(
  mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = cov)
)))

error2 <- as.vector(matrix(t(
  rtmvnorm(
    n = n,
    mean = rep(0, m),
    sigma = cov,
    lower = rep(-.6, m),
    upper = rep(2, m)
  )
)))


# Placebo model
dat <- placeb_model(M = c(0,6,12,18,24), beta = c (2.5,14,5))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()

#Null
dat0 <- dat %>% mutate(y = fixef0 + error) %>%
  mutate(y = plyr::round_any(y, 0.5),
         dgm = "NULL",
         errm = "N")%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
    ungroup() %>%
    mutate(group = as.factor(group), id = as.factor(id))


 #with independent truncated normal
     # dat02 <- dat %>% mutate(y = fixef0 + error2) %>%
     #   mutate(y = plyr::round_any(y, 0.5),
     #          dgm = "NULL",
     #          errm = "TN")%>%
     #   group_by(id) %>%
     #   mutate(chg = y - y[1L]) %>% ungroup() %>%
     #   ungroup() %>%
     #   mutate(group = as.factor(group), id = as.factor(id))


# 30% proportional reduction
  dat1 <- dat %>% mutate(chg = chg*(1-delta1*group),
                         dgm = "PR",
                         errm = "N") %>%
    group_by(id) %>% mutate(fixef0 = chg+fixef0[1L]) %>%
    ungroup() %>%
    mutate(y = fixef0 +error ) %>%
    mutate(y = plyr::round_any(y, 0.5))%>%
    group_by(id) %>%
    mutate(chg = y - y[1L]) %>% ungroup() %>%
    mutate(group = as.factor(group), id = as.factor(id))

      # # 30% proportional reduction with truncated normal error
      # dat12 <- dat %>% mutate(chg = chg*(1-delta1*group),
      #                         dgm = "PR",
      #                         errm = "TN")%>%
      #   group_by(id) %>% mutate(fixef0 = chg+fixef0[1L]) %>%
      #   ungroup() %>%
      #   mutate(y = fixef0 +error2 ) %>%
      #   mutate(y = plyr::round_any(y, 0.5))%>%
      #   group_by(id) %>%
      #   mutate(chg = y - y[1L]) %>% ungroup() %>%
      #   mutate(group = as.factor(group), id = as.factor(id))



# 30% slower progression
  spline_interpolation <- function(df) {
    delta2 = .3
    df= df %>% mutate(x = M * (1 + delta2 * group))%>%
      mutate(y= spline(x = x, y = fixef0, method = "natural", xout = M)$y)
    return(df)
  }


  dat2i <- do.call(rbind, lapply(split(dat, f = dat$id), function(i) spline_interpolation(i)))
  dat2<-dat2i %>% select(-x) %>%
  mutate(y = y + error,
         dgm = "SP",
         errm = "N") %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# 30% slower progression with truncated normal error
      # dat22 <- dat2i%>% select(-x) %>%
      # mutate(y = y + error2,
      #        dgm = "SP",
      #        errm = "TN") %>%
      # mutate(y = plyr::round_any(y, 0.5))%>%
      # group_by(id) %>%
      # mutate(chg = y - y[1L]) %>%
      # ungroup() %>%
      # mutate(group = as.factor(group), id = as.factor(id))


# Linear drug effect and observe delta3 unit absolute reduction from placebo at time 24.
# decrease is linear over time
# (reduction proportional to time ???)
dat3 <- dat %>%
  mutate(chg = chg - (spline(
    x = c(min(M), max(M)),
    y = c(0, delta3),
    method = "natural",
    xout = c(0, 6, 12, 18, 24)
  )$y * group)) %>%
  group_by(id) %>%
  mutate(fixef0 = chg + fixef0[1L]) %>%
  ungroup() %>%
  mutate(y = fixef0 + error) %>%
  mutate(y = plyr::round_any(y, 0.5)) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id),
         dgm = "TP",
         errm = "N")


    # With truncated normal
    #   dat32 <- dat %>%
    #   mutate(chg = chg - (spline(
    #     x = c(min(M), max(M)),
    #     y = c(0, delta3),
    #     method = "natural",
    #     xout = c(0, 6, 12, 18, 24)
    #   )$y * group)) %>%
    #   group_by(id) %>%
    #   mutate(fixef0 = chg + fixef0[1L]) %>%
    #   ungroup() %>%
    #   mutate(y = fixef0 + error2) %>%
    #   mutate(y = plyr::round_any(y, 0.5)) %>%
    #   group_by(id) %>%
    #   mutate(chg = y - y[1L]) %>%
    #   ungroup() %>%
    #   mutate(group = as.factor(group), id = as.factor(id), ,
    #          dgm = "TP",
    #          errm = "TN")

return(list(dat0, dat1, dat2, dat3))
}


