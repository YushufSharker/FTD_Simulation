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
library(tmvtnorm)
library(lme4)

delta1 = 0.3
delta2 = 0.7
delta3 = 4
b0 = 3.4
b1 = 11.3
b2 = 3.7
nhm1 = 3.4
nhm2 = 6.5
nhm3 = 8.9
nhm4 = 9.9
nhm5 = 10.1
n_pbo = 40
n_act = 80
sd1 = 2
sd2 = 3
sd3 = 4
sd4 = 5
sd5 = 6
# r12= 0.65,
# r13= 0.40,
# r14= 0.25,
# r15= 0.15,
# r23= 0.65,
# r24= 0.40,
# r25= 0.25,
# r35= 0.65,
# r45= 0.65,
cor = 0.65
jitter_sd = .8
missingPercentage = .1
mtP1 = 0.05
mtP2 = 0.1
mtP3 = 0.2
mtP4 = 0.3
mtP5 = 0.4

# jitter_sd = 0.8 ### patients visit windows
dgm <- function(delta1 = 0.3,
                delta2 = 0.7,
                delta3 = 4,
                b0 = 3.4,
                b1 = 11.3,
                b2 = 3.7,
                nhm1 = 3.4,
                nhm2 = 6.5,
                nhm3 = 8.9,
                nhm4 = 9.9,
                nhm5 = 10.1,
                n_pbo = 40,
                n_act = 80,
                sd1 = 2,
                sd2 = 3,
                sd3 = 4,
                sd4 = 5,
                sd5 = 6,
                # r12= 0.65,
                # r13= 0.40,
                # r14= 0.25,
                # r15= 0.15,
                # r23= 0.65,
                # r24= 0.40,
                # r25= 0.25,
                # r35= 0.65,
                # r45= 0.65,
                cor = 0.65,
                jitter_sd = .8,
                missingPercentage = .1,
                mtP1 = 0.05,
                mtP2 = 0.1,
                mtP3 = 0.2,
                mtP4 = 0.3,
                mtP5 = 0.4) {

Nhmean = c(nhm1, nhm2, nhm3, nhm4, nhm5)
# paired_correlations = c(r12, r13, r14, r15, r12, r13, r14, r12, r13, r12)
M = c(0,6,12,18,24)
ncs_df = 2
beta = c (b0, b1, b2)
sd = c(sd1, sd2, sd3, sd4, sd5)

#corr = usCorrelation(correlations = paired_correlations)
corr = autocorr.mat(p = 5, rho = cor)
cov = diag(sd) %*% corr %*% diag(sd)
m = length(M)
n = n_pbo + n_act


# error2 <- as.vector(matrix(t(
#   rtmvnorm(
#     n = n,
#     mean = rep(0, m),
#     sigma = cov,
#     lower = rep(-.6, m),
#     upper = rep(2, m)
#   )
# )))


# Placebo model
dat <- placeb_model(M = M, beta = beta,
                    n_pbo = n_pbo, n_act=n_act,
                    ncs_df=ncs_df, jitter_sd =jitter_sd)
  # group_by(id) %>%
  # mutate(chg = fixef0 - fixef0[1L]) %>%
  # ungroup()
set.seed(NULL)
error <- as.vector(matrix(t(
  mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = cov, pre0.9_9994 = TRUE)
)))

#Null by NCS
dat0 <- dat %>% mutate(y = fixef0 + error) %>%
    mutate(y = plyr::round_any(y, 0.5),
           dgm = "NULL",
           errm = "N",
           misrate = 0)%>%
    group_by(id) %>%
    mutate(chg = y - y[1L]) %>% ungroup() %>%
    mutate(group = as.factor(group), id = as.factor(id))

#Null by NCS with x% missing
  dat0m<- introduce_missing(df = dat0, outVariable = "chg",
                            missing_percentage = missingPercentage,
                            prob = c(mtP1, mtP2, mtP3, mtP4, mtP5)) %>%
    mutate(misrate = missingPercentage) %>%
    mutate(y = case_when(is.na(chg) ~ NA_real_, TRUE ~ y ))


# Null by Natural history
dat0N <- dat %>%
  mutate(error = error) %>% group_by(id) %>%
  mutate(y =  Nhmean+error) %>%
  ungroup() %>%
  mutate(y = plyr::round_any(y, 0.5),
         dgm = "NULL_N",
         errm = "N",
         misrate = 0)%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>% select(-error) %>%
  mutate(group = as.factor(group), id = as.factor(id))

#Null by Natural history with x% missing
dat0Nm<- introduce_missing(df = dat0N, outVariable = "chg",
                          missing_percentage = missingPercentage,
                          prob = c(mtP1, mtP2, mtP3, mtP4, mtP5)) %>%
  mutate(misrate = missingPercentage) %>%
  mutate(y = case_when(is.na(chg) ~ NA_real_, TRUE ~ y
  ))


# 30% proportional reduction

  dat1 <- dat %>% mutate(y = fixef0*(1-delta1*group),
                         dgm = "PR",
                         errm = "N",
                         misrate = 0)%>%
    group_by(id) %>%
    mutate(y = if_else(row_number() == 1, fixef0[1L], y))%>%
    ungroup()%>%
    mutate(y = y + error)%>%
    group_by(id) %>%
    mutate(chg = y - y[1L]) %>% ungroup() %>%
    mutate(group = as.factor(group), id = as.factor(id))


#30% proportional reduction with x% missing
  dat1m<- introduce_missing(df = dat1, outVariable = "chg",
                             missing_percentage = missingPercentage,
                             prob = c(mtP1, mtP2, mtP3, mtP4, mtP5)) %>%
    mutate(misrate = missingPercentage)%>%
    mutate(y = case_when(is.na(chg) ~ NA_real_, TRUE ~ y ))


  # 30% slower progression
  spline_interpolation <- function(df, delta2 = delta2) {
    #delta2 = delta2
    df= df %>% mutate(x = M * (1 + delta2 * group))%>%
      mutate(y= spline(x = x, y = fixef0, method = "natural", xout = M)$y)
    return(df)
  }

  dat2i <- do.call(rbind, lapply(split(dat, f = dat$id), function(i) spline_interpolation(i,  delta2 = delta2)))
  dat2<-dat2i %>% select(-x) %>%
  mutate(y = y + error,
         dgm = "SP",
         errm = "N",
         misrate = 0) %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# 30% slower progression with x% missing
  dat2m<- introduce_missing(df = dat2, outVariable = "chg",
                              missing_percentage = missingPercentage,
                              prob = c(mtP1, mtP2, mtP3, mtP4, mtP5)) %>%
      mutate(misrate = missingPercentage) %>%
    mutate(y = case_when(is.na(chg) ~ NA_real_, TRUE ~ y))


# Linear drug effect and observe delta3 unit absolute reduction from placebo at time 24.
# decrease is linear over time
# (reduction proportional to time ???)

dat3 <- dat %>%
  mutate(y = fixef0 - spline(
    x = c(min(M), max(M)),
    y = c(.001, delta3),
    method = "natural",
    xout = c(0, 6, 12, 18, 24)
  )$y * group) %>%
  mutate(y = y + error) %>%
  mutate(y = plyr::round_any(y, 0.5)) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(
    group = as.factor(group),
    id = as.factor(id),
    dgm = "TP",
    errm = "N",
    misrate = 0
  )


# X% missing
dat3m<- introduce_missing(df = dat3, outVariable = "chg",
                          missing_percentage = missingPercentage,
                          prob = c(mtP1, mtP2, mtP3, mtP4, mtP5)) %>%
  mutate(misrate = missingPercentage) %>%
  mutate(y = case_when(is.na(chg) ~ NA_real_, TRUE ~ y ))

return(list(dat0, dat0m, dat0N, dat0Nm, dat1, dat1m, dat2, dat2m, dat3, dat3m))
}


