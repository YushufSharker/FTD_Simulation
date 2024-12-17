# Author: Yushuf Sharker
# Date: 2024-12-16
# These are the functions to supplement or take the lead in the simulation for FTD
# to evaluate the power and type-1 error of CS2 compared to MMRM

# Functions for generating placebo progression

# Packages
library(tidyverse)
library(lme4)
library(emmeans)
library(splines)
library(mice) # Multiple imputation
library(miceadds)
library(pan)
library(longpower) # For power/sample size calculation for longitudinal data

# Function to generate placebo progression
#inputs

M = c(0,6,12,18,24)
n_pbo = 40
n_act = 80
jitter_sd = 0.8 ### patients visit windows
ncs_df = 2

corr = matrix(c(1,    0.65, 0.40, 0.25, 0.15,
                0.65, 1,    0.65, 0.40, 0.25,
                0.40, 0.65, 1,    0.65, 0.40,
                0.25, 0.40, 0.65, 1,    0.65,
                0.15, 0.25, 0.40, 0.65, 1   ),  nrow = 5, byrow = TRUE)
sd = c(2, 3, 4, 5, 6)
beta = c (3.4,11.3,3.7) # ncs time basis coefficients


placeb_model <- function(M, n_pbo = 40, n_act = 80, ncs_df = 2, beta){
  n = n_pbo + n_act
  m = length(M)
  visits <- tibble(visNo = 1:m, M)
   subinfo <- tibble(id = 1:n)
   group <- tibble(cbind(id = 1:n, tibble(group = c(rep(0, n_pbo), rep(1, n_act)))))
  visinfo <- expand_grid(id = 1:n,visNo = 1:m) %>%  group_by(visNo)
  d00 <- subinfo %>% left_join(visinfo, by="id") %>% left_join(visits, by="visNo") %>%
    left_join(group, by="id")
  ns_basis <- ns(d00$M,df=ncs_df)

dd00 <- d00 %>% mutate(Mcat=as.factor(M), month = M )%>%
  mutate(baseline = ifelse(M > 0, 1, 0))
# ns_fun <- lapply(1:2, function(x){function(t){as.numeric(predict(ns(datinput$month,2), t)[,x])}})
dd0 <- model.matrix(~ id + visNo + Mcat+ ns(dd00$month, 2), dd00 ) %>% as_tibble() %>% select(-"(Intercept)") %>%
  rename(ns1 = `ns(dd00$month, 2)1`, ns2 = `ns(dd00$month, 2)2`) %>%
  mutate(fixef0 = beta[1] + ns1*beta[2] + ns2*beta[3]) %>%
  left_join(dd00, by = c('id', 'visNo'))
return(dd0 = dd0)
}

cov = diag(sd) %*% corr %*% diag(sd)
ns_fun_true <- lapply(1:2, function(x){function(t){as.numeric(predict(ns(d00$M,2), t)[,x])}})

#mutate(months_jitter = case_when(visNo == 1 ~ 0, TRUE ~ rnorm(n=n*m, sd=jitter_sd)))%>%
#dd0 = placebo model data without error

ggplot(data = dd0, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) + geom_line()


# Functions to generate treatment progression based on placebo progression.

delta1=0.54
delta2=0.4
delta3=0.4
monthdelay=5

fixef24 = beta[1] + ns_fun_true[[1]](24)*beta[2] + ns_fun_true[[2]](24)*beta[3]
y1true = - delta1 * (ns_fun_true[[1]](24)*beta[2] + ns_fun_true[[2]](24)*beta[3])
y2true = - delta2 * (ns_fun_true[[1]](24-monthdelay)*beta[2] + ns_fun_true[[2]](24-monthdelay)*beta[3])
y3true = - delta3 * (ns_fun_true[[1]](24)*beta[2] + ns_fun_true[[2]](24)*beta[3])
resids_w <- mvtnorm::rmvnorm(n,sigma=cov)
colnames(resids_w) <- 1:ncol(resids_w)
resids <- resids_w %>% as_tibble() %>% mutate(id = 1:n, active = sample(0:1, size = n, replace=TRUE, prob=c(1/3,2/3))) %>%
  pivot_longer('1':'5', names_to = 'visNo', values_to = 'residual') %>%
  mutate(visNo = as.numeric(visNo),  Active = as.factor(active))

#Sim Data matrix
dd1 <- dd0 %>% left_join(resids, by=c('id', 'visNo')) %>% mutate(Active1 = if_else(M==0, "0", Active))%>% # mutate(Act_vis = with(dd0, interaction(Active1, M))) %>%
  mutate(# type I error
    y0 = fixef0,# + residual,
    # stable
    # 20 pct reduction
    y1 = y0 - active * delta1 * (ns_fun[[1]](month)*beta[2] + ns_fun[[2]](month)*beta[3]),
    # 5 month delay
    y2 = case_when(month < monthdelay ~ y0, month >= monthdelay ~ y0 - active * delta2*(ns_fun[[1]](month-monthdelay)*beta[2] + ns_fun[[2]](month-monthdelay)*beta[3])),
    y3 = y0 - active * delta3 * (ns_fun[[1]](month)*beta[2] + ns_fun[[2]](month)*beta[3]))
