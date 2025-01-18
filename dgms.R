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
# ggplot(data = dat, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
# geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()
# dat %>% group_by(M) %>% dplyr::summarize(mean(fixef0))

#Null
dat0 <- dat %>% mutate(y = fixef0 + error) %>%
  mutate(y = plyr::round_any(y, 0.5),
         dgm = "NULL",
         errm = "N")%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
    ungroup() %>%
    mutate(group = as.factor(group), id = as.factor(id))

 # f0<- ggplot(data = dat0, aes(x = month, y = y, colour = as.factor(group), group = id)) +geom_line(lwd = .25)+
 #    geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom")+
 #    scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
 #    xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("Null Scenario")

 #with independent truncated normal
 dat02 <- dat %>% mutate(y = fixef0 + error2) %>%
   mutate(y = plyr::round_any(y, 0.5),
          dgm = "NULL",
          errm = "TN")%>%
   group_by(id) %>%
   mutate(chg = y - y[1L]) %>% ungroup() %>%
   ungroup() %>%
   mutate(group = as.factor(group), id = as.factor(id))


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

  # 30% proportional reduction with truncated normal error
  dat12 <- dat %>% mutate(chg = chg*(1-delta1*group),
                          dgm = "PR",
                          errm = "TN")%>%
    group_by(id) %>% mutate(fixef0 = chg+fixef0[1L]) %>%
    ungroup() %>%
    mutate(y = fixef0 +error2 ) %>%
    mutate(y = plyr::round_any(y, 0.5))%>%
    group_by(id) %>%
    mutate(chg = y - y[1L]) %>% ungroup() %>%
    mutate(group = as.factor(group), id = as.factor(id))

# f1<-ggplot(data = dat1, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd=.25)+
#   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom")+
#   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
#   xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% Proportional reduction")


# 30% slower progression
  spline_interpolation <- function(df) {
    delta2 = .3
    df= df %>% mutate(x = M * (1 + delta2 * group))%>%
      mutate(y= spline(x = x, y = fixef0, method = "natural", xout = M)$y)
    return(df)
  }

  #  x<- dat %>% group_by(id)%>%  nest()
  # dat2 <-1: nrow(x) %>%
  #   map_dfr(~spline_interpolation(x$data[[.x]]))
  dat2i <- do.call(rbind, lapply(split(dat, f = dat$id), function(i) spline_interpolation(i)))
  dat2<-dat2i %>% select(-x) %>%

# dat2 <- dat %>% dplyr::group_by(id) %>%
#
#   mutate(y = spline(
#     x = M * (1 + delta2 * group),
#     y = fixef0 , # c(3.4, 6.63, 9.02, 10.0, 10.1),
#     method = "natural",
#     xout = c(0, 6, 12, 18, 24)
#   )$y) %>%
#
#   dplyr::ungroup() %>%
  mutate(y = y + error,
         dgm = "SP",
         errm = "N") %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# 30% slower progression with truncated normal error
# dat22 <- dat %>% dplyr::group_by(id) %>%
#   mutate(y = spline(
#     x = M * (1 + delta2 * group),
#     y = fixef0     ,# c(3.4, 6.63, 9.02, 10.0, 10.1),
#     method = "natural",
#     xout = c(0, 6, 12, 18, 24)
#   )$y) %>%
#   ungroup() %>%
  dat22 <- dat2i%>% select(-x) %>%
  mutate(y = y + error2,
         dgm = "SP",
         errm = "TN") %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# f2<-ggplot(data = dat22, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd = .25)+
# geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2)  +
#   theme(legend.position = "bottom")+
#   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
#   xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% Slower Progression")


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


dat32 <- dat %>%
  mutate(chg = chg - (spline(
    x = c(min(M), max(M)),
    y = c(0, delta3),
    method = "natural",
    xout = c(0, 6, 12, 18, 24)
  )$y * group)) %>%
  group_by(id) %>%
  mutate(fixef0 = chg + fixef0[1L]) %>%
  ungroup() %>%
  mutate(y = fixef0 + error2) %>%
  mutate(y = plyr::round_any(y, 0.5)) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id), ,
         dgm = "TP",
         errm = "TN")


# f3<-ggplot(data = dat3, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd=.25, alpha = .25)+
# geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) +
#   theme(legend.position = "bottom")+
#   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
#   xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% gradual reduction")



return(list(dat0, dat02, dat1, dat12, dat2, dat22, dat3, dat32))
}


#
# # Model Fits mmrm
# mmrm_fit <- function(dat = dat3)
# fit <- mmrm::mmrm(
#   formula = chg ~ group*Mcat+ ar1(Mcat | id),
#   data = dat32,
#   control = mmrm_control(method = "Kenward-Roger")
# )
#
#
#
# # tests
# # cov2cor((fit)$cov)
# # summary(fit)
# # lsm <-emmeans(fit, ~ group | Mcat, type = "response")
# Overall_prop_change <- tibble(as.data.frame(emmeans(fit, ~ group | Mcat, type = "response")[9:10])) %>% mutate(simno = 1) %>%
# with(., round(1-.[2,3]/.[1,3], 2))
# pvalue <-  round( summary(fit)[["coefficients"]][10, 5], 5 )
# estimate<- round( summary(fit)[["coefficients"]][10, 1], 2 )
#
# #data.frame(confint(pairs(lsm, reverse = TRUE)))
#
# # do experiment#copied
# fit_data<- model.frame(fit)
# LSmeans_object <- emmeans::emmeans(   fit,   data = fit_data,   specs = c("Mcat", "group"),   weights = "proportional" )
# LSmeans_object
#
# library(tern.mmrm)
# get_mmrm_lsmeans(fit,  group * Mcat , conf_level, weights, averages = list())
# tidy(fit) %>% select(term, estimate, `p.value`) %>%filter(term == "group1:Mcat24") %>% select(-term)
# library(sandwich)
# library(rigr)
# test <- c(0,0,0,0,1,0,1,0,0,0)
# lincom(fit, test)

#Ref resource for MMRM https://cran.r-project.org/web/packages/mmrm/vignettes/introduction.html
