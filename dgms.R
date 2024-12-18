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
sd = c(1, 2, 2, 2, 2)
beta = c (3.4,11.3,3.7) # ncs time basis coefficients

delta1 = 0.3 # 30% reduction for times >0
delta2 = 0.3 # 30% slower progression

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
error <-as.vector(matrix(t(mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = cov))))


dat <- placeb_model(M = c(0,6,12,18,24), beta = c (3.4,12,5))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
 ggplot(data = dat, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()

#Null
dat0 <- dat %>% mutate(y = fixef0 + error) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
    ungroup() %>%
    mutate(group = as.factor(group), id = as.factor(id))

  # ggplot(data = dat0, aes(x = month, y = y, colour = as.factor(group), group = id)) +geom_line()+
  #   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()

# 30% reduction for times >0
dat1 <- dat %>% mutate(y = fixef0 + (Mcat6 + Mcat12 + Mcat18 + Mcat24) *
                         (-delta1) * fixef0 * group + error) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# ggplot(data = dat1, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line()+
#   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()


# 30% slower progression
dat2 <- dat %>% group_by(id) %>%
  mutate(y = spline(
    x = M * (1 + delta2 * group),
    y = c(3.4, 6.63, 9.02, 10.0, 10.1),
    method = "natural",
    xout = c(0, 6, 12, 18, 24)
  )$y) %>%
  ungroup() %>%
  mutate(y = y + error) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%  mutate(group = as.factor(group), id = as.factor(id))

ggplot(data = dat2, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line()+
geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()

# Linear drug effect and observe x% reduction from placebo at time 24.
# (reduction proportional to time ???)
dat3 <- dat %>% group_by(id) %>%
  mutate(y = fixef0*(1-spline(x = c(min(M), max(M)), y=c(.001, .3), method = "natural", xout = c(0, 6, 12, 18, 24))$y*group)) %>%
  ungroup() %>%
  mutate(y = y + error) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# ggplot(data = dat3, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line()+
# geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()

library(mmrm)
library(broom)
fit <- mmrm(
  formula = chg ~ group*Mcat+ ar1(Mcat | id),
  data = dat1,
  control = mmrm_control(method = "Kenward-Roger")
)

# tests
cov2cor((fit)$cov)
summary(fit)
lsm <-emmeans(fit, ~ group | Mcat, type = "response")
tibble(as.data.frame(emmeans(fit, ~ group | Mcat, type = "response")[9:10])) %>% mutate(simno = 1)
data.frame(confint(pairs(lsm, reverse = TRUE)))

# do experiment#copied
fit_data<- model.frame(fit)
LSmeans_object <- emmeans::emmeans(   fit,   data = fit_data,   specs = c("Mcat", "group"),   weights = "proportional" )
LSmeans_object

library(tern.mmrm)
get_mmrm_lsmeans(fit,  group * Mcat , conf_level, weights, averages = list())
tidy(fit) %>% select(term, estimate, `p.value`) %>%filter(term == "group1:Mcat24") %>% select(-term)
library(sandwich)
library(rigr)
test <- c(0,0,0,0,1,0,1,0,0,0)
lincom(fit, test)

#Ref resource for MMRM https://cran.r-project.org/web/packages/mmrm/vignettes/introduction.html
