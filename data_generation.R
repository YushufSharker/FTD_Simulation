# Author: Yushuf Sharker
# Data generation for treatment and placebo
# This model uses the PMRM code to generate data

sd = c(2, 3, 4, 5, 6)
beta = c (5,18,8)
delta1 = 0.3 # 30% reduction for times >0
delta2 = 0.3 # 30% slower progression
delta3 = 0.3 # 30% gradual reduction
    nhm1 = 3.4
    nhm2 = 6.5
    nhm3 = 8.9
    nhm4 = 9.9
    nhm5 = 10.1
    n_pbo = 40
    n_act = 80

datgen <- function(sd=c(1, 2, 2, 2, 2), delta1 = .3, delta2 = .3, beta = c (5,13,5)){
    Nhmean = c(nhm1, nhm2, nhm3, nhm4, nhm5)
    n = n_pbo + n_act
corr = matrix(c(1,    0.65, 0.40, 0.25, 0.15,
                0.65, 1,    0.65, 0.40, 0.25,
                0.40, 0.65, 1,    0.65, 0.40,
                0.25, 0.40, 0.65, 1,    0.65,
                0.15, 0.25, 0.40, 0.65, 1   ),  nrow = 5, byrow = TRUE)
cov = diag(sd) %*% corr %*% diag(sd)
error <-as.vector(matrix(t(mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = cov))))

beta = c (5,14,9)
dat <- placeb_model(M = c(0,6,12,18,24), beta = beta, jitter_sd = 0)%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
ggplot(data = dat, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2)

#Null nby NCS
dat0 <- dat %>% mutate(y = fixef0 + error) %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# Null By natural history
dat0N <- dat %>%
  mutate(error = error) %>% group_by(id) %>%
  mutate(y =  Nhmean+error) %>%
  ungroup() %>%
  mutate(y = plyr::round_any(y, 0.5),
         dgm = "NULL",
         errm = "N")

# f0 <- ggplot(data = dat0, aes(
#   x = month,
#   y = y,
#   colour = as.factor(group),
#   group = id
# )) + geom_line(lwd = .25) +
#   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom") +
#   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment")) +
#   xlab("Month") + ylab("FTLD CDR SB score") + ggtitle("Null Scenario")

# 30% reduction for all times >0
dat1 <- dat %>% mutate(chg = chg*(1-delta1*group))%>%
  group_by(id) %>% mutate(fixef0 = chg+fixef0[1L]) %>%
  ungroup() %>%
  mutate(y = fixef0 +error ) %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

dat1 <- dat %>% mutate(y = fixef0*(1-delta1*group))%>%
  group_by(id) %>%
  mutate(y = if_else(row_number() == 1, fixef0[1L], y))%>%
  ungroup()%>%
  mutate(y = y + error)%>%
  #mutate(chg = chg*(1-delta1*group))%>% # this code line was returning same data repeatedly
  # group_by(id) %>% mutate(fixef0 = chg+fixef0[1L]) %>%
  # ungroup() %>%
  #mutate(y = fixef0 +error ) %>%
  #mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

f1 <- ggplot(data = dat1, aes(
  x = month,
  y = chg,
  colour = as.factor(group),
  group = id
)) + geom_line(lwd = .25) +
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom") +
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment")) +
  xlab("Month") + ylab("FTLD CDR SB score") + ggtitle("30% Proportional reduction")

# 30% slower progression
dat2 <- dat %>% dplyr::group_by(id) %>%
  mutate(y = spline(
    x = M * (1 + delta2 * group),
    y = fixef0,
    # c(3.4, 6.63, 9.02, 10.0, 10.1),
    method = "natural",
    xout = c(0, 6, 12, 18, 24)
  )$y) %>%
  ungroup() %>%
  mutate(y = y + error) %>%
  mutate(y = plyr::round_any(y, 0.5)) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# f2 <- ggplot(data = dat2, aes(
#   x = month,
#   y = y,
#   colour = as.factor(group),
#   group = id
# )) + geom_line(lwd = .25) +
#   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2)  +
#   theme(legend.position = "bottom") +
#   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment")) +
#   xlab("Month") + ylab("FTLD CDR SB score") + ggtitle("30% Slower Progression")


# Linear drug effect and observe x% absolute reduction from placebo at time 24.
# (reduction proportional to time ???)

dat3 <- dat %>%
  # mutate(chg = chg * (1 - spline(
  #   x = c(min(M), max(M)),
  #   y = c(.001, delta3),
  #   method = "natural",
  #   xout = c(0, 6, 12, 18, 24)
  # )$y * group)) %>%
  #group_by(id) %>%
  #mutate(fixef0 = chg + fixef0[1L]) %>%
  mutate(y = fixef0-spline(x = c(min(M), max(M)), y=c(.001, delta3),
        method = "natural", xout = c(0, 6, 12, 18, 24))$y*group) %>%
  #ungroup() %>%
  mutate(y = y + error) %>%
  mutate(y = plyr::round_any(y, 0.5)) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id),
         dgm = "TP",
         errm = "N",
         misrate = 0)

f3 <- ggplot(data = dat3, aes(
  x = month,
  y = y,
  colour = as.factor(group),
  group = id
)) + geom_line(lwd = .25) +
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) +
  theme(legend.position = "bottom") +
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment")) +
  xlab("Month") + ylab("FTLD CDR SB score") + ggtitle("30% gradual reduction")
return(list(dat0, dat1, dat2, dat3))
}
