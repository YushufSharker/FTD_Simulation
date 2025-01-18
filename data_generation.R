# Author: Yushuf Sharker
# Data generation for treatment and placebo
# This model uses the PMRM code to generate data

sd = c(1, 1.5, 2, 2.5, 3)
beta = c (5,12,5)
delta1 = 0.3 # 30% reduction for times >0
delta2 = 0.3 # 30% slower progression
delta3 = 0.3 # 30% gradual reduction

datgen <- function(sd=c(1, 2, 2, 2, 2), delta1 = .3, delta2 = .3, beta = c (5,13,5)){

corr = matrix(c(1,    0.65, 0.40, 0.25, 0.15,
                0.65, 1,    0.65, 0.40, 0.25,
                0.40, 0.65, 1,    0.65, 0.40,
                0.25, 0.40, 0.65, 1,    0.65,
                0.15, 0.25, 0.40, 0.65, 1   ),  nrow = 5, byrow = TRUE)
cov = diag(sd) %*% corr %*% diag(sd)
error <-as.vector(matrix(t(mvtnorm::rmvnorm(n, mean = rep(0, m), sigma = cov))))

dat <- placeb_model(M = c(0,6,12,18,24), beta = beta)%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
ggplot(data = dat, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()

#Null
dat0 <- dat %>% mutate(y = fixef0 + error) %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

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

# f1 <- ggplot(data = dat1, aes(
#   x = month,
#   y = y,
#   colour = as.factor(group),
#   group = id
# )) + geom_line(lwd = .25) +
#   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom") +
#   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment")) +
#   xlab("Month") + ylab("FTLD CDR SB score") + ggtitle("30% Proportional reduction")

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
  mutate(chg = chg * (1 - spline(
    x = c(min(M), max(M)),
    y = c(.001, delta3),
    method = "natural",
    xout = c(0, 6, 12, 18, 24)
  )$y * group)) %>%
  group_by(id) %>%
  mutate(fixef0 = chg + fixef0[1L]) %>%
  # mutate(y = fixef0*(1-spline(x = c(min(M), max(M)), y=c(.001, delta2),
  #       method = "natural", xout = c(0, 6, 12, 18, 24))$y*group)) %>%
  ungroup() %>%
  mutate(y = fixef0 + error) %>%
  mutate(y = plyr::round_any(y, 0.5)) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))

# f3 <- ggplot(data = dat3, aes(
#   x = month,
#   y = y,
#   colour = as.factor(group),
#   group = id
# )) + geom_line(lwd = .25) +
#   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) +
#   theme(legend.position = "bottom") +
#   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment")) +
#   xlab("Month") + ylab("FTLD CDR SB score") + ggtitle("30% gradual reduction")
return(list(dat0, dat1, dat2, dat3))
}
