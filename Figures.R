# To create different Figures
# Yushuf Sharker
# 2025-1-17
library(tidyverse)

dat <- placeb_model(M = c(0,6,12,18,24), beta = c (2.5,9,6))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
ggplot(data = dat, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()
dat %>% group_by(M) %>% dplyr::summarize(mean(fixef0))

figdata <- dgm (delta1 = 0.3, delta2 = 0.3, delta3 = 4, b0 = 3.4, b1=11.3, b2 = 3.7,
                 n_pbo = 40, nact=80,
                 sd1 = 2, sd2 = 2, sd3 = 2, sd4 = 2, sd5 = 2  )

# Figure of placebo data origin
Placebo_model <- placeb_model(M = c(0,6,12,18,24), beta = c (2.5,14,5))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
ggplot(data = dat, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()
dat %>% group_by(M) %>% dplyr::summarize(mean(fixef0))

# Under Null Scenario
f0<- ggplot(data = dat0, aes(x = month, y = y, colour = as.factor(group), group = id)) +geom_line(lwd = .25)+
   geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom")+
   scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
   xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("Null Scenario")

# Under 30% Proportional reduction
f1<-ggplot(data = dat1, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd=.25)+
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% Proportional reduction")

# under 30% Slower Progression
f2<-ggplot(data = dat22, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd = .25)+
geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2)  +
  theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% Slower Progression")


# Linear improvement to 30 over time
f3<-ggplot(data = dat3, aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd=.25, alpha = .25)+
geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) +
  theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% at end")

