# To create different Figures
# Yushuf Sharker
# 2025-1-17
library(tidyverse)

dat <- placeb_model(M = c(0,6,12,18,24), beta = c (4,15,5))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
p <-
  ggplot(data = dat, aes(x = month, y = fixef0, group = as.factor(group))) +
  geom_smooth(aes(), se = FALSE, lwd = .5, color = 'black')+theme_bw()+
  xlab("Month") + ylab("FTLD CDR SB score")

for (i in 1:nrow(sim.grid)){
  dat2 <- placeb_model(M = c(0,6,12,18,24),
                       beta = c (sim.grid$b0[i],sim.grid$b1[i],sim.grid$b2[i]))%>%
    mutate(id = as.factor(i))%>%
    group_by(id) %>%
    mutate(chg = fixef0 - fixef0[1L]) %>%
    ungroup()
p = p+geom_smooth(aes(color = id), data = dat2, se = FALSE, lwd = .5)

}

sim.grid <- expand.grid(b0 = 4,
                        b1 = seq(7, 25, by = 2),
                        b2 = seq(7, 9, by = 2))

sim.grid <- sim.grid %>% filter(!(b1 %in% c(17:25) & b2 ==5 )) %>%
  filter(!(b1 %in% c(19:25) & b2 ==7 ))

ggsave("./Outputs/PlaceboModels_test.png", plot = p, width = 5, height = 4, dpi = 300)

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

