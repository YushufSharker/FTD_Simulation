# To create different Figures
# Yushuf Sharker
# 2025-1-17
library(tidyverse)
library(emmeans)
library(broom)
library(splines)
library(mmrm)
library(truncnorm)
library(tmvtnorm)
library(lme4)
library(gridExtra)
dat <- placeb_model(M = c(0,6,12,18,24), beta = c (4,7,7))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
p <-
  ggplot(data = dat, aes(x = month, y = fixef0, group = as.factor(group))) +
  geom_smooth(aes(), se = FALSE, lwd = .5, color = 'black')+theme_bw()+
  xlab("Month") + ylab("FTLD CDR SB score")

for (i in 2:nrow(sim.grid_f)){
  dat2 <- placeb_model(M = c(0,6,12,18,24),
                       beta = c (sim.grid$b0[i],sim.grid$b1[i],sim.grid$b2[i]))%>%
    mutate(id = as.factor(i))%>%
    group_by(id) %>%
    mutate(chg = fixef0 - fixef0[1L]) %>%
    ungroup()
p = p+geom_smooth(aes(color = id), data = dat2, se = FALSE, lwd = .5)

}

sim.grid_f <- expand.grid(b0 = 4,
                        b1 = seq(7, 25, by = 2),
                        b2 = seq(7, 9, by = 2))

sim.grid <- sim.grid %>% filter(!(b1 %in% c(17:25) & b2 ==5 )) %>%
  filter(!(b1 %in% c(19:25) & b2 ==7 ))

ggsave("./Outputs/PlaceboModels_test.png", plot = p, width = 5, height = 4, dpi = 300)

figdata <- dgm (delta1 = 0.3, delta2 = 0.3, delta3 = 4, b0 = 4, b1=25, b2 = 9,
                 n_pbo = 40, n_act=80,
                 sd1 = 2, sd2 = 3, sd3 = 4, sd4 = 5, sd5 = 6  )


# Figure of placebo data origin
Placebo_model <- placeb_model(M = c(0,6,12,18,24), beta = c (4,25,9))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
pb<-ggplot(data = Placebo_model, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) +
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("Placebo model")+theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))

spline_interpolation <- function(df) {
  delta2 = .3
  df= df %>% mutate(x = M * (1 + delta2 * group))%>%
    mutate(y= spline(x = x, y = fixef0, method = "natural", xout = M)$y)
  return(df)
}
dat2i <- do.call(rbind, lapply(split(Placebo_model, f = dat$id), function(i) spline_interpolation(i)))
sp<-dat2i %>% select(-x) %>%
  mutate(y = y + 0,
         dgm = "SP",
         errm = "N",
         misrate = 0) %>%
  mutate(y = plyr::round_any(y, 0.5))%>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id)) %>%
  ggplot(aes(x = month, y = y, colour = as.factor(group), group = id)) +
  geom_line(lwd = .25, alpha = .25)+
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2)  +
  theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% Slower Progression")

tp <- Placebo_model %>%
  mutate(chg = chg - (spline(
    x = c(min(M), max(M)),
    y = c(0, 3),
    method = "natural",
    xout = c(0, 6, 12, 18, 24)
  )$y * group)) %>%
  group_by(id) %>%
  mutate(fixef0 = chg + fixef0[1L]) %>%
  ungroup() %>%
  mutate(y = fixef0 + 0) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>%
  ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id),
         dgm = "TP",
         errm = "N",
         misrate = 0) %>%
  ggplot(aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd=.25, alpha = .25)+
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) +
  theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("3 units at 24 Months")

pr <- Placebo_model%>%
  mutate(chg = chg*(1-.3*group),
                                   dgm = "PR",
                                   errm = "N",
                                   misrate = 0) %>%
  group_by(id) %>% mutate(fixef0 = chg+fixef0[1L]) %>%
  ungroup() %>%
  mutate(y = fixef0 +0 ) %>%
  group_by(id) %>%
  mutate(chg = y - y[1L]) %>% ungroup() %>%
  mutate(group = as.factor(group), id = as.factor(id))%>%
  ggplot(aes(x = month, y = y, colour = as.factor(group), group = id)) +
  geom_line(lwd = .25, alpha = .25)+
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("Null Scenario")

png("./Outputs/combined_plot2.png", width = 500, height = 500)
grid.arrange(pb, pr, sp, tp,  ncol = 2)
dev.off()


# Under Null Scenario
f0<- ggplot(data = figdata[[1]], aes(x = month, y = y, colour = as.factor(group), group = id)) +
  geom_line(lwd = .25, alpha = .15)+
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 1, alpha = 1.5) + theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("Null Scenario")

# Under 30% Proportional reduction
f1<-ggplot(data = figdata[[5]], aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd=.25, alpha = .15)+
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 1, alpha = 1.5) + theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% Proportional reduction")

# under 30% Slower Progression
f2<-ggplot(data = figdata[[7]], aes(x = month, y = y, colour = as.factor(group), group = id)) +
  geom_line(lwd = .25, alpha = .15)+
geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 1, alpha = 1.5)  +
  theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% Slower Progression")


# Linear improvement to 30 over time
f3<-ggplot(data = figdata[[9]], aes(x = month, y = y, colour = as.factor(group), group = id)) + geom_line(lwd=.25, alpha = .15)+
geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 1, alpha = 1.5) +
  theme(legend.position = "bottom")+
  scale_colour_discrete(name = "Group", labels = c("Placebo", "Treatment"))+
  xlab("Month") + ylab("FTLD CDR SB score")+ggtitle("30% at end")

png("./Outputs/combined_dataview.png", width = 500, height = 500)
grid.arrange(f0, f1, f2, f3,  ncol = 2)
dev.off()

