# To create different Figures
# Yushuf Sharker
# 2025-1-17

dat <- placeb_model(M = c(0,6,12,18,24), beta = c (2.5,9,6))%>%
  group_by(id) %>%
  mutate(chg = fixef0 - fixef0[1L]) %>%
  ungroup()
ggplot(data = dat, aes(x = month, y = fixef0, colour = as.factor(group), group = id)) +
  geom_smooth(aes(group = as.factor(group)), se = FALSE, lwd = 2) + theme_minimal()
dat %>% group_by(M) %>% dplyr::summarize(mean(fixef0))
