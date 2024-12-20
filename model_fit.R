# model fitting
# Model Fits mmrm
mmrm_fit <- function(dat = dat3) {
  fit <- mmrm::mmrm(
    formula = chg ~ group * Mcat + ar1(Mcat | id),
    data = dat,
    control = mmrm_control(method = "Kenward-Roger")
  )

  Overall_prop_change <- tibble(as.data.frame(emmeans(fit, ~ group |
                            Mcat, type = "response")[9:10])) %>%
    mutate(simno = 1) %>%
    with(., round(1 - .[2, 3] / .[1, 3], 2))
  pvalue <-  round(summary(fit)[["coefficients"]][10, 5], 5)
  return(c(prop_change = Overall_prop_change$emmean, pvalue = pvalue))
}
