# model fitting
# Model Fits mmrm
mmrm_fit <- function(dat = dat3) {
  fit <- mmrm::mmrm(
    formula = chg ~ group * Mcat + ar1(Mcat | id),
    data = dat,
    control = mmrm_control(method = "Kenward-Roger")
  )

  prop_change_mmrm <- tibble(as.data.frame(emmeans(fit, ~ group |
                            Mcat, type = "response")[9:10])) %>%
    with(., round(1 - .[2, 3] / .[1, 3], 2)) %>% rename(Prop = emmean)
  estimate = round(summary(fit)[["coefficients"]][10, 1], 2)
  pvalue <-  round(summary(fit)[["coefficients"]][10, 5], 5)
  model = 'MMRM'
  return(tibble(prop_change_mmrm, estimate, pvalue, model))
  #return(tibble(c(prop_change = prop_change_mmrm$Prop, estimate = estimate, pvalue = pvalue, mod = 'MMRM')))
}

