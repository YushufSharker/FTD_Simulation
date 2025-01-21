# model fitting
# Model Fits mmrm
mmrm_fit <- function(data = dat3) {
  fit <- mmrm::mmrm(
    formula = chg ~ group * Mcat + ar1(Mcat | id),
    data = data,
    control = mmrm_control(method = "Kenward-Roger")
  )

  prop_change_mmrm <- tibble(as.data.frame(emmeans(fit, ~ group |
                            Mcat, type = "response")[9:10])) %>%
    with(., round(1 - .[2, 3] / .[1, 3], 2)) #%>% dplyr::rename(c("Prop" = "emmean"))
  estimate = round(summary(fit)[["coefficients"]][10, 1], 2)
  pvalue <-  round(summary(fit)[["coefficients"]][10, 5], 5)
  model = 'MMRM'
  dgm = unique(data["dgm"])
  errm = unique(data["errm"])
  return(tibble(prop_change_mmrm, estimate, pvalue, model, dgm, errm))
}

#Ref resource for MMRM https://cran.r-project.org/web/packages/mmrm/vignettes/introduction.html
