# model fitting
# Model Fits mmrm
mmrm_fit <- function(data = dat3) {
  fit <- mmrm::mmrm(
    formula = y ~ group * Mcat + us(Mcat | id),
    data = data,
    control = mmrm_control(method = "Kenward-Roger", accept_singular = TRUE)
  )
  # Calculate the condition number
  condition_number <- kappa(model.matrix(fit))

  sing <- ifelse(condition_number > 1e+04, "Singular", "Non-Singular")
  wtext <- names(warnings())
  w1 = get_phrase(wtext, "unidentifiable")
  w2 = get_phrase(wtext, "failed to converge")

  prop_change_mmrm <- tibble(as.data.frame(emmeans(fit, ~ group |
                                                     Mcat, type = "response")[9:10])) %>%
    with(., round(1 - .[2, 3] / .[1, 3], 2)) #%>% dplyr::rename(c("Prop" = "emmean"))
  estimate = round(summary(fit)[["coefficients"]][10, 1], 2)
  pvalue <-  round(summary(fit)[["coefficients"]][10, 5], 5)
  model = 'MMRM'
  dgm = unique(data["dgm"])
  errm = unique(data["errm"])
  misrate = unique(data["misrate"])
  return(tibble(
    dgm,
    errm,
    prop_change_mmrm,
    estimate,
    pvalue,
    model,
    sing,
    misrate,
    w1,
    w2
  ))
}

#Ref resource for MMRM https://cran.r-project.org/web/packages/mmrm/vignettes/introduction.html
