# Vahe Khachadourian and his intern wrote the code as a part of the intern project
# Yushuf Modified the code to fit the purpose
# Dependency package spline, tidyverse, lme4, lmerTest, emmeans, mmrm
# Date: 12-19-2024
# Ref: https://www.researchgate.net/publication/362759985_Natural_cubic_splines_for_the_analysis_of_Alzheimer's_clinical_trials/fulltext/62fdaf84eb7b135a0e415754/Natural-cubic-splines-for-the-analysis-of-Alzheimers-clinical-trials.pdf?origin=scientificContributions

# Here Month is consisered discrete month

Func_ncs_int <- function(data = dat1, last_visit = 24) {

  ns21 <- function(t){
    as.numeric(predict(splines::ns(data[["M"]], df=2), t)[,1])
  }
  ns22 <- function(t){
    as.numeric(predict(splines::ns(data[["M"]], df=2), t)[,2])
  }

    fit_ncs_int <- lmer(y ~
                          I(ns21(M)) +
                          I(ns22(M))+(I(ns21(M)) +
                         I(ns22(M))):group +   (0 + Mcat| id),
                        data = data,
                        control = lmerControl(check.nobs.vs.nRE = "ignore") )
#0 + Mcat
  sing <- ifelse(isSingular(fit_ncs_int) == TRUE, "Singular", "Non-Singular")
  wtext <- names(warnings())
  w1 = get_phrase(wtext, "unidentifiable")
  w2 = get_phrase(wtext, "failed to converge")
        # # Extract the fitted values and residuals
        # dds <- data %>%
        #   mutate(fitted_values = fitted(fit_ncs_int),
        #          residuals = chg - fitted_values)
        #
        # # Identify the residuals corresponding to M24
        # dd1_M24 <- dds %>% filter(M == last_visit)  # Assuming 'M' is the column indicating the month/time point
        #
        # # Calculate the squared residuals
        # dd1_M24 <- dd1_M24 %>%
        #   mutate(squared_residuals = residuals ^ 2)
        #
        # # Compute the mean square error (MSE) for M54
        # MSE <- mean(dd1_M24$squared_residuals)

  out_ncs_ranslp <- ref_grid(
    fit_ncs_int,
    at = list(M = 24, group = c("0", "1")),
    data = data,
    mode = "kenward-roger") %>%
    emmeans(specs = 'group',
            M=24,
            lmerTest.limit = 2000) #%>%

  prop_change_cs2 <- out_ncs_ranslp %>% as.data.frame() %>% tibble() %>%
    with(., round(1 - .[2, 2] / .[1, 2], 2)) #%>% dplyr::rename(c("Prop"="emmean"))

  dgm = unique(data["dgm"])
  errm = unique(data["errm"])
  misrate = unique(data["misrate"])
  w1 = get_phrase(wtext, "unidentifiable")
  w2 = get_phrase(wtext, "failed to converge")

  output <- bind_cols(
    dgm, errm, misrate, prop_change_cs2,
    out_ncs_ranslp %>% pairs(reverse = TRUE) %>% as.data.frame() %>%
      mutate(
        estimate = round(estimate, 2),
        p.value_ratio = round(p.value, 5),
        mod = 'ncs-ranslp',
        #MSE=MSE,
        # SE.score = SE,
        # df.score = df,
        #t.ratio.score = t.ratio,
        #p.value.score = round(p.value, 5),
        sing = sing,
        w1 = w1,
        w2 = w2
      )%>% select(estimate, pvalue = p.value_ratio, model = mod, sing,
                  w1, w2)

  )
  return(output)
}

