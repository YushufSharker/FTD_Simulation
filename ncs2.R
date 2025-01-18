# Vahe Khachadourian and his intern wrote the code as a part of the intern project
# Yushuf Modified the code to fit the purpose
# Date: 12-19-2024
# Ref: https://www.researchgate.net/publication/362759985_Natural_cubic_splines_for_the_analysis_of_Alzheimer's_clinical_trials/fulltext/62fdaf84eb7b135a0e415754/Natural-cubic-splines-for-the-analysis-of-Alzheimers-clinical-trials.pdf?origin=scientificContributions
require(lme4)
require(tidyverse)
require(lmertest)

ns21 <- function(t){
  as.numeric(predict(splines::ns(dat$month, df=2), t)[,1])
}
ns22 <- function(t){
  as.numeric(predict(splines::ns(dat$month, df=2), t)[,2])
}


Func_ncs_int <- function(data = dat1, last_visit = 24) {
    fit_ncs_int <- lmer(chg ~
                        I(ns21(M)) +
                        I(ns22(M)) +
                        (I(ns21(M))):group +
                        I(ns22(M)):group +
                        (1 | id), data = data)

  # Extract the fitted values and residuals
  dds <- dat1 %>%
    mutate(fitted_values = fitted(fit_ncs_int),
           residuals = chg - fitted_values)

  # Identify the residuals corresponding to M24
  dd1_M24 <- dds %>% filter(M == last_visit)  # Assuming 'M' is the column indicating the month/time point

  # Calculate the squared residuals
  dd1_M24 <- dd1_M24 %>%
    mutate(squared_residuals = residuals ^ 2)

  # Compute the mean square error (MSE) for M54
  MSE <- mean(dd1_M24$squared_residuals)

  out_ncs_ranslp <- ref_grid(
    fit_ncs_int,
    at = list(M = 24, group = c("0", "1")),
    data = data,
    mode = "satterthwaite",
    lmerTest.limit = 15000
  ) %>%
    emmeans(specs = 'group',
            by = 'M',
            lmerTest.limit = 15000) #%>%
  prop_change_cs2 <- out_ncs_ranslp %>% as.data.frame() %>% tibble() %>%
    with(., round(1 - .[2, 3] / .[1, 3], 2))
  names(prop_change_cs2) <- "Prop"

  output <- bind_cols(
    prop_change_cs2,
    out_ncs_ranslp %>% pairs(reverse = TRUE) %>% as.data.frame() %>%
      mutate(
        #p.value_ratio = round(p.value, 5),
        mod = 'ncs-ranslp',
        #MSE=MSE,
        # SE.score = SE,
        # df.score = df,
        #t.ratio.score = t.ratio,
        p.value.score = round(p.value, 2)
      )%>% select(estimate, SE, p.value)
  )
  return(output)
}

