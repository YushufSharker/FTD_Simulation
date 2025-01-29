#' Generate a placebo model dataset
#'
#' This function generates a dataset for a placebo-controlled study with specified parameters.
#'
#' @dpendence tidyverse and splines
#' @param M A numeric vector representing the months or time points.
#' @param n_pbo An integer specifying the number of subjects in the placebo group. Default is 40.
#' @param n_act An integer specifying the number of subjects in the active treatment group. Default is 80.
#' @param ncs_df An integer specifying the degrees of freedom for the natural cubic spline. Default is 2.
#' @param beta A numeric vector of coefficients for the fixed effects.
#' @return A tibble containing the generated dataset with model matrix and additional variables.
#' @export
#'#' @examples
#' M <- 1:12
#' beta <- c(0.5, 1.2, -0.8)
#' placeb_model(M, n_pbo = 40, n_act = 80, ncs_df = 2, beta)

placeb_model <- function(M, n_pbo = 40, n_act = 80, ncs_df = 2, beta, jitter_sd = 0.8) {
  n = n_pbo + n_act
  m = length(M)
  visits <- tibble(visNo = 1:m, M)
  subinfo <- tibble(id = 1:n)
  group <- tibble(cbind(id = 1:n, tibble(group = c(rep(0, n_pbo), rep(1, n_act)))))
  visinfo <- expand_grid(id = 1:n, visNo = 1:m) %>% group_by(visNo)
  d00 <- subinfo %>% left_join(visinfo, by = "id") %>% left_join(visits, by = "visNo") %>%
    left_join(group, by = "id")

  dd00 <- d00 %>% mutate(Mcat = as.factor(M), month = M, monthj = M + rnorm(m * n, mean = 0, jitter_sd)) %>%
    mutate(baseline = ifelse(M > 0, 1, 0))
  dd0 <- model.matrix(~ id + visNo + Mcat + ns(dd00$month, 2), dd00) %>% as_tibble() %>% select(-"(Intercept)") %>%
    dplyr::rename(ns1 = `ns(dd00$month, 2)1`, ns2 = `ns(dd00$month, 2)2`) %>%
    mutate(fixef0 = beta[1] + ns1 * beta[2] + ns2 * beta[3]) %>%
    left_join(dd00, by = c('id', 'visNo'))
  return(dd0 = dd0)
}
