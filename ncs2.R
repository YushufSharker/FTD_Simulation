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
                        (I(ns21(M))) +
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
        p.value_ratio = round(p.value, 5),
        mod = 'ncs-ranslp',
        #MSE=MSE,
        estimate.score = estimate,
        # SE.score = SE,
        # df.score = df,
        #t.ratio.score = t.ratio,
        p.value.score = round(p.value, 2)
      )
  )
  output
}



# Adjusted Scenarios
# Making the difference at the last visit almost the same

results_path <- "sim_results_same_6"

# check for existing results ----
LF <- list.files(file.path(results_path))
completed <- as.numeric(
        gsub('_', '',
             gsub('.rdata', '',
                  gsub('sim', '', LF))))
TODO <- setdiff(1:NSIMS, completed)
## spline functions ----
ns21 <- function(t){
        as.numeric(predict(splines::ns(dd0$Month, df=2), t)[,1])
}
ns22 <- function(t){
        as.numeric(predict(splines::ns(dd0$Month, df=2), t)[,2])
}

# non-linear functions

# run simulations ----
parallel::mclapply(rev(TODO), mc.cores = MC.CORES, FUN = function(i){
        print(i)
        set.seed(SEEDS[i])
        resids_w <- mvtnorm::rmvnorm(N,sigma=Sigma)
        colnames(resids_w) <- 1:ncol(resids_w)
        resids <- resids_w %>%
                as_tibble() %>%
                mutate(id = 1:N, active = sample(0:1, size = N, replace=TRUE)) %>%
                pivot_longer('1':'10', names_to = 'visNo', values_to = 'residual') %>%
                mutate(
                        visNo = as.numeric(visNo),
                        Active = as.factor(active))
        dd1 <- dd0 %>%
                left_join(resids, by=c('id', 'visNo')) %>%
                filter(visNo %in% c(1,3,5,7,9,10)) %>%
                mutate(
                        Y_type_1 = fixef0 + residual,
                        # stable
                        Y_power_stable = Y_type_1 + active*( (yrs>=1.5)*spline(x = c(1.5,2.5,4,4.5), y = c(0,0.06,0.2,0.2), method = "natural",xout=yrs)$y),
                        # fading
                        Y_power_fading = Y_type_1 + active*( (yrs>=1.5)*spline(x = c(1.5,2.5,3.5,4.5), y = c(0,0.08,0.22,0.2), method = "natural",xout=yrs)$y),
                        # 30 pct reduction
                        Y_power_30pctr = Y_type_1 + active * -0.3 * (adni_ns_fun[[1]](yrs)*0.04380665 +
                                                                             adni_ns_fun[[2]](yrs)*-0.4601309 +
                                                                             adni_ns_fun[[3]](yrs)*-2.232262 +
                                                                             adni_ns_fun[[4]](yrs)*-3.509172),
                        # 10 pct progression delay
                        Y_power_10pctd = Y_type_1 + active *  (adni_ns_fun[[1]](yrs*0.9)*0.04380665 +
                                                                       adni_ns_fun[[2]](yrs*0.9)*-0.4601309 +
                                                                       adni_ns_fun[[3]](yrs*0.9)*-2.232262 +
                                                                       adni_ns_fun[[4]](yrs*0.9)*-3.509172 - (adni_ns_fun[[1]](yrs)*0.04380665 +
                                                                                                                      adni_ns_fun[[2]](yrs)*-0.4601309 +
                                                                                                                      adni_ns_fun[[3]](yrs)*-2.232262 +
                                                                                                                      adni_ns_fun[[4]](yrs)*-3.509172) ),

                        # 5.5 month delay
                        Y_power_5.5md = Y_type_1 + active *  (adni_ns_fun[[1]](yrs-5.5/12)*0.04380665 +
                                                                      adni_ns_fun[[2]](yrs-5.5/12)*-0.4601309 +
                                                                      adni_ns_fun[[3]](yrs-5.5/12)*-2.232262 +
                                                                      adni_ns_fun[[4]](yrs-5.5/12)*-3.509172 - (adni_ns_fun[[1]](yrs)*0.04380665 +
                                                                                                                        adni_ns_fun[[2]](yrs)*-0.4601309 +
                                                                                                                        adni_ns_fun[[3]](yrs)*-2.232262 +
                                                                                                                        adni_ns_fun[[4]](yrs)*-3.509172)                                                            )
                )
        dd1$M.n <- as.numeric(levels(dd1$M))[dd1$M]
        # Stable
        dd1$Y <- dd1$Y_power_stable
        Stable = Func_ncs_int(dd1)
        # Fading
        dd1$Y <- dd1$Y_power_fading
        Fading = Func_ncs_int(dd1)
        # 30% reduction
        dd1$Y <- dd1$Y_power_30pctr
        Reduc30pct = Func_ncs_int(dd1)
        # 10% progression delay
        dd1$Y <- dd1$Y_power_10pctd
        Delay10pct = Func_ncs_int(dd1)
        # 5.5 months delay
        dd1$Y <- dd1$Y_power_5.5md
        Delay5.5m = Func_ncs_int(dd1)
        # Type I error --  No trt effect
        dd1$Y <- dd1$Y_type_1
        TYPEI = Func_ncs_int(dd1)

        res <- bind_rows(
                bind_rows(Stable) %>% mutate(scenario='Stable'),
                bind_rows(Fading) %>% mutate(scenario='Fading'),
                bind_rows(Reduc30pct) %>% mutate(scenario='Reduc30pct'),
                bind_rows(Delay10pct) %>% mutate(scenario='Delay10pct'),
                bind_rows(Delay5.5m) %>% mutate(scenario='Delay5.5m'),
                bind_rows(TYPEI) %>% mutate(scenario='type I')) %>%
                mutate(sim = i)
        save(res, file = file.path(results_path, paste0('sim_', i, '.rdata')))
        return(NULL)
})


filenames <- list.files(file.path(results_path), full.names = TRUE)
summary0 <- expand.grid(mod=c("cat", "prop", "time-pmrm", "pstime-pmrm", "ncs-uns", "ncs-ranslp","ps_pncs_int","time_pncs_int"),
                        scenario = c("Stable","Fading","Reduc30pct","Delay10pct", "Delay5.5m", "No Effect"))  %>%
        mutate(true.effect = case_when(
                scenario == "Stable" ~ 0.2,
                scenario == "Fading" ~ 0.2,
                scenario == "Reduc30pct" ~ 0.1958384,
                scenario == "Delay10pct" ~ 0.1976384,
                scenario == "Delay5.5m" ~ 0.200935,
                scenario == "type I" ~ 0))

estimate.score <- c()
Uncovergence <- c()
type2 <- c()
type2.error0 <- c()
for (j in filenames){
        load(j)
        res <- res %>%
                mutate(Uncovergence = as.numeric(is.na(p.value)),
                       p.value.1 =  case_when(
                               mod == "time_pncs_int" ~ 2 * pt(-abs(t.ratio), df = df.score),
                               TRUE ~ p.value
                       ),
                       Type2.error = case_when(
                               Uncovergence == 1 ~ 0,
                               TRUE ~ as.numeric(p.value.1 > .05)),
                       Type2.error0 = as.numeric(p.value.1 > .05)

                )
        Uncovergence <- cbind(Uncovergence, res$Uncovergence)
        type2 <- cbind(type2, res$Type2.error)
        type2.error0 <- cbind(type2.error0, res$Type2.error0)
        estimate.score <- cbind(estimate.score, res$estimate.score)

}

# for (i in filenames){
#         load(i)
#         res <- res %>%
#                 mutate(Uncovergence = as.numeric(is.na(p.value.score)),
#                        Type2.error = case_when(
#                                Uncovergence == 1 ~ 0,
#                                TRUE ~ as.numeric(p.value.score > .05)),
#                        Type2.error0 <- as.numeric(p.value.score > .05)
#
#                 )
#         Uncovergence <- cbind(Uncovergence, res$Uncovergence)
#
#         type2 <- cbind(type2, res$Type2.error)
#         type2.error0 <- cbind(type2.error0, res$Type2.error0)
#         estimate.score <- cbind(estimate.score, res$estimate.score)
#
# }



summary <- summary %>%
        mutate(Sample.Size = N,
               N.visit = 6,
               N.Uncovergence = apply(Uncovergence,1,sum),
               N.type2 = apply(type2,1,sum),
               type2.error = apply(type2.error0,1,function(x) mean(x, na.rm = TRUE)),
               power = 1-type2.error,
               mean.score.est = apply(estimate.score,1,function(x) mean(x, na.rm = TRUE)),
               abs.bias = abs(mean.score.est - true.effect))
