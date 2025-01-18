# Codes for simulation function run
library(tidyverse)
# Conditional power function (Pokok and Mehta modified code for sim)
source("placebo_model.R")
source("dgms.R")
source("ncs2.R")
source("model_fit.R")

dgm <- function(delta1 = 0.3, delta2 = 0.3, delta3 = 4, b0 = 3.4, b1=11.3, b2 = 3.7,
                n_pbo = 40, nact=80,
                sd1 = 2, sd2 = 2, sd3 = 2, sd4 = 2, sd5 = 2  )

# simulation function
  sim <- function(delta1 = 0.3,
                             delta2 = 0.3,
                             delta3 = 4,
                             b0 = 3.4,
                             b1 = 11.3,
                             b2 = 3.7,
                             n_pbo = 40,
                             nact = 80,
                             sd1 = 2,
                             sd2 = 2,
                             sd3 = 2,
                             sd4 = 2,
                             sd5 = 2) {
  #Data generation
  sim_data <- dgm(delta1 = delta1, delta2 = delta2, delta3 = delta3, b0 = b0, b1 = b1, b2 = b2,
      n_pbo = n_pbo, nact = nact, sd1 = sd1, sd2 = sd2, sd3 = sd3, sd4 = sd4, sd5 = sd5)
  do.call(rbind, lapply(dat, function(i) mmrm_fit(i)))
  }













  #delta <- rnorm(1, mean = mdelta, sd = sqrt(sd_t ^ 2/(f*n2)*r + sd_c ^ 2/(f*n2)*(1-r) ))
  mu_t <- mu_c-mu_c*expected_reduction
  mdelta = mu_c*expected_reduction
  n1 <- f*n2
  control <- rnorm(n1*(1-r), mean = mu_c, sd =sd_c)
  treat <- rnorm(n1*r, mean = mu_t, sd = sd_t)



  delta <- mean(control)-mean(treat)
  sd_t_hat <- sd(treat)
  sd_c_hat <- sd(control)

  x <- CP_PZ_sim(
      r = r,
      n1 = n1,
      n2 = n2,
      alpha_1s = alpha_1s,
      eff_est = delta,
      eff_planned = eff_planned,
      eff_null = eff_null,
      SE = NULL,
      p_c = NULL,
      sd_t = sd_t_hat,
      sd_c = sd_c_hat,
      type = "cont",
      pow = pow,
      f = f,
      max = max,
      plot_effect = 0
    )
  nmax = x$n2 * max

  final_n_increase <- dplyr::case_when( # use as_tibble and exclude $ signs
    x$CP_obs < fr ~ 0,
    x$CP_obs >= fr & x$CP_obs < x$CP_ll$CP ~ (x$n2 - x$n1),
    (x$CP_obs >= x$CP_ll$CP & x$CP_obs < x$CP_ul$CP) ~ min(x$n2_inc_new_obs, nmax-x$n1, na.rm = TRUE),
    x$CP_obs >= x$CP_ul$CP ~ (x$n2 -  x$n1)
  )

  Blind_n_increase <- dplyr::case_when( # use as_tibble and exclude $ signs
    x$CP_obs < fr ~ 0,
    x$CP_obs >= fr & x$CP_obs < x$CP_ll$CP ~ (x$n2 - x$n1),
    (x$CP_obs >= x$CP_ll$CP & x$CP_obs < x$CP_ul$CP) ~ (nmax-x$n1),
    x$CP_obs >= x$CP_ul$CP ~ (x$n2 -  x$n1)
  )
  zones <- dplyr::case_when(
    x$CP_obs < fr ~ "Futility",
    x$CP_obs >= fr & x$CP_obs < x$CP_ll$CP ~ "Unfavorable",
    (x$CP_obs >= x$CP_ll$CP & x$CP_obs < x$CP_ul$CP) ~ "Promising",
    x$CP_obs >= x$CP_ul$CP ~ "Favorable"
  )

  fstar <- x$n1 / x$n2
  # simulate z2 for CHW statistic
  if (zones !="Promising") {
    n2_t <- (x$n2 - x$n1) * r
    n2_c <- (x$n2 - x$n1) * (1 - r)
  }
  if (zones =="Promising") {
    n2_t <- final_n_increase * r
    n2_c <- final_n_increase * (1 - r)
  }

  SE2 <- sqrt(sd_t ^ 2 / n2_t + sd_c ^ 2 / n2_c)
  delta2 <- rnorm(1, mean = mdelta, sd = SE2)
  z2     <- (delta2    - eff_null) / SE2
  z <- sqrt(fstar) * x$z1 + sqrt(1 - fstar) * z2
  pvalue_end_CHW <- pnorm(z, 0, 1, lower.tail = FALSE)

  #Simulate z2 conventional
  if (zones !="Promising") {
    n_t <-x$n2 *r
    n_c <-x$n2*(1-r)
  }
  if (zones =="Promising") {
    n_t <-(x$n1+final_n_increase) *r
    n_c <-(x$n1+final_n_increase)*(1-r)
  }
  SE_t <- sqrt(sd_t ^ 2 / n_t + sd_c ^ 2 / n_c)
  delta_t <- rnorm(1, mean = mdelta, sd = SE_t)
  z_t     <- (delta_t  - eff_null) / SE_t
  pvalue_end_conv <- pnorm(z_t, 0, 1, lower.tail = FALSE)

  # Simulate z under blinding condition
  if (zones !="Promising") {
    n_t_b <-x$n2 *r
    n_c_b <-x$n2*(1-r)
  }
  if (zones =="Promising") {
    n_t_b <-(x$n1+Blind_n_increase) *r
    n_c_b <-(x$n1+Blind_n_increase)*(1-r)
  }
  SE_B <- sqrt(sd_t ^ 2 / n_t_b + sd_c ^ 2 / n_c_b)
  delta_B <- rnorm(1, mean = mdelta, sd = SE_B)
  z_B     <- (delta_B - eff_null) / SE_B
  pvalue_B <- pnorm(z_B, 0, 1, lower.tail = FALSE)

  return(
    c(
      eff_o = x$eff_est,
      obs.cp = x$CP_obs,
      zones = zones,
      treat_ratio = r,
      n2star = x$n2_inc_new_obs,
      n1 = x$n1,
      n2 = x$n2,
      nmax = nmax,
      final_n_increase = final_n_increase,
      Blind_n_increase = Blind_n_increase,
      lowermargin_est = x$CP_ll$CP,
      uppermargin_est = x$CP_ul$CP,
      z1 =x$z1,
      delta = delta,
      tdelta = mdelta,
      SE1 = x$se,
      z2 = z2,
      delta2 = delta2,
      SE2 = SE2,
      z_end_chw = z,
      pvalue_end_CHW = pvalue_end_CHW,
      z_conv = z_t,
      pvalue_end_conv = pvalue_end_conv,
      z_B = z_B,
      pvalue_B = pvalue_B
    )
  )
}


sim.grid <- expand.grid(mu_c =c(7, 6, 5, 4, 3, 2 ),
                        expected_reduction = c(0, .3, .35, .4, .45),
                        f = c(0.5), # Fraction of sample at interim
                        fr = c(.05), # futility boundary
                        r = c(.5, 0.6666667), #Allocation Ratio
                        n2 = c(60, 70),
                        max = c(1.666667, 1.428571)
)%>%
  mutate(test = round(max*n2, 0)) %>%
  filter(test %in% c(100))



# Parallel processing function
source("goparallel.R")
goparallel(ncores = 10)

#### Run Parallel processing
parallel::clusterExport(cl = cl,
                        varlist = c("sim_cond.power", "CP_PZ_sim", "sim.grid"))
sim_output4_2<- bind_rows( #future_map_dfr under purrr and furrr lib
  parallel::parApply(
  cl = cl,
  X = matrix(1:nrow(sim.grid)),
  MARGIN = 1,
  FUN = function(x) {
    lapply(1:5000, function(y) { # map_dfr
      # mu_t=sim.grid$mu_c[x]-sim.grid$mu_c[x]*sim.grid$expected_reduction[x]
      bind_rows(c(
        sim_cond.power(
          mu_c = sim.grid$mu_c[x],
          expected_reduction = sim.grid$expected_reduction[x],
          f = sim.grid$f[x],
          fr = sim.grid$fr[x],
          r = sim.grid$r[x],
          n2 = sim.grid$n2[x],
          max = sim.grid$max[x]

        ),
        mu_t=sim.grid$mu_c[x]-sim.grid$mu_c[x]*sim.grid$expected_reduction[x],
        f = sim.grid$f[x],
        fr = sim.grid$fr[x]
      ))
    })
  }
))
parallel::stopCluster(cl)

save.image("./Output/tak594_workspace_sim04_2_05062024.RData") # November 15



