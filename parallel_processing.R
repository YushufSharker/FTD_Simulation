# Codes for simulation function run
start_time <- Sys.time()

# Conditional power function (Pokok and Mehta modified code for sim)
source("placebo_model.R")
source("dgms.R")
source("ncs2.R")
source("model_fit.R")


# simulation function
  sim <- function(delta1 = 0.3,
                             delta2 = 0.3,
                             delta3 = 4,
                             b0 = 3.4,
                             b1 = 11.3,
                             b2 = 3.7,
                             n_pbo = 40,
                             n_act = 80,
                             sd1 = 2,
                             sd2 = 2,
                             sd3 = 2,
                             sd4 = 2,
                             sd5 = 2) {
  #Data generation
  sim_data <- dgm(delta1 = delta1, delta2 = delta2, delta3 = delta3, b0 = b0, b1 = b1, b2 = b2,
      n_pbo = n_pbo, n_act = n_act, sd1 = sd1, sd2 = sd2, sd3 = sd3, sd4 = sd4, sd5 = sd5)
  #Model fitting

  sim_out <- bind_rows(
  do.call(rbind, lapply(sim_data, function(i) mmrm_fit(i))),
  do.call(rbind, lapply(sim_data, function(i) Func_ncs_int(i)))
  )
  return(sim_out)
}

sim.grid <- expand.grid(b0 = 4,
                        b1 = seq(5, 14, length.out = 5),
                        b2 = seq(5, 7, length.out = 3),
                        delta1 = 0.3,
                        delta2 = 0.3,
                        delta3 = 4, # to expected reduction at the end of the time period
                        sd1 = 1,
                        sd2 = 2,
                        sd3 = 2,
                        sd4 = 2,
                        sd5 = 2,
                        npbo = 40,
                        nact = 80
)

# Parallel processing function
source("goparallel.R")
goparallel(ncores = 10)

#### Run Parallel processing
parallel::clusterExport(cl = cl,
                        varlist = c("placeb_model",
                                    "dgm",
                                    "mmrm_fit",
                                    "Func_ncs_int",
                                    "sim.grid",
                                    "sim")
                        )

parallel::clusterEvalQ(cl = cl,
                       expr = {
                         require(mvtnorm)
                         require(lme4)
                         require(lmertest)
                         require(emmeans)
                         require(broom)
                         require(splines)
                         require(mmrm)
                         require(truncnorm)
                         require(plyr)
                         require(tmvtnorm)}
  )
sim_FTD_output<- bind_rows( #future_map_dfr under purrr and furrr lib
  parallel::parApply(
  cl = cl,
  X = matrix(1:nrow(sim.grid)),
  MARGIN = 1,
  FUN = function(x) {
    lapply(1:5000, function(y) { # map_dfr
      # mu_t=sim.grid$mu_c[x]-sim.grid$mu_c[x]*sim.grid$expected_reduction[x]
      bind_rows(c(
        sim(
          b0 = sim.grid$b0[x],
          b1 = sim.grid$b1[x],
          b2 = sim.grid$b2[x],
          delta1 = sim.grid$delta1[x],
          delta2 = sim.grid$delta2[x],
          delta3 = sim.grid$delta3[x],
          sd1 = sim.grid$sd1[x],
          sd2 = sim.grid$sd2[x],
          sd3 = sim.grid$sd3[x],
          sd4 = sim.grid$sd4[x],
          sd5 = sim.grid$sd5[x],
          n_pbo = sim.grid$npbo[x],
          n_act = sim.grid$nact[x]
        )
      ))
    })
  }
))

parallel::stopCluster(cl)
end_time <- Sys.time()
runtime <- (end_time - start_time) # in hours
save.image("./Outputs/tak594_FTD_SIM_01192025.RData")



