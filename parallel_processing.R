# Codes for simulation function run
start_time <- Sys.time()

# Conditional power function (Pokok and Mehta modified code for sim)
source("placebo_model.R")
source("dgms.R")
source("ncs2.R")
source("model_fit.R")
source("warningCapture.R")
source("dropout_function.R")
source("cormatsim.R")


# simulation function
  sim <- function(delta1 = 0.3,
                  delta2 = 0.6,
                  delta3 = 4,
                  b0 = 3.4,
                  b1 = 11.3,
                  b2 = 3.7,
                  nhm1 = 3.4,
                  nhm2 = 6.5,
                  nhm3 = 8.9,
                  nhm4 = 9.9,
                  nhm5 = 10.1,
                  n_pbo = 40,
                  n_act = 80,
                  sd1 = 2,
                  sd2 = 2,
                  sd3 = 2,
                  sd4 = 2,
                  sd5 = 2,
                  # r12= 0.65,
                  # r13= 0.40,
                  # r14= 0.25,
                  # r15= 0.15,
                  # r23= 0.65,
                  # r24= 0.40,
                  # r25= 0.25,
                  # r35= 0.65,
                  # r45= 0.65,
                  cor = 0.65,
                  jitter_sd = 0,
                  missingPercentage = .1,
                  mtP1 = 0.05,
                  mtP2 = 0.1,
                  mtP3 = 0.2,
                  mtP4 = 0.3,
                  mtP5 = 0.4
                  ) {
  #Data generation

  sim_data <- dgm(delta1 = delta1, delta2 = delta2, delta3 = delta3, b0 = b0, b1 = b1, b2 = b2,
      n_pbo = n_pbo, n_act = n_act, sd1 = sd1, sd2 = sd2, sd3 = sd3, sd4 = sd4, sd5 = sd5,
      cor = cor, jitter_sd = jitter_sd, missingPercentage = missingPercentage,
      mtP1 = mtP1, mtP2 = mtP2, mtP3 = mtP3, mtP4 = mtP4, mtP5 = mtP4)
  #Model fitting

  sim_out <- bind_rows(
  do.call(rbind, lapply(sim_data, function(i) mmrm_fit(i))),
  do.call(rbind, lapply(sim_data, function(i) Func_ncs_int(i, last_visit = 24)))
  )
  return(sim_out)
}

sim.grid <- expand.grid(b0 = c(4),
                        b1 = seq(7, 25, by = 2), #seq(5, 25, by = 2),
                        b2 = seq(7, 9, by = 2),# seq(5, 9, by = 2),
                        delta1 = 0.3,
                        delta2 = c(0.3, 0.7),
                        delta3 = c(3, 4), # to expected reduction at the end of the time period
                        sd1 = 2,
                        sd2 = 3,
                        sd3 = 4,
                        sd4 = 5,
                        sd5 = 6,
                        cor = 0.65, #c(0, 0.65),
                        jitter_sd = .8, #c(0, 0.8),
                        missingPercentage = .1,
                        mtP1 = 0.05,
                        mtP2 = 0.05,
                        mtP3 = 0.2,
                        mtP4 = 0.3,
                        mtP5 = 0.4,
                        npbo = 40,
                        nact = 80
) %>% filter(!(b1 %in% c(17:25) & b2 ==5 )) %>%
  filter(!(b1 %in% c(19:25) & b2 ==7 )) %>%
  filter(mtP1 == 0.05, mtP2 == 0.05, mtP3 == 0.2, mtP4 == 0.3, mtP5 == 0.4)

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
                                    "sim",
                                    "autocorr.mat",
                                    "cormat",
                                    "get_phrase",
                                    "get_warning",
                                    "introduce_missing",
                                    "usCorrelation")
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
                         require(tmvtnorm)
                         }
  )

sim_FTD_output<- bind_rows( #future_map_dfr under purrr and furrr lib
  parallel::parApply(
  cl = cl,
  X = matrix(1: nrow(sim.grid)),#10),
  MARGIN = 1,
  FUN = function(x) {
    lapply(1:500, function(y) {
      bind_rows(
        c(
        bind_cols(
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
          cor = sim.grid$cor[x],
          jitter_sd = sim.grid$jitter_sd[x],
          missingPercentage = sim.grid$missingPercentage[x],
          mtP1 = sim.grid$mtP1[x],
          mtP2 = sim.grid$mtP2[x],
          mtP3 = sim.grid$mtP3[x],
          mtP4 = sim.grid$mtP4[x],
          mtP5 = sim.grid$mtP5[x],
          n_pbo = sim.grid$npbo[x],
          n_act = sim.grid$nact[x]
        ), do.call(rbind, (lapply(1:20, function(x) sim.grid[x,])))
        )
      ))
    })
  }
))

parallel::stopCluster(cl)
end_time <- Sys.time()
runtime <- (end_time - start_time) # in hours
save.image("./Outputs/tak594_FTD_SIM_01312025.RData")

# 5000 replication takes 22.40724 hours to run for 15 setup listed in the sim.grid

