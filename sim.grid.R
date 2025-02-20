# Simulation grid setup
#YS 02/12/2025

sim.dfP1 <- expand.grid(b0 = c(4, 5),
            b1 = seq(7, 25, by = 2), #seq(5, 25, by = 2),
            b2 = seq(7, 9, by = 2)) %>% # seq(5, 9, by = 2)
            slice(c(1:15, 26:30)) # controlled so that the highest level is more than 17


sim.dfP2 <- expand.grid(
            delta1 = c(0.1, 0.2, 0.3, .4),
            delta2 = c(0.3, 0.5, 0.7),
            delta3 = c(2, 3, 4, 5), # to expected reduction at the end of the time period
            sd1 = 2,
            sd2 = 3,
            sd3 = 4,
            sd4 = 5,
            sd5 = 5,
            cor = 0.65, #c(0, 0.65),
            jitter_sd = .8, #c(0, 0.8),
            missingPercentage = .1,
            mtP1 = 0.05,
            mtP2 = 0.05,
            mtP3 = 0.2,
            mtP4 = 0.3,
            mtP5 = 0.4,
            npbo = 40,
            nact = 80)


sim.grid_t <- do.call(rbind, lapply(1:nrow(sim.dfP1), function(x)

 do.call(rbind, (lapply(1:nrow(sim.dfP2), function(i) sim.dfP1[x,]))) %>% bind_cols(sim.dfP2)
)
)
# End
