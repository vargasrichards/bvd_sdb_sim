#...............................................................................
### ++++ EBOLA SAFE AND DIGNIFIED BURIALS IN DRC: SCENARIO EXPLORATION +++++ ###
#...............................................................................

#...............................................................................
## ----- SCRIPT TO IMPLEMENT DIFFERENT SCENARIOS AND VISUALISE RESULTS ------ ##
#...............................................................................

                              # Written by Francesco Checchi, LSHTM (May 2026)
                              # francesco.checchi@lshtm.ac.uk 


#...............................................................................                           
### Specifying parameters / setting up scenarios
#...............................................................................

  #...................................       
  ## Read in simulation parameters

    # Read parameter file and create parameter vector
    pars_df <- read_xlsx(paste0(dir_path, "in/evd_sdb_sim_parameters.xlsx"))
    pars <- pars_df$value
    names(pars) <- pars_df$parameter
    
    # Effect of safe and dignified burial (SDB) on reproduction number, based on
      # two methods (Hirano-Imbens [hi], Propensity Weights [pw]) 
      # in Checchi et al. (2025)
    effect_hi <- read.csv(paste0(
      dir_path, "in/out_dose_resp_rn_p_success_hi.csv"))
    effect_pw <- read.csv(paste0(
      dir_path, "in/out_dose_resp_rn_p_success_pw.csv"))
    effect_sdb <- (effect_hi + effect_pw) / 2 # simple average of two methods
    x <- c("mean", "low", "high")
    colnames(effect_sdb) <- c("prop_success", x)
    for (i in 2:nrow(effect_sdb)) {effect_sdb[i, x] <- 
      effect_sdb[i, x] - effect_sdb[1, x]}
    effect_sdb[1, x] <- c(0, 0, 0)
    colnames(effect_sdb) <- c("eff", "mean", "low", "high")
    effect_sdb <- abs(effect_sdb)
    
  #...................................       
  ## Set up scenarios
    
    # Safe and dignified burial (SDB) intervention scenarios    
    scenarios <- expand.grid(
      R0_I = seq(0.50, 1.50, 0.10),
      R0_D = seq(0.40, 0.80, 0.10),
      cov = seq(0, 1, 0.2), 
      eff = unique(effect_sdb$eff)
    )
    scenarios <- merge(scenarios, effect_sdb, by = "eff", all.x = T)    
    scenarios$id <- 1:nrow(scenarios)
    scenarios$R0 <- paste0("R0 = ", scenarios$R0_I + scenarios$R0_D, " (", 
      scenarios$R0_D, " due to burial)")
    
    # Select and name scenarios to highlight
    scenarios$highlight <- F
    scenarios[which(
      scenarios$R0 %in% c(
        "R0 = 1.1 (0.4 due to burial)", 
        "R0 = 1.3 (0.6 due to burial)", 
        "R0 = 1.5 (0.4 due to burial)", 
        "R0 = 1.7 (0.6 due to burial)", 
        "R0 = 1.9 (0.8 due to burial)" 
      ) &
      scenarios$cov %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0) & 
      scenarios$eff %in% c(0.6, 0.8, 1.0)), 
      "highlight"] <- T
    table(scenarios$highlight)
    scenarios_highlight <- subset(scenarios, highlight == T)
    
    
#...............................................................................                           
### Setting up model simulations
#...............................................................................

  #...................................       
  ## Set up stochastic model
    
    # SEIRD time step
    seird_step <- Csnippet("
      double dN_SE = rbinom(S, 1 - exp(-(R0_I/infectious_period) * I/N)) + 
        rbinom(S, 1 - exp(-((R0_D - (cov * eff))/burial_period) * D/N));
      double dN_EI = rbinom(E, 1 - exp(-(1/incubation_period)));
      double dN_IF = rbinom(I, 1 - exp(-(1/infectious_period)));
      double dN_FD = rbinom(dN_IF, cfr);
      double dN_FR = dN_IF - dN_FD;
      double dN_DR = rbinom(D, 1 - exp(-(1/burial_period)));

      S -= dN_SE;
      E += dN_SE - dN_EI;
      I += dN_EI - dN_IF;
      D += dN_FD - dN_DR;
      R += dN_FR + dN_DR;
      new_cases += dN_EI;
    ")

    # Specify initial conditions    
    seird_init <- Csnippet("
      S = N - n_seeds;
      E = round(n_seeds * 
        incubation_period / (incubation_period + infectious_period));
      I = n_seeds - E;
      R = 0;
      D = 0;
      new_cases = 0;
    ")

    # Specify timeline
    timeline <- data.frame(day = 1:pars["n_days"])
    
    # Fold everything into pomp object
    evd_seird <- pomp(
      data = timeline, 
      times = "day", 
      rprocess = euler(seird_step, delta.t = 1), 
      rinit = seird_init,
      t0 = 1,
      paramnames = c("N", "R0_I", "R0_D", "incubation_period", 
        "infectious_period", "burial_period", "cfr", "cov", "eff", "n_seeds"),
      statenames = c("S", "E", "I", "R", "D", "new_cases")
      # accumvars = "new_cases"
    )


#...............................................................................                           
### Running highlight simulations
#...............................................................................
        
  #...................................       
  ## Run simulations
    
    # Initialise output dataframe
    n_sim <- 1:pars["n_sim"]
    out <- merge(scenarios_highlight, n_sim)
    colnames(out)[colnames(out) == "y"] <- "n_sim"  
    out <- merge(out, timeline)
    out <- out[order(out$id, out$n_sim, out$day), ]
    out$new_cases <- NA
    
    # Run highlight simulations
    for (i in 1:nrow(scenarios_highlight)) {
      # run simulation
      sim_i <- suppressWarnings(simulate(evd_seird, 
        params = c(
          N = as.integer(pars["pop"]), 
          R0_I = scenarios_highlight[i, "R0_I"], 
          R0_D = scenarios_highlight[i, "R0_D"], 
          incubation_period = as.numeric(pars["incubation_period"]), 
          infectious_period = as.numeric(pars["infectious_period"]),
          burial_period = as.numeric(pars["burial_period"]),
          cfr = as.numeric( pars["cfr"]),
          cov = scenarios_highlight[i, "cov"], 
          eff = effect_sdb[which(
            effect_sdb$eff == scenarios_highlight[i, "eff"]), "mean"],
          n_seeds = as.integer(pars["n_seeds"])),
        nsim = as.integer(pars["n_sim"]),
        format = "data.frame",
        include.data = F
      ))
      
      # collect results
      colnames(sim_i)[colnames(sim_i) == ".id"] <- "n_sim"
      sim_i <- sim_i[order(sim_i$n_sim, sim_i$day), ]
      out[which(out$id == scenarios_highlight[i, "id"]), 
        c("n_sim", "day", "new_cases")] <- sim_i[,c("n_sim","day", "new_cases")]
    }    


  #...................................       
  ## Visualise final epidemic size
        
    # Compute average results 
    out$eff <- factor(percent(out$eff), levels = c("60%", "80%", "100%"))
    out$R0 <- paste0("R0 = ", out$R0_I + out$R0_D, " (", out$R0_D , 
      " due to burial)")
    out$cov <- factor(paste0("coverage = ", percent(out$cov)),
      levels =paste0("coverage = ",c("0%", "20%", "40%", "60%", "80%", "100%")))
    df <- aggregate(new_cases ~ day + R0 + eff + cov, data = out,
      FUN = function(xx) {c(mean(xx), quantile(xx, c(0.5, 0.1, 0.9)))})
    df <- data.frame(df[, c("day", "R0", "eff", "cov")], unlist(df$new_cases))
    colnames(df) <- c("day", "R0", "eff", "cov", 
      "mean", "median", "quant10", "quant90")
    df <- subset(df, day == max(df$day))
    
    # Visualise
    plot <- ggplot(df, aes(x = eff, y = median, colour = eff, fill = eff)) +
      geom_point(alpha = 0.75, size = 4, shape = 22) +
      geom_errorbar(stat = "identity", aes(ymin = quant10, ymax = quant90),
        width = 0.2) +
      scale_y_continuous("cumulative number of new cases") +
      scale_x_discrete("SDB effectiveness (%)") +
      scale_colour_viridis_d("SDB effectiveness (%)") +
      scale_fill_viridis_d("SDB effectiveness (%)") +
      facet_grid(cov ~ R0) +
      theme_pubclean()+
      theme(legend.position = "none")
    ggsave(paste0(dir_path, "out/scen_highlight_final_size.png"),
      units = "cm", dpi = "print", height = 20, width = 20 * hw)

  #...................................       
  ## Compute and visualise extinction probability
      # extinction = zero new cases during last 14d (edge EVD incubation period)
        
    # Compute extinction probability
    df <- subset(out, day %in% c(max(out$day), max(out$day) - 14))
    df <- aggregate(new_cases ~ n_sim + R0 + eff + cov, data = df, FUN = diff)
    df$extinct <- ifelse(df$new_cases == 0, T, F)
    df <- aggregate(extinct ~ R0 + eff + cov, data = df, FUN = mean)
    
    # Visualise
    plot <- ggplot(df, aes(x = eff, y = extinct, colour = eff, fill = eff)) +
      geom_bar(stat = "identity", alpha = 0.75) +
      scale_y_continuous("probability of outbreak extinction", labels = percent) +
      scale_x_discrete("SDB effectiveness (%)") +
      scale_colour_viridis_d("SDB effectiveness (%)") +
      scale_fill_viridis_d("SDB effectiveness (%)") +
      facet_grid(cov ~ R0) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/scen_highlight_p_extinction.png"),
      units = "cm", dpi = "print", height = 20, width = 20 * hw)

    

#.........................................................................................
### ENDS
#.........................................................................................


