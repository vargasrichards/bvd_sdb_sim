#...............................................................................
### ++++ EBOLA SAFE AND DIGNIFIED BURIALS IN DRC: SCENARIO EXPLORATION +++++ ###
#...............................................................................

#...............................................................................
## --- RUNNING SENSITIVITY ANALYSES: NUMBER OF(SEED) CASES, FUNERAL SIZE ---- ##
#...............................................................................

                              # Written by Francesco Checchi, LSHTM (May 2026)
                              # francesco.checchi@lshtm.ac.uk 



#...............................................................................                           
### Sensitivity analysis 1: Varying the number of seed cases
#...............................................................................

  #...................................       
  ## Select scenarios to highlight
    
    # Select scenarios
    scenarios$highlight <- F
    scenarios[which(
      scenarios$R0 %in% c(
        "R0 = 1.1 (0.4 due to burial)", 
        "R0 = 1.5 (0.4 due to burial)", 
        "R0 = 1.9 (0.8 due to burial)" 
      ) &
      scenarios$cov %in% c(0.4, 0.6, 0.8) & 
      scenarios$eff %in% c(0.6, 0.8, 1.0)), 
      "highlight"] <- T
    
    # Subset highlighted scenarios
    table(scenarios$highlight)
    scenarios_highlight <- subset(scenarios, highlight == T)
    
    
  #...................................       
  ## Run simulations
    
    # Initialise output dataframe
    n_sim <- 1:pars["n_sim"]
    out <- merge(scenarios_highlight, n_sim)
    colnames(out)[colnames(out) == "y"] <- "n_sim"  
    days = c(max(timeline$day - 21), max(timeline$day))
    out <- merge(out, days)
    colnames(out)[ncol(out)] <- "day"
    out <- out[order(out$id, out$n_sim, out$day), ]
    out$cases <- NA
    
    # Run highlight simulations
    pb <- txtProgressBar(min = 1, max = nrow(scenarios_highlight), style = 3)
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
          n_seeds = scenarios_highlight[i, "n_seeds"]),
        nsim = as.integer(length(n_sim)),
        format = "data.frame",
        include.data = F
      ))
      
      # collect results
      sim_i <- subset(sim_i, day %in% days)
      colnames(sim_i)[colnames(sim_i) == ".id"] <- "n_sim"
      sim_i <- sim_i[order(sim_i$n_sim, sim_i$day), ]
      out[which(out$id == scenarios_highlight[i, "id"]), 
        c("n_sim", "day", "cases")] <- sim_i[,c("n_sim","day", "cases")]
      setTxtProgressBar(pb, i)
    }    
    close(pb)

  #...................................       
  ## Visualise epidemic size
        
    # Compute average results 
    df <- subset(out, day == max(out$day))
    df$eff <- factor(percent(df$eff), levels = c("60%", "80%", "100%"))
    df$R0 <- paste0("R0 = ", df$R0_I + df$R0_D, " (", df$R0_D , 
      "\ndue to burial)")
    df$cov <- factor(paste0("coverage = ", percent(df$cov)),
      levels =paste0("coverage = ",c("40%", "60%", "80%")))
    df <- aggregate(cases ~ n_seeds + R0 + eff + cov, data = df,
      FUN = function(xx) {c(mean(xx), quantile(xx, c(0.5, 0.1, 0.9)))})
    df <- data.frame(df[, c("n_seeds", "R0", "eff", "cov")], unlist(df$cases))
    colnames(df) <- c("n_seeds", "R0", "eff", "cov", 
      "mean", "median", "quant10", "quant90")
    df$R0_fr <- df$R0
    df$R0_fr <- gsub("\ndue to burial", "attribuable \naux enterrements",
      df$R0_fr)
    df$R0_fr <- gsub("\\.", ",", df$R0_fr)
    df$cov_fr <- gsub("coverage", "couverture", df$cov)
    df$cov_fr <- factor(df$cov_fr,
      levels = paste0("couverture = ",c("40%", "60%", "80%")))
    write.csv(df, paste0(dir_path, "out/sens_n_seeds_cum_cases.csv"),
      row.names = F)
    
    # Visualise in English
    plot <- ggplot(df, aes(x = n_seeds, y = median, group = eff, 
      colour = eff, fill = eff)) +
      geom_point(alpha = 0.75, size = 2, shape = 1) +
      geom_line(linewidth = 1) +
      # geom_errorbar(stat = "identity", aes(ymin = quant10, ymax = quant90),
      #  width = 0.2) +
      scale_y_continuous("cumulative number of new cases") +
      scale_x_continuous("number of starting cases (epidemic maturity)") +
      scale_colour_manual("SDB effectiveness (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("SDB effectiveness (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0 ~ cov) +
      theme_bw()+
      theme(legend.position = "top", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_n_seeds_cum_cases_en.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)

    # Visualise in French
    plot <- ggplot(df, aes(x = n_seeds, y = median, group = eff, 
      colour = eff, fill = eff)) +
      geom_point(alpha = 0.75, size = 2, shape = 1) +
      geom_line(linewidth = 1) +
      # geom_errorbar(stat = "identity", aes(ymin = quant10, ymax = quant90),
      #  width = 0.2) +
      scale_y_continuous("nombre cumulatif de nouveaux cas") +
      scale_x_continuous("nombre initial de cases (maturité de l'épidémie)") +
      scale_colour_manual("complétude de l'EDS (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("complétude de l'EDS (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0_fr ~ cov_fr) +
      theme_bw()+
      theme(legend.position = "top", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_n_seeds_cum_cases_fr.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)
    
    
  #...................................       
  ## Compute and visualise extinction probability
      # extinction = zero new cases during last 21d (max incubation period)
        
    # Compute extinction probability
    df <- subset(out, day %in% c(max(out$day), max(out$day) - 21))
    df$eff <- factor(percent(df$eff), levels = c("60%", "80%", "100%"))
    df$R0 <- paste0("R0 = ", df$R0_I + df$R0_D, " (", df$R0_D , 
      "\ndue to burial)")
    df$cov <- factor(paste0("coverage = ", percent(df$cov)),
      levels =paste0("coverage = ",c("40%", "60%", "80%")))
    df <- aggregate(cases ~ n_seeds + n_sim + R0 + eff + cov, data = df, 
      FUN = diff)
    df$p_extinction <- ifelse(df$cases == 0, T, F)
    df <- aggregate(p_extinction ~ n_seeds + R0 + eff + cov, data = df, 
      FUN = mean)
    df$R0_fr <- df$R0
    df$R0_fr <- gsub("\ndue to burial", "attribuable \naux enterrements",
      df$R0_fr)
    df$R0_fr <- gsub("\\.", ",", df$R0_fr)
    df$cov_fr <- gsub("coverage", "couverture", df$cov)
    df$cov_fr <- factor(df$cov_fr,
      levels = paste0("couverture = ",c("40%", "60%", "80%")))
    write.csv(df, paste0(dir_path, "out/sens_n_seeds_p_extinction.csv"),
      row.names = F)
    
    # Visualise in English
    plot <- ggplot(df, aes(x = n_seeds, y = p_extinction, group = eff, 
      colour = eff, fill = eff)) +
      geom_point(alpha = 0.75, size = 2, shape = 1) +
      geom_line(linewidth = 1) +
      scale_y_continuous("probability of outbreak extinction by day 90", 
        labels = percent, limits = c(0, 1)) +
      scale_x_continuous("number of starting cases (epidemic maturity)") +
      scale_colour_manual("SDB effectiveness (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("SDB effectiveness (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0 ~ cov) +
      theme_bw()+
      theme(legend.position = "top", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_n_seeds_p_extinction_en.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)

    # Visualise in French
    plot <- ggplot(df, aes(x = n_seeds, y = p_extinction, group = eff, 
      colour = eff, fill = eff)) +
      geom_point(alpha = 0.75, size = 2, shape = 1) +
      geom_line(linewidth = 1) +
      scale_y_continuous("probabilité d'extinction de l'épidémie sur 90 jours", 
        labels = percent, limits = c(0, 1)) +
      scale_x_continuous("nombre initial de cases (maturité de l'épidémie)") +
      scale_colour_manual("complétude de l'EDS (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("complétude de l'EDS (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0_fr ~ cov_fr) +
      theme_bw()+
      theme(legend.position = "top", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_n_seeds_p_extinction_fr.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)

    
#...............................................................................                           
### Sensitivity analysis 2: Two extremes of SDB performance
#...............................................................................

  #...................................       
  ## Select scenarios to highlight
    
    # Select scenarios
    scenarios$highlight <- F
    scenarios[which(
      scenarios$n_seeds == pars["n_seeds"] &
      scenarios$R0 %in% c(
        "R0 = 1.3 (0.4 due to burial)", 
        "R0 = 1.5 (0.6 due to burial)", 
        "R0 = 1.9 (0.8 due to burial)" 
      ) &
      scenarios$cov %in% c(0, 1.0) & 
      scenarios$eff == 1.0), 
      "highlight"] <- T
    
    # Subset highlighted scenarios
    table(scenarios$highlight)
    scenarios_highlight <- subset(scenarios, highlight == T)
    
    
  #...................................       
  ## Run simulations
    
    # Initialise output dataframe
    n_sim <- 1:pars["n_sim"]
    out <- merge(scenarios_highlight, n_sim)
    colnames(out)[colnames(out) == "y"] <- "n_sim"  
    out <- merge(out, timeline)
    out <- out[order(out$id, out$n_sim, out$day), ]
    out$cases <- NA
    
    # Run highlight simulations
    pb <- txtProgressBar(min = 1, max = nrow(scenarios_highlight), style = 3)
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
          n_seeds = scenarios_highlight[i, "n_seeds"]),
        nsim = as.integer(length(n_sim)),
        format = "data.frame",
        include.data = F
      ))
      
      # collect results
      colnames(sim_i)[colnames(sim_i) == ".id"] <- "n_sim"
      sim_i <- sim_i[order(sim_i$n_sim, sim_i$day), ]
      out[which(out$id == scenarios_highlight[i, "id"]), 
        c("n_sim", "day", "cases")] <- sim_i[,c("n_sim","day", "cases")]
      setTxtProgressBar(pb, i)
    }    
    close(pb)

    
  #...................................       
  ## Visualise epidemic size
        
    # Compute average results 
    out$R0 <- paste0("R0 = ", out$R0_I + out$R0_D, " (", out$R0_D , 
      "\ndue to burial)")
    out$cov <- factor(paste0("coverage = ", percent(out$cov)),
      levels = paste0("coverage = ",c("0%", "100%")))
    df <- aggregate(cases ~ day + R0 + eff + cov, data = out,
      FUN = function(xx) {c(mean(xx), quantile(xx, c(0.5, 0.1, 0.9)))})
    df <- data.frame(df[, c("day", "R0", "eff", "cov")], unlist(df$cases))
    colnames(df) <- c("day", "R0", "eff", "cov", 
      "mean", "median", "quant10", "quant90")
    df$R0_fr <- df$R0
    df$R0_fr <- gsub("\ndue to burial", " attribuable \naux enterrements",
      df$R0_fr)
    df$R0_fr <- gsub("\\.", ",", df$R0_fr)
    df$cov_fr <- gsub("coverage", "couverture", df$cov)
    df$cov_fr <- factor(df$cov_fr,
      levels = paste0("couverture = ", c("0%", "100%")))    
    write.csv(df, paste0(dir_path, "out/sens_extreme_cum_cases.csv"),
      row.names = F)
    
    # Visualise in English
    plot <- ggplot(df, aes(x = day, y = median, colour = R0, fill = R0)) +
      geom_line(alpha = 0.75, linewidth = 2) +
      geom_ribbon(aes(ymin = quant10, ymax = quant90), alpha = 0.3) +
      scale_y_continuous("cumulative number of new cases") +
      scale_x_continuous("time since SDB implementation (days)") +
      scale_colour_manual("R0", values = palette_gen[c(13,9,5)]) +
      scale_fill_manual("R0", values = palette_gen[c(13,9,5)]) +
      facet_grid(R0 ~ cov) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_extreme_cum_cases_en.png"),
      units = "cm", dpi = "print", height = 15, width = 20 * hw)

    # Visualise in French
    plot <- ggplot(df, aes(x = day, y = median, colour = R0, fill = R0)) +
      geom_line(alpha = 0.75, linewidth = 2) +
      geom_ribbon(aes(ymin = quant10, ymax = quant90), alpha = 0.3) +
      scale_y_continuous("nombre cumulatif de nouveaux cas") +
      scale_x_continuous("nombre de jours depuis activation de l'EDS") +
      scale_colour_manual("R0", values = palette_gen[c(13,9,5)]) +
      scale_fill_manual("R0", values = palette_gen[c(13,9,5)]) +
      facet_grid(R0_fr ~ cov_fr) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_extreme_cum_cases_fr.png"),
      units = "cm", dpi = "print", height = 15, width = 20 * hw)
    
    
  #...................................       
  ## Compute and visualise extinction probability
      # extinction = zero new cases during last 21d (max incubation period)
        
    # Compute extinction probability
    df <- subset(out, day %in% c(max(out$day), max(out$day) - 21))
    df <- aggregate(cases ~ n_sim + R0 + cov, data = df, FUN = diff)
    df$p_extinction <- ifelse(df$cases == 0, T, F)
    df <- aggregate(p_extinction ~ R0 + cov, data = df, FUN = mean)
    df$R0_fr <- df$R0
    df$R0_fr <- gsub("\ndue to burial", " attribuable \naux enterrements",
      df$R0_fr)
    df$R0_fr <- gsub("\\.", ",", df$R0_fr)
    df$cov <- gsub("coverage = ", "", df$cov)
    df$cov <- factor(df$cov, levels = c("0%", "100%"))
    df$cov_fr <- df$cov
    write.csv(df, paste0(dir_path, "out/sens_extreme_p_extinction.csv"),
      row.names = F)
    
    # Visualise in English
    plot <- ggplot(df, aes(x = cov, y = p_extinction, 
      colour = R0, fill = R0)) +
      geom_bar(stat = "identity", alpha = 0.75) +
      scale_y_continuous("probability of outbreak extinction by day 90", 
        labels = percent, limits = c(0, 1)) +
      scale_x_discrete("SDB coverage (%)") +
      scale_colour_manual("R0", values = palette_gen[c(13,9,5)]) +
      scale_fill_manual("R0", values = palette_gen[c(13,9,5)]) +
      facet_grid(. ~ R0) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_extreme_p_extinction_en.png"),
      units = "cm", dpi = "print", height = 12, width = 20 * hw)

    # Visualise in French
    plot <- ggplot(df, aes(x = cov, y = p_extinction, 
      colour = R0, fill = R0)) +
      geom_bar(stat = "identity", alpha = 0.75) +
      scale_y_continuous("probabilité d'extinction de l'épidémie sur 90 jours", 
        labels = percent, limits = c(0, 1)) +
      scale_x_discrete("couverture de l'EDS (%)") +
      scale_colour_manual("R0", values = palette_gen[c(13,9,5)]) +
      scale_fill_manual("R0", values = palette_gen[c(13,9,5)]) +
      facet_grid(. ~ R0_fr) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/sens_extreme_p_extinction_fr.png"),
      units = "cm", dpi = "print", height = 12, width = 20 * hw)
       

#.........................................................................................
### ENDS
#.........................................................................................


