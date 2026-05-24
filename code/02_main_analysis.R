#...............................................................................
### ++++ EBOLA SAFE AND DIGNIFIED BURIALS IN DRC: SCENARIO EXPLORATION +++++ ###
#...............................................................................

#...............................................................................
## --- SCRIPT TO IMPLEMENT MAIN ANALYSIS SCENARIOS AND VISUALISE RESULTS ---- ##
#...............................................................................

                              # Written by Francesco Checchi, LSHTM (May 2026)
                              # francesco.checchi@lshtm.ac.uk 



#...............................................................................                           
### Running main analysis simulations
#...............................................................................

  #...................................       
  ## Select scenarios to highlight
    
    # Select scenarios
    scenarios$highlight <- F
    scenarios[which(
      scenarios$n_seeds == pars["n_seeds"] &
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
          n_seeds = as.integer(pars["n_seeds"])),
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
    out$eff <- factor(percent(out$eff), levels = c("60%", "80%", "100%"))
    out$R0 <- paste0("R0 = ", out$R0_I + out$R0_D, " (", out$R0_D , 
      "\ndue to burial)")
    out$cov <- factor(paste0("coverage = ", percent(out$cov)),
      levels =paste0("coverage = ",c("0%", "20%", "40%", "60%", "80%", "100%")))
    # saveRDS(out, paste0(dir_path, "out/sims_main_analysis.rds"))
    df <- aggregate(cases ~ day + R0 + eff + cov, data = out,
      FUN = function(xx) {c(mean(xx), quantile(xx, c(0.5, 0.1, 0.9)))})
    df <- data.frame(df[, c("day", "R0", "eff", "cov")], unlist(df$cases))
    colnames(df) <- c("day", "R0", "eff", "cov", 
      "mean", "median", "quant10", "quant90")
    df <- subset(df, day == max(df$day))
    df$R0_fr <- df$R0
    df$R0_fr <- gsub("\ndue to burial", " attribuable \naux enterrements",
      df$R0_fr)
    df$R0_fr <- gsub("\\.", ",", df$R0_fr)
    df$cov_fr <- gsub("coverage", "couverture", df$cov)
    df$cov_fr <- factor(df$cov_fr,
      levels = paste0("couverture = ",
        c("0%", "20%", "40%", "60%", "80%", "100%")))    
    write.csv(df, paste0(dir_path, "out/main_analysis_cum_cases.csv"),
      row.names = F)
    
    # Visualise in English
    plot <- ggplot(df, aes(x = eff, y = median, colour = eff, fill = eff)) +
      geom_point(alpha = 0.75, size = 4, shape = 22) +
      geom_errorbar(stat = "identity", aes(ymin = quant10, ymax = quant90),
        width = 0.2) +
      scale_y_continuous("cumulative number of new cases") +
      scale_x_discrete("SDB effectiveness (%)") +
      scale_colour_manual("SDB effectiveness (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("SDB effectiveness (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0 ~ cov) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/main_analysis_cum_cases_en.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)

    # Visualise in French
    plot <- ggplot(df, aes(x = eff, y = median, colour = eff, fill = eff)) +
      geom_point(alpha = 0.75, size = 4, shape = 22) +
      geom_errorbar(stat = "identity", aes(ymin = quant10, ymax = quant90),
        width = 0.2) +
      scale_y_continuous("nombre cumulatif de nouveaux cas") +
      scale_x_discrete("complétude de l'EDS (%)") +
      scale_colour_manual("complétude de l'EDS (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("complétude de l'EDS (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0_fr ~ cov_fr) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/main_analysis_cum_cases_fr.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)

    
  #...................................       
  ## Compute and visualise extinction probability
      # extinction = zero new cases during last 21d (max incubation period)
        
    # Compute extinction probability
    df <- subset(out, day %in% c(max(out$day), max(out$day) - 21))
    df <- aggregate(cases ~ n_sim + R0 + eff + cov, data = df, FUN = diff)
    df$p_extinction <- ifelse(df$cases == 0, T, F)
    df <- aggregate(p_extinction ~ R0 + eff + cov, data = df, FUN = mean)
    df$R0_fr <- df$R0
    df$R0_fr <- gsub("\ndue to burial", " attribuable \naux enterrements",
      df$R0_fr)
    df$R0_fr <- gsub("\\.", ",", df$R0_fr)
    df$cov_fr <- gsub("coverage", "couverture", df$cov)
    df$cov_fr <- factor(df$cov_fr,
      levels = paste0("couverture = ",
        c("0%", "20%", "40%", "60%", "80%", "100%")))    
    write.csv(df, paste0(dir_path, "out/main_analysis_p_extinction.csv"),
      row.names = F)

    # Visualise in English
    plot <- ggplot(df, aes(x = eff, y = p_extinction, 
      colour = eff, fill = eff)) +
      geom_bar(stat = "identity", alpha = 0.75) +
      scale_y_continuous("probability of outbreak extinction by day 90", 
        labels = percent, limits = c(0, 1)) +
      scale_x_discrete("SDB effectiveness (%)") +
      scale_colour_manual("SDB effectiveness (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("SDB effectiveness (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0 ~ cov) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/main_analysis_p_extinction_en.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)

    # Visualise in French
    plot <- ggplot(df, aes(x = eff, y = p_extinction, 
      colour = eff, fill = eff)) +
      geom_bar(stat = "identity", alpha = 0.75) +
      scale_y_continuous("probabilité d'extinction de l'épidémie sur 90 jours", 
        labels = percent, limits = c(0, 1)) +
      scale_x_discrete("complétude de l'EDS (%)") +
      scale_colour_manual("complétude de l'EDS (%)", 
        values = palette_gen[c(1,9,15)]) +
      scale_fill_manual("complétude de l'EDS (%)",
        values = palette_gen[c(1,9,15)]) +
      facet_grid(R0_fr ~ cov_fr) +
      theme_bw()+
      theme(legend.position = "none", panel.grid.major.x = element_blank())
    ggsave(paste0(dir_path, "out/main_analysis_p_extinction_fr.png"),
      units = "cm", dpi = "print", height = 20, width = 25 * hw)
    

#.........................................................................................
### ENDS
#.........................................................................................


