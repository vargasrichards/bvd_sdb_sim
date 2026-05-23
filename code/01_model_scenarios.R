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
  ## Demographic and social mixing parameters
    # Read health zone population data and compute mean in affected provinces
    pop <- as.data.frame(read_xlsx(paste0(dir_path, 
      "in/drc-hpc-projection-population-2024.xlsx")))  
    x <- pop[which(pop$Province %in% c("Ituri", "Nord-Kivu", "Sud-Kivu")),
      "Population 2024"]
    pop_hz_mean <- round(mean(x), digits = 0)
    pop_hz_sd <- sd(x)
  
    
  #...................................       
  ## Transmission dynamic parameters
    
    # Basic reproduction number
    r0 <- seq(1.0, 2.0, 0.2)
    
    # Infectious period
    infectious_period <- 10 + 1 # symptom onset to death/recovery + 2 days of 
      # infectiousness after death (if CFR ~ 50%, add 1 day on average)
    
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
    
    # Safe and dignified burial (SDB) intervention scenarios    
    scenarios <- expand.grid(coverage = seq(0, 1, 0.1), 
      prop_success = unique(effect_sdb$prop_success))
    scenarios <- merge(scenarios, effect_sdb, by = "prop_success", all.x = T)    
    
    
#...............................................................................                           
### Setting up model simulations
#...............................................................................

  #...................................       
  ## Set up model (no age structure)
    
    # SEIRD time step
    seird_step <- Csnippet("
      double dN_SE = rbinom(S, 1 - exp(-(R0_I/infectious_period) * I/N)) + 
        rbinom(S, 1 - exp(-((R0_D/burial_period) - ((1 - cov) * eff)) * D/N));
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
      Icum += dN_EI;
    ")

    # Specify initial conditions    
    seird_init <- Csnippet("
      S = N - n_seeds;
      E = round(n_seeds * 
        incubation_period / (incubation_period + infectious_period));
      I = n_seeds - E;
      R = 0;
      D = 0;
      Icum = 0;
    ")

    # Specify timeline
    out <- data.frame(day = 1:60)
    
    # Fold everything into pomp object
    evd_seird <- pomp(data = out, times = "day", rprocess = euler(seird_step,delta.t = 1), 
      rinit = seird_init, t0 = 1,
      paramnames = c("N", "R0_I", "R0_D", "incubation_period", 
        "infectious_period", "burial_period", "cfr", "cov", "eff", "n_seeds"),
      statenames = c("S", "E", "I", "R", "D", "Icum"),
      accumvars = "Icum")

    
  #...................................       
  ## Run simulations
    
    # Run simulations
    sim1 <- simulate(evd_seird, 
      params = c(N = pop_hz_mean, R0_I = 1.5, R0_D = 0.5, 
        incubation_period = 3, infectious_period = 10, burial_period = 1,
        cfr = 0.5, cov = 1, eff = 0.3, n_seeds = 10),
      nsim = 5,
      format = "data.frame",
      include.data = F
    )
    
    # Visualise results
    ggplot(sim1, aes(x = day, y = Icum, group = .id)) +
      geom_line()
    
    
    

    


#.........................................................................................
### ENDS
#.........................................................................................


