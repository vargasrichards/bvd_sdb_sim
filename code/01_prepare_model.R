#...............................................................................
### ++++ EBOLA SAFE AND DIGNIFIED BURIALS IN DRC: SCENARIO EXPLORATION +++++ ###
#...............................................................................

#...............................................................................
## ----- SCRIPT TO PREPARE STOCHASTIC MODEL AND DEFINE MODEL SCENARIOS ------ ##
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
      n_seeds = seq(1, 19, 2),
      R0_I = seq(0.50, 1.50, 0.10),
      R0_D = seq(0.40, 0.80, 0.10),
      cov = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), 
      eff = unique(effect_sdb$eff)
    )
    scenarios <- merge(scenarios, effect_sdb, by = "eff", all.x = T)    
    scenarios$id <- 1:nrow(scenarios)
    scenarios$R0 <- paste0("R0 = ", scenarios$R0_I + scenarios$R0_D, " (", 
      scenarios$R0_D, " due to burial)")
    
    
#...............................................................................                           
### Setting up stochastic model
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
      cases += dN_EI;
    ")

    # Specify initial conditions    
    seird_init <- Csnippet("
      S = N - n_seeds;
      E = round(n_seeds * 
        incubation_period / (incubation_period + infectious_period));
      I = n_seeds - E;
      R = 0;
      D = 0;
      cases = n_seeds;
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
      statenames = c("S", "E", "I", "R", "D", "cases")
      # accumvars = "cases" # use this to track daily incidence
    )


#.........................................................................................
### ENDS
#.........................................................................................


