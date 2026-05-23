#...............................................................................
### ++++ EBOLA SAFE AND DIGNIFIED BURIALS IN DRC: SCENARIO EXPLORATION +++++ ###
#...............................................................................

#...............................................................................
## ------- 'CONTROL' R SCRIPT TO SOURCE INPUTS AND CALL OTHER SCRIPTS ------- ##
#...............................................................................

                              # Written by Francesco Checchi, LSHTM (May 2026)
                              # francesco.checchi@lshtm.ac.uk 


#...............................................................................
### Preparatory steps
#...............................................................................

  #...................................      
  ## Install or load required R packages

    # pacman and epidemics
    if (!"pacman" %in% rownames(installed.packages())) {
      install.packages("pacman")}
    # if (!"epidemics" %in% rownames(installed.packages())) {
    #   install.packages('epidemics', 
    #     repos = c('https://epiverse-trace.r-universe.dev',
    #       'https://cloud.r-project.org'))}
    # if (!"epiparameterDB" %in% rownames(installed.packages())) {
    #   install.packages('epiparameterDB', 
    #     repos = c('https://epiverse-trace.r-universe.dev', 
    #       'https://cloud.r-project.org'))}
    # if (!"epiparameter" %in% rownames(installed.packages())) {
    #   install.packages('epiparameter', 
    #     repos = c('https://epiverse-trace.r-universe.dev', 
    #       'https://cloud.r-project.org'))}
    
    
  # Install or load packages from CRAN
    pacman::p_load(
      ggplot2,       # Visualise data
      ggpubr,        # Arrange multiple plots into a single plot
      pomp,          # Implement SEIR model in C++ 
      readxl,        # Read Excel files
      scales,        # Scale and format data for visualisation
      #socialmixr,   # Social mixing matrices
      tidyverse,     # Tidyverse suite of packages
      viridis)       # Colour-blind palette


  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
    # Set font for Windows or Mac
    suppressWarnings(windowsFonts(Arial = windowsFont("Arial")))
    suppressWarnings(par(family = "Arial"))

    # Set working directory to where this file is stored
    dir_path <- paste(dirname(rstudioapi::getActiveDocumentContext()$path  )
      , "/", sep = "")
    setwd(dir_path)
    print( getwd() )
    dir_path <- gsub("/code", "", dir_path)
    suppressWarnings(dir.create(paste0(dir_path, "out")))
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
      # general palette
      palette_gen <- viridis(16)
      show_col(palette_gen)
       
    # height to width ratio for figures (A4 format: should be 1.414, but
      # leave some space for captions)
    hw <- 1.3
    
    
#...............................................................................
### Running scenario simulations
#...............................................................................
    
source(paste0(dir_code, "/01_model_scenarios.R"))    
    
    
#...............................................................................
### ENDS
#...............................................................................


