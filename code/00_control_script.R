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

  # NOTE: need to install RTools for current version of R
    # needed to compile C++ code

  # pacman if not already installed
  if (!"pacman" %in% rownames(installed.packages())) {
    install.packages("pacman", repos = "https://cloud.r-project.org")
    }

  # Install or load packages from CRAN
  pacman::p_load(
    ggplot2,       # Visualise data
    ggpubr,        # Arrange multiple plots into a single plot
    pomp,          # Implement SEIR model in C++ 
    readxl,        # Read Excel files
    scales,        # Scale and format data for visualisation
    tidyverse,     # Tidyverse suite of packages
    viridis,       # Colour-blind palette
    here)          # Construct paths relative to project root directory


#...................................      
## Starting setup

  # Clean up from previous code / runs
  rm(list=ls(all=TRUE))


# Set font for Windows or Mac
if (.Platform$OS.type == "windows") {
  suppressWarnings(windowsFonts(Arial = windowsFont("Arial")))
} else {
  suppressWarnings(par(family = "Arial"))
}

# Set working directory to project root
setwd(here())
print(getwd())

# Derive project root (parent of /code)
dir_path <- here()
suppressWarnings(dir.create(paste0(dir_path, "/out")))

# Initialise random numbers
set.seed(123)

# Colour-blind palette for graphing
palette_gen <- viridis(16)
if (interactive()) show_col(palette_gen)

# Height to width ratio for figures
hw <- 1.3


#...............................................................................
### Running scenario simulations
#...............................................................................
#...................................      
## Set up model and define scenarios
source(here("code", "01_prepare_model.R"))

#...................................      
## Main analysis
source(here("code", "02_main_analysis.R"))

#...................................      
## Sensitivity analyses
source(here("code", "03_sens_analyses.R"))
  
#...............................................................................
### ENDS
#...............................................................................


