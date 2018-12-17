# This file is a script to run a batch of Bayeisan analysis for precipitaiton.
# It calls REA.Gibbs, which uses MCMC methods to estimate the precipitation distribution for the next time period. 
# We repeat this for each virutal observation of temperature and precipitation in each time period


######################################

# Set up to run on a supercomputing cluster that uses a SLURM queueing system

######################################

# Get environment var from Slurm
JOBID = Sys.getenv("SLURM_JOB_ID")

# Check if on Slurm. If not, set working directory
if (JOBID != "") {
  cat("Slurm Job ID Found")
} else {
  setwd("~/Documents/MATLAB/Mombasa_Climate/BMA_code")
  cat("Slurm Job ID Not Found")
}

source("./REA.Gibbs.r")
install.packages("foreach", repos = "http://cran.us.r-project.org")
library(foreach)

# Need to use this package for the dopar commands to work
install.packages("doParallel", repos = "http://cran.us.r-project.org")
library(doParallel)

# Set up the parallelization. If running on slurm, use the number of cores available in the job. If running on desktop, use 2. 
if (JOBID != "") {
  registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
} else {
  registerDoParallel(cores=2)
}

# Get date for file save name
currentDate = Sys.Date()

######################################

# Run Bayesian analysis

######################################


# Load lambda prior information from preprocessing
args(REA.Gibbs)
tmp = read.csv("Input/lambda0.csv",header = FALSE)
lambda0 = as.matrix(tmp)


#  Loop over years and precip, run Bayesian analysis for each
year = c(1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090)
foreach (scen_ii = 1:3) %dopar%{
  foreach (ii = 10) %dopar%{
    yearX = year[ii]    # Historical time period
    yearY = year[ii+1]  # Future time period
    str1 = sprintf("Input/X_%d.csv",yearX)
    str2 = sprintf("Input/X_%d.csv",yearY)
    tmp = read.csv(str1,header = FALSE)
    X = as.matrix(tmp)
    
    tmp = read.csv(str2,header = FALSE)
    Y = as.matrix(tmp)
    
    tmp = read.csv("Input/X0PU.csv",header = FALSE)
    X0P = tmp[c(ii),c(scen_ii)]
    REA.Gibbs(X[c(2),],X0P[c(1)],lambda0[c(2)],Y[c(2),],N=1000)->rg0
    
    mu0str = paste(sprintf("Output/muUP_%d_scen%d",yearY,scen_ii),"job", JOBID, currentDate, ".csv", sep="_")
    nu0str = paste(sprintf("Output/nuUP_%d_scen%d",yearY,scen_ii),"job", JOBID, currentDate, ".csv", sep="_")
    write.table(rg0$mu, file = mu0str,row.names=FALSE, col.names=FALSE)
    write.table(rg0$nu, file = nu0str,row.names=FALSE, col.names=FALSE)
  }
}

