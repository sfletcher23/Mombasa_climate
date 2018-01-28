setwd("~/Desktop/Research/BMA/REA_Bayes")

source("./REA.Gibbs.r")
#source("./REA.GM.r")
#source("./REAMV.GM.r")
install.packages("foreach")
library(foreach)

#[1] 22 21
args(REA.Gibbs)
#function(X, X0, lambda0, Y, N= 15, save.every=50, burn.in= 100*
#save.every)
setwd("~/Desktop/Research/BMA/")
tmp = read.csv("lambda0.csv",header = FALSE)
lambda0 = as.matrix(tmp)

#  Here we could think about running this in parallel
year = c(1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090)
foreach (scen_ii = 1:2) %dopar%{
  foreach (ii = 1:10) %dopar%{
    yearX = year[ii]
    yearY = year[ii+1]
    str1 = sprintf("X_%d.csv",yearX)
    str2 = sprintf("X_%d.csv",yearY)
    str3 = sprintf("output%f.csv",yearY)
    tmp = read.csv(str1,header = FALSE)
    X = as.matrix(tmp)
    
    tmp = read.csv(str2,header = FALSE)
    Y = as.matrix(tmp)
    
    tmp = read.csv("X0TU.csv",header = FALSE)
    X0T = as.matrix(tmp[c(ii),c(scen_ii)])
    
    REA.Gibbs(X[c(1),],X0T[c(1)],lambda0[c(1)],Y[c(1),],N=1000)->rg0
    
    mu0str = sprintf("muUT_output%d_scen%d.csv",yearY,scen_ii)
    nu0str = sprintf("nuUT_output%d_scen%d.csv",yearY,scen_ii)
    write.table(rg0$mu, file = mu0str,row.names=FALSE, col.names=FALSE)
    write.table(rg0$nu, file = nu0str,row.names=FALSE, col.names=FALSE)
    
    #tmp = read.csv("X0PU.csv",header = FALSE)
    #X0P = tmp[c(ii),c(scen_ii)]
    #REA.Gibbs(X[c(2),],X0P[c(1)],lambda0[c(2)],Y[c(2),],N=1000)->rg0

    #mu0str = sprintf("muUP_output%d_scen%d.csv",yearY,scen_ii)
    #nu0str = sprintf("nuUP_output%d_scen%d.csv",yearY,scen_ii)
    #write.table(rg0$mu, file = mu0str,row.names=FALSE, col.names=FALSE)
    #write.table(rg0$nu, file = nu0str,row.names=FALSE, col.names=FALSE)
  }
}

