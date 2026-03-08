# This empirical example estimates the CCAPM model of Hansen and Singleton (1982, ECMA) to illustrate the GMM estimation of 
# nonlinear models

# Data sources: Saint Louis Fed and Yahoo Finance
# URL: https://www.stlouisfed.org/, http://finance.yahoo.com/

#   Description of the data: monthly data from 1959:M01 to 2015:M01 (675 obs) of the following variables: 
# - PCEND: Personal Consumption Expenditures: Nondurable Goods (billions of dollars) from FED
# - PPCEND: Personal Consumption expenditures: Nondurable goods (chain-type price index), DNDGRG3M086SBEA from FED
# - CNP16OV: Civilian Noninstitutional Population (thousands of persons) from FED
# - GS1: annualized 1-Year Treasury Constant Maturity Rate from FED
# - SP500: S&P 500 index at closing price of the last day of the month, ^GSPC from Yahoo Finance
####################

##############################################################################
# Processing Data
##############################################################################
library(AER)
library(sandwich)
library(gmm)
library(quantmod)
library(xtable) 

# WD
setwd("/Users/martinhiti/Desktop/courses/14382/14.382_pset2")

# Reading the data
 raw.data = as.data.frame(read.csv("data/ccapm-long.csv", header=T ))
 attach(raw.data)

# Preparing data
 rCpc        = PCEND/(PPCEND * CNP16OV)
 a.inflation =  PPCEND/Lag(PPCEND, k=12) - 1

 Rb          = (1 + GS1/100 - a.inflation)^(1/12)   #total monthly return to bonds (deflated)
 Rm          = (SP500/Lag(SP500))*(Lag(PPCEND)/PPCEND)  #total monthly return to stocks (deflated)
 Rc          = rCpc/Lag(rCpc)  #total monthly return to per-capita consumption (deflated)
 Ret=        na.omit(cbind(Rc, Rm, Rb))
 colnames(Ret) = c("Rc", "Rm", "Rb")

ts.plot(Ret[500:600,], main="consumption, stock, and bond total returns", col=c(1,4,3))
 detach("raw.data")

pdf(file="output/ccapm_series_plot.PDF", width=6, height=3)
ts.plot(Ret[500:600,], main="consumption, stock, and bond total returns", col=c(1,4,3))
dev.off()

# Write out and Read-in Processed data

write.csv(Ret,file="data/ccapm-ready-to-use.csv",row.names=FALSE)

# We use the same cleaned file as the lecture notes
## Not sure why we created the ready to use file and then don't use it tho
Ret= read.csv(file="data/ccapm-long_mert.csv", header=T)
attach(Ret)
str(Ret)
summary(Ret)


##############################################################################
# Score Functions 
##############################################################################
# the score function here returns the n by m  matrix of scores to be used in gmm package; 
# i.e., the return the m- vector of scores g(X_i, \theta) evaluated at each data point X_i
# for each observation i=1,...,N 

# the score function for estimationg CRRA (power utility) preferences with 2 financials assets

g2.x.theta = function(theta, x)
{
  # input par vector theta
  # input data matrix
  y = x[ ,1:3]    #first three columnds of x are y
  z =  x[, -c(1:3)]  # the rest are instruments
    rho.1 = (theta[1] * y[ ,2] * y[ ,1]^(-theta[2]) - 1)  # structural res 1
  rho.2 =  (theta[1] * y[ ,3] * y[ ,1]^(-theta[2]) - 1)    # structural res 2
  score.1 = rho.1*z #score 1
  score.2 = rho.2*z #score 2
  return (cbind( score.1, score.2))    # output m by n matrix of scores
}

##############################################################################
# Replicate Lecture as Sanity Check 
##############################################################################

# Estimation with  Instruments = 1 lag of consumption and stock returns
y     = na.omit(cbind(Rc, Rb, Rm))
nlags = 1
z    = cbind(1, Lag(y[ ,1], c(1:nlags)), Lag(y[ ,2], c(1:nlags)), Lag(y[ ,3], c(1:nlags)))
x    = na.omit(cbind(y,z))

# GMM fit
gmm2.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1), method = "BFGS" ,  type="iter", vcov="iid"))
print(gmm2.fit)

# CUE fit
gmm2.cue.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1),  type="cue",  vcov="iid"))
print(gmm2.cue.fit)

# Estimation with Instrument = 1 lag of Consumption and Market Returns and Squares and Interactions 
y     = na.omit(cbind(Rc, Rb, Rm))
nlags = 1
z    = cbind(Lag(y[ ,1], c(1:nlags)), Lag(y[ ,2], c(1:nlags)), Lag(y[ ,3], c(1:nlags)))
z=  cbind( 1, z[,1],z[,2],z[,3], z[,1]^2, z[,2]^2, z[,3]^2, z[,1]*z[,2], z[,1]*z[,3], z[,2]*z[,3])
x   = na.omit(cbind(y,z))


# GMM fit
gmm2.zsq.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1),  vcov="iid", type="iter"))
print(gmm2.zsq.fit)

# CUE fit
gmm2.zsq.cue.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1), type="cue", vcov="iid"))
print(gmm2.zsq.cue.fit)

# Printing the results
tableJ= matrix(0, ncol=4, nrow=2)
tableJ[,1]= gmm2.fit$stest[[2]]
tableJ[,2]= gmm2.cue.fit$stest[[2]]
tableJ[,3]= gmm2.zsq.fit$stest[[2]]
tableJ[,4]= gmm2.zsq.cue.fit$stest[[2]]
colnames(tableJ)= c("GMM-1", "CUE-1", "GMM-2", "CUE-2")
rownames(tableJ)= c("J-statistic", "p-value")

tableE= matrix(0, ncol=4, nrow=4);
tableE[c(1:2),1]= gmm2.fit$coef[1,c(1:2)]
tableE[c(3:4),1]= gmm2.fit$coef[2,c(1:2)]
tableE[c(1:2),2]= gmm2.cue.fit$coef[1,c(1:2)]
tableE[c(3:4),2]= gmm2.cue.fit$coef[2,c(1:2)]
tableE[c(1:2),3]= gmm2.zsq.fit$coef[1,c(1:2)]
tableE[c(3:4),3]= gmm2.zsq.fit$coef[2,c(1:2)]
tableE[c(1:2),4]= gmm2.zsq.cue.fit$coef[1,c(1:2)]
tableE[c(3:4),4]= gmm2.zsq.cue.fit$coef[2,c(1:2)]

colnames(tableE)= c("GMM-1", "CUE-1", "GMM-2", "CUE-2")
rownames(tableE)= c("estimated beta", "std error beta",  "estimated alpha", "std error for alpha")

#output tables in latex format

print(xtable(tableE, digits=4, align=c(rep("c", 5))), file="output/ccapm_replication_coefs.tex")
print(xtable(tableJ, digits=3, align=c(rep("c", 5))), file="output/ccapm_replication_J.tex")

##############################################################################
# Homework
# (a) Use only an intercept as an instrument
# (b) Use an intercept and on lag of each return (stock and bond) as instrument
##############################################################################

# (a) Estimation with instruments = intercept
y    = na.omit(cbind(Rc, Rb, Rm)) # outcomes
z    = 1
x    = na.omit(cbind(y,z))

# GMM (a) fit
gmm.a.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1), method = "BFGS" ,  type="iter", vcov="iid"))
print(gmm.a.fit)

# CUE (a) fit
gmm.a.cue.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1),  type="cue",  vcov="iid"))
print(gmm.a.cue.fit)


# (b) Estimation with instrument = 1 lag of stock return and 1 lag of bond returns
y     = na.omit(cbind(Rc, Rb, Rm))
nlags = 1
z    = cbind(1, Lag(y[ ,2], c(1:nlags)), Lag(y[ ,3], c(1:nlags)))
x   = na.omit(cbind(y,z))

# GMM (b) fit
gmm.b.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1),  vcov="iid", type="iter"))
print(gmm.b.fit)

# CUE fit
gmm.b.cue.fit  = summary(gmm(g2.x.theta, x, t0 = c(.99,1), type="cue", vcov="iid"))
print(gmm.b.cue.fit)

# Printing the results

# Only show J-test for (b). With (a) we have m=d so no overidentification
tableJ= matrix(0, ncol=2, nrow=2)
tableJ[,1]= gmm.b.fit$stest[[2]]
tableJ[,2]= gmm.b.cue.fit$stest[[2]]
colnames(tableJ)= c("GMM-(b)", "CUE-(b)")
rownames(tableJ)= c("J-statistic", "p-value")

# Coefficients + SEs for (a) and (b)
tableE= matrix(0, ncol=4, nrow=4);
tableE[c(1:2),1]= gmm.a.fit$coef[1,c(1:2)]
tableE[c(3:4),1]= gmm.a.fit$coef[2,c(1:2)]
tableE[c(1:2),2]= gmm.a.cue.fit$coef[1,c(1:2)]
tableE[c(3:4),2]= gmm.a.cue.fit$coef[2,c(1:2)]
tableE[c(1:2),3]= gmm.b.fit$coef[1,c(1:2)]
tableE[c(3:4),3]= gmm.b.fit$coef[2,c(1:2)]
tableE[c(1:2),4]= gmm.b.cue.fit$coef[1,c(1:2)]
tableE[c(3:4),4]= gmm.b.cue.fit$coef[2,c(1:2)]

colnames(tableE)= c("GMM-(a)", "CUE-(a)", "GMM-(b)", "CUE-(b)")
rownames(tableE)= c("estimated beta", "std error beta",  "estimated alpha", "std error for alpha")

# Output tables in latex format
print(xtable(tableE, digits=4, align=c("l", rep("c", 4))), floating = FALSE, file="output/ccapm_hw_coefs.tex")
print(xtable(tableJ, digits=3, align=c("l", rep("c", 2))), floating = FALSE, file="output/ccapm_hw_J.tex")


# Diagnostics
df <- as.data.frame(x)
colnames(df) <- c("Rc","Rb","Rm","int","Rb_l1", "Rb_l2")
first_stage <- lm(Rc ~ Rb_l1 + Rb_l2, data = df)
summary(first_stage)
