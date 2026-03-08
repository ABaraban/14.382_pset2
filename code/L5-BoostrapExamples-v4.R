
library(boot)

# function to obtain R-Squared from the data 

# Example
#rsq <- function(formula, data, indices) {
#  d <- data[indices,] # allows boot to select sample 
#  fit <- lm(formula, data=d)
#  return(summary(fit)$r.square)
#} 
# bootstrapping with 1000 replications 
#results <- boot(data=mtcars, statistic=rsq, 
#                R=1000, formula=mpg~wt+disp)
# view results
#results 
#plot(results)
# get 95% confidence interval 
#boot.ci(results, type="bca")



### example of EBS working

set.seed(1)
postscript("/Users/chernoz/Dropbox/TEACHING/14.382/Results/Sim_mean.eps",horizontal=F, pointsize=15,width=8.0,height=8.0)
par(mfrow=c(2,1))
simulate.mean<- function(i){ return(mean(rexp(100)))}
mean.draws<- rep(0, length=1000)
for  (i in 1:1000) {  mean.draws[i]<- simulate.mean(i) }
plot(hist(mean.draws,nclass=30,plot=FALSE), col="light blue", main="histogram for true draws", xlab="values of statistic")
se<- round(sqrt(var(mean.draws)), digits=4)
legend(1.2, 60, paste("se=", se))
meanstatistic<- function(data,indices){ return(mean(data[indices])) }
data.fake<- rexp(100)
result<- boot(data=data.fake, statistic=meanstatistic, R=1000)
plot(hist(result$t,nclass=30, plot=FALSE), col="blue", main="histogram for BS draws", xlab="values of statistic")
se<- round(sqrt(var(result$t)), digits=4)
legend(1.1, 60, paste("se=", se))

dev.off()



### example of EBS failing

set.seed(1)
postscript("/Users/chernoz/Dropbox/TEACHING/14.382/Results/Sim_min.eps",horizontal=F, pointsize=15,width=8.0,height=8.0)
par(mfrow=c(2,1))
simulate.min<- function(i){ return(min(rexp(100)))}
min.draws<- rep(0, length=1000)
for  (i in 1:1000) {  min.draws[i]<- simulate.min(i) }
se<- round(sqrt(var(min.draws)),digits=4)
plot(hist(min.draws,nclass=30, plot=FALSE), col="light blue", main="histogram for true draws", xlab="values of statistic")
legend(.06, 100,  paste("se=", se))
minstatistic<- function(data,indices){ return(min(data[indices])) }
data.fake<- rexp(100)
result<- boot(data=data.fake, statistic=minstatistic, R=1000)
plot(hist(result$t,nclass=30, plot=FALSE), col=2, main="histogram for BS draws", xlab="values of statistic")
se<- round(sqrt(var(result$t)),digits=4)
legend(.03, 300,  paste("se=", se))
dev.off()


### BS for Hansen-Singleton

library(AER)
library(sandwich)
library(gmm)
library(quantmod)


Ret<- read.csv(file="data/ccapm-ready-to-use.csv", header=T)
attach(Ret)
y  <- na.omit(cbind(Rc, Rb, Rm))

# bootstrap GMM for technical instruments consisting of one lag



nlasgs=1

z    <- cbind(1, Lag(y[ ,1], c(1:nlags)), Lag(y[ ,2], c(1:nlags)), Lag(y[ ,3], c(1:nlags)))

# also tried: 
# nlags=2
# z<- cbind(rep(1, dim(y)[1]))  Trying this would give m=n, exactly ID-d case or calibration
# to explain "risk premium"

x   <- na.omit(cbind(y,z))


g2b.x.theta <- function(theta, x)
{
  y    <- x[ ,1:3]
  z    <- x[ ,-c(1:3)]
  return (cbind( (theta[1] * y[ ,2] * y[ ,1]^(-theta[2]) - 1)*z, 
                 (theta[1] * y[ ,3] * y[ ,1]^(-theta[2]) - 1)*z))
}

gmm2 <- gmm(g2b.x.theta, x, t0 = c(.99,.22), type="iterative", vcov="iid")
summary(gmm2)
confint(gmm2)


# function for Empirical Bootstap

gmm.bs.statistic<- function(data,indices){ 
  x<- data[indices,]
  gmm2.coef <- gmm(g2b.x.theta, x, t0 = c(.99,.22), type="iterative", vcov="iid")$coef
  return(gmm2.coef)
}

# function for Block Bootstrap


gmm.ts.bs.statistic<- function(data){ 
  x<- data
  gmm2.coef <- gmm(g2b.x.theta, x, t0 = c(.99,.22), type="iterative", vcov="iid")$coef
  return(gmm2.coef)
}


# bootstrap histograms and ses


set.seed(1)
result.gmm.bs<- boot(data=x, statistic=gmm.bs.statistic,R=400)
ses<- round(sqrt(apply(result.gmm.bs$t, 2, var)),digits=4)



# postscript("/Users/chernoz/Dropbox/TEACHING/14.382/Results/BS_GMM_indep.eps",horizontal=F, pointsize=15,width=8.0,height=8.0)
par(mfrow=c(2,1))
plot(hist(result.gmm.bs$t[,1],nclass=23, plot=FALSE), col="light blue", main="histogram for BS draws", xlab="values of statistic")
legend(.999,100, c("se=", paste(ses[1])))
plot(hist(result.gmm.bs$t[,2],nclass=23, plot=FALSE), col="light blue", main="histogram for BS draws", xlab="values of statistic")
legend(.4,100, c("se=", paste(ses[2])))
dev.off()


#confint(gmm(g2.x.theta, x, t0 = c(.99,3), type="iterative", vcov="iid"), level=.90)
quantile(result.gmm.bs$t[,1], c(0.05,.95))
quantile(result.gmm.bs$t[,2], c(.05,.95))


# This employs block bootstap with block size equal to 26

set.seed(1)
result.gmm.bs.ts<- tsboot(tseries=x, statistic=gmm.ts.bs.statistic, l=26, sim="fixed", R=200)
ses.ts<- round(sqrt(apply(result.gmm.bs.ts$t, 2, var)),digits=4)

# postscript("/Users/chernoz/Dropbox/TEACHING/14.382/Results/BS_GMM_tes.eps",horizontal=F, pointsize=15,width=8.0,height=8.0)
par(mfrow=c(2,1))
plot(hist(result.gmm.bs.ts$t[,1],nclass=23, plot=FALSE), col="light blue", main="histogram for BS draws", xlab="values of statistic")
legend(.999,20, c("se=", paste(ses.ts[1])))
plot(hist(result.gmm.bs.ts$t[,2],nclass=23, plot=FALSE), col="light blue", main="histogram for BS draws", xlab="values of statistic")
legend(.5,40, c("se=", paste(ses.ts[2])))
dev.off()

###############################################################################
############################## Homework #######################################
###############################################################################

library(boot)
library(tidyverse)

# the score function for estimationg CRRA (power utility) preferences with 2 financials assets
g2b.x.theta = function(theta, x)
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

gmm.bs.statistic<- function(data,indices){ 
  x<- data[indices,]
  gmm2.coef <- gmm(g2b.x.theta, x, t0 = c(.99,.22), type="iterative", vcov="iid")$coef
  return(gmm2.coef)
}

# Data
Ret = read_csv(file="data/ccapm-long_mert.csv")
attach(Ret)

# (a) Instruments = intercept
y    = na.omit(cbind(Rc, Rb, Rm)) # outcomes
z    = 1
x    = na.omit(cbind(y,z))

# Bootstrap draws
set.seed(1)
gmm_bs_intercept <- boot(data=x, statistic=gmm.bs.statistic,R=400)
gmm_ses_intercept <- summary(gmm_bs_intercept)$bootSE
print(summary(gmm_bs_intercept))

# (b) Instruments = 1 lag of stock return and 1 lag of bond returns
y     = na.omit(cbind(Rc, Rb, Rm))
nlags = 1
z    = cbind(1, Lag(y[ ,2], c(1:nlags)), Lag(y[ ,3], c(1:nlags)))
x   = na.omit(cbind(y,z))

# Bootstrap draws
gmm_bs_returns <- boot(data=x, statistic=gmm.bs.statistic,R=400)
gmm_ses_returns <- summary(gmm_bs_returns)$bootSE
print(summary(gmm_bs_returns))


# Replicate lecture 

# Instruments = 1 lag of consumption, bond, and stock returns
y     = na.omit(cbind(Rc, Rb, Rm))
nlags = 1
z    = cbind(1, Lag(y[ ,1], c(1:nlags)), Lag(y[ ,2], c(1:nlags)), Lag(y[ ,3], c(1:nlags)))
x    = na.omit(cbind(y,z))

# Bootstrap draws
gmm_bs_returns_cons <- boot(data=x, statistic=gmm.bs.statistic,R=400)
gmm_ses_returns_cons <- summary(gmm_bs_returns_cons)$bootSE
print(summary(gmm_bs_returns_cons))


# Display coefficients in a table
table_bootstrap = tibble(
  "Instruments" = c("Intercept only", "Intercept + stock/bond returns"),
  "Beta SE" = c(gmm_ses_intercept[1], gmm_ses_returns[1]),
  "Beta CI" = c(
    paste0("[", as.character(round(gmm_bs_intercept$t0[1] - 1.96 * gmm_ses_intercept[1],3)),", ", as.character(round(gmm_bs_intercept$t0[1] +  1.96 * gmm_ses_intercept[1],3)), "]"),
    paste0("[", as.character(round(gmm_bs_returns$t0[1] - 1.96 * gmm_ses_returns[1],3)),", ", as.character(round(gmm_bs_returns$t0[1] +  1.96 * gmm_ses_returns[1],3)), "]")
    ),
  "Alpha SE" = c(gmm_ses_intercept[2], gmm_ses_returns[2]),
  "Alpha CI" = c(
    paste0("[", as.character(round(gmm_bs_intercept$t0[2] - 1.96 * gmm_ses_intercept[2],3)),", ", as.character(round(gmm_bs_intercept$t0[2] +  1.96 * gmm_ses_intercept[2],3)), "]"),
    paste0("[", as.character(round(gmm_bs_returns$t0[2] - 1.96 * gmm_ses_returns[2],3)),", ", as.character(round(gmm_bs_returns$t0[2] +  1.96 * gmm_ses_returns[2],3)), "]")
  )
)

print(xtable(table_bootstrap, digits=3, align=c("l", "c", "c", "c", "c", "c"), file="output/ccapm_hw_bs.tex"))


# Histograms
ggplot(mapping = aes(x=x)) + 
  geom_histogram(data = tibble(x= gmm_bs_intercept$t[,2])) +
  labs(
    title = "Z = (Constant)",
    x="", y=""
  ) +
  theme_minimal()
ggsave("output/bs_plot1.png", width = 4, height = 3, units = "in")


ggplot(mapping = aes(x=x)) + 
  geom_histogram(data = tibble(x= gmm_bs_returns$t[,2]), fill="red") + 
  labs(
    title = "Z = (Constant, Bond Returns, Stock Returns)",
    x="", y=""
  ) +
  theme_minimal()
ggsave("output/bs_plot2.png", width = 4, height = 3, units = "in")


ggplot(mapping = aes(x=x)) + 
  geom_histogram(data = tibble(x= gmm_bs_returns_cons$t[,2]), fill="blue") +
  labs(
    title = "Z = (Constant, Consumption, Bond Returns, Stock Returns)",
    x="", y=""
  ) +
  theme_minimal()
ggsave("output/bs_plot3.png", width = 4, height = 3, units = "in")
  

  















