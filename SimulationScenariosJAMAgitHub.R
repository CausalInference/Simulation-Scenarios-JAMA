library(simsurv)
library(survival)
library(xtable)


############### SCENARIO 1 #########################

set.seed(1)
#### Input
lambda = 0.0001
lambdacens <- lambda*3 
N = 50000
b0 <- log(1)
b1 <- log(1.0005)
b1cens <- 0 
followUp <- 365*4


#### Construct hazard functions for plotting.
hazardC1 <- function(tt,lamb){return(lamb)}
hazardT1 <- function(tt,lamb, B0, B1){return(lamb*exp(B0+B1*tt))}
timeseq <- seq(0,followUp, length.out = 1e4)
hazardvecC1 <- sapply(timeseq, hazardC1, lambda)
hazardvecT1 <- sapply(timeseq, hazardT1, lambda, b0, b1)
#### Suitible standardisation of hazard vec. 
stdhazardvecC1 <- sapply(timeseq, hazardC1, lambda)/lambda * 0.2 + 0.1
stdhazardvecT1 <- sapply(timeseq, hazardT1, lambda, b0, b1)/lambda * 0.2 +0.1
covariates <- data.frame(z=c(rep(1, N/2),rep(0, N/2)))
covariatesCens <- data.frame(z=c(rep(1, N)))
#### Input done
sim1 <- simsurv(dist="exponential", lambdas = lambda, x=covariates,betas = c(z=b0), tde=c(z=b1),max=followUp, interval=c(1e-8,followUp+1))
surv1 <- survfit(Surv(sim1$eventtime, sim1$status)~covariates$z)
cens1 <- simsurv(dist="exponential",lambdas=lambdacens, x=covariatesCens,betas = c(z=b0), tde=c(z=b1cens),max=followUp, interval=c(1e-8,followUp+1))
censim1 <- sim1
censim1$eventtime[which(sim1$eventtime > cens1$eventtime)] <- cens1$eventtime[which(sim1$eventtime > cens1$eventtime)]
censim1$status[sim1$eventtime > cens1$eventtime] <- rep(0,sum(sim1$eventtime > cens1$eventtime))
censsurv1 <- survfit(Surv(censim1$eventtime, censim1$status)~covariates$z)
censsurv1$time <- censsurv1$time/365
surv1$time <- surv1$time/365

####
cox1 <- coxph(Surv(sim1$eventtime, sim1$status)~covariates$z + cluster(1:N))
summary(cox1)
cox1cens <- coxph(Surv(censim1$eventtime, censim1$status)~covariates$z)
summary(cox1cens)
sum(sim1$eventtime > cens1$eventtime)

############### SCENARIO 1 FINISHED #########################

############### SCENARIO 2 #########################
set.seed(1)
#### Input
lambda = 0.0001
lambdacens <- lambda*2.2 # 4 good choice: lambda/1.2
N = 50000
b0 <- log(2.25)
b1 <- log(0.997)
b1cens <- 0 
followUp <- 365*4

#### Construct hazard functions for plotting.
hazardC2 <- function(tt,lamb){return(lamb)}
hazardT2 <- function(tt,lamb, B0, B1){return(lamb*exp(B0+B1*tt))}
timeseq <- seq(0,followUp, length.out = 1e4)
hazardvecC2 <- sapply(timeseq, hazardC2, lambda)
hazardvecT2 <- sapply(timeseq, hazardT2, lambda, b0, b1)
#### Suitible standardisation of hazard vec. 
stdhazardvecC2 <- sapply(timeseq, hazardC2, lambda)/lambda * 0.2 + 0.2
stdhazardvecT2 <- sapply(timeseq, hazardT2, lambda, b0, b1)/lambda * 0.2 +0.2


covariates <- data.frame(z=c(rep(1, N/2),rep(0, N/2)))
covariatesCens <- data.frame(z=c(rep(1, N)))
#### Input done, now everything is functions of the expressions above.
sim2 <- simsurv(dist="exponential", lambdas = lambda, x=covariates,betas = c(z=b0), tde=c(z=b1),max=followUp, interval=c(1e-8,followUp+1))
surv2 <- survfit(Surv(sim2$eventtime, sim2$status)~covariates$z)
cens2 <- simsurv(dist="exponential",lambdas=lambdacens, x=covariatesCens,betas = c(z=b0), tde=c(z=b1cens),max=followUp, interval=c(1e-8,followUp+1))
censim2 <- sim2
censim2$eventtime[which(sim2$eventtime > cens2$eventtime)] <- cens2$eventtime[which(sim2$eventtime > cens2$eventtime)]
censim2$status[sim2$eventtime > cens2$eventtime] <- rep(0,sum(sim2$eventtime > cens2$eventtime))
censsurv2 <- survfit(Surv(censim2$eventtime, censim2$status)~covariates$z)
censsurv2$time <- censsurv2$time/365
surv2$time <- surv2$time/365

#### Fit Cox models
cox2 <- coxph(Surv(sim2$eventtime, sim2$status)~covariates$z + cluster(1:N))
summary(cox2)
cox2cens <- coxph(Surv(censim2$eventtime, censim2$status)~covariates$z)
summary(cox2cens)

sum(sim2$eventtime > cens2$eventtime)


############### SCENARIO 2 FINISHED #########################

############### SCENARIO 3 #########################

set.seed(1)
#### Input
lambda =  0.00015 
lambdacens <- lambda*2.5 #4
N = 50000
b0 <- log(1)
b1 <- log(6) 
b1cens <- 0 
followUp <- 365*4

scalingFactor <- 0.0015
scalingFactor2 <- 0.3
daysToEquality = 365*2
timeDependentHR <- function(x) scalingFactor2*(exp(scalingFactor*(daysToEquality-x))-1)

#### Construct hazard functions for plotting.
hazardC3 <- function(tt,lamb){return(lamb)}
hazardT3 <- function(tt,lamb, B0, B1){return(lamb*exp(B1*timeDependentHR(tt)))}
timeseq <- seq(0,followUp, length.out = 1e4)
hazardvecC3 <- sapply(timeseq, hazardC3, lambda)
hazardvecT3 <- sapply(timeseq, hazardT3, lambda, b0, b1)
#### Suitible standardisation of hazard vec. 
stdhazardvecC3 <- sapply(timeseq, hazardC3, lambda)/lambda * 0.2 + 0.2
stdhazardvecT3 <- sapply(timeseq, hazardT3, lambda, b0, b1)/lambda * 0.2 +0.2
covariates <- data.frame(z=c(rep(1, N/2),rep(0, N/2)))
covariatesCens <- data.frame(z=c(rep(1, N)))
#### Input done, now everything is functions of the expressions above.
sim3 <- simsurv(dist="exponential", lambdas = lambda, x=covariates,betas = c(z=b0), tde=c(z=b1),
                tdefunction = timeDependentHR,
                max=followUp, interval=c(1e-8,followUp+1))


surv3 <- survfit(Surv(sim3$eventtime, sim3$status)~covariates$z)
cens3 <- simsurv(dist="exponential",lambdas=lambdacens, x=covariatesCens,betas = c(z=b0), tde=c(z=b1cens),max=followUp, interval=c(1e-8,followUp+1))
censim3 <- sim3
censim3$eventtime[which(sim3$eventtime > cens3$eventtime)] <- cens3$eventtime[which(sim3$eventtime > cens3$eventtime)]
censim3$status[sim3$eventtime > cens3$eventtime] <- rep(0,sum(sim3$eventtime > cens3$eventtime))
censsurv3 <- survfit(Surv(censim3$eventtime, censim3$status)~covariates$z)
censsurv3$time <- censsurv3$time/365
surv3$time <- surv3$time/365

#### Fit Cox models
cox3 <- coxph(Surv(sim3$eventtime, sim3$status)~covariates$z + cluster(1:N))
summary(cox3)
cox3cens <- coxph(Surv(censim3$eventtime, censim3$status)~covariates$z)
summary(cox3cens)

############### SCENARIO 3 FINISHED #########################

###############  Output Table (For the Appendix) #########################

#Differences in survival at times t0
getSurvivalDifference <- function(inputSurv, t0){
  t0diff <- diff(summary(inputSurv,times = t0)$surv)
  t0diffSd <- sqrt(sum(summary(inputSurv,times = t0)$std.err^2))
  t0diffU <- t0diff + 1.96 * t0diffSd
  t0diffL <- t0diff - 1.96 * t0diffSd
  return(c(t0diff,t0diffL,t0diffU))
}

outputTable1 <- cbind(rbind(1/summary(cox1)$conf.int[c(1,4,3)], 1/summary(cox1cens)$conf.int[c(1,4,3)]),
                      100*rbind(getSurvivalDifference(surv1,3),getSurvivalDifference(censsurv1,3)))
outputTable2 <- cbind(rbind(summary(cox2)$conf.int[c(1,3,4)],summary(cox2cens)$conf.int[c(1,3,4)]),
                      100*rbind(getSurvivalDifference(surv2,3),getSurvivalDifference(censsurv2,3)))
outputTable3 <- cbind(rbind(summary(cox3)$conf.int[c(1,3,4)],summary(cox3cens)$conf.int[c(1,3,4)]),
                      100*rbind(getSurvivalDifference(surv3,3),getSurvivalDifference(censsurv3,3)))

####
allTables <- rbind(outputTable1[,c(1:4,6,5)],outputTable2,outputTable3)
library(xtable)
xtable(allTables,digits=2,caption = "New Numbers")




