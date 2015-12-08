###################################################################################

#     TAPER MODEL AND BIOMASS ESTIMATION 
#     using a Bayesian hierarchial inference approach to propagate uncertainty 

###################################################################################

###load packages

library(R2jags)
library(ggplot2)
library(dplyr)

###################################################################################
###data: use simulated data by calling 
source("r_code/sim_taper_data.R")

#format data
#1.   number of observations
      Nobs=nrow(census_data)
#2.   vector of diameter measurements
      diam=census_data$dbh
#3.   vector of reference heights for point of measurement for diameter      
      POM=census_data$POM
#4.   species id for data
      species=census_data$species.id
#5.   number of species
      n.sp=length(unique(species))

#store data as a list for jags model
dat<-list(diam=diam, Nobs=Nobs,POM=POM,species=species,n.sp=n.sp)

#parameters to track
params<-c("beta0","beta1","mu","sigma1","mu.alpha","mu.beta","sigma2")

###################################################################################
### jags model structure

cat("model {
    
    for (i in 1:Nobs) {
    
    diam[i] ~ dnorm(mu[i],tau)
    mu[i] <- beta0[species[i]] * exp(-beta1[species[i]] * (POM[i] - 1.3))
      
    }
      #species level random effects
      for (j in 1:n.sp) {
      beta0[j] ~ dnorm(mu.alpha, tau1)
      beta1[j] ~ dnorm(mu.beta,tau2)
      }
      
      #Priors for hyper-parameters    
      mu.alpha ~ dnorm(0,0.0001)
      mu.beta ~ dnorm(0,0.0001)
      
      sigma1 ~ dunif (0,100)
      tau1 <- 1/(sigma1 * sigma1)
    
      sigma2 ~ dunif (0,100)
      tau2 <- 1/(sigma2 * sigma2)

      #residual variance
      sigma ~ dunif (0,100)
      tau <- 1/(sigma * sigma)
      
    }",file="outputs/taper_species_model.txt")

taper.sp.model<-jags(dat=dat, inits=NULL, parameters.to.save=params, "outputs/taper_species_model.txt",n.chains=3, n.iter=3000, n.burnin=1000, n.thin=50)


#if models fails to fit, update the number of iterations
#taper.sp.model2 <- update(taper.sp.model,niter=2000)

#save jags object output
#save(taper.sp.model,file="outputs/taper.sp.model.Rdata")

#load jags object
#load("outputs/taper.sp.model.Rdata")

###################################################################################
#source functions to be used
source(file="r_code/functions.R")


#model convergence assessment using rhat statistic function

Rhat(taper.sp.model) 
#traceplot(taper.sp.model) #look at mcmc chains
plot(taper.sp.model) #display outputs

#outputs for nodes monitored for the MCMC chains using the gipar function 
taper_fit <- gipar(taper.sp.model)
str(taper_fit)
summary(taper_fit)


###################################################################################
### Estimate biomass applying the Metcalf et al. 2009 taper correction

#empty matrix to store species specific intercept and slope parameters
a <-  matrix(NA,nrow=nrow(taper_fit$beta0),ncol=6)
b <-  matrix(NA,nrow=nrow(taper_fit$beta1),ncol=6)

#store jags posterior draws in matrix 

#extract parameter a
for (i in 1:nrow(a)) {
            
      for (j in 1:ncol(taper_fit$beta0)) {
                  
                  a[i,j] <- taper_fit$beta0[i,j]
            }
}

#extract parameter b
for (i in 1:nrow(b)) {
      
      for (j in 1:ncol(taper_fit$beta1)) {
            
            b[i,j] <- taper_fit$beta1[i,j]
      }
}


#empty dataframe to store the diameters estimated at 1.3 m height
dbh_pred <- matrix(NA,nrow=length(POM),ncol=nrow(a)) 

#for-loop to propagate the uncertainty in diameter associated with taper model parameters
for (i in 1:nrow(dbh_pred)) {
            for (j in 1:nrow(a)) {
      dbh_pred[i,j] <- dbh_1.3m(POM[i],a[j,species[i]],b[j,species[i]])
            }
}
dbh_error <- data.frame(t(apply(dbh_pred, 1, quantile, probs=c(.25,.5, .75), na.rm=TRUE)))

#calculate biomass using corrected diameter
biomass_taper <- Chave.dbh(DBH=dbh_pred,0.6,1.3)  
biomass_error <- data.frame(t(apply(biomass_taper, 1, quantile, probs=c(.25,.5, .75), na.rm=TRUE)))
#uncorrected diameter
biomass_notaper <- Chave.dbh(DBH=diam,0.6,1.3)



#Plot species level parameters against simulated values
param_a <- data.frame(a)
colnames(param_a) <- c("species_1","species_2","species_3","species_4","species_5","species_6")
param_a2 <- reshape(param_a, varying = c("species_1","species_2","species_3","species_4","species_5","species_6"),
                    v.names="a", timevar="species", times = c("species_1","species_2","species_3","species_4","species_5","species_6"),
                    direction="long") 

vline_a <- data.frame(a=c(180,160,150,175,180,170), species = c("species_1","species_2","species_3","species_4","species_5","species_6"))
ggplot(data=param_a2,aes(x=a)) +
      geom_density(size=1) +
      geom_vline(data=vline_a, aes(xintercept=a),col="darkgreen") +
      facet_wrap(~species) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_x_continuous("Parameter a")




#Plot biomass against diameter estimates at 1.3 reference height
plot(biomass_error[,1]~dbh_error[,1],col="gray",xlab="Diameter (cm)",ylab="Biomass (Mg)",ylim=c(0,15),main="Diameter estimated at 1.3m height",type="l", lwd=0.6)
lines(biomass_error[,3]~dbh_error[,3],col="gray",lwd=0.6)
lines(biomass_error[,1]~dbh_error[,1],col="black",lwd=0.1)

#plot biomass against measured diameter
plot(biomass_notaper~diam,col="gray",xlab="Diameter (cm)",ylab="Biomass (Mg)",ylim=c(0,30),main="Measured diameter",type="p",cex=0.2)


###################################################################################

