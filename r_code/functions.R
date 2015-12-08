###################################################################################

#FILE CONTAINS FUNCTIONS USED IN THE TAPER MODEL AND BIOMASS ESTIMATION

###################################################################################

#1. Returns MCMC chains for parameters monitored from jags ouput
gipar<-function(model) 
      #Args:
            #model = name of jags model output
{
      return(model[2]$BUGSoutput$sims.list)
      #output: returns a dataframe of posterior draws for nodes monitored from jags model
}

#2. Chave et al. 2014 diameter ~ biomass alometric equation
Chave.dbh <- function(DBH,WSG,E) {
      #Args:
            #DBH = vector of tree diameters
            #WSG = wood specific gravity
            #E = Chave et al. 2014 environment variable
      exp(-1.803 - 0.976*E + 0.976*log(WSG) + 2.673*log(DBH) - 0.0299*(log(DBH)^2))/1000
      #output = biomass in Megagrams (Mg)/~metric tons
}

#3. Function to estimate diameter at 1.3 m reference height; Metcalf et al. 2009
dbh_1.3m <- function (POM,a,b){
      #Args:
            #POM = vector of reference heights where diameter was recorded
            #a = posterior draw for parameter 1
            #b = posterior draw for parameter 2
            a * exp(-b*(POM-1.3))
}
      
#4. Assess model convergence function
Rhat<-function(model) {
      #Args:
            #model = jags mcmc model ouput
      rhat = model[[2]]$summary[,8]
      hist(rhat, breaks=100)
      print(paste("max Rhat =",max(rhat)))
      #outputs:
            #1. histogram of rhat values
            #2. maximum rhat value
}



