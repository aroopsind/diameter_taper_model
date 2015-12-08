###################################################################################

#     SIMULATE DATA FOR TAPER MODEL AND BIOMASS ESTIMATION 
#     Outputs: simulated data has species level and individual effects 

###################################################################################

n.plants=100 #100 trees per species
n.measures=10 #10 diameter measurements per tree
n.sp=6 #6 species

#species level parameters to be used to estimate taper relationship
a.SP=c(180,160,150,175,180,170)
b.SP=c(0.01,0.015,0.02,0.016,0.017,0.009)

#object to store simulated data as a list
species.list=vector("list",n.sp) 

#For loop to simulate data 

for(j in 1:n.sp) {
      
      
      #simulated reference heights (POM) for diameter measurements in meters
      POM=runif(n=n.measures*n.plants,min=1.3,max=8) 
      
      #matrix to store all combinations of trees x measurements
      base.dat.0=expand.grid(c(1:n.plants),c(1:n.measures)) 
      
      #bind trees x measurments with POM
      base.dat=data.frame(POM,base.dat.0)
      colnames(base.dat)=c("POM","tree_id","measurement_id")
      
      #parameter to describe distribution used to sample trees based on POM
      population.mean=0.5
      population.sd=0.1
      
      
      #draw values from a log-transformed normal distribution for random effects for individual trees
      plant.shape=exp(rnorm(n.plants,mean=population.mean,sd=population.sd))
      
      #add individual level effects
      ind.effects=data.frame(c(1:n.plants),plant.shape)
      colnames(ind.effects)=c("tree_id","MEAN.shape")
      
      #binds ind.effects to base file
      base=merge(ind.effects,base.dat,by.x="tree_id",by.y="tree_id")
      
      simulated.data=rep(NA,times=nrow(base)) #vector to store simulated data
      mean.relationship=rep(NA,times=nrow(base)) #vector to store mean values
      for(i in 1:nrow(base)) {
            
            a=a.SP[j]+base$MEAN.shape[i] #adding individual level effect
            b=b.SP[j]
            
            mean.relationship[i]= a * exp(-b*(base$POM[i]-1.3))
            
            simulated.data[i]=rnorm(n=1,mean=mean.relationship[i],sd=22) #fill empty vector with simulated data
            
            
      }  
      
      #final dataframe ouput
      taper_data=cbind(simulated.data,base,rep(j,times=length(simulated.data)))
      colnames(taper_data)=c("dbh","tree_id","mean.shape","POM","Measurement.id","species.id")
      
      #plot of diameter~POM distributions by species
      palette(rainbow(10))
      plot(taper_data[,1]~taper_data[,4],xlab="Point of measurement (m)",ylab="Diameter (cm)",pch=19,col=taper_data[,2],xlim=c(1,8))
      
      species.list[[j]]=taper_data
}


#collapse species level data into a single dataframe
census_data=do.call("rbind",species.list)
#add unique tree id for species
census_data$main_id=paste(census_data$species.id,census_data$tree_id,sep=".") 

###################################################################################