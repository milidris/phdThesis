########################################
# Required Libraries
library(sensitivity)
library(triangle)
library(VineCopula)
library(evd)
library(parallel)

#####################################################
# Simulation function for Flood Case

#Simulate according to bivariate Gaussian Copula
gausCopSim<-function(n, rho){
  BiCopSim(N=n,
           family=1,
           par=rho)
}

#Simulate conditional Gaussian copula
condGausCopSim<-function(n, rho, P_xC){
  BiCopCondSim(N=n,
               family=1,
               cond.val=P_xC,
               cond.var=2,
               par=rho)
}

#Transform uniform sample to marginal
unifToMarg<-function(x, var){
  if(var==1){
    res=qgumbel.trunc(x, loc=1013, scale=558, min=500,max=3000)
  }else if(var==2){
    res=qnorm.trunc(x, mean=30, sd=7, min=15, max=55)
  }else if(var==3){
    res=qtriangle(x, a=49,c=50,b=51)
  }else if(var==4){
    res=qtriangle(x, a=54,c=55,b=56)
  }else if(var==5){
    res=qtriangle(x, a=4990, c=5000, b=5010)
  }else if(var==6){
    res=qtriangle(x, a=295, c=300, b=305)
  }
  
  return(res)
}

#Transform a marginal sample to uniform
margToUnif<-function(x, var){
  if(var==1){
    res=pgumbel.trunc(x, loc=1013, scale=558, min=500,max=3000)
  }else if(var==2){
    res=pnorm.trunc(x, mean=30, sd=7, min=15, max=55)
  }else if(var==3){
    res=ptriangle(x, a=49,c=50,b=51)
  }else if(var==4){
    res=ptriangle(x, a=54,c=55,b=56)
  }else if(var==5){
    res=ptriangle(x, a=4990, c=5000, b=5010)
  }else if(var==6){
    res=ptriangle(x, a=295, c=300, b=305)
  }
  
  return(res)
}


#Simulation according to the joint distribution
sim.RWL<-function(n){
  #First Dependency Bloc
  DepBloc1<-gausCopSim(n, 0.5)
  Q<-qgumbel.trunc(DepBloc1[,1], loc=1013, scale=558, min=500,max=3000)
  Ks<-qnorm.trunc(DepBloc1[,2], mean=30, sd=7, min=15, max=55)
  
  #Second Dependency Bloc
  DepBloc2<-gausCopSim(n, 0.3)
  Zv<-qtriangle(DepBloc2[,1], a=49,c=50,b=51)
  Zm<-qtriangle(DepBloc2[,2], a=54,c=55,b=56)
  
  #Third Dependency Bloc
  DepBloc3<-gausCopSim(n, 0.3)
  L<-qtriangle(DepBloc3[,1], a=4990, c=5000, b=5010)
  B<-qtriangle(DepBloc3[,2], a=295, c=300, b=305)
  
  #Independent Variables
  
  res<-as.data.frame(cbind(Q,Ks,Zv,Zm,L,B))
  colnames(res)=c("Q","Ks","Zv","Zm","L","B")
  return(as.matrix(res))
}

#Conditional simulations
simCond.RWL<-function(n, Sj, Sjc, xjc){
  #Simulation according to everything
  if(is.null(Sjc)){
    return(as.matrix(sim_flood(n)))
  }else{
    xjc=xjc[Sjc]
    xjc=matrix(xjc, ncol=length(Sjc))
    xjc_unif=xjc
    for(i in 1:length(Sjc)){
      xjc_unif[,i]=margToUnif(x=xjc[,i], var=Sjc[i])
    }
    res=matrix(NA, nrow=n, ncol=length(c(Sj, Sjc)))
    #Which dependency bloc ?
    
    #1st Dependency Bloc
    if(!all(c(1,2)%in%Sjc) & any(c(1,2) %in%Sjc)){ 
      condVar=which(Sjc %in% c(1,2))
      givenVar=setdiff(c(1,2), Sjc[condVar])
      unif_B1=condGausCopSim(n,0.5,xjc_unif[condVar])
      res[,givenVar]=sapply(unif_B1, unifToMarg,
                            var=givenVar)
    }else if(all(c(1,2)%in%Sj)){
      unif_B1=gausCopSim(n, 0.5)
      res[,1]=sapply(unif_B1[,1],unifToMarg,
                     var=1)
      res[,2]=sapply(unif_B1[,2],unifToMarg,
                     var=2)
    }
    
    #2nd Dependency Bloc
    if(!all(c(3,4)%in%Sjc) & any(c(3,4) %in%Sjc)){ 
      condVar=which(Sjc %in% c(3,4))
      givenVar=setdiff(c(3,4), Sjc[condVar])
      unif_B1=condGausCopSim(n,0.3,xjc_unif[condVar])
      res[,givenVar]=sapply(unif_B1, unifToMarg,
                            var=givenVar)
    }else if(all(c(3,4)%in%Sj)){
      unif_B1=gausCopSim(n, 0.3)
      res[,3]=sapply(unif_B1[,1],unifToMarg,
                     var=3)
      res[,4]=sapply(unif_B1[,2],unifToMarg,
                     var=4)
    }
    
    #3rd Dependency Bloc
    if(!all(c(5,6)%in%Sjc) & any(c(5,6) %in%Sjc)){ 
      condVar=which(Sjc %in% c(5,6))
      givenVar=setdiff(c(5,6), Sjc[condVar])
      unif_B1=condGausCopSim(n,0.3,xjc_unif[condVar])
      res[,givenVar]=sapply(unif_B1, unifToMarg,
                            var=givenVar)
    }else if(all(c(5,6)%in%Sj)){
      unif_B1=gausCopSim(n, 0.3)
      res[,5]=sapply(unif_B1[,1],unifToMarg,
                     var=5)
      res[,6]=sapply(unif_B1[,2],unifToMarg,
                     var=6)
    }
  }
  
  for(i in 1:length(Sjc)){
    res[,Sjc[i]]=xjc[i]
  }
  return(res)
}

# River Water Level Model
model.RWL<-function(X){
  H<-(X[,1]/(X[,2]*X[,6]*(sqrt((X[,4]-X[,3])/X[,5]))))^(3/5)
  G<-X[,3]+H
  return(G)
}




