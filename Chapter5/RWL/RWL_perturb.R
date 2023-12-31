source(file="useCases/riverWaterLevel/RWL_model.R")
source(file="Chapter3/functions/allocationsMCestimation.R")

# Loading the "interpretability" package : Uncomment if first usage

# install.packages("Chapter5/package/interpretability_0.1.0.tar.gz", 
#                  repos=NULL, 
#                  type='source')

library(interpretability)
library(sensitivity)
library(CVXR)
library(evd)
library(triangle)


#######################################################
# Baseline Shapley Effects (no perturbations)


#######################################################
# Perturbation scheme

# Perturbing L : River width
# - Preserve the median
# - Move application domain from [4990,5010] to [4988,5012]

quantile_L<-function(x){
  return(qtriangle(x, a=4990, c=5000, b=5010))
}

L.pert<-qcwProj(qFun=quantile_L,
                alpha= 0.5,
                b=5000,
                OmegaX=c(4988,5012),
                type="polynomial",
                d=12,
                int.tol=.Machine$double.eps^0.25,
                int.subdiv=1e08,
                normalize=T)

# Perturbing Zm : River width
# - Remove 0.05 to the 0.25 quantile
# - Add +0.1 to the 0.8 and 0.9 quantiles
# - Preserve the median
# - Preserve the application domain [54,56]

quantile_Zm<-function(x){
  return(qtriangle(x, a=54, c=55, b=56))
}

Zm.pert<-qcwProj(qFun=quantile_Zm,
                 alpha= c(0.25,0.5,0.8,0.9),
                 b=c(quantile_Zm(0.25)-0.05, 55, quantile_Zm(0.8)+0.1, quantile_Zm(0.9)+0.1),
                 OmegaX=c(54,56),
                 type="polynomial",
                 d=12,
                 int.tol=.Machine$double.eps^0.25,
                 int.subdiv=1e08,
                 normalize=T)

# Perturbing Q : Annual maximal water flow
# - Application domain : [500,3000] -> [500,3200]
# - Preserve the median
# - Increase by 75 the 0.15 quantile
# - Decrease by 125 the 0.75 quantile


quantile_Q<-function(x){
  return(qgumbel.trunc(x, loc=1013, scale=558, min=500,max=3000))
}

Q.pert<-qcwProj(qFun=quantile_Q,
                alpha= c(0.15,0.75),
                b= c(quantile_Q(0.15)+75, quantile_Q(0.75)-125),
                OmegaX=c(500,3200),
                type="polynomial",
                d=20,
                int.tol=.Machine$double.eps^0.25,
                int.subdiv=1e08,
                normalize=T)

# Perturbing Ks : Strickler Coefficient
# - Application domain dilatation, eta = 2

theta=seq(-1,1,0.01)

appDomainDilat<-function(theta,eta, OmegaX){
  om0=OmegaX[1]
  om1=OmegaX[2]
  
  if(theta < 0){
    b0=0.5*(om0*(2-theta*((1/eta)-1)) + theta*om1*(1/eta -1))
    b1=0.5*(om1*(2-theta*((1/eta)-1)) + theta*om0*(1/eta -1))
  }else if(theta>0){
    b0=0.5*(om0*(2+theta*(eta-1)) - theta*om1*(eta -1))
    b1=0.5*(om1*(2+theta*(eta-1)) - theta*om0*(eta -1))
  }else if(theta==0){
    b0=om0
    b1=om1
  }
  return(c(b0,b1))
}

appDomains=t(sapply(theta, appDomainDilat,
                    OmegaX=c(15,55),
                    eta=1.5))

quantile_Ks<-function(x){
  return(qnorm.trunc(x, mean=30, sd=7, min=15, max=55))
}

perturb.Ks<-function(OmegaX){
  Q.pert<-qcwProj(qFun=quantile_Ks,
                  alpha=NULL,
                  b=NULL,
                  OmegaX=OmegaX,
                  type="polynomial",
                  d=12,
                  int.tol=.Machine$double.eps^0.25,
                  int.subdiv=1e08,
                  normalize=T)
}

cl=parallel::makeCluster(8)
clusterExport(cl, varlist=list("quantile_Ks",
                               "perturb.Ks"))
clusterEvalQ(cl, library(sensitivity))
clusterEvalQ(cl, library(triangle))
clusterEvalQ(cl, library(CVXR))
clusterEvalQ(cl, library(evd))
clusterEvalQ(cl, library(interpretability))

res.pert.Ks<-parApply(cl, appDomains, 1, perturb.Ks)

stopCluster(cl)

#######################################################
# Global Robustness

# Computing Shapley Values for perturbed distributions
res.shap<-matrix(NA, ncol=7, nrow=length(res.pert.Ks))
res.pertVEs<-rep(list(0), length(res.pert.Ks))

colnames(res.shap)=c("theta","Q","Ks","Zv","Zm","L","B")

res.shap[,1]=theta

cl=parallel::makeCluster(8)
clusterExport(cl, varlist=list("condGausCopSim",
                               "gausCopSim",
                               "margToUnif",
                               "model.RWL",
                               "sim.RWL",
                               "simCond.RWL",
                               "unifToMarg",
                               "estim.CE.MC.pert",
                               "L.pert",
                               "Q.pert",
                               "Zm.pert",
                               "res.pert.Ks")
)
clusterEvalQ(cl, library(sensitivity))
clusterEvalQ(cl, library(triangle))
clusterEvalQ(cl, library(VineCopula))
clusterEvalQ(cl, library(evd))
clusterEvalQ(cl, library(interpretability))
doParallel::registerDoParallel(cl)

#Seed for reproducibility
set.seed(12345)

for(theta_idx in 1:length(res.pert.Ks)){
  cat("-------------------------------------------\n")
  cat("Perturbation n°: ",theta_idx, "\n", sep="")
  
  startRep=Sys.time()
  margToPert<-function(x, var){
    if(var==1){
      res=evaluate(Q.pert,margToUnif(x,1))
    }else if(var==2){
      res=evaluate(res.pert.Ks[[theta_idx]], margToUnif(x,2))
    }else if(var==3){
      res=x
    }else if(var==4){
      res=evaluate(Zm.pert, margToUnif(x,4))
    }else if(var==5){
      res=evaluate(L.pert, margToUnif(x,5))
    }else if(var==6){
      res=x
    }
    return(res)
  }
  
  pertMap<-function(X){
    res=X
    for(i in 1:6){
      res[,i]=margToPert(X[,i], i)
    }
    return(res)
  }
  
  clusterExport(cl, varlist=list("margToPert",
                                 "pertMap",
                                 "theta_idx")
  )
  VEs=estim.CEs.MC.pert(model=model.RWL, 
                            jointSim=sim.RWL, 
                            condSim=simCond.RWL,
                            Nv=1e4, 
                            No=1e3,
                            Ni=100,
                            parl=T,
                            cl=cl,
                            pertMapping=pertMap)
  
  res.pertVEs[[theta_idx]]=VEs
  res.shap[theta_idx,-1]=shFromCEs(VEs)$Shap
  endRep=Sys.time()
  
  cat("Computation time:", difftime(endRep, startRep, units="mins"), "mins\n", sep=" ")
}

stopCluster(cl)

#######################################################
# Mean, Var and Quantiles of the output
n_sim=4e5
res.stats<-matrix(NA, ncol=9, nrow=length(res.pert.Ks))
colnames(res.stats)=c("theta","mean","var","min","q025","q05", "q95","q975","max")

res.stats[,1]=theta

#Reproducibility
set.seed(12345)

for(theta_idx in 1:length(res.pert.Ks)){
  cat("-------------------------------------------\n")
  cat("Perturbation n°: ",theta_idx, "\n", sep="")
  startRep=Sys.time()
  
  margToPert<-function(x, var){
    if(var==1){
      res=evaluate(Q.pert,margToUnif(x,1))
    }else if(var==2){
      res=evaluate(res.pert.Ks[[theta_idx]], margToUnif(x,2))
    }else if(var==3){
      res=x
    }else if(var==4){
      res=evaluate(Zm.pert, margToUnif(x,4))
    }else if(var==5){
      res=evaluate(L.pert, margToUnif(x,5))
    }else if(var==6){
      res=x
    }
    return(res)
  }
  
  pertMap<-function(X){
    res=X
    for(i in 1:6){
      res[,i]=margToPert(X[,i], i)
    }
    return(res)
  }
  
  sim_flood.pert<-function(n, pert){
    if(pert){
      return(pertMap(sim.RWL(n)))
    }else{
      return(sim.RWL(n))
    }
  }
  
  Y=model.RWL(sim_flood.pert(n_sim, T))
  res.stats[theta_idx,-1]=c(mean(Y), var(Y), quantile(Y,c(0,0.025, 0.05, 0.95,0.975,1), type=1))
  
  endRep=Sys.time()
  
  cat("Computation time:", difftime(endRep, startRep, units="mins"), "mins\n", sep=" ")
  
}


######################################
# Saving the results
save(res.pertVEs, res.shap, res.stats, L.pert, Q.pert, Zm.pert, res.pert.Ks, file="Chapter5/RWL/results/pertRes_RWL.RData")
