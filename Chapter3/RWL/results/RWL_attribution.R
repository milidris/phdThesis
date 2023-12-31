#######################################################
#Importing functions to simulate and estimate via MC
source(file="useCases/riverWaterLevel/RWL_model.R")
source(file="Chapter3/functions/allocationsMCestimation.R")

#######################################################
# Computing conditional elements for RWL

#Preparing parallel cluster
cl=parallel::makeCluster(8)
clusterExport(cl, varlist=list("condGausCopSim",
                               "gausCopSim",
                               "margToUnif",
                               "model.RWL",
                               "sim.RWL",
                               "simCond.RWL",
                               "unifToMarg",
                               "estim.CE.MC",
                               "estim.CEs.MC")
)
clusterEvalQ(cl, library(sensitivity))
clusterEvalQ(cl, library(triangle))
clusterEvalQ(cl, library(VineCopula))
clusterEvalQ(cl, library(evd))
doParallel::registerDoParallel(cl)

nb_rep=150

model.flood<-function(X){
  Y<-model.RWL(X)
  Y.flood<-as.numeric(Y>54.25)
  return(Y.flood)
}

#Seed for reproducibility
set.seed(12345)

CEs.RWL=rep(list(0), nb_rep)
CEs.flood=rep(list(0), nb_rep)
for(i in 1:nb_rep){
  cat("-------------------------------------------\n")
  cat("Iteration: ",i, "\n", sep="")
  deb=Sys.time()
  CE.RWL=estim.CEs.MC(model=model.RWL, 
                      jointSim=sim.RWL, 
                      condSim=simCond.RWL,
                      Nv=3e5, 
                      No=2e3,
                      Ni=300,
                      parl=T,
                      cl=cl)
  fin=Sys.time()
  cat("RWL iteration:", difftime(fin, deb, units="mins"), "mins\n", sep=" ")
  CEs.RWL[[i]]=CE.RWL
  
  
  deb=Sys.time()
  CE.flood=estim.CEs.MC(model=model.flood, 
                        jointSim=sim.RWL, 
                        condSim=simCond.RWL,
                        Nv=3e5, 
                        No=2e3,
                        Ni=300,
                        parl=T,
                        cl=cl)
  fin=Sys.time()
  cat("Flood iteration:", difftime(fin, deb, units="mins"), "mins\n", sep=" ")
  CEs.flood[[i]]=CE.flood
}

stopCluster(cl)

save(CEs.RWL, file="Chapter3/RWL/results/CE_RWL.RData")
save(CEs.flood, file="Chapter3/RWL/results/CE_flood.RData")


#######################################################
# Computing attribution

load(file="Chapter3/RWL/results/CE_flood.RData")
load(file="Chapter3/RWL/results/CE_RWL.RData")

Shaps.RWL=matrix(ncol=6, nrow=nb_rep)
Shaps.flood=matrix(ncol=6, nrow=nb_rep)

PMEs.RWL=matrix(ncol=6, nrow=nb_rep)
PMEs.flood=matrix(ncol=6, nrow=nb_rep)

for(i in 1:nb_rep){
  
  Shaps.RWL[i,]=shFromCEs(CEs.RWL[[i]])$Shap
  Shaps.flood[i,]=shFromCEs(CEs.flood[[i]])$Shap
  
  PMEs.RWL[i,]=pmeFromCEs(CEs.RWL[[i]])$PME
  PMEs.flood[i,]=pmeFromCEs(CEs.flood[[i]])$PME
}

save(Shaps.RWL, 
     Shaps.flood,
     PMEs.RWL,
     PMEs.flood,
     file="Chapter3/RWL/results/allocations_RWL.RData")

