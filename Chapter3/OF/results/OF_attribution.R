#########################################################################
## Shapley effects and PME computation for the optical filters Use-Case
#########################################################################

#################
# Packages
library(sensitivity)

#################
# Load the Data

load(file="useCases/opticalFilter/OF_data.RData")
X<- data.OptiFilter[,-14]
Y<- data.OptiFilter$Y

#####################
# Indices estimation

#Seed for reproducibility
set.seed(123456)

#Shapley Effects
Shap.OF <- shapleysobol_knn(model=NULL, 
                                    X = X, 
                                    n.knn = 6, 
                                    parl = 8, 
                                    conf=0.9)
Shap.OF<-tell(Shap.OF, Y)

#PME
PME.OF <- pme_knn(model=NULL, 
                          X = X, 
                          n.knn = 6,
                          parl = 8, 
                          conf=0.9)
PME.OF<-tell(PME.OF, Y)

save(Shap.OF, PME.OF, file="Chapter3/OF/results/allocations_OF.RData")

#######################
# Surrogate Modelling

#################
# Packages
library(DiceKriging)

GP1.vars<-Shap.OF$Shap>=0.05 #All inputs
GP2.vars<-PME.OF$PME>=0.05
GP3.vars<-PME.OF$PME>=0.022

#################
#Reproducibility
set.seed(12345)


#################
#GP1
GP1<-km(design=X[,GP1.vars], response=Y, nugget.estim=T)
GP1.loo <- leaveOneOut.km(GP1, type="UK")
GP1.Q2 <- 1 - mean((GP1.loo$mean - Y)^2) / var(Y)
print(paste("Q² (GP1) =", GP1.Q2))

#################
#GP2
GP2<-km(design=X[,GP2.vars], response=Y, nugget.estim=T)
GP2.loo <- leaveOneOut.km(GP2, type="UK")
GP2.Q2 <- 1 - mean((GP2.loo$mean - Y)^2) / var(Y)
print(paste("Q² (GP2) =", GP2.Q2))

#################
#GP3
GP3<-km(design=X[,GP3.vars], response=Y, nugget.estim=T)
GP3.loo <- leaveOneOut.km(GP3, type="UK")
GP3.Q2 <- 1 - mean((GP3.loo$mean - Y)^2) / var(Y)
print(paste("Q² (GP3) =", GP3.Q2))

save(X,Y, GP1, GP1.loo,
     GP2, GP2.loo,
     GP3, GP3.loo, file="Chapter3/OF/results/models_OF.RData")