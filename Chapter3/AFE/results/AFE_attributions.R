########################################
# Required Libraries
library(keras)
library(interpretability)
library(sensitivity)
source(file="Chapter3/functions/allocationsMCestimation.R")

#Loading the NN model & data
model.AFE=load_model_hdf5("useCases/acousticFireExtinguisher/model/model_AFE.hdf5")
load(file="useCases/acousticFireExtinguisher/model/AFE_trainingData.RData")
load(file="useCases/acousticFireExtinguisher/data/AFE_data.RData")

#Prediction black-box model
nn.AFE<-function(X){
  a=model.AFE %>% predict(X)
  res=apply(a,1, function(x) x[2]>x[1])
  return(res)
}

#Vector of predictions
Y<-nn.AFE(as.matrix(data.AFE[,-6]))
X<-as.data.frame(fire[,-7])


#Seed for reproducibility
set.seed(123456)

# Computing the Shapley effects by K-NN using sensitivity
PMEs<-pme_knn(model=NULL,
                        X=X,
                        parl=8,
                        n.knn=6)
PMEs=tell(PMEs, y=Y)

CE.AFE<-list("VEs"=PMEs$VE,
               "indices"=PMEs$indices)

save(CE.AFE, file="Chapter3/AFE/results/CE_AFE.RData")

#######################################################
# Computing attribution

PME.AFE=pmeFromCEs(CE.AFE)
Shap.AFE=shFromCEs(CE.AFE)

save(Shap.AFE, 
     PME.AFE,
     file="Chapter3/AFE/results/allocations_AFE.RData")
