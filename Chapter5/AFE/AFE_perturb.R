########################################
# Required Libraries

# Loading the "interpretability" package : Uncomment if first usage

#install.packages("Chapter5/package/interpretability_0.1.0.tar.gz", 
#                 repos=NULL, 
#                 type='source')

library(keras)
library(interpretability)
library(sensitivity)

#Loading the NN model & data
nn=load_model_hdf5("useCases/acousticFireExtinguisher/model/model_AFE.hdf5")
load(file="useCases/acousticFireExtinguisher/model/AFE_trainingData.RData")
load(file="useCases/acousticFireExtinguisher/data/AFE_data.RData")
fire.nn=data.frame(data.AFE)
fire=data.frame(fire)

#Prediction black-box model
nn.AFE<-function(X){
  a=nn %>% predict(X)
  res=apply(a,1, function(x) x[2]>x[1])
  return(res)
}

###################################
# Modelling quantile constraints

# We preserve the airflow 0.1 to 0.6 empirical quantiles, with 0.05 step
add.const<-cbind(seq(0.1,0.6,0.05),
                 c(quantile(fire.nn[,4], seq(0.1,0.6,0.05), type=1)))

###################################
# Computing the perturbations

# We perturb the 0.8-quantile in eta = [9.5, 14.5], initially at 12
fire.rob<-robustPert(X=as.matrix(fire.nn)[,-6],
                     model=nn.AFE,
                     pert.var=4,
                     theta.prec=0.01,
                     pert.range=c(9.5, 14.5), 
                     init.quant=c(0.8, quantile(fire.nn[,4], 0.8, type=1)),
                     add.const=add.const,
                     OmegaX=c(0,17),
                     pertType="shift",
                     method=c("polynomial",9),
                     par=8, #Computing in parallel on 8 cores
                     standard=F,
                     normalize=T)

##################################
# Global Robustness

#Seed for reproducibility
set.seed(12345)

# Shapley effects with KNN - For each theta
res.shap<-matrix(NA, ncol=ncol(fire), nrow=nrow(fire.rob$theta))

X_=fire[,-7]

for(i in 1:length(fire.rob$pert.data)){
  
  print(paste(i," / ", length(fire.rob$pert.data), sep=""))
  startRep=Sys.time()
  
  
  data.pert=fire.rob$pert.data[[i]]
  Y=data.pert[,8]
  X_[,5]=data.pert[,4]
  X_=data.frame(X_)
  
  res=shapleysobol_knn(model=NULL,
                       X=X_,
                       parl=6,
                       n.knn=6)
  res=tell(res, y=Y)
  
  res.shap[i,]=c(fire.rob$theta[i,1], res$Shap)
  
  
  endRep=Sys.time()
  
  print(paste("Computation Time : ", 
              round(as.numeric(difftime(endRep, startRep, units="secs")),3),
              " seconds."))
  print("-------------------------------------------------------")
}
colnames(res.shap)=c("theta", colnames(X_))


#Computing the baseline shapley effects (no perturbations)
baseline.shap<-shapleysobol_knn(model=NULL,
                                X=fire[,-7],
                                parl=8,
                                n.knn=6)
baseline.shap=tell(baseline.shap, y=fire$Y)

y.pert=fire.rob$pred.pert[,-1] #Predictions on perturbed data
x.dist=fire.rob$dist.pert[,-1] # Magnitude (l2) of the perturbation

for(i in 1:nrow(x.dist)){
  x.dist[i,]=rowSums(fire.rob$pert.data[[i]][,-8] -fire[,-6] )
}

aa=t(apply(y.pert, 1, function(x) x==fire.rob$Y.init))


x.PCP=as.vector(x.dist[aa]) #No change in prediction
x.CP=as.vector(x.dist[!aa]) #Prediction change

save(fire.rob, res.shap, x.PCP, x.CP, baseline.shap, fire.nn,
     file="Chapter5/AFE/results/AFE_pertRes.RData")
