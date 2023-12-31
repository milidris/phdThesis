########################################
# Packages
library(triangle)
library(VineCopula)
library(evd)
library(sensitivity)
library(keras)
library(interpretability)

source(file="useCases/riverWaterLevel/RWL_model.R")

#################################
## Generating Dataset

# Reproducibility
set.seed(12345)

X.train = sim.RWL(1e6)
Y.train<-model.RWL(X.train)

X.val = sim.RWL(1e5)
Y.val<-model.RWL(X.val)

######################################################################
## NN Model

#Reproducibility
set.seed(12345)

#########################
# NN-model architecture

input <- layer_input(shape=6) #6 features as inputs
output <- input %>%
  layer_dense(units = 64, activation = "relu") %>% 
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 1)

nn <- keras_model(input, output)

nn %>%
  compile(
    loss = "mean_squared_error",
    optimizer = optimizer_adam(lr = 0.001), #Learning rate setup
    metrics = c("mean_squared_error")
  )

early_stop <- callback_early_stopping(monitor = "val_loss", patience = 50)
# Stop the training if no improvement on the validation loss for 50 epochs

# Model training
set.seed(54321)
nn %>% fit(
  x = X.train,
  y = Y.train,
  batch_size=1024,
  epochs = 1000,
  verbose = 2,
  validation_data=list(X.val, Y.val),
  callbacks=list(early_stop)
)

# Saving the model
save_model_hdf5(nn, "Chapter5/RWL/surrogate/nn_RWL.hdf5")

############################################################
# Validation

nn=load_model_hdf5("Chapter5/RWL/surrogate/nn_RWL.hdf5")

Y.train.pred<-nn%>%predict(as.matrix(X.train))
Y.val.pred<-nn%>%predict(as.matrix(X.val))

###########################################################
## R^2

1-mean((Y.train.pred - Y.train)^2)/var(Y.train)
1-mean((Y.val.pred - Y.val)^2)/var(Y.val)

# Mean, Var and Quantiles of the output
n_sim=4e5
theta=seq(-1,1,0.01)
res.stats.nn<-matrix(NA, ncol=9, nrow=length(res.pert.Ks))
colnames(res.stats.nn)=c("theta","mean","var","min","q025","q05", "q95","q975","max")

#Loading data and models
load(file="Chapter5/RWL/results/pertRes_RWL.RData")
nn=load_model_hdf5("Chapter5/RWL/surrogate/nn_RWL.hdf5")

nn.predict<-function(X,scores=F){
  a=nn %>% predict(X)
  if(scores){
    res<-a[,2]
  }else{
    res=apply(a,1, function(x) x[2]>x[1])
  }
  return(res)
}

#res.stats[,1]=theta
res.stats.nn[,1]=theta

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
  
  X_pert = sim_flood.pert(n_sim, T)
  
  Y.nn=nn%>%predict(as.matrix(X_pert))
  
  res.stats.nn[theta_idx,-1]=c(mean(Y.nn), var(Y.nn), quantile(Y.nn,c(0,0.025, 0.05, 0.95,0.975,1), type=1))
  
  endRep=Sys.time()
  
  cat("Computation time:", difftime(endRep, startRep, units="mins"), "mins\n", sep=" ")
}


######################################
# Saving the results
save(res.stats.nn, file="Chapter5/RWL/surrogate/pertSurrogate_RWL.RData")
