############################################################################
# Estimation the conditional elements using MC scheme

############################################
# Needed packages
library(parallel)
library(doParallel)
library(foreach)
library(gtools)


############################################
# Estimating a total Sobol' index by MC
estim.CE.MC<-function(condX, condSim, model, subset, Ni){
  ##########################
  # Pre-processing
  No=nrow(condX)
  d=ncol(condX)
  
  varVec=rep(0, No)
  for(i in 1:No){
    X_=condSim(n=Ni,
               Sj=subset,
               Sjc=setdiff(1:d, subset),
               xjc=condX[i,])
    Y_<-model(X_)
    varVec[i]=var(Y_)
  }
  
  #Returning an estimate of St
  St=mean(varVec)
  return(St)
}

############################################
# Estimating every total Sobol' indices

estim.CEs.MC<-function(model, jointSim, condSim, Nv, No, Ni,parl, cl, pert=F){
  #Initial Dataset
  X=jointSim(Nv)
  vy=var(model(X))
  d=ncol(X)
  
  condX=jointSim(No)
  
  ##########################
  #Pre-processing
  
  #List of subsets by order ([[2]]: subsets of order 1...etc)
  indices<-rep(list(0), d+1)
  for(j in 1:d){
    indices[[j+1]]<-t(gtools::combinations(n=d, r=j))
  }
  
  #List of conditional elements (same structure as the indices object)
  VEs<-rep(list(0), d+1)
  
  ####################################
  #Conditional elements estimation
  for(j in 1:(d-1)){
    ##################################################
    # MC Computation of the conditional elements
    if(!is.null(parl)){
      #############################
      #Parallelized VE estimation
      VEs[[j+1]]<-parallel::parApply(cl,indices[[j+1]], 2,estim.CE.MC,
                                     condX = condX, 
                                     condSim = condSim, 
                                     model = model, 
                                     Ni=Ni)/vy
    }else{
      #############################
      # VE estimation
      VEs[[j+1]]<-apply(indices[[j+1]], 2,estim.CE.MC,
                        condX = condX, 
                        condSim = condSim, 
                        model = model, 
                        Ni=Ni)/vy
    }
    VEs[[d+1]]=1
  }
  return(return(list("VEs"=VEs,
                     "indices"=indices,
                     "vy"=vy)))
}

#Returns the Total Sobol' indices from the closed Sobol' indices
get.SobolT<-function(VEs){
  d<-length(VEs)-1
  SobT<-rep(list(0), length(VEs))
  for (ord in 1:d){
    SobT[[ord+1]] = rev(VEs[[d+1]] - VEs[[d+1-ord]])
  }
  return(SobT)
}

########################################################
# Harsanyi Dividends from the conditional elements
harsaDivFromCEs<-function(CEs,dual=F, norm=F){
  # Retrieving data
  VEs=CEs$VEs
  indices=CEs$indices
  d=length(VEs)-1
  
  if(dual){
    VEs=get.SobolT(VEs)
  }
  
  harsa=VEs
  for(ord in 2:d){
    nb_ord<-ncol(indices[[ord+1]])
    res_ord<-rep(0, nb_ord)
    
    for(j in 1:nb_ord){
      indice<-indices[[ord+1]][,j]
      for(k in seq(ord,1,-1)){
        w<-(-1)^(ord-k)
        idx<-which(colSums(matrix(apply(indices[[k+1]], 2, function(x){x %in% indice}),ncol=ncol(indices[[k+1]])))==k)
        res_ord[[j]]=res_ord[[j]]+w*sum(harsa[[k+1]][idx])
      }
    }
    harsa[[ord+1]]=res_ord
  }
  return(harsa)
}

############################################
# Shapley aggregation
shFromCEs<-function(CEs){
    indices=CEs$indices
    VEs=CEs$VEs
    
    d=length(VEs)-1
    
    #Weights vector, each element is the weight for |A| ([1]:|A|=0, [2]:|A|=1, ...)
    comb_weights<-rep(0,d)
    for(j in 1:d){
      comb_weights[j]<-1/choose(n=(d-1), j-1)
    }
  
    ################################
    #Unparallelized weight affectation
    Shaps<-rep(0,d)
    for(var_j in 1:d){
      #For every order, the increments are weighted and summed
      #for the input var_j
      for(ord in 1:d){
        if(ord==1){
          idx<-which(indices[[ord+1]]==var_j)
          Shaps[var_j]=Shaps[var_j]+comb_weights[ord]*VEs[[ord+1]][as.vector(idx)]
        }else{
          #Columns of subsets of order ord+1 containing var_j
          idx_j<-which(indices[[ord+1]]==var_j, arr.ind=T)[,2]
          #Columns of subsets of order ord without j
          idx_woj<-which(apply(indices[[ord]]!=var_j, 2, all))
          #Total increment value for order ord
          tot_incr<-sum(VEs[[ord+1]][as.vector(idx_j)]) - sum(VEs[[ord]][as.vector(idx_woj)])
          Shaps[var_j]=Shaps[var_j]+comb_weights[ord]*tot_incr
        }
      }
    }
  
  Shaps=Shaps/d
  return(list("Shap"=Shaps,
              "VEs"=VEs,
              "indices"=indices))
}

############################################
# PME aggregation

recur.PV<-function(indices, VEs){
  d=length(indices)-1
  Ps<-rep(list(0), d)
  Ps[[1]]=VEs[[2]]
  for(ord in 2:d){
    Ws<-VEs[[ord+1]]
    nbset<-ncol(indices[[ord+1]])#Number of subsets S
    Ps[[ord]]<-rep(0, nbset)#Setting up the Ps results
    for(i in 1:nbset){
      S<-c(indices[[ord+1]][,i]) #For every possible subset S of order ord
      idx_Spi=which(colSums(matrix(indices[[ord]]%in%S, nrow=length(S)-1, ncol=ncol(indices[[ord]])))==length(S)-1)#Find every S\{i} for all i
      Ps[[ord]][i]=Ws[i]/sum(1/(Ps[[ord-1]][idx_Spi]))#Recursively compute Ps
    }
  }
  return(Ps)
}

pmeFromCEs<-function(CEs, norm=F){
  VEs=CEs$VEs
  indices=CEs$indices
  d=length(VEs)-1
  
  #Identifying inputs with 0 Total Sobol (zero players)
  idx.z<-indices[[2]][,which(VEs[[2]]==0)]
  
  if(any(idx.z)){
    #Checking if there are other zero coalitions in each order
    Z.coal.order<-which(sapply(VEs, function(x) any(x==0)))
    
    #Getting the cardinality of the biggest zero coalition
    Z.cardMax<-max(Z.coal.order)-1 #It should always be over or equal to 1 since there are zero players
    
    #Biggest zero coalitions of max cardinality
    Z.coal.cardMax<-indices[[Z.cardMax+1]][,which(VEs[[Z.cardMax+1]]==0), drop=F]
    
    #Boolean of size card(idx.z) indicating if they are in all the zero coalitions of max cardinality
    z.zeroPV<-sapply(idx.z, function(x) sum(colSums(Z.coal.cardMax==x))==1)
    
    #First case: Z.cardMax=d, no one gets anything (should not happen too often)
    if(Z.cardMax==d){
      PV=rep(0,d)
      PV=matrix(PV, ncol=1)
    }else if(Z.cardMax==d-1){
      #Second case: Z.cardMax=d-1, no need to loop
      PV=rep(0,d)
      idx.pv=setdiff(1:d,idx.z[z.zeroPV])
      PV[idx.pv] = VEs[[d+1]]/nrow(Z.coal.cardMax)
      PV=matrix(PV, ncol=1)
    }else{
      #Third case: Recursive computing
      PS.i<-rep(0,d)
      PS.N<-rep(0,d)
      #Setting up VEs and indices by removing every zero coalition
      for(idx.Zcoal in 1:ncol(Z.coal.cardMax)){
        Z.coal=Z.coal.cardMax[,idx.Zcoal]
        indices_<-rep(list(0), d-Z.cardMax+1)
        VEs_<-rep(list(0), d-Z.cardMax+1)
        for(i in 1:(d-Z.cardMax)){
          checkmat<-matrix(is.element(indices[[i+1]], Z.coal), i, ncol(indices[[i+1]]))
          idx_ind_null <- as.vector(which(colSums(checkmat) == 0))
          ind_tmp<-indices[[i+1]][, idx_ind_null, drop=F]
          indices_[[i+1]]<-rbind(ind_tmp, matrix(rep(Z.coal, ncol(ind_tmp)), ncol=ncol(ind_tmp), nrow=length(Z.coal)))
        }
        for(i in 1:(length(indices_)-1)){
          idx.get<-apply(indices_[[i+1]],
                         2,
                         function(x) which(colSums(matrix(indices[[Z.cardMax+i+1]] %in% x,
                                                          nrow=Z.cardMax+i, 
                                                          ncol=ncol(indices[[Z.cardMax+i+1]])))==i+Z.cardMax)
          )
          
          VEs_[[i+1]]<-VEs[[Z.cardMax+i+1]][idx.get]
        }
        PS<-recur.PV(indices=indices_,
                     VEs=VEs_)
        idx.var<-as.vector(indices_[[2]][1,])
        PS.i[idx.var]=PS.i[idx.var]+(1/rev(PS[[length(PS)-1]]))
        PS.N=PS.N+1/PS[[length(PS)]]
      }
      PV=matrix(PS.i/PS.N, ncol=1)
    }
    
  }else{
    PS<-recur.PV(indices=indices,
                 VEs=VEs)
    PV<-matrix(rev(PS[[d]]/PS[[d-1]]), ncol=1)
  }
  return(list("PME"=PV,
              "VEs"=VEs,
              "indices"=indices))
}

############################################################################
# Estimation the conditional elements using MC scheme for perturbations

############################################
# Estimating a total Sobol' index by MC
estim.CE.MC.pert<-function(condX, condSim, model, subset, Ni, pertMapping){
  ##########################
  # Pre-processing
  No=nrow(condX)
  d=ncol(condX)
  
  varVec=rep(0, No)
  for(i in 1:No){
    X_=condSim(n=Ni,
               Sj=subset,
               Sjc=setdiff(1:d, subset),
               xjc=condX[i,])
    X_pert=pertMapping(X_)
    Y_<-model(X_pert)
    varVec[i]=var(Y_)
  }
  
  #Returning an estimate of St
  St=mean(varVec)
  return(St)
}

############################################
# Estimating every total Sobol' indices

estim.CEs.MC.pert<-function(model, jointSim, condSim, Nv, No, Ni,parl, cl, pertMapping){
  #Initial Dataset
  X=pertMapping(jointSim(Nv))
  vy=var(model(X))
  d=ncol(X)
  
  condX=jointSim(No)
  
  ##########################
  #Pre-processing
  
  #List of subsets by order ([[2]]: subsets of order 1...etc)
  indices<-rep(list(0), d+1)
  for(j in 1:d){
    indices[[j+1]]<-t(gtools::combinations(n=d, r=j))
  }
  
  #List of conditional elements (same structure as the indices object)
  VEs<-rep(list(0), d+1)
  
  ####################################
  #Conditional elements estimation
  for(j in 1:(d-1)){
    ##################################################
    # MC Computation of the conditional elements
    if(!is.null(parl)){
      #############################
      #Parallelized VE estimation
      VEs[[j+1]]<-parallel::parApply(cl,indices[[j+1]], 2,estim.CE.MC.pert,
                                     condX = condX, 
                                     condSim = condSim, 
                                     model = model, 
                                     Ni=Ni,
                                     pertMapping=pertMapping)/vy
    }else{
      #############################
      # VE estimation
      VEs[[j+1]]<-apply(indices[[j+1]], 2,estim.CE.MC.pert,
                        condX = condX, 
                        condSim = condSim, 
                        model = model, 
                        Ni=Ni,
                        pertMapping=pertMapping)/vy
    }
    VEs[[d+1]]=1
  }
  return(return(list("VEs"=VEs,
                     "indices"=indices,
                     "vy"=vy)))
}


