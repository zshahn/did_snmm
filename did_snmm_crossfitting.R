library(nleqslv)

#Generate Data
inv_logit = function(x){
  exp(x)/(exp(x)+1)
}

#true blip function used to generate data
blip = function(m,k,l,a,psi){
  a*sum(c(k-m+1,l[m])*psi)
}

#number of time-points including time 0 when treatment is always 0
ntimes=5
#simulate data for one patient 
#ordering L,A,Y at each timepoint (it's L,Y,A in the paper, but L,A,Y in the real applications)
#there's a baseline confounder that influences outcome and treatment at each time
#the psi argument here is the true value of the parameter of the blip function we're trying to estimate, so (1,1) 
sim_pat_coarse_ptT = function(times=ntimes,blip_func=blip,psi=c(1,1)){
  L = rep(NA,times)
  Y = rep(NA,times)
  Y_0 = rep(NA,times)
  A = rep(0,times)
  start=Inf
  U = rnorm(1)
  L[1] = rnorm(1)
  Y[1] = rnorm(1,U+L[1])
  Y_0[1] = Y[1]
  #everyone starts out at 0
  A[1]=0
  #make counterfactual untreated trajectory
  for(t in 2:(times)){
    L[t] = rnorm(1,L[t-1])
    Y_0[t] = rnorm(1,mean(L[1:t])+U)
    A[t] = rbinom(1,1,inv_logit(-1+L[t]+U))
  }
  start = min(which(A!=0))
  #if no treatment, just return untreated trajectory
  if(start==Inf){
    return(list(A=A,Y=Y_0,L=L))
  }else{
    #otherwise blip up the outcomes using the blip function
    A[start:times]=A[start]
    Y= Y_0
    for(k in start:times){
      Y[k] = rnorm(1,Y_0[k] + blip_func(start,k,L,A[start],psi))
    }
  }
  return(list(A=A,Y=Y,L=L))
}
sim_pat_coarse_ptT()
#number of patients to simulate
N = 1000000
#make a dataframe of all the observed values for each patient
#id,A,Y,L,time
data = replicate(N,sim_pat_coarse_ptT())
data = data.frame(cbind(A = unlist(data['A',]),Y=unlist(data['Y',]),L=unlist(data['L',])))
data$id = rep(1:N,each=ntimes)
data$time = rep(0:(ntimes-1),N)

#Add each subject's full history in each row (wide format) to facilitate some later computations
for(t in 0:max(data$time)){
  Lt = data[data$time==t,c('L','id')]
  names(Lt)[1] = paste0('L',t)
  data = merge(data,Lt)
  
  At = data[data$time==t,c('A','id')]
  names(At)[1] = paste0('A',t)
  data = merge(data,At)
  
  Yt = data[data$time==t,c('Y','id')]
  names(Yt)[1] = paste0('Y',t)
  data = merge(data,Yt)
}


#do ML crossfitting

#add time of first treatment for each subject, Inf if never treated
first_treat = aggregate(data$A,by=list(data$id),function(x)min(which(x!=0))-1)
names(first_treat) = c('id','first_treat_time')
data = merge(data,first_treat[,c('id','first_treat_time')])

#do sample splitting
sample1 = sample(N,ceiling(N/2),replace=F)
sample2 = setdiff(1:N,sample1)
data1 = data[data$id %in% sample1,]
data2 = data[data$id %in% sample2,]


#initialize lists to store psi_hat estimates
#as in remark 2, we will compute estimates of variation independent blip functions for each time point
#in the true data generating process, the blip functions are the same at each time (the blip() function in line 9) and have parameter (1,1)
psi_hats1 = vector(mode='list',length=max(data$time))
psi_hats2 = vector(mode='list',length=max(data$time))

#first, estimate the blip function for treatment at the final time-point
#make estimation (data_max) and training (data_max_train) datasets with just subjects not treated prior to the final timepoint
data_max = data2[data2$first_treat_time>=max(data$time) & data2$time==max(data$time),]
data_max_train = data1[data1$first_treat_time>=max(data$time) & data1$time==max(data$time),]

#fit a treatment model to the training dataset
if(max(data$time)>=2){
  vars = c(paste0('A',0:(max(data$time)-1)),paste0('L',0:max(data$time)),paste0('Y',0:(max(data$time)-2)))
}else{
  vars = c(paste0('A',0:(max(data$time)-1)),paste0('L',0:max(data$time)))
}
treat_mod_cv = xgb.cv(data = as.matrix(data_max_train[,vars]),label = data_max_train$A,nrounds = 100, nfold = 5, showsd = T, stratified = T, early_stopping_rounds = 10, maximize = F)
treat_mod = xgboost(data = as.matrix(data_max_train[,vars]),label = data_max_train$A,nrounds = treat_mod_cv$best_iteration)

#add the predictions from the treatment model fit to the training data to the estimation data
data_max$E_A = predict(treat_mod,newdata = as.matrix(data_max[,vars]),type='response')

#Following did_snmm_reboot.pdf, we fit a regression model for the first conditional expectation on the RHS of the equation in Remark 2.
data_max_train$pseudo_Y = data_max_train[,paste0('Y',max(data$time))] - data_max_train[,paste0('Y',max(data$time)-1)]
if(max(data$time)>=2){
  vars2 = c(paste0('A',0:(max(data$time)-1)),paste0('L',0:max(data$time)),paste0('Y',0:(max(data$time)-2)))
}else{
  vars2 = c(paste0('A',0:(max(data$time)-1)),paste0('L',0:max(data$time)))
}# pseudo_Y_mod = ranger(pseudo_Y~.,data = data_max_train[,c('pseudo_Y',vars2)])
pseudo_Y_mod_cv = xgb.cv(data = as.matrix(data_max_train[,vars2]),label = data_max_train$pseudo_Y,nrounds = 100, nfold = 5, showsd = T, stratified = T, early_stopping_rounds = 10, maximize = F)
pseudo_Y_mod = xgboost(data = as.matrix(data_max_train[,vars2]),label = data_max_train$pseudo_Y,nrounds = pseudo_Y_mod_cv$best_iteration)

#we add the predictions from the model for the first conditional expectation on the RHS of the equation in Remark 2 to the estimation dataset
data_max$E_pseudo_Y = predict(pseudo_Y_mod,newdata = as.matrix(data_max[,vars2]),type='response')

#a function to compute the value of the estimating equations for any given candidate value of psi in the estimation data set
compute_est_eq_psi_last = function(psi,H_ta=data_max){
  H_ta$gamma_ta = H_ta$A*psi[1] + H_ta$A*psi[2]*H_ta$L #the blip function
  H_ta$H = H_ta$Y - H_ta$gamma_ta #the blipped down outcomes
  H_ta$lag_H = H_ta[,paste0('Y',max(data$time)-1)] #the lagged blipped down outcomes
  H_ta$H_diffs = H_ta$H - H_ta$lag_H
  #binary treatment
  temp = (H_ta$H_diffs-H_ta$E_pseudo_Y+psi[1]*H_ta$E_A+psi[2]*H_ta[,'L']*H_ta$E_A)*cbind(1,H_ta$L)*(H_ta$A-H_ta$E_A)
  apply(temp,2,sum)
}
#solve estimating equation and store estimate of psi_mk for the last timepoint m
psi_hats2[[max(data$time)]]=nleqslv(x=rep(0,2),fn=compute_est_eq_psi_last,H_ta=data_max)$x


library(tidyverse)
get_lag = function(x){
  c(NA,x[1:(length(x)-1)])
}

#now, backwards recursively estimate the blip function parameters at each earlier timepoint, using previously computed estimates
#from later timepoints as needed
for(m in (max(data$time)-1):1){
  #make training and estimation data sets of subjects not treated prior to time m
  data_m = data2[data2$first_treat_time>=m & data2$time==m,]
  data_m_train = data1[data1$first_treat_time>=m & data1$time==m,]
  
  #fit treatment model in training dataset
  if(m>=2){
    vars = c(paste0('A',0:(m-1)),paste0('L',0:m),paste0('Y',0:(m-2)))
  }else{
    vars = c(paste0('A',0:(m-1)),paste0('L',0:m))
  }
  treat_mod_cv = xgb.cv(data = as.matrix(data_m_train[,vars]),label = data_m_train$A,nrounds = 100, nfold = 5, showsd = T, stratified = T, early_stopping_rounds = 10, maximize = F)
  treat_mod = xgboost(data = as.matrix(data_m_train[,vars]),label = data_m_train$A,nrounds = treat_mod_cv$best_iteration)
  
  #add treatment model prediction to estimation dataset
  data_m$E_A = predict(treat_mod,newdata = as.matrix(data_m[,vars]),type='response')
  
  #for each k>=m, use the training dataset to estimate the regression for the first conditional expectation on the RHS of the equation from Remark 2 and store its
  #prediction in the estimation dataset
  #also, compute blipped down outcomes in the estimation dataset for all subjects treated after m using the previously obtained blip function estimates for the times when they were treated
  for(k in max(data$time):m){
    data_m_train$pseudo_Y = NA
    data_m_train$pseudo_Y[data_m_train$first_treat_time>k|data_m_train$first_treat_time==m] = 
      data_m_train[data_m_train$first_treat_time>k|data_m_train$first_treat_time==m,paste0('Y',k)] - 
      data_m_train[data_m_train$first_treat_time>k|data_m_train$first_treat_time==m,paste0('Y',k-1)]
    if(k>m){
      for(j in k:(m+1)){
        data_m_train$pseudo_Y[data_m_train$first_treat_time==j] = 
          data_m_train[data_m_train$first_treat_time==j,paste0('Y',k)] - psi_hats2[[j]][1]*(k-j+1) - psi_hats2[[j]][2]*data_m_train[data_m_train$first_treat_time==j,paste0('L',j)] - 
          (data_m_train[data_m_train$first_treat_time==j,paste0('Y',k-1)] - ifelse((k-1)>=j,psi_hats2[[j]][1]*(k-1-j+1) + psi_hats2[[j]][2]*data_m_train[data_m_train$first_treat_time==j,paste0('L',j)],0))
      }
    }
    if(m>=2){
      vars2 = c(paste0('A',0:(m-1)),paste0('L',0:m),paste0('Y',0:(m-2)))
    }else{
      vars2 = c(paste0('A',0:(m-1)),paste0('L',0:m))
    }    
    pseudo_Y_mod_cv = xgb.cv(data = as.matrix(data_m_train[,vars2]),label = data_m_train$pseudo_Y,nrounds = 100, nfold = 5, showsd = T, stratified = T, early_stopping_rounds = 10, maximize = F)
    pseudo_Y_mod = xgboost(data = as.matrix(data_m_train[,vars2]),label = data_m_train$pseudo_Y,nrounds = pseudo_Y_mod_cv$best_iteration)
    
    data_m[,paste0('E_pseudo_Y',k)] = predict(pseudo_Y_mod,newdata = as.matrix(data_m[,vars2]),type='response')
    data_m[,paste0('H',k)] = NA
    data_m[data_m$first_treat_time>k|data_m$first_treat_time==m,paste0('H',k)] = data_m[data_m$first_treat_time>k|data_m$first_treat_time==m,paste0('Y',k)] 
    if(k>m){
      for(j in k:(m+1)){
        data_m[data_m$first_treat_time==j,paste0('H',k)] = data_m[data_m$first_treat_time==j,paste0('Y',k)] - 
          psi_hats2[[j]][1]*(k-j+1) - psi_hats2[[j]][2]*data_m[data_m$first_treat_time==j,paste0('L',j)] 
      }
    }
  }
    
    #construct a matrix with all the information that will be needed to compute the estimating equations for the parameters of the blip function for effects of treatment at time m on outcomes at later times
    Hm = cbind(id=rep(1:N,each=max(data$time)-m+2),k=rep((m-1):max(data$time),N))
    
    Hm = merge(Hm,data_m)
    Hm$H = NA
    for(j in m:max(data$time)){
      Hm$H[Hm$k==j]=Hm[Hm$k==j,paste0('H',j)]
    }
    Hm$H[Hm$k==(m-1)]=Hm[Hm$k==(m-1),paste0('Y',m-1)]
    
    Hm$E_pseudo_Y = NA
    for(j in m:max(data$time)){
      Hm$E_pseudo_Y[Hm$k==j]=Hm[Hm$k==j,paste0('E_pseudo_Y',j)]
    }

    Hm = Hm[order(Hm$id,Hm$k),]
    #function to compute estimating equations for given candidate value psi
    compute_est_eq_psi_last = function(psi,H_ta=Hm){
      H_ta$gamma = ifelse(H_ta$first_treat_time==m & H_ta$k>=m,psi[1]*(H_ta$k-m+1) + psi[2]*H_ta$L,0) 
      H_ta$H_g = H_ta$H - H_ta$gamma
      H_ta$H_g_lag = unlist(lapply(split(H_ta$H_g,H_ta$id),get_lag))
      H_ta$H_diffs = H_ta$H_g - H_ta$H_g_lag
      H_ta = H_ta[H_ta$k>=m,]
      #this code is currently for a binary treatment
      temp = (H_ta$H_diffs-H_ta$E_pseudo_Y+psi[1]*(H_ta$k-m+1)*H_ta$E_A+psi[2]*H_ta$L*H_ta$E_A)*cbind(1,H_ta$L)*(H_ta$A-H_ta$E_A)
      apply(temp,2,sum)
    }
    #solve estimating equations and store
  psi_hats2[[m]]=nleqslv(x=rep(0,2),fn=compute_est_eq_psi_last)$x
}
  
#check that the estimated psis are very close to the real values
psi_hats2




