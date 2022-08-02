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
#ordering L,A,Y at each timepoint (it's L,Y,A in the paper, but L,A,Y in the real application)
#there's a baseline confounder that influences outcome and treatment at each time
#(interestingly, seems like if U impacts the time-varying covariates and the time-varying covariates have lagged impacts on the outcome then you can't have parallel trends)
#the psi argument here is the true value of the parameter of the blip function we're trying to estimate, so 1,1 
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
N = 10000
#make a dataframe of all the observed values for each patient
#id,A,Y,L,time
data = replicate(N,sim_pat_coarse_ptT())
data = data.frame(cbind(A = unlist(data['A',]),Y=unlist(data['Y',]),L=unlist(data['L',])))
data$id = rep(1:N,each=ntimes)
data$time = rep(0:(ntimes-1),N)


#Transform standard dataframe into the form we need for our analysis
#one row for each id/m/k triplet where m and k are years with k>=m-1 
#ultimately the important columns will be: 
#time of first treatment
#covariates for blip model at year of first treatment for id if after m (R_*) (L_gamma)
#covariates at year m for blipped down outcome model for time k (D_*)
#estimating equation filler (Q_*)
#estimated treatment probability at time m (E_A)
#treatment value at time m (A_m)
#outcome values at years k and k-1 (Y and lag_Y)
#The below code creates and populates all the above columns, plus some possibly unnecessary ones

mk = expand.grid(m=1:(ntimes-1),k=0:(ntimes-1))
mk = mk[mk$k>=(mk$m-1),]
mk = mk[order(mk$m),]
Hmk = cbind(id=rep(1:N,each=nrow(mk)),m=rep(mk$m,N),k=rep(mk$k,N))
first_treat = aggregate(data$A,by=list(data$id),function(x)min(which(x!=0))-1)
names(first_treat) = c('id','first_treat_time')
Hmk = merge(Hmk,first_treat)

data$first_treat_time = data$time
data$m = data$time
data$k = data$time
data$lag_Y = unlist(lapply(split(data$Y,data$id),function(x)c(NA,x[1:(length(x)-1)])))

Hmk = merge(Hmk,data[,c('id','first_treat_time','A')],all.x = T)
Hmk$A[is.na(Hmk$A)] = 0
data$A_m = data$A
Hmk = merge(Hmk,data[,c('id','m','A_m')],all.x = T)
data$L_m = data$L
Hmk = merge(Hmk,data[,c('id','m','L_m')],all.x = T)


Hmk = Hmk[order(Hmk$id,Hmk$m,Hmk$k),]

Hmk = merge(Hmk,data[,c("id","L","first_treat_time")],all.x = T)
Hmk = merge(Hmk,data[,c("id","Y","k")],all.x = T)
Hmk = merge(Hmk,data[,c("id","lag_Y","k")],all.x = T)


#Fit treatment model
data$past_A = unlist(lapply(split(data$A,data$id),function(x)cumsum(cumsum(x!=0))>1))
A_mod = glm(A~L,data=data[!data$past_A&data$time>0,],family="binomial")

#add estimated treatment expected values to Hmk
data$E_A = NA
data$E_A[data$time>0] = predict(A_mod,newdata = data[data$time>0,],type='response')

#add running mean L covariate to be used in outcome model
cummean = function(x){
  cumsum(x)/(1:length(x))
}
data$L_mean_m = unlist(lapply(split(data$L,data$id),cummean))
Hmk = merge(Hmk,data[,c("id","m",'E_A')],all.x=T)
Hmk = merge(Hmk,data[,c("id","m",'L_mean_m')],all.x=T)

#add treatment and covariate variables to Hmk
data$L_gamma = data$L
Hmk$first_treat_time2 = Hmk$first_treat_time
first_treat$first_treat_time2 = first_treat$first_treat_time
data = merge(data,first_treat[,c('id','first_treat_time2')])
Hmk = merge(Hmk,data[data$time==data$first_treat_time2,c("id","first_treat_time2",'L_gamma')],all.x=T)

Hmk = Hmk[order(Hmk$id,Hmk$m,Hmk$k),]

#Restrict to rows that actually contribute to estimating equations
est_mat = Hmk[Hmk$k>=Hmk$m & Hmk$m<=Hmk$first_treat_time,]

#These estimating equation filler columns are specific to the blip function. 
#because the blip function has 2 dimensional parameter psi, Q_ has 2 dimensions
#Q needs to only be a function of variables knowable at time m 
est_mat$Q_km = est_mat$k - est_mat$m + 1
est_mat$Q_L = est_mat$L_m

#also blip function specific, but has blip function covariates at time start_time instead of m
est_mat$R_km = est_mat$k - est_mat$first_treat_time + 1
est_mat$R_L = est_mat$L_gamma

#i did a flexible outcome model with a coefficient for each mk pair, L_m, and running mean L_mean_m
for(i in 1:(ntimes-1)){
  for(j in i:(ntimes-1)){
    est_mat[,paste0('D_',i,'_',j)] = as.numeric(est_mat$m==i & est_mat$k==j)
  }
}
est_mat$D_Lm = est_mat$L_m
est_mat$D_Lmean = est_mat$L_mean_m

#Now start computing the terms in the closed form linear scenario solution
#first term in closed form linear solution
# A = colSums((est_mat$k_hpi-est_mat$lag_k_hpi)*(est_mat$A_any_m-est_mat$p_any)*est_mat[,grep('Q_',names(est_mat))])
A_dr = colSums((est_mat$Y-est_mat$lag_Y)*cbind((est_mat$A_m-est_mat$E_A)*est_mat[,grep('Q_',names(est_mat))],est_mat[,grep('D_',names(est_mat))]))
#grep('Q_',names(est_mat)) gives filler/garbage columns
#grep('D_',names(est_mat)) gives nuisance outcome model covariate columns

#make second term of closed form linear solution formula
R_inds=grep('R_',names(est_mat))
D_inds=grep('D_',names(est_mat))
Q_inds=grep('Q_',names(est_mat))
p = length(Q_inds)
#Vimk in paper (17)
get_S = function(ind){
  if(est_mat$first_treat_time[ind]>est_mat$k[ind]){
    return(rep(0,p))
  }else{
    return(as.vector(unlist(est_mat[ind,R_inds])))
  }
}

#Vimk-1 in paper (17)
get_S_lag = function(ind){
  if(est_mat$first_treat_time[ind]>(est_mat$k[ind]-1)){
    return(rep(0,p))
  }else{
    return(as.vector(unlist(est_mat[ind-1,R_inds])))
  }
}

get_D = function(ind){
  as.vector(unlist(est_mat[ind,D_inds]))
}

get_Q = function(ind){
  as.vector(unlist(est_mat[ind,Q_inds]))
}

S_diff_list = lapply(1:nrow(est_mat),function(i) get_S(i)-get_S_lag(i))
D_list = lapply(1:nrow(est_mat),get_D)
Q_list = lapply(1:nrow(est_mat),get_Q)

get_B_dr = function(sdiff,d,q,x){
  unlist(c(sdiff,d)) %*% t(unlist(c(q*x,d)))
}
get_B = function(sdiff,q,x){
  unlist(sdiff) %*% t(unlist(q))*x
}
B_list_dr = vector(mode='list',length=length(S_diff_list))
for(i in 1:length(S_diff_list)){
  B_list_dr[[i]] = get_B_dr(S_diff_list[[i]],D_list[[i]],Q_list[[i]],est_mat$A_m[i]-est_mat$E_A[i])
}
B_dr = Reduce('+',B_list_dr)
#this is the estimate of the blip function parameter
#do (17) from paper
psi_hat_dr = as.vector(A_dr%*%solve(B_dr))
#first two elements of psi_hat_dr are the estimates of psi, seems to work correctly and give about 1,1











#do ML crossfitting
library(ranger)
mk = expand.grid(m=1:(ntimes-1),k=0:(ntimes-1))
mk = mk[mk$k>=(mk$m-1),]
mk = mk[order(mk$m),]
Hmk = cbind(id=rep(1:N,each=nrow(mk)),m=rep(mk$m,N),k=rep(mk$k,N))
first_treat = aggregate(data$A,by=list(data$id),function(x)min(which(x!=0))-1)
names(first_treat) = c('id','first_treat_time')
Hmk = merge(Hmk,first_treat)

data$first_treat_time = data$time
data$m = data$time
data$k = data$time
data$lag_Y = unlist(lapply(split(data$Y,data$id),function(x)c(NA,x[1:(length(x)-1)])))

Hmk = merge(Hmk,data[,c('id','first_treat_time','A')],all.x = T)
Hmk$A[is.na(Hmk$A)] = 0
data$A_m = data$A
Hmk = merge(Hmk,data[,c('id','m','A_m')],all.x = T)
data$L_m = data$L
Hmk = merge(Hmk,data[,c('id','m','L_m')],all.x = T)

Hmk = Hmk[order(Hmk$id,Hmk$m,Hmk$k),]

Hmk = merge(Hmk,data[,c("id","L","first_treat_time")],all.x = T)
Hmk = merge(Hmk,data[,c("id","Y","k")],all.x = T)
Hmk = merge(Hmk,data[,c("id","lag_Y","k")],all.x = T)

#add running mean L covariate to be used in outcome model
cummean = function(x){
  cumsum(x)/(1:length(x))
}
data$L_mean_m = unlist(lapply(split(data$L,data$id),cummean))

Hmk = merge(Hmk,data[,c("id","m",'L_mean_m')],all.x=T)

#add treatment and covariate variables to Hmk
data$L_gamma = data$L
Hmk$first_treat_time2 = Hmk$first_treat_time
first_treat$first_treat_time2 = first_treat$first_treat_time
data = merge(data,first_treat[,c('id','first_treat_time2')])
Hmk = merge(Hmk,data[data$time==data$first_treat_time2,c("id","first_treat_time2",'L_gamma')],all.x=T)
Hmk$Q_km = Hmk$k - Hmk$m + 1
Hmk$Q_L = Hmk$L_m
Hmk = Hmk[order(Hmk$id,Hmk$m,Hmk$k),]
data$past_A = unlist(lapply(split(data$A,data$id),function(x)cumsum(cumsum(x!=0))>1))

sample1 = sample(N,ceiling(N/2),replace=F)
sample2 = setdiff(1:N,sample1)
data1 = data[data$id %in% sample1,]
data2 = data[data$id %in% sample2,]
treat_mod1 = ranger(A~time+L+L_mean_m,data=data1[!data1$past_A&data1$time>0,])
treat_mod2 = ranger(A~time+L+L_mean_m,data=data2[!data2$past_A&data2$time>0,])
data1$E_A = predict(treat_mod2,data = data1)$predictions
data2$E_A = predict(treat_mod1,data = data2)$predictions

Hmk1 = Hmk[Hmk$id %in% sample1,]
Hmk2 = Hmk[Hmk$id %in% sample2,]

Hmk1 = merge(Hmk1,data1[,c("id","m",'E_A')],all.x=T)
Hmk2 = merge(Hmk2,data2[,c("id","m",'E_A')],all.x=T)

#Restrict to rows that actually contribute to estimating equations
est_mat1 = Hmk1[Hmk1$k>=Hmk1$m & Hmk1$m<=Hmk1$first_treat_time,]
est_mat2 = Hmk2[Hmk2$k>=Hmk2$m & Hmk2$m<=Hmk2$first_treat_time,]


blip2 = function(m,k,l,a,psi){
  a*sum(c(k-m+1,l)*psi)
}

compute_est_eq_psi_A = function(psi_hat,H_ta=Hmk1){
  gamma_ta = sapply(1:nrow(H_ta),function(i) blip2(H_ta$first_treat_time[i],H_ta$k[i],H_ta$L_gamma[i],H_ta$A[i],psi_hat))
  H_ta$H = ifelse(H_ta$first_treat_time2<=H_ta$k,H_ta$Y - gamma_ta,H_ta$Y)
  H_ta = H_ta %>%
    group_by(id,m) %>%
    dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  H_ta$H_diffs = H_ta$H - H_ta$H_lag
  est_mat_ta = H_ta[H_ta$k>=H_ta$m & H_ta$m<=H_ta$first_treat_time,]
  temp = est_mat_ta$H_diffs*(est_mat_ta[,grep('Q_',names(est_mat_ta))]*(est_mat_ta$A_m-est_mat_ta$E_A))
  apply(temp,2,sum)
}

ss_1_A = nleqslv(x=rep(0,2),fn=compute_est_eq_psi_A,H_ta=Hmk1)
psi_hat_1_A = ss_1_A$x
ss_1_A$termcd

gamma_2 = sapply(1:nrow(Hmk2),function(i) blip2(Hmk2$first_treat_time[i],Hmk2$k[i],Hmk2$L_gamma[i],Hmk2$A[i],psi_hat_1_A))
Hmk2$H = ifelse(Hmk2$first_treat_time<=Hmk2$k,Hmk2$Y - gamma_2,Hmk2$Y)
Hmk2 = Hmk2 %>%
  group_by(id,m) %>%
  dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
Hmk2$H_diffs = Hmk2$H - Hmk2$H_lag
H_mod2 = ranger(H_diffs~m+k+L_m+L_mean_m,data=Hmk2[Hmk2$k>=Hmk2$m & Hmk2$m<=Hmk2$first_treat_time,])

ss_2_A = nleqslv(x=rep(0,2),fn=compute_est_eq_psi_A,H_ta=Hmk2)
psi_hat_2_A = ss_2_A$x
ss_2_A$termcd

gamma_1 = sapply(1:nrow(Hmk1),function(i) blip2(Hmk1$first_treat_time[i],Hmk1$k[i],Hmk1$L_gamma[i],Hmk1$A[i],psi_hat_2_A))
Hmk1$H = ifelse(Hmk1$first_treat_time<=Hmk1$k,Hmk1$Y - gamma_1,Hmk1$Y)
Hmk1 = Hmk1 %>%
  group_by(id,m) %>%
  dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
Hmk1$H_diffs = Hmk1$H - Hmk1$H_lag
H_mod1 = ranger(H_diffs~m+k+L_m+L_mean_m,data=Hmk1[Hmk1$k>=Hmk1$m & Hmk1$m<=Hmk1$first_treat_time,])

compute_est_eq_psi = function(psi_hat,est=1){
  if(est==1){
    H_ta=Hmk1
  }else{
    H_ta=Hmk2
  }
  gamma_ta = sapply(1:nrow(H_ta),function(i) blip2(H_ta$first_treat_time[i],H_ta$k[i],H_ta$L_gamma[i],H_ta$A[i],psi_hat))
  H_ta$H = ifelse(H_ta$first_treat_time2<=H_ta$k,H_ta$Y - gamma_ta,H_ta$Y)
  H_ta = H_ta %>%
    group_by(id,m) %>%
    dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  H_ta$H_diffs = H_ta$H - H_ta$H_lag
  
  # gamma_tr = sapply(1:nrow(H_tr),function(i) blip2(H_tr$first_treat_time[i],H_tr$k[i],H_tr$L_gamma[i],H_tr$A[i],psi_hat))
  # H_tr$H = ifelse(H_tr$first_treat_time<=H_tr$k,H_tr$Y - gamma_tr,H_tr$Y)
  # H_tr = H_tr %>%
  #   group_by(id,m) %>%
  #   dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  # H_tr$H_diffs = H_tr$H - H_tr$H_lag
  # 
  # # H_mod = lm(H_diffs~L_hpi+m+k+m*k+m*L_hpi+k*L_hpi+I(L_hpi^2)+I(m^2)+I(k^2),data=Hmk[Hmk$k>=Hmk$m & Hmk$m>1994 & Hmk$m<=Hmk$first_dereg_year2,c('L_hpi','m','k','H_diffs')])
  # # H_mod = lm(H_diffs~.,data=Hmk[Hmk$k>=Hmk$m & Hmk$m>1994 & Hmk$m<=Hmk$first_dereg_year2,c('L_hpi','m','k','H_diffs')])
  # H_mod = ranger(H_diffs~m+k+L_m+L_mean_m,data=H_tr[H_tr$k>=H_tr$m & H_tr$m<=H_tr$first_treat_time,])
  if(est==1){
    H_mod = H_mod2
  }else{
    H_mod = H_mod1
  }
  est_mat_ta = H_ta[H_ta$k>=H_ta$m & H_ta$m<=H_ta$first_treat_time,]
  est_mat_ta$V = predict(H_mod,data=est_mat_ta)$predictions
  temp = (est_mat_ta$H_diffs-est_mat_ta$V)*(est_mat_ta[,grep('Q_',names(est_mat_ta))]*(est_mat_ta$A_m-est_mat_ta$E_A))
  apply(temp,2,sum)
}

compute_est_eq_psi_2 = function(psi_hat,H_ta=Hmk2,H_tr=Hmk1){
  gamma_ta = sapply(1:nrow(H_ta),function(i) blip2(H_ta$first_treat_time[i],H_ta$k[i],H_ta$L_gamma[i],H_ta$A[i],psi_hat))
  H_ta$H = ifelse(H_ta$first_treat_time2<=H_ta$k,H_ta$Y - gamma_ta,H_ta$Y)
  H_ta = H_ta %>%
    group_by(id,m) %>%
    dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  H_ta$H_diffs = H_ta$H - H_ta$H_lag
  
  gamma_tr = sapply(1:nrow(H_tr),function(i) blip2(H_tr$first_treat_time[i],H_tr$k[i],H_tr$L_gamma[i],H_tr$A[i],psi_hat))
  H_tr$H = ifelse(H_tr$first_treat_time<=H_tr$k,H_tr$Y - gamma_tr,H_tr$Y)
  H_tr = H_tr %>%
    group_by(id,m) %>%
    dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  H_tr$H_diffs = H_tr$H - H_tr$H_lag
  
  # H_mod = lm(H_diffs~L_hpi+m+k+m*k+m*L_hpi+k*L_hpi+I(L_hpi^2)+I(m^2)+I(k^2),data=Hmk[Hmk$k>=Hmk$m & Hmk$m>1994 & Hmk$m<=Hmk$first_dereg_year2,c('L_hpi','m','k','H_diffs')])
  # H_mod = lm(H_diffs~.,data=Hmk[Hmk$k>=Hmk$m & Hmk$m>1994 & Hmk$m<=Hmk$first_dereg_year2,c('L_hpi','m','k','H_diffs')])
  H_mod = ranger(H_diffs~m+k+L_m+L_mean_m,data=H_tr[H_tr$k>=H_tr$m & H_tr$m<=H_tr$first_treat_time,])
  est_mat_ta = H_ta[H_ta$k>=H_ta$m & H_ta$m<=H_ta$first_treat_time,]
  est_mat_ta$V = predict(H_mod,data=est_mat_ta)$predictions
  temp = (est_mat_ta$H_diffs-est_mat_ta$V)*(est_mat_ta[,grep('Q_',names(est_mat_ta))]*(est_mat_ta$A_m-est_mat_ta$E_A))
  apply(temp,2,sum)
}

ss_1 = nleqslv(x=rep(0,2),fn=compute_est_eq_psi,est=1)
psi_hat_1 = ss_1$x
ss_1$termcd
ss_1$fvec


ss_2 = nleqslv(x=rep(0,2),fn=compute_est_eq_psi,est=2)
psi_hat_2 = ss_2$x
ss_2$termcd
ss_2$fvec

psi_hat_cf = (psi_hat_1+psi_hat_2)/2



# #generate estimates of derived quantities using blip function parameter
# cd_estimands = rep(NA,8)
# start_years = c(1995:1998,2000,2001)
# start_psi_inds = c(1,12,22,31,39,45)
# for(i in 0:7){
#   gammas = as.vector(as.matrix(est_mat[,R_inds])%*%psi_hat_dr[1:50])
#   gammas = gammas[est_mat$m==est_mat$first_treat_time & (est_mat$k-est_mat$m)==i]
#   cd_estimands[i+1] = mean(gammas)
# }
# 
# est_mat$gammas=ifelse(est_mat$k<est_mat$first_treat_time,0,as.matrix(est_mat[,R_inds])%*%psi_hat_dr[1:50])
# est_mat$H = est_mat$k_hpi-est_mat$gammas
# no_treat_hat = rep(NA,11)
# years=1995:2005
# for(i in 1:11){
#   no_treat_hat[i] = mean(est_mat$H[est_mat$m==1995 & est_mat$k==years[i]])
# }
# observed_outcome_trajectory = aggregate(data$log_hpi[data$year>=1995],by=list(data$k[data$year>=1995]),mean)[,2]
# observed_outcome_trajectory-no_treat_hat
# 
# B_list = vector(mode='list',length=length(S_diff_list))
# for(i in 1:length(S_diff_list)){
#   B_list[[i]] = get_B(S_diff_list[[i]],Q_list[[i]],est_mat$A_any_m[i]-est_mat$p_any[i])
# }
# B = Reduce('+',B_list)
# psi_hat = as.vector(A%*%solve(B))
# 
# ################
# #BOOTSTRAP
# ################
# 
# nboot = 200
# psi_hat_dr_boot = vector(mode='list',length=nboot)
# no_treat_hat_boot = vector(mode='list',length=nboot)
# observed_outcome_trajectory_boot = vector(mode='list',length=nboot)
# cd_estimands_boot = vector(mode='list',length=nboot)
# counties = unique(Hmk$id)
# data_inds = lapply(1:length(counties),function(i)((i-1)*12+1):(i*12))
# b=1
# count=1
# while(count <= nboot){
#   set.seed(b)
#   samp_counties = sample(1:length(counties),length(counties),replace=T)
#   mod_inds = unlist(data_inds[samp_counties])
#   data_boot = data[mod_inds,]
#   data_boot$id = rep(1:length(counties),each=12)
#   if(sum(data_boot$first_dereg_year3==1995)>0){
#     Hmk_boot = data.frame(cbind(id=rep(1:length(counties),each=77),m=rep(rep(rep(1995:2005,c(12:2)),length(counties))),
#                                       k=rep(c(1994:2005,1995:2005,1996:2005,1997:2005,1998:2005,1999:2005,2000:2005,2001:2005,2002:2005,2003:2005,2004:2005),length(id_keepers))))
#     first_dereg_year_boot = aggregate(data_boot$inter_bra,by=list(data_boot$id),function(x) min(which(x>0)))
#     names(first_dereg_year_boot) = c('id','first_dereg_year')
#     
#     Hmk_boot = merge(Hmk_boot,first_dereg_year_boot)
#     Hmk_boot$first_treat_time = Hmk_boot$first_dereg_year + 1994 - 1
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","first_treat_time","inter_bra")],all.x = T)
#     
#     Hmk_boot$inter_bra[is.na(Hmk_boot$inter_bra)]=0
#     Hmk_boot = Hmk_boot[order(Hmk_boot$id,Hmk_boot$m,Hmk_boot$k),]
#     
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","lag_hpi","first_treat_time")],all.x = T)
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","lag_amt","first_treat_time")],all.x = T)
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","k_amt","k")],all.x = T)
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","k_hpi","k")],all.x = T)
#     
#     # A_mod_boot = glm(inter_bra~year+I(year^2)+I(year^3)+lag_hpi+I(lag_hpi^2),data=data_boot[!data_boot$past_A&data_boot$year>1994,],family="poisson")
#     A_mod_any_boot = glm(any_A~factor(year)*lag_amt+I(lag_amt^2),data=data_boot[!data_boot$past_A&data_boot$year>1994,],family="binomial")
#     data_boot$p_any = NA
#     data_boot$p_any[data_boot$year>1994] = predict(A_mod_any_boot,newdata = data_boot[data_boot$year>1994,],type='response')
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","m",'p_any')],all.x=T)
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","m",'A_m')],all.x=T)
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","m",'L_amt','L_hpi')],all.x=T)
#     
#     Hmk_boot = merge(Hmk_boot,data_boot[data_boot$year==data_boot$first_dereg_year3,c("id","year",'L_gamma_amt','L_gamma_hpi')],all.x=T)
#     Hmk_boot = merge(Hmk_boot,data_boot[,c("id","k",'lag_k_amt','lag_k_hpi')],all.x=T)
#     
#     Hmk_boot = Hmk_boot[order(Hmk_boot$id,Hmk_boot$m,Hmk_boot$k),]
#     Hmk_boot$A_any_m = as.numeric(Hmk_boot$A_m>0)
#     Hmk_boot$inter_bra_any = as.numeric(Hmk_boot$inter_bra>0)
#     
#     est_mat_boot = Hmk_boot[Hmk_boot$k>=Hmk_boot$m & Hmk_boot$m%in%c(1995:1998,2000,2001) & Hmk_boot$m<=Hmk_boot$first_treat_time,]
#     for(i in c(1995:1998,2000,2001)){
#       for(j in i:2005){
#         est_mat_boot[,paste0('Q_',i,'_',j)] = as.numeric(est_mat_boot$m==i & est_mat_boot$k==j)
#       }
#     }
#     est_mat_boot$Q_amt = est_mat_boot$L_amt
#     
#     for(i in c(1995:1998,2000,2001)){
#       for(j in i:2005){
#         est_mat_boot[,paste0('R_',i,'_',j)] = as.numeric(est_mat_boot$first_treat_time==i & est_mat_boot$k==j)
#       }
#     }
#     est_mat_boot$R_amt = est_mat_boot$L_gamma_amt
#     
#     for(i in c(1995:1998,2000,2001)){
#       for(j in i:2005){
#         est_mat_boot[,paste0('D_',i,'_',j,'_amt')] = as.numeric(est_mat_boot$m==i & est_mat_boot$k==j)*est_mat_boot$L_amt
#       }
#     }
#     
#     A_dr_boot = colSums((est_mat_boot$k_hpi-est_mat_boot$lag_k_hpi)*cbind((est_mat_boot$A_any_m-est_mat_boot$p_any)*est_mat_boot[,grep('Q_',names(est_mat_boot))],est_mat_boot[,grep('D_',names(est_mat_boot))]))
#     
#     R_inds=grep('R_',names(est_mat_boot))
#     D_inds=grep('D_',names(est_mat_boot))
#     Q_inds=grep('Q_',names(est_mat_boot))
#     p = length(Q_inds)
#     get_S = function(ind){
#       if(est_mat_boot$first_treat_time[ind]>est_mat_boot$k[ind]){
#         return(rep(0,p))
#       }else{
#         return(as.vector(est_mat_boot[ind,R_inds]))
#       }
#     }
#     
#     get_S_lag = function(ind){
#       if(est_mat_boot$first_treat_time[ind]>(est_mat_boot$k[ind]-1)){
#         return(rep(0,p))
#       }else{
#         return(as.vector(est_mat_boot[ind-1,R_inds]))
#       }
#     }
#     
#     get_D = function(ind){
#       as.vector(est_mat_boot[ind,D_inds])
#     }
#     
#     get_Q = function(ind){
#       as.vector(est_mat_boot[ind,Q_inds])
#     }
#     
#     S_diff_list_boot = lapply(1:nrow(est_mat_boot),function(i) get_S(i)-get_S_lag(i))
#     D_list_boot = lapply(1:nrow(est_mat_boot),get_D)
#     Q_list_boot = lapply(1:nrow(est_mat_boot),get_Q)
#     
#     B_list_dr_boot = vector(mode='list',length=length(S_diff_list_boot))
#     for(i in 1:length(S_diff_list_boot)){
#       B_list_dr_boot[[i]] = get_B_dr(S_diff_list_boot[[i]],D_list_boot[[i]],Q_list_boot[[i]],est_mat_boot$A_any_m[i]-est_mat_boot$p_any[i])
#     }
#     B_dr_boot = Reduce('+',B_list_dr_boot)
#     psi_hat_dr_b = as.vector(A_dr_boot%*%solve(B_dr_boot))
#     psi_hat_dr_boot[[count]] = psi_hat_dr_b
#     
#     cd_estimands_b = rep(NA,8)
#     start_years = c(1995:1998,2000,2001)
#     for(i in 0:7){
#       gammas = as.vector(as.matrix(est_mat_boot[,R_inds])%*%psi_hat_dr_b[1:50])
#       gammas = gammas[est_mat_boot$m==est_mat_boot$first_treat_time & (est_mat_boot$k-est_mat_boot$m)==i]
#       cd_estimands_b[i+1] = mean(gammas)
#     }
#     cd_estimands_boot[[count]] = cd_estimands_b
#     
#     est_mat_boot$gammas=ifelse(est_mat_boot$k<est_mat_boot$first_treat_time,0,as.matrix(est_mat_boot[,R_inds])%*%psi_hat_dr_b[1:50])
#     est_mat_boot$H = est_mat_boot$k_hpi-est_mat_boot$gammas
#     no_treat_hat_b = rep(NA,11)
#     years=1995:2005
#     for(i in 1:11){
#       no_treat_hat_b[i] = mean(est_mat_boot$H[est_mat_boot$m==1995 & est_mat_boot$k==years[i]])
#     }
#     no_treat_hat_boot[[count]] = no_treat_hat_b
#     observed_outcome_trajectory_boot[[count]] = aggregate(data_boot$log_hpi[data_boot$year>=1995],by=list(data_boot$k[data$year>=1995]),mean)[,2]
#     count=count+1
#     b=b+1
#   }else{
#     b=b+1
#   }
# }
# 
# save(psi_hat_dr_boot,no_treat_hat_boot,cd_estimands_boot,observed_outcome_trajectory_boot,
#      psi_hat_dr,no_treat_hat,observed_outcome_trajectory,cd_estimands,
#      file="did_snmm_coarse_hpi_Lamt_results.RData")
# 
# load("did_snmm_coarse_hpi_Lamt_results.RData")
# no_treat_hat_boot_mat = data.frame(no_treat_hat_boot)
# observed_outcome_trajectory_boot_mat = data.frame(observed_outcome_trajectory_boot)
# no_treat_effect_boot_mat = no_treat_hat_boot_mat - observed_outcome_trajectory_boot_mat
# no_treat_effect_sds = apply(no_treat_effect_boot_mat,1,sd)
# uppers = no_treat_hat - observed_outcome_trajectory + 1.96*no_treat_effect_sds
# lowers = no_treat_hat - observed_outcome_trajectory - 1.96*no_treat_effect_sds
# years=1995:2005
# pdf("did_snmm_coarse_hpi_Lamt_plot.pdf")
# plot(years,no_treat_hat - observed_outcome_trajectory,type='l',ylim=c(min(lowers),max(uppers)),main="Effect on log(HPI) of Removing All Deregulation vs Usual Care",ylab="Effect on log(HPI)",xlab="Year")
# points(years,lowers,type='l',lty=2)
# points(years,uppers,type='l',lty=2)
# curve(0*x,from=1990,to=3000,col=2,add=T)
# dev.off()
# 
# psi_hat_dr_mat = data.frame(psi_hat_dr_boot)
# pdf('did_snmm_coarse_psihat_amt.pdf')
# hist(unlist(psi_hat_dr_mat[50,]),main="Bootstrap Distribution of Mortgage Volume Coefficient",xlab=expression(beta))
# points(psi_hat_dr[50],0,col=2)
# dev.off()
# 
# psi_hat_dr[50]
# psi_hat_dr[50] + 1.96*sd(psi_hat_dr_mat[50,])
# psi_hat_dr[50] - 1.96*sd(psi_hat_dr_mat[50,])
# 
# 
# cd_estimands_boot_mat = data.frame(cd_estimands_boot)
# cd_estimands_sds = apply(cd_estimands_boot_mat,1,sd)
# uppers = cd_estimands + 1.96*cd_estimands_sds
# lowers = cd_estimands - 1.96*cd_estimands_sds
# pdf("did_snmm_coarse_cd_estimands_Lamt_plot.pdf")
# plot(0:7,cd_estimands,type='l',ylim=c(min(lowers),max(uppers)),main="Weighted Average ITT-ETTs on log(HPI)",ylab="Effect on log(HPI)",xlab="Year After Treatment")
# points(0:7,lowers,type='l',lty=2)
# points(0:7,uppers,type='l',lty=2)
# curve(0*x,from=0,to=3000,col=2,add=T)
# dev.off()