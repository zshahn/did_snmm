library(haven)
library(nleqslv)
# setwd('C:/Users/zshahn/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
# setwd('~/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
setwd('~/snmm/flood')
data = read_dta('panel_1980_2007.dta')
data$total_floods_pre_1980 = apply(data[,329:350],1,sum)
data$year = as.numeric(data$year)
data$running_flood_rate = (data$total_floods_pre_1980 + unlist(lapply(split(data$disaster_year,data$id),cumsum)))/(data$year-1957)
data$last_flood_year = NA
data$last_flood_year[data$disaster_year==1]=data$year[data$disaster_year==1]
library(tidyverse)
data = data %>% 
  group_by(id) %>%
  fill(last_flood_year,.direction = "down")
data$prev_flood_year = data$last_flood_year
post1980_flood_inds = which(data$disaster_year==1 & data$year>1980)
data$prev_flood_year[post1980_flood_inds] = data$last_flood_year[post1980_flood_inds-1]
data$prev_flood_year[data$year==1980] = NA
flood_hist = data[,329:350]
most_recent_pre1980_floods = apply(flood_hist,1,function(x)1957+max(which(x==1)))
data$prev_flood_year = ifelse(is.na(data$prev_flood_year),most_recent_pre1980_floods,data$prev_flood_year)
data$ts_prev_flood = data$year - data$prev_flood_year
data$ts_prev_flood_short = data$ts_prev_flood<=5
data$ts_prev_flood_medium = data$ts_prev_flood>5 & data$ts_prev_flood<=15
data$ts_prev_flood_long = data$ts_prev_flood>15 

get_last = function(x){
  if(length(x)>1){
    x=c(NA,x[1:(length(x)-1)])
  }else{
    x=NA
  }
  x
}
pre1980_flood_rates = apply(flood_hist,1,sum)/ncol(flood_hist)
data$last_running_flood_rate = unlist(lapply(split(data$running_flood_rate,data$id),get_last))
data$last_running_flood_rate[data$year==1980] = pre1980_flood_rates[data$year==1980]
data$running_flood_rate_5yr = NA
for(yr in 1980:2007){
  yrs = (yr-5):(yr-1)
  data$running_flood_rate_5yr[data$year==yr] = data[data$year==yr,paste0('disyr_cnt',yrs[1])]
}
data$Y = data$ln_takeup_proxy1
# data$Y = unlist(lapply(split(data$ln_takeup_proxy1,data$id),function(x)c(x[2:length(x)],NA)))
data$id = as.numeric(data$id)
ids = unique(data$id)
Hmk_mat_c = data.frame(cbind(id=rep(ids,each=405),m=rep(rep(rep(1981:2007,c(28:2)),length(ids))),
                             k=rep(c(1980:2007,1981:2007,1982:2007,1983:2007,1984:2007,1985:2007,1986:2007,1987:2007,1988:2007,1989:2007,
                                     1990:2007,1991:2007,1992:2007,1993:2007,1994:2007,1995:2007,1996:2007,1997:2007,1998:2007,1999:2007,
                                     2000:2007,2001:2007,2002:2007,2003:2007,2004:2007,2005:2007,2006:2007),length(ids))))

Hmk_mat_c = merge(Hmk_mat_c,data[data$year==1980,c('id',grep('disyr_cnt',names(data),value=T))])
for(yr in 1981:2007){
  last_flood_rate_yr = data[data$year==yr,c('id',"last_running_flood_rate")]
  names(last_flood_rate_yr)[2] = paste0('L_flood_rate_',yr)
  # ts_prev_flood_short_yr = data[data$year==yr,c('id',"ts_prev_flood_short")]
  # names(ts_prev_flood_short_yr)[2] = paste0('L_ts_prev_flood_short_',yr)
  Hmk_mat_c = merge(Hmk_mat_c,last_flood_rate_yr)
  # Hmk_mat_c = merge(Hmk_mat_c,ts_prev_flood_short_yr)
}

data$k = data$year
Hmk_mat_c = merge(Hmk_mat_c,data[,c('id','k','Y')])
Hmk_mat_c = Hmk_mat_c[order(Hmk_mat_c$id,Hmk_mat_c$m,Hmk_mat_c$k),]


# A_mod = glm(disaster_year~factor(year)*last_running_flood_rate*(ts_prev_flood_short+ts_prev_flood_long),data=data[data$year>1980,],family='binomial')
A_mod = glm(disaster_year~factor(year)*last_running_flood_rate,data=data,family='binomial')
data$A_hat = predict(A_mod,newdata=data,type='response')
data$m = data$year
Hmk_mat_c = merge(Hmk_mat_c,data[,c("id","m",'A_hat')],all.x=T)
Hmk_mat_c = merge(Hmk_mat_c,data[,c("id","m",'disaster_year')],all.x=T)

Hmk_mat_c = merge(Hmk_mat_c,data[,c("id","m","last_running_flood_rate")],all.x=T)
Hmk_mat_c = Hmk_mat_c[order(Hmk_mat_c$id,Hmk_mat_c$m,Hmk_mat_c$k),]

Hmk_mat_c$S1 = 1
Hmk_mat_c$S2 = Hmk_mat_c$m-1980
Hmk_mat_c$S3 = Hmk_mat_c$k-Hmk_mat_c$m
Hmk_mat_c$S4 = (Hmk_mat_c$k-Hmk_mat_c$m)^2
Hmk_mat_c$S5 = Hmk_mat_c$last_running_flood_rate

#A*(1,m-1980,k-m,flood rate,short ts prev flood,long ts prev flood)*psi
years = 1981:2007
compute_est_eq_psi = function(psi,Hmk=Hmk_mat_c){
  psi = matrix(rep(psi,each=nrow(Hmk)),nrow=nrow(Hmk))
  gammas = matrix(NA,ncol=length(years),nrow=nrow(Hmk))
  colnames(gammas) = paste0('gamma_',years,'_k')
  for(j in 1:length(years)){
    gammas[,paste0('gamma_',years[j],'_k')] = ifelse((Hmk$k<years[j])|(Hmk$m>years[j]),0,rowSums(Hmk[,paste0("disyr_cnt",years[j])]*cbind(1,years[j]-1980,Hmk$k-years[j],(Hmk$k-years[j])^2,Hmk[,paste0("L_flood_rate_",years[j])])*psi))
  }
  gamma = rowSums(gammas)
  Hmk$H = Hmk$Y - gamma
  Hmk = Hmk %>%
    group_by(id,m) %>%
    dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  Hmk$H_diffs = Hmk$H - Hmk$H_lag
  H_mod = lm(H_diffs~factor(m)*last_running_flood_rate+I(last_running_flood_rate^2)+k+I(k^2)+I(k^3),data=Hmk[Hmk$k>=Hmk$m & Hmk$m>1980,c('last_running_flood_rate','m','k','H_diffs')])
  Hmk$V = predict(H_mod,newdata=Hmk)
  est_mat = Hmk[Hmk$k>=Hmk$m & Hmk$m>1980,]
  temp = (est_mat$H_diffs-est_mat$V)*(est_mat[,c('S1','S2','S3','S4','S5')]*(est_mat$disaster_year-est_mat$A_hat))
  colSums(temp)
}

Sys.time()
ss = nleqslv(x=rep(0,5),fn=compute_est_eq_psi)
Sys.time()
psi_hat = ss$x
ss$termcd
ss$fvec


no_treat_hat = rep(NA,27)
for(i in 1:27){
  no_treat_hat[i] = mean(Hmk_mat_c$H[Hmk_mat_c$m==1981 & Hmk_mat_c$k==years[i]])
}

observed_outcome_trajectory = aggregate(Hmk_mat_c$Y[Hmk_mat_c$k>=1981],by=list(Hmk_mat_c$k[Hmk_mat_c$k>=1981]),mean)[,2]


compute_est_eq_psi_boot = function(psi,Hmk=Hmk_mat_c_boot){
  psi = matrix(rep(psi,each=nrow(Hmk)),nrow=nrow(Hmk))
  gammas = matrix(NA,ncol=length(years),nrow=nrow(Hmk))
  colnames(gammas) = paste0('gamma_',years,'_k')
  for(j in 1:length(years)){
    gammas[,paste0('gamma_',years[j],'_k')] = ifelse((Hmk$k<years[j])|(Hmk$m>years[j]),0,rowSums(Hmk[,paste0("disyr_cnt",years[j])]*cbind(1,years[j]-1980,Hmk$k-years[j],(Hmk$k-years[j])^2,Hmk[,paste0("L_flood_rate_",years[j])])*psi))
  }
  gamma = rowSums(gammas)
  Hmk$H = Hmk$Y - gamma
  Hmk = Hmk %>%
    group_by(id,m) %>%
    dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  Hmk$H_diffs = Hmk$H - Hmk$H_lag
  H_mod = lm(H_diffs~factor(m)*last_running_flood_rate+I(last_running_flood_rate^2)+k+I(k^2)+I(k^3),data=Hmk[Hmk$k>=Hmk$m & Hmk$m>1980,c('last_running_flood_rate','m','k','H_diffs')])
  Hmk$V = predict(H_mod,newdata=Hmk)
  est_mat = Hmk[Hmk$k>=Hmk$m & Hmk$m>1980,]
  temp = (est_mat$H_diffs-est_mat$V)*(est_mat[,c('S1','S2','S3','S4','S5')]*(est_mat$disaster_year-est_mat$A_hat))
  colSums(temp)
}

nboot = 50
psi_hat_boot = vector(mode='list',length=nboot)
ids = unique(Hmk_mat_c$id)
data_inds = lapply(1:length(ids),function(i)((i-1)*28+1):(i*28))
for(b in 1:nboot){
  set.seed(b)
  samp_ids = sample(1:length(ids),length(ids),replace=T)
  mod_inds = unlist(data_inds[samp_ids])
  data_boot = data[mod_inds,]
  data_boot$id = rep(1:length(ids),each=28)
  
  Hmk_mat_c_boot = data.frame(cbind(id=rep(1:length(ids),each=405),m=rep(rep(rep(1981:2007,c(28:2)),length(ids))),
                                    k=rep(c(1980:2007,1981:2007,1982:2007,1983:2007,1984:2007,1985:2007,1986:2007,1987:2007,1988:2007,1989:2007,
                                            1990:2007,1991:2007,1992:2007,1993:2007,1994:2007,1995:2007,1996:2007,1997:2007,1998:2007,1999:2007,
                                            2000:2007,2001:2007,2002:2007,2003:2007,2004:2007,2005:2007,2006:2007),length(ids))))
  
  Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[data_boot$year==1980,c('id',grep('disyr_cnt',names(data_boot),value=T))])
  for(yr in 1981:2007){
    last_flood_rate_yr = data_boot[data_boot$year==yr,c('id',"last_running_flood_rate")]
    names(last_flood_rate_yr)[2] = paste0('L_flood_rate_',yr)
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,last_flood_rate_yr)
  }
  
  Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c('id','k','Y')])
  Hmk_mat_c_boot = Hmk_mat_c_boot[order(Hmk_mat_c_boot$id,Hmk_mat_c_boot$m,Hmk_mat_c_boot$k),]
  
  # A_mod = glm(disaster_year~factor(year)*last_running_flood_rate*(ts_prev_flood_short+ts_prev_flood_long),data=data[data$year>1980,],family='binomial')
  A_mod_boot = glm(disaster_year~factor(year)*last_running_flood_rate,data=data_boot,family='binomial')
  data_boot$A_hat = predict(A_mod_boot,newdata=data_boot,type='response')
  Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("id","m",'A_hat')],all.x=T)
  Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("id","m",'disaster_year')],all.x=T)
  
  Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("id","m","last_running_flood_rate")],all.x=T)
  Hmk_mat_c_boot = Hmk_mat_c_boot[order(Hmk_mat_c_boot$id,Hmk_mat_c_boot$m,Hmk_mat_c_boot$k),]
  
  Hmk_mat_c_boot$S1 = 1
  Hmk_mat_c_boot$S2 = Hmk_mat_c_boot$m-1980
  Hmk_mat_c_boot$S3 = Hmk_mat_c_boot$k-Hmk_mat_c_boot$m
  Hmk_mat_c_boot$S4 = (Hmk_mat_c_boot$k-Hmk_mat_c_boot$m)^2
  Hmk_mat_c_boot$S5 = Hmk_mat_c_boot$last_running_flood_rate
  ss_boot = nleqslv(x=rep(0,5),fn=compute_est_eq_psi_boot)
  psi_hat_boot[[b]] = ss_boot$x
}
save(psi_hat,psi_hat_boot,file='flood_results_boot.RData')

load('flood_results_boot.RData')
psi_hat
psi_hat_boot = data.frame(psi_hat_boot)
point_ests = vector(mode='list',length=27)
sds = vector(mode='list',length=27)
uppers = vector(mode='list',length=27)
lowers = vector(mode='list',length=27)

point_ests_high = vector(mode='list',length=27)
sds_high = vector(mode='list',length=27)
uppers_high = vector(mode='list',length=27)
lowers_high = vector(mode='list',length=27)

est_traj = function(m,b1,b2,b3,b4,b5,L=median(data$last_running_flood_rate)){
  b1+(m-1980)*b2+b3*(m:min(2007,m+15)-m)+b4*(m:min(2007,m+15)-m)^2+b5*L
}

for(i in 1:27){
  m=1980+i
  point_ests[[i]]=psi_hat[1]+(m-1980)*psi_hat[2]+psi_hat[3]*(m:min(2007,m+15)-m)+psi_hat[4]*(m:min(2007,m+15)-m)^2 +psi_hat[5]*median(data$last_running_flood_rate)
  boot_ests = vector(mode='list',length=50)
  for(j in 1:50){
    boot_ests[[j]] = est_traj(m,psi_hat_boot[1,j],psi_hat_boot[2,j],psi_hat_boot[3,j],psi_hat_boot[4,j],psi_hat_boot[5,j])
  }
  boot_ests = data.frame(boot_ests)
  sds[[i]] = apply(boot_ests,1,sd)
  uppers[[i]] = point_ests[[i]] + 1.96*sds[[i]]
  lowers[[i]] = point_ests[[i]] - 1.96*sds[[i]]
  
  # point_ests_high[[i]]=psi_hat[1]+(m-1980)*psi_hat[2]+psi_hat[3]*(m:min(2007,m+15)-m) +psi_hat[4]*.2
  # boot_ests_high = vector(mode='list',length=34)
  # for(j in 1:34){
  #   boot_ests_high[[j]] = est_traj(m,psi_hat_boot[1,j],psi_hat_boot[2,j],psi_hat_boot[3,j],psi_hat_boot[4,j],L=.2)
  # }
  # boot_ests_high = data.frame(boot_ests_high)
  # sds_high[[i]] = apply(boot_ests_high,1,sd)
  # uppers_high[[i]] = point_ests_high[[i]] + 1.96*sds_high[[i]]
  # lowers_high[[i]] = point_ests_high[[i]] - 1.96*sds_high[[i]]
}
pdf('temp.pdf')
plot(1981:1996,point_ests[[1]],type='l',xlim=c(1980,2007),ylim=c(-.05,.1),main='Effects of Floods on Flood Insurance Takeup (Flood Rate=0.1)',xlab='Year',ylab='Effect (log(take-up proxy))')
points(1981:1996,uppers[[1]],type='l',lty='dashed')       
points(1981:1996,lowers[[1]],type='l',lty='dashed')

# points(1981:1991,point_ests_high[[1]],type='l',col=2)
# points(1981:1991,uppers_high[[1]],type='l',lty='dashed',col=2)       
# points(1981:1991,lowers_high[[1]],type='l',lty='dashed',col=2)


for(i in c(5,10,15,20,25)){
  m=1980+i
  points(m:min(2007,m+15),point_ests[[i]],type='l')
  points(m:min(2007,m+15),uppers[[i]],type='l',lty='dashed')       
  points(m:min(2007,m+15),lowers[[i]],type='l',lty='dashed')
}
dev.off()




