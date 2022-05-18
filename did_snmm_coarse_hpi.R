#effect of bank deregulation on home price index
#binary treatment, coarse blip function
#blip function is linear with predictors: indicators for treatment/outcome years, lagged mortgage volume
#doubly robust estimator

#read and clean data
library(haven)
library(nleqslv)
# setwd('C:/Users/zshahn/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
# setwd('~/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
setwd('~/snmm/bank')
tab1 = read_dta("data_tab1.dta")
aggregation = read_dta("data_aggregation.dta")
call = read_dta("call.dta")
distance = read_dta("data_distance.dta")
hmda = read_dta("hmda.dta")
hp_dereg = read_dta("hp_dereg_controls.dta")
data = merge(hp_dereg,hmda)

data = merge(data,aggregation[,c("county","year","hpi","amtoriginated_b")],all.x=T)
data = data[order(data$county,data$year),]
year1995 = which(data$year==1995)
for(i in 1:length(year1995)){
  if(data$year[year1995[i]+1]==1996){
    data$amtoriginated_b[year1995[i]] = data$amtoriginated_b[year1995[i]+1]/exp(data$Dl_vloans_b[year1995[i]+1])
    data$hpi[year1995[i]] = data$hpi[year1995[i]+1]/exp(data$Dl_hpi[year1995[i]+1])
  }
}

year1994 = which(data$year==1994)
for(i in 1:length(year1994)){
  if(data$year[year1994[i]+1]==1995){
    data$amtoriginated_b[year1994[i]] = data$amtoriginated_b[year1994[i]+1]/exp(data$Dl_vloans_b[year1994[i]+1])
    data$hpi[year1994[i]] = data$hpi[year1994[i]+1]/exp(data$Dl_hpi[year1994[i]+1])
  }
}

data$log_hpi = log(data$hpi)
data$log_amt = log(data$amtoriginated_b)

data$log_amt[is.na(data$log_amt)|data$log_amt==-Inf] = median(data$log_amt,na.rm=T)
data$log_hpi[is.na(data$log_hpi)|data$log_hpi==-Inf] = median(data$log_hpi,na.rm=T)


county_years = aggregate(data$year,by=list(data$county),length)
county_keepers = unique(county_years[county_years[,2]==12,1])
data = data[data$county %in% county_keepers,]

#get first year that each county deregulated

first_dereg_year = aggregate(data$inter_bra,by=list(data$county),function(x) min(which(x>0)))
names(first_dereg_year) = c('county','first_dereg_year')
first_dereg_year$first_dereg_year3 = first_dereg_year$first_dereg_year+1994-1
data = merge(data,first_dereg_year[,c("county","first_dereg_year3")])

#Make reformatted data frame Hmk_mat_c
#one row for each county/m/k triplet where m and k are years with k>=m-1 
#ultimately the important columns will be: 
#year of first deregulation
#covariates for blip model at year of first deregulation for county if after m (R_*)
#covariates at year m for blipped down outcome model for year k (D_*)
#estimating equation filler (Q_*)
#estimated treatment probability at time m
#treatment value at time m
#outcome values at years k and k-1
#The below code creates and populates all the above columns, plus some possibly unnecessary ones

Hmk_mat_c = data.frame(cbind(county=rep(county_keepers,each=77),m=rep(rep(rep(1995:2005,c(12:2)),length(county_keepers))),
                             k=rep(c(1994:2005,1995:2005,1996:2005,1997:2005,1998:2005,1999:2005,2000:2005,2001:2005,2002:2005,2003:2005,2004:2005),length(county_keepers))))
Hmk_mat_c = merge(Hmk_mat_c,first_dereg_year)
Hmk_mat_c$first_dereg_year2 = Hmk_mat_c$first_dereg_year + 1994 - 1

data$first_dereg_year2 = data$year
Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","first_dereg_year2","inter_bra")],all.x = T)

data$k = data$year
data$m = data$year

library(tidyverse)
data = data %>%
  group_by(county) %>%
  dplyr::mutate(lag_hpi = lag(log_hpi, n = 1, default = NA))

data = data %>%
  group_by(county) %>%
  dplyr::mutate(lag_amt = lag(log_amt, n = 1, default = NA))

#do some imputation of missing amt and hpi values
data$lag_amt[is.na(data$lag_amt)|data$lag_amt==-Inf] = median(data$lag_amt,na.rm=T)
data$lag_hpi[is.na(data$lag_hpi)|data$lag_hpi==-Inf] = median(data$lag_hpi,na.rm=T)

data$k_hpi = data$log_hpi
data$k_amt = data$log_amt
data$m_hpi = data$log_hpi
data$m_amt = data$log_amt
data = data.frame(data)

Hmk_mat_c$inter_bra[is.na(Hmk_mat_c$inter_bra)]=0
Hmk_mat_c = Hmk_mat_c[order(Hmk_mat_c$county,Hmk_mat_c$m,Hmk_mat_c$k),]

data$past_A = unlist(lapply(split(data$inter_bra,data$county),function(x)cumsum(cumsum(x>0))>1))
data$any_A = data$inter_bra>0

Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","lag_hpi","first_dereg_year2")],all.x = T)
Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","lag_amt","first_dereg_year2")],all.x = T)
Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","k_amt","k")],all.x = T)
Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","k_hpi","k")],all.x = T)

#Fit treatment model
# A_mod = lm(inter_bra~factor(year)*lag_amt+I(lag_amt^2),data=data[!data$past_A&data$year>1994,])
#When taking treatment to be binary, we use the A_mod_any logistic regression with year and the covariate and their interaction as predictors
A_mod_any = glm(any_A~factor(year)*lag_amt+I(lag_amt^2),data=data[!data$past_A&data$year>1994,],family="binomial")

#add predicted treatment probabilities to Hmk_mat_c
# data$A_hat = 0
# data$A_hat[data$year>1994] = predict(A_mod,newdata=data[data$year>1994,])
data$p_any = NA
data$p_any[data$year>1994] = predict(A_mod_any,newdata = data[data$year>1994,],type='response')

Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","m",'p_any')],all.x=T)
# Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","m",'A_hat')],all.x=T)

#add treatment and covariate variables to Hmk_mat_c
data$A_m = data$inter_bra
Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","m",'A_m')],all.x=T)
data$L_amt = data$lag_amt
data$L_hpi = data$lag_hpi
data$lag_k_amt = data$lag_amt
data$lag_k_hpi = data$lag_hpi
data$L_gamma_amt = data$lag_amt
data$L_gamma_hpi = data$lag_hpi

Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","m",'L_amt','L_hpi')],all.x=T)
Hmk_mat_c = merge(Hmk_mat_c,data[,c("county","k",'lag_k_amt','lag_k_hpi')],all.x=T)
Hmk_mat_c = merge(Hmk_mat_c,data[data$year==data$first_dereg_year3,c("county","year",'L_gamma_amt','L_gamma_hpi')],all.x=T)

Hmk_mat_c = Hmk_mat_c[order(Hmk_mat_c$county,Hmk_mat_c$m,Hmk_mat_c$k),]
Hmk_mat_c$A_any_m = as.numeric(Hmk_mat_c$A_m>0)
Hmk_mat_c$inter_bra_any = as.numeric(Hmk_mat_c$inter_bra>0)

#Restrict to rows that actually contribute to estimating equations
est_mat = Hmk_mat_c[Hmk_mat_c$k>=Hmk_mat_c$m & Hmk_mat_c$m%in%c(1995:1998,2000,2001) & Hmk_mat_c$m<=Hmk_mat_c$first_dereg_year2,]

#add Q, R, D columns
for(i in c(1995:1998,2000,2001)){
  for(j in i:2005){
    est_mat[,paste0('Q_',i,'_',j)] = as.numeric(est_mat$m==i & est_mat$k==j)
  }
}
est_mat$Q_amt = est_mat$L_amt

for(i in c(1995:1998,2000,2001)){
  for(j in i:2005){
    est_mat[,paste0('R_',i,'_',j)] = as.numeric(est_mat$first_dereg_year2==i & est_mat$k==j)
  }
}
est_mat$R_amt = est_mat$L_gamma_amt

for(i in c(1995:1998,2000,2001)){
  for(j in i:2005){
    est_mat[,paste0('D_',i,'_',j,'_amt')] = as.numeric(est_mat$m==i & est_mat$k==j)*est_mat$L_amt
  }
}

#first term in closed form linear solution
# A = colSums((est_mat$k_hpi-est_mat$lag_k_hpi)*(est_mat$A_any_m-est_mat$p_any)*est_mat[,grep('Q_',names(est_mat))])
A_dr = colSums((est_mat$k_hpi-est_mat$lag_k_hpi)*cbind((est_mat$A_any_m-est_mat$p_any)*est_mat[,grep('Q_',names(est_mat))],est_mat[,grep('D_',names(est_mat))]))
#k_hpi is Y_k, lag_k_hpi is Y_k-1
#A_any_m is treatment at time m
#p_any is estimated expected treatment at time m output of treatment model
#grep('Q_',names(est_mat)) gives filler/garbage columns
#grep('D_',names(est_mat)) gives nuisance outcome model covariate columns

#make second term of closed form linear solution formula
R_inds=grep('R_',names(est_mat))
D_inds=grep('D_',names(est_mat))
Q_inds=grep('Q_',names(est_mat))
p = length(Q_inds)
#Vimk in paper (17)
get_S = function(ind){
  if(est_mat$first_dereg_year2[ind]>est_mat$k[ind]){
    return(rep(0,p))
  }else{
    return(as.vector(est_mat[ind,R_inds]))
  }
}

#Vimk-1 in paper (17)
get_S_lag = function(ind){
  if(est_mat$first_dereg_year2[ind]>(est_mat$k[ind]-1)){
    return(rep(0,p))
  }else{
    return(as.vector(est_mat[ind-1,R_inds]))
  }
}

get_D = function(ind){
  as.vector(est_mat[ind,D_inds])
}

get_Q = function(ind){
  as.vector(est_mat[ind,Q_inds])
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
  B_list_dr[[i]] = get_B_dr(S_diff_list[[i]],D_list[[i]],Q_list[[i]],est_mat$A_any_m[i]-est_mat$p_any[i])
}
B_dr = Reduce('+',B_list_dr)
#this is the estimate of the blip function parameter
#do (17) from paper
psi_hat_dr = as.vector(A_dr%*%solve(B_dr))

#generate estimates of derived quantities using blip function parameter
cd_estimands = rep(NA,8)
start_years = c(1995:1998,2000,2001)
start_psi_inds = c(1,12,22,31,39,45)
for(i in 0:7){
    gammas = as.vector(as.matrix(est_mat[,R_inds])%*%psi_hat_dr[1:50])
    gammas = gammas[est_mat$m==est_mat$first_dereg_year2 & (est_mat$k-est_mat$m)==i]
    cd_estimands[i+1] = mean(gammas)
}

est_mat$gammas=ifelse(est_mat$k<est_mat$first_dereg_year2,0,as.matrix(est_mat[,R_inds])%*%psi_hat_dr[1:50])
est_mat$H = est_mat$k_hpi-est_mat$gammas
no_treat_hat = rep(NA,11)
years=1995:2005
for(i in 1:11){
  no_treat_hat[i] = mean(est_mat$H[est_mat$m==1995 & est_mat$k==years[i]])
}
observed_outcome_trajectory = aggregate(data$log_hpi[data$year>=1995],by=list(data$k[data$year>=1995]),mean)[,2]
observed_outcome_trajectory-no_treat_hat

B_list = vector(mode='list',length=length(S_diff_list))
for(i in 1:length(S_diff_list)){
  B_list[[i]] = get_B(S_diff_list[[i]],Q_list[[i]],est_mat$A_any_m[i]-est_mat$p_any[i])
}
B = Reduce('+',B_list)
psi_hat = as.vector(A%*%solve(B))

################
#BOOTSTRAP
################

nboot = 200
psi_hat_dr_boot = vector(mode='list',length=nboot)
no_treat_hat_boot = vector(mode='list',length=nboot)
observed_outcome_trajectory_boot = vector(mode='list',length=nboot)
cd_estimands_boot = vector(mode='list',length=nboot)
counties = unique(Hmk_mat_c$county)
data_inds = lapply(1:length(counties),function(i)((i-1)*12+1):(i*12))
b=1
count=1
while(count <= nboot){
  set.seed(b)
  samp_counties = sample(1:length(counties),length(counties),replace=T)
  mod_inds = unlist(data_inds[samp_counties])
  data_boot = data[mod_inds,]
  data_boot$county = rep(1:length(counties),each=12)
  if(sum(data_boot$first_dereg_year3==1995)>0){
    Hmk_mat_c_boot = data.frame(cbind(county=rep(1:length(counties),each=77),m=rep(rep(rep(1995:2005,c(12:2)),length(counties))),
                                      k=rep(c(1994:2005,1995:2005,1996:2005,1997:2005,1998:2005,1999:2005,2000:2005,2001:2005,2002:2005,2003:2005,2004:2005),length(county_keepers))))
    first_dereg_year_boot = aggregate(data_boot$inter_bra,by=list(data_boot$county),function(x) min(which(x>0)))
    names(first_dereg_year_boot) = c('county','first_dereg_year')
    
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,first_dereg_year_boot)
    Hmk_mat_c_boot$first_dereg_year2 = Hmk_mat_c_boot$first_dereg_year + 1994 - 1
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","first_dereg_year2","inter_bra")],all.x = T)
    
    Hmk_mat_c_boot$inter_bra[is.na(Hmk_mat_c_boot$inter_bra)]=0
    Hmk_mat_c_boot = Hmk_mat_c_boot[order(Hmk_mat_c_boot$county,Hmk_mat_c_boot$m,Hmk_mat_c_boot$k),]
    
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","lag_hpi","first_dereg_year2")],all.x = T)
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","lag_amt","first_dereg_year2")],all.x = T)
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","k_amt","k")],all.x = T)
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","k_hpi","k")],all.x = T)
    
    # A_mod_boot = glm(inter_bra~year+I(year^2)+I(year^3)+lag_hpi+I(lag_hpi^2),data=data_boot[!data_boot$past_A&data_boot$year>1994,],family="poisson")
    A_mod_any_boot = glm(any_A~factor(year)*lag_amt+I(lag_amt^2),data=data_boot[!data_boot$past_A&data_boot$year>1994,],family="binomial")
    data_boot$p_any = NA
    data_boot$p_any[data_boot$year>1994] = predict(A_mod_any_boot,newdata = data_boot[data_boot$year>1994,],type='response')
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","m",'p_any')],all.x=T)
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","m",'A_m')],all.x=T)
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","m",'L_amt','L_hpi')],all.x=T)
    
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[data_boot$year==data_boot$first_dereg_year3,c("county","year",'L_gamma_amt','L_gamma_hpi')],all.x=T)
    Hmk_mat_c_boot = merge(Hmk_mat_c_boot,data_boot[,c("county","k",'lag_k_amt','lag_k_hpi')],all.x=T)
    
    Hmk_mat_c_boot = Hmk_mat_c_boot[order(Hmk_mat_c_boot$county,Hmk_mat_c_boot$m,Hmk_mat_c_boot$k),]
    Hmk_mat_c_boot$A_any_m = as.numeric(Hmk_mat_c_boot$A_m>0)
    Hmk_mat_c_boot$inter_bra_any = as.numeric(Hmk_mat_c_boot$inter_bra>0)
    
    est_mat_boot = Hmk_mat_c_boot[Hmk_mat_c_boot$k>=Hmk_mat_c_boot$m & Hmk_mat_c_boot$m%in%c(1995:1998,2000,2001) & Hmk_mat_c_boot$m<=Hmk_mat_c_boot$first_dereg_year2,]
    for(i in c(1995:1998,2000,2001)){
      for(j in i:2005){
        est_mat_boot[,paste0('Q_',i,'_',j)] = as.numeric(est_mat_boot$m==i & est_mat_boot$k==j)
      }
    }
    est_mat_boot$Q_amt = est_mat_boot$L_amt
    
    for(i in c(1995:1998,2000,2001)){
      for(j in i:2005){
        est_mat_boot[,paste0('R_',i,'_',j)] = as.numeric(est_mat_boot$first_dereg_year2==i & est_mat_boot$k==j)
      }
    }
    est_mat_boot$R_amt = est_mat_boot$L_gamma_amt
    
    for(i in c(1995:1998,2000,2001)){
      for(j in i:2005){
        est_mat_boot[,paste0('D_',i,'_',j,'_amt')] = as.numeric(est_mat_boot$m==i & est_mat_boot$k==j)*est_mat_boot$L_amt
      }
    }
    
    A_dr_boot = colSums((est_mat_boot$k_hpi-est_mat_boot$lag_k_hpi)*cbind((est_mat_boot$A_any_m-est_mat_boot$p_any)*est_mat_boot[,grep('Q_',names(est_mat_boot))],est_mat_boot[,grep('D_',names(est_mat_boot))]))
    
    R_inds=grep('R_',names(est_mat_boot))
    D_inds=grep('D_',names(est_mat_boot))
    Q_inds=grep('Q_',names(est_mat_boot))
    p = length(Q_inds)
    get_S = function(ind){
      if(est_mat_boot$first_dereg_year2[ind]>est_mat_boot$k[ind]){
        return(rep(0,p))
      }else{
        return(as.vector(est_mat_boot[ind,R_inds]))
      }
    }
    
    get_S_lag = function(ind){
      if(est_mat_boot$first_dereg_year2[ind]>(est_mat_boot$k[ind]-1)){
        return(rep(0,p))
      }else{
        return(as.vector(est_mat_boot[ind-1,R_inds]))
      }
    }
    
    get_D = function(ind){
      as.vector(est_mat_boot[ind,D_inds])
    }
    
    get_Q = function(ind){
      as.vector(est_mat_boot[ind,Q_inds])
    }
    
    S_diff_list_boot = lapply(1:nrow(est_mat_boot),function(i) get_S(i)-get_S_lag(i))
    D_list_boot = lapply(1:nrow(est_mat_boot),get_D)
    Q_list_boot = lapply(1:nrow(est_mat_boot),get_Q)
    
    B_list_dr_boot = vector(mode='list',length=length(S_diff_list_boot))
    for(i in 1:length(S_diff_list_boot)){
      B_list_dr_boot[[i]] = get_B_dr(S_diff_list_boot[[i]],D_list_boot[[i]],Q_list_boot[[i]],est_mat_boot$A_any_m[i]-est_mat_boot$p_any[i])
    }
    B_dr_boot = Reduce('+',B_list_dr_boot)
    psi_hat_dr_b = as.vector(A_dr_boot%*%solve(B_dr_boot))
    psi_hat_dr_boot[[count]] = psi_hat_dr_b
    
    cd_estimands_b = rep(NA,8)
    start_years = c(1995:1998,2000,2001)
    for(i in 0:7){
      gammas = as.vector(as.matrix(est_mat_boot[,R_inds])%*%psi_hat_dr_b[1:50])
      gammas = gammas[est_mat_boot$m==est_mat_boot$first_dereg_year2 & (est_mat_boot$k-est_mat_boot$m)==i]
      cd_estimands_b[i+1] = mean(gammas)
    }
    cd_estimands_boot[[count]] = cd_estimands_b
    
    est_mat_boot$gammas=ifelse(est_mat_boot$k<est_mat_boot$first_dereg_year2,0,as.matrix(est_mat_boot[,R_inds])%*%psi_hat_dr_b[1:50])
    est_mat_boot$H = est_mat_boot$k_hpi-est_mat_boot$gammas
    no_treat_hat_b = rep(NA,11)
    years=1995:2005
    for(i in 1:11){
      no_treat_hat_b[i] = mean(est_mat_boot$H[est_mat_boot$m==1995 & est_mat_boot$k==years[i]])
    }
    no_treat_hat_boot[[count]] = no_treat_hat_b
    observed_outcome_trajectory_boot[[count]] = aggregate(data_boot$log_hpi[data_boot$year>=1995],by=list(data_boot$k[data$year>=1995]),mean)[,2]
    count=count+1
    b=b+1
  }else{
    b=b+1
  }
}

save(psi_hat_dr_boot,no_treat_hat_boot,cd_estimands_boot,observed_outcome_trajectory_boot,
     psi_hat_dr,no_treat_hat,observed_outcome_trajectory,cd_estimands,
     file="did_snmm_coarse_hpi_Lamt_results.RData")

load("did_snmm_coarse_hpi_Lamt_results.RData")
no_treat_hat_boot_mat = data.frame(no_treat_hat_boot)
observed_outcome_trajectory_boot_mat = data.frame(observed_outcome_trajectory_boot)
no_treat_effect_boot_mat = no_treat_hat_boot_mat - observed_outcome_trajectory_boot_mat
no_treat_effect_sds = apply(no_treat_effect_boot_mat,1,sd)
uppers = no_treat_hat - observed_outcome_trajectory + 1.96*no_treat_effect_sds
lowers = no_treat_hat - observed_outcome_trajectory - 1.96*no_treat_effect_sds
years=1995:2005
pdf("did_snmm_coarse_hpi_Lamt_plot.pdf")
plot(years,no_treat_hat - observed_outcome_trajectory,type='l',ylim=c(min(lowers),max(uppers)),main="Effect on log(HPI) of Removing All Deregulation vs Usual Care",ylab="Effect on log(HPI)",xlab="Year")
points(years,lowers,type='l',lty=2)
points(years,uppers,type='l',lty=2)
curve(0*x,from=1990,to=3000,col=2,add=T)
dev.off()

psi_hat_dr_mat = data.frame(psi_hat_dr_boot)
pdf('did_snmm_coarse_psihat_amt.pdf')
hist(unlist(psi_hat_dr_mat[50,]),main="Bootstrap Distribution of Mortgage Volume Coefficient",xlab=expression(beta))
points(psi_hat_dr[50],0,col=2)
dev.off()

psi_hat_dr[50]
psi_hat_dr[50] + 1.96*sd(psi_hat_dr_mat[50,])
psi_hat_dr[50] - 1.96*sd(psi_hat_dr_mat[50,])


cd_estimands_boot_mat = data.frame(cd_estimands_boot)
cd_estimands_sds = apply(cd_estimands_boot_mat,1,sd)
uppers = cd_estimands + 1.96*cd_estimands_sds
lowers = cd_estimands - 1.96*cd_estimands_sds
pdf("did_snmm_coarse_cd_estimands_Lamt_plot.pdf")
plot(0:7,cd_estimands,type='l',ylim=c(min(lowers),max(uppers)),main="Weighted Average ITT-ETTs on log(HPI)",ylab="Effect on log(HPI)",xlab="Year After Treatment")
points(0:7,lowers,type='l',lty=2)
points(0:7,uppers,type='l',lty=2)
curve(0*x,from=0,to=3000,col=2,add=T)
dev.off()
