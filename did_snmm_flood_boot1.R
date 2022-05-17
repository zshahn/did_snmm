library(nleqslv)
nboot = 25
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
save(psi_hat,psi_hat_boot,file='flood_results_boot_1_25.RData')
