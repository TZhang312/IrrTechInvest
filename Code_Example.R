#------------------------------------------------------------------------------------
#--------------Exmaple of core functions---------------------------------------------
#------------------------------------------------------------------------------------

#--(1) Panel data model analysis
regressionModel <- function (ds, fml) {
  
  ds$water <- ds$irr + ds$prcp
  ds$water2 <- ds$water^2
  ds$nitr2   <- ds$nitr^2
  ds$year   <- ds$year^2
  ds$LogY   <- log(ds$yield)
  reg.model <- lm(formula=fml, na.action = na.omit, data=ds)  
  return (reg.model)
}



#--(2) Estimate future climate and CO2 impacts relative to baseline period
estYld<-function (reg.model, ds.futClim, ssp.i, period.i, ds, CO2_A, CO2_B) {
  
  futClim <- subset (ds.futClim, ds.futClim$ssp==ssp.i, ds.futClim$period==period.i)
  
  yld.baseline <- exp(predict (reg.model))
  yld.futClim <- exp(predict (reg.model, futClim))
  pct.yld.clim <- (yld.futClim - yld.baseline)/yld.baseline * 100
  
  impact.yld.co2.baseline <-(ds$CO2_baseline - 360)/(ds$CO2_baseline - 360 +CO2_A) * CO2_B /100
  impact.yld.co2.futClim <-(ds$CO2_futClim - 360)/(ds$CO2_futClim - 360 +CO2_A) * CO2_B /100
  pct.yld.co2.futClim_relative_co2.baseline <- (impact.yld.co2.futClim - impact.yld.co2.baseline)/(1+impact.yld.co2.baseline) * 100
  pct.yld.climCo2_ssp <- pct.yld.clim + pct.yld.co2.futClim_relative_co2.baseline
  
  return (pct.yld.climCo2_ssp)
}


#---(3) Incoporate technological benifits factors into model and estimate impacts
incorpTechBenifit <- function (ds, ds.futClim, ssp.i, period.i) {
  
  ds <- ds %>%
    left_join(adptFac, by=c("region")) %>% 
    mutate(YLD_ir_Tech=YLD_ir*YLD_chge_fac) %>%
    mutate(YLD_ct_Tech=YLD_ir_Tech*fri_irrA + YLD_rf*pct_rfA) %>%
    mutate(irr_irrA_Tech_field=irr_irrA * IRR_chge_fac)   %>%
    mutate(irr_ct_Tech_field = irr_irrA_Tech_field*fri_irrA)    %>%
    mutate(logIrrP_baseline = log(YLD_ct)/irr)           %>%
    mutate(irr_irrA_Tech_fieldCov = irr_irrA * IRR_chge_fac * COV_chge_fac)   %>%
    mutate(irr_ct_Tech_fieldCov   = irr_irrA_Tech_fieldCov*fri_irrA)    %>%
    mutate(logIrrP_tech_fieldCov    = log(YLD_ct_Tech)/irr_ct_Tech_fieldCov)  %>%
    mutate(IrrP_fieldCov_fac     = logIrrP_tech_fieldCov/logIrrP_baseline) %>%
    mutate (pct_irr_Tech =  (irr_ct_Tech_fieldCov -  irr)/ irr*100 ) 
  
  futClim <- subset (ds.futClim, ds.futClim$ssp==ssp.i, ds.futClim$period==period.i)
  futClim$irr_Tech <- futClim$irr * futClim$IrrP_fieldCov_fac
  futClim$water <- futClim$irr_Tech + futClim$prcp
  yld.futClim_Tech <- exp(predict(reg.model, futClim))
  pct.yld.futClim_Tech <- (yld.futClim_Tech - yld.baseline)/yld.baseline * 100
  pct.yld.climCo2_Tech_ssp <- pct.yld.futClim_Tech + pct.yld.co2.futClim_ssp_relative_co2.baseline
  
  return (pct.yld.climCo2_Tech_ssp)
}


#--- (4) Simulation with irrigation technology upgrades and irrigatgion technology expansion 
simTechEffect <- function (ds.irr, ds.yld) {

  ds.irr.upr.expd <- ds.irr %>%
    mutate (irr_irrA= irr/fri_irrA) %>%
    mutate(fri_rfA_New=fri_rfA - fri_irrA.increase) %>%
    mutate(fri_rfA_New=ifelse (fri_rfA_New < 0, 0, fri_rfA_New)) %>%
    mutate(fri_irrA.increase_adj = fri_rfA - fri_rfA_New) %>%
    mutate (fri_irrA_New= fri_irrA + fri_irrA.increase_adj) %>%
    mutate (irr_ct_Tech_adj = irr_ct_Tech_upr + irr_ct_Tech_expd/fri_irrA *irrA_New ) %>%
    mutate (irr_irrA_Tech_adj = irr_ct_Tech_adj / fri_irrA_New )%>%
    mutate (pct_irr_Tech_New =  (irr_ct_Tech_adj -  irr)/ irr*100 ) %>%
    mutate (pct_irr_Tech_irrA_New =  (irr_irrA_Tech_adj -  irr_irrA)/ irr_irrA*100 ) 
  
  
  ds.yld.upr.expd <- ds.yld %>%
    mutate(fri_rfA_New=pct_rfA - fri_irrA.increase) %>%
    mutate(fri_rfA_New=ifelse (fri_rfA_New < 0, 0, fri_rfA_New)) %>%
    mutate(fri_irrA.increase_adj = fri_rfA - fri_rfA_New) %>%
    mutate(yld_all_adj = yld_rf * fri_rfA_New + yld_irr* fri_irrA + yld_irr_exp * fri_irrA.increase_adj) %>%
    mutate (pct_yld_all =  (yld_all_adj -  yld_baseline)/ yld_baseline*100 )
  
  rsl <- list (ds.irr.upr.expd, 
			   ds.yld.upr.expd)
  return (rsl)
}


#-----(5) Irrigation investment portifilo 
optComb<- function (ds.irr, ds.irr.crp, ds.irrV.crp,  
                    ds.yld, ds.prod.crp, ds.prodV.crp,
                    ds.cost, ds.cost.crp, 
                    cn.prodV.crp, cn.wuV.crp, cn.techV.crp) { 
  

  ds.irrV<-ds.irr
  ds.yldV<-ds.yld
  
  tech.set <- unique(ds.irr$tech_expd_irrA)
  reg.set  <- names(ds.irr[,c(2:ncol(ds.irr))])
  
  tech.n<- length(tech.set)
  reg.n<-  length(reg.set)
  comb.n <- tech.n^reg.n
  

  ds.irr$NC<-ds.irr$N * ds.irr.crp[ds.irr.crp$region3=="NC",]$pct_water
  ds.irr$NE<-ds.irr$NE * ds.irr.crp[ds.irr.crp$region3=="NE",]$pct_water
  ds.irr$NW<-ds.irr$NW * ds.irr.crp[ds.irr.crp$region3=="NW",]$pct_water
  ds.irr$S<-ds.irr$S * ds.irr.crp[ds.irr.crp$region3=="S",]$pct_water
  ds.irr$SW<-ds.irr$SW * ds.irr.crp[ds.irr.crp$region3=="SW",]$pct_water
  ds.irr$YZ<-ds.irr$YZ * ds.irr.crp[ds.irr.crp$region3=="YZ",]$pct_water
    
  ds.irrV$NC<-ds.irrV$N * ds.irrV.crp[ds.irrV.crp$region3=="NC",]$pct_irrV
  ds.irrV$NE<-ds.irrV$NE * ds.irrV.crp[ds.irrV.crp$region3=="NE",]$pct_irrV
  ds.irrV$NW<-ds.irrV$NW * ds.irrV.crp[ds.irrV.crp$region3=="NW",]$pct_irrV
  ds.irrV$S<-ds.irrV$S * ds.irrV.crp[ds.irrV.crp$region3=="S",]$pct_irrV
  ds.irrV$SW<-ds.irrV$SW * ds.irrV.crp[ds.irrV.crp$region3=="SW",]$pct_irrV
  ds.irrV$YZ<-ds.irrV$YZ * ds.irrV.crp[ds.irrV.crp$region3=="YZ",]$pct_irrV
    
  ds.yld$NC<-ds.yld$N * ds.prod.crp[ds.prod.crp$region3=="NC",]$pct_prod
  ds.yld$NE<-ds.yld$NE * ds.prod.crp[ds.prod.crp$region3=="NE",]$pct_prod
  ds.yld$NW<-ds.yld$NW * ds.prod.crp[ds.prod.crp$region3=="NW",]$pct_prod
  ds.yld$S<-ds.yld$S * ds.prod.crp[ds.prod.crp$region3=="S",]$pct_prod
  ds.yld$SW<-ds.yld$SW * ds.prod.crp[ds.prod.crp$region3=="SW",]$pct_prod
  ds.yld$YZ<-ds.yld$YZ * ds.prod.crp[ds.prod.crp$region3=="YZ",]$pct_prod
    
  ds.yldV$NC<-ds.yldV$N * ds.prodV.crp[ds.prodV.crp$region3=="NC",]$pct_prodV
  ds.yldV$NE<-ds.yldV$NE * ds.prodV.crp[ds.prodV.crp$region3=="NE",]$pct_prodV
  ds.yldV$NW<-ds.yldV$NW * ds.prodV.crp[ds.prodV.crp$region3=="NW",]$pct_prodV
  ds.yldV$S<-ds.yldV$S * ds.prodV.crp[ds.prodV.crp$region3=="S",]$pct_prodV
  ds.yldV$SW<-ds.yldV$SW * ds.prodV.crp[ds.prodV.crp$region3=="SW",]$pct_prodV
  ds.yldV$YZ<-ds.yldV$YZ * ds.prodV.crp[ds.prodV.crp$region3=="YZ",]$pct_prodV
    
    
  ds.cost$NC<-ds.cost$N * ds.cost.crp[ds.cost.crp$region3=="NC",]$pct_techCost
  ds.cost$NE<-ds.cost$NE * ds.cost.crp[ds.cost.crp$region3=="NE",]$pct_techCost
  ds.cost$NW<-ds.cost$NW * ds.cost.crp[ds.cost.crp$region3=="NW",]$pct_techCost
  ds.cost$SW<-ds.cost$SW * ds.cost.crp[ds.cost.crp$region3=="SW",]$pct_techCost
  ds.cost$YZ<-ds.cost$YZ * ds.cost.crp[ds.cost.crp$region3=="YZ",]$pct_techCost
  ds.cost$S<-ds.cost$S * ds.cost.crp[ds.cost.crp$region3=="S",]$pct_techCost
  

  ds.vct <- ds.irr[,c(2:ncol(ds.irr))]
  for (j in 1:ncol(ds.vct)) {
    
    if (j==1) {
      vct.irr.1 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      
    } else {
      vct.irr.2 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      vct.irr.1 <- vct.irr.1 + vct.irr.2
      rm(vct.irr.2)
      gc()
    }
  }
  
  ds.vct <- ds.irrV[,c(2:ncol(ds.irrV))]
  for (j in 1:ncol(ds.vct)) {
    
    if (j==1) {
      vct.irrV.1 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      
    } else {
      vct.irrV.2 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      vct.irrV.1 <- vct.irrV.1 + vct.irrV.2
      rm(vct.irrV.2)
      gc()
    }
  }
  
  ds.vct <- ds.yld[,c(2:ncol(ds.yld))]
  for (j in 1:ncol(ds.vct)) {
    
    if (j==1) {
      vct.yld.1 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      
    } else {
      vct.yld.2 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      vct.yld.1 <- vct.yld.1 + vct.yld.2
      rm(vct.yld.2)
      gc()
    }
  }
  
  
  
  ds.vct <- ds.yld[,c(2:ncol(ds.yldV))]
  for (j in 1:ncol(ds.vct)) {
    
    if (j==1) {
      vct.yldV.1 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      
    } else {
      vct.yldV.2 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      vct.yldV.1 <- vct.yldV.1 + vct.yldV.2
      rm(vct.yldV.2)
      gc()
    }
  }
  
  ds.vct <- ds.cost[,c(2:ncol(ds.cost))]
  for (j in 1:ncol(ds.vct)) {
    
    if (j==1) {
      vct.cost.1 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      
    } else {
      vct.cost.2 <- rep(rep(ds.vct[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
      vct.cost.1 <- vct.cost.1 + vct.cost.2
      rm(vct.cost.2)
      gc()
    }
  }
  
  rsl.comb <- data.table(pct_irr=vct.irr.1, pct_yld=vct.yld.1, pct_yldV=vct.yldV.1, pct_irrV=vct.irrV.1, pct_costV=vct.cost.1)
  rsl.comb$roi <- (rsl.comb$pct_yldV *  cn.prodV.crp$value_prod - rsl.comb$pct_irrV *  cn.wuV.crp$value_water )/(rsl.comb$pct_costV * cn.techV.crp$techCost)
  
  ds.target <- rsl.comb %>%
    filter (pct_yld>=0 & pct_irr<=tgt.wu)
  
  # -f there is no any case
  if (nrow(ds.target)==0) {  
    rm(rsl.comb)
    gc()
    return (NA)
  }


  max.roi <- max (ds.target$roi)
  

  pos<-which( rsl.comb$pct_yld>=0
                  & rsl.comb$pct_irr<=tgt.wu
                  & rsl.comb$roi==max.roi)
  
  title.NC<-paste(reg.set[1],tech.set ,sep=":")
  title.NE<-paste(reg.set[2],tech.set ,sep=":")
  title.NW<-paste(reg.set[3],tech.set ,sep=":")
  title.S <-paste(reg.set[4],tech.set ,sep=":")
  title.SW<-paste(reg.set[5],tech.set ,sep=":")
  title.YZ<-paste(reg.set[6],tech.set ,sep=":")
   
  ds.title<-data.frame(cbind (NC=title.NC, NE=title.NE, NW=title.NW, S=title.S,SW=title.SW, YZ=title.YZ) )
  
  best.tech <- c()
  for (j in 1:ncol(ds.title)) {
    
    vct.title <- rep(rep(ds.title[,j], times=comb.n/(tech.n^j)), each=comb.n/(tech.n^(ncol(ds.vct)+1-j)))
    lab <- vct.title[pos]
    rm(vct.title)
    gc()
    best.tech <- c(best.tech, lab)
    
  }
  
  best.comb <- rsl.comb[pos]
  rsl.final <- data.frame(NC=best.tech[1], NE=best.tech[2], NW=best.tech[3],S=best.tech[4], 
                          SW=best.tech[5], YZ=best.tech[6], 
                          pct_irr=best.comb$pct_irr, pct_yld=best.comb$pct_yld,
                          pct_yldV=best.comb$pct_yldV, pct_irrV=best.comb$pct_irrV, 
                          pct_cost=best.comb$pct_cost, roi=best.comb$roi
  )
  
  return (rsl.final)
  
}


optIrrInvest <- function (comb.i, ds.crp) {
  
  tech.comb.i<-unlist(strsplit(comb.i,split="[+]"))
  
  ds.crp.comb <- ds.crp %>%
    filter (tech_expd %in% tech.comb.i)
  
  ds.irr  <- ds.crp.comb %>%
    reshape2::dcast(tech_expd_irrA ~ region, value.var = "pct_irr_Tech_New") 
  
  ds.yld  <- ds.crp.comb %>%
    reshape2::dcast(tech_expd_irrA ~ region, value.var = "pct_yld_clim_co2") 
  
  ds.cost  <- ds.crp.comb %>%
    reshape2::dcast(tech_expd_irrA ~ region, value.var = "pct_techCost") 
  
  
  rsl<-optComb (ds.irr, ds.irr.crp, ds.irrV.crp,  
                ds.yld, ds.prod.crp, ds.prodV.crp,
                ds.cost, ds.cost.crp, 
                cn.prodV.crp, cn.wuV.crp, cn.techV.crp)
				
  return (rsl)
  
}


