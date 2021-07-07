library(ggplot2) # for graph
require(deSolve) # for the "ode" function
library(dplyr)
#library(bbmle) # for ML method  install.packages("bbmle", repos="http://R-Forge.R-project.org")



#Equation#######################
# ver6 and after: D is not death, but severe

sir_1 <- function(rt_ini, rt_rebLength, rt_reb, rt_var, rt_oly, soe_th, soe_eff, soe_time, 
                  gamma, deltaY_temp, deltaA_temp, deltaL, del_var,
                  hosp_rate, hosp_plus_len, 
                  vaceff_inf, vaceff_death_temp, 
                  Sy0, Iy0, I2y0, Ry0, Dy0, Dry0,   # i2 for longer hosp
                  Psyv0, Piyv0, Pi2yv0, Pryv0, Pdyv0, Pdryv0, 
                  Syv0, Iyv0, I2yv0, Ryv0, Dyv0, Dryv0, 
                  Sa0, Ia0, I2a0, Ra0, Da0, Dra0, 
                  Psav0, Piav0, Pi2av0, Prav0, Pdav0, Pdrav0, 
                  Sav0, Iav0, I2av0, Rav0, Dav0, Drav0,
                  Ihy0, Pihyv0, Ihyv0, Iha0, Pihav0, Ihav0,# h for henni-kabu
                  Vy0, Va0, 
                  Soeset0,
                  times) {
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      jinryu_to_rt <- 0.025
      
      if (time>soe_time) {
        rt <- min(rt_ini + soe_time*rt_reb, rt_ini + rt_rebLength*rt_reb)
      } else if (time<rt_rebLength) {
        rt <- rt_ini + time*rt_reb
      } else {
        rt <- rt_ini + rt_rebLength*rt_reb
      }
      

      
      
      
      deltaY <- deltaY_temp
      deltaA <- deltaA_temp
      

      
      #olympic      
      if (44<time && time<60) {  #time0=June09, time44=july23, time60=Aug08
        rt <- (rt_oly*jinryu_to_rt)+rt
      }
      
      
      #soe_th/eff
      soe <- (Iy+Iyv+Ia+Iav)/5
#      
#      if (Soeset>0) {
#        rt <- rt * soe_eff
#      }
      
      
      #soe_time
      if (time>soe_time+10) {
        rt <- rt * soe_eff * soe_eff
      } else if (time>soe_time) {
        rt <- rt * soe_eff
      }
        
      
      

      #var
      rt_varnow <- rt * rt_var
      
      
      #vaccines
      if (time<12+7) {   #time0=June09, time12=june21
        dailyDose <- 40000
      } else {
        dailyDose <- 50000
      }
      
      
      adult_vaccination_max <- 0.8 * (320*(10^5))
      young_vaccination_max <- 0.4 * (1070*(10^5))
      
      adult_vaccinated <- Psav+Piav+Pi2av+Prav+Pdav+Pdrav+Sav+Iav+I2av+Rav+Dav+Drav+Va
      young_vaccinated <- Psyv+Piyv+Pi2yv+Prav+Pdyv+Pdryv+Syv+Iyv+I2yv+Ryv+Dyv+Dryv+Vy

#      adult_vaccination <- 0
#      young_vaccination <- 0
              
      if (time<12+7) {      #time0=June09, time19=june28
        adult_vaccination <- ifelse(adult_vaccinated>=adult_vaccination_max,0,min(dailyDose, (adult_vaccination_max - adult_vaccinated)))
        young_vaccination <- ifelse(young_vaccinated>=young_vaccination_max,0,min(dailyDose - adult_vaccination, (young_vaccination_max - young_vaccinated)))
      } else {
        adult_vaccination <- ifelse(adult_vaccinated>=adult_vaccination_max,0,min(40000, (adult_vaccination_max - adult_vaccinated))) #50000
        young_vaccination <- ifelse(young_vaccinated>=young_vaccination_max,0,min(dailyDose - adult_vaccination, (young_vaccination_max - young_vaccinated)))
      }
      
      
      vaccine_delay <- 1/28
      vaceff_death <- (1-vaceff_death_temp) / (1-vaceff_inf)
      
      
      
      pop <-  Sy+Iy+I2y+Ry+Dy+Dry+Psyv+Piyv+Pi2yv+Pryv+Pdyv+Pdryv+Syv+Iyv+I2yv+Ryv+Dyv+Dryv+Sa+Ia+I2a+Ra+Da+Dra+Psav+Piav+Pi2av+Prav+Pdav+Pdrav+Sav+Iav+I2av+Rav+Dav+Drav+Ihy+Pihyv+Ihyv+Iha+Pihav+Ihav+Vy+Va
      
      infy <- Iy+Piyv+Iyv
      infa <- Ia+Piav+Iav

      infy_var <- Ihy+Pihyv+Ihyv
      infa_var <- Iha+Pihav+Ihav
      
      
            
      dSy <-   -rt / pop * gamma * infy * (1390/1925*1.5) * Sy -rt / pop * gamma * infa * (1390/1925) * Sy -rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Sy -rt_varnow / pop * gamma * infa_var * (1390/1925) * Sy - young_vaccination
      dPsyv <- -rt / pop * gamma * infy * (1390/1925*1.5) * Psyv -rt / pop * gamma * infa * (1390/1925) * Psyv -rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Psyv -rt_varnow / pop * gamma * infa_var * (1390/1925) * Psyv + young_vaccination - Psyv * vaccine_delay
      dSyv <-  -rt / pop * gamma * infy * (1390/1925*1.5) * Syv -rt / pop * gamma * infa * (1390/1925) * Syv -rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Syv -rt_varnow / pop * gamma * infa_var * (1390/1925) * Syv + Psyv * vaccine_delay * (1-vaceff_inf)

      dSa <-  -rt / pop * gamma * infy * (1390/1925*1.5) * Sa -rt / pop * gamma * infa * (1390/1925) * Sa -rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Sa -rt_varnow / pop * gamma * infa_var * (1390/1925) * Sa - adult_vaccination
      dPsav <- -rt / pop * gamma * infy * (1390/1925*1.5) * Psav -rt / pop * gamma * infa * (1390/1925) * Psav -rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Psav -rt_varnow / pop * gamma * infa_var * (1390/1925) * Psav + adult_vaccination - Psav * vaccine_delay
      dSav <- -rt / pop * gamma * infy * (1390/1925*1.5) * Sav -rt / pop * gamma * infa * (1390/1925) * Sav -rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Sav -rt_varnow / pop * gamma * infa_var * (1390/1925) * Sav + Psav * vaccine_delay * (1-vaceff_inf)
      

      dIy <-   rt / pop * gamma * infy * (1390/1925*1.5) * Sy +rt / pop * gamma * infa * (1390/1925) * Sy  - gamma * Iy
      dPiyv <- rt / pop * gamma * infy * (1390/1925*1.5) * Psyv +rt / pop * gamma * infa * (1390/1925) * Psyv - gamma * Piyv 
      dIyv <-  rt / pop * gamma * infy * (1390/1925*1.5) * Syv +rt / pop * gamma * infa * (1390/1925) * Syv - gamma * Iyv

      dIhy <-   rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Sy +rt_varnow / pop * gamma * infa_var * (1390/1925) * Sy  - gamma * Ihy
      dPihyv <- rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Psyv +rt_varnow / pop * gamma * infa_var * (1390/1925) * Psyv - gamma * Pihyv 
      dIhyv <-  rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Syv +rt_varnow / pop * gamma * infa_var * (1390/1925) * Syv - gamma * Ihyv
      
      dIa <-   rt / pop * gamma * infy * (1390/1925*1.5) * Sa +rt / pop * gamma * infa * (1390/1925) * Sa  - gamma * Ia
      dPiav <- rt / pop * gamma * infy * (1390/1925*1.5) * Psav +rt / pop * gamma * infa * (1390/1925) * Psav - gamma * Piav 
      dIav <- rt / pop * gamma * infy * (1390/1925*1.5) * Sav +rt / pop * gamma * infa * (1390/1925) * Sav - gamma * Iav

      dIha <-   rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Sa +rt_varnow / pop * gamma * infa_var * (1390/1925) * Sa  - gamma * Iha
      dPihav <- rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Psav +rt_varnow / pop * gamma * infa_var * (1390/1925) * Psav - gamma * Pihav 
      dIhav <- rt_varnow / pop * gamma * infy_var * (1390/1925*1.5) * Sav +rt_varnow / pop * gamma * infa_var * (1390/1925) * Sav - gamma * Ihav
      
      
      
      dRy <-   gamma * Iy * (1-hosp_rate) + gamma * Ihy * (1-hosp_rate) + hosp_plus_len * I2y
      dPryv <- gamma * Piyv * (1-hosp_rate) + gamma * Pihyv * (1-hosp_rate) + hosp_plus_len * Pi2yv
      dRyv <-  gamma * Iyv * (1-hosp_rate*vaceff_death) + gamma * Ihyv * (1-hosp_rate*vaceff_death) + hosp_plus_len * I2yv
      
      dRa <-   gamma * Ia * (1-hosp_rate) + gamma * Iha * (1-hosp_rate) + hosp_plus_len * I2a
      dPrav <- gamma * Piav * (1-hosp_rate) + gamma * Pihav * (1-hosp_rate) + hosp_plus_len * Pi2av
      dRav <-  gamma * Iav * (1-hosp_rate*vaceff_death) + gamma * Ihav * (1-hosp_rate*vaceff_death) + hosp_plus_len * I2av

      
      
      dI2y <-   gamma * Iy * (hosp_rate - deltaY) + gamma * Ihy * (hosp_rate - deltaY*del_var) - hosp_plus_len * I2y
      dPi2yv <- gamma * Piyv * (hosp_rate - deltaY) + gamma * Pihyv * (hosp_rate - deltaY*del_var) - hosp_plus_len * Pi2yv
      dI2yv <-  gamma * Iyv * (hosp_rate*vaceff_death - deltaY*vaceff_death) + gamma * Ihyv * (hosp_rate*vaceff_death - deltaY*vaceff_death*del_var) - hosp_plus_len * I2yv
      
      dI2a <-   gamma * Ia * (hosp_rate - deltaA) + gamma * Iha * (hosp_rate - deltaA*del_var) - hosp_plus_len * I2a
      dPi2av <- gamma * Piav * (hosp_rate - deltaA) + gamma * Pihav * (hosp_rate - deltaA*del_var) - hosp_plus_len * Pi2av
      dI2av <-  gamma * Iav * (hosp_rate*vaceff_death - deltaA*vaceff_death) + gamma * Ihav * (hosp_rate*vaceff_death - deltaA*vaceff_death*del_var) - hosp_plus_len * I2av

      
      
      dDy <-   gamma * Iy * deltaY + gamma * Ihy * deltaY*del_var - deltaL * Dy
      dPdyv <- gamma * Piyv * deltaY + gamma * Pihyv * deltaY*del_var - deltaL * Pdyv
      dDyv <-  gamma * Iyv * deltaY*vaceff_death + gamma * Ihyv * deltaY*vaceff_death*del_var - deltaL * Dyv
      
      dDa <-   gamma * Ia * deltaA + gamma * Iha * deltaA*del_var - deltaL * Da
      dPdav <- gamma * Piav * deltaA + gamma * Pihav * deltaA*del_var - deltaL * Pdav
      dDav <-  gamma * Iav *deltaA*vaceff_death + gamma * Ihav *deltaA*vaceff_death*del_var - deltaL * Dav

            
      dDry <-  deltaL * Dy
      dPdryv <- deltaL * Pdyv
      dDryv <- deltaL * Dyv
      
      dDra <-  deltaL * Da
      dPdrav <- deltaL * Pdav
      dDrav <- deltaL * Dav      


      
      
      
      

      
      
      
      
      
      
      
      
      dVy <-   Psyv * vaccine_delay * vaceff_inf
      dVa <-   Psav * vaccine_delay * vaceff_inf
      

      
      if (soe > soe_th) {
        dSoeset <- 1
      } else {
        dSoeset <- 0
      }
      
      return(list(c(dSy, dIy, dI2y, dRy, dDy, dDry,
                    dPsyv, dPiyv, dPi2yv, dPryv, dPdyv, dPdryv,
                    dSyv, dIyv, dI2yv, dRyv, dDyv, dDryv,
                    dSa, dIa, dI2a, dRa, dDa, dDra,
                    dPsav, dPiav, dPi2av,dPrav, dPdav, dPdrav,
                    dSav, dIav, dI2av, dRav, dDav, dDrav,
                    dIhy, dPihyv, dIhyv, dIha, dPihav, dIhav,
                    dVy, dVa, 
                    dSoeset)))
      })
  }
  
  # the parameters values:
  parameters_values <- c(rt_ini=rt_ini, rt_rebLength=rt_rebLength, rt_reb=rt_reb, rt_var=rt_var, rt_oly=rt_oly, soe_th=soe_th, soe_eff=soe_eff, soe_time=soe_time, 
                         gamma=gamma, deltaY_temp=deltaY_temp, deltaA_temp=deltaA_temp, deltaL=deltaL, del_var=del_var,
                         hosp_rate=hosp_rate, hosp_plus_len=hosp_plus_len,
                         vaceff_inf=vaceff_inf, vaceff_death_temp=vaceff_death_temp)
  
  # the initial values of variables:
  initial_values <- c(Sy=Sy0, Iy=Iy0, I2y=I2y0, Ry=Ry0, Dy=Dy0, Dry=Dry0,
                      Psyv=Psyv0, Piyv=Piyv0, Pi2yv=Pi2yv0, Pryv=Pryv0, Pdyv=Pdyv0, Pdryv=Pdryv0, 
                      Syv=Syv0, Iyv=Iyv0, I2yv=I2yv0, Ryv=Ryv0, Dyv=Dyv0, Dryv=Dryv0, 
                      Sa=Sa0, Ia=Ia0, I2a=I2a0, Ra=Ra0, Da=Da0, Dra=Dra0, 
                      Psav=Psav0, Piav=Piav0, Pi2av=Pi2av0, Prav=Prav0, Pdav=Pdav0, Pdrav=Pdrav0, 
                      Sav=Sav0, Iav=Iav0, I2av=I2av0, Rav=Rav0, Dav=Dav0, Drav=Drav0, 
                      Ihy=Ihy0, Pihyv=Pihyv0, Ihyv=Ihyv0, Iha=Iha0, Pihav=Pihav0, Ihav=Ihav0,
                      Vy=Vy0, Va=Va0, 
                      Soeset=Soeset0)
  
  # solving
  out <- ode(method="rk4",initial_values, times, sir_equations, parameters_values)
  #out <- lsode(initial_values, times, sir_equations, parameters_values,atol=1e-20,maxsteps=10000)
  
  #method = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
  #           "euler", "rk4", "ode23", "ode45", "radau", 
  #           "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration")

    
  # returning the output:
  simulationresult <- as.data.frame(out)
  return(simulationresult)
}


#Test model
#testresult <- sir_1(rt_ini=1.0, rt_reb=1.3, rt_var=1.5, rt_oly=1.6, gamma = 0.2,
#                    vaceff_inf = 0.8, vaceff_death_temp = 0.9,
#                    deltaY = 0.005, deltaA = 0.05, deltahL = 1/14, dailyDose = 10,
#                    Sy0 = 990, Iy0 = 10, Ry0 = 0, Dy0 = 0, 
#                    Psyv0 = 0, Syv0 = 0, Iyv0 = 0, Ryv0 = 0, Dyv0 = 0, 
#                    Sa0 = 500, Ia0 = 10, Ra0 = 0, Da0 = 0,
#                    Psav0 = 0, Sav0 = 0, Iav0 = 0, Rav0 = 0, Dav0 = 0, 
#                    times = seq(0, 114, by=1))

#ggplot() +
  #  geom_line(data=testresult,aes(x=time,y=Sy, color="Sy")) + 
  #  geom_line(data=testresult,aes(x=time,y=Iy, color="Iy")) + 
  #  geom_line(data=testresult,aes(x=time,y=Ry, color="Ry")) +
  #  geom_line(data=testresult,aes(x=time,y=Dy, color="Dy")) +
  
  #geom_line(data=testresult,aes(x=time,y=Syv, color="Syv")) + 
  #geom_line(data=testresult,aes(x=time,y=Iyv, color="Iyv")) + 
  #geom_line(data=testresult,aes(x=time,y=Ryv, color="Ryv")) +
  #geom_line(data=testresult,aes(x=time,y=Dyv, color="Dyv")) +
  
  #  geom_line(data=testresult,aes(x=time,y=Sa, color="Sa")) + 
  #  geom_line(data=testresult,aes(x=time,y=Ia, color="Ia")) + 
  #  geom_line(data=testresult,aes(x=time,y=Ra, color="Ra")) +
  #  geom_line(data=testresult,aes(x=time,y=Da, color="Da")) +
  
  #  geom_line(data=testresult,aes(x=time,y=Sav, color="Sav")) + 
  #  geom_line(data=testresult,aes(x=time,y=Iav, color="Iav")) + 
  #  geom_line(data=testresult,aes(x=time,y=Rav, color="Rav")) +
  #  geom_line(data=testresult,aes(x=time,y=Dav, color="Dav")) +
  
#xlab('Time') + ylab('Number')




#Simulations

# ver6: D is not death, but severe

#list_rt_reb <- c(1.1,1.2,1.3)
#list_dailyDose <- c(0,25000,50000,100000)


#list_rt_var <- c(1.2)
#list_del_var <- c(1.2)
#list_soe_th <- c(10^10)
#list_soe_time <- c(999)
#list_soe_eff <- c(1)
#list_rt_oly <- c(0)
#list_rt_reb <- c(0.0075)
#list_reb_len <- c(1000)





list_rt_var <- c(1.35)
list_del_var <- c(1.5)
list_soe_th <- c(10^10)

list_soe_time <- c(35,44,49,999)
#1000: 35
#1500: 44
#2000: 49


list_soe_eff <- c(0.8)
list_rt_oly <- c(0)
list_rt_reb <- c(0.006) #delta1.1x0.0095, delta1.2x0.0090,  delta1.3x0.0084, delta1.4x0.0075, delta1.5x0.0063
list_reb_len <- c(31) #time0=June09, time31=july10


i <- 0

for (set_rt_var in list_rt_var) {
  for (set_del_var in list_del_var) {
    for (set_soe_th in list_soe_th) {
      for (set_soe_eff in list_soe_eff) {
        for (set_rt_oly in list_rt_oly) {
          for (set_rt_reb in list_rt_reb) {
            for (set_reb_len in list_reb_len) {
              for (set_soe_time in list_soe_time) {

      i = i+1
      
      vaceff_inf <- 0.8
      vaceff_death_temp <- 0.9
      
      result <- sir_1(rt_ini=1.05, rt_rebLength=set_reb_len, rt_reb=set_rt_reb, rt_var=set_rt_var, rt_oly=set_rt_oly, soe_th=set_soe_th, soe_eff=set_soe_eff, soe_time=set_soe_time, 
                      
                      gamma = 0.2, deltaY_temp = 0.004, deltaA_temp = 0.03, deltaL = 1/14, del_var=set_del_var,
                      
                      hosp_rate = 0.35, hosp_plus_len = 1/5, 
                      
                      vaceff_inf = vaceff_inf, vaceff_death_temp = vaceff_death_temp,

                      
                      Sy0 = 1070*(10^4) - 400*5*0.9 - 50*(10^4) - 37*(10^4)*(1-vaceff_inf) - 37*(10^4)*vaceff_inf - 20 - 500*0.9,
                      Iy0 = 400*5*0.9*0.97, I2y0 = 500*0.9, Ry0 = 0, Dy0 = 20, Dry0 = 0,
                      Psyv0 = 50*(10^4), Piyv0 = 0, Pi2yv0 = 0, Pryv0 = 0, Pdyv0 = 0, Pdryv0 = 0, 
                      Syv0 = 37*(10^4)*(1-vaceff_inf), Iyv0 = 0, I2yv0 = 0, Ryv0 = 0, Dyv0 = 0, Dryv0 =0, 
                      
                      
                      Sa0 = 320*(10^4) - 400*5*0.1 - 76*(10^4) - 7*(10^4)*(1-vaceff_inf) - 7*(10^4)*vaceff_inf - 30 - 500*0.1,
                      Ia0 = 400*5*0.1*0.97, I2a0 = 500*0.1, Ra0 = 0, Da0 = 30, Dra0 = 0, 
                      Psav0 = 76*(10^4), Piav0 = 0, Pi2av0 = 0, Prav0 = 0, Pdav0 = 0, Pdrav0 = 0, 
                      Sav0 = 7*(10^4)*(1-vaceff_inf), Iav0 = 0, I2av0 = 0, Rav0 = 0, Dav0 = 0, Drav0 = 0, 
                      
                      Ihy0 = 400*5*0.9*0.03, Pihyv0 = 0, Ihyv0 = 0, Iha0 = 400*5*0.1*0.03, Pihav0 = 0, Ihav0 = 0,
                      
                      Vy0 =37*(10^4)*vaceff_inf,
                      Va = 7*(10^4)*vaceff_inf, 
                      
                      Soeset0 = 0,
                      
                      times = seq(0, 113, by=1))

      
      
#      sir_1 <- function(rt_ini, rt_reb, rt_var, rt_oly, soe_th, soe_eff,
#                        gamma, deltaY, deltaA, deltaL,
#                        vaceff_inf, vaceff_death_temp, 
#                        Sy0, Iy0, Ry0, Dy0,
#                        Psyv0, Syv0, Iyv0, Ryv0, Dyv0,
#                        Sa0, Ia0, Ra0, Da0,
#                        Psav0, Sav0, Iav0, Rav0, Dav0,
#                        times) {
        

  var_rt_var <- rep(set_rt_var, 114)
  var_del_var <- rep(set_del_var, 114)
  var_soe_th <- rep(set_soe_th, 114)
  var_soe_eff <- rep(set_soe_eff, 114)
  var_rt_oly <- rep(set_rt_oly, 114)
  var_rt_reb <- rep(set_rt_reb, 114)
  var_reb_len <- rep(set_reb_len, 114)
  var_soe_time <- rep(set_soe_time, 114)

  result <- cbind(result,var_rt_var)
  result <- cbind(result,var_del_var)
  result <- cbind(result,var_soe_th)
  result <- cbind(result,var_soe_eff)
  result <- cbind(result,var_rt_oly)
  result <- cbind(result,var_rt_reb)
  result <- cbind(result,var_reb_len)
  result <- cbind(result,var_soe_time)

  
  if (i == 1) {
    comb_result <- result
  } else {
    comb_result <- rbind(comb_result,result)
  }
    
  }}}}}}}}


i_total_per5 <- (comb_result$Iy+comb_result$Piyv+comb_result$Iyv+comb_result$Ia+comb_result$Piav+comb_result$Iav+comb_result$Ihy+comb_result$Pihyv+comb_result$Ihyv+comb_result$Iha+comb_result$Pihav+comb_result$Ihav)/5
severe_total <- comb_result$Dy+comb_result$Pdyv+comb_result$Dyv+comb_result$Da+comb_result$Pdav+comb_result$Dav
delta_prop <- (comb_result$Ihy+comb_result$Pihyv+comb_result$Ihyv+comb_result$Iha+comb_result$Pihav+comb_result$Ihav)/(comb_result$Iy+comb_result$Piyv+comb_result$Iyv+comb_result$Ia+comb_result$Piav+comb_result$Iav+comb_result$Ihy+comb_result$Pihyv+comb_result$Ihyv+comb_result$Iha+comb_result$Pihav+comb_result$Ihav)*100
hosp <- (comb_result$Iy+comb_result$Piyv+comb_result$Ia+comb_result$Piav+comb_result$Ihy+comb_result$Pihyv+comb_result$Iha+comb_result$Pihav+(comb_result$Iyv++comb_result$Iav++comb_result$Ihyv+comb_result$Ihav)/2)*0.35
hosp <- hosp + comb_result$I2y + comb_result$Pi2yv + comb_result$I2yv + comb_result$I2a + comb_result$Pi2av + comb_result$I2av + severe_total

comb_result <- cbind(comb_result, i_total_per5)
comb_result <- cbind(comb_result, severe_total)
comb_result <- cbind(comb_result, delta_prop)
comb_result <- cbind(comb_result, hosp)

