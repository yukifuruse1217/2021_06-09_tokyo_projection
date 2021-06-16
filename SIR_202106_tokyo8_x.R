library(ggplot2) # for graph
require(deSolve) # for the "ode" function
#library(bbmle) # for ML method  install.packages("bbmle", repos="http://R-Forge.R-project.org")



#Equation#######################
# ver6: D is not death, but severe

sir_1 <- function(rt_ini, rt_reb0, rt_reb, rt_var, rt_oly, soe_th, soe_eff,
                  gamma, deltaY_temp, deltaA_temp, deltaL, del_var, var_period,
                  vaceff_inf, vaceff_death_temp, 
                  Sy0, Iy0, Ry0, Dy0,
                  Pyv0, Syv0, Iyv0, Ryv0, Dyv0,
                  Sa0, Ia0, Ra0, Da0,
                  Pav0, Sav0, Iav0, Rav0, Dav0,
                  Dry0, Dryv0, Dra0, Drav0,
                  Vy0, Va0, 
                  Soeset0,
                  times) {
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      jinryu_to_rt <- 0.025
      
      if (time<12) {
        rt <- rt_ini
      } else if (time<(12+28)) {
        rt <- (rt_reb0*jinryu_to_rt)+((rt_reb-rt_reb0)*jinryu_to_rt*(time-12)/28)+rt_ini
      } else {
        rt <- (rt_reb*jinryu_to_rt)+rt_ini
      }
      

      
      
      
      deltaY <- deltaY_temp
      deltaA <- deltaA_temp
      
      #var
      if (time<var_period) {
        rt <- rt * (time/var_period * (rt_var - 1) + 1) * 0.8 + rt * 0.2
      } else {
        rt <- rt * rt_var * 0.8 + rt * 0.2
      }
      
      if (time<var_period) {
        deltaY <- deltaY * (time/var_period * (del_var - 1) + 1) * 0.8 + deltaY * 0.2
        deltaA <- deltaA * (time/var_period * (del_var - 1) + 1) * 0.8 + deltaA * 0.2
      } else {
        deltaY <- deltaY * del_var * 0.8 + deltaY * 0.2
        deltaA <- deltaA * del_var * 0.8 + deltaA * 0.2
      }
      
      
      
      #olympic      
      if (44<time && time<60) {
        rt <- (rt_oly*jinryu_to_rt)+rt
      }
      
      
      #soe_th/eff
      soe <- (Iy+Iyv+Ia+Iav)/5
      
      if (Soeset>0) {
        rt <- rt * soe_eff
      }

      
      
      
      dailyDose <- ifelse(time>13,50000,75000)
      
      adult_vaccination_max <- 0.8 * (320*(10^5))
      young_vaccination_max <- 0.4 * (1070*(10^5))
      
      if (time<12) {      
        adult_vaccination <- ifelse((Pav+Sav+Iav+Rav+Dav+Drav+Va)>=adult_vaccination_max,0,min(dailyDose, (adult_vaccination_max - (Pav+Sav+Iav+Rav+Dav+Drav+Va))))
        young_vaccination <- ifelse((Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy)>=young_vaccination_max,0,min(dailyDose - adult_vaccination, (young_vaccination_max - (Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy))))
      } else if (time<31) {
        adult_vaccination <- ifelse((Pav+Sav+Iav+Rav+Dav+Drav+Va)>=adult_vaccination_max,0,min(50000, (adult_vaccination_max - (Pav+Sav+Iav+Rav+Dav+Drav+Va))))
        young_vaccination <- ifelse((Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy)>=young_vaccination_max,0,min(dailyDose - adult_vaccination, (young_vaccination_max - (Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy))))
      } else if (time<53) {
        adult_vaccination <- ifelse((Pav+Sav+Iav+Rav+Dav+Drav+Va)>=adult_vaccination_max,0,min(25000, (adult_vaccination_max - (Pav+Sav+Iav+Rav+Dav+Drav+Va))))
        young_vaccination <- ifelse((Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy)>=young_vaccination_max,0,min(dailyDose - adult_vaccination, (young_vaccination_max - (Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy))))
      } else {
        adult_vaccination <- ifelse((Pav+Sav+Iav+Rav+Dav+Drav+Va)>=adult_vaccination_max,0,min(dailyDose, (adult_vaccination_max - (Pav+Sav+Iav+Rav+Dav+Drav+Va))))
        young_vaccination <- ifelse((Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy)>=young_vaccination_max,0,min(dailyDose - adult_vaccination, (young_vaccination_max - (Pyv+Syv+Iyv+Ryv+Dyv+Dryv+Vy))))
      }
      
      
      vaccine_delay <- 1/28
      vaceff_death <- (1-vaceff_death_temp) / (1-vaceff_inf)
      
      
      
      
      
      dSy <-  -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Sy -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Sy - young_vaccination
      dPyv <- -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Pyv -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Pyv + young_vaccination - Pyv * vaccine_delay
      dSyv <- -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Syv -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Syv + Pyv * vaccine_delay * (1-vaceff_inf)
      
      dSa <-  -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Sa -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Sa - adult_vaccination
      dPav <- -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Pav -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Pav + adult_vaccination - Pav * vaccine_delay
      dSav <- -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Sav -rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Sav + Pav * vaccine_delay * (1-vaceff_inf)
      
      dIy <-  rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Sy +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Sy +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Pyv +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Pyv - gamma * Iy
      dIyv <- rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Syv +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Syv - gamma * Iyv
      
      dIa <-  rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Sa +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Sa +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Pav +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Pav - gamma * Ia
      dIav <- rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Iy+Iyv) * (1390/1925*1.5) * Sav +rt / (Sy+Iy+Ry+Sa+Ia+Ra+Syv+Iyv+Ryv+Sav+Iav+Rav+Pyv+Pav+Dy+Dyv+Da+Dav+Dry+Dryv+Dra+Drav+Vy+Va) * gamma * (Ia+Iav) * (1390/1925) * Sav - gamma * Iav
      
      dRy <-   gamma * Iy * (1-deltaY)
      dRyv <-  gamma * Iyv * (1-deltaY*vaceff_death)
      dRa <-   gamma * Ia * (1-deltaA)
      dRav <-  gamma * Iav * (1-deltaA*vaceff_death)
      
      dVy <-   Pyv * vaccine_delay * vaceff_inf
      dVa <-   Pav * vaccine_delay * vaceff_inf
      
      dDy <-  gamma * Iy * deltaY - deltaL * Dy
      dDyv <- gamma * Iyv * deltaY*vaceff_death - deltaL * Dyv
      dDa <-  gamma * Ia * deltaA - deltaL * Da
      dDav <- gamma * Iav *deltaA*vaceff_death - deltaL * Dav
      
      dDry <-  deltaL * Dy
      dDryv <- deltaL * Dyv
      dDra <-  deltaL * Da
      dDrav <- deltaL * Dav
      
      if (soe > soe_th) {
        dSoeset <- 1
      } else {
        dSoeset <- 0
      }
      
      return(list(c(dSy, dIy, dRy, dDy, dPyv, dSyv, dIyv, dRyv, dDyv, dSa, dIa, dRa, dDa, dPav, dSav, dIav, dRav, dDav, dDry, dDryv, dDra, dDrav, dVy, dVa, dSoeset)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(rt_ini=rt_ini, rt_reb0=rt_reb0, rt_reb=rt_reb, rt_var=rt_var, rt_oly=rt_oly, soe_th=soe_th, soe_eff=soe_eff,
                         gamma=gamma, deltaY_temp=deltaY_temp, deltaA_temp=deltaA_temp, deltaL=deltaL, del_var=del_var, var_period=var_period, 
                         vaceff_inf=vaceff_inf, vaceff_death_temp=vaceff_death_temp)
  
  # the initial values of variables:
  initial_values <- c(Sy = Sy0, Iy = Iy0, Ry = Ry0, Dy = Dy0,
                      Pyv = Pyv0, Syv = Syv0, Iyv = Iyv0, Ryv = Ryv0, Dyv = Dyv0, 
                      Sa = Sa0, Ia = Ia0, Ra = Ra0, Da = Da0,
                      Pav = Pav0, Sav = Sav0, Iav = Iav0, Rav = Rav0, Dav = Dav0,
                      Dry = Dry0, Dryv = Dryv0, Dra = Dra0, Drav = Drav0,
                      Vy = Vy0, Va = Va0,
                      Soeset = Soeset0)
  
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
#                    Pyv0 = 0, Syv0 = 0, Iyv0 = 0, Ryv0 = 0, Dyv0 = 0, 
#                    Sa0 = 500, Ia0 = 10, Ra0 = 0, Da0 = 0,
#                    Pav0 = 0, Sav0 = 0, Iav0 = 0, Rav0 = 0, Dav0 = 0, 
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


list_rt_var <- c(1,1.2,1.5)
list_del_var <- c(1,1.2,2.0)
list_soe_th <- c(1000,2000,10^8)
list_soe_eff <- c(1,0.7)
list_rt_oly <- c(0,1,4,5,9,10)
list_rt_reb <- c(10,15)
list_val_period <- c(28,56)

i <- 0

for (set_rt_var in list_rt_var) {
  for (set_del_var in list_del_var) {
    for (set_soe_th in list_soe_th) {
      for (set_soe_eff in list_soe_eff) {
        for (set_rt_oly in list_rt_oly) {
          for (set_rt_reb in list_rt_reb) {
            for (set_val_period in list_val_period) {

      i = i+1
      
      vaceff_inf <- 0.8
      vaceff_death_temp <- 0.9
      
      result <- sir_1(rt_ini=0.9, rt_reb0=10, rt_reb=set_rt_reb, rt_var=set_rt_var, rt_oly=set_rt_oly, soe_th=set_soe_th, soe_eff=set_soe_eff,
                      
                      gamma = 0.2, deltaY_temp = 0.004, deltaA_temp = 0.03, deltaL = 1/14, del_var=set_del_var, var_period=set_val_period, 
                      
                      vaceff_inf = vaceff_inf, vaceff_death_temp = vaceff_death_temp,
                      
                      Sy0 = 1070*(10^4) - 390*5*0.9 - 50*(10^4) - 37*(10^4)*(1-vaceff_inf) - 37*(10^4)*vaceff_inf - 20,      
                      Iy0 = 390*5*0.9, Ry0 = 0, Dy0 = 20, 
                      Pyv0 = 50*(10^4), Syv0 = 37*(10^4)*(1-vaceff_inf), Iyv0 = 0, Ryv0 = 0, Dyv0 = 0, Vy0 =37*(10^4)*vaceff_inf, 

                      Sa0 = 320*(10^4) - 390*5*0.1 - 76*(10^4) - 7*(10^4)*(1-vaceff_inf) - 7*(10^4)*vaceff_inf - 30, 
                      Ia0 = 390*5*0.1, Ra0 = 0, Da0 = 30,
                      Pav0 = 76*(10^4), Sav0 = 7*(10^4)*(1-vaceff_inf), Iav0 = 0, Rav0 = 0, Dav0 = 0, Va = 7*(10^4)*vaceff_inf, 
                      
                      Dry0=0, Dryv0=0, Dra0=0, Drav0=0, Soeset0=0,
                      
                      times = seq(0, 113, by=1))

      
      
#      sir_1 <- function(rt_ini, rt_reb, rt_var, rt_oly, soe_th, soe_eff,
#                        gamma, deltaY, deltaA, deltaL,
#                        vaceff_inf, vaceff_death_temp, 
#                        Sy0, Iy0, Ry0, Dy0,
#                        Pyv0, Syv0, Iyv0, Ryv0, Dyv0,
#                        Sa0, Ia0, Ra0, Da0,
#                        Pav0, Sav0, Iav0, Rav0, Dav0,
#                        times) {
        

  var_rt_var <- rep(set_rt_var, 114)
  var_del_var <- rep(set_del_var, 114)
  var_soe_th <- rep(set_soe_th, 114)
  var_soe_eff <- rep(set_soe_eff, 114)
  var_rt_oly <- rep(set_rt_oly, 114)
  var_rt_reb <- rep(set_rt_reb, 114)
  var_val_period <- rep(set_val_period, 114)
  
  result <- cbind(result,var_rt_var)
  result <- cbind(result,var_del_var)
  result <- cbind(result,var_soe_th)
  result <- cbind(result,var_soe_eff)
  result <- cbind(result,var_rt_oly)
  result <- cbind(result,var_rt_reb)
  result <- cbind(result,var_val_period)
  
  
  if (i == 1) {
    comb_result <- result
  } else {
    comb_result <- rbind(comb_result,result)
  }
    
  }}}}}}}









