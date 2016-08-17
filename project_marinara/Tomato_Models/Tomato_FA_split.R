library(deSolve)
rm(list = ls())

time_step <- 0.01
t_run <- seq(0, 25, by = time_step)
small_step <- c(0, time_step)
time_scale = 1
init <- c(0, 0, 0, 0, 0, 0, 1) #SFT, SP, SPGB, INF, FA, #leaves, change in #leaves
inf_info <- data.frame()
inf_first <- data.frame()

#         Protein Stuff
#         X[1] = SFT
#         X[2] = TMF
#         X[3] = INF
#         X[4] = FA
#         X[5] = Change in leaves

#NOT SURE WHY MODEL IS INSENSITIVE TO SP
parms_list <- list(
  sft_split = 0.45, # proportion of SFT that goes to inflorescence meristem
  sft_prod = 0.02, #affects the transition into inflorescence meristem (first segment is inf_threshold/sft_prod)
  fa_prod = 0.01,
  sp_prod = 0.2,
  spgb_prod = 0.1,
  K_1_2 = 0.001,
  K_13 = 9,
  K_23 = 3.3,
  K_13_4 = 0.04, #transition into IM is particularly sensitive to this for values between 0-0.5
  K_23_4 = 3,
  K_5_4 = 0.28,
  h_1_2 = 4,
  h_1_3 = 4,
  h_13_4 = 4,
  h_23_4 = 4,
  h_5_4 = 4,
  delta = c(time_scale*2*c(0.05,0.05,0.05,0.05,0.05),0,0),
  # delta2 = c(time_scale*2*c())
  v_35S = time_scale*c(rep(0, 5)),
  mutants = rep(0,5),
  max_sink = 0.250, #the maximum amount of SFT that can be pooled in the inflorescence meristem/could also be the maximum amount of flux into the inflorescence meristem
  sink_decrease = 0.05 #current model is sensitive to sink decrease
)

parms_list$init <- init

#creates an inflorescence (row) in data frame whenever SFT exceeds the arbitrary threshold

#changing inflorescence threshold (SFT concentration needed) affects flowering time of initial segment and also sympodial segments, which it shouldn't
create_inf <- function(inf,fa,a,split_value) {
  if(inf > 0.075|fa > 0.13) {
    y <- c(0, 0, 0, a, split_value, 0) #SFT flux, SFT, FA, time initiated, sink power, flowering status
    inf_info <<- rbind(inf_info, y)
    return(inf)
  }else{
    return(inf)
  }
}

tomato_inf <- function(t,X,parms = NULL,...){
  with(as.list(parms),{
    p_SFT_FA = X[2]^h_1_2/(K_1_2^h_1_2 + X[2] * h_1_2)
    SFT = X[1]
    FA = (1 - p_SFT_FA) * 0.05 + p_SFT_FA * 0.5
    derivatives_im = rep(0, length(X))
    derivatives_im[1:3] = c(0, SFT, FA)
    derivatives_im[1:3] = derivatives_im[1:3] - 0.2 * X[1:3]
    return(list(
      Derivatives = derivatives_im,
      globals = c(
        p_SFT_FA,
        SFT,
        FA
      )
    ))
  })
}

#currently the thresholds are based on fluxes rather than protein concentrations (gotta see how to change this)
#right now the veg model should only generate SFT, and transfer that into the mersitem
tomato_veg <- function(t, X, parms=NULL,...) {
  
  with(as.list(parms),{
    
    #SFT:SPGB
    x_13 = K_23*X[1]*X[3]/(K_13*K_23 + K_13*X[2]+K_23*X[1])
    #SP:SPGB
    x_23 = K_13*X[2]*X[3]/(K_13*K_23 + K_13*X[2]+K_23*X[1])
    
    #Probabilities of binding to INF protein
    #SFT:SPGB -> INF
    p_13_4 = K_23_4^h_23_4 * x_13^h_13_4/(K_13_4^h_13_4 * K_23_4^h_23_4 + K_23_4^h_23_4*x_13^h_13_4 + K_13_4^h_13_4*x_23^h_23_4)
    
    #SP:SPGB -> INF
    p_23_4 = K_13_4^h_13_4 * x_23^h_23_4/(K_13_4^h_13_4 * K_23_4^h_23_4 + K_23_4^h_23_4*x_13^h_13_4 + K_13_4^h_13_4*x_23^h_23_4)
    
    #FA -> INF
    p_5_4 = X[5]^h_5_4/(K_5_4^h_5_4 + X[5]^h_5_4) #adding this makes the flux thresholds more distinguishable
    
    SFT_tot = sft_prod * X[6] 
    FA_tot = fa_prod * X[6]
    INF_tot = (1 - p_13_4 - p_23_4)*(1 - p_5_4) * 0.01 + 
              (p_13_4*(1 - p_5_4) + (1 - p_13_4 - p_23_4)*p_5_4)*0.05 +
              (p_13_4*p_5_4) * 0.1 
    
    SPGB = spgb_prod

    #check to see if there are any inflorescence meristems; if there are then calculates how much each inflorescence meristem takes away the SFT pool
    ###additional features to add
    ###1. tell when the IM transitions to FM based on LFY concentration
    ###2. add max amount of SFT flux/pool that can go to each inflorescence meristem
    ###3. sink power based on age/maturation of inflorescene meristem !PARTIALLY DONE!
    ###4. output the time of each FM transition
    
    ###issues right now
    ###ODE solver solves for the ODEs twice, so that the inflorescence information does not reset between each "solve" !SOLVED!
    ###need to find out how fast inflorescence meristems mature/produce flowers compared to maturation of sympodial meristems
    
    #partially solved the double solve problem by setting inf_info to be an empty data frame when t == 0
    
    if(t == 0){
      inf_info <<- data.frame()
    }
    
    if(nrow(inf_info) > 0){
      for(i in 1:nrow(inf_info)){
        SFT_im = SFT_tot * inf_info[i,5]
        FA_im = FA_tot * inf_info[i,5]
        INF_im = INF_tot * inf_info[i,5]
        SFT_tot = SFT_tot - SFT_im
        FA_tot = FA_tot - FA_im
        INF_tot = INF_tot - INF_im
        inf_info[i,1] <<- SFT_im
      }
    }
    
    #multiples the SFT flux in each inflorescence mersitem by the time step because running this way does not automatically scale the ODE solving
    if(nrow(inf_info) > 0){
      inf_info[,1] <<- inf_info[,1] * time_step
    }
    
    #calculates SFT pool based on the SFT flux into each inflorescence meristem that will lead to LFY production
    
    if(nrow(inf_info) > 0){
      for(i in 1:nrow(inf_info)){
        s2 <- inflorescence_check(inf_info[i,])
        inf_info[i,] <<- s2[2,][2:7]
        #sink power should follow some curve (beta function? normal distribution?)
        inf_info[i,5] <<- inf_info[i,5] - sink_decrease * time_step
      }
    }
    
    print(t)
    print(inf_info)
    
    if(nrow(inf_info) > 0){
      SP = sp_prod
    }else{
      SP = 0
    }
    
    #trying to get inflorescence maturation information to plot on the same graph as the vegetative growth graph
    # a <<- inf_info[1,]
    # inf_first <<- rbind(inf_first,a)
  
    SFT_veg <- SFT_tot
    FA_veg <- FA_tot
    INF_veg <- INF_tot
    
    create_inf(INF_veg,FA_veg,as.numeric(t),sft_split) * inf_loss
    #check and calculate SFT taken by sink, then subtract from total SFT
    # for(i in 1:nrow(inf_info)){
    #   p = SFT_tot * inf_info[i,4]
    #   if(p > max_sink) {
    #     p = max_sink
    #   }
    #   else{
    #     p = SFT_tot * inf_info[i,4]
    #   }
    #   inf_info[i,1] 
    #   SFT_tot = SFT_tot - p
    # }
    
    #subtract from SFT_total by proportion calculated with sink power
    # for(i in 1:nrow(inf_info)){
    #   p = SFT_tot * inf_info[i,4]
    #   SFT_tot = SFT_tot - p
    # }
    
    
    #inflorescence creator
  
    #inflorescence created based on SFT flux
    #create_inf(X[1], as.numeric(t)) #inflorescence created based on SFT
    
    derivatives = rep(0, length(X))
    derivatives[1:5] = c(SFT_veg, SP, SPGB, INF_veg, FA_veg)
    derivatives[1:5] = derivatives[1:5] * (1 - mutants)
    derivatives[6] = X[7]
    derivatives = derivatives - delta*X
    return(list(
      Derivatives = derivatives,
      globals = c(
        SFT_tot = SFT_tot,
        FA_tot = FA_tot,
        INF_tot = INF_tot,
        SFT_veg = SFT_veg,
        FA_veg = FA_veg,
        INF_veg = INF_veg,
        x_13 = x_13,
        x_23 = x_23,
        p_13_4 = p_13_4,
        p_23_4 = p_23_4,
        p_5_4 = p_5_4
      )
    ))
  })
}

fit_model = function(parms){
  s1 <- ode(y = c(parms$init),
            times = t_run,
            func = tomato_veg,
            parms=parms,
            method='euler')
  # rootfun = root_fun,
  # events = list(func = eventsfun,root=T,terminalroot=terminalroot))
  return(s1)
}

inflorescence_check <- function(values){
  # j <- seq(as.numeric(values[3]),max(t_run), by = 0.1)
  values <- as.numeric(values)
  s2 <- ode(y = values,
            times = small_step,
            func = tomato_inf,
            parms = parms_list,
            method = 'euler')
  return(s2)
}

genotype_parms = function(genotype,parms){
  new_parms <<- parms #parms
  new_parms$init <- init
  #Col 					#check
  if(genotype == '35S:FT'){ #check
    new_parms$v_35S[1] = exp_35S #1.3-1.8
  }
  if(genotype == '35S:LFY'){
    new_parms$v_35S[4] = exp_35S #nothing gets is fast enough (minimum = 5 Ros leaves)
  }
  if(genotype == '35S:TFL1'){ #not with this parameter
    new_parms$v_35S[2] = exp_35S #1
  }
  return(new_parms)
}

inf_info <- data.frame()
inf_first <- data.frame()
s1 <- fit_model(parms_list)


cols = c('red','blue','black','green','gray')
time_scale <- 1
x=t_run
x=x*time_scale
plot(NA,NA,xlim = range(t_run),ylim = c(0,1))#range(s1[,-1]))
lines(s1[,1]*time_scale,s1[,2],col=cols[1])
lines(s1[,1]*time_scale,s1[,3],col=cols[2])
lines(s1[,1]*time_scale,s1[,4],col=cols[3])
lines(s1[,1]*time_scale,s1[,5],col=cols[4])
lines(s1[,1]*time_scale,s1[,6],col=cols[5])
lines(s1[,1]*time_scale,s1[,11],col=cols[5])
lines(s1[,1]*time_scale,s1[,1],col=cols[1], lty = 2)
lines(s1[,1]*time_scale,s1[,13],col=cols[2], lty = 2)
lines(s1[,1]*time_scale,s1[,14],col=cols[3], lty = 2)


for(i in 2:5){
  lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
}

legend('topright',legend=c('SFT','TMF','IM','FA'),col=cols,lty=1, cex = 0.75)


#can optimize model based on time of initiation; they should be 10, 13, 16...etc.