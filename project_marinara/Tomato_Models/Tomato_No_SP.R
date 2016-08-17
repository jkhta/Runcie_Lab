library(deSolve)
rm(list = ls())

time_step <- 0.01
t_run <- seq(0, 30, by = time_step)
small_step <- c(0, 0.01)
time_scale = 1
init <- c(0, 0.1, 0, 0.1, 0, 1)
inf_info <- data.frame()
inf_first <- data.frame()

#         Protein Stuff
#         X[1] = SFT
#         X[2] = TMF
#         X[3] = INF
#         X[4] = FA
#         X[5] = Change in leaves

parms_list <- list(
  sft_split = 0.50, # proportion of SFT that goes to inflorescence meristem
  sft_prod = 0.1, #affects the transition into inflorescence meristem (first segment is inf_threshold/sft_prod)
  fa_prod = 0.075,
  K_1_2 = 0.001,
  K_1_3 = 0.001,
  K_2_3 = 9,
  K_3_4 = 0.001,
  h_1_2 = 4,
  h_1_3 = 4,
  h_2_3 = 4,
  h_3_4 = 4,
  delta = c(time_scale*2*c(0.05,0.05,0.05,0.05),0,0),
  # delta2 = c(time_scale*2*c())
  v_35S = time_scale*c(rep(0, 4)),
  mutants = rep(0,4),
  max_sink = 0.250
)

parms_list$init <- init

#creates an inflorescence (row) in data frame whenever SFT exceeds the arbitrary threshold

#changing inflorescence threshold (SFT concentration needed) affects flowering time of initial segment and also sympodial segments, which it shouldn't
create_inf <- function(sft,fa,a,split_value) {
  if(sft > 0.8|fa > 1) {
    y <- c(0, 0, 0, a, split_value, 0) #SFT flux, SFT, FA, time initiated, sink power, flowering status
    inf_info <<- rbind(inf_info, y)
  }
}

tomato_inf <- function(t,X,parms = NULL,...){
  with(as.list(parms),{
    p_SFT_FA = X[2]^h_1_2/(K_1_2^h_1_2 + X[2] * h_1_2)
    SFT = X[1]
    FA = (1 - p_SFT_FA) * 0.01 + p_SFT_FA * 0.5
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

#right now the veg model should only generate SFT, and transfer that into the mersitem
tomato_veg <- function(t, X, parms=NULL,...) {
  
  with(as.list(parms),{
    
    #state variables
    
    #how much each consecutive inflorescence mersitem is diverting away from the total SFT
    j = nrow(inf_info)
    k = (1 - sft_split)^j #proportion that goes to the vegetative apex
    
    # if(nrow(inf_info) > 1){
    #   l = 1:j-1
    # }else{
    #   l = 0
    # }
    # 
    # m = rep(0.4,j)
    # n = sft_split*m^l
    
    #if inflorescence does not have sink power, assign sink power from vector of calculated sink powers
    # for (i in 1:nrow(inf_info)){
    #   if(is.numeric(inf_info[i,4]) == FALSE){
    #     inf_info[i,4] = n[i]
    #   }
    # }
    
    #subtract from SFT_total by proportion calculated with sink power, up to a maximum of p
    
    SFT_tot = sft_prod * X[5] #total SFT production
    # FA_tot = fa_prod * X[5]
    
    print(t)
    
    #check to see if there are any inflorescence meristems; if there are then calculates how much each inflorescence meristem takes away the SFT pool
    ###additional features to add
    ###1. tell when the IM transitions to FM based on LFY concentration
    ###2. add max amount of SFT flux/pool that can go to each inflorescence meristem
    ###3. sink power based on age/maturation of inflorescene meristem
    ###4. output the time of each FM transition
    
    ###issues right now
    ###ODE solver solves for the ODEs twice, so that the inflorescence information does not reset between each "solve"
    ###need to find out how fast inflorescence meristems mature/produce flowers compared to maturation of sympodial meristems
    
    #partially solved the double solve problem by setting inf_info to be an empty data frame when t == 0
    if(t == 0){
      inf_info <<- data.frame()
    }
    
    if(nrow(inf_info) > 0){
      for(i in 1:nrow(inf_info)){
        SFT_im = SFT_tot * inf_info[i,5]
        SFT_tot <- SFT_tot - SFT_im
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
        m <<- t - inf_info[i,4]
        inf_info[i,5] <<- inf_info[i,5] - 0.02 * time_step
      }
    }
    
    #trying to get inflorescence maturation information to plot on the same graph as the vegetative growth graph
    # a <<- inf_info[1,]
    # inf_first <<- rbind(inf_first,a)
    
    print(inf_info)
    SFT_veg <- SFT_tot
    FA_veg <- fa_prod * X[5]
    
    TMF = 0.1 #TMF production at vegetative apex
    
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
    
    IM = 1
    
    #inflorescence creator
    
    create_inf(SFT_veg,FA_veg,as.numeric(t),sft_split) #inflorescence created based on SFT flux
    #create_inf(X[1], as.numeric(t)) #inflorescence created based on SFT
    
    derivatives = rep(0, length(X))
    derivatives[1:4] = c(SFT_veg, TMF, IM, FA_veg)
    derivatives[1:4] = derivatives[1:4] * (1 - mutants)
    derivatives[5] = X[6]
    derivatives = derivatives - delta*X
    return(list(
      Derivatives = derivatives,
      globals = c(
        SFT_tot,
        SFT_veg
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

inf_info <- data.frame()
inf_first <- data.frame()
s1 <- fit_model(parms_list)

cols = c('red','blue','black','green','gray')
time_scale <- 1
x=t_run
x=x*time_scale
plot(NA,NA,xlim = range(t_run),ylim = c(0,20))#range(s1[,-1]))
lines(s1[,1]*time_scale,s1[,2],col=cols[1])
for(i in 2:5){
  lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
}

legend('topright',legend=c('SFT','TMF','IM','FA'),col=cols,lty=1, cex = 0.75)


#can optimize model based on time of initiation; they should be 10, 13, 16...etc.