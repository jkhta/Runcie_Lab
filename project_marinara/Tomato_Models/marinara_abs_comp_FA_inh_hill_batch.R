run = as.numeric(commandArgs(t=T)[1])

library(deSolve)
library(GenSA)

time_step <- 0.01
t_run <- seq(0, 25, by = time_step)
small_step <- c(0, time_step)
time_scale = 1
init <- c(0, 0, 0, 0, 0, 0, 0, 1) #SFT, TMF, SP, SPGB, FA, INF, #leaves, change in #leaves
inf_info <- data.frame()
inf_first <- data.frame()

parameters_list <- list(
  B_1 = 0.01, #production of SFT per leaf
  B_5 = 0.005, #production of FA per leaf/time unit
  B_3_1 = 0.05, #production of SP before branching
  B_3_2 = 0.1, #production of SP after branching
  B_5_2 = 0.01, #production of TMF without FA inhibition
  B_5_4 = 0.01, #activation of SPGB by FA
  B_14_5 = 0.01, #activation of FA by SFT:SPGB complex
  B_6_5 = 0.01, #activation of FA by INF
  B_14_6 = 0.3, #activation of INF by SFT:SPGB complex
  B_5_6 = 0.05, #activation of INF by FA
  B_1_5 = 0.01, #activation of FA by SFT in the inflorescence meristem
  K_34 = 0.01, 
  K_14 = 0.01,
  K_5_2 = 0.01,
  K_5_4 = 0.01,
  K_14_5 = 0.01,
  K_6_5 = 0.01,
  K_34_5 = 0.01,
  K_5_6 = 0.01,
  K_2_6 = 0.01,
  K_14_6 = 0.01,
  K_34_6 = 0.01,
  K_1_5 = 0.01,
  d1 = 0.3,
  d2 = 0.3,
  d3 = 0.3,
  d4 = 0.3,
  d5 = 0.3,
  d6 = 0.3,
  split = 0.3,
  sink_decrease = 0.01
)

set_parameters <- list(
  v_35S = time_scale*c(rep(0, 6)),
  mutants = rep(0,6)
  #max_sink = 0.250 #the maximum amount of SFT that can be pooled in the inflorescence meristem/could also be the maximum amount of flux into the inflorescence meristem
  #current model is sensitive to sink decrease
)
set_parameters$init <- init

parameter_set <- c(parameters_list, set_parameters)
parameter_set$init <- init

#creates an inflorescence (row) in data frame whenever SFT exceeds the arbitrary threshold

#changing inflorescence threshold (SFT concentration needed) affects flowering time of initial segment and also sympodial segments, which it shouldn't
create_inf <- function(inf,a,split_value) {
  if(inf > 0.1) {
    y <- c(0, 0, 0, a, split_value, 0) #SFT flux, SFT, FA, time initiated, sink power, flowering status
    inf_info <<- rbind(inf_info, y)
    return(inf_info)
  }else{
    return(inf_info)
  }
}

tomato_inf <- function(t,X,parms = NULL,...){
  with(as.list(parms),{
    FA = B_1_5 * X[2]/(K_1_5 + X[2])
    SFT = X[1]
    derivatives_im = rep(0, length(X))
    derivatives_im[1:3] = c(0, SFT, FA)
    derivatives_im[1:3] = derivatives_im[1:3] - 0.1 * X[1:3]
    return(list(
      Derivatives = derivatives_im,
      globals = c(
        SFT,
        FA
      )
    ))
  })
}

#currently the thresholds are based on fluxes rather than protein concentrations (gotta see how to change this)
#right now the veg model should only generate SFT, and transfer that into the mersitem
tomato_veg_valentim <- function(t, X, parms=NULL,...) {
  
  with(as.list(parms),{
    
    #Competitive inhibition/binding between FT/TFL1 and FD
    x_14 = K_34 * X[1] * X[4]/(K_14 * K_34 + K_14 * X[3] + K_34 * X[1])
    x_34 = K_14 * X[3] * X[4]/(K_14 * K_34 + K_14 * X[3] + K_34 * X[1])
    
    #Competitive binding between SFT:SPGB and SP:SPGB for FA and INF binding spots
    p_14_5 = K_34_5^2 * x_14^2/(K_14_5^2 * K_34_5^2 + K_14_5^2 * x_34^2 + K_34_5^2 * x_14^2)
    p_14_6 = K_34_6^2 * x_14^2/(K_14_6^2 * K_34_6^2 + K_14_6^2 * x_34^2 + K_34_6^2 * x_14^2)
    
    #Amount of florigen produced per leaf
    SFT = B_1 * X[7]
    
    #Activation of TMF by FA
    TMF = B_5_2 * K_5_2^2/(K_5_2^2 + X[5]^2)
    
    #When the first inflorescence meristem is created, SP is then produced
    if(nrow(inf_info) > 0){
      SP = B_3_2
    }else{
      SP = B_3_1
    }
    
    #Activation of SPGB by FA
    SPGB = B_5_4 * X[5]^2/(K_5_4^2 + X[5]^2)
    
    #Amount of FA produced over plant age + [activation by SFT:SPGB + activation by INF] * inhibition by SP:SPGB
    FA = B_5 * X[7] + 
         B_14_5 * p_14_5 + 
         B_6_5 * X[6]^2/(K_6_5^2 + X[6]^2)
    
    # * (K_14_5 * x_34/(K_14_5 * K_34_5 + K_14_5 * x_34 + K_34_5 * x_14))
    
    #[Activation of INF by SFT:SPGB + activation by FA]
    INF = B_14_6 * p_14_6 + 
          (B_5_6 * X[5]^2/(K_5_6^2 + X[5]^2)) * K_2_6^2/(K_2_6^2 + X[2]^2)
    
    # * (K_14_5 * x_34/(K_14_5 * K_34_5 + K_14_5 * x_34 + K_34_5 * x_14))
    
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
    
    # This is saying each inflorescence meristem takes a proportion of florigen
    if(nrow(inf_info) > 0){
      for(i in 1:nrow(inf_info)){
        #absolute value of florigen/protein taken by each inflorescence meristem
        SFT_im = inf_info[i,5]
        FA_im = inf_info[i,5]
        INF_im = inf_info[i,5]
        
        #proportion of florigen/protein taken by each inflorescence meristem
        # SFT_im = SFT * inf_info[i,5]
        # FA_im = FA * inf_info[i,5]
        # INF_im = INF * inf_info[i,5]
        SFT <- SFT - SFT_im
        FA <- FA - FA_im
        INF <- INF - INF_im
        inf_info[i,1] <<- SFT_im
      }
    }
    
    SFT_veg = SFT
    FA_veg = FA
    INF_veg = INF
    
    #order of this is really important, it should come after reduction of vegetative partitioning
    create_inf(INF_veg,as.numeric(t),split)
    
    #multiples the SFT flux in each inflorescence mersitem by the time step because running this way does not automatically scale the ODE solving
    # if(nrow(inf_info) > 0){
    #   inf_info[,1] <<- inf_info[,1] * time_step
    # }
    
    #calculates SFT pool based on the SFT flux into each inflorescence meristem that will lead to LFY production
    
    if(nrow(inf_info) > 0){
      for(i in 1:nrow(inf_info)){
        # s2 <- inflorescence_check(inf_info[i,])
        # inf_info[i,] <<- s2[2,][2:7]
        #sink power should follow some curve (beta function? normal distribution?)
        inf_info[i,5] <<- inf_info[i,5] - sink_decrease * time_step
      }
    }
    
    # print(t)
    # print(inf_info)
    
    #trying to get inflorescence maturation information to plot on the same graph as the vegetative growth graph
    # a <<- inf_info[1,]
    # inf_first <<- rbind(inf_first,a)
    
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
    derivatives[1:6] = c(SFT_veg, TMF, SP, SPGB, FA_veg, INF_veg)
    derivatives[1:6] = derivatives[1:6] * (1 - mutants)
    derivatives[7] = X[8]
    delta = c(d1,d2,d3,d4,d5,d6,0,0)
    derivatives = derivatives - delta*X
    
    if(nrow(inf_info) == 4){
      break
    }
    
    return(list(
      Derivatives = derivatives,
      globals = c(
        x_14 = x_14,
        x_34 = x_34,
        p_14_5 = p_14_5,
        p_14_6 = p_14_6
      )
    ))
  })
}

fit_model = function(parms){
  s1 <- tryCatch({
    ode(y = c(parms$init),
            times = t_run,
            func = tomato_veg_valentim,
            parms=parms,
            method='euler')
  },
  error = function(x) return(NA)
  )
            # rootfun = root_fun,
            # events = list(func = events_fun,root=T))
  return(s1)
}

fit_model_old = function(parms){
  s1 <- ode(y = c(parms$init),
            times = t_run,
            func = tomato_veg_valentim,
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
            parms = parameters_list,
            method = 'euler')
  return(s2)
}

genotype_parms = function(genotype,parms){
  new_parms <<- parms #parms
  new_parms$init <- init
  
  if(genotype == 'sft'){
    new_parms$mutants[1] = 0.8
  }
  if(genotype == 'tmf'){
    new_parms$mutants[2] = 0.8 
  }
  if(genotype == 'sp'){
    new_parms$mutants[3] = 0.8
  }
  if(genotype == 'fa'){ 
    new_parms$mutants[5] = 0.8 
  }
  return(new_parms)
}

obj_fun <- function(params) {
  j <- as.list(params)
  names(j) <- names(parameters_list)
  counter <<- counter + 1
  
  data_model <- read.delim('Tomato_Mutant_Data.csv',sep=',')
  data_model$pred_First = NA
  data_model$pred_Second = NA
  data_model$pred_Third = NA
  data_model$pred_Fourth = NA
  
  for(gen in data_model$Genotype){
    i = data_model$Genotype == gen
    new_parms <<- genotype_parms(gen,c(j,set_parameters))
    inf_info <<- data.frame()
    inf_first <<- vector()
    s1 <- fit_model(new_parms)
    inf_first <- inf_info[1:4,4]
    if(length(inf_first) == 0){
      inf_first <- rep(NA,4)
    }
    data_model[i,6:9] = inf_first
  }
  
  data_model[is.na(data_model)] <- 50
  score <- sum((data_model[2:5]-data_model[6:9])^2)
  out_line <- as.vector(unlist(c(j, score, counter)))
  data_out <- as.data.frame(t(unname(out_line)))
  write.table(data_out, file = sprintf("full_run%f.csv", run), append = TRUE, sep = ",", col.names = FALSE)
  print(out_line)
  return(score)
}

low <- c(0.01,0.005,0.05,0.1,rep(0.05,7),rep(0.01,12),rep(0.3,6),0.01,0)
high <- c(1,0.3,0.1,0.5,rep(1,7),rep(10,12),rep(0.7,6),0.9,1)

op_parms <- c(sample(seq(0.01, 1, by = 0.001), 1),
              sample(seq(0.005, 0.3, by = 0.001), 1),
              sample(seq(0.05, 0.1, by = 0.01), 1),
              sample(seq(0.1, 0.5, by = 0.01), 1),
              sample(seq(0.05, 1, by = 0.01), 7),
              sample(seq(0.01, 10, by = 0.01), 12),
              sample(seq(0.3, 0.7, by = 0.01), 6),
              sample(seq(0.01, 0.9, by = 0.01), 1),
              sample(seq(0, 1, by = 0.01), 1))

counter <- 0

out <- GenSA(op_parms, obj_fun, lower = low, upper = high)
dat <- data.frame(out)
data_final <- cbind(names(parameters_list), op_parms, dat)
write.csv(data_final, file = sprintf("run_%f.csv", run))