run = as.numeric(commandArgs(t=T)[1])

library(deSolve)
library(GenSA)
library(Rcpp)

sourceCpp("Silly.cpp")

time_scale = 1  # shifts the timescale by scaling both delta and vs and eta_leaf
time_scale2 = 3
init = c(0,0.6,.1,0.1,0,0,1) # Starts with some FD and LFY, and with a leaf production rate = 1 per unit t.
t = seq(0, 50, by=0.1)

v_lfy = 0.05;v_ap1 = 0.05
op_parms <- c(sample(seq(0.01, 10, by = 0.01), 9), sample(seq(1, 4, by = 0.01), 8))
init_parms <- op_parms
low <- c(rep(0.025, 9), rep(1, 8))
upp <- c(rep(10, 9), rep(4, 8))

exp_35S = 1

baa <- list(
  
  delta    = c(time_scale2*2*c(0.05,0.05,0.05,v_lfy,v_ap1),0,0),
  v_35S    = time_scale2*c(0,rep(0,4)),		#check this. 
  v1       = time_scale2*1*c(rep(0.01,4),0),
  v2       = time_scale2*1*c(0.05,0.05,0.05,v_lfy,v_ap1),
  v3       = time_scale2*2*c(0.05,0.05,0.05,v_lfy,v_ap1),
  eta_leaf = time_scale2*0.01,
  
  T_f = 0.2,
  
  mutants = rep(0,5),
  
  repression = 1
)

baa$init <- init

baa2 <- list(
  
  delta    = c(time_scale*2*c(0.05,0.05,0.05,v_lfy,v_ap1),0,0),
  v_35S    = time_scale*c(0,rep(0,4)),		#check this. 
  v1       = time_scale*1*c(rep(0.01,4),0),
  v2       = time_scale*1*c(0.05,0.05,0.05,v_lfy,v_ap1),
  v3       = time_scale*2*c(0.05,0.05,0.05,v_lfy,v_ap1),
  eta_leaf = time_scale*0.01,
  
  T_f = 0.2,
  
  mutants = rep(0,5),
  
  repression = 1
)

baa2$init <- init
parms_ori1 = list(
  K_13 = 0.39381,
  K_23 = 3.2556,
  K_4_3 = 0.28203,
  K_23_4 = 9.3767,
  K_13_4 = 0.040555,
  K_23_5 = 0.033666,
  K_13_5 = 0.029081,
  K_4_5 = 0.13032,
  K_5_4 = 0.28606,
  h_4_3 = 4.00,
  h_23_4 = 3.8497,
  h_13_4 = 4.00,
  h_23_5 = 4.00,
  h_13_5 = 1.8217,
  h_4_5 = 3.9369,
  h_5_4 = 3.6732,
  h_5_2 = 1.0239
)
parms_ori = list(
  K_13 = 0.39381,
  K_23 = 3.2556,
  K_4_3 = 0.28203,
  K_23_4 = 9.3767,
  K_13_4 = 0.040555,
  K_23_5 = 0.033666,
  K_13_5 = 0.029081,
  K_4_5 = 0.13032,
  K_5_4 = 0.28606,
  h_4_3 = 4.00,
  h_23_4 = 3.8497,
  h_13_4 = 4.00,
  h_23_5 = 4.00,
  h_13_5 = 1.8217,
  h_4_5 = 3.9369,
  h_5_4 = 3.6732,
  h_5_2 = 1.0239,
  
  delta    = c(time_scale*2*c(0.05,0.05,0.05,v_lfy,v_ap1),0,0),
  v_35S    = time_scale*c(0,rep(0,4)),		#check this. 
  v1       = time_scale*1*c(rep(0.01,4),0),
  v2       = time_scale*1*c(0.05,0.05,0.05,v_lfy,v_ap1),
  v3       = time_scale*2*c(0.05,0.05,0.05,v_lfy,v_ap1),
  eta_leaf = time_scale*0.01, #eta_leaf = time_scale*0.01
  
  T_f = 0.2,
  
  mutants = rep(0,5),
  
  repression = 1
)

parms_ori$init <- init

op_parms5 <- list(
  K_13 = 7.018845587,
  K_23 = 1.119881897,
  K_4_3 = 0.228599043,
  K_23_4 = 9.389809234,
  K_13_4 = 0.120358369,
  K_23_5 = 0.400884606,
  K_13_5 = 0.005308391,
  K_4_5 = 2.487942432,
  K_5_4 = 0.157866162,
  h_4_3 = 2.103723153,
  h_23_4 = 3.879031031,
  h_13_4 = 3.692675094,
  h_23_5 = 3.663484422,
  h_13_5 = 2.497241006,
  h_4_5 = 0.824648327,
  h_5_4 = 0.961109439,
  h_5_2 = 2.572238745
)

root_fun = function(t,y,parms,...){
  # This tells the ODE solver to trigger an event. It returns a vector. Events are triggered each time an element = 0.
  return(c(y[5]-0.2,y[5]-0.3)) 
}

# Can use eventsdat to specific changes in pararmeters at specific points in time
# ex. change in FT at certain time
eventsdat = data.frame(var=c(1,2),time=10,value=c(1,0),method='rep')

# eventsfun is called whenever a root is reached 
eventsfun = function(t,y,parms,...){
  if(y[5] > 0.3) y[7] = 0
  y
}

terminalroot = 3 # The 2nd root causes the simulation to stop

fit_model_new = function(parms){
  s1 <- ode(y = c(parms$init),
            times = t,
            func = c_jaeger_model_V3,
            parms=parms,
            method='lsoda',
            rootfun = root_fun,
            events = list(func = eventsfun,root=T,terminalroot=terminalroot))
  
  return(s1)
}

predict_leaves = function(parms){
  s1 = fit_model_new(parms)
  return(attributes(s1)[['troot']])
}

genotype_parms = function(genotype,parms){
  new_parms <- parms #parms
  #Col 					#check
  if(genotype == '35S:FT'){ #check
    new_parms$v_35S[1] = exp_35S #1.3-1.8
    new_parms$init <- c(10, 0.1, 0.1, 0.1, 0, 0, 0)
  }
  if(genotype == '35S:LFY'){
    new_parms$v_35S[4] = exp_35S #nothing gets is fast enough (minimum = 5 Ros leaves)
    new_parms$init <- c(0, 0.1, 0.1, 10.1, 0, 0, 1)
  }
  if(genotype == '35S:TFL1'){ #not with this parameter
    new_parms$v_35S[2] = exp_35S #1
    new_parms$init <- c(0, 10.1, 0.1, 0.1, 0, 0, 1)
  }
  if(genotype == 'lfy-12'){   #check
    new_parms$mutants[4] = 1
    new_parms$init <- c(0, 0.1, 0.1, 0, 0, 0, 1)
  }
  if(genotype == 'ft-10'){	#check
    new_parms$mutants[1] = 1
    new_parms$init <- c(0, 0.1, 0.1, 0.1, 0, 0, 1)
  }
  if(genotype == 'tfl-1'){   #check
    new_parms$mutants[2] = 1
    new_parms$init <- c(0, 0, 0.1, 0.1, 0, 0, 1)
  }
  if(genotype == 'fd-2'){    #check
    new_parms$mutants[3] = 0.75
    new_parms$init <- c(0, 0.1, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == 'fdp-1'){   #check
    new_parms$mutants[3] = 0.2
    new_parms$init <- c(0, 0.1, 0.08, 0.1, 0, 0, 1)
  }
  if(genotype == 'fd-2 fdp-1'){   #check
    new_parms$mutants[3] = 0.95
    new_parms$init <- c(0, 0.1, 0.005, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:TFL1 fd-2'){  #check, exp_35S = 1
    new_parms$v_35S[2] = exp_35S
    new_parms$mutants[3] = 0.75
    new_parms$init <- c(0, 10.1, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == 'tfl1-1 fd-2'){   #check
    new_parms$mutants[2] = 1
    new_parms$mutants[3] = 0.75
    new_parms$init <- c(0, 0, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:FT fd-2'){   #check at exp_35S = 1.0
    new_parms$v_35S[1] = exp_35S
    new_parms$mutants[3] = 0.75
    new_parms$init <- c(10, 0.1, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == 'tfl1-1 fd-2 fdp-1'){  #check
    new_parms$mutants[2] = 1
    new_parms$mutants[3] = .95
    new_parms$init <- c(0, 0, 0.005, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:TFL1 fd-2 fdp-1'){  #check
    new_parms$v_35S[2] = exp_35S
    new_parms$mutants[3] = .95
    new_parms$init <- c(0, 10.1, 0.005, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:FT fd-2 fdp-1'){  #check, regardless of exp_35S
    new_parms$v_35S[1] = exp_35S
    new_parms$mutants[3] = .95
    new_parms$init <- c(10, 0.1, 0.005, 0.1, 0, 0, 1)
  }
  return(new_parms)
}

#Predicts leaf number for various genotypes
predict_genotype = function(genotype,parms){
  new_parms = genotype_parms(genotype,parms)
  new_parms$init[1:5] = new_parms$init[1:5]
  return(predict_leaves(new_parms))
}

#Score for experimental flowering time - predicted flowering time
counter <- 0

obj_fun4_helper <- function(params) {
  op_parms <<- params
  j <- as.list(params)
  names(j) <- names(op_parms5)
  counter <<- counter + 1
  # k <- length(subset(params, c(params < 0.001, params > 10))) * 100
  data_model <- read.delim('Jaeger_data_New.csv',sep=',')
  data_model$pred_R = NA
  data_model$pred_C = NA
  
  for(gen in data_model$Genotype){
    i = data_model$Genotype == gen
    pred = predict_genotype(gen,c(j,baa2))*time_scale
    data_model$pred_R[i] = pred[1]
    data_model$pred_C[i] = pred[2]-pred[1]
  }
  
  data_model[is.na(data_model)] <- 50
  score <- sum((data_model$Ros_Exp - data_model$pred_R)^2 + (data_model$Caul_Exp - data_model$pred_C)^2)
  out_line <- as.vector(unlist(c(j, score, counter)))
  data_out <- as.data.frame(t(unname(out_line)))
  write.table(data_out, file = sprintf("full_run%f.csv", run), append = TRUE, sep = ",", col.names = FALSE)
  print(out_line)
  return(score)
}

out <- GenSA(op_parms, obj_fun4_helper, lower = low, upper = upp, control = list(temperature = 100, maxit = 1000000))
dat <- data.frame(out)
data_final <- cbind(names(op_parms5), op_parms, dat)
write.csv(data_final, file = sprintf("run_%f.csv", run))
