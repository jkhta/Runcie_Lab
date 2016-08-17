library(deSolve)
library(GenSA)
library(Rcpp)

setwd("/Users/jkhta/Runcie_Lab/flowering_model/")

sourceCpp("Silly.cpp")

jaeger_model = function(t,X,parms=NULL,...){
  # X: protein concentrations and developmental state 
  # X[1]: FT
  # X[2]: TFL1
  # X[3]: FD
  # X[4]: LFY
  # X[5]: AP1
  # X[6]: Num. Leaves
  # X[7]: increase in Num. Leaves. Set to 0 once AP1 > 0.2
  
  with(as.list(parms),{
    # recover()
    #Hub protein-Protein Binding
    #FT:FD
    x_13 = K_23*X[1]*X[3]/(K_13*K_23 + K_13*X[2]+K_23*X[1])
    #TFL1:FD
    x_23 = K_13*X[2]*X[3]/(K_13*K_23 + K_13*X[2]+K_23*X[1])
    
    #Hub Gene Activation
    #LFY -> FD
    p_4_3 = X[4]^h_4_3 / (K_4_3^h_4_3 + X[4]^h_4_3)
    #FT:FD -> LFY
    p_13_4 = K_23_4^h_23_4 * x_13^h_13_4 / (K_13_4^h_13_4 * K_23_4^h_23_4 + K_23_4^h_23_4*x_13^h_13_4 + K_13_4^h_13_4*x_23^h_23_4)
    #TFL1:FD -> LFY
    p_23_4 = K_13_4^h_13_4 * x_23^h_23_4 / (K_13_4^h_13_4 * K_23_4^h_23_4 + K_23_4^h_23_4*x_13^h_13_4 + K_13_4^h_13_4*x_23^h_23_4)
    #AP1 -> LFY
    p_5_4 = X[5]^h_5_4 / (K_5_4^h_5_4 + X[5]^h_5_4)
    #FT:FD -> AP1
    p_13_5 = K_23_5^h_23_5 * x_13^h_13_5 / (K_13_5^h_13_5 * K_23_5^h_23_5 + K_23_5^h_23_5*x_13^h_13_5 + K_13_5^h_13_5*x_23^h_23_5)
    #TFL1:FD -> AP1
    p_23_5 = K_13_5^h_13_5 * x_23^h_23_5 / (K_13_5^h_13_5 * K_23_5^h_23_5 + K_23_5^h_23_5*x_13^h_13_5 + K_13_5^h_13_5*x_23^h_23_5)
    #LFY -> AP1
    p_4_5 = X[4]^h_4_5 / (K_4_5^h_4_5 + X[4]^h_4_5)
    
    
    rho = matrix(0,nr=5,nc=3)
    rho[,1] = 1
    
    rho[1,1] = 1; rho[1,2] = 0; rho[1,3] = 0
    rho[2,2] = T_f^h_5_2 / (T_f^h_5_2  + X[5]^h_5_2); rho[2,3] = 0; rho[2,1] = 1 - rho[2,2] # should this be in there? Otherwise it doesn't sum to 1; 
    rho[3,2] = p_4_3; rho[3,1] = 1 - rho[3,2]; rho[3,3] = 0
    rho[4,1] = (1 - p_13_4 - p_23_4) * (1-p_5_4); rho[4,2] = p_13_4*(1-p_5_4) + (1 - p_13_4 - p_23_4)*p_5_4; rho[4,3] = p_13_4*p_5_4
    rho[5,1] = (1 - p_13_5 - p_23_5) * (1-p_4_5); rho[5,2] = p_13_5*(1-p_4_5) + (1 - p_13_5 - p_23_5)*p_4_5; rho[5,3] = p_13_5*p_4_5
    
    v = matrix(0,nr=5,nc=3)
    v[,1] = c(eta_leaf * X[6], v1[2:5])
    v[,2] = v2
    v[,3] = v3  # check this. c(NA,rep(0,1,4))
    
    v = v_35S + rowSums(rho * v)
    
    V = rep(0,length(X))
    V[1:5] = v*(1-mutants) 	# If a gene is mutated, set it's change = 0
    V[6] = X[7]		# change in # leaves
    V = V - delta*X  # subtract degradation
    return(list(
      Derivitaves = V,
      globals = c(
        x_13=x_13,
        x_23=x_23,
        p_13_4 = p_13_4,
        p_13_5 = p_13_5,
        p_23_4 = p_23_4,
        p_23_5 = p_23_5,
        p_4_3 = p_4_3,
        p_4_5 = p_4_5,
        p_5_4 = p_5_4
      )
    ))
  })
}

time_scale = 1  # shifts the timescale by scaling both delta and vs and eta_leaf

init = c(0, 0.6, 0.1, 0.1, 0, 0, 1) # Starts with some FD and LFY, and with a leaf production rate = 1 per unit t.
t = seq(0,200,by=0.1)

v_lfy = 0.05;v_ap1 = 0.05

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
  v_35S    = time_scale*c(rep(0, 5)),		#check this. 
  v1       = time_scale*1*c(0, 0.01, 0.01, 0.01, 0),
  v2       = time_scale*1*c(0.05,0.05,0.05,v_lfy,v_ap1),
  v3       = time_scale*2*c(0.05,0.05,0.05,v_lfy,v_ap1),
  eta_leaf = time_scale*0.01, #eta_leaf = time_scale*0.01
  
  T_f = 0.2,
  
  mutants = rep(0,5)
)

parms_ori$init <- init

parms_optimize = list(
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
  v_35S    = time_scale*c(rep(0, 5)),		#check this. 
  v1       = time_scale*1*c(0, 0.01, 0.01, 0.01, 0),
  v2       = time_scale*1*c(0.05,0.05,0.05,v_lfy,v_ap1),
  v3       = time_scale*2*c(0.05,0.05,0.05,v_lfy,v_ap1),
  eta_leaf = time_scale*0.01, #eta_leaf = time_scale*0.01
  
  T_f = 0.2,
  
  mutants = rep(0,5)
  
)

parms_ori$init <- init

root_fun = function(t,y,parms,...){
  # This tells the ODE solver to trigger an event. It returns a vector. Events are triggered each time an element = 0.
  return(c(y[5]-0.2,y[5]-0.3)) 
}

# Can use eventsdat to specific changes in pararmeters at specific points in time
# ex. change in FT at certain time
eventsdat = data.frame(var=c(1,2),time=10,value=c(1,0),method='rep')

# eventsfun is called whenever a root is reached 
eventsfun = function(time = t,y = parms$init,parms,...){
  if(y[5] > 0.3) y[7] = 0
  y
}

terminalroot = 3 # The 2nd root causes the simulation to stop

fit_model_ori = function(parms){
  s1 <- ode(y = c(parms$init),
            times = t,
            func = jaeger_model,
            parms=parms,
            method='lsoda',
            rootfun = root_fun,
            events = list(func = eventsfun,root=T,terminalroot=terminalroot))
  return(s1)
}

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

exp_35S = 1

#Need to comment out init values for optimizing genotype
genotype_parms = function(genotype,parms){
  new_parms <- parms #parms
  #Col 					#check
  if(genotype == '35S:FT'){ #check
    new_parms$v_35S[1] = exp_35S #1.3-1.8
    # new_parms$init <- c(10, 0.6, 0.1, 0.1, 0, 0, 0)
  }
  if(genotype == '35S:LFY'){
    new_parms$v_35S[4] = exp_35S #nothing gets is fast enough (minimum = 5 Ros leaves)
    # new_parms$init <- c(0, 0.6, 0.1, 10.1, 0, 0, 1)
  }
  if(genotype == '35S:TFL1'){ #not with this parameter
    new_parms$v_35S[2] = exp_35S #1
    # new_parms$init <- c(0, 10.7, 0.1, 0.1, 0, 0, 1)
  }
  if(genotype == 'lfy-12'){   #check
    new_parms$mutants[4] = 1
    # new_parms$init <- c(0, 0.6, 0.1, 0, 0, 0, 1)
  }
  if(genotype == 'ft-10'){	#check
    new_parms$mutants[1] = 1
    # new_parms$init <- c(0, 0.6, 0.1, 0.1, 0, 0, 1)
  }
  if(genotype == 'tfl-1'){   #check
    new_parms$mutants[2] = 1
    # new_parms$init <- c(0, 0, 0.1, 0.1, 0, 0, 1)
  }
  if(genotype == 'fd-2'){    #check
    new_parms$mutants[3] = 0.75
    # new_parms$init <- c(0, 0.6, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == 'fdp-1'){   #check
    new_parms$mutants[3] = 0.2
    # new_parms$init <- c(0, 0.6, 0.08, 0.1, 0, 0, 1)
  }
  if(genotype == 'fd-2 fdp-1'){   #check
    new_parms$mutants[3] = 0.95
    # new_parms$init <- c(0, 0.6, 0.005, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:TFL1 fd-2'){  #check, exp_35S = 1
    new_parms$v_35S[2] = exp_35S
    new_parms$mutants[3] = 0.75
    # new_parms$init <- c(0, 10.7, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == 'tfl1-1 fd-2'){   #check
    new_parms$mutants[2] = 1
    new_parms$mutants[3] = 0.75
    # new_parms$init <- c(0, 0, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:FT fd-2'){   #check at exp_35S = 1.0
    new_parms$v_35S[1] = exp_35S
    new_parms$mutants[3] = 0.75
    # new_parms$init <- c(10, 0.6, 0.025, 0.1, 0, 0, 1)
  }
  if(genotype == 'tfl1-1 fd-2 fdp-1'){  #check
    new_parms$mutants[2] = 1
    new_parms$mutants[3] = .95
    # new_parms$init <- c(0, 0.6, 0.005, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:TFL1 fd-2 fdp-1'){  #check
    new_parms$v_35S[2] = exp_35S
    new_parms$mutants[3] = .95
    # new_parms$init <- c(0, 10.7, 0.005, 0.1, 0, 0, 1)
  }
  if(genotype == '35S:FT fd-2 fdp-1'){  #check, regardless of exp_35S
    new_parms$v_35S[1] = exp_35S
    new_parms$mutants[3] = .95
    # new_parms$init <- c(10, 0.1, 0.005, 0.1, 0, 0, 1)
  }
  return(new_parms)
}

#Predicts leaf number for various genotypes
predict_genotype = function(genotype,parms){
  new_parms = genotype_parms(genotype,parms)
  new_parms$init[1:5] = new_parms$init[1:5] #* (1-new_parms$mutants)
  return(predict_leaves(new_parms))
}

## Run the model.

new_parms <- parms_ori
x=seq(0,50,by=.1)

genotype = 'Col'
new_parms = genotype_parms(genotype,parms_ori)
# new_parms$init = c(0,0.6,.1,0.1,0,0,1) # Modifications for figure SF2

#Plotting protein concentration curves over time
s1 <- fit_model_ori(parms_ori)
cols = c('red','blue','black','green','gray')
x=x*time_scale
plot(NA,NA,xlim = range(x),ylim = c(0,1))#range(s1[,-1]))
for(i in 2:6){
  lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
}
legend('topright',legend=c('FT','TFL1','FD','LFY','AP1'),col=cols,lty=1, cex = 0.75)
abline(h = 0.2, lty = 2)

#Flowering time prediction for various genotypes

data_model = read.delim('Jaeger_data_New.csv',sep=',')
data_model$pred_R = NA
data_model$pred_C = NA

for(gen in data_model$Genotype){
  i = data_model$Genotype == gen
  pred = predict_genotype(gen,parms_ori)*time_scale
  data_model$pred_R[i] = pred[1]
  data_model$pred_C[i] = pred[2]-pred[1]
}

score_model <- function(data){
  i <- with(data, (Ros_Exp - pred_R)^2 + (Caul_Exp - pred_C)^2)
  j <- with(data, (Ros_Model - pred_R)^2 + (Caul_Model - pred_C)^2)
  print(sum(i))
  print(sum(j))
}

score_model(data_model)

data_model[is.na(data_model)] <- 50
plot_ori_R <- plot(data_model$Ros_Exp,data_model$Ros_Mod, xlab = "Ros_Exp", ylab = "Ros_Mod");abline(0,1)
plot(data_model$Ros_Exp,data_model$pred_R, xlab = "Ros_Exp", ylab = "Ros_Mod")
plot_ori_C <- plot(data_model$Caul_Exp,data_model$Caul_Mod, xlab = "Ros_Exp", ylab = "Ros_Mod");abline(0,1)
plot(data_model$Caul_Exp,data_model$pred_C, xlab = "Ros_Exp", ylab = "Ros_Mod");abline(0,1)

score_model <- function(data){
  i <- with(data, (Ros_Exp - pred_R)^2 + (Caul_Exp - pred_C)^2)
  j <- with(data, (Ros_Model - pred_R)^2 + (Caul_Model - pred_C)^2)
  print(sum(i))
  print(sum(j))
}

score_model(data_model)


obj_fun <- function(initial_parms) {
  initial <- initial_parms
  parameters <- parms_optimize
  parameters$init <- initial
  
  pred <- predict_genotype("35S:TFL1 fd-2", parameters)*time_scale
  pred_R <- pred[1]
  pred_C <- pred[2] - pred[1]
  
  score <- (24.4 - pred_R)^2 + (5.4 - pred_C)^2
  print(c(initial, pred_R, pred_C, score))
  return(score)
}

data_model$pred_R[1]
low <- c(0, 0, 0.005, 0, 0, 0, 0.9999999999)
upp <- c(0.0000000001, 1, 0.0050000001, 1, 1e-10, 1e-10, 1)

out_35S_FT <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.0001))
out_35S_LFY <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.0001))
out_35S_TFL1 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.0001))
out_lfy_12 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.0001))
out_ft_10 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.05))
out_tfl_1 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.0001))
out_fd_2 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.05))
out_fdp_1 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.0001))
out_fd_2_fdp_1 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.05))
out_fdp_1 <- GenSA(init, obj_fun, lower = low, upper = upp, control = list(threshold.stop = 0.0001))


dat <- data.frame(out)
data_final <- cbind(names(op_parms5), op_parms, dat)
write.csv(data_final, file = sprintf("run_%f.csv", run))
