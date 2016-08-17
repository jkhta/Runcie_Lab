library(deSolve)
library(GenSA)
library(Rcpp)

setwd("/Users/James/Desktop/Jaeger Model/")

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
    rho[2,2] = T_f^h_5_2 / (T_f^h_5_2  + X[5]^h_5_2); rho[2,3] = 0; rho[2,1] = 1# - rho[2,2] # should this be in there? Otherwise it doesn't sum to 1; 
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
time_scale2 = 3
init = c(0,0,0.1,0.1,0,0,1) # Starts with some FD and LFY, and with a leaf production rate = 1 per unit t.
t_run = seq(0,200,by=0.1)

v_lfy = 0.05;v_ap1 = 0.05

op_parms <- list(
  K_13 = 0.46676042,
  K_23 = 3.25384158,
  K_4_3 = 0.57875631,
  K_23_4 = 9.37669999,
  K_13_4 = 0.27925075,
  K_23_5 = 0.01985732,
  K_13_5 = 0.04917293,
  K_4_5 = 0.15718656,
  K_5_4 = 0.56239711,
  h_4_3 = 4.00,
  h_23_4 = 3.84970009,
  h_13_4 = 4.00,
  h_23_5 = 4.00,
  h_13_5 = 1.88840037,
  h_4_5 = 3.93008488,
  h_5_4 = 3.66735890,
  h_5_2 = 1.02323239
)

op_parms2 <- list(
  K_13 = 0.30332850,
  K_23 = 3.64240056,
  K_4_3 = 0.37486165,
  K_23_4 = 9.37669977,
  K_13_4 = 0.18088952,
  K_23_5 = 0.01289565,
  K_13_5 = 0.03188165,
  K_4_5 = 3.62538877,
  K_5_4 = 0.36426676,
  h_4_3 = 3.29947535,
  h_23_4 = 3.84969937,
  h_13_4 = 3.73071106,
  h_23_5 = 4.00000000,
  h_13_5 = 2.63244114,
  h_4_5 = 2.54532001,
  h_5_4 = 2.37516785,
  h_5_2 = 2.07212274
)

op_parms3 <- list(
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

op_parms4 <- list(
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

op_parms5 <- list(
  K_13 = 8.008541e-01,
  K_23 = 6.013712e+00,
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

parms_ori = list(
  K_13 = 0.482188814,
  K_23 = 2.743194213,
  K_4_3 = 0.311310135,
  K_23_4 = 7.443649784,
  K_13_4 = 0.025,
  K_23_5 = 0.042090169,
  K_13_5 = 0.025,
  K_4_5 = 0.13096201,
  K_5_4 = 0.279180089,
  h_4_3 = 4.00,
  h_23_4 = 3.897157781,
  h_13_4 = 2.89816155,
  h_23_5 = 4,
  h_13_5 = 1.830372619,
  h_4_5 = 3.714213437,
  h_5_4 = 4,
  h_5_2 = 1.019945776,
  
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

parms_high = list(
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
        
        delta    = c(time_scale2*2*c(0.05,0.05,0.05,v_lfy,v_ap1),0,0),
        v_35S    = time_scale2*c(0,rep(0,4)),		#check this. 
        v1       = time_scale2*1*c(rep(0.01,4),0),
        v2       = time_scale2*1*c(0.05,0.05,0.05,v_lfy,v_ap1),
        v3       = time_scale2*2*c(0.05,0.05,0.05,v_lfy,v_ap1),
        eta_leaf = time_scale2*0.01, #eta_leaf = time_scale*0.01
        
        T_f = 0.2,
        
        mutants = rep(0,5),
        
        repression = 1
)

parms_high$init <- init
low_parms <- list(
  K_13 = 0.30332850,
  K_23 = 3.64240056,
  K_4_3 = 0.37486165,
  K_23_4 = 9.37669977,
  K_13_4 = 0.18088952,
  K_23_5 = 0.01289565,
  K_13_5 = 0.03188165,
  K_4_5 = 3.62538877,
  K_5_4 = 0.36426676,
  h_4_3 = 3.29947535,
  h_23_4 = 3.84969937,
  h_13_4 = 3.73071106,
  h_23_5 = 4.00000000,
  h_13_5 = 2.63244114,
  h_4_5 = 2.54532001,
  h_5_4 = 2.37516785,
  h_5_2 = 2.07212274,
  
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

low_parms$init <- init

#Parameters for optimization

low <- c(rep(0.01, 9), rep(1, 8))
upp <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 4, 4, 4, 4, 4, 4, 4, 4)

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

root_fun = function(t,y,parms,...){
  # This tells the ODE solver to trigger an event. It returns a vector. Events are triggered each time an element = 0.
  return(c(y[5]-0.2,y[5]-0.3)) 
}

# Can use eventsdat to specific changes in pararmeters at specific points in time
# ex. change in FT at certain time
eventsdat = data.frame(var=c(1,2),time=10,value=c(1,0),method='rep')

# eventsfun is called whenever a root is reached 
eventsfun = function(t,y,parms,...){
  if(y[5] > 0.2) y[7] = 0
  y
}
terminalroot = 3 # The 2nd root causes the simulation to stop

fit_model_original = function(parms){
  s1 <- ode(y = c(parms$init),
            times = t_run,
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
            method='bdf',
            rootfun = root_fun,
            events = list(func = eventsfun,root=T,terminalroot=terminalroot))
  
  return(s1)
}

predict_leaves = function(parms){
  s1 = fit_model_original(parms)
  return(attributes(s1)[['troot']])
}

exp_35S = 01

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
  if(genotype == 'lfy-12'){   #check
    new_parms$mutants[4] = 1
  }
  if(genotype == 'ft-10'){	#check
    new_parms$mutants[1] = 1
  }
  if(genotype == 'tfl-1'){   #check
    new_parms$mutants[2] = 1
  }
  if(genotype == 'fd-2'){    #check
    new_parms$mutants[3] = 0.81
  }
  if(genotype == 'fdp-1'){   #check
    new_parms$mutants[3] = 0.3
  }
  if(genotype == 'fd-2 fdp-1'){   #check
    new_parms$mutants[3] = 0.95
  }
  if(genotype == '35S:TFL1 fd-2'){  #check, exp_35S = 1
    new_parms$v_35S[2] = exp_35S
    new_parms$mutants[3] = 0.81
  }
  if(genotype == 'tfl1-1 fd-2'){   #check
    new_parms$mutants[2] = 1
    new_parms$mutants[3] = 0.81
  }
  if(genotype == '35S:FT fd-2'){   #check at exp_35S = 1.0
    new_parms$v_35S[1] = exp_35S
    new_parms$mutants[3] = 0.81
  }
  if(genotype == 'tfl1-1 fd-2 fdp-1'){  #check
    new_parms$mutants[2] = 1
    new_parms$mutants[3] = .95
  }
  if(genotype == '35S:TFL1 fd-2 fdp-1'){  #check
    new_parms$v_35S[2] = exp_35S
    new_parms$mutants[3] = .95
  }
  if(genotype == '35S:FT fd-2 fdp-1'){  #check, regardless of exp_35S
    new_parms$v_35S[1] = exp_35S
    new_parms$mutants[3] = .95
  }
  return(new_parms)
}

#Predicts leaf number for various genotypes
predict_genotype = function(genotype,parms){
  new_parms = genotype_parms(genotype,parms)
  new_parms$init[1:5] = new_parms$init[1:5] * (1-new_parms$mutants)
  return(predict_leaves(new_parms))
}

## Run the model.

new_parms = parms
x=seq(0,50,by=.1)
terminalroot = 3

genotype = 'Col'
new_parms = genotype_parms(genotype,parms)
# new_parms$init = c(0,0.6,.1,0.1,0,0,1) # Modifications for figure SF2

#Plotting protein concentration curves over time
s1 <- fit_model_original(parms_ori)
s2 <- fit_model_new(raa)
cols = c('red','blue','black','green','gray')

x = t_run
x=x*time_scale
plot(NA,NA,xlim = range(x),ylim = c(0,1.1))#range(s1[,-1]))
for(i in 2:6){
  lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
}
legend('topright',legend=c('FT','TFL1','FD','LFY','AP1'),col=cols,lty=1, cex = 0.75)

plot(x = "Time (Leaves)", y = "Protein Concentration", xlim = range(x),ylim = c(0,1.1), lwd = 2)#range(s1[,-1]))
for(i in 2:6){
  lines(s2[,1]*time_scale,s2[,i],col=cols[i-1])
}
abline(h = 0.2, lty = 2)
legend('topright',legend=c('FT','TFL1','FD','LFY','AP1'),col=cols,lty=1, cex = 0.75)
title(main = "Degradation Rates: 30%")
#Flowering time prediction for various genotypes
setwd("/Users/jkhta/Runcie_Lab/flowering_model/Experimental:Model Data/")
data_model = read.delim('Jaeger_data_original.csv',sep=',')
data_model$pred_R = NA
data_model$pred_C = NA

for(gen in data_model$Genotype){
  i = data_model$Genotype == gen
  pred = predict_genotype(gen,parms_ori)*time_scale
  data_model$pred_R[i] = pred[1]
  data_model$pred_C[i] = pred[2]-pred[1]
}

data_model[is.na(data_model)] <- 50
plot_ori_R <- plot(data_model$Ros_Exp,data_model$Ros_Mod, xlab = "Ros_Exp", ylab = "Ros_Mod");abline(0,1)
plot(data_model$Ros_Exp,data_model$pred_R, xlab = "Ros_Exp", ylab = "Ros_Mod")
plot_ori_C <- plot(data_model$Caul_Exp,data_model$Caul_Mod, xlab = "Ros_Exp", ylab = "Ros_Mod");abline(0,1)
plot(data_model$Caul_Exp,data_model$pred_C, xlab = "Ros_Exp", ylab = "Ros_Mod");abline(0,1)

#Optimization
counter <- 0
#Gets best parameters for AP1 curve fit
obj_fun = function(params) {
  # params <- unlist(j)
  op_parms <<- params
  counter = counter + 1
  moo <- sum((fit_model(c(params, baa))[,6] - e2)^2)
  print(c(params,moo,counter))
  return(moo)
}

#Gets best parameters for all curves
counter <- 0
obj_fun2 = function(params) {
  # params <- unlist(j)
  counter <<- counter + 1
  op_parms <<- params
  moo <- sum((fit_model_new(c(params, baa))[,2:6] - e2)^2)
  print(c(params,moo,counter))
  return(moo)
}

data_model <- read.delim('Jaeger_data_New.csv',sep=',')
data_model$pred_R = NA
data_model$pred_C = NA

#Score for experimental flowering time - predicted flowering time
counter <- 0
obj_fun3_helper <- function(params) {
  op_parms4 <<- params
  counter <<- counter + 1
  data_model <- read.delim('Jaeger_data_New.csv',sep=',')
  data_model$pred_R = NA
  data_model$pred_C = NA
  for(gen in data_model$Genotype){
    i = data_model$Genotype == gen
    pred = predict_genotype(gen,c(params,baa))*time_scale
    data_model$pred_R[i] = pred[1]
    data_model$pred_C[i] = pred[2]-pred[1]
  }
  
  data_model[is.na(data_model)] <- 50
  score <- sum((data_model$Ros_Exp - data_model$pred_R)^2 + (data_model$Caul_Exp - data_model$pred_C)^2)
  return(score)
}

obj_fun3 <- function(params) {
  moo <- obj_fun3_helper(params)
  print(c(params, moo, counter))
  return(moo)
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
        print(c(params, score, counter))
        return(score)
}

obj_fun4 <- function(params) {
        moo <- obj_fun4_helper(params)
        print(c(params, moo, counter))
        return(moo)
}

e1 <- fit_model(c(op_parms2, baa))[,2:6]
e2 <- fit_model_new(parms_ori)[,2:6]

optim(op_parms5, obj_fun4_helper, method = "SANN", lower = low, upper = upp)#,control = list(maxit= 100))

GenSA(k, obj_fun4_helper, lower = low, upper = upp)

DEoptim(obj_fun2, lower = low, upper = upp)


raa <- c(
  op_parms4,
  baa2
)

op_parms4[]
