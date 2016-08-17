library(data.table)
rm(list = ls())
x <- "/Users/jkhta/Data/30%_Degrade_Batch#2"
setwd(x)
temp <- list.files(pattern = "*.csv")

for(i in 1:length(temp)){
  j <- temp[i]
  assign(j,fread(j))
}

data_names <- c("NA", "K_13", "K_23", "K_4_3", "K_23_4", "K_13_4", "K_23_5", 
                "K_13_5", "K_4_5", "K_5_4", "h_4_3", "h_23_4", "h_13_4", 
                "h_23_4", "h_13_5", "h_4_5", "h_5_4", "h_5_2", "Score", 
                "Iterations")

for(i in ls(pattern = "*_run")) {
  j <- eval(as.symbol(i))
  k <- which.min(unlist(j[,(ncol(j)-1), with=FALSE]))
  l <- j[k,]
  write.table(l, "lowest_scores_.csv", append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
}

find_lowest <- read.csv("lowest_scores_100%.csv", header = TRUE, stringsAsFactors = FALSE)
lowest_index <- which.min(find_lowest[,ncol(find_lowest)-1])
lowest_index <- as.numeric(find_lowest[1,])
lowest_parms <- lowest_index[1:17]
lowest_score <- find_lowest[lowest_index,]
lowest_parms <- as.vector(lowest_score[3:19])

#Add predicted values for the different genotypes----------------------

library(deSolve)
library(GenSA)
library(Rcpp)

setwd("/Users/jkhta/Runcie_Lab/flowering_model")
sourceCpp("Rccp_Files/Silly.cpp")

time_scale = 1  # shifts the timescale by scaling both delta and vs and eta_leaf
time_scale2 = 10
init = c(0,0.1,.1,0.1,0,0,1) # Starts with some FD and LFY, and with a leaf production rate = 1 per unit t.
t_run = seq(0, 50, by=0.01)

v_lfy = 0.05;v_ap1 = 0.05
op_parms <- c(sample(seq(0.01, 10, by = 0.01), 9), sample(seq(1, 4, by = 0.01), 8))
init_parms <- op_parms
low <- c(rep(0.01, 9), rep(1, 8))
upp <- c(rep(10, 9), rep(4, 8))

exp_35S = 1

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

parms_LFY = list(
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
  eta_leaf2 = 0.001
)
parms_ori$init <- init

baa$init <- init

parms_LFY_TFL1 = list(
  K_13 = 0.39381,
  K_23 = 3.2556,
  K_4_3 = 0.28203,
  K_23_4 = 9.3767,
  K_13_4 = 0.040555,
  K_23_5 = 0.033666,
  K_13_5 = 0.029081,
  K_4_5 = 0.13032,
  K_5_4 = 0.28606,
  K_13_2 = 0.1,
  h_4_3 = 4.00,
  h_23_4 = 3.8497,
  h_13_4 = 4.00,
  h_23_5 = 4.00,
  h_13_5 = 1.8217,
  h_4_5 = 3.9369,
  h_5_4 = 3.6732,
  h_5_2 = 1.0239,
  h_13_2 = 2,
  eta_leaf2 = time_scale2*0.01
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
names(lowest_parms) <- names(op_parms5)

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
            times = t_run,
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

setwd("/Users/jkhta/Runcie_Lab/flowering_model/Experimental:Model Data/")

data_model_ori = read.delim('Jaeger_data_New.csv',sep=',')
data_model_ori$pred_R = NA
data_model_ori$pred_C = NA

for(gen in data_model_ori$Genotype){
  i = data_model_ori$Genotype == gen
  pred = predict_genotype(gen,c(lowest_score, baa))*time_scale
  data_model_ori$pred_R[i] = pred[1]
  data_model_ori$pred_C[i] = pred[2]-pred[1]
}

names(lowest_parms) <- names(op_parms5)
lowest_parms <- as.list(lowest_parms)
data_model_hi = read.delim('Jaeger_data_original.csv',sep=',')
data_model_hi$pred_R = NA
data_model_hi$pred_C = NA

for(gen in data_model_hi$Genotype){
  i = data_model_hi$Genotype == gen
  pred = predict_genotype(gen,c(lowest_score, baa))*time_scale
  data_model_hi$pred_R[i] = pred[1]
  data_model_hi$pred_C[i] = pred[2]-pred[1]
}

gen = "Col"
setwd(x)
write.csv(data_model_hi,"lowest_predicted.csv")
for(gen in data_model_hi$Genotype){
  genotype <<- gen
  new_parms <- genotype_parms(gen, c(lowest_parms, baa))
  s1 <- fit_model_new(new_parms)
  cols = c('red','blue','black','green','gray')
  x_time <- t
  tiff(filename = sprintf("Graph_%s.tiff", genotype))
  plot(NA,NA,xlim = range(x_time),ylim = c(0,1))#range(s1[,-1]))
  for(i in 2:6){
    lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
  }
  legend('topright',legend=c('FT','TFL1','FD','LFY','AP1'),col=cols,lty=1, cex = 0.75)
  abline(h = 0.2, lty = 2)
  abline(h = 0.3, lty = 2)
  dev.off()
}
s1 <- fit_model_new(parms_ori)
setwd(x)
write.csv(data_model_hi, file = "lowest_predicted.csv")
new_parms <- c(lowest_parms, baa)
s1 <- fit_model_new(c(lowest_parms,baa))
cols = c('red','blue','black','green','gray')
x_time <- t
x_time=x_time*time_scale2
plot(NA,NA,xlim = range(x_time),ylim = c(0,1))#range(s1[,-1]))
for(i in 2:6){
  lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
}
legend('topright',legend=c('FT','TFL1','FD','LFY','AP1'),col=cols,lty=1, cex = 0.75)
abline(h = 0.2, lty = 2)
abline(h = 0.3, lty = 2)
