#test mode to do odes in a data frame

library(deSolve)
rm(list = ls())

t_run <- seq(0, 300, by = 0.1)
time_scale = 1
init <- c(0.1,0)
inf_info <- data.frame()
y <- c(0.5,0,299,0.6,0) #SFT flux into meristem, SFT concentration, FA concentration, time initiated, sink power, flowering status
inf_info <<- rbind(inf_info, y)

parms_list <- list(
  K_1_2 = 0.02,
  h_1_2 = 2,
  delta = time_scale*2*c(0.05,0.05),
  mutants = rep(0,2)
)
parms_list$init <- init

#creates an inflorescence (row) in data frame whenever SFT exceeds the arbitrary threshold

inf_maturation <- function(t,X, parms = NULL,...){
  with(as.list(parms),{
    p_SFT_FA = X[1]^h_1_2/(K_1_2^h_1_2 + X[1] * h_1_2)
    FA = p_SFT_FA * 0.05
    derivatives = rep(0,length(X))
    derivatives[1:2] = c(SFT, FA)
    derivatives = derivatives - delta * X
    return(list(
      Derivatives = derivatives,
      globals = c(
        p_SFT_FA
      )
    ))
  })
}


inflorescence_check <- function(values){
  j <- seq(as.numeric(values[3]),max(t_run), by = 0.1)
  s2 <- ode(y = as.numeric(values[1:2]),
            times = j,
            func = inf_maturation,
            parms = parms_list,
            method = 'lsoda')
  return(s2)
}

for(i in 1:nrow(inf_info)){
    s1 <- inflorescence_check(inf_info[i,])
    print(head(s1))
}
s1 <- inflorescence_check(inf_info[1,])
