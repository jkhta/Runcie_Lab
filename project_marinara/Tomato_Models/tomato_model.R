library(deSolve)

t <- seq(0, 300, by = 0.01)
init <- c(0, 0.5, 0, 0, 1)
inf_info <- data.frame()

#         Protein Stuff
#         X[1] = SFT
#         X[2] = SP
#         X[3] = FA
#         X[4] = SFT/SP Ratio
#         X[5] = Number of leaves
#         X[6] = Change in number of leaves


parms_list <- list(
        k_1 <- 0.05,
        k_2 <- 0.075,
        k_3 <- 0.075,
        # k_SFT <- 0.1,
        # k_SP <- 0.1,
        # k_FA <- 0.1,
        K_13 <- 0.25,
        h_13 <- 2,
        beta <- 0.33,
        gamma <- 0.2,
        delta <- c(0.1, 0.1, 0.1, 0, 0),
        mutants <- c(0, 0, 0)
)

parms_list$init <- init

#Protein concentrations of SFT, SP, FA, change in leaves, leaf number

SFT_threshold <- function(SFT_conc) {
        flower <- 0
        threshold <- 0
        if(SFT_conc < 0.75) {
                flower <- 0
                vegetative <- 1
        }
        else{
                flower <- 1
                vegetative <- 0
        }
        return(c(flower, vegetative))
}

#creates an inflorescence (row) in data frame whenever SFT exceeds the arbitrary threshold 0.3
create_inf <- function(SFT_conc) {
        if(SFT_conc > 0.3) {
                data <- c(0, "NOT FLOWERED")
                inf_info <<- rbind(inf_info, data)
        }
}

#not used right now but when the flowering meristem gets modeled (individual flowers), this will tell when a inflorescence meristem stop being a sink
inflorescence_check <- function(inf_data) {
        for(i in 1:inf_data[i,]){
                if(inf_data[i,1] > 0.3) {
                        inf_data[i,2] <- "FLOWERED"
                }
        }
}

tomato_model <- function(t, X, parms=NULL,...) {
        
        with(as.list(parms),{
          
        #state variables
        p_13 <- X[1]/(K_13^h_13 + X[1]^h_13)
        SFT <- X[1]
        SFT <- k_1*X[4] - beta*nrow(inf_info)*SFT #needs some work because it doesn't make physical sense with no cap; proportional
        # SFT_meristem <- k_SFT*SFT*(1 - SFT)
        SP <- k_2*X[4]*SFT_threshold(SFT)[2]
        # SP_meristem <- k_SP*SP*(1 - SP)
        FA <- k_3*X[4] + gamma*p_13 - beta*nrow(inf_info)
        # FA_cap <- k_FA*FA*(1 - FA)
        
        # flower_ratio <- SFT/SP
        create_inf(SFT)
        # inflorescence_check(inf_infO)
        
        derivatives <- rep(0, length(X))
        derivatives[1:3] <- c(SFT, SP, FA)#, flower_ratio)
        derivatives[1:3] <- derivatives[1:3] * (1 - mutants)
        derivatives[4] <- X[5]
        derivatives <- derivatives - delta*X
        return(list(
                Derivatives <- derivatives
        ))
        })
}

fit_model = function(parms){
        s1 <- ode(y = c(parms$init),
                  times = t,
                  func = tomato_model,
                  parms=parms,
                  method='lsoda')
                  # rootfun = root_fun,
                  # events = list(func = eventsfun,root=T,terminalroot=terminalroot))
        return(s1)
}

s1 <- fit_model(parms_list)

cols = c('red','blue','black')#'green','gray')
time_scale <- 1
x=t
x=x*time_scale
plot(NA,NA,xlim = range(x),ylim = c(0,5))#range(s1[,-1]))
for(i in 2:4){
        lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
}
