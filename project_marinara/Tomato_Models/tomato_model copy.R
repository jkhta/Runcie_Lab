library(deSolve)
rm(list = ls())
t <- seq(0, 50, by = 0.01)
time_scale = 1
init <- c(0, 0.6, 0, 0.6, 0.1, 0, 1)
inf_info <- data.frame()

#         Protein Stuff
#         X[1] = SFT
#         X[2] = SP
#         X[3] = FA
#         X[4] = SFT/SP Ratio
#         X[5] = Number of leaves
#         X[6] = Change in number of leaves

parms_list <- list(
  alpha <- 0.66, # proportion of SFT that goes to inflorescence meristem
  K_1_3 <- 0.3,
  K_2_3 <- 0.1,
  K_3_5 <- 0.3,
  K_4_5 <- 0.1,
  h_1_3 <- 4,
  h_2_3 <- 4,
  h_3_5 <- 4,
  h_4_5 <- 4,
  delta    = c(time_scale*2*c(0.05,0.05,0.05,0.05,0.05),0,0),
  v_35S    = time_scale*c(rep(0, 5)),
  eta_leaf = time_scale*0.5,
  mutants = rep(0,5)
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
create_inf <- function(SFT_conc){#,time) {
        if(SFT_conc > 1) {
                data <- c(0, 0, "NOT FLOWERED")#, t)
                names(data) <- c("SFT", "FA", "Status")#, "Time Initiated")
                inf_info <<- rbind(inf_info, data)
        }
}

# inflorescence_model <- function(t,X,parms = NULL,...){
#   with(as.list(parms),{
#     
#   }
# }

tomato_model <- function(t, X, parms=NULL,...) {

        with(as.list(parms),{

        #state variables
        p_1_3 = X[1]^h_1_3 / (K_1_3^h_1_3 + X[1]^h_1_3)
        p_2_3 = X[2]^h_2_3 / (K_2_3^h_1_3 + X[2]^h_2_3)
        p_3_5 = X[3]^h_3_5 / (K_3_5^h_3_5 + X[3]^h_3_5)
        p_4_5 = X[4]^h_4_5 / (K_4_5^h_4_5 + X[4]^h_4_5)
        
        #protein concentration rates based on binding probabilities
        SFT <- X[1]
        SFT_split <- alpha*SFT #up to a maximum amount; using carrying capacity
        SFT <- eta_leaf * X[6] - nrow(inf_info) * SFT_split
        TMF <- 0.5
        inf_exp <- (1 - p_1_3*(1 - p_2_3)) * 0.5 + (p_1_3*(1 - p_2_3) * 2.5)
        SP <- eta_leaf * X[6]
        FA <- 0.005 * X[6]
        
        #inflorescence creator
        create_inf(SFT)#,t)
        # for(i in nrow(inf_info)){
        #   j <- inflorescence_check(inf_info[i,1:2])
        #   inf_info[i,1:2] <- j
        # }
        # 
        derivatives <- rep(0, length(X))
        derivatives[1:5] <- c(SFT, TMF, inf_exp, SP, FA)
        derivatives[1:5] <- derivatives[1:5] * (1 - mutants)
        derivatives[6] <- X[7]
        derivatives <- derivatives - delta*X
        return(list(
                Derivatives <- derivatives,
                globals = c(
                  p_1_3,
                  p_2_3,
                  p_3_5,
                  p_4_5
                )
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

# inflorescence_check <- function(parms){
#   s2 <- ode(y = parms,
#             times = t,
#             func = inflorescence_model,
#             parms = parms,
#             method = 'lsoda')
#   return(s2)
# }

s1 <- fit_model(parms_list)


cols = c('red','blue','black','green','gray')
time_scale <- 1
x=t
x=x*time_scale
plot(NA,NA,xlim = range(x),ylim = c(0,5))#range(s1[,-1]))
for(i in 2:5){
        lines(s1[,1]*time_scale,s1[,i],col=cols[i-1])
}

legend('topright',legend=c('SFT','TMF','INF','SP','FA'),col=cols,lty=1, cex = 0.75)

