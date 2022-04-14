
# loading packages --------------------------------------------------------

library(deSolve)
library(parallel)
library(doParallel)
library(openxlsx)
library(tidyr)
library(dplyr)

## change work patch
setwd('1_simulate/')

set.seed(202203)

# function ----------------------------------------------------------------

## estimate beta basing contact matrix (M) and basic reproductive number (R0)
compute.beta <- function(R0,constant,M){
  Eig = max(Re(eigen(M)$values))
  return(R0*constant/Eig)
}

# parameter ---------------------------------------------------------------

## global option -----------------------------------------------------------

## MCMC simulate times
n <- 1000

## simulate times
t_end <- 300 ## time end

## time of intervention
times_lockdown <- c(seq(5, 20, 2), t_end)

## get city list
city_list <- openxlsx::getSheetNames('data/demo.xlsx')[1:2]

## epidemiological parameters
gamma_2 <- 1/7 ## communicability period of symptomatic infections
gamma_1 <- 1/5 ## communicability period of asymptomatic infections
percent_asym <- 0.31 ## proportion of asymptomatic infections
RP_numbers <- c(6, 8, 10, 12)  ## basic reproductive number

para_0_ves <- 0 ## VES of unvaccine
para_1_ves <- 0.5267 ## VES of full vaccine
para_2_ves <- 0.6 ## VES of booster vaccine 

########################## drop #################################
rb_is <- 1   ## relative beta
rb_ia <- 1   ## relative beta
rb_ip <- 1   ## relative beta
#################################################################

df_para <- read.xlsx('data/parameter.xlsx', sheet = 'gamma') |> select(var, value_shape, value_rate)
df_asym <- read.xlsx('data/parameter.xlsx', sheet = 'ratio')
df_asym <- df_asym[,'value']

# model -------------------------------------------------------------------

run.model <- function(i, times_lockdown, p){
  ### loading omega and beta value
  p[52:56] <- as.numeric(mcmc_para[i,])
  ### adjust bate value by relative infectivity of infections
  p[7:9] <- p[55]*p[7:9]
  ### run model by each times of lockdown
  models <- lapply(times_lockdown, FUN = function(t_lockdown){
    p[57] <- t_lockdown
    model <- ode(u0, times, SEIR_extend, p)
    model <- modifiy.data(round(model), age_group)
    #### add identify of each ode progress
    model$n <- i
    model$t_lockdown <- t_lockdown
    model$rp <- p[56]
    return(model)
  })
  models <- do.call('rbind', models)
  return(models)
}

cl <- makeCluster(120)
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(deSolve)
  ## extend SEIR model
  SEIR_extend <- function(time, u, p){
    
    d <- 1
    # beta <- p[55]
    # beta <- ifelse(time<p[56], p[55], p[55]*(1-p[57]))
    if(time >= p[57]){
      if(p[58] == 0){
        d <- 1/0.25
      } else if (p[58] == 1){
        contact_matrix <- contact_matrix_after
      } else {
        contact_matrix <- contact_matrix_after
        d <- 1/0.25
      }
    }
    
    S <- matrix(u[1:age_group], 1, age_group)
    E <- matrix(u[(age_group+1):(age_group*2)], 1, age_group)
    I_a <- matrix(u[(age_group*2+1):(age_group*3)], age_group, 1)
    I_p <- matrix(u[(age_group*3+1):(age_group*4)], age_group, 1)
    I_s <- matrix(u[(age_group*4+1):(age_group*5)], age_group, 1)
    R <- matrix(u[(age_group*5+1):(age_group*6)], 1, age_group)
    
    flow_infect_a <- colSums(I_a %*% (S*p[59:100]/p[10:51]) * contact_matrix)*(1-p[rep(4:6, times = age_group/3)])*p[7]
    flow_infect_p <- colSums(I_p %*% (S*p[59:100]/p[10:51]) * contact_matrix)*(1-p[rep(4:6, times = age_group/3)])*p[8]
    flow_infect_s <- colSums(I_s %*% (S*p[59:100]/p[10:51]) * contact_matrix)*(1-p[rep(4:6, times = age_group/3)])*p[9]
    
    flow_s_e <- matrix(c(flow_infect_a, flow_infect_p, flow_infect_s), nrow = 3, byrow = T) |> colSums()
    flow_e_ia <- E * p[3] * p[52]
    flow_e_ip <- E * (1-p[3]) * p[53]
    flow_ip_is <- I_p * p[54]
    flow_ia_r <- I_a * p[1] * d
    flow_is_r <- I_s * p[2] * d
    
    dS <- -flow_s_e
    dE <- flow_s_e - flow_e_ia - flow_e_ip
    dIa <- flow_e_ia - t(flow_ia_r)
    dIp <- flow_e_ip - t(flow_ip_is)
    dIs <- flow_ip_is - flow_is_r
    dR <- flow_ia_r + flow_is_r
    
    xIa <- sum(flow_e_ia)
    # xIp <- sum(flow_e_ip)
    xIs <- sum(flow_ip_is)
    
    return(list(c(dS, dE, dIa, dIp, dIs, dR, xIa, xIs)))
  }
  
  ## function to merge each age group of model compartment
  modifiy.data <- function(u, age_group){
    time <- u[,1]
    S <- rowSums(u[,2:age_group])
    E <- rowSums(u[,(age_group+2):(age_group*2+1)])
    I_a <- rowSums(u[,(age_group*2+2):(age_group*3+1)])
    I_p <- rowSums(u[,(age_group*3+2):(age_group*4+1)])
    I_s <- rowSums(u[,(age_group*4+2):(age_group*5+1)])
    R <- rowSums(u[,(age_group*5+2):(age_group*6+1)])
    
    # xIa <- rowSums(u[,(age_group*6+2):(age_group*7+1)])
    # xIp <- rowSums(u[,(age_group*7+2):(age_group*8+1)])
    # xIs <- rowSums(u[,(age_group*8+2):(age_group*9+1)])
    xIa <- u[,age_group*6+2]
    # xIp <- u[,age_group*6+3]
    xIs <- u[,age_group*6+3]
    
    outcome <- data.frame(time, S, E, I_a, I_p, I_s, R, xIa, xIs)
    outcome$dIa <- c(0, diff(outcome$xIa))
    # outcome$dIp <- c(0, diff(outcome$xIp))
    outcome$dIs <- c(0, diff(outcome$xIs))
    outcome$dI <- outcome$dIs + outcome$dIa
    outcome$xI <- outcome$xIs + outcome$xIa
    # outcome <- cbind(outcome, u[,(age_group*6+2):(age_group*9+1)])
    return(outcome)
  }
})

for (c in 1:2) {
  ## local option -----------------------------------------------------------
  source('8.1_city_setting.R')
  age_group <- length(pop_s)
  # u0 <- c(pop_s, pop_e, pop_ia, pop_ip, pop_is, pop_r, rep(0, age_group*3))
  u0 <- c(pop_s, pop_e, pop_ia, pop_ip, pop_is, pop_r, rep(0, 2))
  times <- 0:t_end
  contact_matrix <- contact_mt_before[rep(1:(age_group/3), each = 3), rep(1:(age_group/3), each = 3)]
  contact_matrix_after <- contact_mt_after[rep(1:(age_group/3), each = 3), rep(1:(age_group/3), each = 3)]
  
  ## only lockdown -----------------------------------------------------------
  p0 <- c(gamma_1, ## asymptomatic infections
          gamma_2, ## symptomatic infections
          percent_asym,
          para_0_ves, para_1_ves, para_2_ves,
          rb_ia, rb_ip, rb_is, 
          pop_n+1, 
          as.numeric(mcmc_para[1,]), 
          times_lockdown[1],
          1,## plan 1: only lockdown
          rep(df_asym, each = 3)) 
  clusterExport(cl, c('u0', 'p0', 'times', 'mcmc_para', 'age_group', 
                      'contact_matrix_after', 'contact_matrix', 
                      'times_lockdown'), envir = environment())
  
  # system.time({
  #   outcome <- parLapply(cl, 1:nrow(mcmc_para), run.model, times_lockdown = times_lockdown, p = p0)
  #   save(outcome, file = paste0('outcome/', city_list[c], '_0.25_a.RData'))
  # })
  
  ## mass testing and lockdown--------------------------------------------------
  p0[58] <- 2
  clusterExport(cl, c('p0'), envir = environment())
  system.time({
    outcome <- parLapply(cl, 1:nrow(mcmc_para), run.model, times_lockdown = times_lockdown, p = p0)
    save(outcome, file = paste0('outcome/', city_list[c], '_0.25_c.RData'))
  })
  
  ## only mass testing ---------------------------------------------------------
  p0[58] <- 0
  clusterExport(cl, c('p0'), envir = environment())
  system.time({
    outcome <- parLapply(cl, 1:nrow(mcmc_para), run.model, times_lockdown = times_lockdown, p = p0)
    save(outcome, file = paste0('outcome/', city_list[c], '_0.25_b.RData'))
  })
}

stopCluster(cl)


