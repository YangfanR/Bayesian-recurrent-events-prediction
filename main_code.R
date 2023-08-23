library(rjags)
library(ggplot2)
library(dplyr)
library(snowfall)
library(tidyverse)
library(reshape2)

source("func_analysis.R")
source("pred_rec_exp.R")

#########  load data ######### 
load("data.RData")
sim <- sim_tol[[1]] # dataset for analysis and prediction
realtime <- sim_tol[[2]] # real time of reaching the desired number of recurrent events
sim_comp <- sim_tol[[3]] # complete dataset without cut-off
cutoff <- sim_tol[[4]] # cut-off time
lambda_acc <- sim_tol[[5]] # real lambda of different sites

#########  preprocessing data ######### 
unique_id <- unique(sim$subid)
for(i in 1:length(unique_id)){
  id <- unique_id[i]
  sim[sim$subid == id, "subid"] <- i
}
sim_st <- sim[!duplicated(sim$subid),"starttime"]
sim_c <- sim %>% group_by(subid) %>% slice(n())
sim_c$starttime <- sim_st
sim_c$time <- sim_c$stoptime - sim_c$starttime

######### analysis #########

result <- func_analysis(nrow(sim_c), sim, sim_c, n, recruit_prd)
summary(result[,1:(site_num + 5)]) # check estimates

######### predictions #########

max_time <- realtime + 2
pred_save <- pred_rec_exp(iter, n, result, sim, cutoff, max_time, site_num, nmax_site, recruit_prd)
