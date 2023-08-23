pred_rec_exp <- function(iter, n, result, sim, cutoff, max_time, site_num, nmax_site, recruit_prd){
  
  pred <- list()
  n_iter <- nrow(result)
  para_pred <- result[(n_iter-iter+1):n_iter,]
  
  for(i in 1:nrow(para_pred)){
    # if(i %% 100 == 0){print(i)}
    print(i)
    site_cur <- length(unique(sim$site))
    alpha <- para_pred[i,1]
    lambda_acc_tol <- para_pred[i,2:(site_cur+1)]
    lambda_c <- para_pred[i,(site_cur+2)]
    lambda_cr <- para_pred[i,(site_cur+3)]
    lambda_e <- para_pred[i,(site_cur+4)]
    psi <- para_pred[i,(site_cur+5)]
    w_all <- para_pred[i,(site_cur+6):(length(unique(sim$subid))+site_cur+5)]
    sim_copy <- sim
    index_pred <- sim %>% group_by(subid) %>% summarise(cen = sum(Y + D)) %>%
      filter(cen == 0) %>% select(subid) %>% unlist()
    n_cut <- length(unique(sim_copy$subid))
    sim_add <- sim[0,]
    if(n_cut < n){
      for(j in 1:site_num){
        t_current <- cutoff
        pat_num <- 0
        if(j <= site_cur){
          lambda_acc <- lambda_acc_tol[j]
        }else{lambda_acc <- rgamma(1, n, site_num * recruit_prd)}
        while(pat_num < nmax_site){
          samp_acc <- rexp(1, lambda_acc)
          t_current <- t_current + samp_acc
          pat_num <- pat_num + 1
          sim_add <- rbind(sim_add, 
                           data.frame(subid = 1, starttime = t_current, 
                                      stoptime = t_current, delta = 0,
                                      time = 0, D = 0, Y = 0, site = j,
                                      grp = rbinom(1,1,0.5)))
        }
      }
      sim_add <- sim_add %>% arrange(starttime) %>% slice(1:(n - n_cut))
      sim_add$subid <- (n_cut+1):n
      sim_copy <- rbind(sim_copy, sim_add)
    }
    
    if(length(unique(sim$subid)) != length(unique(sim_copy$subid))){
      index_pred <- c(index_pred,
                      (length(unique(sim$subid))+1):length(unique(sim_copy$subid)))
    }
    
    for(ind in index_pred){
      temp_pred <- tail(sim_copy[sim_copy$subid == ind,], n = 1)
      sim_copy <- sim_copy[!(sim_copy$subid == ind & sim_copy$delta == 0 & sim_copy$D == 0 & sim_copy$Y == 0),]
      if(ind %in% 1:length(unique(sim$subid))){
        w <- w_all[ind]
        t_current <- cutoff
      } else{
        w <- rgamma(1, psi, psi)
        t_current <- temp_pred$stoptime
      }
      t_cr <- rexp(1, lambda_cr * w^alpha)
      t_cen <- rexp(1, lambda_cen)
      stoptime_c <- t_cen + t_current
      stoptime_cr <- t_cr + t_current
      starttime <- temp_pred$starttime
      t_cp <- t_cr < t_cen
      stoptime_min <- min(stoptime_cr, stoptime_c)
      
      while(t_current < max_time){
        t_event <- rexp(1, lambda_e * w)
        stoptime_e <- t_event + t_current
        if(stoptime_e >= max_time){break}
        if(stoptime_e < stoptime_min){
          temp_pred$starttime <- starttime
          temp_pred$stoptime <- stoptime_e
          temp_pred$time <- temp_pred$stoptime - temp_pred$starttime
          temp_pred$delta <- 1
          temp_pred$D <- 0
          temp_pred$Y <- 0
          t_current <- temp_pred$stoptime
          starttime <- t_current
        } else if(stoptime_min < max_time){
          temp_pred$starttime <- starttime
          temp_pred$stoptime <- stoptime_min
          temp_pred$time <- temp_pred$stoptime - temp_pred$starttime
          temp_pred$delta <- 0
          temp_pred$D <- as.numeric(t_cp)
          temp_pred$Y <- 1 - t_cp
        }
        sim_copy <- rbind(sim_copy, temp_pred)
        if((temp_pred$Y + temp_pred$D) == 1){break}
      }
      if((temp_pred$Y + temp_pred$D) == 0){
        temp_pred$starttime <- starttime
        temp_pred$stoptime <- max_time
        temp_pred$time <- temp_pred$stoptime - temp_pred$starttime
        temp_pred$delta <- 0
        temp_pred$D <- 0
        temp_pred$Y <- 0
        sim_copy <- rbind(sim_copy, temp_pred)
      }
    }
    pred[[i]] <- sim_copy %>% arrange(subid)
  }
  return(pred)
}