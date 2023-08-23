func_analysis <- function(n, sim, sim_c, n_exp, recruit_prd){
  
  acc <- sim[!duplicated(sim$subid),]
  acc <- acc %>% group_by(site) %>% arrange(starttime) %>% 
    mutate(diff_acc = starttime - lag(starttime))
  acc$diff_acc <- ifelse(is.na(acc$diff_acc),
                         acc$starttime, acc$diff_acc)
  
  ns <- acc %>% group_by(site) %>% summarise(n = n())
  site_cur <- nrow(ns)
  starttime <- matrix(NA, nrow = site_cur, ncol = max(ns$n))
  for(i in 1:site_cur){
    for(j in 1:ns$n[i]){
      temp <- acc[acc$site == i,]
      starttime[i,j] <- temp$diff_acc[j]
    }
  }
  
  J <- sim %>% group_by(subid) %>% tally()
  m_events <- max(J$n)
  time <- delta <- matrix(NA, nrow = n, ncol = m_events)
  for(i in 1:n){
    for(j in 1:J$n[i]){
      unique_id <- unique(sim$subid)
      temp <- sim[sim$subid == unique_id[i],]
      time[i,j] <- temp$time[j] 
      delta[i,j] <- temp$delta[j]
    }
  }
  zeros_e <- matrix(0, nrow = n, ncol = m_events)
  
  d.jags <- list(n = n, J = J$n, D = sim_c$D, Y = sim_c$Y, time = time, delta = delta,
                 time_cr = sim_c$time, time_c = sim_c$time, zeros_e = zeros_e, 
                 zeros_cr = rep(0, n), zeros_c = rep(0, n), starttime = starttime, 
                 s = site_cur, ns = ns$n,
                 n_exp = n_exp, t_exp = recruit_prd)
  
  i.jags <- function() {
    list(lambda_e = runif(1, 0.1), lambda_cr = runif(1, 0.1),
         lambda_c = runif(1, 0.1), psi = runif(1),
         alpha = rnorm(1))}
  
  p.jags <- c("lambda_e", "lambda_cr", "lambda_c", "lambda_acc", "psi", "w", "alpha")
  
  m <- jags.model(data = d.jags, file = "PH_case1.txt", inits = i.jags, n.chains = 1)
  update(m, 1000) 
  res <- coda.samples(m, variable.names = p.jags, n.iter = 20000, thin = 1)
  
  result <- as.mcmc(do.call(rbind, res))
  return(result)
}