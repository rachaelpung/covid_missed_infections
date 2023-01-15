#' create matrix of generation intervals for imported cases by notification status
#' @param param mixture of fixed and unknown model parameters
#' @param obs_data observed imported and local cases
#' @param notified logical input parameter - T if notified, F if missed

matrixGenerationImport <- function(param, obs_data, notified = T){
  
  time_period = param$time_period
  t_max = sum(time_period$duration)
  
  gen_mean_log = param$gen_mean_log
  gen_sd_log = param$gen_sd_log
  
  incub_mean_log = param$incub_mean_log
  incub_sd_log = param$incub_sd_log
  
  # set-up generation matrix 
  generation_matrix =  matrix(0,nrow=t_max,ncol=(t_max+13))
  
  import_N = obs_data$import_N
  import_N = import_N[!is.na(doy_onset)]

  for(p in time_period$period){
    
    doy_start = time_period$doy_start[p]
    doy_end = time_period$doy_end[p]
    subset = import_N[doy_arrival>=doy_start & doy_arrival<=doy_end & is.finite(doy_isolate_from_comm),]
    
    # incubation period distribution
    t = 1:14
    dist_incub = dlnorm(t, meanlog = incub_mean_log, sdlog = incub_sd_log, log = FALSE)
    dist_incub = data.table(day = t, pdf  = dist_incub/sum(dist_incub))
    
    # generation interval distribution
    t = 1:14
    dist_gen = dlnorm(t, meanlog = gen_mean_log, sdlog = gen_sd_log, log = FALSE)
    dist_gen = data.table(day = t, pdf  = dist_gen/sum(dist_gen))
    # or to deconvolute from the serial interval and incubation period
    
    # time from infection to arrival distribution
    subset = subset[rep(seq_len(nrow(subset)), each=14)]
    subset[, doy_infection := doy_onset-rep(1:14, times=.N/14)]
    subset[, duration_incubate := doy_onset-doy_infection] # can only take values more than 1
    
    subset = subset[doy_infection<doy_arrival]  # infection after arrival assumed cannot happen
    subset[dist_incub, prob_infected_days_ago := i.pdf, on=c(duration_incubate='day') ]
    subset[, prob_infected_days_ago:=prob_infected_days_ago/sum(prob_infected_days_ago), by=.(case)]
    
    subset[, duration_infection_to_arrival := doy_arrival - doy_infection]
    dist_infection_to_arrival = subset[,.(sum(prob_infected_days_ago)), by=.(duration_infection_to_arrival)]
    dist_infection_to_arrival = sample(dist_infection_to_arrival$duration_infection_to_arrival, 1000,
                                       prob = dist_infection_to_arrival$ V1,
                                       replace = T)
    
    fit = fitdist(dist_infection_to_arrival, 'exp')
    
    # hist(dist_infection_to_arrival,prob=TRUE,breaks = seq(0,30,1))
    # plot(function(x) dexp(x, rate = fit$estimate[1]), from=0, to=30, col="red", add=TRUE)
    
    t=1:14
    dist_infection_to_arrival = dexp(t, rate = fit$estimate[1], log = FALSE)
    dist_infection_to_arrival = data.table(day = t, pdf_infection_to_arrival  = dist_infection_to_arrival/sum(dist_infection_to_arrival))
    
    
    # work out the generation interval available for infection after left truncation
    dist_gen_imported_missed = dist_infection_to_arrival
    dist_gen_imported_missed[, generation_interval_start := day]
    dist_gen_imported_missed[, generation_interval_end := 14]
    
    dist_gen_imported_missed[, days_at_large := generation_interval_end - generation_interval_start + 1]
    dist_gen_imported_missed=dist_gen_imported_missed[rep(seq_len(nrow(dist_gen_imported_missed)), times = days_at_large)]
    dist_gen_imported_missed[,generation_interval := seq_len(.N), by=.(day)]
    dist_gen_imported_missed[,generation_interval := generation_interval+generation_interval_start-1, by=.(day)]
    dist_gen_imported_missed[dist_gen, pdf_generation_interval := i.pdf, on=c(generation_interval='day')]
    
    dist_gen_imported_missed[, pdf := pdf_infection_to_arrival * pdf_generation_interval]
    dist_gen_imported_missed = dist_gen_imported_missed[,.(sum(pdf)), by=.(generation_interval)]
    setnames(dist_gen_imported_missed, c('day', 'pdf'))
    
    
    # time from infection to isolation
    subset[, duration_infection_to_isolate := doy_isolate_from_comm - doy_infection]
    dist_infection_to_isolate = subset[,.(sum(prob_infected_days_ago)), by=.(duration_infection_to_isolate)]
    dist_infection_to_isolate = sample(dist_infection_to_isolate$duration_infection_to_isolate, 1000, 
                                       prob = dist_infection_to_isolate$V1,
                                       replace = T)
    
    # fit = fitdist(dist_infection_to_isolate, 'lnorm')
    fit = fitdist(dist_infection_to_isolate[dist_infection_to_isolate!=0], 'lnorm')
    
    # hist(dist_infection_to_isolate,prob=TRUE,breaks = c(0:50))
    # plot(function(x) dlnorm(x, meanlog = fit$estimate[1], sdlog =  fit$estimate[2]), from=0, to=50, col="red", add=TRUE)
    
    t=1:14
    dist_infection_to_isolate = dlnorm(t, meanlog = fit$estimate[1], sdlog =  fit$estimate[2], log = FALSE)
    dist_infection_to_isolate = data.table(day = t, pdf_infection_to_isolation  = dist_infection_to_isolate/sum(dist_infection_to_isolate))
    dist_infection_to_isolate[, cdf_infection_to_isolation := cumsum(pdf_infection_to_isolation)]
    
    dist_gen_imported_notified = dist_gen_imported_missed
    dist_gen_imported_notified$pdf = dist_gen_imported_notified$pdf * (1-dist_infection_to_isolate$cdf_infection_to_isolation)
    
    # plot(dist_gen$day, dist_gen$pdf, type='l', col = 'red')
    # lines(dist_gen_imported_missed$day, dist_gen_imported_missed$pdf, col = 'orange')
    # lines(dist_gen_imported_notified$day, dist_gen_imported_notified$pdf, col = 'black')
    # 
    index_start = sum(time_period$duration[1:p]) - time_period$duration[p] + 1
    index_end = sum(time_period$duration[1:p])
    
    for(t in index_start:index_end){
      if(notified == T){
        generation_matrix[t,t:(t+13)] = dist_gen_imported_notified$pdf # from day of arrival
      }else{
        generation_matrix[t,t:(t+13)] = dist_gen_imported_missed$pdf
      }

    }
    
  }

  
  return(generation_matrix)
  
}


#' create matrix of generation intervals for local cases by notification status
#' @param param mixture of fixed and unknown model parameters
#' @param obs_data observed imported and local cases
#' @param notified logical input parameter - T if notified, F if missed

matrixGenerationLocal <- function(param, obs_data, notified = T){
  
  time_period = param$time_period
  t_max = sum(time_period$duration)
  
  gen_mean_log = param$gen_mean_log
  gen_sd_log = param$gen_sd_log
  
  incub_mean_log = param$incub_mean_log
  incub_sd_log = param$incub_sd_log
  
  # set-up generation matrix 
  generation_matrix =  matrix(0,nrow=t_max,ncol=(t_max+13))
  
  local_N = obs_data$local_N
  local_N = local_N[!is.na(doy_onset)]
  
  for(p in time_period$period){
    
    doy_start = time_period$doy_start[p]
    doy_end = time_period$doy_end[p]
    subset = local_N[doy_isolate_from_comm>=doy_start & doy_isolate_from_comm<=doy_end,]
  
    # incubation period distribution
    t = 1:14
    dist_incub = dlnorm(t, meanlog = incub_mean_log, sdlog = incub_sd_log, log = FALSE)
    dist_incub = data.table(day = t, pdf  = dist_incub/sum(dist_incub))
    
    # generation interval distribution
    t = 1:14
    dist_gen = dlnorm(t, meanlog = gen_mean_log, sdlog = gen_sd_log, log = FALSE)
    dist_gen = data.table(day = t, pdf  = dist_gen/sum(dist_gen))
    # or to deconvolute from the serial interval and incubation period
    
    # time from infection to arrival distribution
    subset = subset[rep(seq_len(nrow(subset)), each=14)]
    subset[, doy_infection := doy_onset-rep(1:14, times=.N/14)]
    subset[, duration_incubate := doy_onset-doy_infection] # can only take values more than 1
    
    subset[dist_incub, prob_infected_days_ago := i.pdf, on=c(duration_incubate='day') ]
    subset[, prob_infected_days_ago:=prob_infected_days_ago/sum(prob_infected_days_ago), by=.(case)]
    
    subset[, duration_infection_to_isolate := doy_isolate_from_comm - doy_infection]
    dist_infection_to_isolate = subset[,.(sum(prob_infected_days_ago)), by=.(duration_infection_to_isolate)]
    dist_infection_to_isolate = sample(dist_infection_to_isolate$duration_infection_to_isolate, 1000,
                                       prob = dist_infection_to_isolate$V1,
                                       replace = T)
    
    fit = selm(dist_infection_to_isolate ~ 1) # selm is "skew-elliptic lm"
    # summary(fit)
    # hist(dist_infection_to_isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dsn(x, dp=fit@param$dp), from=-25, to=50, col="red", add=TRUE)
    # 
    # fit = fitdist(dist_infection_to_isolate, 'norm')
    # hist(dist_infection_to_isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dnorm(x, mean = fit$estimate[1], sd =  fit$estimate[2]), from=-25, to=50, col="red", add=TRUE)
    
    t=1:14
    dist_infection_to_isolate = dsn(t, dp=fit@param$dp, log = FALSE)
    # dist_infection_to_isolate = dlnorm(t, meanlog = fit$estimate[1], sdlog =  fit$estimate[2], log = FALSE)
    dist_infection_to_isolate = data.table(day = t, pdf_infection_to_isolation  = dist_infection_to_isolate/sum(dist_infection_to_isolate))
    dist_infection_to_isolate[, cdf_infection_to_isolation := cumsum(pdf_infection_to_isolation)]
    
    dist_gen_local_notified = dist_gen_local_missed = dist_gen
    dist_gen_local_notified$pdf = dist_gen_local_notified$pdf * (1-dist_infection_to_isolate$cdf_infection_to_isolation)
    
    
    # plot(dist_gen_local_missed$day, dist_gen_local_missed$pdf, type='l', col = 'red')
    # lines(dist_gen_local_notified$day, dist_gen_local_notified$pdf, col = 'black')
    # 
    index_start = sum(time_period$duration[1:p]) - time_period$duration[p] + 1
    index_end = sum(time_period$duration[1:p])
    
    for(t in index_start:index_end){
      if(notified == T){
        generation_matrix[t,t:(t+13)] = dist_gen_local_notified$pdf # from day of arrival
      }else{
        generation_matrix[t,t:(t+13)] = dist_gen_local_missed$pdf
      }

    }
    
  }  
  
  return(generation_matrix)
  
}


#' create matrix of infection to isolation intervals for local cases by notification status
#' @param param mixture of fixed and unknown model parameters
#' @param obs_data observed imported and local cases
#' @param notified logical input parameter - T if notified, F if missed

matrixIsolationLocal <- function(param, obs_data, notified = T){
  
  time_period = param$time_period
  t_max = sum(time_period$duration)
  
  gen_mean_log = param$gen_mean_log
  gen_sd_log = param$gen_sd_log
  
  incub_mean_log = param$incub_mean_log
  incub_sd_log = param$incub_sd_log
  
  # set-up generation matrix 
  isolation_matrix =  matrix(0,nrow=t_max,ncol=(t_max+13))
  
  local_N = obs_data$local_N
  local_N = local_N[!is.na(doy_onset)]
  
  
  for(p in time_period$period){
    
    doy_start = time_period$doy_start[p]
    doy_end = time_period$doy_end[p]
    subset = local_N[doy_isolate_from_comm>=doy_start & doy_isolate_from_comm<=doy_end,]
    
    # incubation period distribution
    t = 1:14
    dist_incub = dlnorm(t, meanlog = incub_mean_log, sdlog = incub_sd_log, log = FALSE)
    dist_incub = data.table(day = t, pdf  = dist_incub/sum(dist_incub))
    
    
    # generation interval distribution
    t = 1:14
    dist_gen = dlnorm(t, meanlog = gen_mean_log, sdlog = gen_sd_log, log = FALSE)
    dist_gen = data.table(day = t, pdf  = dist_gen/sum(dist_gen))
    # or to deconvolute from the serial interval and incubation period
    
    # time from infection to arrival distribution
    subset = subset[rep(seq_len(nrow(subset)), each=14)]
    subset[, doy_infection := doy_onset-rep(1:14, times=.N/14)]
    subset[, duration_incubate := doy_onset-doy_infection] # can only take values more than 1
    
    subset[dist_incub, prob_infected_days_ago := i.pdf, on=c(duration_incubate='day') ]
    subset[, prob_infected_days_ago:=prob_infected_days_ago/sum(prob_infected_days_ago), by=.(case)]
    
    subset[, duration_infection_to_isolate := doy_isolate_from_comm - doy_infection]
    dist_infection_to_isolate = subset[,.(sum(prob_infected_days_ago)), by=.(duration_infection_to_isolate)]
    dist_infection_to_isolate = sample(dist_infection_to_isolate$duration_infection_to_isolate, 1000,
                                       prob = dist_infection_to_isolate$ V1,
                                       replace = T)
    
    fit = selm(dist_infection_to_isolate ~ 1) # selm is "skew-elliptic lm"
    # summary(fit)
    # hist(dist_infection_to_isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dsn(x, dp=fit@param$dp), from=-25, to=50, col="red", add=TRUE)
    # 
    # fit = fitdist(dist_infection_to_isolate, 'norm')
    # hist(dist_infection_to_isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dnorm(x, mean = fit$estimate[1], sd =  fit$estimate[2]), from=-25, to=50, col="red", add=TRUE)
    
    t=1:14
    dist_infection_to_isolate = dsn(t, dp=fit@param$dp, log = FALSE)
    # dist_infection_to_isolate = dlnorm(t, meanlog = fit$estimate[1], sdlog =  fit$estimate[2], log = FALSE)
    dist_infection_to_isolate = data.table(day = t, pdf_infection_to_isolation  = dist_infection_to_isolate/sum(dist_infection_to_isolate))
    dist_infection_to_isolate[, cdf_infection_to_isolation := cumsum(pdf_infection_to_isolation)]
    
    dist_gen_local_notified = dist_gen_local_missed = dist_gen
    dist_gen_local_notified$pdf = dist_gen_local_notified$pdf * (1-dist_infection_to_isolate$cdf_infection_to_isolation)
    
  
    index_start = sum(time_period$duration[1:p]) - time_period$duration[p] + 1
    index_end = sum(time_period$duration[1:p])
    
    for(t in index_start:index_end){
      if(notified == T){
        isolation_matrix[t,t:(t+13)] = dist_infection_to_isolate$pdf # from day of arrival
      }else{
        isolation_matrix[t,t:(t+13)] = dist_infection_to_isolate$pdf
      }

    }

    
  }  
  
  return(isolation_matrix)
  
}


#' sample unknown parameters from covariance matrix
#' @param inputTheta current set of model parameters
#' @param covarTheta covariance matrix
#' @param pmask mask parameters that are fixed

sampleTheta <- function(inputTheta, covarTheta, pmask){
  
  mean_vector_theta = inputTheta
  mean_vector_theta[c(3,4)] = log(mean_vector_theta[c(3,4)]/(1-mean_vector_theta[c(3,4)]))
  mean_vector_theta[c(1,2,5)] = log(mean_vector_theta[c(1,2,5)])
 
  if(pmask!=0) diag(covarTheta)[pmask] = 1e-12

  theta_star =  as.numeric(mvtnorm::rmvnorm(1, mean_vector_theta, covarTheta))
  theta_star[c(3,4)] = exp(theta_star[c(3,4)])/(1+exp(theta_star[c(3,4)]))
  theta_star[c(1,2,5)] = exp(theta_star[c(1,2,5)])
  
  if(pmask!=0) theta_star[pmask] = inputTheta[pmask]
  
  
  
  return(theta_star)
}


#' simulate outbreak
#' @param obs_data observed data
#' @param store_previous_period data from previous time period
#' @param theta_star model parameters
#' @param generation_matrix 
#' @param isolation_matrix

simulateOutbreak <- function(obs_data, store_previous_period=NULL, param, theta_star, generation_matrix, isolation_matrix){
  
  period_interest = param$period_interest
  if(period_interest != 1){
    extract = (length(store_previous_period$mod_daily_local_N_linked)-12) : length(store_previous_period$mod_daily_local_N_linked)
  } else{
    extract = 1
  }
  
  
  # format parameter dataframe
  d=length(theta_star)/5
  # d=length(theta_star)/3
  
  time_period = param$time_period
  t_max = sum(time_period$duration)
  theta = data.table(period = time_period$period)
  
  theta[,prop_import_M:= theta_star[1:d]]
  theta[,R:= theta_star[(1+d):(2*d)]]
  theta[,eff_contact_tracing:= theta_star[(1+2*d):(3*d)]]
  theta[,eff_case_finding:= theta_star[(1+3*d):(4*d)]]  
  
  prop_suceptible = param$prop_suceptible[period==period_interest]$prop
  
  # set-up generation matrix
  # generation_matrix dimension of 349 by 369 first index correspondence to doy 18 last index doy 366
  # extract the relevant portion of the generation matrix for time period of interest
  
  # wild
  row = (min(time_period$doy_start) - 17):(max(time_period$doy_end) - 17)
  col = (min(time_period$doy_start) - 17):((max(time_period$doy_end) - 17) + 13 )
  
  # delta
  # row = (min(time_period$doy_start) - 456):(max(time_period$doy_end) - 456)
  # col = (min(time_period$doy_start) - 456):((max(time_period$doy_end) - 456) + 13 )

  
  generation_matrix_imported_N = generation_matrix$generation_matrix_imported_N[row,col]
  generation_matrix_imported_M = generation_matrix$generation_matrix_imported_M[row,col]
  generation_matrix_local_N = generation_matrix$generation_matrix_local_N[row,col]
  generation_matrix_local_M = generation_matrix$generation_matrix_local_M[row,col]
  
  # generate missed imported cases
  # obs_data$daily_arrival_import_N since 1 Jan 2020
  daily_arrival_import_N = obs_data$daily_arrival_import_N_ma[min(time_period$doy_start):max(time_period$doy_end)]
  daily_arrival_import_M = daily_arrival_import_N * rep(theta$prop_import_M, times=time_period$duration)
  
  # calculate the spread of secondary cases arising from imported cases
  sec_cases_infector_N = generation_matrix_imported_N * daily_arrival_import_N * rep(theta$R, times=time_period$duration) * rep(prop_suceptible, times=time_period$duration) 
  sec_cases_infector_M = generation_matrix_imported_M * daily_arrival_import_M * rep(theta$R, times=time_period$duration) * rep(prop_suceptible, times=time_period$duration) 
  
  # daily local notified and missed secondary cases from notified imported infector
  sec_cases_N_infector_N = colSums(sec_cases_infector_N) * rep(theta$eff_contact_tracing, times=time_period$duration+13)
  sec_cases_M_infector_N = colSums(sec_cases_infector_N) * (1-rep(theta$eff_contact_tracing, times=time_period$duration+13))
  # sec_cases_N_infector_N == linked case
  
  # daily local notified and missed secondary cases from missed imported infector
  sec_cases_N_infector_M = colSums(sec_cases_infector_M) * rep(theta$eff_case_finding, times=time_period$duration+13)
  sec_cases_M_infector_M = colSums(sec_cases_infector_M) * (1-rep(theta$eff_case_finding, times=time_period$duration+13))
  # sec_cases_N_infector_M == unlinked case
  
  # add import-related notified secondary cases (i.e. first generation of local notified cases)
  # add import-related missed secondary cases (i.e. first generation of local missed cases)
  daily_local_N = sec_cases_N_infector_N + sec_cases_N_infector_M + 
                  c(store_previous_period$mod_daily_local_N_linked[extract], rep(0,length(sec_cases_N_infector_N)-length(extract))) +
                  c(store_previous_period$mod_daily_local_N_unlinked[extract], rep(0,length(sec_cases_N_infector_N)-length(extract)))
    
  daily_local_M = sec_cases_M_infector_N + sec_cases_M_infector_M + 
                  c(store_previous_period$mod_daily_local_M_linked[extract], rep(0,length(sec_cases_M_infector_N)-length(extract))) +
                  c(store_previous_period$mod_daily_local_M_unlinked[extract], rep(0,length(sec_cases_M_infector_N)-length(extract)))
  
  daily_local_N_linked = sec_cases_N_infector_N + c(store_previous_period$mod_daily_local_N_linked[extract], rep(0,length(sec_cases_N_infector_N)-length(extract))) 
  daily_local_N_unlinked = sec_cases_N_infector_M + c(store_previous_period$mod_daily_local_N_unlinked[extract], rep(0,length(sec_cases_N_infector_M)-length(extract))) 
  
  daily_local_M_linked = sec_cases_M_infector_N + c(store_previous_period$mod_daily_local_M_linked[extract], rep(0,length(sec_cases_M_infector_N)-length(extract))) 
  daily_local_M_unlinked = sec_cases_M_infector_M + c(store_previous_period$mod_daily_local_M_unlinked[extract], rep(0,length(sec_cases_M_infector_M)-length(extract))) 
  
  # Add onwards transmission from subsequent generations 
  for(t in 1:t_max){
    
    sec_cases_infector_N = generation_matrix_local_N[t,] * daily_local_N[t] * rep(theta$R, times=time_period$duration+13) * rep(prop_suceptible, times=time_period$duration+13) 
    sec_cases_infector_M = generation_matrix_local_M[t,] * daily_local_M[t] * rep(theta$R, times=time_period$duration+13) * rep(prop_suceptible, times=time_period$duration+13) 

    sec_cases_N_infector_N = sec_cases_infector_N * rep(theta$eff_contact_tracing, times=time_period$duration+13)
    sec_cases_M_infector_N = sec_cases_infector_N * (1-rep(theta$eff_contact_tracing, times=time_period$duration+13))

    sec_cases_N_infector_M = sec_cases_infector_M * rep(theta$eff_case_finding, times=time_period$duration+13)
    sec_cases_M_infector_M = sec_cases_infector_M * (1-rep(theta$eff_case_finding, times=time_period$duration+13))
    
    daily_local_N = daily_local_N + sec_cases_N_infector_N + sec_cases_N_infector_M
    daily_local_M = daily_local_M + sec_cases_M_infector_N + sec_cases_M_infector_M
    
    daily_local_N_linked = daily_local_N_linked + sec_cases_N_infector_N
    daily_local_N_unlinked = daily_local_N_unlinked + sec_cases_N_infector_M
    
    daily_local_M_linked = daily_local_M_linked + sec_cases_M_infector_N
    daily_local_M_unlinked = daily_local_M_unlinked + sec_cases_M_infector_M
    
  }
  
  # notified cases are based on infection dates 
  # need to convert to date of isolation
  isolation_matrix_local_N = isolation_matrix$isolation_matrix_local_N[row,col]
  
  daily_local_N_linked_by_doy_isolate = isolation_matrix_local_N * daily_local_N_linked[1:length(row)]
  daily_local_N_linked_by_doy_isolate = colSums(daily_local_N_linked_by_doy_isolate) +
                                        c(store_previous_period$mod_daily_local_N_linked_by_doy_isolate[extract], rep(0,length(daily_local_N)-length(extract)))
  
  
  daily_local_N_unlinked_by_doy_isolate = isolation_matrix_local_N * daily_local_N_unlinked[1:length(row)] 
  daily_local_N_unlinked_by_doy_isolate = colSums(daily_local_N_unlinked_by_doy_isolate) +
                                          c(store_previous_period$mod_daily_local_N_unlinked_by_doy_isolate[extract], rep(0,length(daily_local_N)-length(extract)))
  
  list(daily_local_N_unlinked = daily_local_N_unlinked, 
       daily_local_N_linked = daily_local_N_linked,
       daily_local_M_unlinked = daily_local_M_unlinked,
       daily_local_M_linked = daily_local_M_linked,
       daily_local_N_unlinked_by_doy_isolate = daily_local_N_unlinked_by_doy_isolate,
       daily_local_N_linked_by_doy_isolate = daily_local_N_linked_by_doy_isolate)
  
}


#' calculate log likelihood
#' @param obs_data observed data
#' @param output model output

calculateLogLik <- function(store, output, theta_star, pmask){
  
  duration = store$duration
  check_length = length(store$obs_daily_local_N_linked) != length(output$daily_local_N_linked[1:duration]) 
  if(check_length){print('modelled output has different length from observed data')}
  
  # LL = 0
  LL = calculateLogPrior(theta_star, pmask)
  
  LL = LL + sum(dnbinom(store$obs_daily_local_N_unlinked,mu=output$daily_local_N_unlinked_by_doy_isolate[1:duration],size=1/theta_star[5],log=T),
           dnbinom(store$obs_daily_local_N_linked,mu=output$daily_local_N_linked_by_doy_isolate[1:duration],size=1/theta_star[5],log=T))
  
  # LL = sum(dpois((store$obs_daily_local_N_unlinked+
  #                   store$obs_daily_local_N_linked),
  #                lambda=(output$daily_local_N_unlinked_by_doy_isolate[1:duration]+
  #                          output$daily_local_N_linked_by_doy_isolate[1:duration]),log=T))
  
  # LL = LL + sum(dnbinom((store$obs_daily_local_N_unlinked+store$obs_daily_local_N_linked),
  #                mu=(output$daily_local_N_unlinked_by_doy_isolate[1:duration]+
  #                    output$daily_local_N_linked_by_doy_isolate[1:duration]),size=1/theta_star[5],log=T))

  return(LL)
  
}


#' calculate log prior
#' @param obs_data observed data
#' @param output model output
#' 
calculateLogPrior <- function(theta_star, pmask){
  
  # uninformative prior
  # priorProp = ifelse(theta_star[1] > 50 | theta_star[1] < 0, -9999, 0)
  # priorR = ifelse(theta_star[2] > 50 | theta_star[2] < 0, -9999, 0)
  
  # informative prior
  priorProp = dlnorm(theta_star[1], mean = 1, sd=1, log = T)
  priorR = dlnorm(theta_star[2], mean = 1, sd=0.5, log = T)

  # for wild linked unlinked informative prior
  priorEffct = dbeta(theta_star[3], shape1 = 0.1+3, shape2 = 0.1+6-3, log=T)
  priorEffcf = dbeta(theta_star[4], shape1 = 0.1+3, shape2 = 0.1+6-3, log=T)
  
  # for wild and delta total
  # priorEffct = ifelse(theta_star[3] > 0.99 | theta_star[3] <0.01, -9999, 0)
  # priorEffcf = ifelse(theta_star[4] > 0.99 | theta_star[4] <0.01, -9999, 0)
  
  
  prior = c(priorProp, priorR, priorEffct, priorEffcf)
  prior[pmask] = 0
  prior = sum(prior)
  
  return(prior)
  
}



#' metropolis hastings step
#' @param old_theta old model parameter
#' @param theta_star new model parameter 
#' @param old_output old model output 
#' @param new_output new model output
#' @param old_LL old likelihood
#' @param new_LL new likelihood

mh <- function(old_theta, theta_star, old_output, new_output, old_LL, new_LL){
  
  reject = FALSE
  count	= 0
  
  # if(!reject){
    
    differenceLL = new_LL-old_LL 
    
    if(log(runif(1))>differenceLL){
      reject=TRUE
    } else{
      count	= 1
    }
  # }
  
  
  if(reject) return(c(list(theta = old_theta,
                           count = count,
                           LL = old_LL),
                           old_output))
  
  return(c(list(theta = theta_star,
                count = count,
                LL = new_LL),
                new_output))
  
}
