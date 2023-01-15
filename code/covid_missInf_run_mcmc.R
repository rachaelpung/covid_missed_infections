source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_functions.R')
# sourceCpp('codes/covid_missInf_functions_cpp.cpp')

# start MCMC runs
set.seed(123)

# load obs data
load('input/obs_data.RData')
# obs_data = copy(obs.data)
# names(obs_data) = gsub("\\.", "_", names(obs_data))
# names(obs_data$import_N) = gsub("\\.", "_", names(obs_data$import_N))
# names(obs_data$local_N) = gsub("\\.", "_", names(obs_data$local_N))

# load parameter
param = list()

param$variant = 'W' # D:delta, W:wild

param$time_period = read_excel('data/param.xlsx', sheet = 'time_period')
param$time_period = data.table(param$time_period)
setnames(param$time_period, c('period', 'date_start', 'date_end'))

param$time_period = param$time_period[grep(param$variant, period)]
param$time_period[, period:=gsub(param$variant,'',period)]
param$time_period[, period:=as.numeric(period)]

param$time_period[, doy_start := as.numeric(as.Date(date_start)-as.Date('2019-12-31'))]
param$time_period[, doy_end := as.numeric(as.Date(date_end)-as.Date('2019-12-31'))]
param$time_period[, duration := doy_end-doy_start+1]

param$vaccine_coverage = read_excel('data/param.xlsx', sheet = 'vaccine_coverage')
param$vaccine_coverage = data.table(param$vaccine_coverage)
setnames(param$vaccine_coverage, c('period', 've_inf', 'vac_cov'))

param$vaccine_coverage = param$vaccine_coverage[grep(param$variant, period)]
param$vaccine_coverage[, period:=gsub(param$variant,'',period)]
param$vaccine_coverage[, period:=as.numeric(period)]
  
param$prop_suceptible =  data.table(period = param$vaccine_coverag$period,
                                    prop = (1-param$vaccine_coverage$ve_inf)*param$vaccine_coverage$vac_cov +
                                           (1-param$vaccine_coverage$vac_cov))



# initialise model parameter
param$theta = read_excel('data/param.xlsx', sheet = 'theta')
param$theta = data.table(param$theta)
setnames(param$theta, c('chain', 'variable', 'period', 'value'))
param$theta[, variable_period := paste0(variable, '_period_', period)]

# fixed disease transmission parameters
# generation interval 
param$gen_mean_log = log(5.4) 
param$gen_sd_log = 0.4

# incubation period
param$incub_mean_log = 1.63
param$incub_sd_log = 0.50
# param$incub_mean_log = log(4)
# param$incub_sd_log = sqrt(2*(log(4.4)-log(4)))


param$time_period_matrix = read_excel('data/param.xlsx', sheet = 'time_period')
param$time_period_matrix = data.table(param$time_period_matrix)
setnames(param$time_period_matrix, c('period', 'date_start', 'date_end'))

param$time_period_matrix = param$time_period_matrix[grep(param$variant, period)]
param$time_period_matrix[, period:=gsub(param$variant,'',period)]
param$time_period_matrix[, period:=as.numeric(period)]

param$time_period_matrix[, doy_start := as.numeric(as.Date(date_start)-as.Date('2019-12-31'))]
param$time_period_matrix[, doy_end := as.numeric(as.Date(date_end)-as.Date('2019-12-31'))]
param$time_period_matrix[, duration := doy_end-doy_start+1]


# set-up generation matrix
generation_matrix_imported_N = matrixGenerationImport(param, obs_data, notified = T)
generation_matrix_imported_M = matrixGenerationImport(param, obs_data, notified = F)
generation_matrix_local_N = matrixGenerationLocal(param, obs_data, notified = T)
generation_matrix_local_M = matrixGenerationLocal(param, obs_data, notified = F)

generation_matrix = list(generation_matrix_imported_N = generation_matrix_imported_N, 
                         generation_matrix_imported_M = generation_matrix_imported_M, 
                         generation_matrix_local_N = generation_matrix_local_N, 
                         generation_matrix_local_M = generation_matrix_local_M)

# set-up isolation matrix
isolation_matrix_local_N = matrixIsolationLocal(param, obs_data, notified = T)
isolation_matrix = list(isolation_matrix_local_N = isolation_matrix_local_N)



# select parameter for period of interest
period_interest = 6
param$time_period = param$time_period[period == period_interest]
param$theta = param$theta[period == period_interest]
param$period_interest = period_interest

# load data from previous period
if(period_interest == 1){
  store_previous_period = list(mod_daily_local_N_linked = 0,
                               mod_daily_local_N_unlinked = 0,
                               mod_daily_local_M_linked = 0,
                               mod_daily_local_M_unlinked = 0,
                               mod_daily_local_N_linked_by_doy_isolate = 0,
                               mod_daily_local_N_unlinked_by_doy_isolate = 0)
  
} else{
  list_file = dir('output processed/20220521/wild neg binom/', pattern = paste('store_period_', period_interest-1, '_', param$variant, sep = ''))
  load(paste('output processed/20220521/wild neg binom/', list_file, sep=''))
  store_previous_period = store_period
  names(store_previous_period) = gsub("\\.", "_", names(store_previous_period))
  
  
  
  
  for(i in 1:length(store_previous_period)){
    
    names(store_previous_period[[i]]) = gsub("\\.", "_", names(store_previous_period[[i]]))
    store_previous_period[[i]] = apply(store_previous_period[[i]], 2, median)
  }
  
}


# start MCMC runs
seed = scan('input/setSeed.txt')

# set up clusters
cl = makeCluster(detectCores()) #detectCores()
registerDoParallel(cl)

# CHECK NITER AND PMASK OF THETA[1]

foreach(param_set = 1:8, .packages = c("data.table","mvtnorm"), .combine = "c") %dopar%  {

set.seed(seed[param_set])
  
# unknown model parameter
theta = param$theta[chain == param_set]




# set-up output storage
niter = 150000 # 100000
nparam = nrow(theta) 
store = list()

store$theta = matrix(0, niter, nparam)
colnames(store$theta) = theta$variable_period

store$LL = rep(0,niter)
store$accept = rep(0,niter)
store$doy_start = param$time_period$doy_start  
store$doy_end = param$time_period$doy_end
store$doy = min(param$time_period$doy_start):max(param$time_period$doy_end)
store$duration = length(store$doy)

store$obs_daily_local_N_unlinked = obs_data$daily_local_N_unlinked[store$doy]
store$obs_daily_local_N_linked   = obs_data$daily_local_N_linked[store$doy]

# by time of infection
store$mod_daily_local_N_unlinked = matrix(0, niter, store$duration+13) # modelled notified sec cases with missed infector
store$mod_daily_local_N_linked   = matrix(0, niter, store$duration+13) # modelled notified sec cases with notified infector
store$mod_daily_local_M_unlinked = matrix(0, niter, store$duration+13) # modelled missed sec cases with missed infector by day of infection
store$mod_daily_local_M_linked   = matrix(0, niter, store$duration+13) # modelled missed sec cases with notified infector by day of infection

# by time of isolation
store$mod_daily_local_N_unlinked_by_doy_isolate = matrix(0, niter, store$duration+13) 
store$mod_daily_local_N_linked_by_doy_isolate = matrix(0, niter, store$duration+13) 

pmask = c('prop_import_M_period_3') # fix parameters
pmask = match(pmask,theta$variable_period)
if(is.na(pmask)) pmask=0
# pmask=0

# set-up current parameter
theta_star = theta$value
# theta_star[1] = 0 # for delta


# simulate expected daily incidence given current parameter
output = simulateOutbreak(obs_data, store_previous_period, param, theta_star, generation_matrix, isolation_matrix)

# maybe input obs.data instead of daily.import.N





# calculate log likelihood
store$LL[1] = calculateLogLik(store, output, theta_star, pmask)

# store outputs
store$theta[1,] = theta_star
store$mod_daily_local_N_unlinked[1,] = output$daily_local_N_unlinked
store$mod_daily_local_N_linked[1,]   = output$daily_local_N_linked
store$mod_daily_local_M_unlinked[1,] = output$daily_local_M_unlinked
store$mod_daily_local_M_linked[1,]   = output$daily_local_M_linked

store$mod_daily_local_N_unlinked_by_doy_isolate[1,] = output$daily_local_N_unlinked_by_doy_isolate
store$mod_daily_local_N_linked_by_doy_isolate[1,] = output$daily_local_N_linked_by_doy_isolate


# set-up covariance matrix
npc = rep(1,nparam)


cov_matrix_theta0 = diag(npc)
epsilon0=0.001
accept_rate=0.234

old_output = output

for(iter in 2:niter){ # niter
  
  print(sprintf("Iteration %s of %.f. Acceptance rate %.3f", iter, niter, accept_rate))
  
  old_theta  = store$theta[iter-1,]
  # old_output = output
  old_LL = store$LL[iter-1]
  
  # sample new transmission parameters
  epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^iter)),1e-6) # Stop epsilon getting too big or small
  cov_matrix_theta = epsilon0*cov_matrix_theta0
  

  theta_star = sampleTheta(old_theta,cov_matrix_theta,pmask) #sample nearby global parameter space
  # theta_star = sampleTheta_cpp(old_theta,cov_matrix_theta,pmask)
  
  new_output = simulateOutbreak(obs_data, store_previous_period, param, theta_star, generation_matrix, isolation_matrix)
  # new_output = simulateOutbreak_cpp(obs_data, store_previous_period, param, theta_star, generation_matrix, isolation_matrix)
  
  new_LL = calculateLogLik(store, new_output, theta_star, pmask)
  # new_LL = calculateLogLik_cpp(store, new_output, theta_star, pmask)

  # metropolis hasting step
  mh_output = mh(old_theta, theta_star, old_output, new_output, old_LL, new_LL)
  # mh_output = mh_cpp(old_theta, theta_star, old_output, new_output, old_LL, new_LL)
  
  
  
  # store outputs
  store$theta[iter,] = mh_output$theta
  store$LL[iter] = mh_output$LL
  store$accept[iter] = mh_output$count
  store$mod_daily_local_N_unlinked[iter,] = mh_output$daily_local_N_unlinked
  store$mod_daily_local_N_linked[iter,]   = mh_output$daily_local_N_linked
  store$mod_daily_local_M_unlinked[iter,] = mh_output$daily_local_M_unlinked
  store$mod_daily_local_M_linked[iter,]   = mh_output$daily_local_M_linked
  store$mod_daily_local_N_unlinked_by_doy_isolate[iter,] = mh_output$daily_local_N_unlinked_by_doy_isolate
  store$mod_daily_local_N_linked_by_doy_isolate[iter,]   = mh_output$daily_local_N_linked_by_doy_isolate
  
  old_output = mh_output[4:9]
  
  # start adaptation after 100 iterations
  # update acceptance rate
  if(iter>100){
    accept_rate=sum(store$accept[1:iter])/iter
  }
  
  if(iter %% 10000 == 0){ # 10000
    if(param_set<10){
      save(store,file= paste0("~/missed infections/output raw/20220521/output_variant_", param$variant,"_period_", param$period_interest, '_paramset_0', param_set, "_iter_", iter,".rdata"))
    } else{
      save(store,file= paste0("~/missed infections/output raw/20220521/output_variant_", param$variant,"_period_", param$period_interest, '_paramset_', param_set, "_iter_", iter,".rdata"))
    }
  }
  
}

}


# stop clusters
stopCluster(cl)
