source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_functions.R')

load('input/obs_data.RData')
load('input/time_period.RData')

# load parameter
param = list()

param$variant = 'W' # D:delta, W:wild

param$time_period = time_period
param$time_period = param$time_period[variant==param$variant]

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
# param$incub_mean_log = 1.63
# param$incub_sd_log = 0.50
param$incub_mean_log = log(4)
param$incub_sd_log = sqrt(2*(log(4.4)-log(4)))

# set-up isolation matrix
isolation_matrix_local_N = matrixIsolationLocal(param, obs_data, notified = T)
isolation_matrix = list(isolation_matrix_local_N = isolation_matrix_local_N)
isolation_matrix$isolation_matrix_local_N = isolation_matrix$isolation_matrix_local_N[,1:349]

# wild-type icu
icu_wild = isolation_matrix$isolation_matrix_local_N %*% as.matrix(obs_data$daily_local_N_icu[18:366])

# wild-type death
death_wild = isolation_matrix$isolation_matrix_local_N %*% as.matrix(obs_data$daily_local_N_death[18:366])


# load parameter
param = list()

param$variant = 'D' # D:delta, W:wild

param$time_period = time_period
param$time_period = param$time_period[variant==param$variant]

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
# param$incub_mean_log = 1.63
# param$incub_sd_log = 0.50
param$incub_mean_log = log(4)
param$incub_sd_log = sqrt(2*(log(4.4)-log(4)))

# set-up isolation matrix
isolation_matrix_local_N = matrixIsolationLocal(param, obs_data, notified = T)
isolation_matrix = list(isolation_matrix_local_N = isolation_matrix_local_N)
isolation_matrix$isolation_matrix_local_N = isolation_matrix$isolation_matrix_local_N[,1:142]

# Delta icu
obs_data$daily_local_N_icu = c(obs_data$daily_local_N_icu,0,0)
obs_data$daily_local_N_death = c(obs_data$daily_local_N_death,0,0)

icu_delta = isolation_matrix$isolation_matrix_local_N %*% as.matrix(obs_data$daily_local_N_icu[457:598])

# Delta death
death_delta = isolation_matrix$isolation_matrix_local_N %*% as.matrix(obs_data$daily_local_N_death[457:598])


obs_data$daily_local_N_icu_by_inf_doy = rep(0,598)
obs_data$daily_local_N_death_by_inf_doy = rep(0,598)

obs_data$daily_local_N_icu_by_inf_doy[18:366] = icu_wild
obs_data$daily_local_N_icu_by_inf_doy[457:596] = icu_delta

obs_data$daily_local_N_death_by_inf_doy[18:366] = death_wild
obs_data$daily_local_N_death_by_inf_doy[457:596] = death_delta

