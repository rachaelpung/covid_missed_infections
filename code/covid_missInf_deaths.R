source('codes/v6/covid_missInf_load_library.R')
source('codes/v6/covid_missInf_functions.R')

# load obs data
load('input/obs.data.RData')

# load parameter
param = list()

param$variant = 'W' # D:delta, W:wild

param$time.period = read_excel('data/param.xlsx', sheet = 'time.period.param')
param$time.period = data.table(param$time.period)
setnames(param$time.period, c('period', 'date.start', 'date.end'))

param$time.period = param$time.period[grep(param$variant, period)]
param$time.period[, period:=gsub(param$variant,'',period)]
param$time.period[, period:=as.numeric(period)]

param$time.period[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
param$time.period[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
param$time.period[, duration := doy.end-doy.start+1]


# initialise model parameter
param$theta = read_excel('data/param.xlsx', sheet = 'theta')
param$theta = data.table(param$theta)
setnames(param$theta, c('chain', 'variable', 'period', 'value'))
param$theta[, variable.period := paste0(variable, '.period.', period)]

# fixed disease transmission parameters
# generation interval 
param$gen.mean.log = log(5.4) 
param$gen.sd.log = 0.4

# incubation period
param$incub.mean.log = 1.63
param$incub.sd.log = 0.50
# param$incub.mean.log = log(4)
# param$incub.sd.log = sqrt(2*(log(4.4)-log(4)))


param$time.period.matrix = read_excel('data/param.xlsx', sheet = 'time.period.matrix')
param$time.period.matrix = data.table(param$time.period.matrix)
setnames(param$time.period.matrix, c('period', 'date.start', 'date.end'))

param$time.period.matrix = param$time.period.matrix[grep(param$variant, period)]
param$time.period.matrix[, period:=gsub(param$variant,'',period)]
param$time.period.matrix[, period:=as.numeric(period)]

param$time.period.matrix[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
param$time.period.matrix[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
param$time.period.matrix[, duration := doy.end-doy.start+1]

# set-up isolation matrix
isolation.matrix.local.N = matrixIsolationLocal(param, obs.data, notified = T)
isolation.matrix = list(isolation.matrix.local.N = isolation.matrix.local.N)
isolation.matrix$isolation.matrix.local.N = isolation.matrix$isolation.matrix.local.N[,1:349]

# wild-type icu
icu.wild = isolation.matrix$isolation.matrix.local.N %*% as.matrix(obs.data$daily.local.N.icu[18:366])

# wild-type death
death.wild = isolation.matrix$isolation.matrix.local.N %*% as.matrix(obs.data$daily.local.N.death[18:366])


# load parameter
param = list()

param$variant = 'D' # D:delta, W:wild

param$time.period = read_excel('data/param.xlsx', sheet = 'time.period.param')
param$time.period = data.table(param$time.period)
setnames(param$time.period, c('period', 'date.start', 'date.end'))

param$time.period = param$time.period[grep(param$variant, period)]
param$time.period[, period:=gsub(param$variant,'',period)]
param$time.period[, period:=as.numeric(period)]

param$time.period[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
param$time.period[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
param$time.period[, duration := doy.end-doy.start+1]


# initialise model parameter
param$theta = read_excel('data/param.xlsx', sheet = 'theta')
param$theta = data.table(param$theta)
setnames(param$theta, c('chain', 'variable', 'period', 'value'))
param$theta[, variable.period := paste0(variable, '.period.', period)]

# fixed disease transmission parameters
# generation interval 
param$gen.mean.log = log(5.4) 
param$gen.sd.log = 0.4

# incubation period
param$incub.mean.log = 1.63
param$incub.sd.log = 0.50
# param$incub.mean.log = log(4)
# param$incub.sd.log = sqrt(2*(log(4.4)-log(4)))


param$time.period.matrix = read_excel('data/param.xlsx', sheet = 'time.period.matrix')
param$time.period.matrix = data.table(param$time.period.matrix)
setnames(param$time.period.matrix, c('period', 'date.start', 'date.end'))

param$time.period.matrix = param$time.period.matrix[grep(param$variant, period)]
param$time.period.matrix[, period:=gsub(param$variant,'',period)]
param$time.period.matrix[, period:=as.numeric(period)]

param$time.period.matrix[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
param$time.period.matrix[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
param$time.period.matrix[, duration := doy.end-doy.start+1]

# set-up isolation matrix
isolation.matrix.local.N = matrixIsolationLocal(param, obs.data, notified = T)
isolation.matrix = list(isolation.matrix.local.N = isolation.matrix.local.N)
isolation.matrix$isolation.matrix.local.N = isolation.matrix$isolation.matrix.local.N[,1:142]

# Delta icu
obs.data$daily.local.N.icu = c(obs.data$daily.local.N.icu,0,0)
obs.data$daily.local.N.death = c(obs.data$daily.local.N.death,0,0)

icu.delta = isolation.matrix$isolation.matrix.local.N %*% as.matrix(obs.data$daily.local.N.icu[457:598])

# Delta death
death.delta = isolation.matrix$isolation.matrix.local.N %*% as.matrix(obs.data$daily.local.N.death[457:598])


obs.data$daily.local.N.icu.by.inf.doy = rep(0,598)
obs.data$daily.local.N.death.by.inf.doy = rep(0,598)

obs.data$daily.local.N.icu.by.inf.doy[18:366] = icu.wild
obs.data$daily.local.N.icu.by.inf.doy[457:598] = icu.delta

obs.data$daily.local.N.death.by.inf.doy[18:366] = death.wild
obs.data$daily.local.N.death.by.inf.doy[457:598] = death.delta

