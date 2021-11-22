source('codes/v6/covid_missInf_load_library.R')

# simulate observed daily linked and unlinked cases
param = list()
param$time.period = data.table(period = c(1,2),
                               date.start = as.Date(c('2020-01-01','2020-03-01')), 
                               date.end = as.Date(c('2020-02-29','2020-05-31')))

param$time.period[, doy.start := as.numeric(date.start-as.Date('2019-12-31'))]
param$time.period[, doy.end := as.numeric(date.end-as.Date('2019-12-31'))]
param$time.period[, duration := doy.end-doy.start+1]

param$theta = data.table(period = c(1,2),
                         prop.import.M = c(0.2,0.1), 
                         R = c(1.5,1.0), 
                         eff.contact.tracing = c(0.8,0.6),
                         eff.case.finding = c(0.5,0.8))

daily.arrival.import.N = c(rep(20, param$time.period$duration[1]), rep(0, param$time.period$duration[2])) # fixed
daily.arrival.import.M = daily.arrival.import.N * rep(param$theta$prop.import.M, times=param$time.period$duration)

t_max = sum(param$time.period$duration)
generation.matrix.imported.N = generation.matrix.imported.M = matrix(0,nrow=t_max,ncol=(t_max+4))
generation.matrix.local.N = generation.matrix.local.M = matrix(0,nrow=t_max,ncol=(t_max+4))
isolation.matrix.local.N = matrix(0,nrow=t_max,ncol=(t_max+4))
for(t in 1:t_max){
  
  generation.matrix.imported.N[t,t:(t+4)] = c(0,0,0.3,0.4,0)
  generation.matrix.imported.M[t,t:(t+4)] = c(0,0,0.3,0.4,0.2)
  generation.matrix.local.N[t,t:(t+4)] = c(0,0.1,0.3,0.4,0)
  generation.matrix.local.M[t,t:(t+4)] = c(0,0.1,0.3,0.4,0.2)
  isolation.matrix.local.N[t,t:(t+4)] = c(1,0,0,0,0)
  # tmp[t,t:(t+20)] = dlnorm(0:20,theta_f[["serial_mean"]],theta_f[["serial_sd"]])
}

generation.matrix = list(generation.matrix.imported.N = generation.matrix.imported.N,
                         generation.matrix.imported.M = generation.matrix.imported.M,
                         generation.matrix.local.N = generation.matrix.local.N,
                         generation.matrix.local.M = generation.matrix.local.M)

isolation.matrix = list(isolation.matrix.local.N = isolation.matrix.local.N)

# calculate the spread of secondary cases arising from imported cases
sec.cases.infector.N = generation.matrix.imported.N * daily.arrival.import.N * rep(param$theta$R, times=param$time.period$duration)  
sec.cases.infector.M = generation.matrix.imported.M * daily.arrival.import.M * rep(param$theta$R, times=param$time.period$duration)  

# daily local notified and missed secondary cases from notified imported infector
sec.cases.N.infector.N = colSums(sec.cases.infector.N)[1:t_max] * rep(param$theta$eff.contact.tracing, times=param$time.period$duration) 
sec.cases.M.infector.N = colSums(sec.cases.infector.N)[1:t_max] * (1-rep(param$theta$eff.contact.tracing, times=param$time.period$duration))
# sec.cases.N.infector.N == linked case

# daily local notified and missed secondary cases from missed imported infector
sec.cases.N.infector.M = colSums(sec.cases.infector.M)[1:t_max] * rep(param$theta$eff.case.finding, times=param$time.period$duration) 
sec.cases.M.infector.M = colSums(sec.cases.infector.M)[1:t_max] * (1-rep(param$theta$eff.contact.tracing, times=param$time.period$duration))
# sec.cases.N.infector.M == unlinked case

# add import-related notified secondary cases (i.e. first generation of local notified cases)
# add import-related missed secondary cases (i.e. first generation of local missed cases)
daily.local.N = sec.cases.N.infector.N + sec.cases.N.infector.M
daily.local.M = sec.cases.M.infector.N + sec.cases.M.infector.M

daily.local.N.linked = sec.cases.N.infector.N
daily.local.N.unlinked = sec.cases.N.infector.M

daily.local.M.linked = sec.cases.M.infector.N
daily.local.M.unlinked = sec.cases.M.infector.M


# Add onwards transmission from subsequent generations 
for(t in 1:t_max){
  
  sec.cases.infector.N = generation.matrix.local.N[t,1:t_max] * daily.local.N[t] * rep(param$theta$R, times=param$time.period$duration) 
  sec.cases.infector.M = generation.matrix.local.M[t,1:t_max] * daily.local.M[t] * rep(param$theta$R, times=param$time.period$duration) 
  
  sec.cases.N.infector.N = sec.cases.infector.N * rep(param$theta$eff.contact.tracing, times=param$time.period$duration) 
  sec.cases.M.infector.N = sec.cases.infector.M * (1-rep(param$theta$eff.contact.tracing, times=param$time.period$duration))
  
  sec.cases.N.infector.M = sec.cases.infector.M * rep(param$theta$eff.case.finding, times=param$time.period$duration) 
  sec.cases.M.infector.M = sec.cases.infector.M * (1-rep(param$theta$eff.contact.tracing, times=param$time.period$duration))
  
  daily.local.N = daily.local.N + sec.cases.N.infector.N + sec.cases.N.infector.M
  daily.local.M = daily.local.M + sec.cases.M.infector.N + sec.cases.M.infector.M
  
  daily.local.N.linked = daily.local.N.linked + sec.cases.N.infector.N
  daily.local.N.unlinked = daily.local.N.unlinked + sec.cases.N.infector.M
  
  daily.local.M.linked = daily.local.M.linked + sec.cases.M.infector.N
  daily.local.M.unlinked = daily.local.M.unlinked + sec.cases.M.infector.M
}

# store simulated data
obs.data = list(daily.local.N = floor(daily.local.N), 
                daily.local.N.linked = floor(daily.local.N.linked), 
                daily.local.N.unlinked = floor(daily.local.N.unlinked),
                daily.arrival.import.N = daily.arrival.import.N)

import.N = data.table(case = 63521+seq_len(15),
                      doy.onset = c(as.numeric(as.Date('2020-02-03')-as.Date('2019-12-31')) + floor(rnorm(7)),
                                    as.numeric(as.Date('2020-03-03')-as.Date('2019-12-31')) + floor(rnorm(8))))
import.N[, doy.arrival:= doy.onset-1]
import.N[, doy.isolate.from.comm:= doy.onset+1]
import.N[, doy.notification:=doy.isolate.from.comm]
obs.data$import.N = import.N

local.N = data.table(case = 63521+seq_len(15),
                     doy.onset = c(as.numeric(as.Date('2020-02-03')-as.Date('2019-12-31')) + floor(rnorm(7)),
                                   as.numeric(as.Date('2020-03-03')-as.Date('2019-12-31')) + floor(rnorm(8))))
local.N[, doy.isolate.from.comm:= doy.onset+1]
local.N[, doy.notification:=doy.isolate.from.comm]
obs.data$local.N = local.N

save(obs.data, file = 'data/test/obs.data.RData')
save(generation.matrix, file = 'data/test/generation.matrix.RData')
save(isolation.matrix, file = 'data/test/isolation.matrix.RData')


# start MCMC runs
set.seed(123)

# load functions
source('codes/v6/covid_missInf_functions.R')

# load parameter list
param = list()

# time periods of analysis
param$time.period = data.table(period = c(1,2),
                               date.start = as.Date(c('2020-01-01','2020-03-01')), 
                               date.end = as.Date(c('2020-02-29','2020-05-31')))

param$time.period[, doy.start := as.numeric(date.start-as.Date('2019-12-31'))]
param$time.period[, doy.end := as.numeric(date.end-as.Date('2019-12-31'))]
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

# load generation matrix
load('data/test/generation.matrix.RData')

# load isolation matrix
load('data/test/isolation.matrix.RData')

# select parameter for period of interest
period.interest = 1
param$time.period = param$time.period[period == period.interest]
param$theta = param$theta[period == period.interest]
param$period.interest = period.interest

# load data from previous period
if(period.interest == 1){
  store.previous.period = list(mod.daily.local.N.linked = 0,
                               mod.daily.local.N.unlinked = 0,
                               mod.daily.local.M.linked = 0,
                               mod.daily.local.M.unlinked = 0,
                               mod.daily.local.N.linked.by.doy.isolate = 0,
                               mod.daily.local.N.unlinked.by.doy.isolate = 0)
  
} else{
  list.file = dir('output processed', pattern = paste('store.period.', period.interest-1, sep = ''))
  load(paste('output processed/', list.file, sep=''))
  store.previous.period = store.period.1
  
  for(i in 1:length(store.previous.period)){
    store.previous.period[[i]] = apply(store.previous.period[[i]], 2, median)
  }
  
}


# unknown model parameter
theta = param$theta[chain==1 & period %in% c(1,2) & variable != 'rep.k']
theta.limits = data.table(lower = rep(0,5*1),
                          upper = rep(c(5,5,1,1, Inf), each =1))

theta[variable == 'prop.import.M' & period == 2, value:=0]

# reset theta
# theta
# theta$value = c(0.5734058,1.7829113,0.1224401,0.8510073)

# test initial 0.1,0.1, 1.5,1.5, 0.5,0.5, 0.5,0.5
# true 0.2,0.1, 1.5,1.0, 0.8,0.6, 0.5,0.8

theta.limits = data.table(lower = rep(0,4),
                          upper = c(5,5,1,1))

# load obs data
load('data/test/obs.data.RData')
names(obs.data)

# actual
obs.data$daily.arrival.import.N
length(obs.data$daily.arrival.import.N)

# assumed
obs.data$daily.arrival.import.N = c(seq(1,20,length.out = 60), rep(0,92))


# set-up output storage
niter = 10000 # 100000
nparam = nrow(theta) 
store = list()

store$theta = matrix(0, niter, nparam)
colnames(store$theta) = theta$variable.period

store$LL = rep(0,niter)
store$accept = rep(0,niter)
store$doy.start = param$time.period$doy.start  
store$doy.end = param$time.period$doy.end
store$doy = min(param$time.period$doy.start):max(param$time.period$doy.end)
store$duration = length(store$doy)

store$obs.daily.local.N.unlinked = obs.data$daily.local.N.unlinked[store$doy]
store$obs.daily.local.N.linked   = obs.data$daily.local.N.linked[store$doy]

# by time of infection
store$mod.daily.local.N.unlinked = matrix(0, niter, length(store$doy)+4) # modelled notified sec cases with missed infector
store$mod.daily.local.N.linked   = matrix(0, niter, length(store$doy)+4) # modelled notified sec cases with notified infector
store$mod.daily.local.M.unlinked = matrix(0, niter, length(store$doy)+4) # modelled missed sec cases with missed infector by day of infection
store$mod.daily.local.M.linked   = matrix(0, niter, length(store$doy)+4) # modelled missed sec cases with notified infector by day of infection

# by time of isolation
store$mod.daily.local.N.unlinked.by.doy.isolate = matrix(0, niter, length(store$doy)+4) 
store$mod.daily.local.N.linked.by.doy.isolate = matrix(0, niter, length(store$doy)+4) 


# set-up current parameter
theta.star = theta$value



# simulate expected daily incidence given current parameter
output = simulateOutbreak(obs.data, store.previous.period, param, theta.star, generation.matrix, isolation.matrix)

# maybe input obs.data instead of daily.import.N

# calculate log likelihood
store$LL[1] = calculateLogLik(store, output, theta.star)

# store outputs
store$theta[1,] = theta.star
store$mod.daily.local.N.unlinked[1,] = output$daily.local.N.unlinked
store$mod.daily.local.N.linked[1,]   = output$daily.local.N.linked
store$mod.daily.local.M.unlinked[1,] = output$daily.local.M.unlinked
store$mod.daily.local.M.linked[1,]   = output$daily.local.M.linked

store$mod.daily.local.N.unlinked.by.doy.isolate[1,] = output$daily.local.N.unlinked.by.doy.isolate
store$mod.daily.local.N.linked.by.doy.isolate[1,] = output$daily.local.N.linked.by.doy.isolate


# set-up covariance matrix
npc = rep(1,nparam)
pmask = NULL
pmask = c( '') # fix parameters
pmask = match(pmask,theta$variable.period)

# npc[pmask]=0
cov.matrix.theta0 = diag(npc)
epsilon0=0.001


accept.rate=0.234

for(iter in 5001:10000){ # niter
  
  print(sprintf("Iteration %s of %.f. Acceptance rate %.3f", iter, niter, accept.rate))
  
  old.theta  = store$theta[iter-1,]
  old.output = output
  old.LL = store$LL[iter-1]
  
  # sample new transmission parameters
  epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept.rate-0.234)*0.999^iter)),1e-6) # Stop epsilon getting too big or small
  cov.matrix.theta = epsilon0*cov.matrix.theta0
  
  theta.star = sampleTheta(old.theta,cov.matrix.theta,theta.limits,pmask) #sample nearby global parameter space
  # theta.star = sampleTheta(old.theta,covmat.proposal,pmask)
  new.output = simulateOutbreak(obs.data, store.previous.period, param, theta.star, generation.matrix, isolation.matrix)
  new.LL = calculateLogLik(store, new.output, theta.star)
  
  # metropolis hasting step
  mh.output = mh(old.theta, theta.star, old.output, new.output, old.LL, new.LL)
  
  # store outputs
  store$theta[iter,] = mh.output$theta
  store$LL[iter] = mh.output$LL
  store$accept[iter] = mh.output$count
  store$mod.daily.local.N.unlinked[iter,] = mh.output$daily.local.N.unlinked
  store$mod.daily.local.N.linked[iter,]   = mh.output$daily.local.N.linked
  store$mod.daily.local.M.unlinked[iter,] = mh.output$daily.local.M.unlinked
  store$mod.daily.local.M.linked[iter,]   = mh.output$daily.local.M.linked
  store$mod.daily.local.N.unlinked.by.doy.isolate[iter,] = mh.output$daily.local.N.unlinked.by.doy.isolate
  store$mod.daily.local.N.linked.by.doy.isolate[iter,] = mh.output$daily.local.N.linked.by.doy.isolate
  
  output = mh.output[4:9]
  
 
  
  # start adaptation after 100 iterations
  # update acceptance rate
  if(iter>=99){
    accept.rate=sum(store$accept[1:iter])/iter
  }
  
  # if(iter %% 10000 == 0){
  #   save(output,file= paste0("output/output_", paramSet, "_iter_", iter,".rdata"))
  # }
  
}
