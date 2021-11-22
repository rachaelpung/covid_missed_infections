source('codes/v6/covid_missInf_load_library.R')
source('codes/v6/covid_missInf_functions.R')

# start MCMC runs
set.seed(123)

# load obs data
load('input/obs.data.RData')

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


# set-up generation matrix
generation.matrix.imported.N = matrixGenerationImport(param, obs.data, notified = T)
generation.matrix.imported.M = matrixGenerationImport(param, obs.data, notified = F)
generation.matrix.local.N = matrixGenerationLocal(param, obs.data, notified = T)
generation.matrix.local.M = matrixGenerationLocal(param, obs.data, notified = F)

generation.matrix = list(generation.matrix.imported.N = generation.matrix.imported.N, 
                         generation.matrix.imported.M = generation.matrix.imported.M, 
                         generation.matrix.local.N = generation.matrix.local.N, 
                         generation.matrix.local.M = generation.matrix.local.M)

# set-up isolation matrix
isolation.matrix.local.N = matrixIsolationLocal(param, obs.data, notified = T)
isolation.matrix = list(isolation.matrix.local.N = isolation.matrix.local.N)



# select parameter for period of interest
period.interest = 4
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
  list.file = dir('output processed', pattern = paste('store.period.', period.interest-1, '.', param$variant, sep = ''))
  load(paste('output processed/', list.file, sep=''))
  store.previous.period = store.period
  
  for(i in 1:length(store.previous.period)){
    store.previous.period[[i]] = apply(store.previous.period[[i]], 2, median)
  }
  
}


# start MCMC runs
seed = scan('input/setSeed.txt')

# set up clusters
cl = makeCluster(20) #detectCores()
registerDoParallel(cl)

foreach(param.set = 1:20, .packages = c("data.table","tmvtnorm"), .combine = "c") %dopar%  {

set.seed(seed[param.set])
  
# unknown model parameter
theta = param$theta[chain == param.set]
# theta = param$theta[chain == param.set & variable!='rep.k']
# theta=theta[1:3,]

# theta$value[1] = 0.5 ############# CHECK
theta$value[1] = 0
# theta$value[2]= runif(1, min = 0, max=1)
# theta$value[3]= runif(1, min = 0, max=1)
# theta$value[4]= runif(1, min = 0, max=1)
theta.limits = data.table(lower = c(0,0,0,0,0),
                          upper = c(5,5,1,1,Inf))

# theta.limits = data.table(lower = rep(0,3),
#                           upper = c(10^12,5,1))



# set-up output storage
niter = 60000 # 100000
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
store$mod.daily.local.N.unlinked = matrix(0, niter, length(store$doy)+13) # modelled notified sec cases with missed infector
store$mod.daily.local.N.linked   = matrix(0, niter, length(store$doy)+13) # modelled notified sec cases with notified infector
store$mod.daily.local.M.unlinked = matrix(0, niter, length(store$doy)+13) # modelled missed sec cases with missed infector by day of infection
store$mod.daily.local.M.linked   = matrix(0, niter, length(store$doy)+13) # modelled missed sec cases with notified infector by day of infection

# by time of isolation
store$mod.daily.local.N.unlinked.by.doy.isolate = matrix(0, niter, length(store$doy)+13) 
store$mod.daily.local.N.linked.by.doy.isolate = matrix(0, niter, length(store$doy)+13) 


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
# pmask = NULL
pmask = c('prop.import.M.period.5') # fix parameters
pmask = match(pmask,theta$variable.period)

# npc[pmask]=0
cov.matrix.theta0 = diag(npc)
epsilon0=0.001
accept.rate=0.234

old.output = output

for(iter in 2:niter){ # niter
  
  print(sprintf("Iteration %s of %.f. Acceptance rate %.3f", iter, niter, accept.rate))
  
  old.theta  = store$theta[iter-1,]
  # old.output = output
  old.LL = store$LL[iter-1]
  
  # sample new transmission parameters
  epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept.rate-0.234)*0.999^iter)),1e-6) # Stop epsilon getting too big or small
  cov.matrix.theta = epsilon0*cov.matrix.theta0
  
  theta.star = sampleTheta(old.theta,cov.matrix.theta,theta.limits,pmask) #sample nearby global parameter space
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
  store$mod.daily.local.N.linked.by.doy.isolate[iter,]   = mh.output$daily.local.N.linked.by.doy.isolate
  
  old.output = mh.output[4:9]
  
  # start adaptation after 100 iterations
  # update acceptance rate
  if(iter>100){
    accept.rate=sum(store$accept[1:iter])/iter
  }
  
  if(iter %% 10000 == 0){ # 10000
    if(param.set<10){
      save(store,file= paste0("output raw/20210906 test runs/output_variant_", param$variant,"_period_", param$period.interest, '_paramset_0', param.set, "_iter_", iter,".rdata"))
    } else{
      save(store,file= paste0("output raw/20210906 test runs/output_variant_", param$variant,"_period_", param$period.interest, '_paramset_', param.set, "_iter_", iter,".rdata"))
    }
  }
  
}

}


# stop clusters
stopCluster(cl)
