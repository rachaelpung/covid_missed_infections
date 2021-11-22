#' create matrix of generation intervals for imported cases by notification status
#' @param param mixture of fixed and unknown model parameters
#' @param obs.data observed imported and local cases
#' @param notified logical input parameter - T if notified, F if missed

matrixGenerationImport <- function(param, obs.data, notified = T){
  
  time.period = param$time.period.matrix
  t_max = sum(time.period$duration)
  
  gen.mean.log = param$gen.mean.log
  gen.sd.log = param$gen.sd.log
  
  incub.mean.log = param$incub.mean.log
  incub.sd.log = param$incub.sd.log
  
  # set-up generation matrix 
  # generation.matrix =  matrix(0,nrow=t_max,ncol=(t_max+4)) # test only
  generation.matrix =  matrix(0,nrow=t_max,ncol=(t_max+13))
  
  import.N = obs.data$import.N
  import.N = import.N[!is.na(doy.onset)]

  for(p in time.period$period){
    
    doy.start = time.period$doy.start[p]
    doy.end = time.period$doy.end[p]
    subset = import.N[doy.arrival>=doy.start & doy.arrival<=doy.end & is.finite(doy.isolate.from.comm),]
    
    # incubation period distribution
    t = 1:14
    dist.incub = dlnorm(t, meanlog = incub.mean.log, sdlog = incub.sd.log, log = FALSE)
    dist.incub = data.table(day = t, pdf  = dist.incub/sum(dist.incub))
    
    # generation interval distribution
    t = 1:14
    dist.gen = dlnorm(t, meanlog = gen.mean.log, sdlog = gen.sd.log, log = FALSE)
    dist.gen = data.table(day = t, pdf  = dist.gen/sum(dist.gen))
    # or to deconvolute from the serial interval and incubation period
    
    # time from infection to arrival distribution
    subset = subset[rep(seq_len(nrow(subset)), each=14)]
    subset[, doy.infection := doy.onset-rep(1:14, times=.N/14)]
    subset[, duration.incubate := doy.onset-doy.infection] # can only take values more than 1
    
    subset = subset[doy.infection<doy.arrival]  # infection after arrival assumed cannot happen
    subset[dist.incub, prob.infected.days.ago := i.pdf, on=c(duration.incubate='day') ]
    subset[, prob.infected.days.ago:=prob.infected.days.ago/sum(prob.infected.days.ago), by=.(case)]
    
    subset[, duration.infection.to.arrival := doy.arrival - doy.infection]
    dist.infection.to.arrival = subset[,.(sum(prob.infected.days.ago)), by=.(duration.infection.to.arrival)]
    dist.infection.to.arrival = sample(dist.infection.to.arrival$duration.infection.to.arrival, 1000,
                                       prob = dist.infection.to.arrival$ V1,
                                       replace = T)
    
    fit = fitdist(dist.infection.to.arrival, 'exp')
    
    # hist(dist.infection.to.arrival,prob=TRUE,breaks = seq(0,30,1))
    # plot(function(x) dexp(x, rate = fit$estimate[1]), from=0, to=30, col="red", add=TRUE)
    
    t=1:14
    dist.infection.to.arrival = dexp(t, rate = fit$estimate[1], log = FALSE)
    dist.infection.to.arrival = data.table(day = t, pdf.infection.to.arrival  = dist.infection.to.arrival/sum(dist.infection.to.arrival))
    
    
    # work out the generation interval available for infection after left truncation
    dist.gen.imported.missed = dist.infection.to.arrival
    dist.gen.imported.missed[, generation.interval.start := day]
    dist.gen.imported.missed[, generation.interval.end := 14]
    
    dist.gen.imported.missed[, days.at.large := generation.interval.end - generation.interval.start + 1]
    dist.gen.imported.missed=dist.gen.imported.missed[rep(seq_len(nrow(dist.gen.imported.missed)), times = days.at.large)]
    dist.gen.imported.missed[,generation.interval := seq_len(.N), by=.(day)]
    dist.gen.imported.missed[,generation.interval := generation.interval+generation.interval.start-1, by=.(day)]
    dist.gen.imported.missed[dist.gen, pdf.generation.interval := i.pdf, on=c(generation.interval='day')]
    
    dist.gen.imported.missed[, pdf := pdf.infection.to.arrival * pdf.generation.interval]
    dist.gen.imported.missed = dist.gen.imported.missed[,.(sum(pdf)), by=.(generation.interval)]
    setnames(dist.gen.imported.missed, c('day', 'pdf'))
    
    
    # time from infection to isolation
    subset[, duration.infection.to.isolate := doy.isolate.from.comm - doy.infection]
    dist.infection.to.isolate = subset[,.(sum(prob.infected.days.ago)), by=.(duration.infection.to.isolate)]
    dist.infection.to.isolate = sample(dist.infection.to.isolate$duration.infection.to.isolate, 1000, 
                                       prob = dist.infection.to.isolate$V1,
                                       replace = T)
    
    # fit = fitdist(dist.infection.to.isolate, 'lnorm')
    fit = fitdist(dist.infection.to.isolate[dist.infection.to.isolate!=0], 'lnorm')
    
    # hist(dist.infection.to.isolate,prob=TRUE,breaks = c(0:50))
    # plot(function(x) dlnorm(x, meanlog = fit$estimate[1], sdlog =  fit$estimate[2]), from=0, to=50, col="red", add=TRUE)
    
    t=1:14
    dist.infection.to.isolate = dlnorm(t, meanlog = fit$estimate[1], sdlog =  fit$estimate[2], log = FALSE)
    dist.infection.to.isolate = data.table(day = t, pdf.infection.to.isolation  = dist.infection.to.isolate/sum(dist.infection.to.isolate))
    dist.infection.to.isolate[, cdf.infection.to.isolation := cumsum(pdf.infection.to.isolation)]
    
    dist.gen.imported.notified = dist.gen.imported.missed
    dist.gen.imported.notified$pdf = dist.gen.imported.notified$pdf * (1-dist.infection.to.isolate$cdf.infection.to.isolation)
    
    # plot(dist.gen$day, dist.gen$pdf, type='l', col = 'red')
    # lines(dist.gen.imported.missed$day, dist.gen.imported.missed$pdf, col = 'orange')
    # lines(dist.gen.imported.notified$day, dist.gen.imported.notified$pdf, col = 'black')
    # 
    index.start = sum(time.period$duration[1:p]) - time.period$duration[p] + 1
    index.end = sum(time.period$duration[1:p])
    
    for(t in index.start:index.end){
      if(notified == T){
        generation.matrix[t,t:(t+13)] = dist.gen.imported.notified$pdf # from day of arrival
      }else{
        generation.matrix[t,t:(t+13)] = dist.gen.imported.missed$pdf
      }

    }

    # for(t in index.start:index.end){
    #   if(notified == T){
    #     generation.matrix[t,t:(t+4)] = c(0,0,0.3,0.4,0) # from day of arrival
    #   }else{
    #     generation.matrix[t,t:(t+4)] = c(0,0,0.3,0.4,0.2)
    #   }  
    # }
    
  }

  
  return(generation.matrix)
  
}


#' create matrix of generation intervals for local cases by notification status
#' @param param mixture of fixed and unknown model parameters
#' @param obs.data observed imported and local cases
#' @param notified logical input parameter - T if notified, F if missed

matrixGenerationLocal <- function(param, obs.data, notified = T){
  
  time.period = param$time.period.matrix
  t_max = sum(time.period$duration)
  
  gen.mean.log = param$gen.mean.log
  gen.sd.log = param$gen.sd.log
  
  incub.mean.log = param$incub.mean.log
  incub.sd.log = param$incub.sd.log
  
  # set-up generation matrix 
  # generation.matrix =  matrix(0,nrow=t_max,ncol=(t_max+4)) # test only
  generation.matrix =  matrix(0,nrow=t_max,ncol=(t_max+13))
  
  local.N = obs.data$local.N
  local.N = local.N[!is.na(doy.onset)]
  
  for(p in time.period$period){
    
    doy.start = time.period$doy.start[p]
    doy.end = time.period$doy.end[p]
    subset = local.N[doy.isolate.from.comm>=doy.start & doy.isolate.from.comm<=doy.end,]
  
    # incubation period distribution
    t = 1:14
    dist.incub = dlnorm(t, meanlog = incub.mean.log, sdlog = incub.sd.log, log = FALSE)
    dist.incub = data.table(day = t, pdf  = dist.incub/sum(dist.incub))
    
    # generation interval distribution
    t = 1:14
    dist.gen = dlnorm(t, meanlog = gen.mean.log, sdlog = gen.sd.log, log = FALSE)
    dist.gen = data.table(day = t, pdf  = dist.gen/sum(dist.gen))
    # or to deconvolute from the serial interval and incubation period
    
    # time from infection to arrival distribution
    subset = subset[rep(seq_len(nrow(subset)), each=14)]
    subset[, doy.infection := doy.onset-rep(1:14, times=.N/14)]
    subset[, duration.incubate := doy.onset-doy.infection] # can only take values more than 1
    
    subset[dist.incub, prob.infected.days.ago := i.pdf, on=c(duration.incubate='day') ]
    subset[, prob.infected.days.ago:=prob.infected.days.ago/sum(prob.infected.days.ago), by=.(case)]
    
    subset[, duration.infection.to.isolate := doy.isolate.from.comm - doy.infection]
    dist.infection.to.isolate = subset[,.(sum(prob.infected.days.ago)), by=.(duration.infection.to.isolate)]
    dist.infection.to.isolate = sample(dist.infection.to.isolate$duration.infection.to.isolate, 1000,
                                       prob = dist.infection.to.isolate$ V1,
                                       replace = T)
    
    fit = selm(dist.infection.to.isolate ~ 1) # selm is "skew-elliptic lm"
    # summary(fit)
    # hist(dist.infection.to.isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dsn(x, dp=fit@param$dp), from=-25, to=50, col="red", add=TRUE)
    # 
    # fit = fitdist(dist.infection.to.isolate, 'norm')
    # hist(dist.infection.to.isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dnorm(x, mean = fit$estimate[1], sd =  fit$estimate[2]), from=-25, to=50, col="red", add=TRUE)
    
    t=1:14
    dist.infection.to.isolate = dsn(t, dp=fit@param$dp, log = FALSE)
    # dist.infection.to.isolate = dlnorm(t, meanlog = fit$estimate[1], sdlog =  fit$estimate[2], log = FALSE)
    dist.infection.to.isolate = data.table(day = t, pdf.infection.to.isolation  = dist.infection.to.isolate/sum(dist.infection.to.isolate))
    dist.infection.to.isolate[, cdf.infection.to.isolation := cumsum(pdf.infection.to.isolation)]
    
    dist.gen.local.notified = dist.gen.local.missed = dist.gen
    dist.gen.local.notified$pdf = dist.gen.local.notified$pdf * (1-dist.infection.to.isolate$cdf.infection.to.isolation)
    
    
    # plot(dist.gen.local.missed$day, dist.gen.local.missed$pdf, type='l', col = 'red')
    # lines(dist.gen.local.notified$day, dist.gen.local.notified$pdf, col = 'black')
    # 
    index.start = sum(time.period$duration[1:p]) - time.period$duration[p] + 1
    index.end = sum(time.period$duration[1:p])
    
    for(t in index.start:index.end){
      if(notified == T){
        generation.matrix[t,t:(t+13)] = dist.gen.local.notified$pdf # from day of arrival
      }else{
        generation.matrix[t,t:(t+13)] = dist.gen.local.missed$pdf
      }

    }

    # for(t in index.start:index.end){
    #   if(notified == T){
    #     generation.matrix[t,t:(t+4)] = c(0,0,0.3,0.4,0) # from day of arrival
    #   }else{
    #     generation.matrix[t,t:(t+4)] = c(0,0,0.3,0.4,0.2)
    #   }  
    # }
    
  }  
  
  return(generation.matrix)
  
}


#' create matrix of infection to isolation intervals for local cases by notification status
#' @param param mixture of fixed and unknown model parameters
#' @param obs.data observed imported and local cases
#' @param notified logical input parameter - T if notified, F if missed

matrixIsolationLocal <- function(param, obs.data, notified = T){
  
  time.period = param$time.period.matrix
  t_max = sum(time.period$duration)
  
  gen.mean.log = param$gen.mean.log
  gen.sd.log = param$gen.sd.log
  
  incub.mean.log = param$incub.mean.log
  incub.sd.log = param$incub.sd.log
  
  # set-up generation matrix 
  # isolation.matrix =  matrix(0,nrow=t_max,ncol=(t_max+4)) # test only
  isolation.matrix =  matrix(0,nrow=t_max,ncol=(t_max+13))
  
  local.N = obs.data$local.N
  local.N = local.N[!is.na(doy.onset)]
  
  
  for(p in time.period$period){
    
    doy.start = time.period$doy.start[p]
    doy.end = time.period$doy.end[p]
    subset = local.N[doy.isolate.from.comm>=doy.start & doy.isolate.from.comm<=doy.end,]
    
    # incubation period distribution
    t = 1:14
    dist.incub = dlnorm(t, meanlog = incub.mean.log, sdlog = incub.sd.log, log = FALSE)
    dist.incub = data.table(day = t, pdf  = dist.incub/sum(dist.incub))
    
    
    # generation interval distribution
    t = 1:14
    dist.gen = dlnorm(t, meanlog = gen.mean.log, sdlog = gen.sd.log, log = FALSE)
    dist.gen = data.table(day = t, pdf  = dist.gen/sum(dist.gen))
    # or to deconvolute from the serial interval and incubation period
    
    # time from infection to arrival distribution
    subset = subset[rep(seq_len(nrow(subset)), each=14)]
    subset[, doy.infection := doy.onset-rep(1:14, times=.N/14)]
    subset[, duration.incubate := doy.onset-doy.infection] # can only take values more than 1
    
    subset[dist.incub, prob.infected.days.ago := i.pdf, on=c(duration.incubate='day') ]
    subset[, prob.infected.days.ago:=prob.infected.days.ago/sum(prob.infected.days.ago), by=.(case)]
    
    subset[, duration.infection.to.isolate := doy.isolate.from.comm - doy.infection]
    dist.infection.to.isolate = subset[,.(sum(prob.infected.days.ago)), by=.(duration.infection.to.isolate)]
    dist.infection.to.isolate = sample(dist.infection.to.isolate$duration.infection.to.isolate, 1000,
                                       prob = dist.infection.to.isolate$ V1,
                                       replace = T)
    
    fit = selm(dist.infection.to.isolate ~ 1) # selm is "skew-elliptic lm"
    # summary(fit)
    # hist(dist.infection.to.isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dsn(x, dp=fit@param$dp), from=-25, to=50, col="red", add=TRUE)
    # 
    # fit = fitdist(dist.infection.to.isolate, 'norm')
    # hist(dist.infection.to.isolate,prob=TRUE,breaks = seq(-25,50,1))
    # plot(function(x) dnorm(x, mean = fit$estimate[1], sd =  fit$estimate[2]), from=-25, to=50, col="red", add=TRUE)
    
    t=1:14
    dist.infection.to.isolate = dsn(t, dp=fit@param$dp, log = FALSE)
    # dist.infection.to.isolate = dlnorm(t, meanlog = fit$estimate[1], sdlog =  fit$estimate[2], log = FALSE)
    dist.infection.to.isolate = data.table(day = t, pdf.infection.to.isolation  = dist.infection.to.isolate/sum(dist.infection.to.isolate))
    dist.infection.to.isolate[, cdf.infection.to.isolation := cumsum(pdf.infection.to.isolation)]
    
    dist.gen.local.notified = dist.gen.local.missed = dist.gen
    dist.gen.local.notified$pdf = dist.gen.local.notified$pdf * (1-dist.infection.to.isolate$cdf.infection.to.isolation)
    
  
    index.start = sum(time.period$duration[1:p]) - time.period$duration[p] + 1
    index.end = sum(time.period$duration[1:p])
    
    for(t in index.start:index.end){
      if(notified == T){
        isolation.matrix[t,t:(t+13)] = dist.infection.to.isolate$pdf # from day of arrival
      }else{
        isolation.matrix[t,t:(t+13)] = dist.infection.to.isolate$pdf
      }

    }

    # for(t in index.start:index.end){
    #   if(notified == T){
    #     isolation.matrix[t,t:(t+4)] = c(1,0,0,0,0) 
    #   }else{
    #     isolation.matrix[t,t:(t+4)] = c(1,0,0,0,0) 
    #   }  
    # }
    
  }  
  
  return(isolation.matrix)
  
}


#' sample unknown parameters from covariance matrix
#' @param inputTheta current set of model parameters
#' @param covarTheta covariance matrix
#' @param pmask mask parameters that are fixed

sampleTheta <- function(inputTheta, covarTheta, limitsTheta, pmask){
  
  # # sample new parameters from nearby:
  # mean.vector.theta = inputTheta
  # 
  # theta.star = as.numeric(exp(rmvn(1,log(mean.vector.theta), covarTheta)))
  # theta.star[pmask] = mean.vector.theta[pmask]
  # 
  # return(theta.star)
  
  mean.vector.theta = inputTheta
  # covarTheta[,pmask] = 0
  # covarTheta[pmask,] = 0
  diag(covarTheta)[pmask] = 1e-12

  lower.proposal = limitsTheta$lower
  upper.proposal = limitsTheta$upper

  theta.star =  as.vector(rtmvnorm(1,mean = mean.vector.theta, sigma = covarTheta,
                                   lower = lower.proposal, upper = upper.proposal))

  theta.star[pmask] = mean.vector.theta[pmask]
  
  return(theta.star)
}


#' simulate outbreak
#' @param obs.data observed data
#' @param store.previous.period data from previous time period
#' @param theta.star model parameters
#' @param generation.matrix 
#' @param isolation.matrix

simulateOutbreak <- function(obs.data, store.previous.period=NULL, param, theta.star, generation.matrix, isolation.matrix){
  
  period.interest = param$period.interest
  if(period.interest != 1){
    extract = (length(store.previous.period$mod.daily.local.N.linked)-12) : length(store.previous.period$mod.daily.local.N.linked)
  } else{
    extract = 1
  }
  
  
  # format parameter dataframe
  d=length(theta.star)/5
  # d=length(theta.star)/3
  
  time.period = param$time.period
  t_max = sum(time.period$duration)
  theta = data.table(period = time.period$period)
  
  theta[,prop.import.M:= theta.star[1:d]]
  theta[,R:= theta.star[(1+d):(2*d)]]
  theta[,eff.contact.tracing:= theta.star[(1+2*d):(3*d)]]
  theta[,eff.case.finding:= theta.star[(1+3*d):(4*d)]]  #### check this!!!
  # theta[,eff.case.finding:= theta.star[(1+2*d):(3*d)]]
  
  # set-up generation matrix
  # generation.matrix dimension of 349 by 369 first index correspondence to doy 18 last index doy 366
  # extract the relevant portion of the generation matrix for time period of interest
  
  # wild
  # row = (min(time.period$doy.start) - 17):(max(time.period$doy.end) - 17)
  # col = (min(time.period$doy.start) - 17):((max(time.period$doy.end) - 17) + 13 )
  
  # delta
  row = (min(time.period$doy.start) - 456):(max(time.period$doy.end) - 456)
  col = (min(time.period$doy.start) - 456):((max(time.period$doy.end) - 456) + 13 )
  
  # # test
  # row = min(time.period$doy.start):max(time.period$doy.end)
  # col = min(time.period$doy.start):(max(time.period$doy.end) + 4 )
  
  generation.matrix.imported.N = generation.matrix$generation.matrix.imported.N[row,col]
  generation.matrix.imported.M = generation.matrix$generation.matrix.imported.M[row,col]
  generation.matrix.local.N = generation.matrix$generation.matrix.local.N[row,col]
  generation.matrix.local.M = generation.matrix$generation.matrix.local.M[row,col]
  
  # generate missed imported cases
  # obs.data$daily.arrival.import.N since 1 Jan 2020
  daily.arrival.import.N = obs.data$daily.arrival.import.N.ma[min(time.period$doy.start):max(time.period$doy.end)]
  daily.arrival.import.M = daily.arrival.import.N * rep(theta$prop.import.M, times=time.period$duration)
  # daily.arrival.import.M = rep(theta$prop.import.M, times=time.period$duration)
  
  # calculate the spread of secondary cases arising from imported cases
  sec.cases.infector.N = generation.matrix.imported.N * daily.arrival.import.N * rep(theta$R, times=time.period$duration)  
  sec.cases.infector.M = generation.matrix.imported.M * daily.arrival.import.M * rep(theta$R, times=time.period$duration)  
  
  # daily local notified and missed secondary cases from notified imported infector
  sec.cases.N.infector.N = colSums(sec.cases.infector.N) * rep(theta$eff.contact.tracing, times=time.period$duration+13)
  sec.cases.M.infector.N = colSums(sec.cases.infector.N) * (1-rep(theta$eff.contact.tracing, times=time.period$duration+13))
  # sec.cases.N.infector.N == linked case
  
  # test 
  # sec.cases.N.infector.N = colSums(sec.cases.infector.N) * rep(theta$eff.contact.tracing, times=param$time.period$duration+4)
  # sec.cases.M.infector.N = colSums(sec.cases.infector.N) * (1-rep(theta$eff.contact.tracing, times=param$time.period$duration+4))
  
  # daily local notified and missed secondary cases from missed imported infector
  sec.cases.N.infector.M = colSums(sec.cases.infector.M) * rep(theta$eff.case.finding, times=time.period$duration+13)
  sec.cases.M.infector.M = colSums(sec.cases.infector.M) * (1-rep(theta$eff.contact.tracing, times=time.period$duration+13))
  # sec.cases.N.infector.M == unlinked case
  
  # test
  # sec.cases.N.infector.M = colSums(sec.cases.infector.M) * rep(theta$eff.case.finding, times=param$time.period$duration+4)
  # sec.cases.M.infector.M = colSums(sec.cases.infector.M) * (1-rep(theta$eff.contact.tracing, times=param$time.period$duration+4))
  
  
  # add import-related notified secondary cases (i.e. first generation of local notified cases)
  # add import-related missed secondary cases (i.e. first generation of local missed cases)
  daily.local.N = sec.cases.N.infector.N + sec.cases.N.infector.M + 
                  c(store.previous.period$mod.daily.local.N.linked[extract], rep(0,length(sec.cases.N.infector.N)-length(extract))) +
                  c(store.previous.period$mod.daily.local.N.unlinked[extract], rep(0,length(sec.cases.N.infector.N)-length(extract)))
    
  daily.local.M = sec.cases.M.infector.N + sec.cases.M.infector.M + 
                  c(store.previous.period$mod.daily.local.M.linked[extract], rep(0,length(sec.cases.M.infector.N)-length(extract))) +
                  c(store.previous.period$mod.daily.local.M.unlinked[extract], rep(0,length(sec.cases.M.infector.N)-length(extract)))
  
  daily.local.N.linked = sec.cases.N.infector.N + c(store.previous.period$mod.daily.local.N.linked[extract], rep(0,length(sec.cases.N.infector.N)-length(extract))) 
  daily.local.N.unlinked = sec.cases.N.infector.M + c(store.previous.period$mod.daily.local.N.unlinked[extract], rep(0,length(sec.cases.N.infector.M)-length(extract))) 
  
  daily.local.M.linked = sec.cases.M.infector.N + c(store.previous.period$mod.daily.local.M.linked[extract], rep(0,length(sec.cases.M.infector.N)-length(extract))) 
  daily.local.M.unlinked = sec.cases.M.infector.M + c(store.previous.period$mod.daily.local.M.unlinked[extract], rep(0,length(sec.cases.M.infector.M)-length(extract))) 
  
  # Add onwards transmission from subsequent generations 
  for(t in 1:t_max){
    
    sec.cases.infector.N = generation.matrix.local.N[t,] * daily.local.N[t] * rep(theta$R, times=time.period$duration+13)
    sec.cases.infector.M = generation.matrix.local.M[t,] * daily.local.M[t] * rep(theta$R, times=time.period$duration+13)

    sec.cases.N.infector.N = sec.cases.infector.N * rep(theta$eff.contact.tracing, times=time.period$duration+13)
    sec.cases.M.infector.N = sec.cases.infector.M * (1-rep(theta$eff.contact.tracing, times=time.period$duration+13))

    sec.cases.N.infector.M = sec.cases.infector.M * rep(theta$eff.case.finding, times=time.period$duration+13)
    sec.cases.M.infector.M = sec.cases.infector.M * (1-rep(theta$eff.contact.tracing, times=time.period$duration+13))
    
    # test
    # sec.cases.infector.N = generation.matrix.local.N[t,] * daily.local.N[t] * rep(theta$R, times=param$time.period$duration+4)
    # sec.cases.infector.M = generation.matrix.local.M[t,] * daily.local.M[t] * rep(theta$R, times=param$time.period$duration+4)
    # 
    # sec.cases.N.infector.N = sec.cases.infector.N * rep(theta$eff.contact.tracing, times=param$time.period$duration+4)
    # sec.cases.M.infector.N = sec.cases.infector.M * (1-rep(theta$eff.contact.tracing, times=param$time.period$duration+4))
    # 
    # sec.cases.N.infector.M = sec.cases.infector.M * rep(theta$eff.case.finding, times=param$time.period$duration+4)
    # sec.cases.M.infector.M = sec.cases.infector.M * (1-rep(theta$eff.contact.tracing, times=param$time.period$duration+4))

    daily.local.N = daily.local.N + sec.cases.N.infector.N + sec.cases.N.infector.M
    daily.local.M = daily.local.M + sec.cases.M.infector.N + sec.cases.M.infector.M
    
    daily.local.N.linked = daily.local.N.linked + sec.cases.N.infector.N
    daily.local.N.unlinked = daily.local.N.unlinked + sec.cases.N.infector.M
    
    daily.local.M.linked = daily.local.M.linked + sec.cases.M.infector.N
    daily.local.M.unlinked = daily.local.M.unlinked + sec.cases.M.infector.M
    
  }
  
  # notified cases are based on infection dates 
  # need to convert to date of isolation
  isolation.matrix.local.N = isolation.matrix$isolation.matrix.local.N[row,col]
  
  daily.local.N.linked.by.doy.isolate = isolation.matrix.local.N * daily.local.N.linked[1:length(row)]
  daily.local.N.linked.by.doy.isolate = colSums(daily.local.N.linked.by.doy.isolate) +
                                        c(store.previous.period$mod.daily.local.N.linked.by.doy.isolate[extract], rep(0,length(daily.local.N)-length(extract)))
  
  
  daily.local.N.unlinked.by.doy.isolate = isolation.matrix.local.N * daily.local.N.unlinked[1:length(row)] 
  daily.local.N.unlinked.by.doy.isolate = colSums(daily.local.N.unlinked.by.doy.isolate) +
                                          c(store.previous.period$mod.daily.local.N.unlinked.by.doy.isolate[extract], rep(0,length(daily.local.N)-length(extract)))
  
  list(daily.local.N.unlinked = daily.local.N.unlinked, 
       daily.local.N.linked = daily.local.N.linked,
       daily.local.M.unlinked = daily.local.M.unlinked,
       daily.local.M.linked = daily.local.M.linked,
       daily.local.N.unlinked.by.doy.isolate = daily.local.N.unlinked.by.doy.isolate,
       daily.local.N.linked.by.doy.isolate = daily.local.N.linked.by.doy.isolate)
  
}


#' calculate log likelihood
#' @param obs.data observed data
#' @param output model output

calculateLogLik <- function(store, output, theta.star){
  
  duration = store$duration
  check.length = length(store$obs.daily.local.N.linked) != length(output$daily.local.N.linked[1:duration]) 
  if(check.length){print('modelled output has different length from observed data')}
  
  
  # LL = sum(dnbinom(store$obs.daily.local.N.unlinked,mu=output$daily.local.N.unlinked.by.doy.isolate[1:duration],size=1/theta.star[5],log=T),
  #          dnbinom(store$obs.daily.local.N.linked,mu=output$daily.local.N.linked.by.doy.isolate[1:duration],size=1/theta.star[5],log=T))
  
  # LL = sum(dpois((store$obs.daily.local.N.unlinked+
  #                   store$obs.daily.local.N.linked),
  #                lambda=(output$daily.local.N.unlinked.by.doy.isolate[1:duration]+
  #                          output$daily.local.N.linked.by.doy.isolate[1:duration]),log=T))
  
  
  LL = sum(dnbinom((store$obs.daily.local.N.unlinked+
                    store$obs.daily.local.N.linked),
                 mu=(output$daily.local.N.unlinked.by.doy.isolate[1:duration]+
                           output$daily.local.N.linked.by.doy.isolate[1:duration]),size=1/theta.star[5],log=T))
  
  # LL = sum(dpois(store$obs.daily.local.N.unlinked, lambda=output$daily.local.N.unlinked.by.doy.isolate[1:duration],log=T),
  #          dpois(store$obs.daily.local.N.linked, lambda=output$daily.local.N.linked.by.doy.isolate[1:duration],log=T))
  
  return(LL)
  
}


#' metropolis hastings step
#' @param old.theta old model parameter
#' @param theta.star new model parameter 
#' @param old.output old model output 
#' @param new.output new model output
#' @param old.LL old likelihood
#' @param new.LL new likelihood

mh <- function(old.theta, theta.star, old.output, new.output, old.LL, new.LL){
  
  reject = FALSE
  count	= 0
  
  # d=length(theta.star)/4
  # 
  # if(any(theta.star < 0)) reject =TRUE
  # if(any(theta.star[1:d] > 1 )) reject =TRUE # reject if proportion of missed daily imported > 1
  # if(any(theta.star[(1+d):(2*d)] > 5)) reject =TRUE # reject if R > 5
  # if(any(theta.star[(1+2*d):(3*d)] > 1 )) reject =TRUE # reject if effectiveness of contact tracing > 1
  # if(any(theta.star[(1+3*d):(4*d)] > 1)) reject =TRUE # reject if effectiveness of case finding > 1
  
  # if(!reject){
    
    differenceLL = new.LL-old.LL 
    
    if(log(runif(1))>differenceLL){
      reject=TRUE
    } else{
      count	= 1
    }
  # }
  
  
  if(reject) return(c(list(theta = old.theta,
                           count = count,
                           LL = old.LL),
                           old.output))
  
  return(c(list(theta = theta.star,
                count = count,
                LL = new.LL),
                new.output))
  
}
