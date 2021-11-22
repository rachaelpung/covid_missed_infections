source('codes/v6/covid_missInf_load_library.R')
load('output processed/inc.all.RData')

# prob of sero conversion
prob.sero.convert = 0.871

# pdf incubation period
dist.incub = data.table(time.since.inf = 1:21)
dist.incub[, prob.onset.since.inf := dlnorm(time.since.inf, meanlog = 1.63, sdlog = 0.50, log = FALSE)]
dist.incub[, prob.onset.since.inf := prob.onset.since.inf/sum(prob.onset.since.inf)]

# prob dist of being S+ since time of onset
dist.prob.sero.since.onset = read_csv('input/dist_Sero_Positivity.csv')
fit.spline = gam(PERCENTAGE_POSITIVE~s(TIME), data = dist.prob.sero.since.onset)
predict.spline = predict(fit.spline,newdata=data.frame(TIME=0:60))/100
predict.spline[predict.spline>1] = 1
dist.prob.sero.since.onset = data.table(time.since.onset = 0:60, prob.sero.since.onset = as.vector(predict.spline))

# assume antibody level will decline after 2 months (~8 weeks)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244126
# https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00354-6/fulltext  

# Not used
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0255208


# interpolate between points
dist.prob.sero.since.onset.interpolate = data.table(time.since.onset = c(60,90,120,150,180,210,240,270,330),
                                                    prob.sero.since.onset = c(1,0.58,0.58,0.58,0.58, 0.56,0.56,0.56,0.19))


fit.spline = smooth.spline(dist.prob.sero.since.onset.interpolate$time.since.onset, 
                           dist.prob.sero.since.onset.interpolate$prob.sero.since.onset,
                           df=5)

predict.spline = predict(fit.spline,61:330)$y
dist.prob.sero.since.onset.interpolate = data.table(time.since.onset = 61:330, prob.sero.since.onset = as.vector(predict.spline))
dist.prob.sero.since.onset = rbind(dist.prob.sero.since.onset,dist.prob.sero.since.onset.interpolate)

# convert time since onset to time since infection
dist.prob.sero.since.onset = dist.prob.sero.since.onset[rep(seq_len(.N), times=21),]
dist.prob.sero.since.onset[, time.inf.to.onset:=rep(1:21, each=.N/21)]
dist.prob.sero.since.onset[, time.since.inf := time.inf.to.onset+time.since.onset]
dist.prob.sero.since.onset[, prob.onset.since.inf := rep(dist.incub[,prob.onset.since.inf], each = .N/21)]
dist.prob.sero.since.onset[, prob.sero.since.inf:=prob.onset.since.inf*prob.sero.since.onset]

# prob dist of being S+ since time of infection if successfully seroconvert
dist.prob.sero.since.inf = dist.prob.sero.since.onset[,.(sum(prob.sero.since.inf)), by=.(time.since.inf)]
dist.prob.sero.since.inf = dist.prob.sero.since.inf[time.since.inf <= 330]
setnames(dist.prob.sero.since.inf, old = 'V1', new = 'prob.sero.since.inf')
dist.prob.sero.since.inf = rbind(data.table(time.since.inf=-14:0, prob.sero.since.inf = 0), dist.prob.sero.since.inf)


# validate against population sero survey
doy.test.start = as.numeric(as.Date("2020-09-07") - as.Date('2019-12-31')) 
doy.test.end   = as.numeric(as.Date("2020-10-31") - as.Date('2019-12-31'))
duration.test = doy.test.end - doy.test.start + 1

# prob of being S+ days of testing
prob.sero.pos = data.table(time.inf = rep(1:doy.test.end, each=duration.test),
                           doy.test = rep(doy.test.start:doy.test.end, times=doy.test.end),
                           prob.doy.test = 1/117,
                           prob.sero.convert = prob.sero.convert)

prob.sero.pos[,time.since.inf:= doy.test-time.inf]
prob.sero.pos[dist.prob.sero.since.inf, prob.sero.since.inf:=prob.sero.since.inf, on=c(time.since.inf='time.since.inf')]
prob.sero.pos[is.na(prob.sero.since.inf),prob.sero.since.inf:=0 ]

prob.sero.pos[,detected:=prob.doy.test*prob.sero.convert*prob.sero.since.inf]


prob.sero.pos = prob.sero.pos[,.(sum(detected)), by=.(time.inf)]$V1


# confidence intervals of observed seroprevalence
binom.CI <- function(events, #events = outcomes
                     trials, #number of individuals, test, etc
                     alpha = 0.05){
  
  n <- trials
  x <- events
  p.hat <- x/n
  
  # Calculate upper and lower limit
  upper.lim <- (p.hat + (qnorm(1-(alpha/2))^2/(2*n)) + qnorm(1-(alpha/2)) * sqrt(((p.hat*(1-p.hat))/n) + (qnorm(1-(alpha/2))^2/(4*n^2))))/(1 + (qnorm(1-(alpha/2))^2/(n)))
  
  lower.lim <- (p.hat + (qnorm(alpha/2)^2/(2*n)) + qnorm(alpha/2) * sqrt(((p.hat*(1-p.hat))/n) + (qnorm(alpha/2)^2/(4*n^2))))/(1 + (qnorm(alpha/2)^2/(n)))
  
  
  # Modification for probabilities close to boundaries
  
  if ((n <= 50 & x %in% c(1, 2)) | (n >= 51 & n <= 100 & x %in% c(1:3))) {
    lower.lim <- 0.5 * qchisq(alpha, 2 * x)/n
  }
  
  if ((n <= 50 & x %in% c(n - 1, n - 2)) | (n >= 51 & n <= 100 & x %in% c(n - (1:3)))) {
    upper.lim <- 1 - 0.5 * qchisq(alpha, 2 * (n - x))/n
  }
  
  out <- c(lower.lim,upper.lim)
  names(out) <- c("lower.CI","upper.CI")
  return(out)
  
}

(2/1578)*100
binom.CI(events = 2, trials = 1578)*100

# 0.1267427
# 
# lower.CI   upper.CI 
# 0.03476431 0.46095281 





# validate against model using linked and unlinked cases
model = inc.all$wild.neg.binom.mod.daily.local.N.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.mod.daily.local.N.linked[,1:doy.test.end] +
  inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.mod.daily.local.M.linked[,1:doy.test.end]

model.sero.pos = model %*% prob.sero.pos
quantile(model.sero.pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%       2.5%      97.5% 
# 0.01985097 0.01463559 0.05258105 

quantile(model.sero.pos*100/5300000 - (2/1578)*100, probs = c(0.5,0.025,0.975))


# validate against model using all cases
model = inc.all$wild.neg.binom.total.mod.daily.local.N.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.total.mod.daily.local.N.linked[,1:doy.test.end] +
  inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,1:doy.test.end]

model.sero.pos = model %*% prob.sero.pos
quantile(model.sero.pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%       2.5%      97.5% 
# 0.06741919 0.03610059 0.64898756 

quantile(model.sero.pos*100/5300000 - (2/1578)*100, probs = c(0.5,0.025,0.975))

























### ASSUME NO WANING IMMUNITY 330 DAYS SINCE INFECTION

# prob of sero conversion
prob.sero.convert = 0.871

# pdf incubation period
dist.incub = data.table(time.since.inf = 1:21)
dist.incub[, prob.onset.since.inf := dlnorm(time.since.inf, meanlog = 1.63, sdlog = 0.50, log = FALSE)]
dist.incub[, prob.onset.since.inf := prob.onset.since.inf/sum(prob.onset.since.inf)]

# prob dist of being S+ since time of onset
dist.prob.sero.since.onset = read_csv('input/dist_Sero_Positivity.csv')
fit.spline = gam(PERCENTAGE_POSITIVE~s(TIME), data = dist.prob.sero.since.onset)
predict.spline = predict(fit.spline,newdata=data.frame(TIME=0:60))/100
predict.spline[predict.spline>1] = 1
dist.prob.sero.since.onset = data.table(time.since.onset = 0:60, prob.sero.since.onset = as.vector(predict.spline))

# assume antibody level will decline after 2 months (~8 weeks)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244126
# https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00354-6/fulltext  

# Not used
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0255208


# interpolate between points
dist.prob.sero.since.onset.interpolate = data.table(time.since.onset = c(60,90,120,150,180,210,240,270,330),
                                                    prob.sero.since.onset = c(1,1,1,1,1,1,1,1,1))


fit.spline = smooth.spline(dist.prob.sero.since.onset.interpolate$time.since.onset, 
                           dist.prob.sero.since.onset.interpolate$prob.sero.since.onset,
                           df=5)

predict.spline = predict(fit.spline,61:330)$y
dist.prob.sero.since.onset.interpolate = data.table(time.since.onset = 61:330, prob.sero.since.onset = as.vector(predict.spline))
dist.prob.sero.since.onset = rbind(dist.prob.sero.since.onset,dist.prob.sero.since.onset.interpolate)

# convert time since onset to time since infection
dist.prob.sero.since.onset = dist.prob.sero.since.onset[rep(seq_len(.N), times=21),]
dist.prob.sero.since.onset[, time.inf.to.onset:=rep(1:21, each=.N/21)]
dist.prob.sero.since.onset[, time.since.inf := time.inf.to.onset+time.since.onset]
dist.prob.sero.since.onset[, prob.onset.since.inf := rep(dist.incub[,prob.onset.since.inf], each = .N/21)]
dist.prob.sero.since.onset[, prob.sero.since.inf:=prob.onset.since.inf*prob.sero.since.onset]

# prob dist of being S+ since time of infection if successfully seroconvert
dist.prob.sero.since.inf = dist.prob.sero.since.onset[,.(sum(prob.sero.since.inf)), by=.(time.since.inf)]
dist.prob.sero.since.inf = dist.prob.sero.since.inf[time.since.inf <= 330]
setnames(dist.prob.sero.since.inf, old = 'V1', new = 'prob.sero.since.inf')
dist.prob.sero.since.inf = rbind(data.table(time.since.inf=-14:0, prob.sero.since.inf = 0), dist.prob.sero.since.inf)


# validate against population sero survey
doy.test.start = as.numeric(as.Date("2020-09-07") - as.Date('2019-12-31')) 
doy.test.end   = as.numeric(as.Date("2020-10-31") - as.Date('2019-12-31'))
duration.test = doy.test.end - doy.test.start + 1

# prob of being S+ days of testing
prob.sero.pos = data.table(time.inf = rep(1:doy.test.end, each=duration.test),
                           doy.test = rep(doy.test.start:doy.test.end, times=doy.test.end),
                           prob.doy.test = 1/117,
                           prob.sero.convert = prob.sero.convert)

prob.sero.pos[,time.since.inf:= doy.test-time.inf]
prob.sero.pos[dist.prob.sero.since.inf, prob.sero.since.inf:=prob.sero.since.inf, on=c(time.since.inf='time.since.inf')]
prob.sero.pos[is.na(prob.sero.since.inf),prob.sero.since.inf:=0 ]

prob.sero.pos[,detected:=prob.doy.test*prob.sero.convert*prob.sero.since.inf]


prob.sero.pos = prob.sero.pos[,.(sum(detected)), by=.(time.inf)]$V1

# validate against model using linked and unlinked cases
model = inc.all$wild.neg.binom.mod.daily.local.N.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.mod.daily.local.N.linked[,1:doy.test.end] +
  inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.mod.daily.local.M.linked[,1:doy.test.end]

model.sero.pos = model %*% prob.sero.pos
quantile(model.sero.pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%       2.5%      97.5% 
# 0.03199490 0.02360734 0.08869945

quantile(model.sero.pos*100/5300000 - (2/1578)*100, probs = c(0.5,0.025,0.975))


# validate against model using all cases
model = inc.all$wild.neg.binom.total.mod.daily.local.N.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.total.mod.daily.local.N.linked[,1:doy.test.end] +
  inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,1:doy.test.end] +
  inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,1:doy.test.end]

model.sero.pos = model %*% prob.sero.pos
quantile(model.sero.pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%      2.5%     97.5% 
# 0.11125363 0.05918522 1.07295108  

quantile(model.sero.pos*100/5300000 - (2/1578)*100, probs = c(0.5,0.025,0.975))


binom.CI <- function(events, #events = outcomes
                     trials, #number of individuals, test, etc
                     alpha = 0.05){
  
  n <- trials
  x <- events
  p.hat <- x/n
  
  # Calculate upper and lower limit
  upper.lim <- (p.hat + (qnorm(1-(alpha/2))^2/(2*n)) + qnorm(1-(alpha/2)) * sqrt(((p.hat*(1-p.hat))/n) + (qnorm(1-(alpha/2))^2/(4*n^2))))/(1 + (qnorm(1-(alpha/2))^2/(n)))
  
  lower.lim <- (p.hat + (qnorm(alpha/2)^2/(2*n)) + qnorm(alpha/2) * sqrt(((p.hat*(1-p.hat))/n) + (qnorm(alpha/2)^2/(4*n^2))))/(1 + (qnorm(alpha/2)^2/(n)))
  
  
  # Modification for probabilities close to boundaries
  
  if ((n <= 50 & x %in% c(1, 2)) | (n >= 51 & n <= 100 & x %in% c(1:3))) {
    lower.lim <- 0.5 * qchisq(alpha, 2 * x)/n
  }
  
  if ((n <= 50 & x %in% c(n - 1, n - 2)) | (n >= 51 & n <= 100 & x %in% c(n - (1:3)))) {
    upper.lim <- 1 - 0.5 * qchisq(alpha, 2 * (n - x))/n
  }
  
  out <- c(lower.lim,upper.lim)
  names(out) <- c("lower.CI","upper.CI")
  return(out)
  
}

(2/1578)*100
binom.CI(events = 2, trials = 1578)*100

# 0.1267427
# 
# lower.CI   upper.CI 
# 0.03476431 0.46095281 
