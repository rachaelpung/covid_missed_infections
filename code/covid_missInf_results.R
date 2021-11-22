source('codes/v6/covid_missInf_load_library.R')

load('input/obs.data.RData')
load('output processed/inc.all.RData')
load('output processed/theta.all.RData')

# load time period
time.period = read_excel('data/param.xlsx', sheet = 'time.period.param')
time.period = data.table(time.period)
setnames(time.period, c('period', 'date.start', 'date.end'))

time.period[, variant:=gsub('(.*)\\d', '\\1',period)]
time.period[, period:=gsub(variant,'',period), by=seq_len(nrow(time.period))]
time.period[, period:=as.numeric(period)]

time.period[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
time.period[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
time.period[, duration := doy.end-doy.start+1]

# tabulate cases by time period
period.start = c(time.period$doy.start[1], time.period$doy.start[1:5],
                 time.period$doy.start[1], time.period$doy.start[1:5],
                 time.period$doy.start[6], time.period$doy.start[6:10])
period.end = c(time.period$doy.end[5], time.period$doy.end[1:5],
               time.period$doy.end[5], time.period$doy.end[1:5],
               time.period$doy.end[10], time.period$doy.end[6:10])

case.table = data.table(period.variant = rep(c('wild.neg.binom.mod','wild.neg.binom.total','delta.neg.binom.total.mod'), each=6), period.start = period.start, period.end = period.end)
case.table[, import.N.shn := sum(obs.data$daily.arrival.import.N.shn[period.start:period.end]), by=seq_len(nrow(case.table))]
case.table[, import.N.no.shn := sum(obs.data$daily.arrival.import.N.no.shn[period.start:period.end]), by=seq_len(nrow(case.table))]
case.table[, local.N.linked := sum(obs.data$daily.local.N.linked[period.start:period.end]), by=seq_len(nrow(case.table))]
case.table[, local.N.unlinked := sum(obs.data$daily.local.N.unlinked[period.start:period.end]), by=seq_len(nrow(case.table))]

case.table[, local.N.death.by.inf := sum(obs.data$daily.local.N.death.by.inf.doy[period.start:period.end]), by=seq_len(nrow(case.table))]
case.table[, local.N.icu.by.inf := sum(obs.data$daily.local.N.icu.by.inf.doy[period.start:period.end]), by=seq_len(nrow(case.table))]
case.table[, local.N.death := sum(obs.data$daily.local.N.death[period.start:period.end]), by=seq_len(nrow(case.table))]
case.table[, local.N.icu := sum(obs.data$daily.local.N.icu[period.start:period.end]), by=seq_len(nrow(case.table))]



# wild type missed
case.table[, period.missed.start := period.start]
case.table[, period.missed.end := period.end]
case.table[period.variant == 'delta.neg.binom.total.mod', period.missed.start := period.missed.start-456]
case.table[period.variant == 'delta.neg.binom.total.mod', period.missed.end := period.missed.end-456]


for(v in case.table[,unique(period.variant)]){
  
  subset.inc.all = inc.all[grep(v, names(inc.all))]
  
  # daily.local.M.unlinked + daily.local.M.linked
  case.table[period.variant == v, local.M.med := quantile(apply(subset.inc.all[[3]][,period.missed.start:period.missed.end],1,sum)+
                                                          apply(subset.inc.all[[4]][,period.missed.start:period.missed.end],1,sum),
                                                          probs = c(0.5)),
             by=seq_len(nrow(case.table[period.variant == v,]))]
  
  case.table[period.variant == v, local.M.lwr := quantile(apply(subset.inc.all[[3]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[4]][,period.missed.start:period.missed.end],1,sum),
                                                          probs = c(0.025)),
             by=seq_len(nrow(case.table[period.variant == v,]))]
  
  case.table[period.variant == v, local.M.upp := quantile(apply(subset.inc.all[[3]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[4]][,period.missed.start:period.missed.end],1,sum),
                                                          probs = c(0.975)),
             by=seq_len(nrow(case.table[period.variant == v,]))]
  
  # daily.local.M.unlinked + daily.local.M.linked + daily.local.N.unlinked + daily.local.N.linked
  case.table[period.variant == v, local.T.med := quantile(apply(subset.inc.all[[3]][,period.missed.start:period.missed.end],1,sum)+
                                                          apply(subset.inc.all[[4]][,period.missed.start:period.missed.end],1,sum)+
                                                          apply(subset.inc.all[[1]][,period.missed.start:period.missed.end],1,sum)+
                                                          apply(subset.inc.all[[2]][,period.missed.start:period.missed.end],1,sum),
                                                          probs = c(0.5)),
             by=seq_len(nrow(case.table[period.variant == v,]))]
  
  case.table[period.variant == v, local.T.lwr := quantile(apply(subset.inc.all[[3]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[4]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[1]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[2]][,period.missed.start:period.missed.end],1,sum),
                                                          probs = c(0.025)),
             by=seq_len(nrow(case.table[period.variant == v,]))]
  
  case.table[period.variant == v, local.T.upp := quantile(apply(subset.inc.all[[3]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[4]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[1]][,period.missed.start:period.missed.end],1,sum)+
                                                            apply(subset.inc.all[[2]][,period.missed.start:period.missed.end],1,sum),
                                                          probs = c(0.975)),
             by=seq_len(nrow(case.table[period.variant == v,]))]
  
  
}


case.table[,CFR.med := local.N.death.by.inf*100/(local.T.med-local.M.med)]
case.table[,CFR.lwr := local.N.death.by.inf*100/(local.T.lwr-local.M.lwr)]
case.table[,CFR.upp := local.N.death.by.inf*100/(local.T.upp-local.M.upp)]

case.table[,IFR.med := local.N.death.by.inf*100/local.T.med]
case.table[,IFR.lwr := local.N.death.by.inf*100/local.T.lwr]
case.table[,IFR.upp := local.N.death.by.inf*100/local.T.upp]

case.table[,CIR.med := local.N.icu.by.inf*100/(local.T.med-local.M.med)]
case.table[,CIR.lwr := local.N.icu.by.inf*100/(local.T.lwr-local.M.lwr)]
case.table[,CIR.upp := local.N.icu.by.inf*100/(local.T.upp-local.M.upp)]

case.table[,IIR.med := local.N.icu.by.inf*100/local.T.med]
case.table[,IIR.lwr := local.N.icu.by.inf*100/local.T.lwr]
case.table[,IIR.upp := local.N.icu.by.inf*100/local.T.upp]


View(case.table)


# incidence for wild-type SARS-CoV-2

# notified cases before lockdown
sum(obs.data$daily.local.N[1:97])

# using linked and unlinked cases for model fitting

# missed cases before lockdown
quantile(apply(inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,61:97],1,sum)+
         apply(inc.all$wild.neg.binom.mod.daily.local.M.linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))

# daily missed cases one week before lockdown
quantile(apply(inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,91:97],1,sum) + 
           apply(inc.all$wild.neg.binom.mod.daily.local.M.linked[,91:97],1,sum), 
         probs = c(0.5,0.025,0.975))/7

# daily missed cases before lockdown
quantile(apply(inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,61:97],1,sum)+
         apply(inc.all$wild.neg.binom.mod.daily.local.M.linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))/37

# daily missed cases during partial lockdown and phase 1 reopening
quantile(apply(inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,98:170],1,sum)+
        apply(inc.all$wild.neg.binom.mod.daily.local.M.linked[,98:170],1,sum),
        probs=c(0.5,0.025, 0.975))/73

# daily missed cases during partial lockdown and phase 1 reopening
quantile(apply(inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,171:366],1,sum)+
           apply(inc.all$wild.neg.binom.mod.daily.local.M.linked[,171:366],1,sum),
         probs=c(0.5,0.025, 0.975))/196

# missed cases in 2020
quantile(apply(inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,1:366],1,sum)+
           apply(inc.all$wild.neg.binom.mod.daily.local.M.linked[,1:366],1,sum),
         probs=c(0.5,0.025, 0.975))


# using cases with no info on linkage for model fitting

# missed cases before lockdown
quantile(apply(inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,61:97],1,sum)+
            apply(inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))

# daily missed cases before lockdown
quantile(apply(inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,61:97],1,sum)+
            apply(inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))/37

# daily missed cases one week before lockdown
quantile(apply(inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,91:97],1,sum) + 
           apply(inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,91:97],1,sum), 
         probs = c(0.5,0.025,0.975))/7

# daily missed cases during partial lockdown and phase 1 reopening
quantile(apply(inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,98:170],1,sum)+
           apply(inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,98:170],1,sum),
         probs=c(0.5,0.025, 0.975))/73

# difference between models in estimating missed infections during lockdown period
quantile(((apply(inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,98:170],1,sum)+
           apply(inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,98:170],1,sum))/73 )/
           ((apply(inc.all$wild.neg.binom.mod.daily.local.M.unlinked[,98:170],1,sum)+
              apply(inc.all$wild.neg.binom.mod.daily.local.M.linked[,98:170],1,sum))/73),
         probs=c(0.5,0.025, 0.975))

# missed cases in 2020
quantile(apply(inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked[,1:366],1,sum)+
            apply(inc.all$wild.neg.binom.total.mod.daily.local.M.linked[,1:366],1,sum),
         probs=c(0.5,0.025, 0.975))


# estimating parameters

# start of the pandemic
apply(theta.all$wild.neg.binom.period.1, 2, quantile, probs = c(0.5,0.025,0.975))

# March to April
apply(theta.all$wild.neg.binom.period.2, 2, quantile, probs = c(0.5,0.025,0.975))

apply(theta.all$wild.neg.binom.period.2[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs.data$daily.arrival.import.N[61:97])

# Partial lockdown
apply(theta.all$wild.neg.binom.period.3, 2, quantile, probs = c(0.5,0.025,0.975))

# reopening
apply(theta.all$wild.neg.binom.period.4, 2, quantile, probs = c(0.5,0.025,0.975))

apply(theta.all$wild.neg.binom.period.4[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs.data$daily.arrival.import.N[171:194])

apply(theta.all$wild.neg.binom.total.period.4, 2, quantile, probs = c(0.5,0.025,0.975))



# incidence for delta-variant SARS-CoV-2
sum(obs.data$daily.local.N[457:501])/45

# early importation of delta cases
apply(theta.all$delta.neg.binom.total.period.1[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs.data$daily.arrival.import.N[457:501])
apply(theta.all$delta.neg.binom.total.period.1, 2, quantile, probs = c(0.5,0.025,0.975))

# missed cases in early delta outbreak
quantile(apply(inc.all$delta.neg.binom.total.mod.daily.local.M.unlinked[,1:45],1,sum)+
           apply(inc.all$delta.neg.binom.total.mod.daily.local.M.linked[,1:45],1,sum),
         probs=c(0.5,0.025, 0.975))

# decline phase
apply(theta.all$delta.neg.binom.total.period.3, 2, quantile, probs = c(0.5,0.025,0.975))

# surge phase
apply(theta.all$delta.neg.binom.total.period.4, 2, quantile, probs = c(0.5,0.025,0.975))

quantile(apply(inc.all$delta.neg.binom.total.mod.daily.local.M.unlinked[,103:140],1,sum)+
           apply(inc.all$delta.neg.binom.total.mod.daily.local.M.linked[,103:140],1,sum),
         probs=c(0.5,0.025, 0.975))/38
