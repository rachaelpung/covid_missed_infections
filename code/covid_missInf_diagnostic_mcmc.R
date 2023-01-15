# load storage
list.file = dir('output raw/20220521/wild neg binom/', pattern = 'period_4')
chain = list()
list.file

for(f in 1:8){
  load(paste('output raw/20220521/wild neg binom/', list.file[f], sep=''))
  # load(paste('output raw/test/', list.file[f], sep=''))
  chain[[f]] = store$theta
}


# trace plot
niter = 150000
xmax = 150000
nparam = ncol(chain[[1]])
par(mfrow=c(2,3))

source('codes/covid_missInf_plot_colour.R')
chain.colour = rep(c(CONFIG$colsLight2[4],CONFIG$colsLight2[3],CONFIG$colsLight2[8]), length.out = 20)

for(p in 1:nparam){ 
  
  iteration = seq_len(niter)
  
  col.name = colnames(chain[[1]])[p]
  
  if(grepl('prop_import_M', col.name)) {ymin = 0; ymax = 50}
  if(grepl('R', col.name)) {ymin = 0; ymax = 50}
  if(grepl('eff_contact_tracing', col.name)) {ymin = 0; ymax = 1}
  if(grepl('eff_case_finding', col.name)) {ymin = 0; ymax = 1}
  if(grepl('rep_k', col.name)) {ymin = 0; ymax = 50}
  
  for(c in 1:length(chain)){ # length(chain)
    
    if(c==1){
      plot(iteration,  chain[[c]][,p], type = 'l', col=chain.colour[c], 
           xlab='iteration', ylab=col.name, ylim = c(ymin,ymax), xlim=c(0,xmax))
    }else{
      lines(iteration, chain[[c]][,p], col=chain.colour[c])  
    }
  }
  
}

# incidence plot
par(mfrow=c(2,2))
# plot
xmin=5000
plot(store$doy, store$obs_daily_local_N_linked, type = 'l' )
lines(store$doy, apply(store$mod_daily_local_N_linked_by_doy_isolate[xmin:xmax,1:length(store$doy)], 2, median), col='red')
lines(store$doy, apply(store$mod_daily_local_N_linked_by_doy_isolate[xmin:xmax,1:length(store$doy)], 2, quantile, 0.95), col='red')

plot(store$doy, store$obs_daily_local_N_unlinked, type = 'l' )
lines(store$doy, apply(store$mod_daily_local_N_unlinked_by_doy_isolate[xmin:xmax,1:length(store$doy)], 2, median), col='red')
lines(store$doy, apply(store$mod_daily_local_N_unlinked_by_doy_isolate[xmin:xmax,1:length(store$doy)], 2, quantile, 0.95), col='red')

plot(store$doy, store$obs_daily_local_N_unlinked+store$obs_daily_local_N_linked, type = 'l' )
lines(store$doy, apply(store$mod_daily_local_N_unlinked_by_doy_isolate[xmin:xmax,1:length(store$doy)]+
                         store$mod_daily_local_N_linked_by_doy_isolate[xmin:xmax,1:length(store$doy)], 2, median), col='red')
lines(store$doy, apply(store$mod_daily_local_N_unlinked_by_doy_isolate[xmin:xmax,1:length(store$doy)]+
                         store$mod_daily_local_N_linked_by_doy_isolate[xmin:xmax,1:length(store$doy)], 2, quantile, 0.95), col='red')
lines(store$doy, apply(store$mod_daily_local_N_unlinked_by_doy_isolate[xmin:xmax,1:length(store$doy)]+
                         store$mod_daily_local_N_linked_by_doy_isolate[xmin:xmax,1:length(store$doy)], 2, quantile, 0.05), col='red')

store.M = store$mod_daily_local_M_unlinked[xmin:xmax,1:length(store$doy)]+store$mod_daily_local_M_linked[xmin:xmax,1:length(store$doy)]
plot(store$doy, apply(store.M, 2, quantile, 0.75), col='red', type = 'l')
lines(store$doy, apply(store.M, 2, quantile, 0.5), col='red')
lines(store$doy, apply(store.M, 2, quantile, 0.25), col='red')

test=store$mod_daily_local_M_unlinked[xmin:xmax,1:length(store$doy)]+store$mod_daily_local_M_linked[xmin:xmax,1:length(store$doy)]
apply(test, 2, quantile, 0.5)


# autocorrelation plot
par(mfrow=c(2,2))

for(p in 1:nparam){ 
  for(c in 1:1){ # length(chain)
    
    col.name = colnames(chain[[c]])[p]
    lag = acf(chain[[c]][,p], lag.max = 1000, plot = FALSE)$lag
    autocorrelation = acf(chain[[c]][,p], lag.max = 1000, plot = FALSE)$acf
    
    if(any(is.na(autocorrelation))){
      plot(lag, rep(1, length(lag)), col=chain.colour[c], type = 'l', main=col.name) 
    } else{
      plot(lag, autocorrelation, col=chain.colour[c], type = 'l', main=col.name) 
    }
  }
}


# effective sample size
eff.sample.size = matrix(0, nrow = length(chain), ncol = nparam)

for(p in 1:nparam){
  for(c in 1:length(chain)){
    
    n = length(chain[[c]][,p])
    acf = acf(chain[[c]][,p], lag.max = 1000, plot = FALSE)$acf
    acf.sums = acf[1:999] + acf[2:1000]
    
    if(any(is.na(acf.sums))){
      eff.sample.size[c,p] = NA
    } else{
      denominator = 1 + 2*sum(acf[2:(which(acf.sums < 0)[1]+1)])
      n.eff = n/denominator
      eff.sample.size[c,p] = n.eff
    }
  }
}

eff.sample.size = data.table(eff.sample.size)
colnames(eff.sample.size) = colnames(chain[[1]])

# Gelman-Rubin convergence diagnostic
compute.Rhat <- function(chain){
  
  chain.split = rbind(chain[1:(length(chain)/2)],
                      chain[((length(chain)/2)+1):length(chain)])
  
  n = dim(chain.split)[2]
  m = dim(chain.split)[1]
  
  Ws = sapply(1:m, function(x){ var(chain.split[x,]) })
  W = mean(Ws)
  
  chain.means <- sapply(1:m, function(x){ mean(chain.split[x,]) })
  
  B = var(chain.means)
  V = (n-1)/n*W + B
  Rhat = sqrt(V/W)
  
  return(Rhat)
}


total.n = dim(chain[[1]])[1]
terminal.n = seq(10, total.n, 50)
r.hats = matrix(nrow = nparam*length(chain), ncol = length(terminal.n))

for(p in 1:nparam){
  for(c in 1:length(chain)){
    
    r.hats[p+(c-1)*(nparam),] = sapply(1:length(terminal.n), function(x){
      
      truncated.chain = chain[[c]][,p][1:terminal.n[x]]
      compute.Rhat(truncated.chain)
      
    })
  }
}

# acceptance rate
accept.rate = matrix(0, nrow = length(chain), ncol = 1)

compute.accept.rate <- function(chain){
  
  n = length(chain)
  acceptance.rate = (uniqueN(chain)-1)/(n-1)
  return(acceptance.rate)
  
}

for(p in 1:1){
  for(c in 1:length(chain)){
    accept.rate[c,p] = compute.accept.rate(chain[[c]][,p])
  }
}





