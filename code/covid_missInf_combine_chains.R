# create variable list
v = c('theta', 'LL', 'mod.daily.local.N.unlinked','mod.daily.local.N.linked',
      'mod.daily.local.M.unlinked', 'mod.daily.local.M.linked', 
      'mod.daily.local.N.unlinked.by.doy.isolate', 'mod.daily.local.N.linked.by.doy.isolate')
for(n in 1:length(v)){ assign(v[n], list())}

period.num=5
period.variant='D'

# load all chains outputs
list.file = dir('output raw/20210906 test runs/', pattern = paste(period.variant, '_period_', period.num, sep=''))
list.file

for(f in 1:length(list.file)){
  load(paste('output raw/20210906 test runs/', list.file[f], sep=''))
  
  # burn-in first 5000 and thin every 10th sample
  retain = seq(5001,60000, 50) #60000, 50
  # retain = 2000:5000
  
  theta[[f]] = as.data.table(store$theta[retain,])
  LL[[f]] = as.data.table(store$LL[retain])
  mod.daily.local.N.unlinked[[f]] = as.data.table(store$mod.daily.local.N.unlinked[retain,])
  mod.daily.local.N.linked[[f]] = as.data.table(store$mod.daily.local.N.linked[retain,])
  mod.daily.local.M.unlinked[[f]] = as.data.table(store$mod.daily.local.M.unlinked[retain,])
  mod.daily.local.M.linked[[f]] = as.data.table(store$mod.daily.local.M.linked[retain,])
  mod.daily.local.N.unlinked.by.doy.isolate[[f]] = as.data.table(store$mod.daily.local.N.unlinked.by.doy.isolate[retain,])
  mod.daily.local.N.linked.by.doy.isolate[[f]] = as.data.table(store$mod.daily.local.N.linked.by.doy.isolate[retain,])

  
}

theta = rbindlist(theta, use.names = T)
LL = rbindlist(LL, use.names = T)
mod.daily.local.N.unlinked = rbindlist(mod.daily.local.N.unlinked, use.names = T)
mod.daily.local.N.linked = rbindlist(mod.daily.local.N.linked, use.names = T)
mod.daily.local.M.unlinked = rbindlist(mod.daily.local.M.unlinked, use.names = T)
mod.daily.local.M.linked = rbindlist(mod.daily.local.M.linked, use.names = T)
mod.daily.local.N.unlinked.by.doy.isolate = rbindlist(mod.daily.local.N.unlinked.by.doy.isolate, use.names = T)
mod.daily.local.N.linked.by.doy.isolate = rbindlist(mod.daily.local.N.linked.by.doy.isolate, use.names = T)

store.period = list(theta = theta,
                      LL = LL,
                      mod.daily.local.N.unlinked = mod.daily.local.N.unlinked,
                      mod.daily.local.N.linked = mod.daily.local.N.linked,
                      mod.daily.local.M.unlinked = mod.daily.local.M.unlinked,
                      mod.daily.local.M.linked = mod.daily.local.M.linked,
                      mod.daily.local.N.unlinked.by.doy.isolate = mod.daily.local.N.unlinked.by.doy.isolate,
                      mod.daily.local.N.linked.by.doy.isolate = mod.daily.local.N.linked.by.doy.isolate)


save(store.period, file=paste('output processed/store.period.', period.num, '.', period.variant, '.RData', sep=''))
