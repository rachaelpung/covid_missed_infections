# combine and thin chains for one time period
# create variable list
v = c('theta', 'LL', 'mod_daily_local_N_unlinked','mod_daily_local_N_linked',
      'mod_daily_local_M_unlinked', 'mod_daily_local_M_linked', 
      'mod_daily_local_N_unlinked_by_doy_isolate', 'mod_daily_local_N_linked_by_doy_isolate')
for(n in 1:length(v)){ assign(v[n], list())}

period_num=1
period_variant='W'

# load all chains outputs
list.file = dir('output raw/20220521/wild neg binom/', pattern = paste(period_variant, '_period_', period_num, sep=''))
list.file

for(f in 1:length(list.file)){
  load(paste('output raw/20220521/wild neg binom/', list.file[f], sep=''))
  
  # burn-in first 5000 and thin every 50th sample
  retain = seq(5001,150000, 50) #60000, 50
  # retain = seq(2001,10000, 50)
  
  theta[[f]] = as.data.table(store$theta[retain,])
  LL[[f]] = as.data.table(store$LL[retain])
  mod_daily_local_N_unlinked[[f]] = as.data.table(store$mod_daily_local_N_unlinked[retain,])
  mod_daily_local_N_linked[[f]] = as.data.table(store$mod_daily_local_N_linked[retain,])
  mod_daily_local_M_unlinked[[f]] = as.data.table(store$mod_daily_local_M_unlinked[retain,])
  mod_daily_local_M_linked[[f]] = as.data.table(store$mod_daily_local_M_linked[retain,])
  mod_daily_local_N_unlinked_by_doy_isolate[[f]] = as.data.table(store$mod_daily_local_N_unlinked_by_doy_isolate[retain,])
  mod_daily_local_N_linked_by_doy_isolate[[f]] = as.data.table(store$mod_daily_local_N_linked_by_doy_isolate[retain,])

  
}

theta = rbindlist(theta, use.names = T)
LL = rbindlist(LL, use.names = T)
mod_daily_local_N_unlinked = rbindlist(mod_daily_local_N_unlinked, use.names = T)
mod_daily_local_N_linked = rbindlist(mod_daily_local_N_linked, use.names = T)
mod_daily_local_M_unlinked = rbindlist(mod_daily_local_M_unlinked, use.names = T)
mod_daily_local_M_linked = rbindlist(mod_daily_local_M_linked, use.names = T)
mod_daily_local_N_unlinked_by_doy_isolate = rbindlist(mod_daily_local_N_unlinked_by_doy_isolate, use.names = T)
mod_daily_local_N_linked_by_doy_isolate = rbindlist(mod_daily_local_N_linked_by_doy_isolate, use.names = T)

store_period = list(theta = theta,
                      LL = LL,
                      mod_daily_local_N_unlinked = mod_daily_local_N_unlinked,
                      mod_daily_local_N_linked = mod_daily_local_N_linked,
                      mod_daily_local_M_unlinked = mod_daily_local_M_unlinked,
                      mod_daily_local_M_linked = mod_daily_local_M_linked,
                      mod_daily_local_N_unlinked_by_doy_isolate = mod_daily_local_N_unlinked_by_doy_isolate,
                      mod_daily_local_N_linked_by_doy_isolate = mod_daily_local_N_linked_by_doy_isolate)


save(store_period, file=paste('output processed/20220521/wild neg binom/store_period_', period_num, '_', period_variant, '.RData', sep=''))



