source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_functions.R')

load('input/obs_data.RData')
load('input/time_period.RData')

# theta posteriors
period=c(1,2,3,4,5)

folder = c('wild neg binom', 'wild neg binom total', 'delta neg binom total')

theta_all = list()
theta_period = list()

for(folder_num in 1:length(folder)){

  print(folder[folder_num])

  for(period_num in period){

    list_file = dir(paste('output raw/20220521/', folder[folder_num], '/', sep =''), pattern = paste('_period_', period_num, sep=''))
    list_file

    if(length(list_file) !=0){
      for(f in 1:length(list_file)){
        
        load(paste('output raw/20220521/', folder[folder_num], '/', list_file[f], sep=''))
        
        # burn-in first 5000 and thinned every 50
        retain = seq(5001,150000,50) 
        theta_period[[f]] = as.data.table(store$theta[retain,])
        theta_period[[f]]$chain = ceiling(f/2)
        theta_period[[f]]$period = period_num
        
        print(f)
      }  
      
      # theta_all[[period_num+(folder_num-1)*5]] = rbindlist(theta_period, use.names = T)
      name_list = gsub(' ','_',paste(folder[folder_num], '_period_', period_num, sep=''))
      theta_all[[name_list]] = rbindlist(theta_period, use.names = T)
    }
    print(period_num)
  }
}

rm(period, period_num, f, folder, folder_num, list_file, retain, name_list, theta_period)


save(theta_all, file='output processed/20220818/theta_all.RData')
save(theta_all, file='output processed/20220521/theta_wild.RData')

# incidence posteriors
extract_wild = seq_len(sum(time_period[variant =='W']$duration))
extract_delta = seq_len(sum(time_period[variant =='D']$duration))

period=c(1,2,3,4,5)

folder = c('wild neg binom', 'wild neg binom total', 'delta neg binom total')

inc_all = list() # incidence of all wild and delta cases
inc_variant = list() # temporary variable, variant specific incidence
inc_variant_period = list() # temporary variable, variant and period specific incidence

for(folder_num in 1:length(folder)){

  print(folder[folder_num])

  for(period_num in period){


    list_file = dir(paste('output processed/20220521/', folder[folder_num], '/', sep =''), pattern = as.character(period_num))
    list_file

    if(length(list_file) !=0){
      load(paste('output processed/20220521/', folder[folder_num], '/', list_file, sep=''))
      
      # remove theta and likelihood estimates
      store_period[[1]] = store_period[[2]] = NULL
      
      if(period_num == 1){
        inc_variant = store_period
      } else{
        
        inc_variant_period = store_period
        
        for(i in 1:length(inc_variant_period)){
          
          
          inc_variant_period[[i]][,1:13] = sweep(inc_variant_period[[i]][,1:13], 2, inc_variant_period_median[[i]], FUN='-')
          extract = (length(inc_variant[[i]])-12) : length(inc_variant[[i]])
          
          
          inc_variant[[i]][,extract] = inc_variant[[i]][,extract,with=FALSE] + inc_variant_period[[i]][,1:13]
          inc_variant[[i]] = cbind(inc_variant[[i]], inc_variant_period[[i]][,14:length(inc_variant_period[[i]])])
          
          if (any(inc_variant[[i]] < 0)){
            print("Negative values detected")
            break
          }
        }
      }
      
      
      
      
      # extract median estimates to substract from the stored outputs of next time period
      inc_variant_period_median = store_period
      
      for(i in 1:length(inc_variant_period_median)){
        inc_variant_period_median[[i]] = apply(inc_variant_period_median[[i]], 2, median)
        extract = (length(inc_variant_period_median[[i]])-12) : length(inc_variant_period_median[[i]])
        inc_variant_period_median[[i]] = inc_variant_period_median[[i]][extract]
      }
      
      
      
      print(period_num)
    }
    
  }

  name_list = paste(gsub(' ','_',folder[folder_num]), names(inc_variant), sep="_")

  names(inc_variant) = name_list

  inc_all = c(inc_all, inc_variant)

}

rm(inc_variant, inc_variant_period, inc_variant_period_median, period, period_num, folder, folder_num,
   extract, name_list, list_file, store, store_period)



# add zeros to wild type incidence
zero_mat = matrix(0, nrow = dim(inc_all[[1]])[1], ncol = time_period[variant =='W']$doy_start[1]-1)

for(i in 1:length(inc_all)){

  inc_all[[i]] = as.matrix(inc_all[[i]])

  if(grepl('wild',names(inc_all)[i])){

    inc_all[[i]] = inc_all[[i]][,extract_wild]
    inc_all[[i]] = cbind(zero_mat,inc_all[[i]])

  } else if(grepl('delta',names(inc_all)[i])){

    inc_all[[i]] = inc_all[[i]][,extract_delta]

  }

  colnames(inc_all[[i]]) = NULL
}

rm(i, zero_mat, extract_wild, extract_delta)

save(inc_all, file='output processed/20220818/inc_all.RData')
save(inc_all, file='output processed/20220521/inc_wild.RData')


# convert R missed into R notified, 
# R effective, R notified import, R missed import
param_W = list(time_period_matrix = time_period[variant =='W'])
param_D = list(time_period_matrix = time_period[variant =='D'])

# fixed disease transmission parameters
# generation interval
param_W$gen_mean_log = param_D$gen_mean_log = log(5.4)
param_W$gen_sd_log = param_D$gen_sd_log = 0.4

# incubation period
param_W$incub_mean_log = 1.63
param_W$incub_sd_log =  0.50

param_D$incub_mean_log = log(4)
param_D$incub_sd_log = sqrt(2*(log(4.4)-log(4)))

# load the distribution of duration from infection to isolation
generation_matrix_local_N_wild = matrixGenerationLocal(param_W, obs_data, notified = T)
generation_matrix_local_N_delta = matrixGenerationLocal(param_D, obs_data, notified = T)

# load the distribution of duration from infection to travel to isolation
generation_matrix_imported_M_wild = matrixGenerationImport(param_W, obs_data, notified = F)
generation_matrix_imported_M_delta = matrixGenerationImport(param_D, obs_data, notified = F)
generation_matrix_imported_N_wild = matrixGenerationImport(param_W, obs_data, notified = T)
generation_matrix_imported_N_delta = matrixGenerationImport(param_D, obs_data, notified = T)

# reformat data to extract the distribution of duration from infection to isolation in respective time periods
wild_period_start = param_W$time_period_matrix$doy_start-17
wild_period_end = (param_W$time_period_matrix$doy_start-17 +13)
delta_period_start = param_D$time_period_matrix$doy_start-456
delta_period_end = (param_D$time_period_matrix$doy_start-456 + 13)

period = 1:5
folder = c('wild neg binom', 'wild neg binom total', 'delta neg binom total')
prop_generation_at_large_local_N = list()
prop_generation_at_large_imported_N = list()
prop_generation_at_large_imported_M = list()

for(folder_num in 1:length(folder)){
  print(folder[folder_num])

  for(period_num in period){
    name_list = gsub(' ','_',paste(folder[folder_num], '_period_', period_num, sep=''))

    if(grepl('wild',folder[folder_num])){

      prop_generation_at_large_local_N[[name_list]] = sum(generation_matrix_local_N_wild[wild_period_start[period_num],wild_period_start[period_num]:wild_period_end[period_num]])
      prop_generation_at_large_imported_N[[name_list]] = sum(generation_matrix_imported_N_wild[wild_period_start[period_num],wild_period_start[period_num]:wild_period_end[period_num]])
      prop_generation_at_large_imported_M[[name_list]] = sum(generation_matrix_imported_M_wild[wild_period_start[period_num],wild_period_start[period_num]:wild_period_end[period_num]])

    } else if(grepl('delta',folder[folder_num])){

      prop_generation_at_large_local_N[[name_list]] = sum(generation_matrix_local_N_delta[delta_period_start[period_num],delta_period_start[period_num]:delta_period_end[period_num]])
      prop_generation_at_large_imported_N[[name_list]] = sum(generation_matrix_imported_N_delta[delta_period_start[period_num],delta_period_start[period_num]:delta_period_end[period_num]])
      prop_generation_at_large_imported_M[[name_list]] = sum(generation_matrix_imported_M_delta[delta_period_start[period_num],delta_period_start[period_num]:delta_period_end[period_num]])

    }
  }
}

# compute the R of notified cases and effective R
eigenNGM <- function(theta){

  sapply(1:nrow(theta), function(x){

    NGMComm = matrix(c((1-theta[[4]][x])*theta[[2]][x],  (1-theta[[3]][x])*theta[[8]][x],
                          theta[[4]][x]*theta[[2]][x],      theta[[3]][x]*theta[[8]][x]),
                     nrow = 2, ncol = 2, byrow = T)

    # missed offspring missed infector    missed offspring notified infector
    # notified offspring missed infector    notified offspring notified infector

    max(eigen(NGMComm)$values)

  })

}

for(i in 1:length(theta_all)){

  col_add_name =  paste(c('R_notified_period_', 'R_eff_period_', 'R_notified_imported_period_', 'R_missed_imported_period_'), unique(theta_all[[i]]$period), sep='')

  theta_all[[i]][[col_add_name[1]]] = prop_generation_at_large_local_N[[i]]*theta_all[[i]][[2]]
  theta_all[[i]][[col_add_name[2]]] = eigenNGM(theta_all[[i]])

  theta_all[[i]][[col_add_name[3]]] = prop_generation_at_large_imported_N[[i]]*theta_all[[i]][[2]]
  theta_all[[i]][[col_add_name[4]]] = prop_generation_at_large_imported_M[[i]]*theta_all[[i]][[2]]

  print(i)
}

save(theta_all, file='output processed/20220818/theta_all.RData')
save(theta_all, file='output processed/20220521/theta_wild.RData')

