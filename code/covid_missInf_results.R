source('codes/covid_missInf_load_library.R')

load('input/obs_data.RData')
load('output processed/20220521/inc_all.RData')  #0521 informative prior, 0818 non informative
load('output processed/20220521/theta_all.RData')

# load time period
time_period = read_excel('data/param.xlsx', sheet = 'time_period')
time_period = data.table(time_period)
setnames(time_period, c('period', 'date_start', 'date_end'))

time_period[, variant:=gsub('(.*)\\d', '\\1',period)]
time_period[, period:=gsub(variant,'',period), by=seq_len(nrow(time_period))]
time_period[, period:=as.numeric(period)]

time_period[, doy_start := as.numeric(as.Date(date_start)-as.Date('2019-12-31'))]
time_period[, doy_end := as.numeric(as.Date(date_end)-as.Date('2019-12-31'))]
time_period[, duration := doy_end-doy_start+1]

# tabulate cases by time period
period_start = c(time_period$doy_start[1], time_period$doy_start[1:5],
                 time_period$doy_start[1], time_period$doy_start[1:5],
                 time_period$doy_start[7], time_period$doy_start[7:10])
period_end = c(time_period$doy_end[5], time_period$doy_end[1:5],
               time_period$doy_end[5], time_period$doy_end[1:5],
               time_period$doy_end[10], time_period$doy_end[7:10])

case_table = data.table(period_variant = c(rep(c('wild_neg_binom_mod','wild_neg_binom_total'), each=6),
                                           rep(c('delta_neg_binom_total_mod'), each=5)), 
                        period_start = period_start, period_end = period_end)
case_table[, import_N := sum(obs_data$daily_arrival_import_N_shn[period_start:period_end]), by=seq_len(nrow(case_table))]
case_table[, import_N_no_shn := sum(obs_data$daily_arrival_import_N_no_shn[period_start:period_end]), by=seq_len(nrow(case_table))]
case_table[, local_N_linked := sum(obs_data$daily_local_N_linked[period_start:period_end]), by=seq_len(nrow(case_table))]
case_table[, local_N_unlinked := sum(obs_data$daily_local_N_unlinked[period_start:period_end]), by=seq_len(nrow(case_table))]

case_table[, local_N_death_by_inf := sum(obs_data$daily_local_N_death_by_inf_doy[period_start:period_end]), by=seq_len(nrow(case_table))]
case_table[, local_N_icu_by_inf := sum(obs_data$daily_local_N_icu_by_inf_doy[period_start:period_end]), by=seq_len(nrow(case_table))]
case_table[, local_N_death := sum(obs_data$daily_local_N_death[period_start:period_end]), by=seq_len(nrow(case_table))]
case_table[, local_N_icu := sum(obs_data$daily_local_N_icu[period_start:period_end]), by=seq_len(nrow(case_table))]



# wild type missed
case_table[, period_missed_start := period_start]
case_table[, period_missed_end := period_end]
case_table[period_variant == 'delta_neg_binom_total_mod', period_missed_start := period_missed_start-456] 
case_table[period_variant == 'delta_neg_binom_total_mod', period_missed_end := period_missed_end-456]


for(v in case_table[,unique(period_variant)]){
  
  subset_inc_all = inc_all[grep(v, names(inc_all))]
  
  # daily_local_M_unlinked + daily_local_M_linked
  case_table[period_variant == v, local_M_med := quantile(apply(subset_inc_all[[3]][,period_missed_start:period_missed_end],1,sum)+
                                                          apply(subset_inc_all[[4]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.5)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  case_table[period_variant == v, local_M_lwr := quantile(apply(subset_inc_all[[3]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[4]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.025)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  case_table[period_variant == v, local_M_upp := quantile(apply(subset_inc_all[[3]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[4]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.975)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  # daily.local_M.unlinked + daily.local_M.linked + daily_local_N_unlinked + daily_local_N_linked
  case_table[period_variant == v, local_T_med := quantile(apply(subset_inc_all[[3]][,period_missed_start:period_missed_end],1,sum)+
                                                          apply(subset_inc_all[[4]][,period_missed_start:period_missed_end],1,sum)+
                                                          apply(subset_inc_all[[1]][,period_missed_start:period_missed_end],1,sum)+
                                                          apply(subset_inc_all[[2]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.5)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  case_table[period_variant == v, local_T_lwr := quantile(apply(subset_inc_all[[3]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[4]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[1]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[2]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.025)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  case_table[period_variant == v, local_T_upp := quantile(apply(subset_inc_all[[3]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[4]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[1]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[2]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.975)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  case_table[period_variant == v, local_N_med := quantile(apply(subset_inc_all[[1]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[2]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.5)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  case_table[period_variant == v, local_N_lwr := quantile(apply(subset_inc_all[[1]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[2]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.025)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  case_table[period_variant == v, local_N_upp := quantile(apply(subset_inc_all[[1]][,period_missed_start:period_missed_end],1,sum)+
                                                            apply(subset_inc_all[[2]][,period_missed_start:period_missed_end],1,sum),
                                                          probs = c(0.975)),
             by=seq_len(nrow(case_table[period_variant == v,]))]
  
  print(v)
  
}




case_table[,CFR_med := local_N_death_by_inf*100/(local_N_med)]
case_table[,CFR_lwr := local_N_death_by_inf*100/(local_N_lwr)]
case_table[,CFR_upp := local_N_death_by_inf*100/(local_N_upp)]

case_table[,IFR_med := local_N_death_by_inf*100/local_T_med]
case_table[,IFR_lwr := local_N_death_by_inf*100/local_T_lwr]
case_table[,IFR_upp := local_N_death_by_inf*100/local_T_upp]

case_table[,CIR_med := local_N_icu_by_inf*100/(local_N_med)]
case_table[,CIR_lwr := local_N_icu_by_inf*100/(local_N_lwr)]
case_table[,CIR_upp := local_N_icu_by_inf*100/(local_N_upp)]

case_table[,IIR_med := local_N_icu_by_inf*100/local_T_med]
case_table[,IIR_lwr := local_N_icu_by_inf*100/local_T_lwr]
case_table[,IIR_upp := local_N_icu_by_inf*100/local_T_upp]


View(case_table)


# incidence for wild-type SARS-CoV-2

# notified cases before lockdown
sum(obs_data$daily_local_N[1:97])

# using linked and unlinked cases for model fitting

# missed cases before lockdown
quantile(apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,61:97],1,sum)+
           apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))

# daily missed cases one week before lockdown
quantile(apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,91:97],1,sum) + 
           apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,91:97],1,sum), 
         probs = c(0.5,0.025,0.975))/7

# daily missed cases before lockdown
quantile(apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,61:97],1,sum)+
         apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))/37

# daily missed cases during partial lockdown and phase 1 reopening
quantile(apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,98:170],1,sum)+
        apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,98:170],1,sum),
        probs=c(0.5,0.025, 0.975))/73

# daily missed cases after partial lockdown 
quantile(apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,171:366],1,sum)+
           apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,171:366],1,sum),
         probs=c(0.5,0.025, 0.975))/196

# missed cases in 2020
quantile(apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,1:366],1,sum)+
           apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,1:366],1,sum),
         probs=c(0.5,0.025, 0.975))

# missed imported cases 
quantile(sum(obs_data$daily_arrival_import_N_ma[18:60])*theta_all$wild_neg_binom_period_1$prop_import_M_period_1,
         probs=c(0.5,0.025, 0.975))


# using cases with no info on linkage for model fitting

# missed cases before lockdown
quantile(apply(inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,61:97],1,sum)+
            apply(inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))

# daily missed cases before lockdown
quantile(apply(inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,61:97],1,sum)+
            apply(inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,61:97],1,sum),
         probs=c(0.5,0.025, 0.975))/37

# daily missed cases one week before lockdown
quantile(apply(inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,91:97],1,sum) + 
           apply(inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,91:97],1,sum), 
         probs = c(0.5,0.025,0.975))/7

# daily missed cases during partial lockdown and phase 1 reopening
quantile(apply(inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,98:170],1,sum)+
           apply(inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,98:170],1,sum),
         probs=c(0.5,0.025, 0.975))/73

# difference between models in estimating missed infections during lockdown period
quantile(((apply(inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,98:170],1,sum)+
           apply(inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,98:170],1,sum))/73 )/
           ((apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,98:170],1,sum)+
              apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,98:170],1,sum))/73),
         probs=c(0.5,0.025, 0.975))

# missed cases in 2020
quantile(apply(inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,1:366],1,sum)+
            apply(inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,1:366],1,sum),
         probs=c(0.5,0.025, 0.975))

# proportion of missed cases in 2020
missed_2020 = (apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,1:366],1,sum)+
                 apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,1:366],1,sum))
total_2020 = (apply(inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,1:366],1,sum)+
                apply(inc_all$wild_neg_binom_mod_daily_local_M_linked[,1:366],1,sum)+
                apply(inc_all$wild_neg_binom_mod_daily_local_N_unlinked[,1:366],1,sum)+
                apply(inc_all$wild_neg_binom_mod_daily_local_N_linked[,1:366],1,sum))


quantile(missed_2020/total_2020,
         probs=c(0.5,0.025, 0.975))

# proportion of missed cases in 2021
missed_2021 = (apply(inc_all$delta_neg_binom_total_mod_daily_local_M_unlinked[,1:140],1,sum)+
                 apply(inc_all$delta_neg_binom_total_mod_daily_local_M_linked[,1:140],1,sum))
total_2021 = (apply(inc_all$delta_neg_binom_total_mod_daily_local_M_unlinked[,1:140],1,sum)+
                apply(inc_all$delta_neg_binom_total_mod_daily_local_M_linked[,1:140],1,sum)+
                apply(inc_all$delta_neg_binom_total_mod_daily_local_N_unlinked[,1:140],1,sum)+
                apply(inc_all$delta_neg_binom_total_mod_daily_local_N_linked[,1:140],1,sum))


quantile(missed_2021/total_2021,
         probs=c(0.5,0.025, 0.975))


# estimating parameters

# start of the pandemic
apply(theta_all$wild_neg_binom_period_1, 2, quantile, probs = c(0.5,0.025,0.975))
apply(theta_all$wild_neg_binom_period_1[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs_data$daily_arrival_import_N[18:60])


# March to April
apply(theta_all$wild_neg_binom_period_2, 2, quantile, probs = c(0.5,0.025,0.975))

apply(theta_all$wild_neg_binom_period_2[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs_data$daily_arrival_import_N[61:97])

# Partial lockdown
apply(theta_all$wild_neg_binom_period_3, 2, quantile, probs = c(0.5,0.025,0.975))

# reopening
apply(theta_all$wild_neg_binom_period_4, 2, quantile, probs = c(0.5,0.025,0.975))

apply(theta_all$wild_neg_binom_period_4[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs_data$daily_arrival_import_N[171:194])

apply(theta_all$wild_neg_binom_total_period_4, 2, quantile, probs = c(0.5,0.025,0.975))

apply(theta_all$wild_neg_binom_period_5, 2, quantile, probs = c(0.5,0.025,0.975))
apply(theta_all$wild_neg_binom_period_5[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs_data$daily_arrival_import_N[195:366])


# incidence for delta-variant SARS-CoV-2
sum(obs_data$daily_local_N[457:501])/45

# early importation of delta cases
apply(theta_all$delta_neg_binom_total_period_1[,1], 2, quantile, probs = c(0.5,0.025,0.975)) * mean(obs_data$daily_arrival_import_N[457:501])
apply(theta_all$delta_neg_binom_total_period_1, 2, quantile, probs = c(0.5,0.025,0.975))

# missed cases in early delta outbreak
quantile(apply(inc_all$delta_neg_binom_total_mod_daily_local_M_unlinked[,1:45],1,sum)+
           apply(inc_all$delta_neg_binom_total_mod_daily_local_M_linked[,1:45],1,sum),
         probs=c(0.5,0.025, 0.975))

# tighten and relaxation
apply(theta_all$delta_neg_binom_total_period_2, 2, quantile, probs = c(0.5,0.025,0.975))


# surge phase
apply(theta_all$delta_neg_binom_total_period_3, 2, quantile, probs = c(0.5,0.025,0.975))

# decline phase
apply(theta_all$delta_neg_binom_total_period_4, 2, quantile, probs = c(0.5,0.025,0.975))

quantile(apply(inc_all$delta_neg_binom_total_mod_daily_local_M_unlinked[,109:140],1,sum)+
           apply(inc_all$delta_neg_binom_total_mod_daily_local_M_linked[,109:140],1,sum),
         probs=c(0.5,0.025, 0.975))/32
