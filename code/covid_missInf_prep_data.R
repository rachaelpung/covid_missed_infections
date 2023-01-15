source('codes/covid_missInf_load_library.R')

# load data
data_case = read_excel('data/dataCase.xlsx',
                      col_types = c("text", "numeric", "text", "text", "text", "numeric",  "date", "date", 
                                    "date", "date", "date", "date", "date", "date", "numeric", "text", 
                                    "numeric", "text", "numeric", "text", "date", "numeric", "text", "text", 
                                    "date", "date", "text", "text", "text", "text", "date", "text", "text", 
                                    "date", "date", "date", "date"))

# covert to data table
data_case = data.table(data_case)

# format column names
colnames(data_case) = gsub("_", '.', colnames(data_case))
colnames(data_case) = tolower(colnames(data_case))
col_with_dates = grep('date', colnames(data_case))

# convert date columns into day of year (i.e. doy)
for(c in col_with_dates){
  
  col_name = colnames(data_case)[c]
  col_name = gsub("date", 'doy', col_name)
  
  data_case$new = as.numeric(as.Date(data_case[[c]]) - as.Date('2019-12-31'))
  colnames(data_case)[length(data_case)] = col_name

}


# remove expunge cases
data_case[source_of_infection==2, rm:=1]

# remove old infections (i.e. case with positive serology 8 days or more before notification date)
# with effect from 31 Jul, cases with S+ then C+ (>7 days but <=90 days) will not be assigned a case number
data_case[case == 520, rm := 1]
data_case[case ==  91, rm := 1]
data_case[doy_notification - doy_serology_pos >7 & 
          doy_notification <= as.numeric(as.Date("2020-07-31") - as.Date('2019-12-31')), rm := 1]

# remove cases that did not acquire infection from the community (e.g. dorm ops cases, persons caring for COVID-19 patients)
data_case[remarks == "COVID OPS TRANSMISSION", rm := 1]

# remove dorm cases notified between 7 Apr (start of lockdown) and 31 Oct  
data_case[doy_notification >=  as.numeric(as.Date("2020-04-07") - as.Date('2019-12-31')) &
          doy_notification <=  as.numeric(as.Date("2020-10-31") - as.Date('2019-12-31')) & 
          dorm_dweller == 1, rm := 1]

data_case = data_case[is.na(rm)]
data_case[,rm := NULL]

# load dormitory gazette date
data_gazette = read_excel('data/dataCase.xlsx', sheet = "Sheet3")
data_gazette = data.table(data_gazette)

# format column names
colnames(data_gazette) = gsub("_", '.', colnames(data_gazette))
colnames(data_gazette) = tolower(colnames(data_gazette))

# convert to doy when movement restrictions started for the whole day
data_gazette[,doy_restrict_start := as.Date(date_restrict_start, "%d %b %Y") ]
data_gazette[1:7,doy_restrict_start := doy_restrict_start + 1]
data_gazette[,doy_restrict_start := as.numeric(doy_restrict_start - as.Date('2019-12-31')) ]

# convert to doy when movement restrictions lasted for the whole day 
data_gazette[,doy_restrict_end := as.Date(date_restrict_end, "%d %b %Y") - 1 ]
data_gazette[,doy_restrict_end := as.numeric(doy_restrict_end - as.Date('2019-12-31')) ]

# assign doy movement restrictions to dorm cases with no AccessCode 
# case in gazetted dorms and notified before gazette end date
for(i in 1:data_gazette[,.N]){
  data_case[dorm_dweller == 1 & 
            grepl(data_gazette$dorm.name[i], cluster.remarks, fixed = T) & 
            doy_notification >= data_gazette$doy_restrict_start[i] &
            doy_notification <= data_gazette$doy_restrict_end[i],  
            doy_restrict_start := data_gazette$doy_restrict_start[i]]
}

# case in non gazetted dorms and notified between 21 Apr (workers not allowed to move out of their dorm) and 11 Aug (dorms declared to be cleared)
data_case[dorm_dweller==1 & is.na(doy_restrict_start) &
          doy_notification >= as.numeric(as.Date("2020-04-21") - as.Date('2019-12-31')) &
          doy_notification <= as.numeric(as.Date('2020-08-31') - as.Date('2019-12-31')),
          doy_restrict_start := as.numeric(as.Date("2020-04-21") - as.Date('2019-12-31'))]

# if case notified in Sep onwards set doy_accesscode to NA, assume that these cases that were not quarantine would have exposure to the community  
data_case[doy_notification >=  as.numeric(as.Date("2020-09-01") - as.Date('2019-12-31')), doy_accesscode := NA]

# if case under HOME SHN prior to issuance of tracker (10 Aug 2020) assumed potential community exposure 
# else assumed no community exposure
data_case[shn_location == "HOME" &
          doy_arrival <= as.numeric(as.Date('2020-08-10') - as.Date('2019-12-31'))
          , doy_shn_or_isolated_for_testing := NA]

# find earliest DATE_ISOLATE_FROM_COMM (i.e. min of doy_admission, doy_cif, doy_shn_or_isolated_for_testing,  
#                                                   doy_quarantine, doy_accesscode, 
#                                                   doy_reported_isolation, doy.movement.restriction.start)
extract.col = c('doy_admission', 'doy_cif', 'doy_shn_or_isolated_for_testing','doy_quarantine', 
                'doy_accesscode', 'doy_reported_isolation', 'doy_restrict_start')
data_case[, doy_isolate_from_comm := apply(data_case[,extract.col, with=F], 1, min, na.rm = T)] 

# find which case whose date of isolation cannot be determined
# majority are dorm cases notified before 21 Apr
View(data_case[doy_isolate_from_comm == "Inf",])

# must have isolation date for all acute infection
data_case[doy_isolate_from_comm == "Inf" & doy_notification<=366, doy_isolate_from_comm := doy_notification]
data_case[doy_isolate_from_comm == "Inf" & local.overseas==0, doy_isolate_from_comm := doy_notification]


# set-up storage
t_max = max(data_case$doy_isolate_from_comm[is.finite(data_case$doy_isolate_from_comm)])
obs.data = list()
obs.data$doy = 1:t_max

# tabulate local cases by day of isolation
obs.data$daily_local_N = c(0,hist(data_case[local.overseas==0, doy_isolate_from_comm], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily_local_N_linked = c(0,hist(data_case[local.overseas==0 & source_of_infection==1, doy_isolate_from_comm], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily_local_N_unlinked = c(0,hist(data_case[local.overseas==0 & source_of_infection==0, doy_isolate_from_comm], breaks = obs.data$doy,plot = F)$counts)

# assume only imported cases with community exposure can generate local infections
# obs.data$daily_arrival_import_N = c(0,hist(data_case[local.overseas==1 & is.na(doy_shn_or_isolated_for_testing), doy_arrival], breaks = obs.data$doy,plot = F)$counts)
# obs.data$daily_arrival_import_N = gam(obs.data$daily_arrival_import_N~s(obs.data$doy), family='poisson')$fitted.values

obs.data$daily_arrival_import_N = c(0,hist(data_case[local.overseas==1, doy_arrival], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily_arrival_import_N_no_shn = c(0,hist(data_case[local.overseas==1 & is.na(doy_shn_or_isolated_for_testing), doy_arrival], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily_arrival_import_N_shn = c(0,hist(data_case[local.overseas==1 & !is.na(doy_shn_or_isolated_for_testing), doy_arrival], breaks = obs.data$doy,plot = F)$counts)



obs.data$daily_arrival_import_N_ma = ma(obs.data$daily_arrival_import_N,7)
obs.data$daily_arrival_import_N_ma[is.na(obs.data$daily_arrival_import_N_ma )] = 0
obs.data$daily_arrival_import_N_ma = as.numeric(obs.data$daily_arrival_import_N_ma)


# line list of imported cases
obs.data$import_N = data_case[local.overseas==1, .(case,doy_notification,doy_arrival,doy.onset,doy_isolate_from_comm)]

# line list of local cases
obs.data$local_N = data_case[local.overseas==0, .(case,doy_notification,doy.onset,doy_isolate_from_comm)]

# tabulate local icu cases by day of isolation
obs.data$daily_local_N_icu = c(0,hist(data_case[local.overseas==0 & !is.na(date.icu), doy_isolate_from_comm], breaks = obs.data$doy,plot = F)$counts)

# tabulate local death cases by day of isolation
obs.data$daily_local_N_death = c(0,hist(data_case[local.overseas==0 & !is.na(date.death), doy_isolate_from_comm], breaks = obs.data$doy,plot = F)$counts)


save(obs.data, file = 'input/obs.data.RData')

