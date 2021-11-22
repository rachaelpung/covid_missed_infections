source('codes/v6/covid_missInf_load_library.R')

# load data
data.case = read_excel('data/dataCase.xlsx',
                      col_types = c("text", "numeric", "text", "text", "text", "numeric",  "date", "date", 
                                    "date", "date", "date", "date", "date", "date", "numeric", "text", 
                                    "numeric", "text", "numeric", "text", "date", "numeric", "text", "text", 
                                    "date", "date", "text", "text", "text", "text", "date", "text", "text", 
                                    "date", "date", "date", "date"))

# covert to data table
data.case = data.table(data.case)

# format column names
colnames(data.case) = gsub("_", '.', colnames(data.case))
colnames(data.case) = tolower(colnames(data.case))
col.with.dates = grep('date', colnames(data.case))

# convert date columns into day of year (i.e. doy)
for(c in col.with.dates){
  
  col.name = colnames(data.case)[c]
  col.name = gsub("date", 'doy', col.name)
  
  data.case$new = as.numeric(as.Date(data.case[[c]]) - as.Date('2019-12-31'))
  colnames(data.case)[length(data.case)] = col.name

}


# remove expunge cases
data.case[source.of.infection==2, rm:=1]

# remove old infections (i.e. case with positive serology 8 days or more before notification date)
# with effect from 31 Jul, cases with S+ then C+ (>7 days but <=90 days) will not be assigned a case number
data.case[case == 520, rm := 1]
data.case[case ==  91, rm := 1]
data.case[doy.notification - doy.serology.pos >7 & 
          doy.notification <= as.numeric(as.Date("2020-07-31") - as.Date('2019-12-31')), rm := 1]

# remove cases that did not acquire infection from the community (e.g. dorm ops cases, persons caring for COVID-19 patients)
data.case[remarks == "COVID OPS TRANSMISSION", rm := 1]

# remove dorm cases notified between 7 Apr (start of lockdown) and 31 Oct  
data.case[doy.notification >=  as.numeric(as.Date("2020-04-07") - as.Date('2019-12-31')) &
          doy.notification <=  as.numeric(as.Date("2020-10-31") - as.Date('2019-12-31')) & 
          dorm.dweller == 1, rm := 1]

data.case = data.case[is.na(rm)]
data.case[,rm := NULL]

# load dormitory gazette date
data.gazette = read_excel('data/dataCase.xlsx', sheet = "Sheet3")
data.gazette = data.table(data.gazette)

# format column names
colnames(data.gazette) = gsub("_", '.', colnames(data.gazette))
colnames(data.gazette) = tolower(colnames(data.gazette))

# convert to doy when movement restrictions started for the whole day
data.gazette[,doy.restrict.start := as.Date(date.restrict.start, "%d %b %Y") ]
data.gazette[1:7,doy.restrict.start := doy.restrict.start + 1]
data.gazette[,doy.restrict.start := as.numeric(doy.restrict.start - as.Date('2019-12-31')) ]

# convert to doy when movement restrictions lasted for the whole day 
data.gazette[,doy.restrict.end := as.Date(date.restrict.end, "%d %b %Y") - 1 ]
data.gazette[,doy.restrict.end := as.numeric(doy.restrict.end - as.Date('2019-12-31')) ]

# assign doy movement restrictions to dorm cases with no AccessCode 
# case in gazetted dorms and notified before gazette end date
for(i in 1:data.gazette[,.N]){
  data.case[dorm.dweller == 1 & 
            grepl(data.gazette$dorm.name[i], cluster.remarks, fixed = T) & 
            doy.notification >= data.gazette$doy.restrict.start[i] &
            doy.notification <= data.gazette$doy.restrict.end[i],  
            doy.restrict.start := data.gazette$doy.restrict.start[i]]
}

# case in non gazetted dorms and notified between 21 Apr (workers not allowed to move out of their dorm) and 11 Aug (dorms declared to be cleared)
data.case[dorm.dweller==1 & is.na(doy.restrict.start) &
          doy.notification >= as.numeric(as.Date("2020-04-21") - as.Date('2019-12-31')) &
          doy.notification <= as.numeric(as.Date('2020-08-31') - as.Date('2019-12-31')),
          doy.restrict.start := as.numeric(as.Date("2020-04-21") - as.Date('2019-12-31'))]

# if case notified in Sep onwards set doy.accesscode to NA, assume that these cases that were not quarantine would have exposure to the community  
data.case[doy.notification >=  as.numeric(as.Date("2020-09-01") - as.Date('2019-12-31')), doy.accesscode := NA]

# if case under HOME SHN prior to issuance of tracker (10 Aug 2020) assumed potential community exposure 
# else assumed no community exposure
data.case[shn.location == "HOME" &
          doy.arrival <= as.numeric(as.Date('2020-08-10') - as.Date('2019-12-31'))
          , doy.shn.or.isolated.for.testing := NA]

# find earliest DATE_ISOLATE_FROM_COMM (i.e. min of doy.admission, doy.cif, doy.shn.or.isolated.for.testing,  
#                                                   doy.quarantine, doy.accesscode, 
#                                                   doy.reported.isolation, doy.movement.restriction.start)
extract.col = c('doy.admission', 'doy.cif', 'doy.shn.or.isolated.for.testing','doy.quarantine', 
                'doy.accesscode', 'doy.reported.isolation', 'doy.restrict.start')
data.case[, doy.isolate.from.comm := apply(data.case[,extract.col, with=F], 1, min, na.rm = T)] 

# find which case whose date of isolation cannot be determined
# majority are dorm cases notified before 21 Apr
View(data.case[doy.isolate.from.comm == "Inf",])

# must have isolation date for all acute infection
data.case[doy.isolate.from.comm == "Inf" & doy.notification<=366, doy.isolate.from.comm := doy.notification]
data.case[doy.isolate.from.comm == "Inf" & local.overseas==0, doy.isolate.from.comm := doy.notification]


# set-up storage
t_max = max(data.case$doy.isolate.from.comm[is.finite(data.case$doy.isolate.from.comm)])
obs.data = list()
obs.data$doy = 1:t_max

# tabulate local cases by day of isolation
obs.data$daily.local.N = c(0,hist(data.case[local.overseas==0, doy.isolate.from.comm], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily.local.N.linked = c(0,hist(data.case[local.overseas==0 & source.of.infection==1, doy.isolate.from.comm], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily.local.N.unlinked = c(0,hist(data.case[local.overseas==0 & source.of.infection==0, doy.isolate.from.comm], breaks = obs.data$doy,plot = F)$counts)

# assume only imported cases with community exposure can generate local infections
# obs.data$daily.arrival.import.N = c(0,hist(data.case[local.overseas==1 & is.na(doy.shn.or.isolated.for.testing), doy.arrival], breaks = obs.data$doy,plot = F)$counts)
# obs.data$daily.arrival.import.N = gam(obs.data$daily.arrival.import.N~s(obs.data$doy), family='poisson')$fitted.values

obs.data$daily.arrival.import.N = c(0,hist(data.case[local.overseas==1, doy.arrival], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily.arrival.import.N.no.shn = c(0,hist(data.case[local.overseas==1 & is.na(doy.shn.or.isolated.for.testing), doy.arrival], breaks = obs.data$doy,plot = F)$counts)
obs.data$daily.arrival.import.N.shn = c(0,hist(data.case[local.overseas==1 & !is.na(doy.shn.or.isolated.for.testing), doy.arrival], breaks = obs.data$doy,plot = F)$counts)



obs.data$daily.arrival.import.N.ma = ma(obs.data$daily.arrival.import.N,7)
obs.data$daily.arrival.import.N.ma[is.na(obs.data$daily.arrival.import.N.ma )] = 0
obs.data$daily.arrival.import.N.ma = as.numeric(obs.data$daily.arrival.import.N.ma)


# line list of imported cases
obs.data$import.N = data.case[local.overseas==1, .(case,doy.notification,doy.arrival,doy.onset,doy.isolate.from.comm)]

# line list of local cases
obs.data$local.N = data.case[local.overseas==0, .(case,doy.notification,doy.onset,doy.isolate.from.comm)]

# tabulate local icu cases by day of isolation
obs.data$daily.local.N.icu = c(0,hist(data.case[local.overseas==0 & !is.na(date.icu), doy.isolate.from.comm], breaks = obs.data$doy,plot = F)$counts)

# tabulate local death cases by day of isolation
obs.data$daily.local.N.death = c(0,hist(data.case[local.overseas==0 & !is.na(date.death), doy.isolate.from.comm], breaks = obs.data$doy,plot = F)$counts)


save(obs.data, file = 'input/obs.data.RData')

