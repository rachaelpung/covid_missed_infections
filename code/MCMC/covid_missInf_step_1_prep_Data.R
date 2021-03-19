library(mgcv)
library(readxl)
library(data.table)

source('codes/v5/MCMC/covid_missInf_functions_distributions.r')

# load incidence data
dataCase = read_excel('data/dataCase.xlsx',
                      col_types = c("text", "numeric", "text", "text", "date", "date", "date", "date", "date", 
                                    "date", "date", "date", "numeric", "text", "numeric", "text", "numeric", 
                                    "text", "date", "numeric", "text", "text", "date", "date"))
dataCase = data.table(dataCase)

# convert dates to DOY of year in dataCase
for(i in c(5:12,19,23,24)) {dataCase[,i] = as.numeric(as.Date(dataCase[[i]]) - as.Date('2019-12-31'))} 

# remove expunge cases
dataCase[SOURCE_OF_INFECTION == 2, rm:= 1]

# remove old infections (i.e. case with positive serology 8 days or more before notification date)
# with effect from 31 Jul, cases with S+ then C+ (>7 days but <=90 days) will not be assigned a case number
dataCase[CASE == 520, rm := 1]
dataCase[CASE ==  91, rm := 1]
dataCase[DATE_NOTIFICATION - DATE_SEROLOGY_POS >7 & 
         DATE_NOTIFICATION <= as.numeric(as.Date("2020-07-31") - as.Date('2019-12-31')), rm := 1]

# remove cases that did not acquire infection from the community (e.g. dorm ops cases, persons caring for COVID-19 patients)
dataCase[REMARKS == "COVID OPS TRANSMISSION", rm := 1]

# remove dorm cases notified between 7 Apr (start of lockdown) and 31 Oct  
dataCase[DATE_NOTIFICATION >=  as.numeric(as.Date("2020-04-07") - as.Date('2019-12-31')) &
         DATE_NOTIFICATION <=  as.numeric(as.Date("2020-10-31") - as.Date('2019-12-31')) & 
         DORM_DWELLER == 1, rm := 1]

dataCase = dataCase[is.na(rm)]
dataCase[,rm := NULL]

# load dormitory gazette date
dataGaz = read_excel('data/dataCase.xlsx', sheet = "Sheet3")
dataGaz = data.table(dataGaz)

# DOY when movement restrictions started for the whole day
dataGaz[,DATE_MOVEMENT_RESTRICTION_START := as.Date(DATE_MOVEMENT_RESTRICTION_START, "%d %b %Y") ]
dataGaz[,DATE_MOVEMENT_RESTRICTION_START := dataGaz[,2] + c(rep(1, times = 7), rep(0, times = 9)) ]
dataGaz[,DATE_MOVEMENT_RESTRICTION_START := as.numeric(dataGaz[[2]] - as.Date('2019-12-31')) ]

# DOY when movement restrictions lasted for the whole day 
dataGaz[,DATE_MOVEMENT_RESTRICTION_END := as.Date(DATE_MOVEMENT_RESTRICTION_END, "%d %b %Y") - 1 ]
dataGaz[,DATE_MOVEMENT_RESTRICTION_END := as.numeric(dataGaz[[3]] - as.Date('2019-12-31')) ]

# determine dates of movement restrictions in dorm cases with no AccessCode 
# case in gazetted dorms and notified before gazette end date
for(i in 1:dataGaz[,.N]){
  dataCase[DORM_DWELLER == 1 & 
           grepl(dataGaz$DORM_NAME[i], CLUSTER_REMARKS, fixed = T) & 
           DATE_NOTIFICATION >= dataGaz$DATE_MOVEMENT_RESTRICTION_START[i] &
           DATE_NOTIFICATION <= dataGaz$DATE_MOVEMENT_RESTRICTION_END[i],  
           DATE_MOVEMENT_RESTRICTION_START := dataGaz$DATE_MOVEMENT_RESTRICTION_START[i]]
}

# case in non gazetted dorms and notified between 21 Apr (workers not allowed to move out of their dorm) and 11 Aug (dorms declared to be cleared)
dataCase[DORM_DWELLER == 1 & is.na(DATE_MOVEMENT_RESTRICTION_START) &
         DATE_NOTIFICATION >= as.numeric(as.Date("2020-04-21") - as.Date('2019-12-31')) &
         DATE_NOTIFICATION <= as.numeric(as.Date('2020-08-31') - as.Date('2019-12-31')),    
         DATE_MOVEMENT_RESTRICTION_START := as.numeric(as.Date("2020-04-21") - as.Date('2019-12-31'))]

# case notified in Sep onwards set DATE_ACCESSCODE to NA, assume to all cases that were not quarantine would have exposure to the community  
dataCase[DATE_NOTIFICATION >=  as.numeric(as.Date("2020-09-01") - as.Date('2019-12-31')), DATE_ACCESSCODE := NA]

# case under HOME SHN assumed potential community exposure, else assumed no community exposure
dataCase[SHN_LOCATION == "HOME", DATE_SHN_OR_ISOLATED_FOR_TESTING := NA]

# find earliest DATE_ISOLATE_FROM_COMM (i.e. min of DATE_ADMISSION, DATE_CIF, DATE_SHN_OR_ISOLATED_FOR_TESTING,  
#                                                   DATE_QUARANTINE, DATE_ACCESSCODE, 
#                                                   DATE_REPORTED_ISOLATION, DATE_MOVEMENT_RESTRICTION_START)
colnames(dataCase)[c(7,8,10,11,23,24,25)]
dataCase[, DATE_ISOLATE_FROM_COMM := apply(dataCase[,c(7,8,10,11,23,24,25)], 1, min, na.rm = T)] 

# find which case whose date of isolation cannot be determined
# majority are dorm cases notified before 21 Apr
View(dataCase[DATE_ISOLATE_FROM_COMM == "Inf" & SOURCE_OF_INFECTION != 2,])

# must have isolation date for all acute infection
dataCase[DATE_ISOLATE_FROM_COMM == "Inf", DATE_ISOLATE_FROM_COMM := DATE_NOTIFICATION]

# remove local cases that were notified more than 14 days since isolation from community, assume to have acquired infection while isolated from community
# dataCase[DATE_NOTIFICATION - DATE_ISOLATE_FROM_COMM > 14 &
#          LOCAL_OVERSEAS == 0, rm := 1]

# remove cases with onset more than 21 days since isolation from community, assume to have acquired infection while isolated from community
# dataCase[!is.na(DATE_ONSET) & DATE_ONSET - DATE_ISOLATE_FROM_COMM > 21, rm := 1]

# dataCase = dataCase[is.na(rm)]
# dataCase[,rm := NULL]


# ----------------------------------------------------------------------------------------------------------------------
# tabulate daily incidence of linked and unlinked community cases who were infected before isolation
# 18 Jan 2020: DOY 18, 1 Jan 2021: DOY 367
dataCaseCommNotified = dataCase[LOCAL_OVERSEAS == 0 & DATE_ISOLATE_FROM_COMM <= 367]

dataInc = data.table(
  TIME =  19:367, 
  NOTIFIED_COMM_UNLINKED = hist(dataCaseCommNotified[SOURCE_OF_INFECTION == 0, DATE_ISOLATE_FROM_COMM], 
                                breaks = 18:367,plot = F)$counts,
  NOTIFIED_COMM_LINKED   = hist(dataCaseCommNotified[SOURCE_OF_INFECTION == 1, DATE_ISOLATE_FROM_COMM], 
                                breaks = 18:367, plot = F)$counts
)

dataInc = rbind(data.frame(TIME = 18, NOTIFIED_COMM_UNLINKED = 0, NOTIFIED_COMM_LINKED = 0), dataInc)

# save 
write.csv(dataInc, 'input/dataInc.csv', row.names = F)

# tabulate daily incidence of all imported cases and those with community exposure
dataIncImport   = data.table(
  TIME =  19:367, 
  NOTIFIED_IMPORT_COMM_EXPOSURE = hist(dataCase[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING), 
                                                DATE_ARRIVAL], breaks = 18:367, plot = F)$counts,
  NOTIFIED_IMPORT_ISOLATED = hist(dataCase[LOCAL_OVERSEAS == 1 & !is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING), 
                                           DATE_ARRIVAL], breaks = 18:367, plot = F)$counts
)

dataIncImport[1, NOTIFIED_IMPORT_COMM_EXPOSURE := dataIncImport[1, NOTIFIED_IMPORT_COMM_EXPOSURE] - 1]
dataIncImport = rbindlist(list(data.table(TIME = 18, NOTIFIED_IMPORT_COMM_EXPOSURE = 1, NOTIFIED_IMPORT_ISOLATED =0), 
                               dataIncImport), use.names = T, fill = T)

# fit smoothing spline for all notified imported cases with community exposure
dataIncImport[, SPLINE := gam(NOTIFIED_IMPORT_COMM_EXPOSURE~s(TIME), data = dataIncImport, family='poisson')$fitted.values]

# save 
write.csv(dataIncImport, 'input/dataIncImport.csv', row.names = F)

# tabulate imported cases with community exposure 
dataCaseImportNotified = dataCase[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING)]
dataCaseImportNotified = dataCaseImportNotified[, .(NOTIFIED, DATE_ARRIVAL, SOURCE_OF_INFECTION, 
                                                    LOCAL_OVERSEAS, DATE_ISOLATE_FROM_COMM)]
dataCaseImportNotified[, CASE_COUNT := 1]

# save 
write.csv(dataCaseImportNotified, 'input/dataCaseImportNotified.csv', row.names = F)

# generate distribution of generation interval over time
currentParam = list()
currentParam$mu_genInterval = 7.5
currentParam$sigma_genInterval = 3.4 
distGen = genIntDensity(currentParam)

# save 
write.csv(distGen, 'input/dist_Gen_Interval.csv', row.names = F)

# generate distribution of duration from infection to isolation
distIsolate = rbindlist(list(isolateDensity(dataCase[DATE_NOTIFICATION <= 60]),
                             isolateDensity(dataCase[DATE_NOTIFICATION >= 61 & DATE_NOTIFICATION <= 97]),
                             isolateDensity(dataCase[DATE_NOTIFICATION >= 98 & DATE_NOTIFICATION <= 170]),
                             isolateDensity(dataCase[DATE_NOTIFICATION >= 171 & DATE_NOTIFICATION <= 244]),
                             isolateDensity(dataCase[DATE_NOTIFICATION >= 245])),
                        use.names = T, fill = T)
  
distIsolate[, PERIOD := rep(1:5, each =22)] 

# save 
write.csv(distIsolate, 'input/dist_Isolate.csv', row.names = F)

# generate distribution of time of infection for imported cases
distTravel = travelDensity(dataCase)

# save 
write.csv(distTravel, 'input/dist_Travel.csv', row.names = F)

# generate distribution of time of arrival to isolation for imported cases
distArrivalIsolate = dataCase[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING), DATE_ISOLATE_FROM_COMM - DATE_ARRIVAL]
distArrivalIsolate = table(distArrivalIsolate)
distArrivalIsolate = prop.table(distArrivalIsolate)
distArrivalIsolate = data.table(distArrivalIsolate)
setnames(distArrivalIsolate, c('DAY_SINCE_ARRIVAL', 'PROB_ISOLATE_DAY_SINCE_ARRIVAL'))

distArrivalIsolate[, DAY_SINCE_ARRIVAL := as.numeric(DAY_SINCE_ARRIVAL)]
distArrivalIsolate = distArrivalIsolate[DAY_SINCE_ARRIVAL <= 21]
distArrivalIsolate[, PROB_ISOLATE_DAY_SINCE_ARRIVAL := PROB_ISOLATE_DAY_SINCE_ARRIVAL/sum(PROB_ISOLATE_DAY_SINCE_ARRIVAL)]

# save 
write.csv(distArrivalIsolate, 'input/dist_Arrival_Isolate.csv', row.names = F)
