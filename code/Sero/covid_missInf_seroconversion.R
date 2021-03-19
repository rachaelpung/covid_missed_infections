library(mgcv)
library(readr)
library(readxl)
library(data.table)

load('output processed/20210223/outputTotalChain_10000_thinned_iter_period_all.rdata')
outputTotal = copy(outputTotalChain_period_all)

# PCR and sero testing in preschool staff
startTest = as.numeric(as.Date("2020-05-15") - as.Date('2019-12-31')) 
endTest   = as.numeric(as.Date("2020-05-29") - as.Date('2019-12-31'))

startTime = startTest - 30
endTime   = endTest - 1

# extract missed infections on and since 31 Mar
colStart = which(outputTotal$time == startTime)
colEnd   = which(outputTotal$time == endTime) # assumes missed infections at least 1 days ago can be PCR and sero pos
dataIncMissed = outputTotal$IncMissedSim_DOI[,colStart:colEnd]
  
# reformat 
dataIncMissed = data.table(TIME = rep(startTime:endTime, times = 10000),
                           MISSED = as.vector(t(dataIncMissed)),
                           ITER = rep(1:10000, each = 44))

# prob of sero conversion
probSeroConvert = 0.961

# pdf incubation period
distIncub = data.table(TIME_SINCE_INFECTION = 1:21)
distIncub[, PROB_ONSET := dlnorm(TIME_SINCE_INFECTION, meanlog = 1.63, sdlog = 0.50, log = FALSE)]
distIncub[, PROB_ONSET := PROB_ONSET/sum(PROB_ONSET)]

# cdf of igg detection probability since onset given sero given seroconversion
distSero =  read_csv('input/dist_Sero_Positivity.csv')
spline = gam(PERCENTAGE_POSITIVE~s(TIME), data = distSero)
predict_spline = predict(spline,newdata=data.frame(TIME=0:60))/100
distSero = data.table(TIME_SINCE_ONSET = 0:60, PROB_SERO_DETECTION = as.vector(predict_spline))
distSero[PROB_SERO_DETECTION>1, PROB_SERO_DETECTION := 1]
# distSero[PROB_SERO_DETECTION<0, PROB_SERO_DETECTION := 0]

# convert time since onset to time since infection
distSero = distSero[rep(seq_len(.N), times = 21),]
distSero[, TIME_INFECTION_TO_ONSET := rep(1:21, each = .N/21)]
distSero[, TIME_SINCE_INFECTION := TIME_INFECTION_TO_ONSET+TIME_SINCE_ONSET]
distSero[, PROB_ONSET := rep(distIncub[,PROB_ONSET], each = .N/21)]
distSero[, PROB_SERO_DETECTION := PROB_SERO_DETECTION*PROB_ONSET]
distSero = distSero[, .(TIME_SINCE_INFECTION, PROB_SERO_DETECTION)]
distSero = distSero[, .(sum(PROB_SERO_DETECTION)), by = 'TIME_SINCE_INFECTION']
setnames(distSero, c('TIME_SINCE_INFECTION', 'PROB_SERO_DETECTION'))
distSero = distSero[TIME_SINCE_INFECTION <=61]

# distribution of probability of detecting pcr pos curves since time of infection
distPCR = data.table(read_csv('input/dist_PCR_Positivity.csv'))
spline = gam(PERCENTAGE_POSITIVE~s(TIME), data = distPCR)
predict_spline = predict(spline,newdata=data.frame(TIME=0:60))/100
distPCR = data.table(TIME_SINCE_ONSET = 0:60, PROB_PCR_DECTECTION = as.vector(predict_spline))
distPCR[PROB_PCR_DECTECTION>1, PROB_PCR_DECTECTION := 1]
distPCR[PROB_PCR_DECTECTION<0, PROB_PCR_DECTECTION := 0]

# convert time since onset to time since infection
distPCR = distPCR[rep(seq_len(.N), times = 21),]
distPCR[, TIME_INFECTION_TO_ONSET := rep(1:21, each = .N/21)]
distPCR[, TIME_SINCE_INFECTION := TIME_INFECTION_TO_ONSET+TIME_SINCE_ONSET]
distPCR[, PROB_ONSET := rep(distIncub[,PROB_ONSET], each = .N/21)]
distPCR[, PROB_PCR_DECTECTION := PROB_PCR_DECTECTION*PROB_ONSET]
distPCR = distPCR[, .(TIME_SINCE_INFECTION, PROB_PCR_DECTECTION)]
distPCR = distPCR[, .(sum(PROB_PCR_DECTECTION)), by = 'TIME_SINCE_INFECTION']
setnames(distPCR, c('TIME_SINCE_INFECTION', 'PROB_PCR_DETECTION'))
distPCR = distPCR[TIME_SINCE_INFECTION <=61]

# read in swab results
HPB_Masterlist = data.table(read_excel("Z:/Wuhan pneumonia/Lab results/Labs/HPB results/HPB_Masterlist.xlsx"))
HPB_Masterlist = HPB_Masterlist[`REGISTRATION TIME` <= '2020-05-31', ]
HPB_Masterlist = HPB_Masterlist[`SITE OF TESTING` != 'WESTLITE JUNIPER']

HPB_Masterlist$date = as.Date(as.character(HPB_Masterlist$`DATE OF TESTING`), "%Y%m%d")
HPB_Masterlist = HPB_Masterlist[c(which(duplicated(HPB_Masterlist[,c(1,5)]) == FALSE)),]

# probability of being tested on a particular day from 15-29 May
probTestDay = data.frame(table(HPB_Masterlist$date))$Freq/41852 # 1/(29-15+1)

# for every missed infection infected on day startTime = 91, 
# how many will be PCR pos & sero pos on startTest = 136 till endTest = 150
dataIncMissed = dataIncMissed[rep(seq_len(.N), each = length(startTest:endTest) )]
dataIncMissed[, TIME_TEST := rep(startTest:endTest, times = .N/15)]
dataIncMissed[, PROB_TIME_TEST := rep(probTestDay, times = .N/15)]

dataIncMissed[, TIME_TEST_SINCE_INFECTION := TIME_TEST - TIME]
dataIncMissed[, PROB_SERO_CONVERT := probSeroConvert]
dataIncMissed = merge(dataIncMissed, distSero, by.x = 'TIME_TEST_SINCE_INFECTION', by.y = 'TIME_SINCE_INFECTION')
dataIncMissed = merge(dataIncMissed, distPCR[,.(TIME_SINCE_INFECTION, PROB_PCR_DETECTION)], by.x = 'TIME_TEST_SINCE_INFECTION', by.y = 'TIME_SINCE_INFECTION')
dataIncMissed[, PROB_IDENTIFIED_TEST_PERIOD :=  PROB_TIME_TEST*PROB_SERO_CONVERT*PROB_SERO_DETECTION*PROB_PCR_DETECTION]

dataIncMissed[, MISSED_CASES_IDENTIFIED_TEST_PERIOD := MISSED*PROB_IDENTIFIED_TEST_PERIOD]
dataIncMissed[, MISSED_CASES_IDENTIFIED_TEST_PERIOD_INC_RATE := MISSED_CASES_IDENTIFIED_TEST_PERIOD/(5.381)]


summary(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD), by = ITER]$V1)
summary(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD_INC_RATE), by = ITER]$V1)
quantile(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD_INC_RATE), by = ITER]$V1, probs = 0.025)
quantile(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD_INC_RATE), by = ITER]$V1, probs = 0.975)

# 13/41000*1000000

# Sero testing in the population from Oct 7-31 
startTime = 18
endTime   = 274 # as.numeric(as.Date("2020-09-30") - as.Date('2019-12-31'))

# extract missed infections on and since 31 Mar
colStart = which(outputTotal$time == startTime)
colEnd   = which(outputTotal$time == endTime) # assumes missed infections at least 1 days ago can be PCR and sero pos
dataIncMissed = outputTotal$IncMissedSim_DOI[,colStart:colEnd]

# reformat 
dataIncMissed = data.table(TIME = rep(startTime:endTime, times = 10000),
                           MISSED = as.vector(t(dataIncMissed)),
                           ITER = rep(1:10000, each = 257))

dataIncMissed[, PROB_SERO_CONVERT := probSeroConvert]
dataIncMissed[, MISSED_CASES_IDENTIFIED_TEST_PERIOD := MISSED*PROB_SERO_CONVERT]
dataIncMissed[, MISSED_CASES_IDENTIFIED_TEST_PERIOD_SEROPREVALENCE_RATE := MISSED_CASES_IDENTIFIED_TEST_PERIOD/(5381000)*100]


summary(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD), by = ITER]$V1)
summary(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD_SEROPREVALENCE_RATE), by = ITER]$V1)

quantile(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD_SEROPREVALENCE_RATE), by = ITER]$V1, probs = 0.025)
quantile(dataIncMissed[, sum(MISSED_CASES_IDENTIFIED_TEST_PERIOD_SEROPREVALENCE_RATE), by = ITER]$V1, probs = 0.975)


