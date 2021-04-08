library(data.table)

load('output processed/20210223/outputTotalChain_100000_iter_period_5.RData')

outputTotal = copy(outputTotalChain_period_5)
timePeriod = 5

dataSpilloverDateIsolateSummary = apply(outputTotal$SpilloverDateIsolate[seq(1,100000,100),],2, mean, na.rm = T)
dataSpilloverDateIsolateSummary = data.table(TIME = dataSpilloverDateIsolateSummary[1:21],
                                             NOTIFIED_COMM_UNLINKED = dataSpilloverDateIsolateSummary[22:42],
                                             NOTIFIED_COMM_LINKED = dataSpilloverDateIsolateSummary[43:63],
                                             TIMEPERIOD = timePeriod)

dataSpilloverOffspringSummary = apply(outputTotal$SpilloverOffspring[seq(1,100000,100),],2, mean, na.rm = T)  
dataSpilloverOffspringSummary = data.table(DAY_OF_INFECTION = rep(dataSpilloverDateIsolateSummary[,TIME], 2),
                                           NOTIFIED_INFECTOR = rep(c('Y', 'N'), each = 21),
                                           OFFSPRING = dataSpilloverOffspringSummary,
                                           TIMEPERIOD = timePeriod)

write.csv(dataSpilloverDateIsolateSummary, file = 'output processed/20210223/dataSpilloverDateIsolateSummary_tmp.csv', row.names = F)
write.csv(dataSpilloverOffspringSummary, file = 'output processed/20210223/dataSpilloverOffspringSummary_tmp.csv', row.names = F)
