library(readr)
library(data.table)

# load chains in each period
listOutput = dir(path = 'output raw/20210223', pattern = '150000_timeperiod_4')
length(listOutput)

for(file in 1:length(listOutput)){
  load(paste('output raw/20210223/', listOutput[file], sep = ''))
  if(file<10) assign(paste('outputChain_0', file, sep = '') , copy(output))
  if(file>=10) assign(paste('outputChain_', file, sep = '') , copy(output))
}

rm(output)

# trace plot to check for burnin
iterN = xmax = 15000
par(mfrow=c(2,2))

for(p in 1:4){ 
  
  iteration = seq_len(iterN)
  
  chain1 = outputChain_01$param[1:xmax, p]; chain2  = outputChain_02$param[1:xmax, p]
  chain3 = outputChain_03$param[1:xmax, p]; chain4  = outputChain_04$param[1:xmax, p]
  chain5 = outputChain_05$param[1:xmax, p]; chain6  = outputChain_06$param[1:xmax, p]
  chain7 = outputChain_07$param[1:xmax, p]; chain8  = outputChain_08$param[1:xmax, p]
  chain9 = outputChain_09$param[1:xmax, p]; chain10 = outputChain_10$param[1:xmax, p]
  
  if(p %in% c(1,5, 9,13,17)) {ymin = 0; ymax = 5}
  if(p %in% c(2,6,10,14,18, 3,7,11,15,19)) {ymin = 0; ymax = 1}
  if(p %in% c(4,8,12,16,20)) {ymin = 0; ymax = 20}
  
  plot(iteration,  chain1, type = 'l', col='#00468BFF', xlab='iteration', 
       ylab=colnames(outputChain_01$param)[p], ylim = c(ymin,ymax), xlim=c(0,xmax))
  lines(iteration, chain2, col=rgb(1,165/255,0,alpha=0.5))  # "#FFA500FF"
  lines(iteration, chain3, col=rgb(237/255,0,0,alpha=0.5))
  lines(iteration, chain4, col='#00468BFF')
  lines(iteration, chain5, col=rgb(1,165/255,0,alpha=0.5))
  lines(iteration, chain6, col=rgb(237/255,0,0,alpha=0.5))
  lines(iteration, chain7, col='#00468BFF')
  lines(iteration, chain8, col=rgb(1,165/255,0,alpha=0.5))
  lines(iteration, chain9, col=rgb(237/255,0,0,alpha=0.5))
  lines(iteration, chain10, col=rgb(237/255,0,0,alpha=0.5))
  # for test outputs
  # lines(iteration, rep(trueValue, max(iteration)), col="#ED0000FF", lty = "dashed")
  
}

# trace plot for LL
iteration = seq_len(iterN)
chain1 = outputChain_01$LL[1:xmax]; chain2  = outputChain_02$LL[1:xmax]
chain3 = outputChain_03$LL[1:xmax]; chain4  = outputChain_04$LL[1:xmax]
chain5 = outputChain_05$LL[1:xmax]; chain6  = outputChain_06$LL[1:xmax]
chain7 = outputChain_07$LL[1:xmax]; chain8  = outputChain_08$LL[1:xmax]
chain9 = outputChain_09$LL[1:xmax]; chain10 = outputChain_10$LL[1:xmax]
plot(iteration, chain1, type = 'l', col='#00468BFF', xlab='iteration', ylab='LL', ylim = c(-120,-100))
lines(iteration, chain2, col=rgb(1,165/255,0,alpha=0.5))  
lines(iteration, chain3, col=rgb(237/255,0,0,alpha=0.5))  
lines(iteration, chain4, col='#00468BFF')
lines(iteration, chain5, col=rgb(1,165/255,0,alpha=0.5)) 
lines(iteration, chain6, col=rgb(237/255,0,0,alpha=0.5))
lines(iteration, chain7, col='#00468BFF')
lines(iteration, chain8, col=rgb(1,165/255,0,alpha=0.5)) 
lines(iteration, chain9, col=rgb(237/255,0,0,alpha=0.5))
lines(iteration, chain10, col=rgb(237/255,0,0,alpha=0.5))
#legend('topright', legend=c('Chain 1', 'Chain 2'), col=c('#00468BFF', rgb(1,165/255,0,alpha=0.5)), lty=1, cex=0.8)


# combine chains
timePeriod = 5
nParam = 4
outputTotal = list()

nameParam = c(paste('mu_R_period_', timePeriod, sep = ''),
              paste('effOffspringParentNotified_period_', timePeriod, sep = ''),
              paste('effOffspringParentMissed_period_', timePeriod, sep = ''),
              paste('ratioImportMissed_period_', timePeriod, sep = ''))

outputTotal = list()
outputTotal$param = matrix(NA, ncol = nParam, nrow = 1)
colnames(outputTotal$param) = nameParam

outputTotal$LL = NA
outputTotal$startTimePeriod = outputChain_01$startTimePeriod
outputTotal$endTimePeriod = outputChain_01$endTimePeriod
outputTotal$time = outputChain_01$time

outputTotal$IncNotifiedUnlinkedObs = as.numeric(outputChain_01$IncNotifiedUnlinkedObs)
outputTotal$IncNotifiedLinkedObs   = as.numeric(outputChain_01$IncNotifiedLinkedObs)

outputTotal$IncNotifiedUnlinkedSim = NA 
outputTotal$IncNotifiedLinkedSim   = NA 

outputTotal$IncNotifiedUnlinkedSim_DOI = NA
outputTotal$IncNotifiedLinkedSim_DOI   = NA
outputTotal$IncMissedUnlinkedSim_DOI   = NA
outputTotal$IncMissedLinkedSim_DOI     = NA

outputTotal$SpilloverOffspring = NA
outputTotal$SpilloverDateIsolate = NA

burnin = c(1:5000)
ls(pat = 'outputChain')  

# combine into one chain
for(c in 1:length(ls(pat = 'outputChain'))){ 
  chain = get(ls(pat = 'outputChain')[c])
  for(d in c(1, 8:15)){ outputTotal[[d]] = rbind(outputTotal[[d]], chain[[d]][-burnin,]) }
}

for(d in c(1, 8:15)){ outputTotal[[d]] = outputTotal[[d]][-1,] }

rm(chain)
outputTotalChain_period_5 = copy(outputTotal)
save(outputTotalChain_period_5, file = 'output processed/20210223/outputTotalChain_100000_iter_period_5.RData')

# combine into four chains for params only 
outputTotalChain_1_period_4 = list(param = NA)
outputTotalChain_2_period_4 = list(param = NA)
outputTotalChain_3_period_4 = list(param = NA)
outputTotalChain_4_period_4 = list(param = NA)

outputTotalChain_1_period_4$param = rbind(outputChain_01$param,outputChain_05$param[-c(1:5000),],outputChain_09$param[c(5001:77500),])   # c(5001:10000) c(5001:77500) 
outputTotalChain_2_period_4$param = rbind(outputChain_02$param,outputChain_06$param[-c(1:5000),],outputChain_09$param[c(77501:150000),]) # c(10001:15000) c(77501:150000)
outputTotalChain_3_period_4$param = rbind(outputChain_03$param,outputChain_07$param[-c(1:5000),],outputChain_10$param[c(5001:77500),])
outputTotalChain_4_period_4$param = rbind(outputChain_04$param,outputChain_08$param[-c(1:5000),],outputChain_10$param[c(77501:150000),])

# 367500 iter for chain 2-4, # 30000 iter for chain 1 & 5
# save(outputTotalChain_1_period_4, file = 'output processed/20210223/outputChain_1_period_4.RData')
# save(outputTotalChain_2_period_4, file = 'output processed/20210223/outputChain_2_period_4.RData')
# save(outputTotalChain_3_period_4, file = 'output processed/20210223/outputChain_3_period_4.RData')
# save(outputTotalChain_4_period_4, file = 'output processed/20210223/outputChain_4_period_4.RData')

# column bind different periods 
outputFourChain_period_5 = rbind(outputTotalChain_1_period_5$param,
                                 outputTotalChain_2_period_5$param,
                                 outputTotalChain_3_period_5$param,
                                 outputTotalChain_4_period_5$param)

outputFourChain_period_5 = cbind(outputFourChain_period_5, chain = rep(1:4, each = 30000)) # 30000  367500
save(outputFourChain_period_5, file = 'output processed/20210223/outputFourChain_period_5.RData')




# combine chains from all time periods
# load full chains
listOutput = dir(path = 'output processed/20210223', pattern = 'outputTotalChain')
length(listOutput)
for(file in 1:length(listOutput)){ load(paste('output processed/20210223/', listOutput[file], sep = '')) }

# load summary spillover
dataSpilloverDateIsolate = data.table(read_csv('output processed/20210223/dataSpilloverDateIsolateSummary.csv'))
dataSpilloverOffspring   = data.table(read_csv('output processed/20210223/dataSpilloverOffspringSummary.csv'))

# thin the chains
for(d in c(1, 8:15)){
  
  outputTotalChain_period_1[[d]] = outputTotalChain_period_1[[d]][seq(1,100000,10),]
  outputTotalChain_period_2[[d]] = outputTotalChain_period_2[[d]][seq(1,1450000,145),]
  outputTotalChain_period_3[[d]] = outputTotalChain_period_3[[d]][seq(1,1450000,145),]
  outputTotalChain_period_4[[d]] = outputTotalChain_period_4[[d]][seq(1,1450000,145),]
  outputTotalChain_period_5[[d]] = outputTotalChain_period_5[[d]][seq(1,100000,10),]

}


# from period 2 onwards, at the start of each period, subtract the mean spillover notified infections
# replace with actual spillover notified infections from previous time period
for(i in 1:10000){
  
  chainAfter = outputTotalChain_period_2
  chainBefore = outputTotalChain_period_1
  timeBefore = 1
  
  # notified by day of notification
  chainAfter$IncNotifiedUnlinkedSim[i,c(1:21)] = chainAfter$IncNotifiedUnlinkedSim[i,c(1:21)] - 
                                                 dataSpilloverDateIsolate[TIMEPERIOD == timeBefore, NOTIFIED_COMM_UNLINKED]
  
  chainAfter$IncNotifiedUnlinkedSim[i,c(1:21)] = chainAfter$IncNotifiedUnlinkedSim[i,c(1:21)] + 
                                                 chainBefore$SpilloverDateIsolate[i, 22:42]
  
  chainAfter$IncNotifiedLinkedSim[i,c(1:21)]   = chainAfter$IncNotifiedLinkedSim[i,c(1:21)] - 
                                                 dataSpilloverDateIsolate[TIMEPERIOD == timeBefore, NOTIFIED_COMM_LINKED]
  
  chainAfter$IncNotifiedLinkedSim[i,c(1:21)]   = chainAfter$IncNotifiedLinkedSim[i,c(1:21)] + 
                                                 chainBefore$SpilloverDateIsolate[i, 43:63]
  
  # missed by day of infection
  chainAfter$IncMissedLinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncMissedLinkedSim_DOI[i,c(1:21)] -
    (dataSpilloverOffspring[TIMEPERIOD == timeBefore & NOTIFIED_INFECTOR == 'Y', OFFSPRING]*(1-chainAfter$param[i,2]))
  
  chainAfter$IncMissedLinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncMissedLinkedSim_DOI[i,c(1:21)] +
    (chainBefore$SpilloverOffspring[i, 1:21]*(1-chainAfter$param[i,2]))
  
  chainAfter$IncMissedUnlinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncMissedUnlinkedSim_DOI[i,c(1:21)] -
    (dataSpilloverOffspring[TIMEPERIOD == timeBefore & NOTIFIED_INFECTOR == 'N', OFFSPRING]*(1-chainAfter$param[i,3]))
  
  chainAfter$IncMissedUnlinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncMissedUnlinkedSim_DOI[i,c(1:21)] +
    (chainBefore$SpilloverOffspring[i, 22:42]*(1-chainAfter$param[i,3]))
  
  # notified by day of infection
  chainAfter$IncNotifiedLinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncNotifiedLinkedSim_DOI[i,c(1:21)] -
    (dataSpilloverOffspring[TIMEPERIOD == timeBefore & NOTIFIED_INFECTOR == 'Y', OFFSPRING]*chainAfter$param[i,2])
  
  chainAfter$IncNotifiedLinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncNotifiedLinkedSim_DOI[i,c(1:21)] +
    (chainBefore$SpilloverOffspring[i, 1:21]*chainAfter$param[i,2])
  
  chainAfter$IncNotifiedUnlinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncNotifiedUnlinkedSim_DOI[i,c(1:21)] -
    (dataSpilloverOffspring[TIMEPERIOD == timeBefore & NOTIFIED_INFECTOR == 'N', OFFSPRING]*chainAfter$param[i,3])
  
  chainAfter$IncNotifiedUnlinkedSim_DOI[i,c(1:21)] = 
    chainAfter$IncNotifiedUnlinkedSim_DOI[i,c(1:21)] +
    (chainBefore$SpilloverOffspring[i, 22:42]*chainAfter$param[i,3])
  
  outputTotalChain_period_2 = chainAfter
  
  
}
  
rm(chainAfter, chainBefore, timeBefore)


# bind all time periods together
outputTotalChain_period_all = list()

outputTotalChain_period_all$param = cbind(outputTotalChain_period_1$param,
                                          outputTotalChain_period_2$param,
                                          outputTotalChain_period_3$param,
                                          outputTotalChain_period_4$param,
                                          outputTotalChain_period_5$param)

outputTotalChain_period_all$startTimePeriod = c(outputTotalChain_period_1$startTimePeriod,
                                                outputTotalChain_period_2$startTimePeriod,
                                                outputTotalChain_period_3$startTimePeriod,
                                                outputTotalChain_period_4$startTimePeriod,
                                                outputTotalChain_period_5$startTimePeriod)

outputTotalChain_period_all$endTimePeriod = c(outputTotalChain_period_1$endTimePeriod,
                                              outputTotalChain_period_2$endTimePeriod,
                                              outputTotalChain_period_3$endTimePeriod,
                                              outputTotalChain_period_4$endTimePeriod,
                                              outputTotalChain_period_5$endTimePeriod)

outputTotalChain_period_all$time = c(outputTotalChain_period_1$time,
                                     outputTotalChain_period_2$time,
                                     outputTotalChain_period_3$time,
                                     outputTotalChain_period_4$time,
                                     outputTotalChain_period_5$time)

outputTotalChain_period_all$IncNotifiedUnlinkedObs = c(outputTotalChain_period_1$IncNotifiedUnlinkedObs,
                                                       outputTotalChain_period_2$IncNotifiedUnlinkedObs,
                                                       outputTotalChain_period_3$IncNotifiedUnlinkedObs,
                                                       outputTotalChain_period_4$IncNotifiedUnlinkedObs,
                                                       outputTotalChain_period_5$IncNotifiedUnlinkedObs)

outputTotalChain_period_all$IncNotifiedLinkedObs = c(outputTotalChain_period_1$IncNotifiedLinkedObs,
                                                     outputTotalChain_period_2$IncNotifiedLinkedObs,
                                                     outputTotalChain_period_3$IncNotifiedLinkedObs,
                                                     outputTotalChain_period_4$IncNotifiedLinkedObs,
                                                     outputTotalChain_period_5$IncNotifiedLinkedObs)

outputTotalChain_period_all$IncNotifiedUnlinkedSim = cbind(outputTotalChain_period_1$IncNotifiedUnlinkedSim,
                                                           outputTotalChain_period_2$IncNotifiedUnlinkedSim,
                                                           outputTotalChain_period_3$IncNotifiedUnlinkedSim,
                                                           outputTotalChain_period_4$IncNotifiedUnlinkedSim,
                                                           outputTotalChain_period_5$IncNotifiedUnlinkedSim)

outputTotalChain_period_all$IncNotifiedLinkedSim = cbind(outputTotalChain_period_1$IncNotifiedLinkedSim,
                                                         outputTotalChain_period_2$IncNotifiedLinkedSim,
                                                         outputTotalChain_period_3$IncNotifiedLinkedSim,
                                                         outputTotalChain_period_4$IncNotifiedLinkedSim,
                                                         outputTotalChain_period_5$IncNotifiedLinkedSim)

outputTotalChain_period_all$IncNotifiedUnlinkedSim_DOI = cbind(outputTotalChain_period_1$IncNotifiedUnlinkedSim_DOI,
                                                               outputTotalChain_period_2$IncNotifiedUnlinkedSim_DOI,
                                                               outputTotalChain_period_3$IncNotifiedUnlinkedSim_DOI,
                                                               outputTotalChain_period_4$IncNotifiedUnlinkedSim_DOI,
                                                               outputTotalChain_period_5$IncNotifiedUnlinkedSim_DOI)

outputTotalChain_period_all$IncNotifiedLinkedSim_DOI = cbind(outputTotalChain_period_1$IncNotifiedLinkedSim_DOI,
                                                             outputTotalChain_period_2$IncNotifiedLinkedSim_DOI,
                                                             outputTotalChain_period_3$IncNotifiedLinkedSim_DOI,
                                                             outputTotalChain_period_4$IncNotifiedLinkedSim_DOI,
                                                             outputTotalChain_period_5$IncNotifiedLinkedSim_DOI)

outputTotalChain_period_all$IncMissedUnlinkedSim_DOI = cbind(outputTotalChain_period_1$IncMissedUnlinkedSim_DOI,
                                                             outputTotalChain_period_2$IncMissedUnlinkedSim_DOI,
                                                             outputTotalChain_period_3$IncMissedUnlinkedSim_DOI,
                                                             outputTotalChain_period_4$IncMissedUnlinkedSim_DOI,
                                                             outputTotalChain_period_5$IncMissedUnlinkedSim_DOI)

outputTotalChain_period_all$IncMissedLinkedSim_DOI = cbind(outputTotalChain_period_1$IncMissedLinkedSim_DOI,
                                                           outputTotalChain_period_2$IncMissedLinkedSim_DOI,
                                                           outputTotalChain_period_3$IncMissedLinkedSim_DOI,
                                                           outputTotalChain_period_4$IncMissedLinkedSim_DOI,
                                                           outputTotalChain_period_5$IncMissedLinkedSim_DOI)

outputTotalChain_period_all$IncMissedSim_DOI = outputTotalChain_period_all$IncMissedUnlinkedSim_DOI + outputTotalChain_period_all$IncMissedLinkedSim_DOI 


outputTotalChain_period_all$SpilloverOffspring = cbind(outputTotalChain_period_1$SpilloverOffspring,
                                                       outputTotalChain_period_2$SpilloverOffspring,
                                                       outputTotalChain_period_3$SpilloverOffspring,
                                                       outputTotalChain_period_4$SpilloverOffspring,
                                                       outputTotalChain_period_5$SpilloverOffspring)

outputTotalChain_period_all$SpilloverDateIsolate = cbind(outputTotalChain_period_1$SpilloverDateIsolate,
                                                         outputTotalChain_period_2$SpilloverDateIsolate,
                                                         outputTotalChain_period_3$SpilloverDateIsolate,
                                                         outputTotalChain_period_4$SpilloverDateIsolate,
                                                         outputTotalChain_period_5$SpilloverDateIsolate)

save(outputTotalChain_period_all, file = 'output processed/20210223/outputTotalChain_10000_thinned_iter_period_all.RData')


