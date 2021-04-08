library(mgcv)
library(MASS)
library(readr)
library(data.table)
library(doParallel)

source('codes/v5/MCMC/covid_missInf_functions_simulations.r')

# generate imported, isolated cases over two time periods, each lasting 50 days
# in period 1, there were 50 cases isolated upon arrival, 100 cases isolated few days after arrival 
# in period 2, there were 150 cases isolated upon arrival 
dataTest = data.table(NOTIFIED = rep('Y', times = 300),
                      LOCAL_OVERSEAS = 1,
                      SOURCE_OF_INFECTION = 1,
                      DATE_ARRIVAL = c(rep(seq(1,50,1), 3),rep(seq(51,100,1), 3)) )

dataTest[c(1:50, 151:300), DATE_ISOLATE_FROM_COMM := DATE_ARRIVAL]
dataTest[c(1:50, 151:300), DATE_SHN_OR_ISOLATED_FOR_TESTING := DATE_ARRIVAL]
dataTest[51:150, DATE_ISOLATE_FROM_COMM := DATE_ARRIVAL + 5]

# subset imported test cases with community exposure 
testCaseImportNotified = dataTest[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING)]
testCaseImportNotified = dataCaseImportNotified[, .(NOTIFIED, DATE_ARRIVAL, SOURCE_OF_INFECTION, 
                                                    LOCAL_OVERSEAS, DATE_ISOLATE_FROM_COMM)]
testCaseImportNotified[, CASE_COUNT := 1]

# save
write.csv(testCaseImportNotified, file = 'test/input/testCaseImportNotified.csv', row.names = F)

# tabulate imported cases by day of arrival
testIncImport   = data.table(
  TIME =  2:100, 
  NOTIFIED_IMPORT_COMM_EXPOSURE = hist(dataTest[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING), 
                                                DATE_ARRIVAL], breaks = 1:100, plot = F)$counts,
  NOTIFIED_IMPORT_ISOLATED      = hist(dataTest[LOCAL_OVERSEAS == 1 & !is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING), 
                                                DATE_ARRIVAL], breaks = 1:100, plot = F)$counts
)

dataTest[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING) & DATE_ARRIVAL == 1,.N]
dataTest[LOCAL_OVERSEAS == 1 & !is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING) & DATE_ARRIVAL == 1,.N]

testIncImport[1, NOTIFIED_IMPORT_COMM_EXPOSURE := testIncImport[1, NOTIFIED_IMPORT_COMM_EXPOSURE] - 2]
testIncImport[1, NOTIFIED_IMPORT_ISOLATED := testIncImport[1, NOTIFIED_IMPORT_ISOLATED] - 1]
testIncImport = rbindlist(list(data.table(TIME = 1, NOTIFIED_IMPORT_COMM_EXPOSURE = 0, NOTIFIED_IMPORT_ISOLATED = 1), 
                               testIncImport), use.names = T, fill = T)

# fit smoothing spline for all notified imported cases with community exposure
testIncImport[, SPLINE := gam(NOTIFIED_IMPORT_COMM_EXPOSURE~s(TIME), data = testIncImport, family='poisson')$fitted.values]

# save
write.csv(testIncImport, file = 'test/input/testIncImport.csv', row.names = F)

## generate test outbreak 

# load distribution of generation interval over time
distGen = data.table(read_csv('test/input/dist_Gen_Interval.csv')) # gamma distributed, mu 7.5, sigma 3.4

# load distribution of duration from infection to isolation
distIsolate = data.table(read_csv('test/input/dist_Isolate.csv'))

# load distribution of time of infection for imported cases
distTravel = data.table(read_csv('test/input/dist_Travel.csv')) 

# load initalParam
initialParam = read_csv('test/input/initialParamTest.csv')
testSpilloverOffspring = data.table()
testSpilloverDateIsolate = data.table()

# load observed data
testCaseImportNotified = data.table(read_csv('test/input/testCaseImportNotified.csv')) 
testIncImport = data.table(read_csv('test/input/testIncImport.csv')) 

testInc = data.table()

# set up time periods
for(t in 1:2){
timePeriod = t
startTimePeriod = c(1,51)  
endTimePeriod =  c(50,100)

# set up parameters
nParam = 4 # no of param for the period of interest
paramSet = 3
testParam = list()
testParam$param = as.matrix(initialParam[paramSet,((timePeriod-1)*4+1):((timePeriod-1)*4+nParam)])
testParam$LL = NA
testParam$time = startTimePeriod[timePeriod]:endTimePeriod[timePeriod]
testParam$lenTime = length(testParam$time)
testParam$startTimePeriod = startTimePeriod[timePeriod]
testParam$endTimePeriod = endTimePeriod[timePeriod]

# imported cases
testCaseImportNotifiedPeriod = testCaseImportNotified[DATE_ARRIVAL %in% testParam$time]
testIncImportPeriod = testIncImport[TIME %in% testParam$time]

distIsolatePeriod = distIsolate[PERIOD == timePeriod,]

# simulate daily incidence, run the simCommInc function line by line
testTmpInc = simCommInc(testCaseImportNotifiedPeriod, testIncImportPeriod, testSpilloverOffspring, testSpilloverDateIsolate,
                        testParam, distGen, distIsolatePeriod, distTravel)

testSpilloverOffspring = testTmpInc$spilloverOffspring
testSpilloverDateIsolate = testTmpInc$spilloverDateIsolate
  
testTmpInc = copy(testTmpInc$currentIncIsolate)

testTmpInc[, `:=` (NOTIFIED_COMM_UNLINKED = round(NOTIFIED_COMM_UNLINKED, digits = 0),
                   NOTIFIED_COMM_LINKED = round(NOTIFIED_COMM_LINKED, digits = 0))]

testInc = rbind(testInc, testTmpInc)

}

write.csv(testInc, 'test/input/testInc.csv', row.names = F)

# attempt to recover parameter values in respective time period
# load distribution of generation interval over time
distGen = data.table(read_csv('test/input/dist_Gen_Interval.csv')) # gamma distributed, mu 7.5, sigma 3.4

# load distribution of duration from infection to isolation
distIsolate = data.table(read_csv('test/input/dist_Isolate.csv'))

# load distribution of time of infection for imported cases
distTravel = data.table(read_csv('test/input/dist_Travel.csv')) 

# load observed data
testCaseImportNotified = data.table(read_csv('test/input/testCaseImportNotified.csv')) 
testIncImport = data.table(read_csv('test/input/testIncImport.csv')) 
testInc = data.table(read_csv('test/input/testInc.csv')) 

# spillover cases
initialParam = read_csv('test/input/initialParamTest.csv')
testSpilloverOffspring = data.table() # for time period 1
testSpilloverDateIsolate = data.table() # for time period 1

# start MCMC runs
setSeed     = scan('input/setSeed.txt')

paramSet = 1

## set up clusters
# cl = makeCluster(2) #detectCores()
# registerDoParallel(cl)

# foreach(paramSet = 1:2, .packages = c("data.table", "mgcv"), .combine = "c") %dopar%  {
#nrow(initalParam)
set.seed(setSeed[paramSet])
  
iterN = 1000 # 100000
timePeriod = 1 # no. of time intervals in the dataset
nParam = 4

startTimePeriod = c(1,51)  
endTimePeriod =  c(50,100)


nameParam = character()
for(period in 1:timePeriod){
  
  name = c(paste('mu_R_period_', period, sep = ''),
           paste('effOffspringParentNotified_period_', period, sep = ''),
           paste('effOffspringParentMissed_period_', period, sep = ''),
           paste('ratioImportMissed_period_', period, sep = ''))
  
  nameParam = c(nameParam, name)
  
}

output = list()
output$param = matrix(0, iterN, nParam)
colnames(output$param) = nameParam

output$LL = rep(0,iterN)
output$startTimePeriod = startTimePeriod[timePeriod]  
output$endTimePeriod = endTimePeriod[timePeriod]
output$time = startTimePeriod[timePeriod]:endTimePeriod[timePeriod]

output$IncNotifiedUnlinkedObs = matrix(0, 1, length(output$time)) 
output$IncNotifiedLinkedObs   = matrix(0, 1, length(output$time)) 

# store output by day of isolation
output$IncNotifiedUnlinkedSim = matrix(0, iterN, length(output$time)) 
output$IncNotifiedLinkedSim   = matrix(0, iterN, length(output$time)) 

# store output by day of infection
output$IncNotifiedUnlinkedSim_DOI = matrix(0, iterN, length(output$time)) 
output$IncNotifiedLinkedSim_DOI   = matrix(0, iterN, length(output$time))
output$IncMissedUnlinkedSim_DOI   = matrix(0, iterN, length(output$time))
output$IncMissedLinkedSim_DOI     = matrix(0, iterN, length(output$time))

output$SpilloverOffspring   = matrix(0, iterN, 2*21)
output$SpilloverDateIsolate = matrix(0, iterN, 3*21)

# initialise parameters
currentParam = list()
currentParam$param = as.matrix(initialParam[paramSet,((timePeriod-1)*4+1):((timePeriod-1)*4+nParam)])
currentParam$LL = NA
currentParam$time = output$time
currentParam$lenTime = length(output$time)
currentParam$startTimePeriod = output$startTimePeriod
currentParam$endTimePeriod = output$endTimePeriod


testCaseImportNotified = testCaseImportNotified[DATE_ARRIVAL %in% currentParam$time]
testIncImport = testIncImport[TIME %in% currentParam$time]
distIsolate = distIsolate[PERIOD == timePeriod,]

testInc = testInc[TIME %in% currentParam$time]

# simulate expected daily incidence given current parameter
currentInc = simCommInc(testCaseImportNotified, testIncImport, testSpilloverOffspring, testSpilloverDateIsolate,
                        currentParam, distGen, distIsolate, distTravel)



# calculate log likelihood
currentParam$LL = calculateLogLik(testInc, currentInc$currentIncIsolate, currentParam$lenTime)

# store currentParam and timeSeriesInc into output
output$param[1,] = currentParam$param 
output$LL[1]    = currentParam$LL 

output$IncNotifiedUnlinkedObs[1,] = testInc$NOTIFIED_COMM_UNLINKED
output$IncNotifiedLinkedObs[1,]   = testInc$NOTIFIED_COMM_LINKED

output$IncNotifiedUnlinkedSim[1,] = currentInc$currentIncIsolate$NOTIFIED_COMM_UNLINKED
output$IncNotifiedLinkedSim[1,]   = currentInc$currentIncIsolate$NOTIFIED_COMM_LINKED

output$IncNotifiedUnlinkedSim_DOI[1,] = currentInc$currentIncInfection$NOTIFIED_COMM_UNLINKED
output$IncNotifiedLinkedSim_DOI[1,]   = currentInc$currentIncInfection$NOTIFIED_COMM_LINKED
output$IncMissedUnlinkedSim_DOI[1,]   = currentInc$currentIncInfection$MISSED_COMM_UNLINKED
output$IncMissedLinkedSim_DOI[1,]     = currentInc$currentIncInfection$MISSED_COMM_LINKED

output$SpilloverOffspring[1,]     = currentInc$spilloverOffspring$OFFSPRING
output$SpilloverDateIsolate[1,]   = c(currentInc$spilloverDateIsolate$TIME, currentInc$spilloverDateIsolate$NOTIFIED_COMM_UNLINKED, currentInc$spilloverDateIsolate$NOTIFIED_COMM_LINKED)


lambda = 0.01
sd = rep(1,nParam); covMatrix = diag(sd)
acceptRate = 0.234 ;  acceptedN = 0

# start MCMC
for(iter in 501:iterN){ #iterN
  
  print(sprintf("Iteration %s of %.f. Acceptance rate %.03f", iter, iterN, acceptRate))
  
  # MH for current transmission parameter
  oldParam               = currentParam
  oldInc                 = currentInc
  
  currentParam$param[1,1:nParam] = mvrnorm(1, currentParam$param[1,1:nParam], Sigma=lambda*covMatrix)

  iterOutput             = mh(oldParam,currentParam, 
                              testInc, oldInc,
                              testCaseImportNotified, testIncImport, testSpilloverOffspring, testSpilloverDateIsolate)
  
  currentParam      =  iterOutput[[1]]
  acceptedN         =  acceptedN + iterOutput[[2]]
  currentInc        =  iterOutput[[3]]  
    
    
  if(iter >= 100){
    # start adaptation after 100 iterations
    acceptRate = acceptedN/(iter-1)
    lambda = max(0.01,min(1,exp(log(lambda)+(acceptRate-0.234)*0.999^iter))) 
    covMatrix  = cov(output$param[2:iter,1:nParam])
    
  }
  
    
  output$param[iter,] = currentParam$param 
  output$LL[iter]    = currentParam$LL 
  
  output$IncNotifiedUnlinkedSim[iter,] = currentInc$currentIncIsolate$NOTIFIED_COMM_UNLINKED
  output$IncNotifiedLinkedSim[iter,]   = currentInc$currentIncIsolate$NOTIFIED_COMM_LINKED
  
  output$IncNotifiedUnlinkedSim_DOI[iter,] = currentInc$currentIncInfection$NOTIFIED_COMM_UNLINKED
  output$IncNotifiedLinkedSim_DOI[iter,]   = currentInc$currentIncInfection$NOTIFIED_COMM_LINKED
  output$IncMissedUnlinkedSim_DOI[iter,]   = currentInc$currentIncInfection$MISSED_COMM_UNLINKED
  output$IncMissedLinkedSim_DOI[iter,]     = currentInc$currentIncInfection$MISSED_COMM_LINKED
  
  output$SpilloverOffspring[iter,]     = currentInc$spilloverOffspring$OFFSPRING
  output$SpilloverDateIsolate[iter,]   = c(currentInc$spilloverDateIsolate$TIME, currentInc$spilloverDateIsolate$NOTIFIED_COMM_UNLINKED, currentInc$spilloverDateIsolate$NOTIFIED_COMM_LINKED)
  
  
  if(iter %% 1000 == 0){
    save(output,file= paste0("test/output raw/outputTest_", paramSet, "_iter_", iter, "_timeperiod_", timePeriod,  ".rdata"))
  }
  

}

#acceptedN/iterN
#accept_rateT


#save(output,file= paste0("test/output raw/outputTest_", paramSet, ".rdata"))

paramSet

# }

## stop clusters
# stopCluster(cl)
