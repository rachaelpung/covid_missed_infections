library(mgcv)
library(MASS)
library(readr)
library(data.table)
library(doParallel)

setwd('C:/Users/rachaelpung/Desktop/Secured thumb/Projects/COVID-19/Missed infections/')
source('codes/v5/covid_missInf_functions_simulations.r')

# load distribution of generation interval over time
distGen = data.table(read_csv('input/dist_Gen_Interval.csv')) # gamma distributed, mu 7.5, sigma 3.4

# load distribution of duration from infection to isolation
# distIsolate = data.table(read_csv('input/dist_Isolate.csv')) 

# generate distribution of duration from infection to isolation
t = 0:21
distIsolate = dgamma(t,shape=25,rate=2.5)
distIsolate = distIsolate/sum(distIsolate)
distIsolate = data.table(DAY = t, PROB_ISOLATE = distIsolate)
distIsolate[, CUM_PROB_ISOLATE := cumsum(distIsolate$PROB_ISOLATE)]
distIsolate[, PROB_AT_LARGE := 1-CUM_PROB_ISOLATE]
distIsolate = rbind(distIsolate,distIsolate)
distIsolate[, PERIOD := rep(seq_len(2), each = 22)]

# load distribution of time of infection for imported cases
distTravel = data.table(read_csv('input/dist_Travel.csv')) 

# load initalParam
initialParam = read_csv('input/initialParamTest.csv')

# generate linked and unlinked cases over time based on known parameters
dataTest = data.table(CASE = 1:300,
                      NOTIFIED = rep('Y', times = 300),
                      LOCAL_OVERSEAS = 1,
                      SOURCE_OF_INFECTION = 1,
                      DATE_SHN_OR_ISOLATED_FOR_TESTING = c(seq_len(100), rep(NA, times = 50), seq(51,100,1), rep(NA, times = 50), seq(51,100,1)),
                      DATE_ARRIVAL = rep(seq_len(100), times = 3),
                      DATE_NOTIFICATION = rep(seq_len(100), times = 3) + 5,
                      DATE_ONSET = c(rep(seq_len(100)+3, times = 2), rep(NA, times = 100)))

dataTest[, DATE_ISOLATE_FROM_COMM := DATE_NOTIFICATION]

# set up time periods
timePeriod = 1 # period of interest
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

# load distribution of spillover infections in previous time period
dataSpilloverOffspring = data.table()
dataSpilloverDateIsolate = data.table()

# subset imported test cases with community exposure 
dataCaseImportNotified = dataTest[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING)]
dataCaseImportNotified = dataCaseImportNotified[, .(NOTIFIED, DATE_ARRIVAL, SOURCE_OF_INFECTION, 
                                                    LOCAL_OVERSEAS, DATE_ISOLATE_FROM_COMM)]
dataCaseImportNotified[, CASE_COUNT := 1]
dataCaseImportNotified = dataCaseImportNotified[DATE_ARRIVAL %in% testParam$time]

# imported cases
dataIncImport   = data.table(
  TIME =  2:100, 
  NOTIFIED_IMPORT_COMM_EXPOSURE = hist(dataTest[LOCAL_OVERSEAS == 1 & is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING), 
                                                DATE_ARRIVAL], breaks = 1:100, plot = F)$counts,
  NOTIFIED_IMPORT_ISOLATED      = hist(dataTest[LOCAL_OVERSEAS == 1 & !is.na(DATE_SHN_OR_ISOLATED_FOR_TESTING), 
                                                DATE_ARRIVAL], breaks = 1:100, plot = F)$counts
)

dataIncImport[1, NOTIFIED_IMPORT_COMM_EXPOSURE := dataIncImport[1, NOTIFIED_IMPORT_COMM_EXPOSURE] - 2]
dataIncImport[1, NOTIFIED_IMPORT_ISOLATED := dataIncImport[1, NOTIFIED_IMPORT_ISOLATED] - 1]

dataIncImport = rbindlist(list(data.table(TIME = 1, NOTIFIED_IMPORT_COMM_EXPOSURE = 2, NOTIFIED_IMPORT_ISOLATED = 1), 
                               dataIncImport), use.names = T, fill = T)

# fit smoothing spline for all notified imported cases with community exposure
dataIncImport[, SPLINE := gam(NOTIFIED_IMPORT_COMM_EXPOSURE~s(TIME), data = dataIncImport, family='poisson')$fitted.values]
dataIncImport = dataIncImport[TIME %in% testParam$time]

distIsolate = distIsolate[PERIOD == 1,]

# simulate daily incidence, run the simCommInc function line by line
dataInc = simCommInc(dataCaseImportNotified, dataIncImport, dataSpilloverOffspring, dataSpilloverDateIsolate,
                     testParam, distGen, distIsolate, distTravel)

dataInc = copy(dataInc$currentIncIsolate)

dataInc[, `:=` (NOTIFIED_COMM_UNLINKED = round(NOTIFIED_COMM_UNLINKED, digits = 0),
                NOTIFIED_COMM_LINKED = round(NOTIFIED_COMM_LINKED, digits = 0))]

# start MCMC runs
setSeed     = scan('input/setSeed.txt')

paramSet = 1

## set up clusters
# cl = makeCluster(2) #detectCores()
# registerDoParallel(cl)

# foreach(paramSet = 1:2, .packages = c("data.table", "mgcv"), .combine = "c") %dopar%  {
#nrow(initalParam)

# attempt to recover parameters used to simulate the above dataset
set.seed(setSeed[paramSet])
  
iterN = 5000 # 100000
timePeriod = 1 # no. of time intervals in the dataset
nParam = 4

startTimePeriod = c(1,51)  
endTimePeriod =  c(50,100)


nameParam = character()
for(period in 1:timePeriod){
  
  name = c(paste('mu_R_period_', period, sep = ''),
           paste('effOffspringParentNotified_period_', period, sep = ''),
           paste('effOffspringParentMissed_period_', period, sep = ''),
           paste('ratioImportMissed_period_', period, sep = ''))#,
           #paste('mu_R_missed_period_', period, sep = ''))
  
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

# distribution of generation interval over time
distGen 

# distribution of duration from infection to isolation
distIsolate 

# distribution of time of infection for imported cases
distTravel 

# imported cases
dataCaseImportNotified
dataIncImport

# simulate expected daily incidence given current parameter
currentInc = simCommInc(dataCaseImportNotified, dataIncImport, dataSpilloverOffspring, dataSpilloverDateIsolate,
                        currentParam, distGen, distIsolate, distTravel)

# calculate log likelihood
currentParam$LL = calculateLogLik(dataInc, currentInc$currentIncIsolate, currentParam$lenTime)

# store currentParam and timeSeriesInc into output
output$param[1,] = currentParam$param 
output$LL[1]    = currentParam$LL 

output$IncNotifiedUnlinkedObs[1,] = dataInc$NOTIFIED_COMM_UNLINKED
output$IncNotifiedLinkedObs[1,]   = dataInc$NOTIFIED_COMM_LINKED

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
for(iter in 1001:5000){ #iterN
  
  print(sprintf("Iteration %s of %.f. Acceptance rate %.03f", iter, iterN, acceptRate))
  
  # MH for current transmission parameter
  oldParam               = currentParam
  oldInc                 = currentInc
  
  currentParam$param[1,1:nParam] = mvrnorm(1, currentParam$param[1,1:nParam], Sigma=lambda*covMatrix)

  iterOutput             = mh(oldParam,currentParam, 
                              dataInc, oldInc,
                              dataCaseImportNotified, dataIncImport, dataSpilloverOffspring, dataSpilloverDateIsolate)
  
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
    save(output,file= paste0("output raw/outputTest_", paramSet, "_iter_", iter, "_timeperiod_", timePeriod,  ".rdata"))
  }
  

}

#acceptedN/iterN
#accept_rateT


#save(output,file= paste0("output/outputTest_", paramSet, ".rdata"))

paramSet

# }

## stop clusters
# stopCluster(cl)
